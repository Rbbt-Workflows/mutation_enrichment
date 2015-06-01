module MutationEnrichment
  dep do |jobname, inputs| job(inputs[:baseline] || :pathway_base_counts, inputs[:database].to_s, inputs) end
  input :database, :select, "Database code", nil, :select_options => DATABASES
  input :baseline, :select, "Type of baseline to use", :pathway_base_counts, :select_options => [:pathway_base_counts, :pathway_gene_counts]
  input :mutations, :array, "Genomic Mutation"
  input :fdr, :boolean, "BH FDR corrections", true
  input :masked_genes, :array, "Ensembl Gene ID list of genes to mask", []
  input :background, :array, "Enrichment background", nil
  input :invert_background, :boolean, "Restrict to elements NOT in background"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :watson, :boolean, "Alleles reported in the watson (forward) strand", true
  task :mutation_pathway_enrichment => :tsv do |database,baseline,mutations,fdr,masked_genes,background,invert_background,organism, watson|
    counts        = step(baseline).load
    total_covered = step(baseline).info[:total_size] || step(baseline).info[:total_genes]
    GenomicMutation.setup(mutations, "MutationEnrichment", organism, watson)


    masked_genes = background if background and invert_background

    affected_genes = mutations.genes.compact.flatten.uniq
    affected_genes = affected_genes.remove(masked_genes) if masked_genes and masked_genes.any?
    affected_genes = affected_genes.subset(background) if background and not background.empty? and not invert_background

    # Get database tsv and native ids

    database_tsv, all_db_genes, db_gene_field, db_pathway_field = database_info database, organism

    all_db_genes = all_db_genes.ensembl
    all_db_genes = all_db_genes.remove(masked_genes) if masked_genes and masked_genes.any?
    all_db_genes = all_db_genes.subset(background) if background and not background.empty? and not invert_background
    all_db_genes.compact!
    all_db_genes.sort!


    database_tsv = database_tsv.reorder db_gene_field, [db_pathway_field]


    # Annotate each pathway with the affected genes that are involved in it

    log :pathway_matches, "Finding affected genes per pathway"
    affected_genes_per_pathway = {}
    database_tsv.with_unnamed do
      affected_genes.each do |gene|
        next if gene.nil? or gene.empty?
        next unless database_tsv[gene]
        pathways = database_tsv[gene].flatten
        next if pathways.nil? or pathways.empty?
        pathways.uniq.each do |pathway|
          next if pathway.nil? or pathway.empty?
          affected_genes_per_pathway[pathway] ||= []
          affected_genes_per_pathway[pathway] << gene
        end
      end
    end

    log :mutation_genes, "Finding genes overlapping mutations"
    mutation_genes = {}
    gene_mutations = {}
    mutations.genes.zip(mutations.clean_annotations).each do |genes, mutation|
      next if genes.nil?
      mutation_genes[mutation] = genes.sort
      genes.each do |gene|
        gene_mutations[gene] ||= []
        gene_mutations[gene] << mutation
      end
    end
    mutations = mutations.clean_annotations

    log :covered_mutations, "Finding mutations overlapping genes in pathway"
    all_db_genes = all_db_genes.ensembl

    covered_mutations = mutations.select{|mutation| Misc.intersect_sorted_arrays((mutation_genes[mutation] || []).dup, all_db_genes.dup).any? }.length
    set_info :covered_mutations, covered_mutations

    log :pvalue, "Calculating binomial pvalues"
    pvalues = TSV.setup({}, :key_field => database_tsv.fields.first, :fields => ["Matches", "Pathway total", "Ensembl Gene ID"], :namespace => organism, :type => :double)
    counts.unnamed = true
    affected_genes_per_pathway.each do |pathway, genes|
      next if pathway.nil? or pathway.empty?
      pathway_total = counts[pathway]
      matches = gene_mutations.values_at(*genes).compact.flatten.length
      next if matches == 0

      common_genes = affected_genes.subset(genes).uniq
      pvalues[pathway] = [[matches], [pathway_total], common_genes.sort_by{|g| g.name || g}]
    end

    if pvalues.size == 0
      pvalues
    else
      pvalues = pvalues.R("pvalues = apply(data, 1, function(v){ binom.test(as.numeric(v[1]), #{ covered_mutations }, as.numeric(v[2]) /  #{total_covered.to_f}, 'greater')$p.value });
      data = cbind(data, p.value = pvalues);
      data = data[names(data)[c(1,2,4,3)]];", :key => pvalues.key_field) unless pvalues.empty?

      pvalues.process "p.value" do |v|
        v.to_f
      end

      pvalues = FDR.adjust_hash! pvalues, 2 if fdr and not pvalues.empty?

      pvalues.type = :double

      set_info :total_covered, total_covered

      pvalues
    end
  end
  export_asynchronous :mutation_pathway_enrichment
end
