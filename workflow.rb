require 'rbbt'
require 'rbbt/util/misc'
require 'rbbt/tsv/change_id'
require 'rbbt/workflow'
require 'rbbt/statistics/hypergeometric'
require 'rbbt/statistics/random_walk'

Workflow.require_workflow 'Genomics'
require 'rbbt/knowledge_base/Genomics'
require 'rbbt/entity/gene'
require 'rbbt/entity/genomic_mutation'

module MutationEnrichment
  extend Workflow
  extend Resource
  
  MASKED_TERMS = %w(cancer melanoma carcinoma glioma hepatitis leukemia leukaemia disease infection opathy hepatitis sclerosis hepatatis glioma Shigellosis plasmosis maniasis)
  MASKED_IDS = {}

  self.subdir = "var/MutationEnrichment"

  class << self
    attr_accessor :knowledge_base_dir

    def knowledge_base_dir
      @knowledge_base_dir ||= MutationEnrichment.var.knowledge_base
    end
  end

  DATABASES = Genomics.knowledge_base.registry.keys

  helper :database_info do |database, organism|
    @organism_kb ||= {}
    @organism_kb[organism] ||= begin
                                 dir = MutationEnrichment.knowledge_base_dir

                                 kb = KnowledgeBase.new dir, organism
                                 kb.format["Gene"] = "Ensembl Gene ID"
                                 kb.registry = Genomics.knowledge_base.registry
                                 kb
                               end

    db = @organism_kb[organism].get_database(database, :persist => true)

    tsv, total_keys, source_field, target_field = [db, db.keys, db.key_field, db.fields.first]

    if target_field == "Ensembl Gene ID"
      pathway_field, gene_field = source_field, target_field
      total_genes = Gene.setup(tsv.values.flatten.compact.uniq, "Ensembl Gene ID", organism)
    else
      pathway_field, gene_field = target_field, source_field
      total_genes = Gene.setup(total_keys, gene_field, organism)
    end

    tsv.namespace = organism

    [tsv, total_genes, gene_field, pathway_field]
  end
  
  helper :database_info do |database, organism|
    Persist.memory([database, organism] * ": ") do
      @organism_kb ||= {}
      @organism_kb[organism] ||= begin
                                  dir = MutationEnrichment.knowledge_base_dir

                                  kb = KnowledgeBase.new dir, organism
                                  kb.format["Gene"] = "Ensembl Gene ID"
                                  kb.registry = Genomics.knowledge_base.registry
                                  kb
                                end

      db = @organism_kb[organism].get_database(database, :persist => true, :target => "Gene=>Ensembl Gene ID" )
      db = Association.add_reciprocal db if @organism_kb[organism].registry[database][1][:undirected]

      tsv, total_keys, source_field, target_field = [db, db.keys, db.key_field, db.fields.first]

      if target_field == "Ensembl Gene ID"
        pathway_field, gene_field = source_field, target_field
        total_genes = Gene.setup(tsv.values.flatten.compact.uniq, "Ensembl Gene ID", organism)
      else
        pathway_field, gene_field = target_field, source_field
        total_genes = total_keys
      end

      tsv.namespace = organism

      [tsv, total_genes, gene_field, pathway_field]
    end
  end

  #{{{ BASE AND GENE COUNTS

  input :masked_genes, :array, "Ensembl Gene ID list of genes to mask", []
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :pathway_base_counts => :tsv do |masked_genes, organism|
    database = clean_name
    log :loading_genes, "Loading genes from #{ database } #{ organism }"

    tsv, total_genes, gene_field, pathway_field = database_info database, organism

    if pathway_field != gene_field
      tsv = tsv.reorder pathway_field, [gene_field]
    else
      tsv = tsv.reorder 0, [:key]
    end

    counts = TSV.setup({}, :key_field => tsv.key_field, :fields => ["Bases"], :type => :single, :cast => :to_i, :namespace => organism)

    log :processing_database, "Processing database #{database}"
    tsv.with_monitor :desc => "Computing exon bases for pathways" do
      tsv.through do |pathway, values|
        next if values.empty?
        genes = values.first.flatten
        size = Gene.gene_list_exon_bases((genes.compact.uniq - masked_genes))
        counts[pathway] = size
      end
    end

    log :computing_exome_size, "Computing number of exome bases covered by pathway annotations"
    total_size = Gene.gene_list_exon_bases(total_genes - masked_genes)

    set_info :total_size, total_size
    set_info :total_gene_list, (total_genes - masked_genes)

    counts
  end

  input :masked_genes, :array, "Ensembl Gene ID list of genes to mask", []
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :pathway_gene_counts => :tsv do |masked_genes,organism|
    database = clean_name

    tsv, total_genes, gene_field, pathway_field = database_info database, organism

    tsv = tsv.reorder pathway_field, [gene_field] unless tsv.key_field == pathway_field and tsv.fields.first == gene_field

    counts = TSV.setup({}, :key_field => tsv.key_field, :fields => ["Genes"], :type => :single, :cast => :to_i, :namespace => organism)

    tsv.through do |pathway, values|
      genes = values[0]
      next if genes.nil? or genes.empty? 
      genes = genes.remove(masked_genes)
      num = genes.length
      counts[pathway] = num
    end

    set_info :total_genes, total_genes.remove(masked_genes).length
    set_info :total_gene_list, total_genes.remove(masked_genes).clean_annotations

    counts
  end

  #{{{ Mutation enrichment
   
  dep do |jobname, inputs| job(inputs[:baseline] || :pathway_base_counts, inputs[:database].to_s, inputs) end
  input :database, :select, "Database code", nil, :select_options => DATABASES
  input :baseline, :select, "Type of baseline to use", :pathway_base_counts, :select_options => [:pathway_base_counts, :pathway_gene_counts]
  input :mutations, :array, "Genomic Mutation"
  input :fdr, :boolean, "BH FDR corrections", true
  input :masked_genes, :array, "Ensembl Gene ID list of genes to mask", []
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :watson, :boolean, "Alleles reported in the watson (forward) strand", true
  task :mutation_pathway_enrichment => :tsv do |database,baseline,mutations,fdr,masked_genes,organism, watson|
    counts        = step(baseline).load
    total_covered = step(baseline).info[:total_size] || step(baseline).info[:total_genes]
    GenomicMutation.setup(mutations, "MutationEnrichment", organism, watson)


    affected_genes = mutations.genes.compact.flatten.uniq
    affected_genes = affected_genes.remove(masked_genes)

    # Get database tsv and native ids

    database_tsv, all_db_genes, db_gene_field, db_pathway_field = database_info database, organism

    all_db_genes = all_db_genes.ensembl.remove(masked_genes).compact.sort

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

  #{{{ Sample enrichment
  
  dep do |jobname, inputs| inputs[:baseline] ||= :pathway_base_counts; job(inputs[:baseline], inputs[:database].to_s, inputs) end
  input :database, :select, "Database code", nil, :select_options => DATABASES
  input :baseline, :select, "Type of baseline to use", :pathway_base_counts, :select_options => [:pathway_base_counts, :pathway_gene_counts]
  input :mutations, :tsv, "Genomic Mutation and Sample. Example row '10:12345678:A{TAB}Sample01{TAB}Sample02'"
  input :permutations, :integer, "Number of permutations in test", 10000
  input :fdr, :boolean, "BH FDR corrections", true
  input :masked_genes, :array, "Ensembl Gene ID list of genes to mask", []
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :watson, :boolean, "Alleles reported in the watson (forward) strand", true
  task :sample_pathway_enrichment => :tsv do |database,baseline,mutations,permutations,fdr,masked_genes,organism,watson|
    pathway_counts                         = step(baseline).load
    total_covered                          = step(baseline).info[:total_size] || step(baseline).info[:total_genes]
    total_pathway_genes_list               = step(baseline).info[:total_gene_list]

    mutations.extend TSV unless TSV === mutations

    if mutations.fields.nil?
      mutations.key_field = "Genomic Mutation"
      mutations.fields = ["Sample"]
      mutations.type = :double
    end

    database, all_db_genes, gene_field, pathway_field = database_info database, organism

    all_db_genes = all_db_genes.ensembl

    if database.key_field == pathway_field
      database_p2g = database
      database_g2p = database_p2g.reorder 0, [:key]
    else
      database_g2p = database
      database_p2g = database_p2g.reorder 0, [:key]
    end

    all_mutations = GenomicMutation.setup(mutations.keys, "MutationEnrichment", organism, watson)
    mutation_genes = Misc.process_to_hash(all_mutations){|all_mutations| all_mutations.genes}
    mutation_genes = Sequence.job(:affected_genes, clean_name, :mutations => all_mutations).run

    affected_samples_per_pathway = TSV.setup({}, :key_field => pathway_field, :fields => ["Sample"], :type => :flat)
    covered_genes_per_samples = {}
    all_samples = []
    sample_mutation_tokens = []
    covered_mutations = []
    log :classify, "Classifying mutations by pathway"
    mutations.slice("Sample").each do |mutation,samples|
      samples = [samples] unless Array === samples
      samples = samples.flatten

      next if mutation_genes[mutation].nil? or mutation_genes[mutation].empty?
      _genes = mutation_genes[mutation]
      _genes = _genes.to(gene_field)
      pathways = database_g2p.values_at(*_genes).compact.flatten.compact
      next if pathways.empty?
      pathways.each do |pathway|
        affected_samples_per_pathway[pathway] ||= []
        affected_samples_per_pathway[pathway].concat samples
      end
      samples.each do |sample|
        covered_genes_per_samples[sample] ||= []
        covered_genes_per_samples[sample].concat mutation_genes[mutation] unless mutation_genes[mutation].nil?
      end
      all_samples.concat samples

      if (mutation_genes[mutation] & all_db_genes).any?
        sample_mutation_tokens.concat samples
        covered_mutations << mutation
      end
    end

    affected_genes = mutation_genes.values.compact.flatten.uniq

    set_info :covered_mutations, covered_mutations.length

    pathways = pathway_counts.keys

    pathway_expected_counts = {}
    log :expected_counts, "Calculating expected counts"
    pathway_counts.with_unnamed do
      pathway_counts.with_monitor :desc => "Calculating expected counts" do
        affected_samples_per_pathway.with_unnamed do
          pathway_counts.through do |pathway, count|
            next unless affected_samples_per_pathway.include?(pathway) and affected_samples_per_pathway[pathway].any?
            ratio = count.to_f / total_covered
            num_token_list = R.eval_a "rbinom(#{ permutations }, #{ sample_mutation_tokens.length }, #{ ratio })"
            pathway_expected_counts[pathway] = num_token_list.collect{|num_tokens|
              # Add 1 to estabilize estimates
              Misc.sample(sample_mutation_tokens, num_tokens.to_i).uniq.length + 1
            }
          end
        end
      end
    end

    tsv = TSV.setup({}, :key_field => affected_samples_per_pathway.key_field, :fields => ["Sample", "Matches", "Expected", "Ratio", "Pathway total", "p-value", "Ensembl Gene ID"], :namespace => organism, :type => :double)
    log :pvalues, "Comparing observed vs expected counts"

    affected_samples_per_pathway.through do |pathway, samples|
      next unless samples.any?
      next unless pathway_expected_counts.include? pathway
      pathway_genes = database_p2g[pathway][0]

      samples = samples.uniq.select{|sample| (covered_genes_per_samples[sample] & pathway_genes).any?}
      # Add 1 to estabilize estimates
      count = samples.length
      expected = Misc.mean(pathway_expected_counts[pathway]).floor
      pvalue = pathway_expected_counts[pathway].select{|exp_c| exp_c > count}.length.to_f / permutations
      tsv[pathway] = [samples.sort, [count], [expected], [count.to_f / expected], [pathway_counts[pathway]], [pvalue], (pathway_genes & affected_genes)]
    end

    FDR.adjust_hash! tsv, 5 if fdr

    set_info :total_covered, total_covered

    tsv
  end
  export_asynchronous :sample_pathway_enrichment


  dep do |jobname, inputs| job(inputs[:baseline] || :pathway_base_counts, inputs[:database].to_s, inputs) end
  input :database, :select, "Database code", nil, :select_options => DATABASES
  input :baseline, :select, "Type of baseline to use", :pathway_base_counts, :select_options => [:pathway_base_counts, :pathway_gene_counts]
  input :gene_counts, :tsv, "Number of mutations per gene"
  input :mask_diseases, :boolean, "Mask disease related terms", true
  input :fdr, :boolean, "BH FDR corrections", true
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :gene_count_enrichment => :tsv do |database,baseline,gene_counts,mask_diseases,fdr,organism|
    counts        = step(baseline).load
    total_covered = step(baseline).info[:total_size] || step(baseline).info[:total_genes]

    gene_counts.identifiers = Organism.identifiers(Organism.default_code("Hsa"))
    gene_counts = gene_counts.change_key "Ensembl Gene ID"
    gene_counts = gene_counts.to_single

    affected_genes = gene_counts.keys
    total_mutations = gene_counts.values.flatten.compact.inject(0){|acc,e| acc += e.to_i}

    # Get database tsv and native ids

    database_tsv, all_db_genes, db_gene_field, db_pathway_field = database_info database, organism

    unless database_tsv.key_field == db_gene_field and database_tsv.fields.first == db_pathway_field
      log :reorder, "Reordering database"
      database_tsv = database_tsv.reorder db_gene_field, [db_pathway_field] 
    end

    # Annotate each pathway with the affected genes that are involved in it

    log :pathway_matches, "Finding affected genes per pathway"
    affected_genes_per_pathway = {}
    database_tsv.with_unnamed do
      affected_genes.each do |gene|
        next unless database_tsv[gene]
        pathways = database_tsv[gene].first.flatten
        next if pathways.nil? or pathways.empty?
        pathways.uniq.each do |pathway|
          affected_genes_per_pathway[pathway] ||= []
          affected_genes_per_pathway[pathway] << gene
        end
      end
    end

    if mask_diseases and not Gene == Entity.formats[db_pathway_field]
      Log.debug("Masking #{MASKED_TERMS * ", "}")
      masked = MASKED_IDS[database] ||= database_tsv.with_unnamed do
        terms = database_tsv.values.flatten.uniq
        terms = Misc.prepare_entity(terms, db_pathway_field)
        if terms.respond_to? :name
          terms.select{|t| t.name =~ /#{MASKED_TERMS * "|"}/i}
        else
          masked = nil
        end
      end
    else
      masked = nil
    end

    Gene.setup(affected_genes, "Ensembl Gene ID", Organism.default_code("Hsa"))

    log :covered_mutations, "Finding mutations overlapping genes in a pathway"
    all_db_genes = all_db_genes.ensembl

    log :pvalue, "Calculating binomial pvalues"
    pvalues = TSV.setup({}, :key_field => database_tsv.fields.first, :fields => ["Matches", "Pathway total", "Ensembl Gene ID"], :namespace => organism, :type => :double)
    counts.unnamed = true
    affected_genes_per_pathway.each do |pathway, genes|
      next if genes.empty?
      next if pathway.nil? or pathway.empty?
      next if masked and masked.include? pathway
      next unless counts.include? pathway
      pathway_total = counts[pathway]
      matches = gene_counts.values_at(*genes).inject(0){|acc,e| acc += e.nil? ? 0 : e.to_i }
      next if matches == 0

      common_genes = affected_genes.subset(genes).uniq
      next if common_genes.length < 3
      pvalues[pathway] = [[matches], [pathway_total], common_genes.sort_by{|g| g.name || g}]
    end
    
    pvalues = pvalues.R(<<-EOF, :key => pvalues.key_field)
      pvalues = apply(data, 1, function(v){ binom.test(as.numeric(v[1]), #{ total_mutations }, as.numeric(v[2]) /  #{total_covered.to_f}, 'greater')$p.value });
      data = cbind(data, p.value = pvalues);
      data = data[names(data)[c(1,2,4,3)]];
    EOF

    pvalues.process "p.value" do |v|
      v.to_f
    end

    pvalues = FDR.adjust_hash! pvalues, 2 if fdr

    pvalues.type = :double

    set_info :total_covered, total_covered

    pvalues
  end
  export_asynchronous :gene_count_enrichment

end
