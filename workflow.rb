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
end

require 'mutation_enrichment/tasks/databases'
require 'mutation_enrichment/tasks/pathway_counts'
require 'mutation_enrichment/tasks/mutation_pathway_enrichment'
require 'mutation_enrichment/tasks/sample_pathway_enrichment'
  
module MutationEnrichment

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

