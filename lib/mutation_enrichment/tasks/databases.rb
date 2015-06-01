module MutationEnrichment
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
end
