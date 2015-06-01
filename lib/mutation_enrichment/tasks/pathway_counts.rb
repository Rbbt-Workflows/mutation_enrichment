
module MutationEnrichment
  #{{{ BASE AND GENE COUNTS

  input :masked_genes, :array, "Ensembl Gene ID list of genes to mask", []
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :pathway_base_counts => :tsv do |masked_genes, organism|
    database = clean_name
    log :loading_genes, "Loading genes from #{ database } #{ organism }"

    tsv, total_genes, gene_field, pathway_field = database_info database, organism

    if pathway_field != gene_field
      tsv = tsv.reorder pathway_field, [gene_field], :persist => true
    else
      tsv = tsv.reorder 0, [:key], :persist => true
    end

    counts = TSV.setup({}, :key_field => tsv.key_field, :fields => ["Bases"], :type => :single, :cast => :to_i, :namespace => organism)

    log :processing_database, "Processing database #{database}"
    tsv.unnamed = true
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
end
