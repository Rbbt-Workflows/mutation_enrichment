module MutationEnrichment

  #dep do |jobname, inputs| inputs[:baseline] ||= :pathway_base_counts; job(inputs[:baseline], inputs[:database].to_s, inputs) end
  #input :database, :select, "Database code", nil, :select_options => DATABASES
  #input :baseline, :select, "Type of baseline to use", :pathway_base_counts, :select_options => [:pathway_base_counts, :pathway_gene_counts]
  #input :mutations, :tsv, "Genomic Mutation and Sample. Example row '10:12345678:A{TAB}Sample01{TAB}Sample02'"
  #input :permutations, :integer, "Number of permutations in test", 10000
  #input :fdr, :boolean, "BH FDR corrections", true
  #input :masked_genes, :array, "Ensembl Gene ID list of genes to mask", []
  #input :organism, :string, "Organism code", Organism.default_code("Hsa")
  #input :watson, :boolean, "Alleles reported in the watson (forward) strand", true
  #task :sample_pathway_enrichment_old => :tsv do |database,baseline,mutations,permutations,fdr,masked_genes,organism,watson|
  #  pathway_counts                         = step(baseline).load
  #  total_covered                          = step(baseline).info[:total_size] || step(baseline).info[:total_genes]
  #  total_pathway_genes_list               = step(baseline).info[:total_gene_list]

  #  mutations.extend TSV unless TSV === mutations

  #  if mutations.fields.nil?
  #    mutations.key_field = "Genomic Mutation"
  #    mutations.fields = ["Sample"]
  #    mutations.type = :double
  #  end

  #  database, all_db_genes, gene_field, pathway_field = database_info database, organism

  #  all_db_genes = all_db_genes.ensembl

  #  log :reordering, "Reordering database"
  #  if database.key_field == pathway_field
  #    database_p2g = database
  #    database_g2p = database_p2g.reorder 0, [:key]
  #  else
  #    database_g2p = database
  #    database_p2g = database_p2g.reorder 0, [:key]
  #  end

  #  log :affected_genes, "Getting affected genes"
  #  all_mutations = GenomicMutation.setup(mutations.keys, "MutationEnrichment", organism, watson)
  #  mutation_genes = Sequence.job(:affected_genes, clean_name, :mutations => all_mutations).run

  #  affected_samples_per_pathway = TSV.setup({}, :key_field => pathway_field, :fields => ["Sample"], :type => :flat)
  #  covered_genes_per_samples = {}
  #  all_samples = []
  #  sample_mutation_tokens = []
  #  covered_mutations = []
  #  log :classify, "Classifying mutations by pathway"
  #  TSV.traverse mutations, :fields => ["Sample"], :type => :flat do |mutation,samples|
  #    mutation = mutation.first if Array === mutation
  #    samples = samples.split("|") unless Array === samples
  #    samples = samples.flatten

  #    next if mutation_genes[mutation].nil? or mutation_genes[mutation].empty?
  #    _genes = mutation_genes[mutation]
  #    _genes = _genes.to(gene_field)
  #    pathways = database_g2p.values_at(*_genes).compact.flatten.compact
  #    next if pathways.empty?
  #    pathways.each do |pathway|
  #      affected_samples_per_pathway[pathway] ||= []
  #      affected_samples_per_pathway[pathway].concat samples
  #    end
  #    samples.each do |sample|
  #      covered_genes_per_samples[sample] ||= []
  #      covered_genes_per_samples[sample].concat mutation_genes[mutation] unless mutation_genes[mutation].nil?
  #    end
  #    all_samples.concat samples

  #    if (mutation_genes[mutation] & all_db_genes).any?
  #      sample_mutation_tokens.concat samples
  #      covered_mutations << mutation
  #    end
  #  end

  #  affected_genes = mutation_genes.values.compact.flatten.uniq

  #  set_info :covered_mutations, covered_mutations.length

  #  pathways = pathway_counts.keys

  #  pathway_expected_counts = {}
  #  log :expected_counts, "Calculating expected counts"
  #  pathway_counts.with_unnamed do
  #    pathway_counts.with_monitor :desc => "Calculating expected counts" do
  #      affected_samples_per_pathway.with_unnamed do
  #        pathway_counts.through do |pathway, count|
  #          next unless affected_samples_per_pathway.include?(pathway) and affected_samples_per_pathway[pathway].length > 1
  #          ratio = count.to_f / total_covered
  #          num_token_list = R.eval_a "rbinom(#{ permutations }, #{ sample_mutation_tokens.length }, #{ ratio })"
  #          pathway_expected_counts[pathway] = num_token_list.collect{|num_tokens|
  #            # Add 1 to estabilize estimates
  #            Misc.sample(sample_mutation_tokens, num_tokens.to_i).uniq.length + 1
  #          }
  #        end
  #      end
  #    end
  #  end

  #  tsv = TSV.setup({}, :key_field => affected_samples_per_pathway.key_field, :fields => ["Sample", "Matches", "Expected", "Ratio", "Pathway total", "p-value", "Ensembl Gene ID"], :namespace => organism, :type => :double)
  #  log :pvalues, "Comparing observed vs expected counts"

  #  affected_samples_per_pathway.through do |pathway, samples|
  #    next unless samples.any?
  #    next unless pathway_expected_counts.include? pathway
  #    pathway_genes = database_p2g[pathway][0]

  #    samples = samples.uniq.select{|sample| (covered_genes_per_samples[sample] & pathway_genes).any?}
  #    # Add 1 to estabilize estimates
  #    count = samples.length
  #    expected = Misc.mean(pathway_expected_counts[pathway]).floor
  #    pvalue = pathway_expected_counts[pathway].select{|exp_c| exp_c > count}.length.to_f / permutations
  #    tsv[pathway] = [samples.sort, [count], [expected], [count.to_f / expected], [pathway_counts[pathway]], [pvalue], (pathway_genes & affected_genes)]
  #  end

  #  FDR.adjust_hash! tsv, 5 if fdr

  #  set_info :total_covered, total_covered

  #  tsv
  #end









  #{{{




  dep do |jobname, inputs| inputs[:baseline] ||= :pathway_base_counts; job(inputs[:baseline], inputs[:database].to_s, inputs) end
  input :database, :select, "Database code", nil, :select_options => DATABASES
  input :baseline, :select, "Type of baseline to use", :pathway_base_counts, :select_options => [:pathway_base_counts, :pathway_gene_counts]
  input :mutations, :tsv, "Genomic Mutation and Sample. Example row '10:12345678:A{TAB}Sample01{TAB}Sample02'", nil, :stream => false
  input :permutations, :integer, "Number of permutations in test", 10000
  input :fdr, :boolean, "BH FDR corrections", true
  input :masked_genes, :array, "Ensembl Gene ID list of genes to mask", []
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :watson, :boolean, "Alleles reported in the watson (forward) strand", true
  dep Sequence, :affected_genes, :compute => :produce, :mutations => :mutations, :organism => :organism, :watson => :watson do |jobname,options|
    all_mutations = TSV.traverse options[:mutations], :into => :stream, :type => :keys do |mutation|
      mutation = mutation.first if Array === mutation
      mutation
    end
    Sequence.job(:affected_genes, jobname, options.merge(:mutations => all_mutations))
  end
  task :sample_pathway_enrichment => :tsv do |database,baseline,mutations,permutations,fdr,masked_genes,organism,watson|
    pathway_counts                         = step(baseline).load
    total_covered                          = step(baseline).info[:total_size] || step(baseline).info[:total_genes]
    total_pathway_genes_list               = step(baseline).info[:total_gene_list]

    database, all_db_genes, gene_field, pathway_field = database_info database, organism

    all_db_genes = all_db_genes.ensembl

    log :reordering, "Reordering database"
    if database.key_field == pathway_field
      database_p2g = database
      database_g2p = database_p2g.reorder 0, [:key], :persist => true
    else
      database_g2p = database
      database_p2g = database_p2g.reorder 0, [:key], :persist => true
    end

    affected_samples_per_pathway = TSV.setup({}, :key_field => pathway_field, :fields => ["Sample"], :type => :flat)
    covered_genes_per_samples = {}
    all_samples = []
    sample_mutation_tokens = []
    covered_mutations = []
    mutation_genes = nil

    TmpFile.with_file do |tmp_file|
      log :saving_orig_file
      Open.open(tmp_file, :mode => 'w') do |tmp_file_io|
        TSV.traverse mutations, :type => :array, :into => tmp_file_io do |*parts|
          parts.flatten * "\t"
        end
      end

      mutations = tmp_file

      #log :affected_genes, "Getting affected genes"
      #all_mutations = TSV.traverse mutations, :into => :stream, :type => :keys do |mutation|
      #  mutation = mutation.first if Array === mutation
      #  mutation
      #end

      #mutation_genes_job = Sequence.job(:affected_genes, clean_name, :mutations => all_mutations, :organism => organism)
      #mutation_genes_job.recursive_clean.produce
      mutation_genes_job = step(:affected_genes)
      mutation_genes = mutation_genes_job.path.tsv :persist => true, :persist_file => file('mutation_genes.tch')

      TSV.traverse mutations, :bar => self.progress_bar("Classifying mutations by pathway"), :type => :flat do |mutation,samples|
        mutation = mutation.first if Array === mutation
        samples = samples.split("|") unless Array === samples
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
    end

    affected_genes = mutation_genes.values.compact.flatten.uniq

    set_info :covered_mutations, covered_mutations.length

    pathways = pathway_counts.keys

    affected_pathways =  affected_samples_per_pathway.keys
    pathway_expected_counts = {}
    affected_samples_per_pathway.with_unnamed do
      TSV.traverse pathway_counts, :bar => self.progress_bar("Calculating expected counts") do |pathway, count|
        pathway = pathway.first if Array === pathway
        next unless affected_samples_per_pathway.include?(pathway) and affected_samples_per_pathway[pathway].length > 1
        ratio = count.to_f / total_covered
        num_token_list = R.eval_a "rbinom(#{ permutations }, #{ sample_mutation_tokens.length }, #{ ratio })"
        pathway_expected_counts[pathway] = num_token_list.collect{|num_tokens|
          # Add 1 to estabilize estimates
          Misc.sample(sample_mutation_tokens, num_tokens.to_i).uniq.length + 1
        }
      end
    end

    tsv = TSV.setup({}, :key_field => affected_samples_per_pathway.key_field, :fields => ["Sample", "Matches", "Expected", "Ratio", "Pathway total", "p-value", "Ensembl Gene ID"], :namespace => organism, :type => :double)
    log :pvalues, "Comparing observed vs expected counts"

    TSV.traverse affected_samples_per_pathway, :bar => self.progress_bar("Comparing observed vs expected counts") do |pathway, samples|
      pathway = pathway.first if Array === pathway
      next unless samples.any?
      next unless pathway_expected_counts.include? pathway
      pathway_genes = database_p2g[pathway][0]
      matched_genes = (pathway_genes & affected_genes)
      next unless matched_genes.length > 1

      samples = samples.uniq.select{|sample| (covered_genes_per_samples[sample] & pathway_genes).any?}
      # Add 1 to estabilize estimates
      count = samples.length
      expected = Misc.mean(pathway_expected_counts[pathway]).floor
      next if count <= expected

      better_count = pathway_expected_counts[pathway].select{|exp_c| exp_c >= count}.length
      pvalue = (1.0 + better_count) / (1.0 + permutations)

      tsv[pathway] = [samples.sort, [count], [expected], [count.to_f / expected], [pathway_counts[pathway]], [pvalue], matched_genes]
    end

    FDR.adjust_hash! tsv, 5 if fdr

    set_info :total_covered, total_covered

    tsv
  end
  export_asynchronous :sample_pathway_enrichment
end
