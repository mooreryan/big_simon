module BigSimon
  class Pipeline
    # @param collated_results_table { virus => host => score_type => program => score }
    def self.map_taxa collated_results_table, virus_name_map, host_name_map
      new_results_table = {}

      collated_results_table.each do |virus_name, host_table|
        if virus_name_map.include? virus_name
          new_virus_name = virus_name_map[virus_name]
        else
          new_virus_name = virus_name
        end

        new_results_table[new_virus_name] = {}

        host_table.each do |host_name, score_table|
          if host_name_map.include? host_name
            new_host_name = host_name_map[host_name]
          else
            new_host_name = host_name
          end

          new_results_table[new_virus_name][new_host_name] = score_table
        end
      end

      new_results_table
    end

    # @param [Array<Hash>] results_table host info hash tables.  See functions in Parsers class.
    # @param [Array<String>] programs names of programs generating hash tables (in same order as host_data)
    def self.collate_host_results results_table, programs
      Rya::AbortIf.assert results_table.count == programs.count

      virus_host_scores = {}
      all_viruses       = results_table.reduce(Set.new) { |acc, ht| acc + ht.keys }

      all_viruses.each do |virus|
        virus_host_scores[virus] = {}
      end

      results_table.each_with_index do |ht, idx|
        program = programs[idx]

        ht.each do |virus, host_scores|
          host_scores.each do |ht|
            host = ht[:host]
            score = ht[:score]
            scaled_score = ht[:scaled_score]

            unless virus_host_scores[virus].has_key? host
              virus_host_scores[virus][host] = { scores: {}, scaled_scores: {}}
            end

            virus_host_scores[virus][host][:scores][program] = score
            virus_host_scores[virus][host][:scaled_scores][program] = scaled_score
          end
        end
      end

      virus_host_scores
    end
  end
end