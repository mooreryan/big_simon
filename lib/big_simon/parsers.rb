module BigSimon
  # Methods for parsing output files
  class Parsers

    # @note VirHostMatcher returns true distances that run from 0 to 1, so it doesn't need scaling.
    # @note VirHostMatcher includes the whole file name as the id of the organism, so we chop off some common endings.
    def self.vir_host_matcher fname
      hosts = nil

      host_info = {}
      File.open(fname, "rt").each_line.with_index do |line, idx|
        line.chomp!
        line.sub! /,$/, "" # git rid of trailing commas

        if idx.zero?
          stat, *hosts = line.split ","

          hosts.map! { |str| BigSimon::Utils.strip_suffix str }
        else
          ary   = line.split ","
          virus = BigSimon::Utils.strip_suffix ary.shift

          # In this case the best value is the lowest distance.
          dists = ary.map.
            with_index do |dist, idx|
            { host: hosts[idx], score: dist.to_f, scaled_score: dist.to_f }
          end.sort_by { |ht| ht[:scaled_score] }


          host_info[virus] = dists
        end
      end

      host_info
    end

    # @note WIsH gives log likelihoods so the scaled value is actually scaled.
    # @note The viruses and hosts will have the ID rather than the file name.
    def self.wish fname
      viruses = nil

      host_info = {}

      hosts = nil
      File.open(fname, "rt").each_line.map.with_index do |line, idx|
        line.chomp!

        if idx.zero?
          ary = line.split("\t")
          ary.unshift("")
        else
          ary = line.split("\t")
        end
      end.transpose.each_with_index do |line_ary, idx|
        if idx.zero?
          hosts = line_ary.drop(1)
        else
          virus = line_ary.shift

          scores      = line_ary.map(&:to_f)

          host_vals = scores.map.with_index do |score, idx|
            { host: hosts[idx], score: score, scaled_score: 1 - Math.exp(score) }
          end

          host_info[virus] = host_vals
        end

        host_info.each do |virus, hosts|
          hosts.sort_by! { |ht| ht[:scaled_score] }
        end
      end

      host_info
    end
  end

end