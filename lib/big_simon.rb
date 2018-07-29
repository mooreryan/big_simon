require "rya"
require "set"
require "pp"

require "big_simon/version"

Time.extend Rya::CoreExtensions::Time
Process.extend Rya::CoreExtensions::Process
Array.include Rya::CoreExtensions::Array
Math.extend Rya::CoreExtensions::Math

# @todo Does this only work with negative numbers?
# https://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
def normalize_log_likelihood_scores scores
  max = scores.max.to_f
  tmp = scores.map { |score| Math.exp(score - max) }

  min = tmp.min
  tmp.map { |val| val / (1 + min) }
end

module BigSimon
  # Project directories
  ROOT       = File.join __dir__, ".."
  BIN        = File.join ROOT, "vendor", "bin", "mac"
  SPEC       = File.join ROOT, "spec"
  TEST_FILES = File.join SPEC, "test_files"
  WISH       = File.join BIN, "WIsH"
  VHM        = File.join BIN, "vhm.py"

  # @todo These don't have unit tests yet.
  # @note Skips any duplicate IDs.  Only keeps the first one.
  class Utils
    def self.scale_log_likelihood ll
      1 - Math.exp(ll)
    end

    def self.set_up_tmp_dirs fastas, tmpdir, which
      Object::FileUtils.mkdir_p tmpdir

      name_map = {}
      all_ids  = Set.new

      seq_num = -1
      fastas.each do |fname|
        ParseFasta::SeqFile.open(fname).each_record do |rec|
          if all_ids.include? rec.id
            Rya::AbortIf.logger.warn { "#{rec.id} was seen more than one time!  Duplicate organism IDs are not allowed, so we will only keep the first one." }
          else
            all_ids << rec.id

            seq_num += 1

            new_id           = "#{which}_#{seq_num}"
            name_map[new_id] = rec.id

            outfname = File.join tmpdir, "#{new_id}.fa"

            File.open(outfname, "w") do |f|
              f.puts rec
            end
          end
        end
      end

      [name_map, all_ids]
    end

    def self.strip_suffix fname
      fname.sub /.fasta$|.fa$/, ""
    end

    def self.check_file! fname
      Rya::AbortIf.abort_if fname && !File.exist?(fname),
                            "#{fname} doesn't exist!  Try big_simon --help for help."
    end

    def self.check_opt! opts, arg
      Rya::AbortIf.abort_unless opts.send(:fetch, "#{arg}_given".to_sym),
                                "You must specify --#{arg.to_s.tr('_', '-')}.  Try big_simon --help for help."
    end
  end

  class Pipeline
    # @param [Array<Hash>] host_data host info hash tables.  See functions in Parsers class.
    # @param [Array<String>] programs names of programs generating hash tables (in same order as host_data)
    def self.collate_host_results host_data, programs
      Rya::AbortIf.assert host_data.count == programs.count

      virus_host_scores = {}
      all_viruses       = host_data.reduce(Set.new) { |acc, ht| acc + ht.keys }

      all_viruses.each do |virus|
        virus_host_scores[virus] = {}
      end

      host_data.each_with_index do |ht, idx|
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
          scores_norm = normalize_log_likelihood_scores scores

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


  class Runners

    # Runs the WIsH program
    #
    # @raise [AbortIf::Exit] if commands fail
    def self.wish exe, vir_dir, host_dir, outdir, threads
      model_dir = File.join outdir, "model"

      FileUtils.mkdir_p model_dir

      build_model = "#{exe} " \
      "-t #{threads} " \
      "-c build " \
      "-g #{host_dir} " \
      "-m #{model_dir}"

      predict = "#{exe} " \
      "-t #{threads} " \
      "-c predict " \
      "-g #{vir_dir} " \
      "-m #{model_dir} " \
      "-r #{outdir}"

      Process.run_and_time_it! "Building model", build_model
      Process.run_and_time_it! "Predicting host", predict

      FileUtils.rm_r model_dir if Dir.exist? model_dir

      outf     = File.join outdir, "llikelihood.matrix"
      new_outf = File.join outdir, "wish.txt"
      FileUtils.mv outf, new_outf

      new_outf
    end

    def self.vir_host_matcher exe, vir_dir, host_dir, outdir
      FileUtils.mkdir_p outdir

      cmd = "python #{exe} " \
      "-v #{vir_dir} " \
      "-b #{host_dir} " \
      "-o #{outdir} " \
      "-d 1" # only compute d2star dissimilarity

      Process.run_and_time_it! "Computing d2star dissimilarity", cmd

      tmp_dir = File.join outdir, "tmp"
      FileUtils.rm_r tmp_dir if Dir.exist? tmp_dir

      bad_files = %w[d2star_k6_main.html hostTaxa.txt_new.txt]
      bad_files.each do |fname|
        path = File.join outdir, fname

        FileUtils.rm path if File.exist? path
      end

      outf     = File.join outdir, "d2star_k6.csv"
      new_outf = File.join outdir, "vir_host_matcher.txt"
      FileUtils.mv outf, new_outf

      new_outf
    end
  end
end
