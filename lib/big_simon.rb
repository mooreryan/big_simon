require "rya"
require "set"

require "big_simon/version"

Time.extend Rya::CoreExtensions::Time
Process.extend Rya::CoreExtensions::Process


module BigSimon
  class Utils
    def self.strip_suffix fname
      fname.sub /.fasta$|.fa$/, ""
    end

    # @param [Array<Hash>] host_data host info hash tables.  See functions in Parsers class.
    # @param [Array<String>] programs names of programs generating hash tables (in same order as host_data)
    def self.collate_host_results host_data, programs
      Rya::AbortIf.assert host_data.count == programs.count

      hosts = {}
      all_viruses = host_data.reduce(Set.new) { |acc, ht| acc + ht.keys }

      all_viruses.each do |virus|
        hosts[virus] = {}
      end

      host_data.each_with_index do |ht, idx|
        program = programs[idx]

        ht.each do |virus, host_scores|
          hosts[virus][program] = host_scores
        end
      end

      hosts
    end
  end

  # Project directories
  ROOT       = File.join __dir__, ".."
  BIN        = File.join ROOT, "vendor", "bin", "mac"
  SPEC       = File.join ROOT, "spec"
  TEST_FILES = File.join SPEC, "test_files"



  # Methods for parsing output files
  class Parsers

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
              with_index { |dist, idx| [hosts[idx], dist.to_f] }.
              sort_by { |_, dist| dist }

          host_info[virus] = dists
        end
      end

      host_info
    end

    # @note The viruses and hosts will have the ID rather than the file name.
    def self.wish fname
      viruses = nil

      host_info = {}
      File.open(fname, "rt").each_line.with_index do |line, idx|
        line.chomp!

        if idx.zero?
          viruses = line.split "\t"

          # Set up the hash table
          viruses.each do |virus|
            host_info[virus] = []
          end
        else
          ary  = line.split "\t"
          host = ary.shift

          ary.each_with_index do |val, idx|
            host_info[viruses[idx]] << [host, val.to_f]
          end
        end
      end

      host_info.each do |virus, hosts|
        hosts.sort_by! { |host, val| val }
        hosts.reverse!
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
      "-r #{outdir} -b"

      Process.run_and_time_it! "Building model", build_model
      Process.run_and_time_it! "Predicting host", predict

      FileUtils.rm_r model_dir if Dir.exist? model_dir
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
    end
  end
end
