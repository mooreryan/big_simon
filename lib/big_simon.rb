require "rya"

require "big_simon/version"

Time.extend Rya::CoreExtensions::Time
Process.extend Rya::CoreExtensions::Process

module BigSimon

  # Project directories
  ROOT       = File.join __dir__, ".."
  BIN        = File.join ROOT, "vendor", "bin", "mac"
  SPEC       = File.join ROOT, "spec"
  TEST_FILES = File.join SPEC, "test_files"

  class Parsers

    def self.vir_host_matcher fname
      hosts = nil

      host_info = {}
      File.open(fname, "rt").each_line.with_index do |line, idx|
        line.chomp!
        line.sub! /,$/, "" # git rid of trailing commas

        if idx.zero?
          stat, *hosts = line.split ","
        else
          ary   = line.split ","
          virus = ary.shift

          dists = ary.map.
              with_index { |dist, idx| [hosts[idx], dist.to_f] }.
              sort_by { |_, dist| dist }

          best_host = dists[0][0]

          host_info[virus] = {
              best: best_host,
              all:  dists
          }
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
