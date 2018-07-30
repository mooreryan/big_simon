require "tempfile"
require "parallel"

module BigSimon
  class Runners

    # This one's a bit different as it parses as well and returns original names.
    # @todo Also do the reverse of each genome in case it's a contig.
    def self.mummer exe, vir_dir, host_dir, outdir, threads
      klass = Class.new.extend Rya::CoreExtensions::Math
      FileUtils.mkdir_p outdir

      # TODO put these all in one file then do it?

      results = {}

      # Takes names in files and puts them to the file names
      name_map = {}

      Dir.glob(vir_dir + "/*").each do |vir_fname|
        host_fnames = Dir.glob(host_dir + "/*")

        return_value = Parallel.map(host_fnames, in_processes: threads) do |host_fname|

          # end

          # Dir.glob(host_dir + "/*").each do |host_fname|

          vir_base  = File.basename vir_fname
          host_base = File.basename host_fname
          outfname  = File.join outdir, "#{vir_base}___#{host_base}.mummer"

          virus = nil

          # -l is min length of a match TODO pull this into a const
          # -F to force 4 columns
          cmd = "#{exe} -F " \
                "-maxmatch " \
                "-l 15 " \
                "#{host_fname} " \
                "#{vir_fname} " \
                "> #{outfname}"

          Process.run_and_time_it! "Calculating matches", cmd

          # Note there should only be one '>' per file here.
          host  = nil
          score = 0
          File.open(outfname, "rt").each_line.with_index do |line, idx|
            if idx.zero?
              this_virus = line.chomp.sub(/^>/, "").sub(/___reverse$/, "").strip

              Rya::AbortIf::abort_unless(this_virus == virus, "OOPS") if virus

              virus ||= this_virus
            else
              ary = line.chomp.strip.split(" ")
              Rya::AbortIf.abort_unless ary.count == 4, "Problem parsing #{outfname} (mummer output)"

              host  = ary[0].sub(/___reverse$/, "").strip
              len   = ary[3].to_i

              score = len if len > score
            end
          end

          # unless results.has_key? virus
          #   results[virus] = []
          # end
          #
          # results[virus] << { host: host, score: score, scaled_score: nil }

          FileUtils.rm outfname

          host_table = { host: host, score: score, scaled_score: nil }

          [virus, host_table]
        end

        virus = return_value.first.first
        host_tables = return_value.map(&:last)

        all_scores = host_tables.map { |table| table[:score] }
        min        = all_scores.min
        max        = all_scores.max
        from       = 1
        to         = 0

        host_tables.each do |host_table|
          host_table[:scaled_score] = klass.scale host_table[:score], min, max, from, to
        end

        results[virus] = host_tables
      end

      results
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

    def self.heatmaps exe, indir, outdir
      FileUtils.mkdir_p outdir

      fnames = Dir.glob("#{indir}/scores*.txt").map do |in_fname|
        extname  = File.extname in_fname
        basename = File.basename in_fname, extname

        out_fname = File.join outdir, "#{basename}.heatmap.pdf"

        [in_fname, out_fname]
      end


      rcode_str = BigSimon::Utils.rcode fnames

      Object::Tempfile.open do |f|
        f.puts rcode_str
        f.fsync # ensure no data is buffered


        cmd = "#{exe} #{f.path}"
        Process.run_and_time_it! "Drawing heatmaps", cmd
      end

      out_fnames = fnames.map(&:last)
    end
  end
end


