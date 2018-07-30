require "tempfile"
require "parse_fasta"

module BigSimon
  class Runners

    # @note To match the other things, you'd like them to be key'd on the file name.
    def self.mummer exe, vir_dir, host_dir, outdir, threads
      klass = Class.new.extend Rya::CoreExtensions::Math
      FileUtils.mkdir_p outdir

      mummer_outfname = File.join outdir, "mummer_out.txt"

      virus_fnames = Dir.glob(vir_dir + "/*")
      host_fnames  = Dir.glob(host_dir + "/*")

      hit_table = {}

      Tempfile.open do |vir_f|
        Tempfile.open do |host_f|
          virus_fnames.each do |fname|
            Rya::AbortIf.assert fname.match(/.fa$/), "bad fname: #{fname}"

            Object::ParseFasta::SeqFile.open(fname).each_record do |rec|
              # id needs to be the file name
              new_id = File.basename fname.sub(/.fa$/, "")

              hit_table[new_id] = {}

              vir_f.puts ">#{new_id}\n#{rec.seq}"

              vir_f.puts ">#{new_id}___reverse\n#{rec.seq.reverse}"
            end
          end

          host_fnames.each do |fname|
            Rya::AbortIf.assert fname.match(/.fa$/), "bad fname: #{fname}"

            Object::ParseFasta::SeqFile.open(fname).each_record do |rec|
              new_id = File.basename fname.sub(/.fa$/, "")

              # Add this host to each virus in the hit_table
              hit_table.each do |virus, host_table|
                host_table[new_id] = 0 # set it to defualt score of 0
              end

              host_f.puts ">#{new_id}\n#{rec.seq}"
              host_f.puts ">#{new_id}___reverse\n#{rec.seq.reverse}"
            end
          end

          vir_f.fsync
          host_f.fsync

          cmd = "mummer -threads #{threads} -qthreads #{threads} -maxmatch -l 15 #{host_f.path} #{vir_f.path} > #{mummer_outfname}"
          Process.run_and_time_it! "MUMMER", cmd
        end
      end

      virus = nil
      overall_max_score = 0
      File.open(mummer_outfname, "rt").each_line.with_index do |line, idx|
        line.chomp!

        unless line.empty?
          if line.start_with? ">"
            virus = line.chomp.sub(/^>/, "").sub(/___reverse$/, "").strip

            # It can be duplicated as there are forward and reverse for each sequence (in case they're contigs.)

            Rya::AbortIf.assert hit_table.has_key?(virus)
            # unless hit_table.has_key? virus
            #   hit_table[virus] = {}
            # end
          else
            ary = line.strip.split " "

            host = ary[0].sub(/___reverse$/, "").strip
            score = ary[3].to_i

            Rya::AbortIf.assert hit_table[virus].has_key?(host)

            # unless hit_table[virus].has_key? host
            #   hit_table[virus][host] = -1
            # end

            # We only want the longest hit.
            hit_table[virus][host] = score if score > hit_table[virus][host]

            # Track the overall max for scaling.
            overall_max_score = score if score > overall_max_score
          end
        end
      end

      results_table = {}

      min, max, from, to = 0, overall_max_score, 1, 0

      hit_table.each do |virus, host_table|
        results_table[virus] = []

        host_table.each do |host, score|
          scaled_score = klass.scale score, min, max, from, to

          results_table[virus] << { host: host, score: score, scaled_score: scaled_score }
        end
      end

      results_table
    end

    # This one's a bit different as it parses as well and returns original names.
    # @todo Also do the reverse of each genome in case it's a contig.
    def self.mummer2 exe, vir_dir, host_dir, outdir, threads
      klass = Class.new.extend Rya::CoreExtensions::Math
      FileUtils.mkdir_p outdir

      # TODO put these all in one file then do it?

      results = {}

      # Takes names in files and puts them to the file names
      name_map = {}

      Dir.glob(vir_dir + "/*").each do |vir_fname|
        this_virus_scores = []
        virus             = nil

        Dir.glob(host_dir + "/*").each do |host_fname|
          vir_base  = File.basename vir_fname
          host_base = File.basename host_fname
          outfname  = File.join outdir, "#{vir_base}___#{host_base}.mummer"

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

          this_virus_scores << score

          unless results.has_key? virus
            results[virus] = []
          end

          results[virus] << { host: host, score: score, scaled_score: nil }

          FileUtils.rm outfname
        end

        # This was the original scaling, i.e. per virus
        # min = 0 # this_virus_scores.min # Technically, this should range from 0 to 15.  Any data missing from this table would give a zero.  TODO we don't actually account for this though.
        # max = this_virus_scores.max
        # from = 1
        # to = 0
        #
        # results[virus].each do |host_table|
        #   host_table[:scaled_score] = klass.scale host_table[:score], min, max, from, to
        # end
      end

      all_scores = []
      results.each do |virus, host_tables|
        all_scores << host_tables.map { |table| table[:score] }
      end

      all_scores.flatten!
      max = all_scores.max

      results.each do |virus, host_tables|
        host_tables.each do |host_table|
          host_table[:scaled_score] = klass.scale host_table[:score], 0, max, 1, 0
        end
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


