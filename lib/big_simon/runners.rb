require "tempfile"
require "parse_fasta"
require "parallel"

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

          # -k 3 index every third position in reference (broken now, bug in mummer)
          # -n -k 3 -threads 3
          # -n match only A C T G
          cmd = "mummer -n -threads #{threads} -qthreads #{threads} -maxmatch -l 15 #{host_f.path} #{vir_f.path} > #{mummer_outfname}"
          Process.run_and_time_it! "MUMMER", cmd
        end
      end

      virus             = nil
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

            host  = ary[0].sub(/___reverse$/, "").strip
            score = ary[3].to_i

            Rya::AbortIf.assert hit_table[virus].has_key?(host)

            # unless hit_table[virus].has_key? host
            #   hit_table[virus][host] = -1
            # end

            # We only want the longest hit.
            hit_table[virus][host] = score if score > hit_table[virus][host]

            # Track the overall max for scaling.
            overall_max_score      = score if score > overall_max_score
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

    # For scoring homology-ness, I just sum the bitscore for all significant hits for all genomes.
    #
    # @note I will make the specified outdir if it doesn't exist.
    # @note Assumes that the files end with *.fa
    # @note Assumes that the file names match the IDs.  This SHOULD be taken care of by the big_simon program.
    # @todo assert that fname thing matches sequence ID name.
    def self.homology vir_dir, host_dir, outdir, threads
      FileUtils.mkdir_p outdir

      host_orfs          = File.join outdir, "host_orfs.homology"
      host_orfs_blast_db = host_orfs + ".blast_db.homology"

      # Call ORFs on Hosts
      cmd = "cat #{host_dir}/*.fa | #{BigSimon::PRODIGAL} " \
      "-d #{host_orfs} " \
      "> /dev/null"

      Process.run_and_time_it! "Predicting host ORFs", cmd

      # Make blast db's for the host genes.
      cmd = "#{BigSimon::MAKEBLASTDB} " \
      "-in #{host_orfs} " \
      "-out #{host_orfs_blast_db} " \
      "-dbtype nucl"

      Process.run_and_time_it! "Making host blast db", cmd

      vir_genome_fnames = Dir.glob(vir_dir + "/*.fa")

      blast_info = Parallel.map(vir_genome_fnames, in_processes: threads) do |vir_genome_fname|
        vir_orfs      = File.join outdir, File.basename(vir_genome_fname) + ".vir_orfs.homology"
        blast_results = File.join outdir, File.basename(vir_genome_fname) + ".blast_results.homology"

        # this will be used as a viral ID.
        vir_simple_fname              = File.basename vir_genome_fname, ".fa"
        blast_table                   = {}
        blast_table[vir_simple_fname] = Hash.new 0

        # Call ORFs on the virus.
        cmd = "#{BigSimon::PRODIGAL} " \
        "-d #{vir_orfs} -p meta -i #{vir_genome_fname} " \
        "> /dev/null"

        Process.run_and_time_it! "Predicting ORFs for #{File.basename vir_genome_fname}", cmd

        # Blast the ORFs against genomes.
        cmd = "#{BigSimon::BLASTN} -query #{vir_orfs} -db #{host_orfs_blast_db} -outfmt 6 -evalue 0.01 -word_size 11 -out #{blast_results}"
        Process.run_and_time_it! "Blasting ORFs for #{File.basename vir_genome_fname}", cmd

        # Remove ORFs file
        FileUtils.rm vir_orfs if File.exist? vir_orfs

        Rya::AbortIf.logger.info { "Parsing #{blast_results}" }
        # Parse the blast.

        File.open(blast_results, "rt").each_line do |line|
          ary = line.chomp.split "\t"

          # The .sub() is to remove the annotation that prodigal gives.
          vir_id  = ary[0].sub(/_[0-9]+$/, "")
          host_id = ary[1].sub(/_[0-9]+$/, "")
          score   = ary[11].to_f

          Rya::AbortIf.assert blast_table.has_key?(vir_id), "blast_table: got #{vir_id} should have been #{vir_simple_fname}"

          blast_table[vir_id][host_id] += score
        end

        # Remove blast file
        # FileUtils.rm_r blast_results if File.exist? blast_results

        # Again, we're assuming the input is .fa, which the big_simon program SHOULD ensure.  TODO check these things with assertions.
        simple_vir_name = File.basename vir_genome_fname.sub(/.fa$/, "")

        [simple_vir_name, blast_table]
      end

      collated_blast_table = {}
      host_simple_names    = Dir.glob(host_dir + "/*.fa").map { |fname| File.basename(fname, ".fa") }

      Rya::AbortIf.assert host_simple_names.length == host_simple_names.uniq.length, "host simple names are not unique"

      Rya::AbortIf.logger.info { "Collating blast results" }

      # Get max score
      max_score = -1
      blast_info.each do |_, blast_table|
        blast_table.each do |vir_id, host_scores|
          this_max  = host_scores.values.max || -1 # sometimes there are no hits at all

          max_score = this_max if this_max > max_score
        end
      end
      Rya::AbortIf.assert max_score > -1, "didn't get any scores"


      klass = Class.new.extend Rya::CoreExtensions::Math
      blast_info.each do |simple_vir_name, blast_table|
        blast_table.each do |vir_id, host_scores|
          collated_blast_table[vir_id] = []

          host_simple_names.each do |host_id|
            scaled_score = klass.scale host_scores[host_id], 0, max_score, 1, 0

            host_table = { host: host_id, score: host_scores[host_id], scaled_score: scaled_score }
            collated_blast_table[vir_id] << host_table
          end
        end
      end

      pp collated_blast_table

      collated_blast_table
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

        [File.absolute_path(in_fname), File.absolute_path(out_fname)]
      end


      rcode_str = BigSimon::Utils.rcode fnames

      Object::File.open(File.join(outdir, "RCODE.r"), "w") do |f|
        f.puts rcode_str
        f.fsync # ensure no data is buffered


        cmd = "#{exe} #{f.path}"
        Process.run_and_time_it! "Drawing heatmaps", cmd
      end

      out_fnames = fnames.map(&:last)
    end
  end
end


