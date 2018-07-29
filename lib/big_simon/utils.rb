module BigSimon
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

end