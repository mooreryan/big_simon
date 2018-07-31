module BigSimon
  # @todo These don't have unit tests yet.
  # @note Skips any duplicate IDs.  Only keeps the first one.
  class Utils
    def self.check_file! fname
      Rya::AbortIf.abort_if fname && !File.exist?(fname),
                            "#{fname} doesn't exist!  Try big_simon --help for help."
    end

    def self.check_opt! opts, arg
      Rya::AbortIf.abort_unless opts.send(:fetch, "#{arg}_given".to_sym),
                                "You must specify --#{arg.to_s.tr('_', '-')}.  Try big_simon --help for help."
    end

    def self.rcode fnames
      functions = %Q|
  library(reshape2)
library(gplots)
library(RColorBrewer)

file.join <- function(...) {
    paste(..., sep="/")
}

draw.heatmap <- function(infname, outfname) {
    dat <- read.table(infname, header=T, sep="\t")

    wide.dat <- dcast(dat, host ~ virus, value.var="score")

    hosts <- wide.dat[, 1]
    scores <- wide.dat[, 2:ncol(wide.dat)]
    scores.numeric <- apply(scores, 2, as.numeric)

    scores.matrix <- as.matrix(scores.numeric)

    rownames(scores.matrix) <- hosts

    palette <- "YlOrBr"
    col <- colorRampPalette(brewer.pal(n=9, palette))(n = 25)
    size <- 0.75

    pdf(outfname, height=5, width=8)

    heatmap.2(scores.matrix,
              trace="none", ## Disable those wonky lines.
              col=col, ## Set the color.

              ## Size opts
              margins=c(11, 11), cexRow=size, cexCol=size,

              ## Key labeling
              key.xlab="Score")

    invisible(dev.off())
}

|

      drawing = fnames.map do |in_fname, out_fname|
        %Q{

draw.heatmap("#{in_fname}", "#{out_fname}")
}
      end.join

      [functions, drawing].join "\n"
    end

    def self.scale_log_likelihood ll
      1 - Math.exp(ll)
    end

    # @note I also rename all the sequences in the tmp fasta files with the new ID.
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
              f.puts ">#{new_id}\n#{rec.seq}" # TODO HERE
            end
          end
        end
      end

      [name_map, all_ids]
    end

    def self.strip_suffix fname
      fname.sub /.fasta$|.fa$/, ""
    end
  end
end