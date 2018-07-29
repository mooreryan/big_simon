module BigSimon
  class Runners

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

      outfiles = []

      Dir.glob("#{indir}/scores*.txt").each do |fname|
        extname = File.extname fname
        basename = File.basename fname, extname

        out_fname = File.join outdir, "#{basename}.heatmap.pdf"
        outfiles << out_fname

        rcode_str = rcode fname, out_fname

        Tempfile.open do |f|
          f.puts rcode_str
          f.fsync # ensure no data is buffered


          cmd = "#{exe} #{f.path}"
          Process.run_and_time_it! "Drawing heatmaps", cmd
        end
      end

      outfiles
    end
  end
end

def rcode infname, out_fname
  %Q{
library(reshape2)
library(gplots)
library(RColorBrewer)

file.join <- function(...) {
    paste(..., sep="/")
}

dat <- read.table("#{infname}", header=T, sep="\t")

wide.dat <- dcast(dat, host ~ virus, value.var="score")

hosts <- wide.dat[, 1]
scores <- wide.dat[, 2:ncol(wide.dat)]
scores.numeric <- apply(scores, 2, as.numeric)

scores.matrix <- as.matrix(scores.numeric)

rownames(scores.matrix) <- hosts

palette <- "YlOrBr"
col <- colorRampPalette(brewer.pal(n=9, palette))(n = 25)
size <- 0.75

pdf("#{out_fname}", height=5, width=8)

heatmap.2(scores.matrix,
          trace="none", ## Disable those wonky lines.
          col=col, ## Set the color.

          ## Size opts
          margins=c(11, 11), cexRow=size, cexCol=size,

          ## Key labeling
          key.xlab="Mean score (lower is better)")

invisible(dev.off())
}
end

