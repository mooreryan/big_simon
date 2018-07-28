require "rya"

require "big_simon/version"

Time.extend Rya::CoreExtensions::Time
Process.extend Rya::CoreExtensions::Process

module BigSimon

  # Project directories
  ROOT       = File.join __dir__, ".."
  BIN        = File.join ROOT, "vendor", "bin"
  SPEC       = File.join ROOT, "spec"
  TEST_FILES = File.join SPEC, "test_files"

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
  end
end
