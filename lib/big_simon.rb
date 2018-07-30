require "rya"
require "set"
require "pp"

require "big_simon/version"

require "big_simon/utils"

require "big_simon/runners"
require "big_simon/parsers"
require "big_simon/pipeline"

Time.extend Rya::CoreExtensions::Time
Process.extend Rya::CoreExtensions::Process
Array.include Rya::CoreExtensions::Array
Math.extend Rya::CoreExtensions::Math

module BigSimon
  ROOT       = File.join __dir__, ".."
  BIN        = File.join ROOT, "vendor", "bin", "mac"
  SPEC       = File.join ROOT, "spec"
  TEST_FILES = File.join SPEC, "test_files"

  # Programs
  WISH       = File.join BIN, "WIsH"
  VHM        = File.join BIN, "vhm.py"
  MUMMER     = File.join BIN, "mummer"
  RSCRIPT = "Rscript"

  BLASTN = File.join BIN, "blastn"
  MAKEBLASTDB = File.join BIN, "makeblastdb"
  PRODIGAL = File.join BIN, "prodigal"
end
