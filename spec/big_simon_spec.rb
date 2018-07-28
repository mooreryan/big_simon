RSpec.describe BigSimon do
  it "has a version number" do
    expect(BigSimon::VERSION).not_to be nil
  end

  describe BigSimon::Parsers do
    describe "::vir_host_matcher" do
      let(:fname) { File.join BigSimon::TEST_FILES, "d2star_k6.csv" }

      it "parses vir_host_matcher output" do
        result   = BigSimon::Parsers.vir_host_matcher fname
        expected = {
            "AJ609634.fasta" => {
                best: "NC_002971.fasta",
                all:  [["NC_002971.fasta", 0.365838], ["NC_002163.fasta", 0.380042]]
            },
            "DQ113772.fasta" => {
                best: "NC_002971.fasta",
                all:  [["NC_002971.fasta", 0.360301], ["NC_002163.fasta", 0.404101]]
            },
        }

        expect(result).to eq expected
      end

      context "bad input data" do
        it "handles duplicate virus names"
        it "handles duplicate host names"
      end
    end
  end

  describe BigSimon::Runners do
    let(:exe) { File.join BigSimon::BIN, program_name }
    let(:outdir) { File.join BigSimon::TEST_FILES, "test_output", program_name }
    let(:vir_dir) { File.join BigSimon::TEST_FILES, "virus" }
    let(:host_dir) { File.join BigSimon::TEST_FILES, "host" }
    let(:threads) { 3 }

    describe "::wish" do
      let(:program_name) { "WIsH" }

      it "removes model dir" do
        FileUtils.rm_r outdir if Dir.exist? outdir

        expect {
          BigSimon::Runners.wish exe, vir_dir, host_dir, outdir, threads
        }.not_to raise_error

        expect(File).not_to exist(File.join outdir, "model")
      end
    end

    describe "::vir_host_matcher" do
      let(:program_name) { "vhm.py" }

      it "deletes the extra stuff" do
        expect {
          BigSimon::Runners.vir_host_matcher exe, vir_dir, host_dir, outdir
        }.not_to raise_error

        expect(File).not_to exist(File.join outdir, "tmp")
        expect(File).not_to exist(File.join outdir, "d2star_k6_main.html")
        expect(File).not_to exist(File.join outdir, "hostTaxa.txt_new.txt")
      end
    end
  end
end
