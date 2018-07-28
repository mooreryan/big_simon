RSpec.describe BigSimon do
  it "has a version number" do
    expect(BigSimon::VERSION).not_to be nil
  end

  describe BigSimon::Runners do
    let(:exe) { File.join BigSimon::BIN, program_name }
    let(:outdir) { File.join BigSimon::TEST_FILES, "test_output", program_name }
    let(:vir_dir) { File.join BigSimon::TEST_FILES, program_name, "virus" }
    let(:host_dir) { File.join BigSimon::TEST_FILES, program_name, "host" }
    let(:threads) { 3 }

    describe "#wish" do
      let(:program_name) { "WIsH" }

      it "removes model dir" do
        FileUtils.rm_r outdir if Dir.exist? outdir

        expect {
          BigSimon::Runners.wish exe, vir_dir, host_dir, outdir, threads
        }.not_to raise_error

        expect(File).not_to exist(File.join outdir, "model")
      end
    end
  end
end
