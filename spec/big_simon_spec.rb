module SpecConst
  VIRUS_1 = "AJ609634.snazzy.lala.long_name"
  VIRUS_2 = "DQ113772"

  HOST_1 = "NC_002971"
  HOST_2 = "NC_002163.snazzy.lala.long_name"

  PROGRAM_1 = "VirHostMatcher"
  PROGRAM_2 = "WIsH"
end

RSpec.describe BigSimon do
  it "has a version number" do
    expect(BigSimon::VERSION).not_to be nil
  end

  let(:parsed_output_vhm) do
    # These values checked by hand
    a1, a2 = 0.365838, 0.380042
    b1, b2 = 0.360301, 0.404101

    {
      SpecConst::VIRUS_1 => [
        { host: SpecConst::HOST_1, score: a1, scaled_score: a1 },
        { host: SpecConst::HOST_2, score: a2, scaled_score: a2 }
      ],
      SpecConst::VIRUS_2 => [
        { host: SpecConst::HOST_1, score: b1, scaled_score: b1 },
        { host: SpecConst::HOST_2, score: b2, scaled_score: b2 },
      ],
    }
  end
  let(:parsed_output_wish) do
    # These numbers varified by hand
    a1, a2 = -1.38868, -1.39021
    b1, b2 = -1.37943, -1.38371

    a1_scaled, a2_scaled, b1_scaled, b2_scaled = [a1, a2, b1, b2].map { |val| BigSimon::Utils.scale_log_likelihood val }

    {
      SpecConst::VIRUS_1 => [
        { host: SpecConst::HOST_1, score: a1, scaled_score: a1_scaled },
        { host: SpecConst::HOST_2, score: a2, scaled_score: a2_scaled },
      ],
      SpecConst::VIRUS_2 => [
        { host: SpecConst::HOST_2, score: b1, scaled_score: b1_scaled },
        { host: SpecConst::HOST_1, score: b2, scaled_score: b2_scaled },
      ],
    }
  end


  describe BigSimon::Pipeline do
    describe "::collate_host_results" do
      it "raises if length of inputs doesn't match" do
        expect {
          BigSimon::Pipeline::collate_host_results [Hash.new], ["", ""]
        }.to raise_error Rya::AbortIf::Assert::AssertionFailureError
      end

      it "collates host results" do
        expected = {
          SpecConst::VIRUS_1 => {
            SpecConst::HOST_1 => {
              scores:        {
                SpecConst::PROGRAM_1 => 0.365838,
                SpecConst::PROGRAM_2 => -1.38868,
              },
              scaled_scores: {
                SpecConst::PROGRAM_1 => 0.365838,
                SpecConst::PROGRAM_2 => BigSimon::Utils.scale_log_likelihood(-1.38868),
              },
            },
            SpecConst::HOST_2 => {
              scores:        {
                SpecConst::PROGRAM_1 => 0.380042,
                SpecConst::PROGRAM_2 => -1.39021,
              },
              scaled_scores: {
                SpecConst::PROGRAM_1 => 0.380042,
                SpecConst::PROGRAM_2 => BigSimon::Utils.scale_log_likelihood(-1.39021),
              },
            },
          },
          SpecConst::VIRUS_2 => {
            SpecConst::HOST_1 => {
              scores:        {
                SpecConst::PROGRAM_1 => 0.360301,
                SpecConst::PROGRAM_2 => -1.38371,
              },
              scaled_scores: {
                SpecConst::PROGRAM_1 => 0.360301,
                SpecConst::PROGRAM_2 => BigSimon::Utils.scale_log_likelihood(-1.38371),
              },
            },
            SpecConst::HOST_2 => {
              scores:        {
                SpecConst::PROGRAM_1 => 0.404101,
                SpecConst::PROGRAM_2 => -1.37943,
              },
              scaled_scores: {
                SpecConst::PROGRAM_1 => 0.404101,
                SpecConst::PROGRAM_2 => BigSimon::Utils.scale_log_likelihood(-1.37943),
              },
            },
          },
        }

        actual = BigSimon::Pipeline.collate_host_results [parsed_output_vhm, parsed_output_wish], [SpecConst::PROGRAM_1, SpecConst::PROGRAM_2]

        expect(actual).to eq expected
      end
    end
  end

  describe BigSimon::Parsers do
    describe "::vir_host_matcher" do
      let(:fname) { File.join BigSimon::TEST_FILES, "d2star_k6.csv" }

      it "parses vir_host_matcher output" do
        result = BigSimon::Parsers.vir_host_matcher fname

        expect(result).to eq parsed_output_vhm
      end

      context "bad input data" do
        it "handles duplicate virus names"
        it "handles duplicate host names"
      end
    end

    describe "::wish" do
      let(:fname) { File.join BigSimon::TEST_FILES, "llikelihood.matrix" }

      it "parses wish output" do
        result = BigSimon::Parsers.wish fname

        expect(result).to eq parsed_output_wish
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
      let(:program_name) { SpecConst::PROGRAM_2 }

      it "has the right output files" do
        FileUtils.rm_r outdir if Dir.exist? outdir

        expected_outf = File.join outdir, "wish.txt"
        actual_outf   = nil

        expect {
          actual_outf = BigSimon::Runners.wish exe, vir_dir, host_dir, outdir, threads
        }.not_to raise_error

        expect(File).not_to exist(File.join outdir, "model")

        expect(File).to exist(expected_outf)
        expect(actual_outf).to eq expected_outf
      end
    end

    describe "::vir_host_matcher" do
      let(:program_name) { "vhm.py" }

      # Do all the output tests in one run so that we don't rerun the pipeline a lot.
      it "has the right output files" do
        FileUtils.rm_r outdir if Dir.exist? outdir

        expected_outf = File.join outdir, "vir_host_matcher.txt"
        actual_outf   = nil

        expect {
          actual_outf = BigSimon::Runners.vir_host_matcher exe, vir_dir, host_dir, outdir
        }.not_to raise_error

        expect(File).not_to exist(File.join outdir, "tmp")
        expect(File).not_to exist(File.join outdir, "d2star_k6_main.html")
        expect(File).not_to exist(File.join outdir, "hostTaxa.txt_new.txt")

        expect(File).to exist(expected_outf)
        expect(actual_outf).to eq expected_outf
      end
    end
  end
end
