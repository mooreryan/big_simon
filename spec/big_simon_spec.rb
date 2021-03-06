module SpecConst
  VIRUS_1 = "AJ609634.snazzy.lala.long_name"
  VIRUS_2 = "DQ113772"

  VIRUS_1_NEW_NAME = "apple pie good"
  VIRUS_2_NEW_NAME = "and reaaly tasty"

  HOST_1 = "NC_002971"
  HOST_2 = "NC_002163.snazzy.lala.long_name"

  HOST_1_NEW_NAME = "ryan moore"
  HOST_2_NEW_NAME = "amelia harrison"

  PROGRAM_1 = "VirHostMatcher"
  PROGRAM_2 = "WIsH"
  PROGRAM_3 = "mummer"
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

  # This one is for version2 of the function.
  let(:parsed_output_mummer2) do
    # Note that this one doesn't have them ordered by anything in particular.
    # These numbers varified by hand
    klass = Class.new.extend Rya::CoreExtensions::Math

    # Scale to the overall max, and 0, the theoretical min.

    { # H1 is 2971, H2 is 2163
      "gi|38638842|emb|AJ609634.1|" => [# AJ
        { host: "gi|15791399|ref|NC_002163.1|", score: 18, scaled_score: klass.scale(18, 0, 22, 1, 0) },
        { host: "gi|77358712|ref|NC_002971.3|", score: 19, scaled_score: klass.scale(19, 0, 22, 1, 0) },
      ],
      "gi|73747015|gb|DQ113772.1|"  => [# DQ
        { host: "gi|15791399|ref|NC_002163.1|", score: 22, scaled_score: klass.scale(22, 0, 22, 1, 0) },
        { host: "gi|77358712|ref|NC_002971.3|", score: 18, scaled_score: klass.scale(18, 0, 22, 1, 0) },
      ],
    }
  end
  let(:parsed_output_mummer) do
    # Note that this one doesn't have them ordered by anything in particular.
    # These numbers varified by hand
    klass = Class.new.extend Rya::CoreExtensions::Math

    # Scale to the overall max, and 0, the theoretical min.
    #
    # host_3 won't show up in the mummer output file.
    #

    multiplier = 1000

    { # H1 is 2971, H2 is 2163
      "virus_1" => [# AJ
        { host: "host_1", score: 19 * multiplier, scaled_score: klass.scale(19 * multiplier, 0, 22 * multiplier, 1, 0) }, # virus_1 to host_1 and others
        { host: "host_2", score: 19 * multiplier, scaled_score: klass.scale(19 * multiplier, 0, 22 * multiplier, 1, 0) }, # this one is from virus_1_reverse to host_2
        { host: "host_3", score: 0, scaled_score: 1.0 },
      ],
      "virus_2" => [# DQ
        { host: "host_1", score: 18 * multiplier, scaled_score: klass.scale(18 * multiplier, 0, 22 * multiplier, 1, 0) },
        { host: "host_2", score: 22 * multiplier, scaled_score: klass.scale(22 * multiplier, 0, 22 * multiplier, 1, 0) },
        { host: "host_3", score: 0, scaled_score: 1.0 },
      ],
      "virus_3" => [ # A fake one with no hosts and also not even in the mummer output
        { host: "host_1", score: 0, scaled_score: 1.0 },
        { host: "host_2", score: 0, scaled_score: 1.0 },
        { host: "host_3", score: 0, scaled_score: 1.0 },
      ]
    }
  end
  let(:parsed_output_homology) do
    klass = Class.new.extend Rya::CoreExtensions::Math
    multiplier = 1000

    min, max, from, to = 0, 234.0 * multiplier, 1, 0

    {
      "virus_1" => [# AJ
        { host: "host_1", score: 0, scaled_score: klass.scale(0, min, max, from, to) },
        { host: "host_2", score: 234.0 * multiplier, scaled_score: klass.scale(234.0 * multiplier, min, max, from, to) },
        { host: "host_3", score: 0, scaled_score: 1.0 },
      ],
      "virus_2" => [# DQ
        { host: "host_1", score: 39.9 * multiplier, scaled_score: klass.scale(39.9 * multiplier, min, max, from, to) },
        { host: "host_2", score: 201.4 * multiplier, scaled_score: klass.scale(201.4 * multiplier, min, max, from, to) },
        { host: "host_3", score: 0, scaled_score: 1.0 },
      ],
      "virus_3" => [ # A fake one with no hosts and also not even in the mummer output
        { host: "host_1", score: 0, scaled_score: 1.0 },
        { host: "host_2", score: 0, scaled_score: 1.0 },
        { host: "host_3", score: 0, scaled_score: 1.0 },
      ]
    }
  end
  let(:fake_seq_lengths) do
    {
      "virus_1" => 1,
      "virus_2" => 1,
      "virus_3" => 1,
      "host_1" => 0,
      "host_2" => 0,
      "host_3" => 0,
    }
  end


  let(:collated_host_table) do
    {
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
  end

  let :collated_host_table_mapped do
    {
      SpecConst::VIRUS_1_NEW_NAME => {
        SpecConst::HOST_1_NEW_NAME => {
          scores:        {
            SpecConst::PROGRAM_1 => 0.365838,
            SpecConst::PROGRAM_2 => -1.38868,
          },
          scaled_scores: {
            SpecConst::PROGRAM_1 => 0.365838,
            SpecConst::PROGRAM_2 => BigSimon::Utils.scale_log_likelihood(-1.38868),
          },
        },
        SpecConst::HOST_2_NEW_NAME => {
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
      SpecConst::VIRUS_2_NEW_NAME => {
        SpecConst::HOST_1_NEW_NAME => {
          scores:        {
            SpecConst::PROGRAM_1 => 0.360301,
            SpecConst::PROGRAM_2 => -1.38371,
          },
          scaled_scores: {
            SpecConst::PROGRAM_1 => 0.360301,
            SpecConst::PROGRAM_2 => BigSimon::Utils.scale_log_likelihood(-1.38371),
          },
        },
        SpecConst::HOST_2_NEW_NAME => {
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
  end


  describe BigSimon::Pipeline do
    describe "::collate_host_results" do
      it "raises if length of inputs doesn't match" do
        expect {
          BigSimon::Pipeline::collate_host_results [Hash.new], ["", ""]
        }.to raise_error Rya::AbortIf::Assert::AssertionFailureError
      end

      it "collates host results" do
        actual = BigSimon::Pipeline.collate_host_results [parsed_output_vhm, parsed_output_wish], [SpecConst::PROGRAM_1, SpecConst::PROGRAM_2]

        expect(actual).to eq collated_host_table
      end
    end

    describe "::map_taxa" do
      it "maps viral and host taxa in collated results to original name" do
        virus_name_map = {
          SpecConst::VIRUS_1 => SpecConst::VIRUS_1_NEW_NAME,
          SpecConst::VIRUS_2 => SpecConst::VIRUS_2_NEW_NAME,
        }
        host_name_map  = {
          SpecConst::HOST_1 => SpecConst::HOST_1_NEW_NAME,
          SpecConst::HOST_2 => SpecConst::HOST_2_NEW_NAME,
        }

        actual = BigSimon::Pipeline.map_taxa collated_host_table, virus_name_map, host_name_map

        expect(actual).to eq collated_host_table_mapped
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

    # describe "::mummer2" do
    #   let(:program_name) { SpecConst::PROGRAM_3 }
    #
    #   it "calculates matches" do
    #     actual_results = nil
    #     expect {
    #       actual_results = BigSimon::Runners.mummer exe, vir_dir, host_dir, outdir, threads
    #     }.not_to raise_error
    #
    #     expect(actual_results).to eq parsed_output_mummer2
    #   end
    # end

    describe "::homology" do
      # Need this for some of the output variables.
      let(:program_name) { "homology_seach" }

      vdir = File.join BigSimon::TEST_FILES, "homology_files", "virus"
      hdir = File.join BigSimon::TEST_FILES, "homology_files", "host"
      odir = File.join BigSimon::TEST_FILES, "homology_files", "output"


      it "runs the homology pipeline" do
        actual_output = nil
        expected_output = parsed_output_homology

        expect {
          actual_output = BigSimon::Runners.homology vdir, hdir, odir, threads, fake_seq_lengths
        }.not_to raise_error

        expect(actual_output).to eq expected_output
      end
    end

    describe "::mummer" do
      let(:program_name) { SpecConst::PROGRAM_3 }

      vdir = File.join BigSimon::TEST_FILES, "mummer_files", "virus"
      hdir = File.join BigSimon::TEST_FILES, "mummer_files", "host"
      odir = File.join BigSimon::TEST_FILES, "mummer_files", "output"


      it "calculates matches" do
        actual_results = nil
        expect {
          actual_results = BigSimon::Runners.mummer exe, vdir, hdir, odir, threads, fake_seq_lengths
        }.not_to raise_error

        expect(actual_results).to eq parsed_output_mummer
      end
    end


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

    describe "::heatmaps" do
      it "makes heatmaps" do
        indir  = File.join BigSimon::TEST_FILES, "outdir_for_heatmaps"
        outdir = File.join indir, "outdir"

        FileUtils.rm_r outdir if Dir.exist? outdir

        actual_outfiles = nil
        expect {
          actual_outfiles = BigSimon::Runners.heatmaps BigSimon::RSCRIPT, indir, outdir
        }.not_to raise_error

        expected_outfiles = []
        Dir.glob("#{indir}/*.txt").each do |infile|
          ext  = File.extname infile
          base = File.basename infile, ext

          outfile = File.join outdir, "#{base}.heatmap.pdf"
          expected_outfiles << File.absolute_path(outfile)

          expect(File).to exist outfile
        end

        expect(actual_outfiles).to eq expected_outfiles
      end
    end
  end
end
