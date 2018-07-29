test_small:
	rm -r 0000TEST/; exe/big_simon -v spec/test_files/virus/* -h spec/test_files/host/* -o 0000TEST && tree 0000TEST

test_small_install:
	rm -r 0000TEST/; rake install && exe/big_simon -v spec/test_files/virus/* -h spec/test_files/host/* -o 0000TEST && tree 0000TEST

test_toy:
	rm -r toyexample_out; time exe/big_simon -v vendor/repos/VirHostMatcher/test/toyexample/virus/* -h vendor/repos/VirHostMatcher/test/toyexample/host/* -o toyexample_out -t 3
