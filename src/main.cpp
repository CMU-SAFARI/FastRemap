// C++ Implementation of CrossMap 
// For speed and improved accuracy
// by: Jeremie Kim

#include "utils.h" 
#include "mapbam.h" 
#include "IntervalTree.h" 
#include "common.h" 

#include<stdio.h> 
#include<iostream> 
#include<tuple> 
#include<cstring> 
#include<ctype.h> 
#include<stdlib.h> 
#include<unistd.h>
#include<getopt.h> 

int main(int argc, char **argv) {
	// if enough arguments, check if the second == "bam" 
 	// chain file, in_file, unmapped_file, outfile 

	std::string chain_file, in_file, unmapped_file, out_file;  

    int insert_size = 200; 
    int insert_size_stdev = 30; 
    int insert_size_fold = 3; 
    bool addtags = false; 

	const char* const short_opts = "s:t:m:ha";
    const option long_opts[] = {
			{"mean", required_argument, nullptr, 'm'},
			{"stdev", required_argument, nullptr, 's'},
			{"times", required_argument, nullptr, 't'},
			{"append-tags", no_argument, nullptr, 'a'}, 
            {"help", no_argument, nullptr, 'h'},
            {nullptr, no_argument, nullptr, 0}
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (opt == -1) { 
            break;
		} 

        switch (opt) {
			case 'm':
				insert_size = std::stoi(optarg); 
				break; 
			case 's':
				insert_size_stdev = std::stoi(optarg); 
				break; 
			case 't': 
				insert_size_fold = std::stoi(optarg); 
				break; 
			case 'a': 
				addtags = true; 
				break; 

			case 'h': // -h or --help
			case '?': // Unrecognized option
			default:
				break;
        }
    }

    int num_args = argc - optind; 
	if (num_args >= 1) { 
		if (strcmp(argv[optind], "bam") == 0) { 

			if (argc >= 3) { 
				// TODO: add optional arguments 
				chain_file = std::string(argv[optind+1]);
				in_file = std::string(argv[optind+2]); 
				unmapped_file = std::string(argv[optind+3]);  
				if (argc >= 4) { 
					out_file = std::string(argv[optind+4]); 
				} 

				std::map<std::string, int> target_chrom_size; 
				std::map<std::string, int> source_chrom_size; 
				std::map<std::string, ITree> mapTree; 
				read_chain_file(chain_file, target_chrom_size, source_chrom_size, mapTree); 

				std::cout << "Input File:    " << in_file << "\n";  
				std::cout << "Unmapped File: " << unmapped_file << "\n";  
				std::cout << "Output File:   " << out_file << "\n"; 
				std::cout << "mean : " << insert_size << "\n"; 
				std::cout << "addtags: " << addtags << "\n"; 



				crossmap_bam_file(mapTree, chain_file, in_file, unmapped_file, out_file, target_chrom_size, insert_size, insert_size_stdev, insert_size_fold, addtags); 
			} 
		} 
	} 

	return 0; 
} 
