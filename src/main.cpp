// C++ Implementation of CrossMap 
// For speed and improved accuracy
// by: Jeremie Kim

#include "utils.h" 
#include "mapbam.h" 
#include "mapsam.h" 
#include "mapbed.h"
#include "IntervalTree.h" 
#include "common.h" 
#include "CommandLineParser.h"

#include<stdio.h> 
#include<iostream> 
#include<tuple> 
#include<cstring> 
#include<ctype.h> 
#include<stdlib.h> 
#include<unistd.h>
#include<getopt.h> 

int main(int argc, const char **argv) {
	// parse arguments / parameters
	CommandLineParser options;
    seqan::ArgumentParser::ParseResult parseRes = parseCommandOptions(options, argc, argv);

    if(parseRes != seqan::ArgumentParser::PARSE_OK){
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    std::cout << "Chain File:    " << std::string(toCString(options.chain_file)) << "\n";  
	std::cout << "Input File:    " << std::string(toCString(options.in_file)) << "\n";  
	std::cout << "Unmapped File: " << std::string(toCString(options.unmapped_file)) << "\n";  
	std::cout << "Output File:   " << std::string(toCString(options.out_file)) << "\n"; 
	if(strcmp(toCString(options.type), "bed") != 0){
		std::cout << "mean : " << options.insert_size << "\n";
		std::cout << "stdev : " << options.insert_size_stdev << "\n";
		std::cout << "times : " << options.insert_size_fold << "\n";
		std::cout << "addtags: " << options.addtags << "\n"; 
	}

    std::map<std::string, int> target_chrom_size; 
	std::map<std::string, int> source_chrom_size; 
	std::map<std::string, ITree> mapTree; 
	read_chain_file(std::string(toCString(options.chain_file)), target_chrom_size, source_chrom_size, mapTree); 

	// parsing BAM file 
	if(strcmp(toCString(options.type), "bam") == 0){  
		crossmap_bam_file(mapTree, std::string(toCString(options.chain_file)), std::string(toCString(options.in_file)), std::string(toCString(options.unmapped_file)), std::string(toCString(options.out_file)), target_chrom_size, options.insert_size, options.insert_size_stdev, options.insert_size_fold, options.addtags);
	}
	// parsing SAM file 
	else if(strcmp(toCString(options.type), "sam") == 0){ 
		crossmap_sam_file(mapTree, std::string(toCString(options.chain_file)), std::string(toCString(options.in_file)), std::string(toCString(options.unmapped_file)), std::string(toCString(options.out_file)), target_chrom_size, options.insert_size, options.insert_size_stdev, options.insert_size_fold, options.addtags);
	}
	// parsing BED file 
	else if(strcmp(toCString(options.type), "bed") == 0){  
		crossmap_bed_file(mapTree, std::string(toCString(options.chain_file)), std::string(toCString(options.in_file)), std::string(toCString(options.unmapped_file)), std::string(toCString(options.out_file)));
	}

	return 0; 
} 
