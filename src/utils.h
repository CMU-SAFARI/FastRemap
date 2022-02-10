#ifndef __UTILS_H 
#define __UTILS_H 

#include "common.h" 
#include<iostream>
#include<string> 
#include<map> 
//#include<seqan3/alphabet/cigar/cigar.hpp> 
#include<seqan/bam_io.h> 


int read_chain_file(std::string, std::map<std::string, int>&, std::map<std::string, int>&, std::map<std::string, ITree>&); 
int map_coordinates(std::map<std::string, ITree>&, std::string, int, int, std::vector<std::tuple<std::string, int, int, std::string>>&, std::string, bool); 
bool check_bed12(std::string bedline);
int intersectBed(std::string, int, int, std::string, int, int, std::string&, int&, int&); 
std::string update_chromID(std::string c_temp, std::string c_target); 
std::string revcomp_DNA(std::string dna, bool extended);
int get_reference_length(seqan::String<seqan::CigarElement<>>); 

// string trimming functions
// https://www.techiedelight.com/trim-string-cpp-remove-leading-trailing-spaces/ 
std::string righttrim(const std::string&);
std::string lefttrim(const std::string&);
std::string trim(const std::string&); 

#endif 
