#ifndef __MAPBED_H 
#define __MAPBED_H

#include "common.h" 

#include<map> 
#include<string> 
#include<stdio.h> 

int crossmap_bed_file(std::map<std::string, ITree>& mapping, std::string chainfile, std::string infile, std::string unmapped_file, std::string outfile_prefix); 

#endif 
