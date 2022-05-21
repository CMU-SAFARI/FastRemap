#ifndef __MAPVCF_H 
#define __MAPVCF_H

#include "common.h" 

#include<map> 
#include<string> 
#include<stdio.h> 

void crossmap_vcf_file(std::map<std::string, ITree>& mapping, std::string infile, std::string outfile, std::string refgenome, bool noCompAllele = false, bool compress = false);

#endif 
