#ifndef __MAPBAM_H 
#define __MAPBAM_H

#include "common.h" 

#include<map> 
#include<string> 
#include<stdio.h> 

int crossmap_bam_file(std::map<std::string, ITree>&, std::string, std::string, std::string, std::string, std::map<std::string, int>, int, int, int, bool); 

#endif 
