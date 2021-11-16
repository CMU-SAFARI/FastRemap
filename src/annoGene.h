#ifndef __ANNOGENE_H 
#define __ANNOGENE_H

#include "common.h" 

#include<map> 
#include<string> 
#include<stdio.h> 
 
void getExonFromLine(std::vector<std::string> fields, std::vector<std::tuple<std::string, int, int>>& list);
#endif 
