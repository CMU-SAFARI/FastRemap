#include "mapbed.h" 
#include "common.h" 
#include "utils.h" 

#include<algorithm> 
#include<string> 
#include<stdio.h>
#include<iostream>
#include<sstream> 
#include<map> 
#include<vector> 
#include<filesystem> 


void getExonFromLine(std::vector<std::string> fields, std::vector<std::tuple<std::string, int, int>>& list){
    // Extract ALL exon regions from input bed line (must be 12-column). return list of [chrom st end] in list param
    list.clear();

    /**
     * not needed as we take the fields vector as input not the actual line
    std::string line = bedline;
    trim(line);

    std::vector<std::string> fields; 
    std::stringstream ss(line); 
    std::string token; 
    while(ss >> token) { 
        fields.push_back(token); 
    }
    */

    // no need to check for the correctness of start value, already checked
    int txStart = std::stoi(fields[1]); 
    std::string chrom = fields[0];
    /** Not needed params as not utilized later
     * std::string strand = fields[5];
     * std::string geneName = fields[3];
     * score=fields[4]
     * */

    std::stringstream blockStarts(fields[11]);
    std::stringstream blockSizes(fields[10]);
    std::vector<int> exon_start;
    std::vector<int> exon_end;
    std::string item;
    int counter = 0;

    while(std::getline(blockStarts, item, ',')){
        exon_start.push_back(std::stoi(item) + txStart);
    }
    while(std::getline(blockSizes, item, ',')){
        exon_end.push_back(std::stoi(item) + exon_start[counter]);
        counter++;
    }
    // #chrom = chrom + ':' + strand

    for(int i = 0; i < exon_start.size(); i++){
        list.push_back(std::make_tuple(chrom, exon_start[i], exon_end[i])); 
    }
}

