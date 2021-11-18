#include "mapbed.h" 
#include "common.h" 
#include "utils.h" 
#include "annoGene.h"

#include<algorithm> 
#include<string> 
#include<stdio.h>
#include<iostream>
#include<sstream> 
#include<map> 
#include<vector> 
#include <set>
#include<filesystem> 

int crossmap_bed_file(std::map<std::string, ITree>& mapping, std::string chainfile, std::string infile, std::string unmapped_file, std::string outfile_prefix) { 
    // determine input file format (BAM, CRAM or SAM) 
    std::string file_type; 
    std::string file_extension = infile.substr(infile.find_last_of(".") + 1); 
    std::vector<std::string> comments; 
    if (file_extension == "bed") { 
        file_type = "BED"; 
        comments.push_back("ORIGINAL_BED_FILE=" + std::string(std::filesystem::current_path()) + "/" + infile); 
        // TODO: check header. 
    } 
    else {
        std::cout << "File Extension: " << file_extension << " not supported or unknown\n"; 
        exit(1); 
    } 
    comments.push_back("CHAIN_FILE=" + std::string(std::filesystem::current_path()) + "/" + chainfile); 

    // get output file names and open the output file  
    std::ofstream OUTFILE;
    std::string outfile_name = outfile_prefix;
    if (outfile_prefix != "") { 
        if (file_type == "BED") {
            outfile_name = outfile_name + ".bed"; 
            OUTFILE.open(outfile_name);
        } 
    } 
    else { 
        std::cout << "Output to Screen not yet implemented\n"; 
        exit(1); 
    } 

	/////////////////////////////////////////////////////////////////////////
    // open the unmapped file
    std::ofstream UNMAP; 
    if (unmapped_file != "") { 
        UNMAP.open(unmapped_file);
    }
    else {
        std::cout << "unmapped file does not exist\n"; 
        exit(1); 
    }


    // open and read the input bed file
    std::ifstream INFILE(infile); 
    std::string line; 
    std::vector<std::tuple<std::string, int, int, std::string>> matches;

    if (INFILE.is_open()) {    
        while (std::getline(INFILE, line)) {
            trim(line); 
            if (line.empty()) continue; 
            if (line.at(0) == '#') continue; 
            if (line.substr(0, 5) == "track") continue; 
            if (line.substr(0, 7) == "browser") continue; 

            std::vector<std::string> fields; 
            std::stringstream ss(line); 
            std::string token; 
            while(ss >> token) { 
                fields.push_back(token); 
            }

            // filter out line less than 3 columns
            if(fields.size() < 3){
                // std::cout << "ERROR: Less than 3 fields. skip (Line: " << line << ")\n"; 
                UNMAP << line << "\tInvalidBedFormat.\n"; 
			    continue;
            }

            // try to extract chromosome number, start & end locations
            std::string chrom;
			int start;
			int end;
            std::string strand = "+";
            
            chrom = fields[0];
            try{
                start = std::stoi(fields[1]);	           
            }
            catch(std::exception &err){
                UNMAP << line << "\tInvalidStartPosition.\n";
                continue;
            }

            try{
                end = std::stoi(fields[2]);	           
            }
            catch(std::exception &err){
                UNMAP << line << "\tInvalidEndPosition.\n";
                continue;
            }

            if(start > end){
                UNMAP << line << "\tStart>End.\n";
                continue;
            }
                    
            // deal with bed less than 12 columns
            if(fields.size() < 12){
                // try to reset strand
                for(int vector_index = 0; vector_index < fields.size(); vector_index++){
                    if(fields[vector_index] == "+")
                        strand = "+";
                    else if(fields[vector_index] == "-")
                        strand = "-";
                    
                }
          //      std::cout << "chr: " << chrom << " start: " << start << " end: " << end << std::endl; 

                // find the chains this bed interval is intersecting
                int retval = map_coordinates(mapping, chrom, start, end, matches, strand, false); 
                // std::cout << "return value of map_coordinates : " << retval << std::endl;
            //    std::cout << "matches.size(): " << matches.size() << std::endl; 

                // analyze the intersecting chains
                // no match or incorrect number of matches
                if ((matches.size() == 0) || ((matches.size() % 2) != 0)) {
                    //std::cout << line << "\tinterval failed to liftover\n"; 
                    UNMAP << line << "\tUnmap.\n"; 
                    continue; 
                } 

                // interval uniquely mapped 
                else if (matches.size() == 2) { 
					// reset fields
					fields[0] = std::get<0>(matches[1]);
					fields[1] = std::to_string(std::get<1>(matches[1]));
					fields[2] = std::to_string(std::get<2>(matches[1]));
              //      std::cout << "1. fields[0]: " << fields[0] << " fields[1]: " << fields[1] << " fields[2]: " << fields[2] << std::endl; 

                    // update the strand information
                    for(int index = 0; index < fields.size(); index++){
                        if((fields[index] == "+") || (fields[index] == "-")){
                            fields[index] = std::get<3>(matches[1]);
                        }
                    }
                //    std::cout << "2. fields[0]: " << fields[0] << " fields[1]: " << fields[1] << " fields[2]: " << fields[2] << std::endl; 

                    // write the updated interval to out file
                    for(int index = 0; index < fields.size(); index++){
                        OUTFILE << fields[index] << "\t"; 
                    }
                    OUTFILE << "\n"; 
                }

                // interval multi mapped 
                else if (matches.size() > 2){
                    int count = 0;
                    for(int match_index = 1; match_index < matches.size(); match_index = match_index + 2){
                        count += 1;
						fields[0] = std::get<0>(matches[match_index]);
						fields[1] = std::to_string(std::get<1>(matches[match_index]));
						fields[2] = std::to_string(std::get<2>(matches[match_index]));

                        // update the strand information
                        for(int index = 0; index < fields.size(); index++){
                            if((fields[index] == "+") || (fields[index] == "-")){
                                fields[index] = std::get<3>(matches[match_index]);
                            }
                        }
                        
                        // TODO: maybe we can report the number of matches as well, or the mapping number etc
                        // write the updated interval to out file
                        for(int index = 0; index < fields.size(); index++){
                            OUTFILE << fields[index] << "\t"; 
                        }
                        OUTFILE << "\n"; 

                    }
                }
				    
            }

            // WARNING: there might be an issue regarding int-string conversion, check it later
            // deal with bed12 and bed12+8 (genePred format), THE SECOND ONE ACTUALLY NOT SUPPORTED
            else if((fields.size() == 12) || (fields.size() == 20)){
                bool fail_flag = false;
                strand = fields[5];
                if((strand != "+") && (strand != "-")){
                    UNMAP << line << "\tUnknown strand.\n";
                    continue;
                }

                std::vector<std::tuple<std::string, int, int>> exons_old_pos;
                std::vector<std::tuple<std::string, int, int, std::string>> exons_new_pos;
                getExonFromLine(fields, exons_old_pos);

                for(int i = 0; i < exons_old_pos.size(); i++){
                    // matches has two elements, first is query, 2nd is target. 
                    map_coordinates(mapping, std::get<0>(exons_old_pos[i]), std::get<1>(exons_old_pos[i]), 
                                    std::get<2>(exons_old_pos[i]), matches, strand, false);
                    
                    if ((matches.size() == 0) || ((matches.size() % 2) != 0)) {
                        fail_flag = true; 
                        break; 
                    }
                    else if (matches.size() == 2) { 
                        exons_new_pos.push_back(std::make_tuple(std::get<0>(matches[1]), 
                        std::get<1>(matches[1]), std::get<2>(matches[1]), std::get<3>(matches[1])));
                    }
                    else{
                        fail_flag = true;
                        break;
                    }
                }

                if(!fail_flag){
                    // check if all exons were mapped to the same chromosome and the same strand
                    std::set<std::string> chr_id;
                    std::set<int> exon_strand;

                    for(int i = 0; i < exons_new_pos.size(); i++){
                        chr_id.insert(std::get<0>(exons_new_pos[i]));
                        exon_strand.insert(std::get<1>(exons_new_pos[i]));
                    }
                    
                    if((chr_id.size() != 1) || (exon_strand.size() != 1))
                        fail_flag = true;

                    if(!fail_flag){
                        // build new bed
                        // TODO: make sure these parts are correctly transformed to c++
                        int cds_start_offset = std::stoi(fields[6]) - std::stoi(fields[1]);
                        int cds_end_offset = std::stoi(fields[2]) - std::stoi(fields[7]);
                        std::string new_chrom = std::get<0>(exons_new_pos[0]);
                        int new_chrom_st = std::get<1>(exons_new_pos[0]);
                        int new_chrom_end = std::get<2>(exons_new_pos.back());
                        std::string new_chrom_st_str = std::to_string(new_chrom_st);
                        std::string new_chrom_end_str = std::to_string(new_chrom_end);
                        std::string new_name = fields[3];
                        std::string new_score = fields[4];
                        std::string new_strand = std::get<3>(exons_new_pos[0]);
                        std::string new_thickStart = std::to_string(new_chrom_st + cds_start_offset);
                        std::string new_thickEnd = std::to_string(new_chrom_end - cds_end_offset);
                        std::string new_ittemRgb = fields[8];
                        std::string new_blockCount = std::to_string(exons_new_pos.size());
                        
                        std::string new_blockSizes = "";
                        for(int i = 0; i < exons_new_pos.size(); i++){
                            new_blockSizes = new_blockSizes + std::to_string(std::get<2>(exons_new_pos[i]) - std::get<1>(exons_new_pos[i])) + ",";
                        }

                        std::string new_blockStarts = "";
                        for(int i = 0; i < exons_new_pos.size(); i++){
                            new_blockStarts = new_blockStarts + std::to_string(std::get<1>(exons_new_pos[i]) - new_chrom_st) + ",";
                        }

                        std::vector<std::string> new_fields_list {new_chrom,new_chrom_st_str,new_chrom_end_str,new_name,new_score,new_strand,
                                            new_thickStart,new_thickEnd,new_ittemRgb,new_blockCount,new_blockSizes,new_blockStarts};

                        std::string new_bedline = "";
                        for(int i = 0; i < new_fields_list.size(); i++){
                            new_bedline = new_bedline + new_fields_list[i] + "\t";
                        }

                        if(check_bed12(new_bedline) == false){
                            fail_flag = true;
                        }
                        else{
                            OUTFILE << new_bedline << "\n"; 
                        }
                                
                    }
                }

                if(fail_flag){ 
                    UNMAP << line << "\tUnmap.\n"; 
                }       
            }
        }
    }

    return 0; 
}
