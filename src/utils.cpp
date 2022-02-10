#include "utils.h"
#include "common.h" 

#include<iostream> 
#include<map> 
#include<fstream> 
#include<sstream> 
#include<stdio.h> 
#include<string> 
#include<algorithm> 
//#include<seqan3/alphabet/views/all.hpp>
//#include<seqan3/alphabet/cigar/cigar.hpp> 
#include<stdio.h> 
#include<regex> 

using namespace std; 
const std::string WHITESPACE = " \n\r\t\f\v";


int read_chain_file(std::string chain_file, std::map<std::string, int>& target_chromSize, std::map<std::string, int>& source_chromSize, std::map<std::string, ITree>& maps) { 
    std::string source_name,   target_name;
    std::string source_strand, target_strand; 
    int source_size,           target_size; 
    int source_start,          target_start; 
    int sfrom,                 tfrom; 
    int sgap,                  tgap; 
    int size; 

    std::map<std::string, ITree::interval_vector> interval_vector_map; 

    std::cout << "Reading the chain file: " << chain_file << "\n"; 

    std::ifstream infile(chain_file); 
    if (infile.is_open()) { 
        std::string line; 
        int line_number = 0; 
        while (std::getline(infile, line)) {
            trim(line); 

            if (line.empty()) continue; 
            if (line.at(0) == '#') continue; 

            vector<string> fields; 
            std::stringstream ss(line); 
            std::string token; 
            while(ss >> token) { 
                fields.push_back(token); 
            } 

            // found a chainfile header line 
            if (fields[0] == "chain" && (fields.size() == 12 || fields.size() == 13)) { 
                source_name = fields[2]; 
                source_size = stoi(fields[3]); 
                source_strand = fields[4]; 
                source_start = stoi(fields[5]); 
                target_name = fields[7]; 
                target_size = stoi(fields[8]); 
                target_strand = fields[9]; 
                target_start = stoi(fields[10]); 
                target_chromSize[target_name] = target_size; 
                source_chromSize[source_name] = source_size; 

                if (source_strand != "+") {
                    std::cout << "ERROR: source_strand must be + (Line: " << line_number << ")\n"; 
                    exit(0); 
                }
                if (target_strand != "+" and target_strand != "-") {
                    std::cout << "ERROR: target_strand must be + or - (Line: " << line_number << ")\n"; 
                    exit(0); 
                }

                if (interval_vector_map.find(source_name) == interval_vector_map.end()) { 
                    interval_vector_map[source_name] = ITree::interval_vector();
                } 

                sfrom = source_start;
                tfrom = target_start; 
            }
            else if (fields[0] != "chain" and fields.size() == 3) {
                size = stoi(fields[0]); 
                sgap = stoi(fields[1]); 
                tgap = stoi(fields[2]); 
                
                if (target_strand == "+") { 
                    interval_vector_map[source_name].push_back(ITree::interval(sfrom, sfrom+size, make_tuple(target_name, tfrom, tfrom+size, target_strand))); 
                } 
                else if (target_strand == "-") { 
                    interval_vector_map[source_name].push_back(ITree::interval(sfrom, sfrom+size, make_tuple(target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))); 
                } 
                sfrom += size + sgap;
                tfrom += size + tgap; 
            } 
            else if (fields[0] != "chain" and fields.size() == 1) { 
                size = stoi(fields[0]); 
                
                if (target_strand == "+") { 
                    interval_vector_map[source_name].push_back(ITree::interval(sfrom, sfrom+size, make_tuple(target_name, tfrom, tfrom+size, target_strand))); 
                }
                else if (target_strand == "-") { 
                    interval_vector_map[source_name].push_back(ITree::interval(sfrom, sfrom+size, make_tuple(target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))); 
                } 
            } 
            else {
                std::cout << "ERROR: invalid chain file format (Line: " << line_number << ")\n"; 
                exit(0); 
            } 

            line_number++; 
        } 
        infile.close();  
    } 

    map<std::string, ITree::interval_vector>::iterator it; 
    for (it = interval_vector_map.begin(); it != interval_vector_map.end(); it++) { 
        maps[it->first] = ITree(std::move(it->second), 16, 1); 
    } 

    return 0; 
} 


int map_coordinates(std::map<std::string, ITree>& mapping, std::string q_chr, int q_start, int q_end, std::vector<std::tuple<std::string, int, int, std::string>>& matches, std::string q_strand = "+", bool print_match = false) { 
    // initialize matches by clearing everything 
    matches.clear();     
    std::map<std::string, std::string> complement = {{"+" , "-"}, {"-" , "+"}};
    
    ITree::interval_vector targets; 
    std::string tmp_q_chr = std::regex_replace(q_chr, std::regex("chr"), ""); 
    std::string tmp2_q_chr = "chr" + q_chr; 
//    std::cout << "q_chr: " << q_chr; 
//    std::cout << "tmp_q_chr: " << tmp_q_chr; 
//    std::cout << "tmp2_q_chr: " << tmp2_q_chr; 
    std::string mapping_query; 
    if (mapping.find(q_chr) != mapping.end()) { 
        mapping_query = q_chr; 
    } 
    else if (mapping.find(tmp_q_chr) != mapping.end()) { 
        mapping_query = tmp_q_chr; 
    } 
    else if (mapping.find(tmp2_q_chr) != mapping.end()) { 
        mapping_query = tmp2_q_chr; 
    }
    else {
//        std::cout << "INSIDE map_coordinates() no matching query name found\n";
        return -1; 
    } 

	targets = mapping[mapping_query].findOverlapping(q_start, q_end); 

    if (targets.size() == 0) { 
  //      std::cout << "INSIDE map_coordinates() no overlapping interval found\n";
        return -1; 
    } 
    else if (targets.size() >= 1) { 
        
        for (ITree::interval_vector::iterator it = targets.begin(); it != targets.end(); it++) {
            //std::cout << "checking\n"; 
            int s_start = it->start; 
            int s_end = it->stop; 
            std::string t_chrom = std::get<0>(it->value); 

            t_chrom = update_chromID(q_chr, t_chrom); 
            int t_start = std::get<1>(it->value); 
            int t_end = std::get<2>(it->value); 
            std::string t_strand = std::get<3>(it->value); 
            
            std::string chr; 
            int real_start, real_end; 
            intersectBed(q_chr, q_start, q_end, q_chr, s_start, s_end, chr, real_start, real_end); 

            int l_offset = abs(real_start - s_start); 
            int size = abs(real_end - real_start); 

            //std::cout << "NEW\n"; 
            //std::cout << s_start << "\n" << s_end << "\n" << t_chrom << "\n" << t_start << "\n" << t_end << "\n" << t_strand << "\n\n"; 

            
            matches.push_back(std::make_tuple(chr, real_start, real_end, q_strand)); 
            if (t_strand == "+") { 
                int i_start = t_start + l_offset; 
                if (q_strand == "+") {
                    matches.push_back(std::make_tuple(t_chrom, i_start, i_start + size, t_strand)); 
                }
                else {
                    matches.push_back(std::make_tuple(t_chrom, i_start, i_start + size, complement[t_strand])); 
                } 
            } 
            else if (t_strand == "-") {  // TODO: CHECK VALIDITY HERE? why complement when both - - ? 
                int i_start = t_end - l_offset - size; 
                if (q_strand == "+") {
                    matches.push_back(std::make_tuple(t_chrom, i_start, i_start + size, t_strand)); 
                }
                else {
                    matches.push_back(std::make_tuple(t_chrom, i_start, i_start + size, complement[t_strand])); 
                } 
            } 
            else {
                std::cout << "Unknown strand " << q_strand << ". Can only be + or -.\n";  
                exit(1); 
            }
        } 
    } 
    //std::cout << "matches \n"; 
    //for (int it = 0; it < matches.size(); it++) { 
    //    std::cout << std::get<0>(matches[it]) << " " << std::get<1>(matches[it])  << " " << std::get<2>(matches[it])  << " " << std::get<3>(matches[it]) << "\n"; 
    //} 
    return 0; 
} 

std::string update_chromID(std::string c_temp, std::string c_target) { 
    if (c_temp.find("chr") == 0) { 
        if (c_target.find("chr") == 0) { 
            return c_target; 
        }
        else { 
            return "chr" + c_target; 
        } 
    } 
    else {
        if (c_target.find("chr") == 0) { 
            return std::regex_replace(c_target, std::regex("chr"), ""); 
        } 
        else {
            return c_target; 
        } 
    } 
    return ""; 
} 

bool check_bed12(std::string bedline){
    // Check if bed12 format is correct or not.

    std::vector<std::string> fields; 
    std::stringstream ss(bedline); 
    std::string token; 
    while(ss >> token) { 
        fields.push_back(token); 
    }

    if(fields.size() != 12)
        return false;
    
    if((fields[5] != "+") && (fields[5] != "-") && (fields[5] != ".")){
        return false;
    }

    int chromStart;
    int chromEnd;
	int thickStart;
	int thickEnd;
	int blockCount;

    std::stringstream blockStartsStream(fields[11]);
    std::stringstream blockSizesStream(fields[10]);
    std::vector<int> blockStarts;
    std::vector<int> blockSizes;
    std::string item;

    try{	 
        chromStart = std::stoi(fields[1]);
		chromEnd = std::stoi(fields[2]);
		thickStart = std::stoi(fields[6]);
		thickEnd = std::stoi(fields[7]);
		blockCount = std::stoi(fields[9]);

        while(std::getline(blockStartsStream, item, ','))
            blockStarts.push_back(std::stoi(item));
        while(std::getline(blockSizesStream, item, ','))
            blockSizes.push_back(std::stoi(item));         
    }
    catch(std::exception &err){
        return false;
    }

    if((chromStart > chromEnd) || (thickStart > thickEnd))
        return false;
    if((thickStart < chromStart) || (thickEnd > chromEnd))
        return false;
    if(blockSizes.size() != blockCount)
        return false;
    if(blockStarts.size() != blockCount)
        return false;
	if(blockCount < 1)
        return false;

    for(int i = 0; i < blockSizes.size(); i++){
        if(blockSizes[i] < 0)
            return false;
    }

    for(int i = 0; i < blockStarts.size(); i++){
        if(blockStarts[i] < 0)
            return false;
    }
    return true;
}

int intersectBed(std::string chr1, int st1, int end1, std::string chr2, int st2, int end2, std::string& ret_chr, int& ret_st, int& ret_end) { 
    
    if (st1 > end1 || st2 > end2) { 
        std::cout << "Start cannot be larger than end\n";  
        exit(1); 
    } 
    if (chr1 != chr2) {
        return -1; 
    } 
    if (st1 > end2 || end1 < st2) { 
        return -1; 
    } 

    ret_chr = chr1; 
    ret_st = std::max(st1, st2); 
    ret_end = std::min(end1, end2); 
    
    return 0; 
}

std::string revcomp_DNA(std::string dna, bool extended){
    std::map<std::string, std::string> extended_map = {
                { "A", "T"}, { "C", "G"}, { "G", "C"}, { "T", "A"},
                { "Y", "R"}, { "R", "Y"}, { "S", "W"}, { "W", "S"},
                { "K", "M"}, { "M", "K"}, { "B", "V"}, { "V", "B"},
                { "D", "H"}, { "H", "D"}, { "N", "N"}, { ".", "."}, { "*", "*"}};
    std::map<std::string, std::string> not_extended_map = {
                { "A", "T"}, { "C", "G"}, { "G", "C"}, { "T", "A"},
                { "N", "N"}, { "X", "X"}};

	// seq = dna.replace(' ','').upper()

    std::string reversed_dna = "";
    std::string current_rev = "";
    if(extended){
        for(int i = 0; i < dna.length(); i++){
            if(dna.substr(i, 1) == " ")
                ;
            else if(dna.substr(i, 1) == ","){
                reversed_dna = reversed_dna + "," + current_rev;
                current_rev = "";
            }
            else
                current_rev = extended_map.at(dna.substr(i, 1));
        }
    }
    else{
        for(int i = 0; i < dna.length(); i++){
            if(dna.substr(i, 1) == " ")
                ;
            else if(dna.substr(i, 1) == ","){
                reversed_dna = reversed_dna + "," + current_rev;
                current_rev = "";
            }
            else
                current_rev = not_extended_map.at(dna.substr(i, 1));
        }   
    }
    reversed_dna = reversed_dna + "," + current_rev;
    reversed_dna = reversed_dna.substr(1); // to remove the , in the beginning

    return reversed_dna;
}

int get_reference_length(seqan::String<seqan::CigarElement<>> cigar) { 
    int len = 0; 

    for (auto it : cigar) { 
        auto count = it.count; 
        auto operation = it.operation; 
        if (operation == (char)('D') || operation == (char)('M')) { 
            len += count; 
        }  
    } 
    return len; 
} 

 
// string trimming functions
// https://www.techiedelight.com/trim-string-cpp-remove-leading-trailing-spaces/ 
std::string lefttrim(const std::string &s) {
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string righttrim(const std::string &s) {
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}
 
std::string trim(const std::string &s) {
    return righttrim(lefttrim(s));
}



