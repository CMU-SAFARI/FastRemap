#include "mapbam.h" 
#include "common.h" 
#include "utils.h" 

#include<algorithm> 
#include<string> 
#include<stdio.h>
#include<iostream>
#include<sstream> 
#include<map> 
#include<span> 
#include<vector> 
#include<filesystem> 
#include<seqan/bam_io.h>
#include<seqan/basic.h>
#include<seqan/sequence.h> 
#include<seqan/modifier.h> 

typedef typename seqan::BamHeaderRecord::TTag TTag; 

int crossmap_bam_file(std::map<std::string, ITree>& mapping, std::string chainfile, std::string infile, std::string unmapped_file, std::string outfile_prefix, std::map<std::string, int> chrom_size, int IS_size, int IS_std, int fold, bool addtag) { 
    // determine input file format (BAM, CRAM or SAM) 
    std::string file_type; 
    std::string file_extension = infile.substr(infile.find_last_of(".") + 1); 
    std::vector<std::string> comments; 
    if (file_extension == "bam") { 
        file_type = "BAM"; 
        comments.push_back("ORIGINAL_BAM_FILE=" + std::string(std::filesystem::current_path()) + "/" + infile); 
        // TODO: check header. 
    } 
    else {
        std::cout << "File Extension: " << file_extension << " not supported or unknown\n"; 
        exit(1); 
    } 
    comments.push_back("CHAIN_FILE=" + std::string(std::filesystem::current_path()) + "/" + chainfile); 

    // get output file names.  
    std::string outfile_name = outfile_prefix;
    if (outfile_prefix != "") { 
        if (file_type == "BAM") {
            outfile_name = outfile_name + ".bam"; 
        } 
    } 
    else { 
        std::cout << "Output to Screen not yet implemented\n"; 
        exit(1); 
    } 
 	
	seqan::CharString bamFileName = infile;

    auto test_outfile_name2 = std::filesystem::current_path() / outfile_name; 
	seqan::BamFileIn bamFileIn(seqan::toCString(bamFileName)); 

	std::ofstream sam_out;
	sam_out.open(test_outfile_name2);  
	seqan::BamFileOut bamFileOut(seqan::context(bamFileIn), sam_out, seqan::Bam()); 


	/////////////////////////////////////////////////////////
	// Updating Headers 
	/////////////////////////////////////////////////////////

	seqan::BamHeader header;
	seqan::BamHeader new_header;
	seqan::readHeader(header, bamFileIn);  

	seqan::BamAlignmentRecord record;
	seqan::BamHeaderRecord header_record; 
	seqan::CharString tmp; 

    std::map<std::string, int> name_to_id;
    std::map<int, std::string> id_to_name;

	seqan::StringSet<seqan::CharString> contigNameStore; 
	seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > contigNameStoreCache(contigNameStore); 
	seqan::BamIOContext<seqan::StringSet<seqan::CharString>> bamIOContext(contigNameStore, contigNameStoreCache); 

	for (unsigned recIdx = 0; searchRecord(recIdx, header, seqan::BAM_HEADER_REFERENCE, recIdx); recIdx++) { 
		unsigned this_idx = 0; 
		
		std::string chr_name;
		int chr_length; 
		while(getTagValue(tmp, this_idx, header[recIdx])) { 
			switch (this_idx) {
				case 0: 
					chr_name = seqan::toCString(tmp); 
					break;
				case 1: 
					chr_length = atoi(seqan::toCString(tmp)); 
					break; 
			}
			this_idx++; 
		}

		seqan::assignValueById(contigLengths(context(bamFileOut)), nameToId(contigNamesCache(context(bamFileOut)), header[recIdx].tags[0].i2), chrom_size[seqan::toCString(header[recIdx].tags[0].i2)]);
		
		// TODO: This needs to be implemented in seqan to get rid of useless chromosomes in header 
		// seqan::removeValueById(contigLengths(context(bamFileOut)), nameToId(contigNamesCache(context(bamFileOut)), header[recIdx].tags[0].i2)); 
	}

	// TODO: or implement this in seqan 
	// seqan::clear(contigLengths(context(bamFileOut))); 
	
	for (auto it : chrom_size) { 
		seqan::assignValueById(contigLengths(context(bamFileOut)), nameToId(contigNamesCache(context(bamFileOut)), it.first), it.second); 
		name_to_id[it.first] = nameToId(contigNamesCache(context(bamFileOut)), it.first); 
		id_to_name[nameToId(contigNamesCache(context(bamFileOut)), it.first)] = it.first; 
	} 
	seqan::BamHeaderRecord seqRecord; 
	seqRecord.type = seqan::BAM_HEADER_PROGRAM; 
	seqan::appendValue(seqRecord.tags, TTag("ID", "FastRemap")); 
	seqan::appendValue(seqRecord.tags, TTag("VN", "1.0")); 
	seqan::appendValue(header, seqRecord); 

	seqan::BamHeaderRecord seqRecord1; 
	seqRecord1.type = seqan::BAM_HEADER_COMMENT; 
	seqan::appendValue(seqRecord1.tags, TTag("ID", comments[0]), seqan::Exact()); 
	seqan::appendValue(header, seqRecord1); 

	seqan::BamHeaderRecord seqRecord2; 
	seqRecord2.type = seqan::BAM_HEADER_COMMENT; 
	seqan::appendValue(seqRecord2.tags, TTag("ID", comments[1]), seqan::Exact()); 
	seqan::appendValue(header, seqRecord2); 

	seqan::writeHeader(bamFileOut, header);
	
	/////////////////////////////////////////////////////////////////////////

    std::ofstream UNMAP; 
    if (unmapped_file != "") { 
        UNMAP.open(unmapped_file);
    }
    else {
        std::cout << "unmapped file does not exist\n"; 
        exit(1); 
    }

    int QF = 0, NN = 0, NU = 0, NM = 0, UN = 0, UU = 0, UM = 0, MN = 0, MU = 0, MM = 0, SN = 0, SM = 0, SU = 0, total_item = 0; 
    std::string unmap_queryname; 

    int read1_start, read1_end; 
    int read2_start, read2_end; 
    std::string read1_chr;
    std::string read2_chr; 
    std::string read1_strand, read2_strand;
    std::vector<std::tuple<std::string, int, int, std::string>> read1_maps; 
    std::vector<std::tuple<std::string, int, int, std::string>> read2_maps; 
    std::string new_reference_id; 


	seqan::BamAlignmentRecord old_alignment; 
	while (!seqan::atEnd(bamFileIn)) { 

        total_item++; 
		seqan::readRecord(old_alignment, bamFileIn);
		
		seqan::BamAlignmentRecord new_alignment = old_alignment; 
		seqan::BamTagsDict tagsDict(new_alignment.tags); 
        seqan::eraseTag(tagsDict, "RG"); 
       
        auto old_flag = old_alignment.flag; 
        unmap_queryname = seqan::toCString(old_alignment.qName); 
        int old_reference_length = get_reference_length(old_alignment.cigar); 
		//std::cout << "old_reference_length: " << old_reference_length << "\n"; 
		
        
        int old_reference_start = old_alignment.beginPos; 
        int old_reference_end = old_reference_start + old_reference_length; 

        // TODO: something about rg tags.. maybe unecessary for seqan3

        //////////////////////////////
        // Pair-end sequencing
        //////////////////////////////
        if (old_flag & seqan::BAM_FLAG_MULTIPLE) { 
            //std::cout << "Read paired\n"; 
            new_alignment.flag = seqan::BAM_FLAG_MULTIPLE; 
            if (old_flag & seqan::BAM_FLAG_FIRST) { 
                new_alignment.flag |= seqan::BAM_FLAG_FIRST;
                unmap_queryname += ".1"; 
            } 
            else if (old_flag & seqan::BAM_FLAG_LAST) { 
                new_alignment.flag |= seqan::BAM_FLAG_LAST; 
                unmap_queryname += ".2"; 
            } 
            
            if (old_alignment.flag & 0x800) { 
                new_alignment.flag |= 0x800; 
            } 

            if (old_flag & seqan::BAM_FLAG_QC_NO_PASS) { 
                //std::cout << "Failed filter\n"; 
                new_alignment.flag |= seqan::BAM_FLAG_QC_NO_PASS; 
                QF++; 
                if (addtag) { 
					seqan::setTagValue(tagsDict, "QF", 0); 
                } 
                UNMAP << "chr\t" << old_reference_start << "\t" << old_reference_end << "\t" << unmap_queryname << "\n"; 
                continue; 
            } 
            //////////////////////////////////////////
            // R1 originally unmapped 
            //////////////////////////////////////////
            else if (old_flag & seqan::BAM_FLAG_UNMAPPED) { 
                //std::cout << "R1 originally unmapped\n"; 
                NU++; // not accurate. could be NU, NN, NM 
                UNMAP << "chr\t" << old_reference_start << "\t" << old_reference_end << "\t" << unmap_queryname << "\n"; 
				continue; 
            } 
            /////////////////////////////////////////
            // R1 is originally mapped 
            /////////////////////////////////////////
            else {
                //std::cout << "R1 originally mapped\n"; 
                read1_chr = id_to_name[old_alignment.rID]; 

                if (old_flag & seqan::BAM_FLAG_RC) { 
                    read1_strand = "-";    
                } 
                else {
                    read1_strand = "+"; 
                }
                read1_start = old_reference_start; 
                read1_end = old_reference_end; 
                if (map_coordinates(mapping, read1_chr, read1_start, read1_end, read1_maps, read1_strand, false) == -1) {
                    //std::cout << "failed map_coordinates\n"; 
                }    
                
                if (!(old_flag & seqan::BAM_FLAG_NEXT_UNMAPPED)) { 
                    read2_chr = id_to_name[old_alignment.rNextId];
                    if (old_flag & seqan::BAM_FLAG_NEXT_RC) { 
                        read2_strand = "-";    
                    } 
                    else {
                        read2_strand = "+"; 
                    }
                    read2_start = old_alignment.pNext; 
                    read2_end = read2_start + 1; // TODO: double check this?? 
                    map_coordinates(mapping, read2_chr, read2_start, read2_end, read2_maps, read2_strand, false); 
                }

                /////////////////////////////////////
                // R1 failed to liftover 
                ///////////////////////////////////// 
                if (read1_maps.size() == 0) {
                    //std::cout << "R1 failed to liftover\n"; 
                    UNMAP << read1_chr << "\t" << read1_start << "\t" << read1_end << "\t" << unmap_queryname << "\n"; 
                    continue; 
                } 
                /////////////////////////////////////
                // R1 uniquely mapped 
                ///////////////////////////////////// 
                else if (read1_maps.size() == 2) {
                    //std::cout << "R1 uniquely mapped\n"; 
                    if (std::get<3>(read1_maps[1]) == "-") { 
                        new_alignment.flag |= seqan::BAM_FLAG_RC; 
                    }
					// REQUIRED TO AVOID MARKDUP ISSUES when remapping across chromosomes
					new_alignment.rID = name_to_id[std::get<0>(read1_maps[1])]; 
                    new_alignment.beginPos = std::get<1>(read1_maps[1]); 

                    // opposite strand 
                    if (std::get<3>(read1_maps[0]) != std::get<3>(read1_maps[1])) { 
						seqan::reverse(new_alignment.cigar); 
						seqan::Dna5String tmp = new_alignment.seq; 
						seqan::Dna5StringReverseComplement rc(tmp); 
						new_alignment.seq = seqan::CharString(rc); 

						tmp = new_alignment.qual; 
						seqan::Dna5StringReverse r(tmp); 
						new_alignment.qual = seqan::CharString(r); 
                    } 
                    
                    // R2 unmapped before or after conversion 
                    if ((old_flag & seqan::BAM_FLAG_NEXT_UNMAPPED) || (read2_maps.size() == 0)) { 
                        new_alignment.flag |= seqan::BAM_FLAG_NEXT_UNMAPPED; 
                        new_alignment.rNextId = name_to_id[std::get<0>(read1_maps[1])];  // TODO double check. 
                        new_alignment.pNext = std::get<1>(read1_maps[1]);
                        new_alignment.tLen = 0; 
                        UN += 1; 
                        if (addtag) { 
							seqan::setTagValue(tagsDict, "UN", 0); 
                        }
                    } 
                    // R2 is unique mapped 
                    else if (read2_maps.size() == 2) { 
                        if (std::get<3>(read2_maps[1]) == "-") {
                            new_alignment.flag |= seqan::BAM_FLAG_NEXT_RC; 
                        } 
                        new_alignment.rNextId = name_to_id[std::get<0>(read2_maps[1])]; 
                        new_alignment.pNext = std::get<1>(read2_maps[1]); 
                        
                        new_alignment.tLen = abs(new_alignment.beginPos - int(new_alignment.pNext)) + old_reference_length; //(old_alignment.template_length() - old_alignment.mate_position().value_or(0)); // TODO: double check  old reference_length 
                        if (std::get<3>(read2_maps[1]) != std::get<3>(read1_maps[1]) && (new_alignment.tLen <= IS_size + fold * IS_std) && (new_alignment.tLen >= IS_size - fold * IS_std)) { 
                            new_alignment.flag |= seqan::BAM_FLAG_ALL_PROPER; 
                        } 
                        UU++; 
                        if (addtag) { 
							seqan::setTagValue(tagsDict, "UU", 0); 
                        }
                    } 
                    // R2 is multiple mapped 
                    else { 
                        if (std::get<3>(read2_maps[read2_maps.size()-1]) == "-") { 
                            new_alignment.flag |= seqan::BAM_FLAG_NEXT_UNMAPPED; 
                        }
                        new_alignment.flag |= seqan::BAM_FLAG_SECONDARY; 
                        new_alignment.rNextId = name_to_id[std::get<0>(read2_maps[read2_maps.size()-1])]; 
                        new_alignment.pNext = std::get<1>(read2_maps[read2_maps.size()-1]);
                        new_alignment.tLen = 0; 
                        UM++; 
                        if (addtag) { 
							seqan::setTagValue(tagsDict, "UM", 0); 
                        }
                    } 
                    
					seqan::writeRecord(bamFileOut, new_alignment); 
                    continue; 
                } 
        
                /////////////////////////////////////
                // R1 multiple mapped 
                ///////////////////////////////////// 
                else if (read1_maps.size() > 2 && read1_maps.size() % 2 == 0) { 
                    //std::cout << "R1 Multiple Mapped\n"; 
                    new_alignment.flag |= seqan::BAM_FLAG_SECONDARY; 
                    if (std::get<3>(read1_maps[1]) == "-") {
                        new_alignment.flag |= seqan::BAM_FLAG_RC; 
                    }
                    new_alignment.rID = name_to_id[std::get<0>(read1_maps[1])]; 
                    new_alignment.beginPos = std::get<1>(read1_maps[1]);
                    //new_alignment.beginPos = std::get<1>(read1_maps[1]);
                    new_alignment.mapQ = 255; 

                    // opposite strand 
                    if (std::get<3>(read1_maps[0]) != std::get<3>(read1_maps[1])) { 
//                        //std::cout << "Opposite strand\n"; 
						seqan::reverse(new_alignment.cigar); 
						seqan::Dna5String tmp = new_alignment.seq; 
						seqan::Dna5StringReverseComplement rc(tmp); 
						new_alignment.seq = seqan::CharString(rc); 

						tmp = new_alignment.qual; 
						seqan::Dna5StringReverse r(tmp); 
						new_alignment.qual = seqan::CharString(r); 
                    } 
                    
                    // R2 is unmapped 
                    if ((old_flag & seqan::BAM_FLAG_NEXT_UNMAPPED) || read2_maps.size() == 0) { 
                        new_alignment.flag |= seqan::BAM_FLAG_NEXT_UNMAPPED; 
                        new_alignment.rNextId = name_to_id[std::get<0>(read1_maps[1])]; 
                        new_alignment.pNext = std::get<1>(read1_maps[1]); 
                        new_alignment.tLen = 0; 
                        MN++; 
                        if (addtag) { 
							seqan::setTagValue(tagsDict, "MN", 0); 
                        }
                    } 
                    // R2 is unique mapped 
                    else if (read2_maps.size() == 2) { 
                        if (std::get<3>(read2_maps[1]) == "-") { 
                            new_alignment.flag |= seqan::BAM_FLAG_NEXT_RC; 
                        } 
                        new_alignment.rNextId = name_to_id[std::get<0>(read2_maps[1])]; 
                        new_alignment.pNext = std::get<1>(read2_maps[1]); 
                        new_alignment.tLen = 0; 
                        MU++; 
                        if (addtag) { 
							seqan::setTagValue(tagsDict, "MU", 0); 
                        } 
                    }
                    // R2 is multiple mapped. 
                    else {
                        if (std::get<3>(read2_maps[read2_maps.size()-1]) == "-") { 
                            new_alignment.flag |= seqan::BAM_FLAG_NEXT_RC; 
                        }
                        new_alignment.flag |= seqan::BAM_FLAG_SECONDARY; 
                        new_alignment.rNextId = name_to_id[std::get<0>(read2_maps[read2_maps.size()-1])]; 
                        new_alignment.pNext = std::get<1>(read2_maps[read2_maps.size()-1]); 
                        new_alignment.tLen = 0; 
                        MM++; 
                        if (addtag) { 
							seqan::setTagValue(tagsDict, "MM", 0); 
                        }
                    } 

                    seqan::writeRecord(bamFileOut, new_alignment); 
					
                    continue; 
                } 
            } 
        } 
        // single end sequencing 
        else {
            new_alignment.rNextId = -1;
            new_alignment.pNext = 0; 
            new_alignment.tLen = 0; 
        
            // originally unmapped 
            if (static_cast<bool>(old_flag & seqan::BAM_FLAG_UNMAPPED)) { 
                UNMAP << "chr\t" << old_reference_start << "\t" << old_reference_end << "\t" << unmap_queryname << "\n"; 
                continue; 
            }
            else {
                new_alignment.flag = 0x0; // clear flag 
                read1_chr = id_to_name[old_alignment.rID]; 

                if (static_cast<bool>(old_flag & seqan::BAM_FLAG_RC)) { 
                    read1_strand = "-"; 
                }
                else { 
                    read1_strand = "+"; 
                }
                read1_start = old_reference_start; 
                read1_end = old_reference_end;
                if (map_coordinates(mapping, read1_chr, read1_start, read1_end, read1_maps, read1_strand, false) == -1) {
                    //std::cout << "failed map_coordinates\n"; 
                }
                
                // unmapped after liftover 
                if (read1_maps.size() == 0) { 
                    UNMAP << read1_chr << "\t" << old_reference_start << "\t" << old_reference_end << "\t" << unmap_queryname << "\n"; 
                    continue; 
                } 

                // unique mapped 
                if (read1_maps.size() == 2) { 
                    if (std::get<3>(read1_maps[1]) == "-") { 
                        new_alignment.flag |= seqan::BAM_FLAG_RC; 
                    } 
                    if (std::get<3>(read1_maps[0]) != std::get<3>(read1_maps[1])) { 
						seqan::reverse(new_alignment.cigar); 
						seqan::Dna5String tmp = new_alignment.seq; 
						seqan::Dna5StringReverseComplement rc(tmp); 
						new_alignment.seq = seqan::CharString(rc); 

						tmp = new_alignment.qual; 
						seqan::Dna5StringReverse r(tmp); 
						new_alignment.qual = seqan::CharString(r); 
                    }

                    new_alignment.rID = name_to_id[std::get<0>(read1_maps[1])]; 
                    new_alignment.beginPos = std::get<1>(read1_maps[1]); 
                    SU++; 
                    if (addtag) { 
						seqan::setTagValue(tagsDict, "SU", 0); 
                    } 
                    seqan::writeRecord(bamFileOut, new_alignment); 
                    continue; 
                } 
                
                // multiple mapped 
                if (read1_maps.size() > 2 && read1_maps.size() % 2 == 0) { 
                    new_alignment.flag |= seqan::BAM_FLAG_SECONDARY; 
                    if (std::get<3>(read1_maps[1]) == "-") { 
                        new_alignment.flag |= seqan::BAM_FLAG_RC; 
                    } 
                    if (std::get<3>(read1_maps[0]) != std::get<3>(read1_maps[1])) { 
						seqan::reverse(new_alignment.cigar); 
						seqan::Dna5String tmp = new_alignment.seq; 
						seqan::Dna5StringReverseComplement rc(tmp); 
						new_alignment.seq = seqan::CharString(rc); 

						tmp = new_alignment.qual; 
						seqan::Dna5StringReverse r(tmp); 
						new_alignment.qual = seqan::CharString(r); 
                    } 

                    new_alignment.rID = name_to_id[std::get<0>(read1_maps[1])]; 
                    new_alignment.beginPos = std::get<1>(read1_maps[1]); 
                    SM++; 
                    if (addtag) {
						seqan::setTagValue(tagsDict, "SM", 0); 
                    } 
                    seqan::writeRecord(bamFileOut, new_alignment); 
                    continue; 
                } 
            }
        } 
 	} 
    
    UNMAP.close(); 

    // TODO: sort the entries.     

    std::cout << "Total alignments: " << total_item << "\n"; 
    std::cout << "       QC failed: " << QF << "\n"; 

    if (NN+NU+NM+UN+UU+UM+MN+MU+MM > 0) { 
        std::cout << "  Paired-end reads:\n";
        std::cout << "\tR1 unique, R2 unique (UU):     " << UU << "\n";  
        std::cout << "\tR1 unique, R2 unmap  (UN):     " << UN << "\n"; 
        std::cout << "\tR1 unique, R2 multiple (UM):   " << UM << "\n"; 
        
        std::cout << "\tR1 multiple, R2 unique (MU):   " << MU << "\n";  
        std::cout << "\tR1 multiple, R2 unmap  (MN):   " << MN << "\n"; 
        std::cout << "\tR1 multiple, R2 multiple (MM): " << MM << "\n"; 

        std::cout << "\tR1 unmap, R2 unique (NU):      " << NU << "\n";  
        std::cout << "\tR1 unmap, R2 unmap  (NN):      " << NN << "\n"; 
        std::cout << "\tR1 unmap, R2 multiple (NM):    " << NM << "\n"; 
    } 
    if (SN+SU+SM > 0) { 
        std::cout << "  Single-end reads:\n"; 
        std::cout << "\tUniquely mapped (SU): " << SU << "\n";  
        std::cout << "\tUnmapped (SN):        " << SN << "\n"; 
        std::cout << "\tMultiple mapped (SM): " << SM << "\n"; 
    }

    return 0; 
}
