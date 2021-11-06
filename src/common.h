#ifndef __COMMON_H 
#define __COMMON_H 

#include "IntervalTree.h"

#include<string> 

typedef IntervalTree<int, std::tuple<std::string, int, int, std::string>> ITree;

// same format as seqan3::sam_record; 
typedef struct my_sam_record {
    std::string id;                 // string QNAME  COL 1    qName      id() 
    int flag;                       // int    FLAG   COL 2    flag       flag() 
    std::string reference_id;       // string RNAME  COL 3    rID        reference_id() 
    int reference_position;         // int    POS    COL 4    beginPos   reference_position()/sequence_position() 
    int mapping_quality;            // int    MAPQ   COL 5    mapQ       mapping_quality() 
    std::string cigar_sequence;     // string CIGAR  COL 6    cigar      cigar_sequence() 
    std::string mate_reference_id;  // string RNEXT  COL 7    rNextId    mate_reference_id() 
    std::string mate_position;      // string PNEXT  COL 8    pNext      mate_position() 
    std::string template_length;    // string TLEN   COL 9    tLen       template_length() 
    std::string sequence;           // string SEQ    COL 10   seq        sequence() 
    std::string base_qualities;     // string QUAL   COL 11   qual       base_qualities() 
    std::string tags;               // string tags   COL 12   tags       tags() 
} mysam_record; 

#endif 
