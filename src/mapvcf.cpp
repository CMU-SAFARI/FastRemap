#include "mapvcf.h" 
#include "common.h" 
#include "utils.h" 

#include<algorithm> 
#include<string> 
#include<stdio.h>
#include<iostream>
#include<sstream> 
#include<map> 
#include<vector> 
#include <set>
#include<filesystem> 
#include <regex>
#include <iterator>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;
void crossmap_vcf_file(std::map<std::string, ITree>& mapping, std::string infile, std::string outfile, std::string liftoverfile, std::string refgenome, bool noCompAllele = false, bool compress = false, std::string cstyle = "a"){
    /**
	Convert genome coordinates in VCF format.
	Parameters
	----------
	mapping : dict
		Dictionary with source chrom name as key, IntervalTree object as value.
	infile : file
		Input file in VCF format. Can be a regular or compressed (*.gz, *.Z,*.z, *.bz,
		*.bz2, *.bzip2) file, local file or URL (http://, https://, ftp://) pointing to
		remote file.
	outfile : str
		prefix of output files.
	liftoverfile : file
		Chain (https://genome.ucsc.edu/goldenPath/help/chain.html) format file. Can be a
		regular or compressed (*.gz, *.Z,*.z, *.bz, *.bz2, *.bzip2) file, local file or
		URL (http://, https://, ftp://) pointing to remote file.
	refgenome : file
		The genome sequence file of 'target' assembly in FASTA format.
	noCompAllele : bool
		A logical value indicates whether to compare ref_allele to alt_allele after
		liftover. If True, the variant will be marked as "unmap" if
		ref_allele == alt_allele.
	cstyle : str, optional
		Chromosome ID style. Must be one of ['a', 's', 'l'], where
		'a' : as-is. The chromosome ID of the output file is in the same style of the input file.
		's' : short ID, such as "1", "2", "X.
		'l' : long ID, such as "chr1", "chr2", "chrX.
	*/

    if(noCompAllele)
        std::cout << "Keep variants [reference_allele == alternative_allele] ..." << std::endl;
    else    
        std::cout << "Filter out variants [reference_allele == alternative_allele] ..." << std::endl;


    // TODO: convert these parts into seqan code
	// index refegenome file if it hasn't been done
    /**
	if not os.path.exists(refgenome + '.fai'):
		logging.info("Creating index for: %s" % refgenome)
		pysam.faidx(refgenome)
	if os.path.getmtime(refgenome + '.fai') < os.path.getmtime(refgenome):
		logging.info("Index file is older than reference genome. Re-creating index for: %s" % refgenome)
		pysam.faidx(refgenome)
    
	refFasta = pysam.Fastafile(refgenome)
    */

    std::ofstream FILE_OUT;
    FILE_OUT.open(outfile);
    std::ofstream UNMAP;
    UNMAP.open(outfile + ".unmap");

	int total = 0;
	int fail = 0;

    std::ifstream INFILE(infile); 
    std::string line; 
    std::vector<std::tuple<std::string, int, int, std::string>> matches;
    std::string chr_template;

    if (INFILE.is_open()) {    
        while (std::getline(INFILE, line)) {
            trim(line); 
            if (line.empty()) 
                continue; 

            // deal with meta-information lines.
            // meta-information lines needed in both mapped and unmapped files
            if (line.substr(0, 12) == "##fileformat"){
                FILE_OUT << line << "\n"; 
                UNMAP << line << "\n"; 
            }
            else if (line.substr(0, 6) == "##INFO"){
                FILE_OUT << line << "\n"; 
                UNMAP << line << "\n"; 
            }
            else if (line.substr(0, 8) == "##FILTER"){
                FILE_OUT << line << "\n"; 
                UNMAP << line << "\n"; 
            }
            else if (line.substr(0, 8) == "##FORMAT"){
                FILE_OUT << line << "\n"; 
                UNMAP << line << "\n"; 
            }
            else if (line.substr(0, 5) == "##ALT"){
                FILE_OUT << line << "\n"; 
                UNMAP << line << "\n"; 
            }
            else if (line.substr(0, 8) == "##SAMPLE"){
                FILE_OUT << line << "\n"; 
                UNMAP << line << "\n"; 
            }
            else if (line.substr(0, 10) == "##PEDIGREE"){
                FILE_OUT << line << "\n"; 
                UNMAP << line << "\n"; 
            }   


            // meta-information lines needed in unmapped files
            else if (line.substr(0, 10) == "##assembly"){
                UNMAP << line << "\n"; 
            }  
            else if (line.substr(0, 8) == "##contig"){
                UNMAP << line << "\n"; 
                std::size_t found = line.find("ID=chr");
                if (found != std::string::npos)
                    chr_template = "chr1";
                else
                    chr_template = "1";
            }  

            // update contig information
            else if (line.substr(0, 6) == "#CHROM"){
                /**
                 * This part is not important I think
                logging.info("Updating contig field ... ")
                target_gsize = dict(list(zip(refFasta.references, refFasta.lengths)))
                for chr_id in sorted(target_gsize):
                    print("##contig=<ID=%s,length=%d,assembly=%s>" % (update_chromID(chr_template, chr_id, cstyle), target_gsize[chr_id], os.path.basename(refgenome)), file=FILE_OUT)

                print("##liftOverProgram=<CrossMap,version=%s,website=https://sourceforge.net/projects/crossmap>" % __version__, file=FILE_OUT)
                print("##liftOverChainFile=<%s>" % liftoverfile, file=FILE_OUT)
                print("##originalFile=<%s>" % infile, file=FILE_OUT)
                print("##targetRefGenome=<%s>" % refgenome, file=FILE_OUT)
                print("##liftOverDate=<%s>" % datetime.date.today().strftime("%B%d,%Y"), file=FILE_OUT)
                */
                FILE_OUT << line << "\n"; 
                UNMAP << line << "\n"; 
                std::cout << "Lifting over ... "  << std::endl;
            }  
            else{
                if (line.substr(0, 1) == "#")
                    continue;
                
                std::vector<std::string> fields; 
                std::stringstream ss(line); 
                std::string token; 
                while(ss >> token) { 
                    fields.push_back(token); 
                }

                total = total + 1;

                std::string chrom = fields[0];
                int start = std::stoi(fields[1]) - 1;	 // 0 based
                int end = start + (fields[3]).length();

                map_coordinates(mapping, chrom, start, end, matches, "+", false);

                if(matches.size() == 0){
                    UNMAP << line << "\tFail(Unmap)\n"; 
                    fail = fail + 1;
                    continue;
                }

                if(matches.size() == 2){
                    // update chrom
                    std::string target_chr = std::get<0>(matches[1]);	// target_chr is from chain file, could be 'chr1' or '1'
                    int target_start = std::get<1>(matches[1]);
                    int target_end = std::get<2>(matches[1]);
                    fields[0] = target_chr;

                    // update start coordinate
                    fields[1] = std::to_string(target_start + 1);
                    
                    // TODO: convert these parts into seqan code
                    // update ref allele
                    /**
                    target_chr = update_chromID(refFasta.references[0], target_chr)
                    try:
                        fields[3] = refFasta.fetch(target_chr,target_start,target_end).upper()
                    except:
                        print (line + "\tFail(KeyError)", file=UNMAP)
                        fail += 1
                        continue
                    */

                    // update END if any
                    std::regex regex_end ("END=[0-9]+");
                    std::string replacement = "END=" + std::to_string(target_end);
                    fields[7] = std::regex_replace (fields[7], regex_end,replacement);

                    if(std::get<3>(matches[1]) == "-")
                        fields[4] = revcomp_DNA(fields[4], true);

                    // check if ref_allele is the same as alt_allele
                    if(noCompAllele){
                        for(int index = 0; index < fields.size(); index++){
                            FILE_OUT << fields[index] << "\t"; 
                        }
                        FILE_OUT << "\n"; 
                    }   
                    else{
                        if(fields[3] != fields[4]){
                            for(int index = 0; index < fields.size(); index++){
                                FILE_OUT << fields[index] << "\t"; 
                            }
                            FILE_OUT << "\n"; 
                        }
                        else{
                            UNMAP << line << "\tFail(REF==ALT)\n"; 
                            fail = fail + 1;
                        }
                    }
                }
                else{
                    UNMAP << line << "\tFail(Multiple_hits)\n"; 
                    fail = fail + 1;
                    continue;
                }
            }
        }
    }

	FILE_OUT.close();
	UNMAP.close();

    std::cout << "Total entries: " << total << std::endl;
    std::cout << "Failed to map: " << fail << std::endl;

	if(compress){
        /**
		try:
			logging.info("Compressing \"%s\" ..." % outfile)
			subprocess.call("gzip " + outfile, shell=True)
		except:
			pass
        */
    }

}