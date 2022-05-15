/** @file CommandLineParser.h
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/

#ifndef COMMAND_LINE_PARSER_H_
#define COMMAND_LINE_PARSER_H_

#include <stdio.h>
#include <iostream>
#include <vector>
#include <seqan/arg_parse.h>

/** @brief Struct for holding command line options and their values.
 *
 *  @see parseCommandOptions()
 *  @see parseInitialCommand()
 */
struct CommandLineParser{
    
    CommandLineParser():insert_size(200), insert_size_stdev(30), insert_size_fold(3), addtags(false){}

    seqan::CharString type;
    seqan::CharString chain_file;
    seqan::CharString in_file;
    seqan::CharString out_file;
    seqan::CharString unmapped_file;

    int insert_size;
    int insert_size_stdev;
    int insert_size_fold;
    bool addtags;
};

/** @brief Parse values in order to run either preprocessing step or correction step.
 *
 *  @param options Stores parsed values
 *  @param argc Number of arguments specified while running Apollo
 *  @param argv Argument values array
 *  @return seqan::ArgumentParser::PARSE_OK if everything went well
 */
seqan::ArgumentParser::ParseResult
parseCommandOptions(CommandLineParser& options, int argc, char const **argv){

    using namespace std;
    seqan::ArgumentParser parser("FastRemap: A Tool for Quickly Remapping Reads between Genome Assemblies");

    setVersion(parser, "1.0");
    setDate(parser, "May 2022");

    addOption(parser, seqan::ArgParseOption("f", "file-type", "'bam', 'sam', or 'bed' file depending on input file",
                                            seqan::ArgParseArgument::STRING, "STR", false));
    setValidValues(parser, "file-type", "bam sam bed");
    setRequired(parser, "file-type");

    addOption(parser, seqan::ArgParseOption("c", "chain-file", "A chain file (https://genome.ucsc.edu/goldenPath/help/chain.html) that describes regions of similarity between references.", seqan::ArgParseArgument::INPUT_FILE, "FILE", false));
    setRequired(parser, "chain-file");

    addOption(parser, seqan::ArgParseOption("i", "input", "{s,b}am or bed file containing elements to be remapped based on chain file", seqan::ArgParseArgument::INPUT_FILE, "FILE", false));
    setRequired(parser, "input");

    addOption(parser, seqan::ArgParseOption("o", "output", "File containing all the remapped elements from the input file", seqan::ArgParseArgument::OUTPUT_FILE, "FILE", false));
    setRequired(parser, "output");

    addOption(parser, seqan::ArgParseOption("u", "output-unmapped", "File containing all the elements that could not be remapped from the input file based on the provided chain file", seqan::ArgParseArgument::OUTPUT_FILE, "FILE", false));
    setRequired(parser, "output-unmapped");

    addOption(parser, seqan::ArgParseOption("a", "append-tags", "Add tag to each alignment in BAM file. Tags for pair-end alignments include: QF = QC failed, NN = both read1 and read2 unmapped, NU = read1 unmapped, read2 unique mapped, NM = read1 unmapped, multiple mapped, UN = read1 uniquely mapped, read2 unmap, UU = both read1 and read2 uniquely mapped, UM = read1 uniquely mapped, read2 multiple mapped, MN = read1 multiple mapped, read2 unmapped, MU = read1 multiple mapped, read2 unique mapped, MM = both read1 and read2 multiple mapped. Tags for single-end alignments include: QF = QC failed, SN = unmaped, SM = multiple mapped, SU = uniquely mapped."));

    addOption(parser, seqan::ArgParseOption("m", "mean", "Average insert size of pair-end sequencing (bp).", 
                                            seqan::ArgParseArgument::INTEGER, "INT", false));
    setDefaultValue(getOption(parser, "mean"), options.insert_size);
    seqan::setMinValue(parser, "mean", "0");

    addOption(parser, seqan::ArgParseOption("s", "stdev", "Stanadard deviation of insert size.", 
                                            seqan::ArgParseArgument::INTEGER, "INT", false));
    setDefaultValue(getOption(parser, "stdev"), options.insert_size_stdev);
    seqan::setMinValue(parser, "stdev", "0");

    addOption(parser, seqan::ArgParseOption("t", "times", "A mapped pair is considered as proper pair if both ends mapped to different strand and the distance between them is less then '-t' * stdev from the mean.", 
                                            seqan::ArgParseArgument::INTEGER, "INT", false));
    setDefaultValue(getOption(parser, "times"), options.insert_size_fold);
    seqan::setMinValue(parser, "times", "0");

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    if (res == seqan::ArgumentParser::PARSE_OK){

        getOptionValue(options.type, parser, "f");
        getOptionValue(options.chain_file, parser, "c");
        getOptionValue(options.in_file, parser, "i");
        getOptionValue(options.out_file, parser, "o");
        getOptionValue(options.unmapped_file, parser, "u");
        options.addtags = isSet(parser, "a");
        getOptionValue(options.insert_size, parser, "m");
        getOptionValue(options.insert_size_stdev, parser, "s");
        getOptionValue(options.insert_size_fold, parser, "t");
    }

    return res;
}

#endif
