#include "optionsParser.h"

/***************************************************************************************
    Parse the command line to retrieve the options
    as seen in http://seqan.readthedocs.org/en/seqan-v2.0.0/Tutorial/ParsingCommandLineArguments.html
***************************************************************************************/
seqan::ArgumentParser::ParseResult
parseCommandLine(AlphaOptions &options, int argc, char const **argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("ovgraphbuild");
    // Set short description, version, and date.
    seqan::setShortDescription(parser, "Build a read overlap graph from a mapping file");
    seqan::setVersion(parser, "1.0.0");
    seqan::setDate(parser, "February 2016");

    // Define Options -- Section Input Files
    seqan::addSection(parser, "Input");
    seqan::addOption(parser, seqan::ArgParseOption( "r", "reference", "Reference fasta file.",
        seqan::ArgParseArgument::INPUT_FILE, "IN_FILE"));
    seqan::addOption(parser, seqan::ArgParseOption( "p", "pairwise", "Reference pairwise alignment file, fasta format.",
        seqan::ArgParseArgument::INPUT_FILE, "IN_FILE"));
    seqan::addOption(parser, seqan::ArgParseOption( "s", "sam", "SAM/BAM input file.",
        seqan::ArgParseArgument::INPUT_FILE, "IN_FILE"));
    // Define Options -- Section Output Files
    seqan::addSection(parser, "Output");
    seqan::addOption(parser, seqan::ArgParseOption( "o", "output_basename", "Output files basename.",
        seqan::ArgParseArgument::STRING, "BASENAME"));
    seqan::addOption(parser, seqan::ArgParseOption( "", "asqg", "Output overlap graph in ASQG format."));
    seqan::addOption(parser, seqan::ArgParseOption( "", "csv", "Output overlap graph in 2 CSV files (nodes and edges)."));
    // Define Options -- Section Parameters
    seqan::addSection(parser, "Parameters");
    seqan::addOption(parser, seqan::ArgParseOption( "", "min_ref_pairwise_pid", "Minimum percent id for the ref pairwise alignments.",
        seqan::ArgParseArgument::DOUBLE, "MIN"));
    seqan::addOption(parser, seqan::ArgParseOption( "m", "min_overlap", "Minimum overlap size.",
        seqan::ArgParseArgument::INTEGER, "MIN"));
    seqan::addOption(parser, seqan::ArgParseOption( "i", "id_threshold", "Sequence identity to keep an overlap between [0.0:1.0].",
        seqan::ArgParseArgument::DOUBLE, "ID"));
    seqan::addOption(parser, seqan::ArgParseOption( "", "min_trail_matches", "Minimum trailing matches needed to anchor the pairwise overlaps (correction).",
        seqan::ArgParseArgument::INTEGER, "MIN"));
    seqan::addOption(parser, seqan::ArgParseOption( "", "no_indel", "Indels are not allowed in reads overlap."));
    seqan::addOption(parser, seqan::ArgParseOption( "", "illumina", "Optimize the parameters for illumina sequencing."));

    // Define Options -- Section Misc.
    seqan::addSection(parser, "Misc.");
    seqan::addOption(parser, seqan::ArgParseOption( "v", "verbose", "Turn on verbose output."));
    seqan::addOption(parser, seqan::ArgParseOption( "", "debug", "Turn on debug infos."));
    seqan::addOption(parser, seqan::ArgParseOption( "", "test", "Turn on test infos."));

    // Set required options
    setRequired(parser, "reference");
//    setRequired(parser, "pairwise");
    setRequired(parser, "sam");

    // Set format restrictions
    seqan::setValidValues(parser, "reference", "fasta fa");
    seqan::setValidValues(parser, "pairwise", "fasta fa");
    seqan::setValidValues(parser, "sam", "sam bam");

    // Set default values
    seqan::setDefaultValue(parser, "output_basename", "ovgraphbuild_output");

    seqan::setDefaultValue(parser, "min_ref_pairwise_pid", 90);
    seqan::setDefaultValue(parser, "min_overlap", 50);
    seqan::setDefaultValue(parser, "id_threshold", 1);
    seqan::setDefaultValue(parser, "min_trail_matches", 3);

    // Parse command line
    auto res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract and print the options.
    seqan::getOptionValue(options.myRefFastaFile, parser, "reference");
    seqan::getOptionValue(options.myRefPairwiseAlignFile, parser, "pairwise");
    seqan::getOptionValue(options.mySamFile, parser, "sam");

    seqan::getOptionValue(options.outputBasename, parser, "output_basename");
    options.outputASQG = seqan::isSet(parser, "asqg");
//    options.outputASQG = true;
    options.outputCSV = seqan::isSet(parser, "csv");
//    options.outputCSV = true;

    seqan::getOptionValue(options.minRefPairwisePercentId, parser, "min_ref_pairwise_pid");
    seqan::getOptionValue(options.minOverlapLength, parser, "min_overlap");
    seqan::getOptionValue(options.idRateThreshold, parser, "id_threshold");
    seqan::getOptionValue(options.minNumTrailingMatches, parser, "min_trail_matches");
    options.noIndel = seqan::isSet(parser, "no_indel");

    if (seqan::isSet(parser, "illumina"))
    {
        options.noIndel = true;
    }

    options.verbose = seqan::isSet(parser, "verbose");
    options.debug = seqan::isSet(parser, "debug");
    options.test = seqan::isSet(parser, "test");

    // TO DO: Add parameters check (e.g. 0 <= pid <= 100)

    return seqan::ArgumentParser::PARSE_OK;
}

/***************************************************************************************
    Print all the debug and verbose info at the start of the program, after
    arguments parsing
***************************************************************************************/
int printStartingDebugAndVerboseInfo(AlphaOptions &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()" << "\n" << "\n";

        std::cout << "DEBUG: size of bool      = " << sizeof(bool) << "\n"
                  << "DEBUG: size of char      = " << sizeof(char) << "\n"
                  << "DEBUG: size of short     = " << sizeof(short) << "\n"
                  << "DEBUG: size of int       = " << sizeof(int) << "\n"
                  << "DEBUG: size of int8_t    = " << sizeof(int8_t) << "\n"
                  << "DEBUG: size of int16_t   = " << sizeof(int16_t) << "\n"
                  << "DEBUG: size of int32_t   = " << sizeof(int32_t) << "\n"
                  << "DEBUG: size of int64_t   = " << sizeof(int64_t) << "\n"
                  << "DEBUG: size of long      = " << sizeof(long) << "\n"
                  << "DEBUG: size of long long = " << sizeof(long long) << "\n"
                  << "DEBUG: size of float     = " << sizeof(float) << "\n"
                  << "DEBUG: size of double    = " << sizeof(double) << "\n"
                  << "DEBUG: size of char*             = " << sizeof(char*) << "\n"
                  << "DEBUG: size of std::string       = " << sizeof(std::string) << "\n"
                  << "DEBUG: size of seqan::CharString = " << sizeof(seqan::CharString)
                  << "\n" << "\n";
    }

    if (options.verbose)
    {
        std::cout << "PARAM: References:         \t" << options.myRefFastaFile << "\n"
                  << "PARAM: Pairwise:           \t" << options.myRefPairwiseAlignFile << "\n"
                  << "PARAM: Sam file:           \t" << options.mySamFile << "\n"
                  << "PARAM: Output basename:    \t" << options.outputBasename << "\n"
                  << "PARAM: ASQG output:        \t" << options.outputASQG << "\n"
                  << "PARAM: CSV output:         \t" << options.outputCSV << "\n"
                  << "PARAM: Min Ref Pairw. pid: \t" << options.minRefPairwisePercentId << "\n"
                  << "PARAM: Min Overlap:        \t" << options.minOverlapLength << "\n"
                  << "PARAM: Id Threshold:       \t" << options.idRateThreshold << "\n"
                  << "PARAM: NoIndel:            \t" << options.noIndel << "\n"
                  << "PARAM: Debug:              \t" << options.debug << "\n"
                  << "PARAM: Verbose:            \t" << options.verbose << "\n"
                  << "PARAM: Test:               \t" << options.test << "\n" << "\n";
    }

    return 0;
}

/***************************************************************************************
    Main procedure
***************************************************************************************/
int getOptions(AlphaOptions &options, int argc, char const **argv)
{
    // Parse the command line.
    auto res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
    {
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }

    std::cout << "\n";

    printStartingDebugAndVerboseInfo(options);

    return 0;
}
