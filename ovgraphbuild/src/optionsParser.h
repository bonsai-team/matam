#ifndef GETOPTIONS_H
#define GETOPTIONS_H

#include <stdlib.h>
#include <stdint.h>

#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

/***************************************************************************************
    Options structure
    as seen in http://seqan.readthedocs.org/en/seqan-v2.0.0/Tutorial/ParsingCommandLineArguments.html
***************************************************************************************/
struct AlphaOptions
{
    seqan::CharString myRefFastaFile;
    seqan::CharString myRefPairwiseAlignFile;
    seqan::CharString mySamFile;
    seqan::CharString outputBasename;
    bool outputASQG;
    bool outputCSV;
    double minRefPairwisePercentId;
    int minOverlapLength;
    double idRateThreshold;
    int minNumTrailingMatches;
    bool verbose;
    bool debug;
    bool test;
    bool multiRef;
    bool noIndel;

    // Initialize the variables with initializer lists
    AlphaOptions() :
        myRefFastaFile(""), myRefPairwiseAlignFile(""), mySamFile(""),
        outputBasename(""), outputASQG(false), outputCSV(false),
        minRefPairwisePercentId(0),
        minOverlapLength(0), idRateThreshold(100.0), minNumTrailingMatches(0),
        verbose(false), debug(false), test(false), multiRef(false),
        noIndel(false)
    {}
};

/***************************************************************************************
    Parse the command line to retrieve the options
    as seen in http://seqan.readthedocs.org/en/seqan-v2.0.0/Tutorial/ParsingCommandLineArguments.html
***************************************************************************************/
auto parseCommandLine(AlphaOptions &, int, char const **);

/***************************************************************************************
    Print all the debug and verbose info at the start of the program, after
    arguments parsing
***************************************************************************************/
int printStartingDebugAndVerboseInfo(AlphaOptions &);

/***************************************************************************************
    Main procedure
***************************************************************************************/
int getOptions(AlphaOptions &, int, char const **);

#endif // GETOPTIONS_H
