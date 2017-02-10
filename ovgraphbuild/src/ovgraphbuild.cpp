/*
 * ovgraphbuild.cpp
 *
 * Created on: July 25, 2014
 * Modified on: February 24, 2016
 * Author: Pierre Pericard
 * Version: 1.0.0
 */

#define SEQAN_PROFILE // enable time measurements

#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <regex>
#include <typeinfo>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/map.h>

#include "optionsParser.h"
#include "assemblyStatistics.h"
#include "compatibilityGraphBuilding.h"
#include "compatibilityGraph.h"
#include "seqanObjectsFormating.h"

/***************************************************************************************
    Read the reference fasta file, store the sequences, the ids, and the lengths
***************************************************************************************/
void readRefFastaFile(seqan::StringSet<seqan::CharString> &refIds,
                      seqan::StringSet<seqan::Dna5String> &refSeqs,
                      seqan::Map<seqan::Pair<seqan::CharString, int> > &refLengthMap,
                      char* myRefFastaFile,
                      AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, myRefFastaFile))
    {
        std::cerr << "ERROR: Could not open " << myRefFastaFile << "!" << "\n";
        exit(1);
    }

    try
    {
        seqan::readRecords(refIds, refSeqs, seqFileIn);

        // This is the regex to retrieve the sequence identifier
        std::regex re("^([\\w\\.]+)\\s+.*");

        // Store reference seqs lengths in a map
        auto totalRefIdsNum = static_cast<int>(length(refIds));
        for (int i = 0; i < totalRefIdsNum; ++i)
        {
            std::cmatch match;

            try
            {
                std::regex_search(seqan::toCString(refIds[i]), match, re);
//                std::cout << seqan::toCString(refIds[i]) << "| |" << match.str(1)
//                          << "| |" << length(refSeqs[i]) << "\n";
                add(refLengthMap, match.str(1), length(refSeqs[i]));
            }
            catch (std::regex_error& e)
            {
                // Syntax error in the regular expression
                std::cerr << "ERROR: Pb with regex while parsing refId. " << e.what() << "\n";
            }
        }
    }
    catch (seqan::Exception const &e)
    {
        std::cerr << "ERROR: " << e.what() << "\n";
    }
}

/***************************************************************************************
    Read the reference name table from the input SAM file
***************************************************************************************/
void readRefNameFromSam(seqan::Map<seqan::Pair<seqan::CharString, int> > &refNameToIndexMap,
                        seqan::String<seqan::CharString> &refNameList,
                        char* mySamFile,
                        AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    std::ifstream inSamFile(mySamFile);
    std::string line;

    while(getline(inSamFile, line))
    {
        std::stringstream linestream(line);
        std::string readName;
        int tag;
        std::string refName;

        // read the spaces/tab separated line
        linestream >> readName >> tag >> refName;

        if (! hasKey(refNameToIndexMap, refName))
        {
            seqan::add(refNameToIndexMap, refName, -1);
            seqan::appendValue(refNameList, refName);
        }
    }

    // We MUST NOT sort the reference name list, to keep the same order
    // when reading bamRecord later (where rID are stored as int)

    // Store each reference name position in the list as a value
    //   in the map
    auto refSeqNum = static_cast<int>(length(refNameList));
    for (int i = 0; i < refSeqNum; ++i)
    {
//        std::cout << refNameList[i] << " index=" << i << "\n";
        refNameToIndexMap[refNameList[i]] = i;
    }

//    for (auto const &it : refNameList)
//    {
//        std::cout << it << "\t";
//    }
//    std::cout << "\n" << "\n";

//    for (auto const &it : refNameToIndexMap)
//    {
////        std::cout << seqan::getValueI1(it) << "\t";
//        std::cout << it << "\t";
//    }
//    std::cout << "\n" << "\n";
}

/***************************************************************************************
    Read the SAM file and store all BamRecords in a vector
***************************************************************************************/
void readSamFile(seqan::BamHeader &header,
                 std::vector<std::vector<seqan::BamAlignmentRecord> > &bamRecordBuffer,
                 GlobalStatistics &globalStats,
                 char* mySamFile,
                 AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    // Open input sam/bam file
    seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, mySamFile))
    {
        std::cerr << "ERROR: Could not open " << mySamFile << "!" << "\n";
        exit(1);
    }

    // Copy header
    seqan::readHeader(header, bamFileIn);

    try
    {
        // Read the bam records
        seqan::BamAlignmentRecord bamRecord;
        seqan::CharString oldQName = "";
        std::vector<seqan::BamAlignmentRecord> readBamRecordBuffer;

        while (!atEnd(bamFileIn))
        {
            seqan::readRecord(bamRecord, bamFileIn);
            ++globalStats.totalBamRecords;

//            if (options.debug)
//            {
//                std::cerr << formatBamRecordString(bamRecord) << "\n";
//            }

            if (bamRecord.qName != oldQName)
                ++globalStats.totalReadsNum;

            if (hasFlagUnmapped(bamRecord))
            {
                ++globalStats.unmappedAlignmentsNum;
            }
            else
            {
                if (bamRecord.qName != oldQName)
                {
                    ++globalStats.mappedReadsNum;

                    if (!readBamRecordBuffer.empty())
                        bamRecordBuffer.push_back(readBamRecordBuffer);

                    readBamRecordBuffer.clear();
                }

                readBamRecordBuffer.push_back(std::move(bamRecord));
            }

            oldQName = bamRecord.qName;
        }

        if (!readBamRecordBuffer.empty())
        {
            bamRecordBuffer.push_back(readBamRecordBuffer);

//            if (options.debug)
//            {
//                std::cerr << "\n";
//            }
        }

        readBamRecordBuffer.clear();
    }
    catch (seqan::Exception const &e)
    {
        std::cerr << "ERROR: " << e.what() << "\n";
    }
}

/******************************************************************************
    Print all compatibility stats following graph building
******************************************************************************/
void printCompatibilityStats(GlobalStatistics &globalStats)
{
    std::cout << "INFO: " << globalStats.readPairsWithCommonRefNum
              << " reads pairs have at least one common reference" << "\n"
              << "INFO: " << globalStats.numConsideredAlignmentsPairs
              << " alignments pairs were considered" << "\n"
              << "INFO: Of which " << globalStats.alignmentsPairsWithCommonRefNum
              << " pairs were mapped on the same reference"
              << "\n" << "\n";

    std::cout << std::setprecision(4)
              << "INFO: Number of potential overlaps: " << globalStats.numPotentialOverlaps << "\n"
              << "INFO: Number of true overlaps:      " << globalStats.numOverlaps
              << " (" << (globalStats.numOverlaps * 100.0 / globalStats.numPotentialOverlaps) << "% potOv)" << "\n"
//              << "INFO: Number of similar pairs:      " << globalStats.compatibleReadPairsNum
//              << " (" << (globalStats.compatibleReadPairsNum * 100.0 / globalStats.numPotentialOverlaps) << "% potOv)"
              << "INFO: Number of enclosed overlaps:  " << globalStats.enclosedOverlapsNum << "\n"
              << "\n";

    std::cout << "INFO: True Positives:  " << globalStats.truePositive << "\n"
              << "INFO: False Positives: " << globalStats.falsePositive << "\n"
              << "INFO: False Negatives: " << globalStats.falseNegative << "\n"
              << "INFO: True Negatives:  " << globalStats.trueNegative
              << "\n" << "\n";

    float sensitivity = 100.0 * globalStats.truePositive / (globalStats.truePositive + globalStats.falseNegative);
    float specificity = 100.0 * globalStats.trueNegative / (globalStats.falsePositive + globalStats.trueNegative);
    float precision = 100.0 * globalStats.truePositive / (globalStats.truePositive + globalStats.falsePositive);

    std::cout << std::setprecision(4)
              << "INFO: Sensitivity: " << sensitivity << "%" << "\n"
              << "INFO: Specificity: " << specificity << "%" << "\n"
              << "INFO: Precision  : " << precision << "%"
              << "\n" << "\n";
}

/******************************************************************************
    Main
******************************************************************************/
int main(int argc, char const ** argv)
{
    // Start a clock for the function
    clock_t begin_fct = clock();
    // Start an temp clock for intermediary blocks
    clock_t begin_tmp = clock();

    /// ////////////
    /// Get options

    AlphaOptions options;

    getOptions(options, argc, argv);

    GlobalStatistics globalStats;

    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    /// ///////////////////////////
    /// Read references fasta file

    seqan::StringSet<seqan::CharString> refIds;
    seqan::StringSet<seqan::Dna5String> refSeqs;
    seqan::Map<seqan::Pair<seqan::CharString, int> > refLengthMap;

    readRefFastaFile(refIds, refSeqs, refLengthMap, seqan::toCString(options.myRefFastaFile), options);

    if (options.verbose)
    {
        std::cout << "TIME: Reference fasta file read in " << double(clock() - begin_tmp) / CLOCKS_PER_SEC
                  << " seconds." << "\n";
        std::cout << "INFO: " << length(refIds) << " reference sequences were loaded" << "\n" << "\n";
    }

    // Reset tmp clock
    begin_tmp = clock();

    /// ////////////////////////////////////////////////
    /// Read SAM file to construct a reference name map

    seqan::Map<seqan::Pair<seqan::CharString, int> > refNameToIndexMap;
    seqan::String<seqan::CharString> refNameList;

    readRefNameFromSam(refNameToIndexMap, refNameList, seqan::toCString(options.mySamFile), options);

    int refSeqNum = length(refNameToIndexMap);

    if (options.verbose)
    {
        std::cout << "TIME: References names loaded from the SAM file in " << double(clock() - begin_tmp) / CLOCKS_PER_SEC
                  << " seconds." << "\n";
        std::cout << "INFO: " << refSeqNum << " references are present in the SAM file"
                  << "\n" << "\n";
    }

    // Reset tmp clock
    begin_tmp = clock();

    /// /////////////////
    /// SAM file reading

    seqan::BamHeader header;
    std::vector<std::vector<seqan::BamAlignmentRecord> > bamRecordBuffer;

    readSamFile(header, bamRecordBuffer, globalStats, seqan::toCString(options.mySamFile), options);

    auto mappedAlignmentsNum = globalStats.totalBamRecords - globalStats.unmappedAlignmentsNum;

    if (options.verbose)
    {
        std::cout << "TIME: SAM file reading finished in " << double(clock() - begin_tmp) / CLOCKS_PER_SEC
                  << " seconds." << "\n";
        std::cout << "INFO: " << globalStats.totalBamRecords << " bam records were read"
                  << ", representing " << globalStats.totalReadsNum << " reads"
                  << "\n"
                  << "INFO: " << mappedAlignmentsNum << " bam record were mapped on a reference"
                  << ", representing " << globalStats.mappedReadsNum << " mapped reads"
                  << "\n" << "\n";
    }

    // Reset the tmp clock
    begin_tmp = clock();

    /// ////////////////
    /// Graph building

    TGraph graph;
    std::vector<TVertexDescriptor > vertices;
    TProperties readNames;

    buildCompatibilityGraph(graph,
                            vertices,
                            readNames,
                            globalStats,
                            bamRecordBuffer,
                            options);

    if (options.verbose)
    {
        std::cout << "TIME: Compatibility graph was build in " << double(clock() - begin_tmp) / CLOCKS_PER_SEC
                  << " seconds." << "\n";
        std::cout << "INFO: " << globalStats.numConsideredReadsPairs << " reads pairs were considered" << "\n"
                  << "INFO: Of which " << globalStats.compatibleReadPairsNum
                  << " were found compatible" << "\n"
                  << "INFO:          " << globalStats.incompatibleReadPairsNum
                  << " were found incompatible" << "\n"
                  << "INFO: and      " << globalStats.neitherCompatNorIncompatReadPairsNum
                  << " were neither compatible nor incompatible"
                  << "\n" << "\n";
    }

    if (options.debug)
    {
        printCompatibilityStats(globalStats);
    }

    // Reset the tmp clock
    begin_tmp = clock();

    /// /////////////
    /// Main ending

    if (options.debug)
    {
        std::cout << "TIME: Function " << __FUNCTION__ << "() finished in "
                  << double(clock() - begin_fct) / CLOCKS_PER_SEC << " Seconds." << "\n" << "\n";
    }

    return 0;
} //~main()
