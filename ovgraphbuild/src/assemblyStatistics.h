#ifndef ASSEMBLYSTATISTICS_H
#define ASSEMBLYSTATISTICS_H

#include <stdlib.h>
#include <stdint.h>
#include <iostream>

/***************************************************************************************
    Global Statistics structure
***************************************************************************************/
struct GlobalStatistics
{
    int64_t totalReadsNum;
    int64_t totalBamRecords;
    int64_t mappedReadsNum;
    int64_t unmappedAlignmentsNum;
    int64_t numConsideredReadsPairs;
    int64_t readPairsWithCommonRefNum;
    int64_t numConsideredAlignmentsPairs;
    int64_t alignmentsPairsWithCommonRefNum;
    int64_t numPotentialOverlaps;
    int64_t numOverlaps;
    int64_t enclosedOverlapsNum;
    int64_t compatibleReadPairsNum;
    int64_t incompatibleReadPairsNum;
    int64_t neitherCompatNorIncompatReadPairsNum;
    int64_t truePositive;
    int64_t falsePositive;
    int64_t trueNegative;
    int64_t falseNegative;

    // Initialize the variables with initializer lists
    GlobalStatistics() :
        totalReadsNum(0), totalBamRecords(0), mappedReadsNum(0),
        unmappedAlignmentsNum(0), numConsideredReadsPairs(0), readPairsWithCommonRefNum(0),
        numConsideredAlignmentsPairs(0), alignmentsPairsWithCommonRefNum(0), numPotentialOverlaps(0),
        numOverlaps(0), enclosedOverlapsNum(0),
        compatibleReadPairsNum(0), incompatibleReadPairsNum(0), neitherCompatNorIncompatReadPairsNum(0),
        truePositive(0), falsePositive(0), trueNegative(0), falseNegative(0)
    {}
};

/***************************************************************************************
    Overlap Statistics structure // ATTENTION: shoudl be in alignmentsComparison and should be merged with the OverlapDescription object
***************************************************************************************/
struct OverlapStatistics
{
    int32_t totalOverlapPositions;
    int32_t coherentPositions;
    int32_t indelNum;
    int32_t matchesNum;
    int32_t mismatchesNum;

    // Initialize the variables with initializer lists
    OverlapStatistics() :
        totalOverlapPositions(0), coherentPositions(0), indelNum(0),
        matchesNum(0), mismatchesNum(0)
    {}
};

/******************************************************************************
    Add the stats of a local OverlapStatistics object to a global one
******************************************************************************/
void addOverlapStatsInto(OverlapStatistics &,
                         OverlapStatistics const &);

/******************************************************************************
    Print the info from an OverlapStatistics object
******************************************************************************/
void printOverlapStats(OverlapStatistics &);

#endif // ASSEMBLYSTATISTICS_H
