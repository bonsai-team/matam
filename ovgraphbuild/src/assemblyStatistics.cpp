#include "assemblyStatistics.h"

/******************************************************************************
    Add the stats of a local OverlapStatistics object to a global one
******************************************************************************/
void addOverlapStatsInto(OverlapStatistics &totalOvStats,
                         OverlapStatistics const &localOvStats)
{
    totalOvStats.totalOverlapPositions += localOvStats.totalOverlapPositions;
    totalOvStats.coherentPositions += localOvStats.coherentPositions;
    totalOvStats.indelNum += localOvStats.indelNum;
    totalOvStats.matchesNum += localOvStats.matchesNum;
    totalOvStats.mismatchesNum += localOvStats.mismatchesNum;
}

/******************************************************************************
    Print the info from an OverlapStatistics object
******************************************************************************/
void printOverlapStats(OverlapStatistics &ovStats)
{
    std::cout << "totalOverlapPositions = " << ovStats.totalOverlapPositions << "\n"
              << "coherentPositions = " << ovStats.coherentPositions << "\n"
              << "indelNum = " << ovStats.indelNum << "\n"
              << "matchesNum = " << ovStats.matchesNum << "\n"
              << "mismatchesNum = " << ovStats.mismatchesNum << "\n";
}
