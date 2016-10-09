#ifndef COMPATIBILITYGRAPHBUILDING_H
#define COMPATIBILITYGRAPHBUILDING_H

#include <stdlib.h>
#include <stdint.h>
#include <ostream>
#include <algorithm>
#include <set>

#include <seqan/bam_io.h>
#include <seqan/align.h>

#include "compatibilityGraph.h"
#include "alignmentsComparison.h"
#include "optionsParser.h"
#include "assemblyStatistics.h"
#include "seqanObjectsFormating.h"

using TSequence = seqan::String<char>;                 // sequence type
using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;      // align type
//using TRow = seqan::Row<TAlign>::Type; // gappend sequence type

/******************************************************************************
    Build the compatibility graph
******************************************************************************/
void buildCompatibilityGraph(TGraph &,
                             std::vector<TVertexDescriptor > &,
                             TProperties &,
                             GlobalStatistics &,
                             std::vector<std::vector<seqan::BamAlignmentRecord> > const &,
                             AlphaOptions const &);

/******************************************************************************
    Initialize the graph vertices and store read names
******************************************************************************/
void initializeGraph(TGraph &,
                     std::vector<TVertexDescriptor > &,
                     TProperties &,
                     std::vector<std::vector<seqan::BamAlignmentRecord> > const &,
                     AlphaOptions const &);

/******************************************************************************
    Write all reads sequences to ASQG output file
******************************************************************************/
void writeReadsToASQG(std::ostream &,
                      std::vector<std::vector<seqan::BamAlignmentRecord> > const &,
                      AlphaOptions const &);

/******************************************************************************
    Write all reads sequences to CSV output file
******************************************************************************/
void writeReadsToCSV(std::ostream &,
                     TProperties const &,
                     std::vector<std::vector<seqan::BamAlignmentRecord> > const &,
                     AlphaOptions const &);

/******************************************************************************
    Compute the compatibility between 2 reads given all their bam records
******************************************************************************/
void computeReadsPairCompatibility(GlobalStatistics &,
                                   std::ofstream &,
                                   std::ofstream &,
                                   int64_t,
                                   int64_t,
                                   std::vector<seqan::BamAlignmentRecord> const &,
                                   std::vector<seqan::BamAlignmentRecord> const &,
                                   TProperties const &,
                                   AlphaOptions const &);

/******************************************************************************
    Write an overlap to the ASQG output file
******************************************************************************/
void writeOverlapToASQG(std::ostream &,
                        OverlapDescription const &);

/******************************************************************************
    Write an overlap to the CSV Edges output file
******************************************************************************/
void writeOverlapToCSV(std::ostream &,
                       OverlapDescription const &,
                       bool const,
                       bool const,
                       OverlapStatistics const &);

/******************************************************************************
    Align 2 read sequences
******************************************************************************/
int align2ReadSequences(TAlign &,
                        TSequence const &,
                        TSequence const &);

/******************************************************************************
    Compute statistics for a pairwise alignment from Seqan
******************************************************************************/
void computeAlignmentStats(OverlapStatistics &,
                           OverlapDescription &,
                           TAlign const &,
                           TSequence const &,
                           TSequence const &,
                           bool const);

/******************************************************************************
    Print the progress of the nested loop
******************************************************************************/
void printProgress(std::ostream &, int32_t, int64_t, int64_t);

#endif // COMPATIBILITYGRAPHBUILDING_H
