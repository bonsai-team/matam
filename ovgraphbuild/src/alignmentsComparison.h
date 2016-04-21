#ifndef ALIGNMENTSCOMPARISON_H
#define ALIGNMENTSCOMPARISON_H

#include <stdlib.h>
#include <stdint.h>

#include <seqan/bam_io.h>

#include "assemblyStatistics.h"
#include "TransitionMatrix.h"
#include "optionsParser.h"

/***************************************************************************************
    Overlap Description structure
    Positions of both reads are given on the sense strand, as in the ASQG format
***************************************************************************************/
struct OverlapDescription
{
    int64_t readIDI;
    int64_t readIDJ;
    int32_t readLengthI;
    int32_t readLengthJ;
    int32_t ovBeginPosI;
    int32_t ovEndPosI;
    int32_t ovBeginPosJ;
    int32_t ovEndPosJ;
    bool isEnclosedOverlap;
    bool isRead2RC;

    // Initialize the variables with initializer lists
    OverlapDescription() :
        readIDI(0), readIDJ(0),
        readLengthI(0), readLengthJ(0),
        ovBeginPosI(0), ovEndPosI(0),
        ovBeginPosJ(0), ovEndPosJ(0),
        isEnclosedOverlap(false),
        isRead2RC(false)
    {}
};

/******************************************************************************
    Reverse Complement an OverlapDescription
******************************************************************************/
void reverseComplementOverlapDescription(OverlapDescription &);

/******************************************************************************
    Reverse Complement the read i overlap positions
******************************************************************************/
void reverseComplementReadI(OverlapDescription &);

/******************************************************************************
    Reverse Complement the read j overlap positions
******************************************************************************/
void reverseComplementReadJ(OverlapDescription &);

/******************************************************************************
    Compute the compatibility between 2 bam records
******************************************************************************/
void computeAlignmentsPairCompatibility(bool &,
                                        bool &,
                                        bool &,
                                        GlobalStatistics &,
                                        seqan::BamAlignmentRecord const &,
                                        seqan::BamAlignmentRecord const &,
                                        TransitionMatrix &,
                                        AlphaOptions const &);

/******************************************************************************
    Compute the compatibility between 2 bam records mapped on the same ref
******************************************************************************/
void computeCompatOnSameRef(bool &,
                            bool &,
                            bool &,
                            GlobalStatistics &,
                            seqan::BamAlignmentRecord const &,
                            seqan::BamAlignmentRecord const &,
                            AlphaOptions const &);

/******************************************************************************
    Compute detailled stats on an overlap between two reads mapped on the
    same ref.
******************************************************************************/
void computeOverlapOnSameRef(OverlapStatistics &,
                             seqan::String<seqan::Iupac> const &,
                             seqan::String<seqan::Iupac> const &,
                             seqan::String<seqan::CigarElement<> > const &,
                             seqan::String<seqan::CigarElement<> > const &,
                             int32_t const, // beginPosI
                             int32_t const, // beginPosJ
                             int32_t const,
                             int32_t const,
                             AlphaOptions const &);

/******************************************************************************
    Unfold a cigar and return a string (ex: 3M --> MMM)
******************************************************************************/
std::string unfoldCigar(seqan::String<seqan::CigarElement<> > const &);

/******************************************************************************

******************************************************************************/
void initializeUnfoldedCigarIt(std::string::iterator &,
                               seqan::Iterator<const seqan::String<seqan::Iupac> >::Type &);

/******************************************************************************

******************************************************************************/
int32_t simpleParcours(int32_t &,
//                       seqan::Iter<const seqan::String<seqan::SimpleType<unsigned char, seqan::Iupac_>, seqan::Alloc<> >, seqan::AdaptorIterator<const seqan::SimpleType<unsigned char, seqan::Iupac_>*, seqan::Tag<seqan::Default_> > > &,
                       seqan::Iterator<const seqan::String<seqan::Iupac> >::Type &,
                       std::string::iterator &,
                       int32_t const);

/******************************************************************************

******************************************************************************/
void doubleParcours(OverlapStatistics &,
                    int32_t &,
                    seqan::Iterator<const seqan::String<seqan::Iupac> >::Type &,
                    seqan::Iterator<const seqan::String<seqan::Iupac> >::Type &,
                    std::string::iterator &,
                    std::string::iterator &,
                    int32_t const,
                    AlphaOptions const &);

/******************************************************************************
    Compute the compatibility between 2 bam records mapped on different refs
    Must be used with bamRecordI.rID < bamRecordJ.rID
******************************************************************************/
void computeCompatOnDifferentRefs(bool &,
                                  bool &,
                                  bool &,
                                  GlobalStatistics &,
                                  seqan::BamAlignmentRecord const &, //
                                  seqan::BamAlignmentRecord const &,
                                  TransitionMatrix &,
                                  AlphaOptions const &);

/******************************************************************************
    Compute transition blocs overlap between two reads mapped on
    different refs.
******************************************************************************/
void computeBlocOverlapOnDifferentRefs(OverlapStatistics &,
                                       int32_t &,
                                       bool &,
                                       bool &,
                                       seqan::BamAlignmentRecord const &,
                                       seqan::BamAlignmentRecord const &,
                                       int32_t const,
                                       int32_t const,
                                       std::vector<ConservedBloc> const &,
                                       AlphaOptions const &);

/******************************************************************************
    Compute detailled stats on an overlap between two reads mapped on
    different refs. nly to be called when the overlap is determined.
******************************************************************************/
void computeOverlapOnDifferentRefs(OverlapStatistics &,
                                   int32_t &,
                                   seqan::BamAlignmentRecord const &,
                                   seqan::BamAlignmentRecord const &,
                                   int32_t const,
                                   int32_t const,
                                   std::vector<ConservedBloc> const &,
                                   AlphaOptions const &);

/******************************************************************************
    Print all the info about a pair of alignments
******************************************************************************/
void printAlignmentsPairInfo(seqan::BamAlignmentRecord const &,
                             seqan::BamAlignmentRecord const &,
                             std::ostream &);

#endif // ALIGNMENTSCOMPARISON_H
