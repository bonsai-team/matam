#include "alignmentsComparison.h"

template<typename T> class TD;

/******************************************************************************
    Reverse-complement an OverlapDescription
******************************************************************************/
void reverseComplementOverlapDescription(OverlapDescription &ovDescription)
{
    int32_t tmp = ovDescription.ovEndPosI;
    ovDescription.ovEndPosI = ovDescription.readLengthI - 1 - ovDescription.ovBeginPosI;
    ovDescription.ovBeginPosI = ovDescription.readLengthI - 1 - tmp;

    tmp = ovDescription.ovEndPosJ;
    ovDescription.ovEndPosJ = ovDescription.readLengthJ - 1 - ovDescription.ovBeginPosJ;
    ovDescription.ovBeginPosJ = ovDescription.readLengthJ - 1 - tmp;
}

/******************************************************************************
    Reverse Complement the read i overlap positions
******************************************************************************/
void reverseComplementReadI(OverlapDescription &ovDescription)
{
    int32_t tmp = ovDescription.ovEndPosI;
    ovDescription.ovEndPosI = ovDescription.readLengthI - 1 - ovDescription.ovBeginPosI;
    ovDescription.ovBeginPosI = ovDescription.readLengthI - 1 - tmp;
}

/******************************************************************************
    Reverse Complement the read j overlap positions
******************************************************************************/
void reverseComplementReadJ(OverlapDescription &ovDescription)
{
    int32_t tmp = ovDescription.ovEndPosJ;
    ovDescription.ovEndPosJ = ovDescription.readLengthJ - 1 - ovDescription.ovBeginPosJ;
    ovDescription.ovBeginPosJ = ovDescription.readLengthJ - 1 - tmp;
}

/******************************************************************************
    Compute the compatibility between 2 bam records
******************************************************************************/
void computeAlignmentsPairCompatibility(bool &areReadsCompatible,
                                        bool &areReadsIncompatible,
                                        bool &areReadsOverlapping,
                                        GlobalStatistics &globalStats,
                                        seqan::BamAlignmentRecord const &bamRecordI,
                                        seqan::BamAlignmentRecord const &bamRecordJ,
                                        AlphaOptions const &options)
{
    if (options.test)
    {
        printAlignmentsPairInfo(bamRecordI, bamRecordJ, std::cerr);
    }

    computeCompatOnSameRef(areReadsCompatible,
                           areReadsIncompatible,
                           areReadsOverlapping,
                           globalStats,
                           bamRecordI,
                           bamRecordJ,
                           options);
}

/******************************************************************************
    Compute the compatibility between 2 bam records mapped on the same ref
******************************************************************************/
void computeCompatOnSameRef(bool &areReadsCompatible,
                            bool &areReadsIncompatible,
                            bool &areReadsOverlapping,
                            GlobalStatistics &globalStats,
                            seqan::BamAlignmentRecord const &bamRecordI,
                            seqan::BamAlignmentRecord const &bamRecordJ,
                            AlphaOptions const &options)
{
    if (options.test)
    {
        std::cerr << "I'm in " << __FUNCTION__ << "\n";
    }

    ++globalStats.alignmentsPairsWithCommonRefNum;

    auto endPosI = bamRecordI.beginPos + seqan::getAlignmentLengthInRef(bamRecordI);
    auto endPosJ = bamRecordJ.beginPos + seqan::getAlignmentLengthInRef(bamRecordJ);

    auto overlapOnRefBeginPos = std::max(bamRecordI.beginPos, bamRecordJ.beginPos);
    auto overlapOnRefEndPos = std::min(endPosI, endPosJ);

    auto overlapOnRefSize = static_cast<int32_t>(overlapOnRefEndPos - overlapOnRefBeginPos);

    if (options.test)
    {
        std::cerr << overlapOnRefBeginPos << " " << overlapOnRefEndPos << " "
                  << overlapOnRefSize;
    }

    // First quick test, to determine if there is at least a one nucleotide
    // overlap between the reads on the reference
    if (overlapOnRefSize >= 1)
    {
        ++globalStats.numPotentialOverlaps;

        OverlapStatistics ovStats;

        computeOverlapOnSameRef(ovStats,
                                bamRecordI.seq,
                                bamRecordJ.seq,
                                bamRecordI.cigar,
                                bamRecordJ.cigar,
                                bamRecordI.beginPos,
                                bamRecordJ.beginPos,
                                overlapOnRefBeginPos,
                                overlapOnRefEndPos,
                                options);

        if (options.test)
        {
            printOverlapStats(ovStats);
        }

        // Test for real overlap size
        areReadsOverlapping = (ovStats.totalOverlapPositions >= options.minOverlapLength);

        if (areReadsOverlapping)
        {
            ++globalStats.numOverlaps;

            auto seqIdentityPercent = static_cast<double>(ovStats.matchesNum) / ovStats.totalOverlapPositions;

            if (options.test)
            {
                std::cout << seqIdentityPercent << " " << options.idRateThreshold << "\n";
            }

            if (seqIdentityPercent >= options.idRateThreshold)
            {
                areReadsCompatible = true;
//                printAlignmentsPairInfo(bamRecordI, bamRecordJ, std::cout);
            }
            else
            {
                areReadsIncompatible = true;
            }
        }
        else
        {
            areReadsIncompatible = true;
        }
    }
    else
    {
        areReadsIncompatible = true;
    }
}

/******************************************************************************
    Compute detailled stats on an overlap between two reads mapped on the
    same ref.
******************************************************************************/
void computeOverlapOnSameRef(OverlapStatistics &ovStats,
                             seqan::String<seqan::Iupac> const &seqI,
                             seqan::String<seqan::Iupac> const &seqJ,
                             seqan::String<seqan::CigarElement<> > const &cigarI,
                             seqan::String<seqan::CigarElement<> > const &cigarJ,
                             int32_t const beginPosI,
                             int32_t const beginPosJ,
                             int32_t const overlapOnRefBeginPos,
                             int32_t const overlapOnRefEndPos,
                             AlphaOptions const &options)
{
    if (options.test)
    {
        std::cerr << "I'm in " << __FUNCTION__ << "\n";
    }

    seqan::Iterator<const seqan::String<seqan::Iupac> >::Type seqItI = seqan::begin(seqI);
    seqan::Iterator<const seqan::String<seqan::Iupac> >::Type seqItJ = seqan::begin(seqJ);

    auto unfoldedCigarI = unfoldCigar(cigarI);
    auto unfoldedCigarJ = unfoldCigar(cigarJ);

    auto unfoldedCigarItI = std::begin(unfoldedCigarI);
    auto unfoldedCigarItJ = std::begin(unfoldedCigarJ);

    initializeUnfoldedCigarIt(unfoldedCigarItI, seqItI);
    initializeUnfoldedCigarIt(unfoldedCigarItJ, seqItJ);

    int32_t referencePosition;

    if (beginPosI <= beginPosJ)
    {
        referencePosition = beginPosI;

//        TD<decltype(unfoldedCigarItI)> cigarIt;

        if (options.test)
        {
            std::cerr << "simpleParcours i" << "\n";
        }

        simpleParcours(referencePosition,
                       seqItI,
                       unfoldedCigarItI,
                       overlapOnRefBeginPos);
    }
    else
    {
        referencePosition = beginPosJ;

        if (options.test)
        {
            std::cerr << "simpleParcours j" << "\n";
        }

        simpleParcours(referencePosition,
                       seqItJ,
                       unfoldedCigarItJ,
                       overlapOnRefBeginPos);
    }

    if (options.test)
    {
        std::cerr << "doubleParcours" << "\n";
    }

    doubleParcours(ovStats,
                   referencePosition,
                   seqItI,
                   seqItJ,
                   unfoldedCigarItI,
                   unfoldedCigarItJ,
                   overlapOnRefEndPos);

//    std::cerr << "I'm out " << __FUNCTION__ << "\n";
}

/******************************************************************************
    Unfold a cigar and return a string (ex: 3M --> MMM)
******************************************************************************/
std::string unfoldCigar(seqan::String<seqan::CigarElement<> > const &cigar)
{
    std::string unfoldedCigar;

    for (auto const &cigarElement : cigar)
    {
        unfoldedCigar.append(cigarElement.count, cigarElement.operation);
    }

    return unfoldedCigar;
}

/******************************************************************************

******************************************************************************/
void initializeUnfoldedCigarIt(std::string::iterator &unfoldedCigarIt,
                               seqan::Iterator<const seqan::String<seqan::Iupac> >::Type &seqIt)
{
    int softClippedNuclNum = 0;

    while (*unfoldedCigarIt == 'S')
    {
        ++unfoldedCigarIt;
        ++softClippedNuclNum;
    }

    seqIt += softClippedNuclNum;
}

/******************************************************************************

******************************************************************************/
int32_t simpleParcours(int32_t &refPosition,
//                       seqan::Iter<const seqan::String<seqan::SimpleType<unsigned char, seqan::Iupac_>, seqan::Alloc<> >, seqan::AdaptorIterator<const seqan::SimpleType<unsigned char, seqan::Iupac_>*, seqan::Tag<seqan::Default_> > > &seqIt,
                       seqan::Iterator<const seqan::String<seqan::Iupac> >::Type &seqIt,
                       std::string::iterator &unfoldedCigarIt,
                       int32_t const refPositionEnd)
{
//    std::cerr << "I'm in " << __FUNCTION__ << "\n";
//
//    std::cerr << refPosition << " " << refPositionEnd << "\n";

    int32_t readNuclNum = 0;

    while (refPosition < refPositionEnd)
    {
//        std::cerr << *seqIt << " " << *unfoldedCigarIt << " " << refPosition << "\n";

        if (*unfoldedCigarIt == 'M')
        {
            ++seqIt;
            ++readNuclNum;
            ++unfoldedCigarIt;
            ++refPosition;
        }
        else if (*unfoldedCigarIt == 'D')
        {
            ++unfoldedCigarIt;
            ++refPosition;
        }
        else if (*unfoldedCigarIt == 'I')
        {
            ++seqIt;
            ++readNuclNum;
            ++unfoldedCigarIt;
        }
    }

    return readNuclNum;
}

/******************************************************************************

******************************************************************************/
void doubleParcours(OverlapStatistics &ovStats,
                    int32_t &refPosition,
                    seqan::Iterator<const seqan::String<seqan::Iupac> >::Type &seqItI,
                    seqan::Iterator<const seqan::String<seqan::Iupac> >::Type &seqItJ,
                    std::string::iterator &unfoldedCigarItI,
                    std::string::iterator &unfoldedCigarItJ,
                    int32_t const refPositionEnd)
{
    while (refPosition < refPositionEnd)
    {
//        std::cerr << refPosition << " "
//                  << *seqItI << " " << *unfoldedCigarItI << " "
//                  << *seqItJ << " " << *unfoldedCigarItJ << " "
//                  << refPositionEnd << "\n";

        if (*unfoldedCigarItI == 'M')
        {
            ++ovStats.totalOverlapPositions;

            if (*unfoldedCigarItJ == 'M')
            {
                ++ovStats.coherentPositions;

                if (*seqItI == *seqItJ)
                {
                    //Match
                    ++ovStats.matchesNum;
                }
                else
                {
                    //Mismatch
                    ++ovStats.mismatchesNum;
                }

                ++seqItI;
                ++seqItJ;
                ++unfoldedCigarItI;
                ++unfoldedCigarItJ;
                ++refPosition;
            }
            else if (*unfoldedCigarItJ == 'D')
            {
                ++ovStats.indelNum;
                //
                ++seqItI;
                ++unfoldedCigarItI;
                ++unfoldedCigarItJ;
                ++refPosition;
            }
            else if (*unfoldedCigarItJ == 'I')
            {
                ++ovStats.indelNum;
                //
                ++seqItJ;
                ++unfoldedCigarItJ;
            }
        }
        else if (*unfoldedCigarItI == 'D')
        {
            if (*unfoldedCigarItJ == 'M')
            {
                ++ovStats.indelNum;
                ++ovStats.totalOverlapPositions;
                //
                ++seqItJ;
                ++unfoldedCigarItI;
                ++unfoldedCigarItJ;
                ++refPosition;
            }
            else if (*unfoldedCigarItJ == 'D')
            {
                //Nothing
                ++unfoldedCigarItI;
                ++unfoldedCigarItJ;
                ++refPosition;
            }
            else if (*unfoldedCigarItJ== 'I')
            {
                ++ovStats.indelNum;
                ++ovStats.totalOverlapPositions;
                //
                ++seqItJ;
                ++unfoldedCigarItJ;
            }
        }
        else if (*unfoldedCigarItI == 'I')
        {
            ++ovStats.totalOverlapPositions;

            if (*unfoldedCigarItJ == 'M' || *unfoldedCigarItJ == 'D')
            {
                ++ovStats.indelNum;
                //
                ++seqItI;
                ++unfoldedCigarItI;
            }
            else if (*unfoldedCigarItJ == 'I')
            {
                ++ovStats.coherentPositions;

                if (*seqItI == *seqItJ)
                {
                    //Match
                    ++ovStats.matchesNum;
                }
                else
                {
                    //Mismatch
                    ++ovStats.mismatchesNum;
                }

                ++seqItI;
                ++seqItJ;
                ++unfoldedCigarItI;
                ++unfoldedCigarItJ;
            }
        }
    }
}

/******************************************************************************
    Print all the info about a pair of alignments
******************************************************************************/
void printAlignmentsPairInfo(seqan::BamAlignmentRecord const &bamRecordI,
                             seqan::BamAlignmentRecord const &bamRecordJ,
                             std::ostream &stream)
{
    auto endPosI = bamRecordI.beginPos + seqan::getAlignmentLengthInRef(bamRecordI);
    auto endPosJ = bamRecordJ.beginPos + seqan::getAlignmentLengthInRef(bamRecordJ);

    stream << "###############################################" << "\n";

    stream << "read i: " << bamRecordI.qName << ", ref ID: " << bamRecordI.rID
           << ", beginPosI: " <<  bamRecordI.beginPos << ", endPosI: " << endPosI
           << ", revCompI: " << seqan::hasFlagRC(bamRecordI) << ", cigarI: ";

    for (auto const &cigarElement: bamRecordI.cigar)
    {
        stream << cigarElement.count << cigarElement.operation;
    }

    stream << ", cigarI: " << unfoldCigar(bamRecordI.cigar);

    stream << "\n";

    stream << "read j: " << bamRecordJ.qName << ", ref ID: " << bamRecordJ.rID
           << ", beginPosJ: " <<  bamRecordJ.beginPos << ", endPosJ: " << endPosJ
           << ", revCompJ: " << seqan::hasFlagRC(bamRecordJ) << ", cigarJ: ";

    for (auto const &cigarElement: bamRecordJ.cigar)
    {
        stream << cigarElement.count << cigarElement.operation;
    }

    stream << ", cigarJ: " << unfoldCigar(bamRecordJ.cigar) << "\n";
}
