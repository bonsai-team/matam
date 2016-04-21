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
                                        TransitionMatrix &transitionMatrix,
                                        AlphaOptions const &options)
{
//    TD<decltype(bamRecordI.seq)> bamRecordSeq;
//    TD<decltype(bamRecordI.cigar)> bamRecordCigar;

    if (options.test)
    {
        printAlignmentsPairInfo(bamRecordI, bamRecordJ, std::cerr);
    }

    if (bamRecordI.rID == bamRecordJ.rID)
    {
        computeCompatOnSameRef(areReadsCompatible,
                               areReadsIncompatible,
                               areReadsOverlapping,
                               globalStats,
                               bamRecordI,
                               bamRecordJ,
                               options);
    }
    // Only transition informations with referenceAId <= referenceBId were stored to prevent redundancy
    else if (bamRecordI.rID < bamRecordJ.rID)
    {
        computeCompatOnDifferentRefs(areReadsCompatible,
                                     areReadsIncompatible,
                                     areReadsOverlapping,
                                     globalStats,
                                     bamRecordI,
                                     bamRecordJ,
                                     transitionMatrix,
                                     options);
    }
    else
    {
        computeCompatOnDifferentRefs(areReadsCompatible,
                                     areReadsIncompatible,
                                     areReadsOverlapping,
                                     globalStats,
                                     bamRecordJ,
                                     bamRecordI,
                                     transitionMatrix,
                                     options);
    }
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
                   overlapOnRefEndPos,
                   options);

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
                    int32_t const refPositionEnd,
                    AlphaOptions const &options)
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
    Compute the compatibility between 2 bam records mapped on different refs
    Must be used with bamRecordI.rID < bamRecordJ.rID
******************************************************************************/
void computeCompatOnDifferentRefs(bool &areReadsCompatible,
                                  bool &areReadsIncompatible,
                                  bool &areReadsOverlapping,
                                  GlobalStatistics &globalStats,
                                  seqan::BamAlignmentRecord const &bamRecordI,
                                  seqan::BamAlignmentRecord const &bamRecordJ,
                                  TransitionMatrix &transitionMatrix,
                                  AlphaOptions const &options)
{
    if (options.test)
    {
        std::cerr << "I'm in " << __FUNCTION__ << "\n";
    }

    if (transitionMatrix.isThereTransitionVector(bamRecordI.rID, bamRecordJ.rID))
    {
        auto const p_transitionVector = transitionMatrix._getData(bamRecordI.rID, bamRecordJ.rID);

        auto endPosOnRefI = bamRecordI.beginPos + seqan::getAlignmentLengthInRef(bamRecordI);
        auto endPosOnRefJ = bamRecordJ.beginPos + seqan::getAlignmentLengthInRef(bamRecordJ);

        auto isOverlapPossible = true; // If the 2 reads doesnt share at least one common matching region, this will be set to false
        auto isOverlapDetermined = false; // If one read ends in a non-defined region and the other read begins in the same region, this will be set to false

        auto readsOverlapMinSize = static_cast<int32_t>(0);

        OverlapStatistics ovStats;

        computeBlocOverlapOnDifferentRefs(ovStats,
                                          readsOverlapMinSize,
                                          isOverlapPossible,
                                          isOverlapDetermined,
                                          bamRecordI,
                                          bamRecordJ,
                                          endPosOnRefI,
                                          endPosOnRefJ,
                                          *p_transitionVector,
                                          options);

        if (isOverlapPossible)
        {
            if (isOverlapDetermined)
            {
                ++globalStats.numPotentialOverlaps;

//                areReadsOverlapping = (readsOverlapMinSize >= options.minOverlapLength);
                areReadsOverlapping = (ovStats.totalOverlapPositions >= options.minOverlapLength);

                if (areReadsOverlapping)
                {
                    ++globalStats.numOverlaps;

                    auto seqIdentityPercent = static_cast<double>(ovStats.matchesNum) / ovStats.totalOverlapPositions;

                    if (seqIdentityPercent >= options.idRateThreshold)
                    {
                        areReadsCompatible = true;
//                        printAlignmentsPairInfo(bamRecordI, bamRecordJ, std::cout);
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
                // We cannot conclude so we go to the next alignments pair
                if (options.test)
                {
                    std::cout << "On ne peut pas conclure, donc on regarde la paire d'alignement suivant"
                              << "\n";
                }
            }
        }
        else
        {
            areReadsIncompatible = true;
        }
    }
    else
    {
        // No transition info between those 2 ref, so we go to the next alignments pair
    }
}

/******************************************************************************
    Compute transition blocs overlap between two reads mapped on
    different refs.
******************************************************************************/
void computeBlocOverlapOnDifferentRefs(OverlapStatistics &ovStats,
                                       int32_t &readsOverlapMinSize,
                                       bool &isOverlapPossible,
                                       bool &isOverlapDetermined,
                                       seqan::BamAlignmentRecord const &bamRecordI,
                                       seqan::BamAlignmentRecord const &bamRecordJ,
                                       int32_t const endPosOnRefI,
                                       int32_t const endPosOnRefJ,
                                       std::vector<ConservedBloc> const &transitionVector,
                                       AlphaOptions const &options)
{
    if (options.test)
    {
        std::cerr << "I'm in " << __FUNCTION__ << "\n";
    }

    int32_t blocNumBeginI = 0, blocNumEndI = 0, blocNumBeginJ = 0, blocNumEndJ = 0;

//    int i=0;

    auto transitionVectorIt = std::begin(transitionVector);
    auto const transitionVectorItEnd = std::end(transitionVector);

    while (transitionVectorIt != transitionVectorItEnd
           && ((endPosOnRefI >= transitionVectorIt->startPosA)
               && (endPosOnRefJ >= transitionVectorIt->startPosB)))
    {
        auto conservedBloc = *transitionVectorIt;

        int32_t blocEndPosA = conservedBloc.startPosA + conservedBloc.length;
        int32_t blocEndPosB = conservedBloc.startPosB + conservedBloc.length;

        if (bamRecordI.beginPos > blocEndPosA) blocNumBeginI += 2;
        else if (bamRecordI.beginPos >= conservedBloc.startPosA) ++blocNumBeginI;

        if (endPosOnRefI > blocEndPosA) blocNumEndI += 2;
        else if (endPosOnRefI >= conservedBloc.startPosA) ++blocNumEndI;

        if (bamRecordJ.beginPos > blocEndPosB) blocNumBeginJ += 2;
        else if (bamRecordJ.beginPos >= conservedBloc.startPosB) ++blocNumBeginJ;

        if (endPosOnRefJ > blocEndPosB) blocNumEndJ += 2;
        else if (endPosOnRefJ >= conservedBloc.startPosB) ++blocNumEndJ;

//        if (i>0)
//        {
//            printAlignmentsPairInfo(bamRecordI, bamRecordJ, std::cerr);
//
//            std::cerr << "\n" << "\n" << "\n"
//                      << conservedBloc.startPosA << " " << blocEndPosA << " "
//                      << conservedBloc.startPosB << " " << blocEndPosB << " "
//                      << blocNumBeginI << " " << blocNumEndI << " "
//                      << blocNumBeginJ << " " << blocNumEndJ << " "
//                      << "\n" << "\n" << "\n";
//        }

        ++transitionVectorIt;
//        ++i;
    }

    int32_t overlapBeginBlocNum = std::max(blocNumBeginI, blocNumBeginJ);
    int32_t overlapEndBlocNum = std::min(blocNumEndI, blocNumEndJ);
    int32_t overlapBlocSize = overlapEndBlocNum - overlapBeginBlocNum;

    if (overlapBlocSize >= 0)
    {
        isOverlapPossible = true;

        if (overlapBlocSize >= 1)
        {
            isOverlapDetermined = true;
        }
        else
        {
            if (overlapBeginBlocNum%2==1) isOverlapDetermined = true;
        }
    }

    if (isOverlapDetermined)
    {
        computeOverlapOnDifferentRefs(ovStats,
                                      readsOverlapMinSize,
                                      bamRecordI,
                                      bamRecordJ,
                                      endPosOnRefI,
                                      endPosOnRefJ,
                                      transitionVector,
                                      options);
    }
}

/******************************************************************************
    Compute detailled stats on an overlap between two reads mapped on
    different refs. Only to be called when the overlap is determined.
******************************************************************************/
void computeOverlapOnDifferentRefs(OverlapStatistics &ovStats,
                                   int32_t &readsOverlapMinSize,
                                   seqan::BamAlignmentRecord const &bamRecordI,
                                   seqan::BamAlignmentRecord const &bamRecordJ,
                                   int32_t const endPosOnRefI,
                                   int32_t const endPosOnRefJ,
                                   std::vector<ConservedBloc> const &transitionVector,
                                   AlphaOptions const &options)
{
    if (options.test)
    {
        std::cerr << "I'm in " << __FUNCTION__ << "\n";
    }

//    printAlignmentsPairInfo(bamRecordI, bamRecordJ, std::cerr);

    bool overlapExists = true;

    seqan::Iterator<const seqan::String<seqan::Iupac> >::Type seqItI = seqan::begin(bamRecordI.seq);
    seqan::Iterator<const seqan::String<seqan::Iupac> >::Type seqItJ = seqan::begin(bamRecordJ.seq);

    auto unfoldedCigarI = unfoldCigar(bamRecordI.cigar);
    auto unfoldedCigarJ = unfoldCigar(bamRecordJ.cigar);

    auto unfoldedCigarItI = std::begin(unfoldedCigarI);
    auto unfoldedCigarItJ = std::begin(unfoldedCigarJ);

    initializeUnfoldedCigarIt(unfoldedCigarItI, seqItI);
    initializeUnfoldedCigarIt(unfoldedCigarItJ, seqItJ);

    auto transitionVectorIt = std::begin(transitionVector);
    auto const transitionVectorItEnd = std::end(transitionVector);

    auto conservedBloc = *transitionVectorIt;

    int32_t blocEndPosA = conservedBloc.startPosA + conservedBloc.length;
    int32_t blocEndPosB = conservedBloc.startPosB + conservedBloc.length;

    // Go to the first bloc involved in the overlap
    while (bamRecordI.beginPos > blocEndPosA || bamRecordJ.beginPos > blocEndPosB)
    {
        if (options.test)
        {
            std::cerr << conservedBloc.startPosA << " " << blocEndPosA << " "
                      << conservedBloc.startPosB << " " << blocEndPosB << " "
                      << "\n";
        }

        ++transitionVectorIt;

        conservedBloc = *transitionVectorIt;
        blocEndPosA = conservedBloc.startPosA + conservedBloc.length;
        blocEndPosB = conservedBloc.startPosB + conservedBloc.length;
    }

    int32_t referencePositionA = bamRecordI.beginPos;
    int32_t referencePositionB = bamRecordJ.beginPos;

    if (options.test)
    {
        std::cerr << conservedBloc.startPosA << " " << blocEndPosA << " "
                  << conservedBloc.startPosB << " " << blocEndPosB << " "
                  << referencePositionA << " " << referencePositionB << " "
                  << "\n";
    }

    // Set all iterators to the start of the overlap between reads
    if (bamRecordI.beginPos < conservedBloc.startPosA
        && bamRecordJ.beginPos < conservedBloc.startPosB)
    {
        if (options.test)
        {
            std::cerr << "toto 1" << "\n";
        }

        simpleParcours(referencePositionA, seqItI, unfoldedCigarItI, conservedBloc.startPosA);
        simpleParcours(referencePositionB, seqItJ, unfoldedCigarItJ, conservedBloc.startPosB);
    }
    else if (bamRecordI.beginPos < conservedBloc.startPosA)
    {
        if (options.test)
        {
            std::cerr << "toto 2" << "\n";
        }

        int32_t beginPosJOnRefA = conservedBloc.startPosA + (bamRecordJ.beginPos - conservedBloc.startPosB);

        if (endPosOnRefI < beginPosJOnRefA)
        {
            overlapExists = false;
        }
        else
        {
            simpleParcours(referencePositionA, seqItI, unfoldedCigarItI, beginPosJOnRefA);
        }
    }
    else if (bamRecordJ.beginPos < conservedBloc.startPosB)
    {
        if (options.test)
        {
            std::cerr << "toto 3" << "\n";
        }

        int32_t beginPosIOnRefB = conservedBloc.startPosB + (bamRecordI.beginPos - conservedBloc.startPosA);

        if (endPosOnRefJ < beginPosIOnRefB)
        {
            overlapExists = false;
        }
        else
        {
            simpleParcours(referencePositionB, seqItJ, unfoldedCigarItJ, beginPosIOnRefB);
        }
    }
    else // the 2 alignements start in the bloc
    {
        int32_t beginPosIOnRefB = conservedBloc.startPosB + (bamRecordI.beginPos - conservedBloc.startPosA);

        if (bamRecordJ.beginPos < beginPosIOnRefB)
        {
            if (options.test)
            {
                std::cerr << "toto 4" << "\n";
            }

            if (endPosOnRefJ < beginPosIOnRefB)
            {
                overlapExists = false;
            }
            else
            {
                simpleParcours(referencePositionB, seqItJ, unfoldedCigarItJ, beginPosIOnRefB);
            }
        }
        else if (beginPosIOnRefB < bamRecordJ.beginPos)
        {
            if (options.test)
            {
                std::cerr << "toto 5" << "\n";
            }

            int32_t beginPosJOnRefA = conservedBloc.startPosA + (bamRecordJ.beginPos - conservedBloc.startPosB);

            if (endPosOnRefI < beginPosJOnRefA)
            {
                overlapExists = false;
            }
            else
            {
                simpleParcours(referencePositionA, seqItI, unfoldedCigarItI, beginPosJOnRefA);
            }
        }
        else
        {
            // bamRecordJ.beginPos == beginPosIOnRefB,
            // so all iterators are already synchronized

            if (options.test)
            {
                std::cerr << "toto 6" << "\n";
            }
        }
    }

//    std::cerr << conservedBloc.startPosA << " " << blocEndPosA << " "
//              << conservedBloc.startPosB << " " << blocEndPosB << " "
//              << referencePositionA << " " << endPosOnRefI << " "
//              << referencePositionB << " " << endPosOnRefJ << " "
//              << overlapExists << "\n";

    if (overlapExists)
    {
//        std::cerr << "Bloc suivant: " << (transitionVectorIt+1)->startPosA << " "
//                      << (transitionVectorIt+1)->startPosB << "\n";

        // Start overlap computing
        while (transitionVectorIt+1 != transitionVectorItEnd
               && endPosOnRefI >= (transitionVectorIt+1)->startPosA
               && endPosOnRefJ >= (transitionVectorIt+1)->startPosB)
        {
//            std::cerr << "Bloc suivant: " << (transitionVectorIt+1)->startPosA << " "
//                      << (transitionVectorIt+1)->startPosB << "\n";

//            std::cerr << "Nouveau bloc" << "\n";

            OverlapStatistics blocOvStats;

            int32_t oldReferencePositionA = referencePositionA;

            doubleParcours(blocOvStats,
                           referencePositionA,
                           seqItI,
                           seqItJ,
                           unfoldedCigarItI,
                           unfoldedCigarItJ,
                           blocEndPosA,
                           options);

//            std::cerr << "Out doubleParcours" << "\n";

            addOverlapStatsInto(ovStats, blocOvStats);

            referencePositionB += referencePositionA - oldReferencePositionA;

//            std::cerr << conservedBloc.startPosA << " " << blocEndPosA << " "
//                      << conservedBloc.startPosB << " " << blocEndPosB << " "
//                      << referencePositionA << " " << endPosOnRefI << " "
//                      << referencePositionB << " " << endPosOnRefJ << " "
//                      << ovStats.totalOverlapPositions << "\n";

            ++transitionVectorIt;

            conservedBloc = *transitionVectorIt;
            blocEndPosA = conservedBloc.startPosA + conservedBloc.length;
            blocEndPosB = conservedBloc.startPosB + conservedBloc.length;

//            std::cerr << conservedBloc.startPosA << " " << blocEndPosA << " "
//                      << conservedBloc.startPosB << " " << blocEndPosB << " "
//                      << referencePositionA << " " << endPosOnRefI << " "
//                      << referencePositionB << " " << endPosOnRefJ << " "
//                      << "\n";

            auto readNuclCrossedNumI = simpleParcours(referencePositionA,
                                                      seqItI,
                                                      unfoldedCigarItI,
                                                      conservedBloc.startPosA);

            auto readNuclCrossedNumJ = simpleParcours(referencePositionB,
                                                      seqItJ,
                                                      unfoldedCigarItJ,
                                                      conservedBloc.startPosB);

//            std::cerr << conservedBloc.startPosA << " " << blocEndPosA << " "
//                      << conservedBloc.startPosB << " " << blocEndPosB << " "
//                      << referencePositionA << " " << endPosOnRefI << " "
//                      << referencePositionB << " " << endPosOnRefJ << " "
//                      << "\n";
        }

        OverlapStatistics blocOvStats;

        int32_t oldReferencePositionA = referencePositionA;
        int32_t oldReferencePositionB = referencePositionB;

        // Last doubleParcours to reach the overlap end
//        std::cerr << "Commence à parcourir la fin de l'overlap" << "\n";

        if (endPosOnRefI > blocEndPosA && endPosOnRefJ > blocEndPosB)
        {
//            std::cerr << "tata 1" << "\n";

            doubleParcours(blocOvStats,
                           referencePositionA,
                           seqItI,
                           seqItJ,
                           unfoldedCigarItI,
                           unfoldedCigarItJ,
                           blocEndPosA,
                           options);

            referencePositionB += referencePositionA - oldReferencePositionA;
        }
        else if (endPosOnRefI > blocEndPosA)
        {
//            std::cerr << "tata 2" << "\n";

            doubleParcours(blocOvStats,
                           referencePositionB,
                           seqItI,
                           seqItJ,
                           unfoldedCigarItI,
                           unfoldedCigarItJ,
                           endPosOnRefJ,
                           options);

            referencePositionA += referencePositionB - oldReferencePositionB;
        }
        else if (endPosOnRefJ > blocEndPosB)
        {
//            std::cerr << "tata 3" << "\n";

            doubleParcours(blocOvStats,
                           referencePositionA,
                           seqItI,
                           seqItJ,
                           unfoldedCigarItI,
                           unfoldedCigarItJ,
                           endPosOnRefI,
                           options);

            referencePositionB += referencePositionA - oldReferencePositionA;
        }
        else
        {
            int32_t endPosIOnRefB = conservedBloc.startPosB + (endPosOnRefI - conservedBloc.startPosA);

            if (endPosIOnRefB <= endPosOnRefJ)
            {
//                std::cerr << "tata 4" << "\n";

                doubleParcours(blocOvStats,
                               referencePositionB,
                               seqItI,
                               seqItJ,
                               unfoldedCigarItI,
                               unfoldedCigarItJ,
                               endPosIOnRefB,
                               options);
            }
            else
            {
//                std::cerr << "tata 5" << "\n";

                doubleParcours(blocOvStats,
                               referencePositionB,
                               seqItI,
                               seqItJ,
                               unfoldedCigarItI,
                               unfoldedCigarItJ,
                               endPosOnRefJ,
                               options);
            }

            referencePositionA += referencePositionB - oldReferencePositionB;
        }

        addOverlapStatsInto(ovStats, blocOvStats);

//        std::cerr << conservedBloc.startPosA << " " << blocEndPosA << " "
//                  << conservedBloc.startPosB << " " << blocEndPosB << " "
//                  << referencePositionA << " " << referencePositionB << " "
//                  << ovStats.totalOverlapPositions << "\n";
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
