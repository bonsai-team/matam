#include "TransitionMatrix.h"

TransitionMatrix::TransitionMatrix(seqan::Map<seqan::Pair<seqan::CharString, int> > &refNameToIndexMap,
                                   seqan::Map<seqan::Pair<seqan::CharString, int> > &refLengthMap,
                                   double pIdThreshold,
                                   bool verbose,
                                   bool debug):
    _refNameToIndexMap(refNameToIndexMap),
    _refLengthMap(refLengthMap),
    _pIdThreshold(pIdThreshold),
    _verbose(verbose),
    _debug(debug),
    // >18a|score=294,expect=1.8e-161,pid=84.45%|AXOR01000023.32966.34489: [206-1437]+
    _reRefId("^(\\d+)([ab])\\|(?:score=(\\d+),expect=([\\d\\.e\\-]+?),pid=([\\d\\.]+?)%)?\\|([\\w\\.]+):\\s+\\[(\\d+)-(\\d+)\\]([+-])\\s*", std::regex_constants::optimize),
    _referenceNum(length(refNameToIndexMap)),
    _data(((_referenceNum-1)*_referenceNum/2), nullptr),
    _dataExistVector(((_referenceNum-1)*_referenceNum/2), false)
{
    if (_debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }
}

TransitionMatrix::~TransitionMatrix()
{

}

int64_t TransitionMatrix::_coordinatesToIndice(int i, int j)
{
    return i*(_referenceNum-1) - i*(i+1)/2 + j-1;
}

void TransitionMatrix::_indexToCoordinates(int &i, int &j, int64_t k)
{
    i = std::floor(2*_referenceNum - 1 - std::sqrt((2*_referenceNum - 1)*(2*_referenceNum - 1) -8*k)/2);
    j = k - i*(_referenceNum-1) - i*(i-1)/2;
}

void TransitionMatrix::_setData(std::shared_ptr<std::vector<ConservedBloc> > p_transitionVector, int i, int j)
{
    _data[_coordinatesToIndice(i, j)] = p_transitionVector;
}

std::shared_ptr<std::vector<ConservedBloc> > TransitionMatrix::_getData(int i, int j)
{
    return _data[_coordinatesToIndice(i, j)];
}

/******************************************************************************
    Store all pairwise alignments from a buffer
******************************************************************************/
int TransitionMatrix::addPairwiseAlignments(seqan::StringSet<seqan::CharString> &pairwiseIds,
                                            seqan::StringSet<seqan::CharString> &pairwiseSeqs)
{
    int storedAlignmentNum = 0;

    auto itIds = seqan::begin(pairwiseIds);
    auto itIdsEnd = seqan::end(pairwiseIds);
    auto itSeqs = seqan::begin(pairwiseSeqs);

    while(itIds != itIdsEnd)
    {
        if (_storePairwiseAlignment(*itIds, *(itIds+1), *itSeqs, *(itSeqs+1)))
            ++storedAlignmentNum;

        itIds += 2;
        itSeqs += 2;
    }

    // BUG: storedAlignmentNum need to be computed again after sorting and filtering
    return storedAlignmentNum;
}

/******************************************************************************
    Given a pairwise alignment (fasta format), store it in the data vector
    if the 2 references are present in the reference map and if the alignment
    passes all filters. If information about A->B has already been stored,
    B->A info is not stored.
******************************************************************************/
bool TransitionMatrix::_storePairwiseAlignment(seqan::CharString const &idStringA,
                                               seqan::CharString const &idStringB,
                                               seqan::CharString &seqA,
                                               seqan::CharString &seqB)
{
    bool isAlignmentStored = false;

    FastaAlignmentHeader alignHeaderA, alignHeaderB;

    _parsePairwiseId(alignHeaderA, idStringA);
    _parsePairwiseId(alignHeaderB, idStringB);

    //Beware refNameToIndexMap does not contain all references, only the ones in the SAM file
    if (seqan::hasKey(_refNameToIndexMap, alignHeaderA.refId)
        && seqan::hasKey(_refNameToIndexMap, alignHeaderB.refId))
    {
        if (alignHeaderA.pairNum == alignHeaderB.pairNum
            && alignHeaderA.pairChar == 'a'
            && alignHeaderB.pairChar == 'b')
        {
            int indexRefA = _refNameToIndexMap[alignHeaderA.refId];
            int indexRefB = _refNameToIndexMap[alignHeaderB.refId];

            // To insure each ref pair is treated only once, and remove self matches
            //      Potential bug if matches are generated only once, then we cannot insure
            //      all pairs are treated.
            if (indexRefA < indexRefB)
            {
                // Filter based on the reference pairwise alignment pid threshold and
                // the orientation of the alignment
                if (alignHeaderA.percentId >= _pIdThreshold
                    && alignHeaderA.senseStrand)
                {
                    auto p_transitionVector = _getData(indexRefA, indexRefB);

                    // Create the transition vector on the Heap if it doesnt already exist
                    if (p_transitionVector == nullptr)
                    {
                        p_transitionVector = std::make_shared<std::vector<ConservedBloc> >();
                        _setData(p_transitionVector, indexRefA, indexRefB);

                        _dataExistVector[_coordinatesToIndice(indexRefA, indexRefB)] = true;
                    }

                    _fillTransitionVectors(p_transitionVector,
                                           indexRefA,
                                           indexRefB,
                                           seqA,
                                           seqB,
                                           alignHeaderA,
                                           alignHeaderB);

                    isAlignmentStored = true;
                }
            }
        }
        else
        {
            std::cerr << "ERROR: Pb with pairwise alignments parsing, please sort your file." << "\n";
        }
    }

    return isAlignmentStored;
}

/******************************************************************************
    Parse a pairwise alignment header
******************************************************************************/
void TransitionMatrix::_parsePairwiseId(FastaAlignmentHeader &alignHeader,
                                        seqan::CharString const &idString)
{
    try
    {
        std::cmatch match;
        std::regex_search(toCString(idString), match, _reRefId);

        alignHeader.pairNum = std::stoi(match.str(1));
        alignHeader.pairChar = match.str(2)[0];

        //BUG: score, evalue and percentId cannot be optional
        alignHeader.score = std::stoi(match.str(3));
        alignHeader.evalue = std::atof(match.str(4).c_str());
        alignHeader.percentId = std::atof(match.str(5).c_str());

        alignHeader.refId = match.str(6);

        alignHeader.startPos = std::stoi(match.str(7));
        alignHeader.stopPos = std::stoi(match.str(8));
        alignHeader.senseStrand = (match.str(9)[0] == '+');

//        if (_debug)
//        {
//            std::cout << "DEBUG: pairwiseId: " << idString << "\n"
//                      << "DEBUG: \t" << alignHeader.pairNum << " " << alignHeader.pairChar << " "
//                      << alignHeader.score << " " << alignHeader.evalue << " " << alignHeader.percentId << " "
//                      << alignHeader.refId << " " << alignHeader.startPos << " " << alignHeader.stopPos << " "
//                      << alignHeader.senseStrand << "\n";
//        }
    }
    catch (std::regex_error& e)
    {
        // Syntax error in the regular expression
        std::cerr << "ERROR: Pb with regex while parsing refId in TransitionMatrix. " << e.what() << "\n";
    }
}

/******************************************************************************
    Given a good pairwise alignment, convert it in conserved blocs and store
    them in the transition vector
******************************************************************************/
void TransitionMatrix::_fillTransitionVectors(std::shared_ptr<std::vector<ConservedBloc> > p_transitionVector,
                                              int const indexRefA,
                                              int const indexRefB,
                                              seqan::CharString const &seqA,
                                              seqan::CharString const &seqB,
                                              FastaAlignmentHeader const &alignHeaderA,
                                              FastaAlignmentHeader const &alignHeaderB)
{
//    if (_debug)
//    {
//        std::cout << alignHeaderA.refId << " index=" << indexRefA << "\n";
//        std::cout << alignHeaderB.refId << " index=" << indexRefB << "\n";
//
//        for (auto const &conservedBloc : *p_transitionVector)
//        {
//            std::cout << "{" << conservedBloc.startPosA << ", "
//                      << conservedBloc.startPosB << ", "
//                      << conservedBloc.length << "}"
//                      << "\n" << std::flush;
//        }
//
//        std::cout << alignHeaderA.pairNum << alignHeaderA.pairChar << "||"
//                  << alignHeaderA.refId << ": [" << alignHeaderA.startPos << "-"
//                  << alignHeaderA.stopPos << "]" << alignHeaderA.senseStrand << "\n"
//                  << seqA << "\n" << std::flush;
//        std::cout << alignHeaderB.pairNum << alignHeaderB.pairChar << "||"
//                  << alignHeaderB.refId << ": [" << alignHeaderB.startPos << "-"
//                  << alignHeaderB.stopPos << "]" << alignHeaderB.senseStrand << "\n"
//                  << seqB << "\n" << std::flush;
//    }

    int32_t seqAPos = alignHeaderA.startPos;
    int32_t seqBPos = alignHeaderB.startPos;

    int32_t blocStartA = -1;
    int32_t blocStartB = -1;
    int32_t blocLength = 0;

    bool isInABloc = false;

    auto seqAIt = seqan::begin(seqA);
    auto seqBIt = seqan::begin(seqB);

    for (; seqAIt != seqan::end(seqA); ++seqAIt, ++seqBIt)
    {
        if (*seqAIt == '-' || *seqBIt == '-')
        {
            if (isInABloc)
            {
                p_transitionVector->push_back(ConservedBloc(blocStartA, blocStartB, blocLength));
            }

            isInABloc = false;

            if (*seqAIt != '-')
            {
                ++seqAPos;
            }
            else if (*seqBIt != '-')
            {
                ++seqBPos;
            }
        }
        else
        {
            if (!isInABloc)
            {
                blocStartA = seqAPos;
                blocStartB = seqBPos;
                blocLength = 0;
            }

            isInABloc = true;
            ++blocLength;

            ++seqAPos;
            ++seqBPos;
        }
    }

    if (isInABloc)
    {
        p_transitionVector->push_back(ConservedBloc(blocStartA, blocStartB, blocLength));
    }

//    if (_debug)
//    {
//        for (auto const &conservedBloc : *p_transitionVector)
//        {
//            std::cout << "{" << conservedBloc.startPosA << ", "
//                      << conservedBloc.startPosB << ", "
//                      << conservedBloc.length << "}"
//                      << "\n" << std::flush;
//        }
//        std::cout << "\n" << std::flush;
//    }
}

/******************************************************************************
    Sort and filter out overlapping conserved blocs for each transition vector.
    When all transition vectors are filled, this function should be called
    to sort conserved blocs by position of ref A and remove chains of
    overlapping conserved blocs
******************************************************************************/
void TransitionMatrix::sortAndFilter()
{
    if (_debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    int64_t k = 0;

    for (auto &p_transitionVector : _data)
    {
        if (p_transitionVector != nullptr)
        {
            std::sort(std::begin(*p_transitionVector),
                      std::end(*p_transitionVector),
                      compareConservedBlocs);

            int indexRefA, indexRefB;
            _indexToCoordinates(indexRefA, indexRefB, k);

//            std::cout << "k=" << k << "\tA=" << indexRefA << "\tB=" << indexRefB << "\n";
//
//            for (auto const &conservedBloc : *p_transitionVector)
//            {
//                std::cout << "{" << conservedBloc.startPosA << ", "
//                          << conservedBloc.startPosB << ", "
//                          << conservedBloc.length << "}"
//                          << "\n" << std::flush;
//            }

            _filterOutOverlappingConservedBlocs(*p_transitionVector);

//            std::cout << "====" << "\n";
//
//            for (auto const &conservedBloc : *p_transitionVector)
//            {
//                std::cout << "{" << conservedBloc.startPosA << ", "
//                          << conservedBloc.startPosB << ", "
//                          << conservedBloc.length << "}"
//                          << "\n" << std::flush;
//            }
//
//            std::cout << "\n" << std::flush;

            if ((*p_transitionVector).size() == 0)
            {
                p_transitionVector = nullptr;
                _dataExistVector[k] = false;
            }
        }

        ++k;
    }
}

/******************************************************************************
    Comparison function between 2 conserved bloc (on their startPosA)
******************************************************************************/
bool compareConservedBlocs(ConservedBloc const &blocI, ConservedBloc const &blocJ)
{
    return (blocI.startPosA < blocJ.startPosA);
}

/******************************************************************************
    Given a sorted vector of ConservedBloc, remove all overlapping blocs
    The algorithm will identify chains of overlapping blocs (either on
    ref A or ref B) and will remove all blocs involved in each chain
******************************************************************************/
void TransitionMatrix::_filterOutOverlappingConservedBlocs(std::vector<ConservedBloc> &transitionVector)
{
    auto lastConservedBloc = ConservedBloc(-1, -1, 0);

    auto startOverlapIt = std::begin(transitionVector); // Will point to the beginning of a chain of overlapping blocs
    auto endOverlapIt = std::begin(transitionVector); // Will point to the end of a chain of overlapping blocs

    for (auto itCurrent = std::begin(transitionVector); itCurrent != std::end(transitionVector); ++itCurrent)
    {
        // Test to identify if the current bloc is overlapping with the previous one,
        // on either of each ref
        if ((*itCurrent).startPosA < (lastConservedBloc.startPosA + lastConservedBloc.length) ||
            (*itCurrent).startPosB < (lastConservedBloc.startPosB + lastConservedBloc.length))
        {
            // If the blocs are overlapping, the chain is extended
            endOverlapIt = itCurrent;
        }
        else
        {
            // If the chain length is > 1, then the blocs from the chain are overlapping
            // and should be removed from the vector
            if (endOverlapIt != startOverlapIt)
            {
                itCurrent = transitionVector.erase(startOverlapIt, endOverlapIt+1);
            }

            startOverlapIt = itCurrent;
            endOverlapIt = itCurrent;
        }

        // We store the bloc for the next turn
        lastConservedBloc = *itCurrent;
    }

    if (endOverlapIt != startOverlapIt)
    {
        transitionVector.erase(startOverlapIt, endOverlapIt+1);
    }
}

/******************************************************************************
    Get the number of references pairs with transition info
******************************************************************************/
int64_t TransitionMatrix::getFilledBoxesNum()
{
    int64_t filledBoxesCount = 0;

    for (auto const &p_transitionVector : _data)
    {
        if (p_transitionVector != nullptr) { ++filledBoxesCount; }
    }

    return filledBoxesCount;
}

/******************************************************************************
    Dump all transition vectors
******************************************************************************/
void TransitionMatrix::printData()
{
    int64_t k = 0;

    for (auto const &p_transitionVector : _data)
    {
        if (p_transitionVector != nullptr)
        {
            int indexRefA, indexRefB;
            _indexToCoordinates(indexRefA, indexRefB, k);

            std::cout << "k=" << k << "\tA=" << indexRefA << "\tB=" << indexRefB << "\n";

            for (auto const &conservedBloc : *p_transitionVector)
            {
                std::cout << "{" << conservedBloc.startPosA << ", "
                          << conservedBloc.startPosB << ", "
                          << conservedBloc.length << "}"
                          << "\n" << std::flush;
            }

            std::cout << "\n" << std::flush;
        }

        ++k;
    }
}

/******************************************************************************
    Return true if there is transition information between 2 given ref
******************************************************************************/
bool TransitionMatrix::isThereTransitionVector(int rIDI, int rIDJ)
{
    if (_dataExistVector[_coordinatesToIndice(rIDI, rIDJ)])
        return true;
    else
        return false;
}

///******************************************************************************
//    Return true if there is transition information between 2 given ref
//******************************************************************************/
//bool TransitionMatrix::isThereTransitionVector(int rIDI, int rIDJ)
//{
//    if (rIDI < rIDJ)
//    {
//        if (_dataExistVector[_coordinatesToIndice(rIDI, rIDJ)])
//            return true;
//    }
//    else
//    {
//        if (_dataExistVector[_coordinatesToIndice(rIDJ, rIDI)])
//            return true;
//    }
//
//    return false;
//}
