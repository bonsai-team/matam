#ifndef TRANSITIONMATRIX_H
#define TRANSITIONMATRIX_H

#include <stdlib.h>
#include <stdint.h>
#include <memory>
#include <vector>
#include <iostream>
#include <regex>
#include <string>
#include <cmath>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/map.h>

// Conserved blocs are assumed to represent alignments between
// the positive strand of 2 sequences. A bloc can contain only
// matches and/or mismatches, an indel leads to the creation
// of a new bloc.
struct ConservedBloc
{
    int32_t startPosA;
    int32_t startPosB;
    int32_t length;

    ConservedBloc():
        startPosA(-1),
        startPosB(-1),
        length(-1)
    {}

    ConservedBloc(int32_t startA,
                  int32_t startB,
                  int32_t lgth):
        startPosA(startA),
        startPosB(startB),
        length(lgth)
    {}
};

bool compareConservedBlocs(ConservedBloc const &, ConservedBloc const &);

struct FastaAlignmentHeader
{
    int pairNum;                // 4 byte
    int score;                  // 4 byte
    int startPos;               // 4 byte
    int stopPos;                // 4 byte
    double evalue;              // 8 byte
    double percentId;           // 8 byte
    seqan::CharString refId;    // 24 bytes
    bool senseStrand;           // 1 byte
    char pairChar;              // 1 byte

    FastaAlignmentHeader():
        pairNum(-1),
        score(0),
        startPos(-1),
        stopPos(-1),
        evalue(10e3),
        percentId(0),
        refId(""),
        senseStrand(true),
        pairChar('\0')
    {}
};

class TransitionMatrix
{
public:
    TransitionMatrix(seqan::Map<seqan::Pair<seqan::CharString, int> > &,
                     seqan::Map<seqan::Pair<seqan::CharString, int> > &,
                     double,
                     bool,
                     bool);
    virtual ~TransitionMatrix();
    int addPairwiseAlignments(seqan::StringSet<seqan::CharString> &,
                              seqan::StringSet<seqan::CharString> &);
    void sortAndFilter();
    int getCorrespondingPosition(int, seqan::CharString&, seqan::CharString &);
    bool isThereTransitionVector(int, int);
    int64_t getFilledBoxesNum();
    void printData();
    std::shared_ptr<std::vector<ConservedBloc> > _getData(int, int);
protected:
private:
    seqan::Map<seqan::Pair<seqan::CharString, int> > &_refNameToIndexMap;
    seqan::Map<seqan::Pair<seqan::CharString, int> > &_refLengthMap;
    double _pIdThreshold;
    bool _verbose, _debug;
    std::regex _reRefId;
    int _referenceNum;
    std::vector<std::shared_ptr<std::vector<ConservedBloc> > > _data;
    std::vector<bool> _dataExistVector;
    bool _storePairwiseAlignment(seqan::CharString const &,
                                 seqan::CharString const &,
                                 seqan::CharString &,
                                 seqan::CharString &);
    void _fillTransitionVectors(std::shared_ptr<std::vector<ConservedBloc> >,
                                int const,
                                int const,
                                seqan::CharString const &,
                                seqan::CharString const &,
                                FastaAlignmentHeader const &,
                                FastaAlignmentHeader const &);
    void _parsePairwiseId(FastaAlignmentHeader&, seqan::CharString const &);
    void _setData(std::shared_ptr<std::vector<ConservedBloc> >, int, int);
    void _filterOutOverlappingConservedBlocs(std::vector<ConservedBloc> &);
    int64_t _coordinatesToIndice(int, int);
    void _indexToCoordinates(int &, int &, int64_t);
};

#endif // TRANSITIONMATRIX_H
