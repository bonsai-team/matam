#include "seqanObjectsFormating.h"

/***************************************************************************************
    Format a bamRecord and return it as a string
***************************************************************************************/
std::string formatBamRecordString(seqan::BamAlignmentRecord const &bamRecord)
{
    std::stringstream buffer;

    buffer << "qName=" << bamRecord.qName << "\t"
           << "rID=" << bamRecord.rID << "\t"
           << "beginPos=" << bamRecord.beginPos << "\t"
           << "cigar=";

    for (auto const &cigar : bamRecord.cigar)
        buffer << cigar.count << cigar.operation;

    buffer << "\t" << "seq=" << bamRecord.seq << "\t";

    return buffer.str();
}
