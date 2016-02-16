#ifndef SEQANOBJECTSFORMATING_H
#define SEQANOBJECTSFORMATING_H

#include <stdlib.h>
#include <stdint.h>

#include <seqan/bam_io.h>

/***************************************************************************************
    Format a bamRecord and return it as a string
***************************************************************************************/
std::string formatBamRecordString(seqan::BamAlignmentRecord const &bamRecord);


#endif // SEQANOBJECTSFORMATING_H
