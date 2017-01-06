#include "compatibilityGraphBuilding.h"

//template<typename T> class TD;

/******************************************************************************
    Build the compatibility graph
******************************************************************************/
void buildCompatibilityGraph(TGraph &graph,
                             std::vector<TVertexDescriptor > &vertices,
                             TProperties &readNames,
                             GlobalStatistics &globalStats,
                             std::vector<std::vector<seqan::BamAlignmentRecord> > const &bamRecordBuffer,
                             AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    // Initialize graph vertices and store read names
    initializeGraph(graph, vertices, readNames, bamRecordBuffer, options);

    // Compute the compatibility for each read pair
    int64_t mappedReadsNum = bamRecordBuffer.size();
    int64_t maxReadsPairs = mappedReadsNum * (mappedReadsNum - 1) / 2;
    int64_t numConsideredReadsPairs = 0;

    if (options.debug)
    {
        std::cout << "DEBUG: Computing compatibility graph"
                  << "\n" << "\n" << std::flush;
    }

    // for performance reason to improve stream speed
    std::ios::sync_with_stdio(false);

    // Declare output files for the overlap graph
    std::ofstream asqgFile;
    std::ofstream csvNodesFile;
    std::ofstream csvEdgesFile;

    // Initialise the ASQG output file if needed
    if (options.outputASQG)
    {
        std::string asqgFilename(seqan::toCString(options.outputBasename));
        asqgFilename += ".asqg";

        // Opening ASQG output file
        asqgFile.open(asqgFilename, std::ofstream::out | std::ofstream::trunc);

        // Write ASQG header
        asqgFile << "HT\tVN:i:1\tER:f:0\t"
                 << "OL:i:" << options.minOverlapLength << "\t"
                 << "IN:Z:" << options.mySamFile << "\t"
                 << "CN:i:0\tTE:i:0\n";

        // Write ASQG nodes (aka reads)
        writeReadsToASQG(asqgFile, bamRecordBuffer, options);
    }

    // Initialise the CSV output files if needed
    if (options.outputCSV)
    {
        std::string csvNodesFilename(seqan::toCString(options.outputBasename));
        csvNodesFilename += ".nodes.csv";
        std::string csvEdgesFilename(seqan::toCString(options.outputBasename));
        csvEdgesFilename += ".edges.csv";

        // Opening CSV output files
        csvNodesFile.open(csvNodesFilename, std::ofstream::out | std::ofstream::trunc);
        csvEdgesFile.open(csvEdgesFilename, std::ofstream::out | std::ofstream::trunc);

        // Write CSV nodes
        writeReadsToCSV(csvNodesFile, readNames, bamRecordBuffer, options);
        csvNodesFile.close();

        // Write CSV Edges header
        csvEdgesFile << "Source;Target;Type;Coemitted;Weight;Pid;MultiRef" << "\n";
    }

    // Start the nested loop
    auto bamRecordBufferItI = std::begin(bamRecordBuffer);
    auto const bamRecordBufferEndItI = std::end(bamRecordBuffer)-1;
    auto bamRecordBufferItJ = bamRecordBufferItI+1;
    auto const bamRecordBufferEndItJ = std::end(bamRecordBuffer);

    int64_t i=0, j=1;

    for (; bamRecordBufferItI != bamRecordBufferEndItI; ++bamRecordBufferItI)
    {
//        std::cerr << i << "\n" << std::flush;

        for (bamRecordBufferItJ = bamRecordBufferItI+1;
             bamRecordBufferItJ != bamRecordBufferEndItJ;
             ++bamRecordBufferItJ)
        {
//            std::cerr << j << "\n" << std::flush;

            ++numConsideredReadsPairs;

            computeReadsPairCompatibility(globalStats,
                                          asqgFile,
                                          csvEdgesFile,
                                          i,
                                          j,
                                          *bamRecordBufferItI,
                                          *bamRecordBufferItJ,
                                          readNames,
                                          options);

            if (options.verbose)
            {
                if (numConsideredReadsPairs <= (float)maxReadsPairs / 100.0)
                    printProgress(std::cout, maxReadsPairs/100000, numConsideredReadsPairs, maxReadsPairs);
                else
                    printProgress(std::cout, maxReadsPairs/1000, numConsideredReadsPairs, maxReadsPairs);
            }

            ++j;
        }

        ++i;
        j=i+1;
    }

    if (options.verbose)
    {
        printProgress(std::cout, 1, numConsideredReadsPairs, maxReadsPairs);
        std::cout << "\n" << "\n";
    }

    asqgFile.close();
    csvEdgesFile.close();

    globalStats.numConsideredReadsPairs += numConsideredReadsPairs;
}

/******************************************************************************
    Initialize the graph vertices and store read names
******************************************************************************/
void initializeGraph(TGraph &graph,
                     std::vector<TVertexDescriptor > &vertices,
                     TProperties &readNames,
                     std::vector<std::vector<seqan::BamAlignmentRecord> > const &bamRecordBuffer,
                     AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    seqan::CharString oldQName = "";

    // Create a vertex for each read
    for (auto i=static_cast<unsigned>(0); i<bamRecordBuffer.size(); ++i)
    {
        vertices.push_back(addVertex(graph));
    }

    // Resize the property map containing the read names to fit all reads from the graph
    seqan::resizeVertexMap(readNames, graph);

    // Assign each read name to the corresponding vertex
    for (auto i=static_cast<unsigned>(0); i<bamRecordBuffer.size(); ++i)
    {
        seqan::assignProperty(readNames, vertices[i], bamRecordBuffer[i][0].qName);
    }
}

/******************************************************************************
    Write all reads sequences to ASQG output file
******************************************************************************/
void writeReadsToASQG(std::ostream &asqgFile,
                      std::vector<std::vector<seqan::BamAlignmentRecord> > const &bamRecordBuffer,
                      AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    int64_t i = 0;

    for (auto const &readBamRecordBuffer : bamRecordBuffer)
    {
        auto &readBamRecord = readBamRecordBuffer[0];

        seqan::CharString readSeq(readBamRecord.seq);

        if (seqan::hasFlagRC(readBamRecord))
        {
            seqan::reverseComplement(readSeq);
        }

        asqgFile << "VT\t" << i << "\t"
                 << readSeq << "\n";

        ++i;
    }
}

/******************************************************************************
    Write all reads sequences to CSV output file
******************************************************************************/
void writeReadsToCSV(std::ostream &csvNodesFile,
                     TProperties const &readNames,
                     std::vector<std::vector<seqan::BamAlignmentRecord> > const &bamRecordBuffer,
                     AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    csvNodesFile << "Id;Label;Specie" << "\n";

    int64_t i = 0;

    for (auto const &readBamRecordBuffer : bamRecordBuffer)
    {
        auto &readBamRecord = readBamRecordBuffer[0];

        auto specieID = seqan::prefix(readNames[i], 3);

        csvNodesFile << i << ";" << readBamRecord.qName << ";" << specieID << "\n";

        ++i;
    }
}

/******************************************************************************
    Compute the compatibility between 2 reads given all their bam records
******************************************************************************/
void computeReadsPairCompatibility(GlobalStatistics &globalStats,
                                   std::ofstream &asqgFile,
                                   std::ofstream &csvEdgesFile,
                                   int64_t i,
                                   int64_t j,
                                   std::vector<seqan::BamAlignmentRecord> const &readBamRecordBufferI,
                                   std::vector<seqan::BamAlignmentRecord> const &readBamRecordBufferJ,
                                   TProperties const &readNames,
                                   AlphaOptions const &options)
{
    bool atLeastOneCommonRef = false;
    bool areReadsCompatible = false;
    bool areReadsIncompatible = false;
    bool areReadsOverlapping = false;
    bool wasFoundWithMultiRef = false;

    bool trueCoEmittedReadsSpecie = (seqan::prefix(readNames[i], 3) == seqan::prefix(readNames[j], 3));
//    areReadsOverlapping = alignOvStats.totalOverlapPositions >= options.minOverlapLength;

    OverlapDescription alignOvDescription;

    // Test each pair of bamRecords from read_i vs. read_j
    // and exit the double loop as soon as a compatible alignment
    // has been found

    // TO DO: Envisager de parcourir la double boucle en proposant les paires
    // d'alignements de bonne qualité en premier si possible
    // int32_t min_align_num = std::min(readBamRecordBufferI.size(), readBamRecordBufferJ.size());

    for (auto const &bamRecordI : readBamRecordBufferI)
    {
        for (auto const &bamRecordJ : readBamRecordBufferJ)
        {
            if (bamRecordI.rID == bamRecordJ.rID)
            {
                ++globalStats.numConsideredAlignmentsPairs;

                alignOvDescription.isRead2RC = !(seqan::hasFlagRC(bamRecordI)==seqan::hasFlagRC(bamRecordJ));

                atLeastOneCommonRef = true;

                computeAlignmentsPairCompatibility(areReadsCompatible,
                                                   areReadsIncompatible,
                                                   areReadsOverlapping,
                                                   globalStats,
                                                   bamRecordI,
                                                   bamRecordJ,
                                                   options);
            }

            // Exit the nested loop if there is enough information
            // to decide whether the reads are compatible or not
            if (areReadsCompatible || areReadsIncompatible)
            {
                goto fastExit;
            }
        }
    }

    // There we never went through the fastExit
    // these reads are neither compatible neither incompatible
    // TO DO: Understand what is happening. WTF ???
    // Si ni compatible ni incompatoble --> statut particulier
    // les mettre dans un paquet pour les traiter plus tard ?
    ++globalStats.neitherCompatNorIncompatReadPairsNum;
    //    std::cerr << "Reads are neither compatibles nor incompatibles..." << "\n";

    areReadsIncompatible = true;

// goto tag to exit the nested loop just above
fastExit:

    if (atLeastOneCommonRef)
        ++globalStats.readPairsWithCommonRefNum;

    // If at least one alignment pair was compatible
    // we add an edge between the 2 reads

    OverlapStatistics alignOvStats;

    if (areReadsCompatible)
    {
        if (areReadsIncompatible)
        {
            std::cerr << "ERROR: 2 reads cannot be both compatible and incompatible, go check your algo..."
                      << "\n" << "\n" << std::flush;
        }

        //
        auto const &bamRecordI = readBamRecordBufferI[0];
        auto const &bamRecordJ = readBamRecordBufferJ[0];

        // De novo overlap alignment computing
        auto readSeqI(bamRecordI.seq);
        if (seqan::hasFlagRC(bamRecordI)) seqan::reverseComplement(readSeqI);

        auto readSeqJ(bamRecordJ.seq);
        if (seqan::hasFlagRC(bamRecordJ) && alignOvDescription.isRead2RC)
        {
            //nothing
        }
        else if (seqan::hasFlagRC(bamRecordJ) || alignOvDescription.isRead2RC)
            seqan::reverseComplement(readSeqJ);

        TAlign align;
        // int score = align2ReadSequences(align, readSeqI, readSeqJ);
        align2ReadSequences(align, readSeqI, readSeqJ);

        alignOvDescription.readIDI = i;
        alignOvDescription.readIDJ = j;

        computeAlignmentStats(alignOvStats, alignOvDescription, align, readSeqI, readSeqJ, false);

        auto seqIdentityPercent = static_cast<double>(alignOvStats.matchesNum) / alignOvStats.totalOverlapPositions;

        if (alignOvDescription.isEnclosedOverlap)
        {
            ++globalStats.enclosedOverlapsNum;
        }

        areReadsOverlapping = alignOvStats.totalOverlapPositions >= options.minOverlapLength;
        areReadsCompatible = (areReadsOverlapping && (alignOvStats.indelNum == 0)
                              && (seqIdentityPercent >= options.idRateThreshold));

        if (alignOvDescription.isRead2RC) reverseComplementReadJ(alignOvDescription);
    }

    if (areReadsCompatible)
    {
//        std::cerr << "Les 2 reads sont compatibles" << "\n";

        ++globalStats.compatibleReadPairsNum;

        // Add a new edge between the two reads
        if (options.outputASQG)
        {
            writeOverlapToASQG(asqgFile, alignOvDescription);
        }

        if (options.outputCSV)
        {
            writeOverlapToCSV(csvEdgesFile,
                              alignOvDescription,
                              trueCoEmittedReadsSpecie,
                              wasFoundWithMultiRef,
                              alignOvStats);
        }

        // Evaluate TP and FP
        if (trueCoEmittedReadsSpecie) ++globalStats.truePositive;
        else ++globalStats.falsePositive;
    }
    else // reads are considered incompatible by default
    {
        if (areReadsIncompatible)
        {
            ++globalStats.incompatibleReadPairsNum;
        }

//        std::cerr << "Les 2 reads ne sont pas compatibles" << "\n";
        // Evaluate FN et TN
        if(areReadsOverlapping)
        {
            if (trueCoEmittedReadsSpecie) ++globalStats.falseNegative;
            else ++globalStats.trueNegative;
        }
    }
}

/******************************************************************************
    Write an overlap to the ASQG output file
******************************************************************************/
void writeOverlapToASQG(std::ostream &asqgFile,
                        OverlapDescription const &aOvDcr)
{
    asqgFile << "ED\t" << aOvDcr.readIDI << " " << aOvDcr.readIDJ
             << " " << aOvDcr.ovBeginPosI << " " << aOvDcr.ovEndPosI
             << " " << aOvDcr.readLengthI
             << " " << aOvDcr.ovBeginPosJ << " " << aOvDcr.ovEndPosJ
             << " " << aOvDcr.readLengthJ << " " << aOvDcr.isRead2RC
             << " -1\n";
}

/******************************************************************************
    Write an overlap to the CSV Edges output file
******************************************************************************/
void writeOverlapToCSV(std::ostream &csvEdgesFile,
                       OverlapDescription const &aOvDcr,
                       bool const trueCoEmittedReadsSpecie,
                       bool const wasFoundWithMultiRef,
                       OverlapStatistics const &aOvStat)
{
    auto seqIdentityPercent = static_cast<double>(aOvStat.matchesNum) / aOvStat.totalOverlapPositions;

    csvEdgesFile << std::setprecision(4)
                 << aOvDcr.readIDI << ";" << aOvDcr.readIDJ << ";"
                 << "Undirected;"
                 << trueCoEmittedReadsSpecie << ";"
                 << (trueCoEmittedReadsSpecie ? 1 : 2) << ";"
                 << seqIdentityPercent << ";"
                 << wasFoundWithMultiRef << "\n";
}

/******************************************************************************
    Align 2 read sequences
******************************************************************************/
int align2ReadSequences(TAlign &align,
                        TSequence const &readSeqI,
                        TSequence const &readSeqJ)
{
    seqan::resize(seqan::rows(align), 2);
    seqan::assignSource(seqan::row(align, 0), readSeqI);
    seqan::assignSource(seqan::row(align, 1), readSeqJ);

//    auto scoringScheme = seqan::Score<int, seqan::Simple>(2, -4, -16, -3);
//    auto scoringScheme = seqan::Score<int, seqan::Simple>(1, -1, -2);
    auto scoringScheme = seqan::Score<int, seqan::Simple>(1, -1, -1);
    auto alignConfig = seqan::AlignConfig<true, true, true, true>();

//    int maxError = (int) (1 + std::min(length(readSeqI), length(readSeqJ)) * 0.05);

    int score = seqan::globalAlignment(align, scoringScheme, alignConfig);
//    int score = seqan::globalAlignment(align, scoringScheme, alignConfig, -maxError, maxError);

    return score;
}

/******************************************************************************
    Compute statistics for a pairwise alignment from Seqan
******************************************************************************/
void computeAlignmentStats(OverlapStatistics &alignOvStats,
                           OverlapDescription &alignOvDescription,
                           TAlign const &align,
                           TSequence const &readSeqI,
                           TSequence const &readSeqJ,
                           bool const toPrint)
{
    if (toPrint)
        std::cout << "\n" << align;

    auto &rowI = row(align, 0);
    auto &rowJ = row(align, 1);

    auto beginPosInAlignI = toViewPosition(rowI, 0);
    auto beginPosInAlignJ = toViewPosition(rowJ, 0);

    alignOvDescription.readLengthI = seqan::length(readSeqI);
    alignOvDescription.readLengthJ = seqan::length(readSeqJ);

    int lastPosI = seqan::length(readSeqI)-1;
    int lastPosJ = seqan::length(readSeqJ)-1;

    auto endPosInAlignI = toViewPosition(rowI, lastPosI);
    auto endPosInAlignJ = toViewPosition(rowJ, lastPosJ);

    int overlapBeginPos = std::max(beginPosInAlignI, beginPosInAlignJ);
    int overlapEndPos = std::min(endPosInAlignI, endPosInAlignJ);
    alignOvStats.totalOverlapPositions = overlapEndPos - overlapBeginPos + 1;

    alignOvDescription.ovBeginPosI = std::max(0, (int)toSourcePosition(rowI, beginPosInAlignJ));
    alignOvDescription.ovBeginPosJ = std::max(0, (int)toSourcePosition(rowJ, beginPosInAlignI));

    alignOvDescription.ovEndPosI = std::min(lastPosI, (int)toSourcePosition(rowI, endPosInAlignJ));
    alignOvDescription.ovEndPosJ = std::min(lastPosJ, (int)toSourcePosition(rowJ, endPosInAlignI));

    if ((alignOvDescription.ovBeginPosI == 0 && alignOvDescription.ovEndPosI == lastPosI)
        || (alignOvDescription.ovBeginPosJ == 0 && alignOvDescription.ovEndPosJ == lastPosJ))
    {
        alignOvDescription.isEnclosedOverlap = true;
    }

    if (!alignOvDescription.isEnclosedOverlap)
    {
        for (int i = overlapBeginPos; i <= overlapEndPos; ++i)
        {
            if (seqan::isGap(rowI, i) || seqan::isGap(rowJ, i))
            {
                ++alignOvStats.indelNum;
            }
            else
            {
                if (rowI[i] == rowJ[i])
                {
                    ++alignOvStats.matchesNum;
                }
                else
                {
                    ++alignOvStats.mismatchesNum;
                }
            }

            if (toPrint)
            {
                std::cout << i << " " << rowI[i] << " " << rowJ[i] << " "
                          << alignOvStats.matchesNum << " "
                          << alignOvStats.mismatchesNum << " "
                          << alignOvStats.indelNum << "\n";
            }
        }
    }

    if (toPrint)
    {
        std::cout << "\n" << alignOvStats.totalOverlapPositions << "\n"
                  << alignOvStats.matchesNum << "\n"
                  << alignOvStats.indelNum << "\n";
    }
}

/******************************************************************************
    Print the progress of the nested loop
******************************************************************************/
void printProgress(std::ostream &stream,
                   int32_t step,
                   int64_t numConsideredReadsPairs,
                   int64_t maxReadsPairs)
{
    if (numConsideredReadsPairs % step == 0)
    {
        double progress = (double) numConsideredReadsPairs * 100.0 / maxReadsPairs;
        stream << std::setprecision(5)
               << "\rINFO: Compatibility graph building progress at "
               << progress << "%         " << std::flush;
    }
}