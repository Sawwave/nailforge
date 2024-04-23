#include "StringTree.hpp"
#include "../SeedExtension/SeedExtension.hpp"
#include "MaxExtensionTable/MaxExtensionTable.hpp"
#include "../PhmmProcessor/PhmmProcessor.hpp"
#include <vector>
#include <string>
#include <cstring>



namespace NailForge::StringTree {


    struct StringTreeSparseScoreVectorEntry {
        StringTreeSparseScoreVectorEntry(const uint32_t modelPosition, const float score) :
            modelPosition(modelPosition), score(score), maxScore(score) {}
        StringTreeSparseScoreVectorEntry(const uint32_t modelPosition, const float score, const float maxScore) :
            modelPosition(modelPosition), score(score), maxScore(maxScore) {}
        uint32_t modelPosition;
        float score;
        float maxScore;
    };

    struct StringTreeStackEntry {
        AwFmSearchRange searchRange;
        uint8_t symbol;
        std::vector<StringTreeSparseScoreVectorEntry> diagonalEntries;
    };

    void verifyDiagonalsPassingThreshold(const StringTree::Context& context, const StringTreeStackEntry& stackEntry,
        const uint32_t& thisSymbolModelPosition, const float maxScoreAlongDiagonal, std::vector<AlignmentSeed>& seedList,
        const uint8_t hitLength);


    void findSeeds(const StringTree::Context& context, std::vector<AlignmentSeed>& seedList) {
        
        const auto awfmBackwardStepFunction = context.phmm.header.alphabet == P7HmmReaderAlphabetAmino ?
            awFmAminoIterativeStepBackwardSearch : awFmNucleotideIterativeStepBackwardSearch;

        const uint32_t nextDiagonalCellOffset = context.isReverseCompliment ? 1 : -1;
        const uint8_t modelSymbolEncodingBitmask = context.isReverseCompliment ? 0b11 : 0x00;
        const uint8_t alphabetSize = context.phmm.header.alphabet == P7HmmReaderAlphabetAmino ? 20 : 4;
        const uint32_t modelLength = context.phmm.header.modelLength;

        //set up the stacks that allow us to perform the recursion iteratively
        std::vector<StringTreeStackEntry> stack;
        stack.resize(context.searchParams.maximumHitLength);

        Table::MaxExtensionTable maxExtensionTable =
            Table::MaxExtensionTable(modelLength, context.searchParams.maximumHitLength, alphabetSize);
        maxExtensionTable.populateTable(context.matchScores, context.isReverseCompliment);


        //string tree loop begins here
        int_fast8_t currentDepth = 0;
        stack[0].symbol = 0;
        while (currentDepth >= 0) {
            const uint8_t maxExtensionPositionsRemaining = (context.searchParams.maximumHitLength - 1) - currentDepth;
            const auto& currentSymbol = stack[currentDepth].symbol;
            //character encoding is the character after complimenting (if doing reverse compliment search)
            const uint8_t phmmSymbolIdx = currentSymbol ^ modelSymbolEncodingBitmask;
            stack[currentDepth].diagonalEntries.clear();

            if (__builtin_expect(currentDepth == 0, false)) {
                //create search range from first symbol
                stack[0].searchRange.startPtr = context.fmIndex.prefixSums[currentSymbol];
                stack[0].searchRange.endPtr = context.fmIndex.prefixSums[currentSymbol + 1] - 1;
                for (uint32_t modelPosition = 0; modelPosition < modelLength; modelPosition++) {
                    //get the match score. note that we look up the non-complimented score, but we'll write down the complimented character later
                    const float matchScore = context.matchScores[(modelPosition * alphabetSize) + phmmSymbolIdx];

                    //only try extending this position if the match score was positive.
                    if (__builtin_expect(matchScore > 0, false)) {
                        const float bestExtension = maxExtensionTable.getMaxExtensionScore(modelPosition, context.searchParams.maximumHitLength - 1);
                        const float bestScorePossible = matchScore + bestExtension;
                        if (bestScorePossible >= context.searchParams.mainDiagonalThresholdScore) {
                            stack[0].diagonalEntries.emplace_back(modelPosition, matchScore);
                        }
                    }
                }
            }

            else {
                const auto& prevDiagonalEntries = stack[currentDepth - 1].diagonalEntries;

                bool hasResolvedSearchRange = false;
                for (const auto& diagonalEntry : prevDiagonalEntries) {

                    const uint32_t thisSymbolModelPosition = diagonalEntry.modelPosition + nextDiagonalCellOffset;
                    const float matchScore = context.matchScores[(thisSymbolModelPosition * alphabetSize) + phmmSymbolIdx];
                    const float accumulatedScore = matchScore + diagonalEntry.score;
                    const float maxScoreSeen = std::max(accumulatedScore, diagonalEntry.maxScore);
                    if (accumulatedScore > 0.0f) {
                        const float bestExtension = maxExtensionTable.getMaxExtensionScore(thisSymbolModelPosition, maxExtensionPositionsRemaining);
                        const float bestScorePossible = accumulatedScore + bestExtension;

                        const bool accumulatedScorePassesThreshold = accumulatedScore >= context.searchParams.mainDiagonalThresholdScore;
                        const bool bestExtensionPossiblePassesThreshold = bestScorePossible >= context.searchParams.mainDiagonalThresholdScore;

                        //resolve the search range if we haven't yet.
                        if (__builtin_expect(!hasResolvedSearchRange && bestExtensionPossiblePassesThreshold, false)) {
                            stack[currentDepth].searchRange = stack[currentDepth - 1].searchRange;
                            awfmBackwardStepFunction(&context.fmIndex, &stack[currentDepth].searchRange, currentSymbol);
                            hasResolvedSearchRange = true;
                        }

                        if (__builtin_expect(accumulatedScorePassesThreshold && (currentDepth < (context.searchParams.maximumHitLength)), false)) {
                            verifyDiagonalsPassingThreshold(context, stack[currentDepth], thisSymbolModelPosition, maxScoreSeen,
                                seedList, currentDepth);

                        }
                        //if we didn't pass the threshold but if fully extended this kmer might, add it to the
                        //next layer diagonal list
                        else if (bestExtensionPossiblePassesThreshold) {
                            stack[currentDepth].diagonalEntries.emplace_back(thisSymbolModelPosition, accumulatedScore,
                                maxScoreSeen);
                        }
                    }
                }
            }

            //roll over the current symbol, and current depth if necessary
            bool exploreFurtherDownTree =
                ((currentDepth) < (context.searchParams.maximumHitLength - 1)) &&
                stack[currentDepth].diagonalEntries.size() != 0 &&
                (awFmSearchRangeLength(&stack[currentDepth].searchRange) != 0);

            if (exploreFurtherDownTree) {
                stack[++currentDepth].symbol = 0;
            }
            else {
                while (stack[currentDepth].symbol == alphabetSize - 1 && currentDepth >= 0) {
                    stack[currentDepth--].symbol = 0;
                }
                if (__builtin_expect(currentDepth >= 0, true)) {
                    stack[currentDepth].symbol++;
                }
            }

        }
    }

    void verifyDiagonalsPassingThreshold(const StringTree::Context& context, const StringTreeStackEntry& stackEntry,
        const uint32_t& thisSymbolModelPosition, const float maxScoreAlongDiagonal, std::vector<AlignmentSeed>& seedList,
        const uint8_t hitLength) {

        uint64_t searchRangeLength = awFmSearchRangeLength(&stackEntry.searchRange);
        //if there are way too many sequence hits, it's probably in some repetitive data, and we shouldn't report it.
        const double hitsPerMillionPositions = (double)searchRangeLength / (((double)context.fastaVector.sequence.count / (double)1e6));
        if (hitsPerMillionPositions > context.searchParams.maxSeqHitsPerMillion) {
            return;
        }
        for (uint64_t bwtPosition = stackEntry.searchRange.startPtr;
            bwtPosition <= stackEntry.searchRange.endPtr;bwtPosition++) {
            uint64_t globalSequencePosition, localSequencePosition, sequenceIdx;
            AwFmReturnCode awfmrc;

            globalSequencePosition = awFmFindDatabaseHitPositionSingle(&context.fmIndex,
                bwtPosition, &awfmrc);
            if (__builtin_expect(awfmrc != AwFmFileReadOkay, false)) {
                std::cout << "Error: could not resolve bwt position into global sequence position" << std::endl;
            }

            //get the global position for this bwt position
            awfmrc = awFmGetLocalSequencePositionFromIndexPosition(&context.fmIndex,
                globalSequencePosition, &sequenceIdx, &localSequencePosition);
            if (__builtin_expect(awfmrc != AwFmSuccess, false)) {
                std::cout << "Error: could not resolve local sequence position from fm index" << std::endl;
            }

            const StringTree::HitPosition hitPosition(sequenceIdx, localSequencePosition, thisSymbolModelPosition, hitLength);
            bool seedIsVerified = SeedExtension::verifySeedViaExtension(context, hitPosition, maxScoreAlongDiagonal);
            if (seedIsVerified) {
                seedList.emplace_back(localSequencePosition, thisSymbolModelPosition, sequenceIdx);

            }
        }
    }
}