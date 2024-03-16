#include "StringTree.hpp"
#include "../SeedExtension/SeedExtension.hpp"
#include "MaxExtensionTable/MaxExtensionTable.hpp"
#include "../PhmmProcessor/PhmmProcessor.hpp"
#include <vector>

namespace NailForge::StringTree {


    struct StringTreeSparseScoreVectorEntry {
        StringTreeSparseScoreVectorEntry(const uint32_t modelPosition, const float score) :
            modelPosition(modelPosition), score(score) {}
        uint32_t modelPosition;
        float score;
    };

    struct StringTreeStackEntry {
        AwFmSearchRange searchRange;
        uint8_t symbol;
        std::vector<StringTreeSparseScoreVectorEntry> diagonalEntries;
    };

    void verifyDiagonalsPassingThreshold(const AwFmIndex* const fmIndex, const StringTreeStackEntry& stackEntry,
        const uint32_t& thisSymbolModelPosition, const NailForge::SearchParams& params, const bool isReverseCompliment,
        const FastaVector& fastaVector, const P7Hmm& phmm);


    void findSeeds(const AwFmIndex* const fmIndex, const FastaVector& fastaVector, const P7Hmm& phmm,
        const uint32_t modelIdx, const std::vector<float>& matchScores, const std::vector<AlignmentSeed>& seedList,
        const NailForge::SearchParams& params, const bool isReverseComplimentSearch) {

        const uint32_t nextDiagonalCellOffset = isReverseComplimentSearch ? 1 : -1;
        const uint8_t symbolEncodingComplimentBitmask = isReverseComplimentSearch ? 0x03 : 0x00;
        const uint8_t alphabetSize = phmm.header.alphabet == P7HmmReaderAlphabetAmino ? 20 : 4;
        const uint32_t modelLength = phmm.header.modelLength;

        //set up the stacks that allow us to perform the recursion iteratively
        std::vector<StringTreeStackEntry> stack;
        stack.resize(params.maximumHitLength);
        for (auto& entry : stack) {
            entry.diagonalEntries.reserve(modelLength);
        }

        Table::MaxExtensionTable maxExtensionTable =
            Table::MaxExtensionTable(modelLength, params.maximumHitLength, alphabetSize);
        maxExtensionTable.populateTable(matchScores, isReverseComplimentSearch);


        //string tree loop begins here
        int8_t currentDepth = 0;
        stack[0].symbol = 0;
        while (currentDepth >= 0) {
            const uint8_t maxExtensionPositionsRemaining = (params.maximumHitLength - 1) - currentDepth;
            auto& currentSymbol = stack[currentDepth].symbol;
            //character encoding is the character after complimenting (if doing reverse compliment search)
            const uint8_t phmmSymbolIdx = currentSymbol ^ symbolEncodingComplimentBitmask;
            stack[currentDepth].diagonalEntries.clear();

            if (__builtin_expect(currentDepth == 0, false)) {
                //create search range from first symbol
                stack[0].searchRange.startPtr = fmIndex->prefixSums[currentSymbol];
                stack[0].searchRange.endPtr = fmIndex->prefixSums[currentSymbol + 1] - 1;
                for (uint32_t modelPosition = 0; modelPosition < modelLength;modelPosition++) {
                    //get the match score. note that we look up the non-complimented score, but we'll write down the complimented character later
                    const float matchScore = matchScores[(modelPosition * alphabetSize) + phmmSymbolIdx];

                    //only try extending this position if the match score was positive.
                    if (__builtin_expect(matchScore > 0, false)) {
                        const float bestExtension = maxExtensionTable.getMaxExtensionScore(modelPosition, params.maximumHitLength - 1);
                        const float bestScorePossible = matchScore + bestExtension;
                        if (bestScorePossible >= params.mainDiagonalThresholdScore) {
                            stack[0].diagonalEntries.emplace_back(modelPosition, matchScore);
                        }
                    }
                }
            }

            else {
                const auto& prevDiagonalEntries = stack[currentDepth - 1].diagonalEntries;

                for (const auto& diagonalEntry : prevDiagonalEntries) {
                    bool hasResolvedSearchRange = false;

                    const uint32_t thisSymbolModelPosition = diagonalEntry.modelPosition + nextDiagonalCellOffset;
                    const float matchScore = matchScores[(thisSymbolModelPosition * alphabetSize) + phmmSymbolIdx];
                    const float accumulatedScore = matchScore + diagonalEntry.score;
                    if (accumulatedScore > 0.0f) {
                        const float bestExtension = maxExtensionTable.getMaxExtensionScore(thisSymbolModelPosition, maxExtensionPositionsRemaining);
                        const float bestScorePossible = accumulatedScore + bestExtension;
                        const bool accumulatedScorePassesThreshold = accumulatedScore >= params.mainDiagonalThresholdScore;
                        const bool bestExtensionPossiblePassesThreshold = bestScorePossible >= params.mainDiagonalThresholdScore;

                        //resolve the search range if we haven't yet.
                        if (__builtin_expect(!hasResolvedSearchRange && bestExtensionPossiblePassesThreshold, false)) {
                            stack[currentDepth].searchRange = stack[currentDepth - 1].searchRange;
                            if (phmm.header.alphabet == P7HmmReaderAlphabetAmino) {
                                awFmAminoIterativeStepBackwardSearch(fmIndex, &stack[currentDepth].searchRange, currentSymbol);
                            }
                            else {
                                awFmNucleotideIterativeStepBackwardSearch(fmIndex, &stack[currentDepth].searchRange, currentSymbol);
                            }
                            hasResolvedSearchRange = true;
                        }

                        if (__builtin_expect(accumulatedScorePassesThreshold, false)) {
                            verifyDiagonalsPassingThreshold(fmIndex, stack[currentDepth], thisSymbolModelPosition,
                                params, isReverseComplimentSearch, fastaVector, phmm);
                        }
                        //if we didn't pass the threshold but if fully extended this kmer might, add it to the
                        //next layer diagonal list
                        else if (bestExtensionPossiblePassesThreshold) {
                            stack[currentDepth].diagonalEntries.emplace_back(thisSymbolModelPosition, accumulatedScore);
                        }
                    }
                }
            }

            //roll over the current symbol, and current depth if necessary
            bool exploreFurtherDownTree =
                ((currentDepth + 1) < params.maximumHitLength) &&
                stack[currentDepth].diagonalEntries.size() != 0 &&
                (awFmSearchRangeLength(&stack[currentDepth].searchRange) != 0);

            if (exploreFurtherDownTree) {
                stack[++currentDepth].symbol = 0;
            }
            else if (currentSymbol == alphabetSize - 1) {
                currentSymbol = 0;
                currentDepth--;
            }
            else {
                currentSymbol++;
            }
        }
    }

    void verifyDiagonalsPassingThreshold(const AwFmIndex* const fmIndex, const StringTreeStackEntry& stackEntry,
        const uint32_t& thisSymbolModelPosition, const NailForge::SearchParams& params, const bool isReverseCompliment,
        const FastaVector& fastaVector, const P7Hmm& phmm) {

        uint64_t searchRangeLength = awFmSearchRangeLength(&stackEntry.searchRange);
        //if there are way too many sequence hits, it's probably in some repetitive data, and we shouldn't report it.
        const float hitsPerMillionPositions = searchRangeLength / ((float)fastaVector.sequence.count / 1e6);
        if (hitsPerMillionPositions > params.maxSeqHitsPerMillion) {
            return;
        }

        for (uint64_t bwtPosition = stackEntry.searchRange.startPtr;
            bwtPosition <= stackEntry.searchRange.endPtr;bwtPosition++) {
            uint64_t globalSequencePosition, localSequencePosition, sequenceIdx;
            AwFmReturnCode awfmrc;

            globalSequencePosition = awFmFindDatabaseHitPositionSingle(fmIndex,
                bwtPosition, &awfmrc);
            if (__builtin_expect(awfmrc != AwFmFileReadOkay, false)) {
                std::cout << "Error: could not resolve bwt position into global sequence position" << std::endl;
            }


            //get the global position for this bwt position
            awfmrc = awFmGetLocalSequencePositionFromIndexPosition(fmIndex,
                globalSequencePosition, &sequenceIdx, &localSequencePosition);
            if (__builtin_expect(awfmrc != AwFmSuccess, false)) {
                std::cout << "Error: could not resolve local sequence position from fm index" << std::endl;
            }

            SeedExtension::verifySeedViaExtension(localSequencePosition,
                thisSymbolModelPosition, sequenceIdx, params, isReverseCompliment, fastaVector, phmm);
        }
    }

}