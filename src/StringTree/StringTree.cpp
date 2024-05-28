#include "StringTree.hpp"
#include "../SeedExtension/SeedExtension.hpp"
#include "MaxExtensionTable/MaxExtensionTable.hpp"
#include "../PhmmProcessor/PhmmProcessor.hpp"
#include <vector>
#include <string>
#include <cstring>
#include <optional>


namespace NailForge::StringTree {

    /// @brief Intermediary data structure detailing a position that currently passes the string tree, along with the score at that position.
    struct StringTreeSparseScoreVectorEntry {
        StringTreeSparseScoreVectorEntry(const uint32_t modelPosition, const float score) :
            modelPosition(modelPosition), score(score), maxScore(score) {}
        StringTreeSparseScoreVectorEntry(const uint32_t modelPosition, const float score, const float maxScore) :
            modelPosition(modelPosition), score(score), maxScore(maxScore) {}
        uint32_t modelPosition;
        float score;
        float maxScore;
    };

    /// @brief List of the positions passing the current node at the string tree.
    struct StringTreeStackEntry {
        AwFmSearchRange searchRange;
        uint8_t symbol;
        std::vector<StringTreeSparseScoreVectorEntry> diagonalEntries;
    };

    /// @brief struct storing the resolved FM-index position from BWT-space, now in sequence-space
    struct SequencePosition {
        SequencePosition(const uint64_t bwtPosition, const AwFmIndex& fmIndex)noexcept {
            uint64_t globalSequencePosition;
            AwFmReturnCode awfmrc;

            globalSequencePosition = awFmFindDatabaseHitPositionSingle(&fmIndex, bwtPosition, &awfmrc);
            if (__builtin_expect(awfmrc != AwFmFileReadOkay, false)) {
                std::cout << "Error: could not resolve bwt position into global sequence position" << std::endl;
            }
            //get the global position for this bwt position
            awfmrc = awFmGetLocalSequencePositionFromIndexPosition(&fmIndex,
                globalSequencePosition, &sequenceIdx, &localSequencePosition);
            if (__builtin_expect(awfmrc != AwFmSuccess, false)) {
                std::cout << "Error: could not resolve local sequence position from fm index" << std::endl;
            }

        }
        uint64_t localSequencePosition;
        uint64_t sequenceIdx;
    };


    /// @brief Prepares a node that has positions passing the score threshold for extension via the SeedExtension module
    /// @param context Search context for this model
    /// @param searchRange FM-index search range representing positions in the target for the string represented by the current path through the string tree.  
    /// @param resolvedSequencePositionList optional list of resolved positions. if this list is nullopt, this function will use the fm index to find the positions.
    /// @param thisSymbolModelPosition position in the model that this node represents
    /// @param maxScoreAlongDiagonal best score seen along the entire diagonal
    /// @param seedList list of seeds to append any verified extensions to
    /// @param hitLength length of the hit from the string tree, expected to be currentDepth
    void verifyDiagonalsPassingThreshold(const StringTree::Context& context, const AwFmSearchRange& searchRange,
        std::optional<std::vector<SequencePosition>>& resolvedSequencePositionList, const uint32_t& thisSymbolModelPosition,
        const float maxScoreAlongDiagonal, std::vector<AlignmentSeed>& seedList, const uint8_t hitLength) noexcept;


    void findSeeds(const StringTree::Context& context, std::vector<AlignmentSeed>& seedList) noexcept {

        const auto awfmBackwardStepFunction = context.phmm.header.alphabet == P7HmmReaderAlphabetAmino ?
            awFmAminoIterativeStepBackwardSearch : awFmNucleotideIterativeStepBackwardSearch;

        const uint32_t nextDiagonalCellOffset = context.isReverseComplement ? 1 : -1;
        const uint8_t modelSymbolEncodingBitmask = context.isReverseComplement ? 0b11 : 0x00;
        const uint8_t alphabetSize = context.phmm.header.alphabet == P7HmmReaderAlphabetAmino ? 20 : 4;
        const uint32_t modelLength = context.phmm.header.modelLength;

        //set up the stacks that allow us to perform the recursion iteratively
        std::vector<StringTreeStackEntry> stack;
        stack.resize(context.searchParams.maximumHitLength);

        Table::MaxExtensionTable maxExtensionTable =
            Table::MaxExtensionTable(modelLength, context.searchParams.maximumHitLength, alphabetSize);
        maxExtensionTable.populateTable(context.matchScores, context.isReverseComplement);


        //string tree loop begins here
        int_fast8_t currentDepth = 0;
        stack[0].symbol = 0;
        while (currentDepth >= 0) {
            auto& thisSearchRange = stack[currentDepth].searchRange;

            //default initialized to nullopt for every new string tree node
            std::optional<std::vector<SequencePosition>> resolvedSequencePositionList;

            const uint8_t maxExtensionPositionsRemaining = (context.searchParams.maximumHitLength - 1) - currentDepth;
            const auto& currentSymbol = stack[currentDepth].symbol;
            //character encoding is the character after complementing (if doing reverse complement search)
            const uint8_t phmmSymbolIdx = currentSymbol ^ modelSymbolEncodingBitmask;
            stack[currentDepth].diagonalEntries.clear();

            if (__builtin_expect(currentDepth == 0, false)) {
                //create search range from first symbol
                stack[0].searchRange.startPtr = context.fmIndex.prefixSums[currentSymbol];
                stack[0].searchRange.endPtr = context.fmIndex.prefixSums[currentSymbol + 1] - 1;
                for (uint32_t modelPosition = 0; modelPosition < modelLength; modelPosition++) {
                    //get the match score. note that we look up the non-complemented score, but we'll write down the complemented character later
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
                const auto& prevSearchRange = stack[currentDepth - 1].searchRange;
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
                            thisSearchRange = prevSearchRange;
                            awfmBackwardStepFunction(&context.fmIndex, &thisSearchRange, currentSymbol);
                            hasResolvedSearchRange = true;
                        }

                        if (__builtin_expect(accumulatedScorePassesThreshold && 
                        (currentDepth < static_cast<int32_t>(context.searchParams.maximumHitLength)), false)) {
                            verifyDiagonalsPassingThreshold(context, thisSearchRange, resolvedSequencePositionList,
                                thisSymbolModelPosition, maxScoreSeen, seedList, currentDepth);

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
                ((currentDepth) < (static_cast<int32_t>(context.searchParams.maximumHitLength) - 1)) &&
                stack[currentDepth].diagonalEntries.size() != 0 &&
                (awFmSearchRangeLength(&thisSearchRange) != 0);

            if (exploreFurtherDownTree) {
                stack[++currentDepth].symbol = 0;
            }
            else {
                while (currentDepth >= 0 && stack[currentDepth].symbol == alphabetSize - 1) {
                    stack[currentDepth--].symbol = 0;
                }
                if (__builtin_expect(currentDepth >= 0, true)) {
                    stack[currentDepth].symbol++;
                }
            }

        }
    }

    void verifyDiagonalsPassingThreshold(const StringTree::Context& context, const AwFmSearchRange& searchRange,
        std::optional<std::vector<SequencePosition>>& resolvedSequencePositionList, const uint32_t& thisSymbolModelPosition,
        const float maxScoreAlongDiagonal, std::vector<AlignmentSeed>& seedList, const uint8_t hitLength) noexcept {

        const uint64_t searchRangeLength = awFmSearchRangeLength(&searchRange);
        if (searchRangeLength == 0) {
            return;
        }
        //if there are way too many sequence hits, it's probably in some repetitive data, and we shouldn't report it.
        const double hitsPerMillionPositions = (double)searchRangeLength / (((double)context.fastaVector.sequence.count / (double)1e6));
        if (hitsPerMillionPositions > context.searchParams.maxSeqHitsPerMillion) {
            return;
        }
        if (!resolvedSequencePositionList.has_value()) {
            resolvedSequencePositionList = std::vector<SequencePosition>();
            resolvedSequencePositionList->reserve(searchRangeLength);

            //resolve the search range into a list of actual positions in the target sequence.
            for (uint64_t bwtPosition = searchRange.startPtr;
                bwtPosition <= searchRange.endPtr; bwtPosition++) {
                resolvedSequencePositionList->emplace_back(bwtPosition, context.fmIndex);
            }
        }
        for (const auto& resolvedSequencePosition : *resolvedSequencePositionList) {
            const StringTree::HitPosition hitPosition(resolvedSequencePosition.sequenceIdx,
                resolvedSequencePosition.localSequencePosition, thisSymbolModelPosition, hitLength);
            const auto verificationResult = SeedExtension::verifySeedViaExtension(context, hitPosition, maxScoreAlongDiagonal);
            if (verificationResult.isVerified) {
                seedList.emplace_back(resolvedSequencePosition.localSequencePosition,
                    thisSymbolModelPosition, resolvedSequencePosition.sequenceIdx, verificationResult.maximumScore);
            }
        }
    }
}