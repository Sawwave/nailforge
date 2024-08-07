#ifndef NAIL_FORGE_SEED_EXTENSION_HPP
#define NAIL_FORGE_SEED_EXTENSION_HPP

#include "../nailforge.hpp"
#include "../Alphabet/LetterConversion.hpp"
#include "../StringTree/StringTree.hpp"
#include <iostream>

namespace NailForge::SeedExtension {


    struct ExtensionResult {
        ExtensionResult(const float maximumScore, const bool isVerified) :
            maximumScore(maximumScore), isVerified(isVerified) {}
        float maximumScore;
        bool isVerified;
    };

    struct ScoreLengthTuple {
        ScoreLengthTuple() {}
        ScoreLengthTuple(const float score, const uint32_t length) :
            score(score), length(length) {}
        float score;
        uint32_t length;
    };

    /// @brief Verifies a hit from the String Tree by extending the flanking regions of the diagonal to attempt to pass the final threshold score
    /// @param context String Tree search context data
    /// @param hitPosition position of the hit
    /// @param maxScoreAlongDiagonal maximum score seen along the initial diagonal
    /// @return struct containing the final score, and whether it was verified
    ExtensionResult verifySeedViaExtension(const StringTree::Context& context,
        const StringTree::HitPosition& hitPosition, const float maxScoreAlongDiagonal);


    /// @brief 
    /// @tparam isReverseComplement determines if the implementation should be aligning as the reverse complement strand 
    /// @param sequencePtr pointer to the first character of the sequence that the hit occurs in 
    /// @param sequenceLength length of the sequence the hit occurs in
    /// @param context string tree search context data
    /// @param hitPosition position of the hit in model/sequence space
    /// @param isPriorFlank determines if the computation should be looking at the prior flank or post flank
    /// @return score of the flanking region, and how long the flank was. This does not include the original diagonal's score.
    template<bool isReverseComplement>
    ScoreLengthTuple findFlankingDiag(const char* sequencePtr, const uint64_t sequenceLength, const StringTree::Context& context,
        const StringTree::HitPosition& hitPosition, const float thresholdScore, const bool isPriorFlank) noexcept;


    template<bool isReverseComplement>
    ScoreLengthTuple findFlankingDiag(const char* sequencePtr, const uint64_t sequenceLength, const StringTree::Context& context,
        const StringTree::HitPosition& hitPosition, const float thresholdScore, const bool isPriorFlank) noexcept {

        const uint8_t symbolEncodingComplementBitmask = isReverseComplement ? 0x03 : 0x00;
        const NailForge::Alphabet alphabet = context.phmm.header.alphabet == P7HmmReaderAlphabetAmino ?
            NailForge::Alphabet::Amino : NailForge::Alphabet::Dna;
        const auto alphabetSize = NailForge::getAlphabetSize(alphabet);

        float accumulatedScore = 0.0f;
        float maxAccumulatedScore = 0.0f;
        uint64_t flankLength;

        if constexpr (isReverseComplement) {
            flankLength = isPriorFlank ?
                std::min(hitPosition.localSequencePosition,
                    static_cast<uint64_t>(context.phmm.header.modelLength - hitPosition.modelPosition - 1)) :
                std::min(sequenceLength - hitPosition.localSequencePosition - hitPosition.hitLength,
                    static_cast<uint64_t>(hitPosition.modelPosition - hitPosition.hitLength));
        }
        else {
            flankLength = isPriorFlank ?
                std::min(hitPosition.localSequencePosition, static_cast<uint64_t>(hitPosition.modelPosition)) :
                std::min(sequenceLength - hitPosition.localSequencePosition - hitPosition.hitLength,
                    static_cast<uint64_t>(context.phmm.header.modelLength - hitPosition.modelPosition - hitPosition.hitLength));

        }
        flankLength = std::min(flankLength, static_cast<uint64_t>(context.searchParams.flankExtensionLength));

        for (uint32_t flankPosition = 0; flankPosition < flankLength; flankPosition++) {
            uint32_t flankModelPosition;
            uint32_t flankSequencePosition;

            flankSequencePosition = isPriorFlank ?
                hitPosition.localSequencePosition - flankPosition - 1 :
                hitPosition.localSequencePosition + flankPosition + hitPosition.hitLength;

            if constexpr (isReverseComplement) {
                flankModelPosition = isPriorFlank ?
                    hitPosition.modelPosition + flankPosition + 1 :
                    hitPosition.modelPosition - flankPosition - hitPosition.hitLength;
            }
            else {
                flankModelPosition = isPriorFlank ?
                    hitPosition.modelPosition - flankPosition - 1 :
                    hitPosition.modelPosition + flankPosition + hitPosition.hitLength;
            }

            const char sequenceSymbol = sequencePtr[flankSequencePosition];
            const uint8_t symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);
            const uint8_t phmmSymbolIdx = symbolIdx ^ symbolEncodingComplementBitmask;
            const float matchScore = context.matchScores[(flankModelPosition * alphabetSize) + phmmSymbolIdx];

            accumulatedScore += matchScore;
            maxAccumulatedScore = std::max(maxAccumulatedScore, accumulatedScore);
            //here, we are killing the flank accumulation if the flank's accumulated score goes below -1*threshold.
            //the assumption here is that a flank this bad would require the threshold to be hit by the main diag and the
            //other flank, without any contribution from this flank. I think that flanks could probably be killed much earlier,
            //but I haven't done any analysis to find a good place to kill it.
            //note that the *max* accumulated score is returned, so this can be positive even if the flank is killed in this way.
            if (__builtin_expect(accumulatedScore < (-thresholdScore), false)) {
                return ScoreLengthTuple(maxAccumulatedScore, flankPosition);
            }
        }
        return ScoreLengthTuple(maxAccumulatedScore, flankLength);
    }
}

#endif