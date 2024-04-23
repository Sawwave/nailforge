#ifndef NAIL_FORGE_STRING_TREE_HPP
#define NAIL_FORGE_STRING_TREE_HPP

#include "../nailforge.hpp"

namespace NailForge::StringTree {

    struct Context {
        Context(const AwFmIndex& fmIndex, const FastaVector& fastaVector,
            const P7Hmm& phmm, const std::vector<float>& matchScores, const NailForge::SearchParams& searchParams,
            const float extensionThresholdScore, const bool isReverseCompliment)noexcept :
            fmIndex(fmIndex), fastaVector(fastaVector), phmm(phmm), matchScores(matchScores),
            searchParams(searchParams), extensionThresholdScore(extensionThresholdScore),
            isReverseCompliment(isReverseCompliment) {}
        const AwFmIndex& fmIndex;
        const FastaVector& fastaVector;
        const P7Hmm& phmm;
        const std::vector<float>& matchScores;
        const NailForge::SearchParams& searchParams;
        const float extensionThresholdScore;
        const bool isReverseCompliment;
    };

    struct HitPosition {
        HitPosition(const uint64_t sequenceIdx, const uint64_t localSequencePosition,
            const uint32_t modelPosition, const uint8_t hitLength) noexcept :
            sequenceIdx(sequenceIdx), localSequencePosition(localSequencePosition),
            modelPosition(modelPosition), hitLength(hitLength) {}
        const uint64_t sequenceIdx;
        const uint64_t localSequencePosition;
        const uint32_t modelPosition;
        const uint8_t hitLength;
    };


    void findSeeds(const StringTree::Context& context, std::vector<AlignmentSeed>& seedList);
}


#endif