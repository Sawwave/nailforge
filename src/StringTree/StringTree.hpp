#ifndef NAIL_FORGE_STRING_TREE_HPP
#define NAIL_FORGE_STRING_TREE_HPP

#include "../nailforge.hpp"

namespace NailForge::StringTree {

    /// @brief Data bundle for constant data for a full findSeeds invocation
    struct Context {
        Context(const AwFmIndex& fmIndex, const FastaVector& fastaVector,
            const P7Hmm& phmm, const std::vector<float>& matchScores, const NailForge::SearchParams& searchParams,
            const bool isReverseComplement)noexcept :
            fmIndex(fmIndex), fastaVector(fastaVector), phmm(phmm), matchScores(matchScores),
            searchParams(searchParams), isReverseComplement(isReverseComplement) {}
        const AwFmIndex& fmIndex;
        const FastaVector& fastaVector;
        const P7Hmm& phmm;
        const std::vector<float>& matchScores;
        const NailForge::SearchParams& searchParams;
        const bool isReverseComplement;
    };


    /// @brief Describes a position in model/sequence space where a StringTree hit occurrs.
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

    /// @brief finds seeds, first with the string tree, and then validates them using the SeedExtension module, and adds verified seeds to the seed list.
    /// @param context Context of the search
    /// @param seedList data structure to append hits to.
    void findSeeds(const StringTree::Context& context, std::vector<AlignmentSeed>& seedList) noexcept;
}


#endif