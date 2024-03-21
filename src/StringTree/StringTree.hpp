#ifndef NAIL_FORGE_STRING_TREE_HPP
#define NAIL_FORGE_STRING_TREE_HPP

#include "../nailforge.hpp"

namespace NailForge::StringTree {

    void findSeeds(const AwFmIndex* const fmIndex, const FastaVector& fastaVector, const P7Hmm& phmm,
        const uint32_t modelIdx, const std::vector<float>& matchScores, std::vector<AlignmentSeed>& seedList,
        const NailForge::SearchParams& params, const bool isReverseComplimentSearch);
}


#endif