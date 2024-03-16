#ifndef NAIL_FORGE_SEED_EXTENSION_HPP
#define NAIL_FORGE_SEED_EXTENSION_HPP

#include "../nailforge.hpp"

namespace NailForge::SeedExtension {

    bool verifySeedViaExtension(const uint64_t localSequencePosition, const uint32_t modelPosition,
        const uint64_t sequenceIdx, const NailForge::SearchParams& searchParams, const bool isReverseCompliment,
        const FastaVector& fastaVector, const P7Hmm& phmm);
}

#endif