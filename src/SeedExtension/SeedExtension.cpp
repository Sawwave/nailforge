#include "SeedExtension.hpp"
#include "../PhmmProcessor/PhmmProcessor.hpp"
#include "../Alphabet/LetterConversion.hpp"
#include "iostream"
#include <array>
#include <algorithm>
#include <immintrin.h>

namespace NailForge::SeedExtension {

    bool verifySeedViaExtension(const StringTree::Context& context,
        const StringTree::HitPosition& hitPosition, const float maxScoreAlongDiagonal) {

        char* sequencePtr;
        size_t sequenceLength;
        fastaVectorGetSequence(&context.fastaVector, hitPosition.sequenceIdx, &sequencePtr, &sequenceLength);
        if (sequencePtr == NULL) {
            std::cerr << "Error: could not read sequence from fasta at given sequence idx" << std::endl;
            return false;
        }

        float priorFlankScore, postFlankScore;
        if (context.isReverseCompliment) {
            priorFlankScore = findFlankingDiag<true>(sequencePtr, sequenceLength, context, hitPosition, true);
            postFlankScore = findFlankingDiag<false>(sequencePtr, sequenceLength, context, hitPosition, false);
        }
        else {
            priorFlankScore = findFlankingDiag<true>(sequencePtr, sequenceLength, context, hitPosition, true);
            postFlankScore = findFlankingDiag<false>(sequencePtr, sequenceLength, context, hitPosition, false);
        }

        return (priorFlankScore + postFlankScore + maxScoreAlongDiagonal) >= context.extensionThresholdScore;
    }
}
