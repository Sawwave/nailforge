#include "SeedExtension.hpp"
#include "../PhmmProcessor/PhmmProcessor.hpp"
#include "../Alphabet/LetterConversion.hpp"
#include "iostream"
#include <array>
#include <algorithm>
#include <immintrin.h>

namespace NailForge::SeedExtension {

    ExtensionResult verifySeedViaExtension(const StringTree::Context& context,
        const StringTree::HitPosition& hitPosition, const float maxScoreAlongDiagonal) {

        char* sequencePtr;
        size_t sequenceLength;
        fastaVectorGetSequence(&context.fastaVector, hitPosition.sequenceIdx, &sequencePtr, &sequenceLength);
        if (sequencePtr == NULL) {
            std::cerr << "Error: could not read sequence from fasta at given sequence idx" << std::endl;
            return ExtensionResult(0, false);
        }

        const float extensionThresholdScore = PhmmProcessor::findThreshold(context.phmm,
            context.searchParams.extensionPValue, sequenceLength);

        float priorFlankScore, postFlankScore;
        if (context.isReverseComplement) {
            priorFlankScore = findFlankingDiag<true>(sequencePtr, sequenceLength, context, hitPosition, true);
            postFlankScore = findFlankingDiag<true>(sequencePtr, sequenceLength, context, hitPosition, false);
        }
        else {
            priorFlankScore = findFlankingDiag<false>(sequencePtr, sequenceLength, context, hitPosition, true);
            postFlankScore = findFlankingDiag<false>(sequencePtr, sequenceLength, context, hitPosition, false);
        }

        const float finalScore = priorFlankScore + postFlankScore + maxScoreAlongDiagonal;
        return ExtensionResult(finalScore, finalScore >= extensionThresholdScore);
    }
}
