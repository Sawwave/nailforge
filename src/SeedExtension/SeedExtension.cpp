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

        const float thresholdScore = PhmmProcessor::findThreshold(context.phmm,
            context.searchParams.extensionPValue, sequenceLength);

        //The scoreLengthTuples also contain the flank lengths, but they aren't being used. Future developers
        //might just change this back to returning the score from findFlankingDiag()
        ScoreLengthTuple priorFlankScoreLength, postFlankScoreLength;
        if (context.isReverseComplement) {
            priorFlankScoreLength = findFlankingDiag<true>(sequencePtr, sequenceLength, context, hitPosition, thresholdScore, true);
            postFlankScoreLength = findFlankingDiag<true>(sequencePtr, sequenceLength, context, hitPosition, thresholdScore, false);
        }
        else {
            priorFlankScoreLength = findFlankingDiag<false>(sequencePtr, sequenceLength, context, hitPosition, thresholdScore, true);
            postFlankScoreLength = findFlankingDiag<false>(sequencePtr, sequenceLength, context, hitPosition, thresholdScore, false);
        }

        const float finalScore = priorFlankScoreLength.score + postFlankScoreLength.score + maxScoreAlongDiagonal;

        return ExtensionResult(finalScore, finalScore >= thresholdScore);
    }
}
