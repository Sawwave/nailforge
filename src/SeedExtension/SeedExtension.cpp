#include "SeedExtension.hpp"
#include "../Alphabet/LetterConversion.hpp"
#include "iostream"
#include <array>
#include <algorithm>

namespace NailForge::SeedExtension {


    const int32_t numDiagonals = 9;


    float findFlankingDiagMax(const char* sequencePtr, const uint64_t sequenceLength,
        const P7Hmm& phmm, const NailForge::SearchParams& searchParams, const uint64_t localSequencePosition,
        const uint32_t modelPosition, const bool isReverseCompliment, const bool isPriorFlank);

    bool verifySeedViaExtension(const uint64_t localSequencePosition, const uint32_t modelPosition,
        const uint64_t sequenceIdx, const NailForge::SearchParams& searchParams, const bool isReverseCompliment,
        const FastaVector& fastaVector, const P7Hmm& phmm) {
        char* sequencePtr;
        size_t sequenceLength;
        fastaVectorGetSequence(&fastaVector, sequenceIdx, &sequencePtr, &sequenceLength);
        if (sequencePtr == NULL) {
            std::cerr << "Error: could not read sequence from fasta at given sequence idx" << std::endl;
            return false;
        }

        float diagMax = findFlankingDiagMax(sequencePtr, sequenceLength,
            phmm, searchParams, modelPosition, localSequencePosition, isReverseCompliment, true);

        if (__builtin_expect(diagMax >= searchParams.extensionThresholdScore, false)) {
            return true;
        }
        diagMax += findFlankingDiagMax(sequencePtr, sequenceLength,
            phmm, searchParams, modelPosition, localSequencePosition, isReverseCompliment, false);
        return diagMax >= searchParams.extensionThresholdScore;
    }


    float findFlankingDiagMax(const char* sequencePtr, const uint64_t sequenceLength,
        const P7Hmm& phmm, const NailForge::SearchParams& searchParams, const uint64_t localSequencePosition,
        const uint32_t modelPosition, const bool isReverseCompliment, const bool isPriorFlank) {


        const uint8_t symbolEncodingComplimentBitmask = isReverseCompliment ? 0x03 : 0x00;
        const NailForge::Alphabet alphabet = phmm.header.alphabet == P7HmmReaderAlphabetAmino ?
            NailForge::Alphabet::Amino : NailForge::Alphabet::Dna;
        const auto alphabetSize = NailForge::getAlphabetSize(alphabet);

        std::array<float, 8> diagonalMaxes;    //will be initialized to all zeros

        for (int32_t diagonalIdx = 0; diagonalIdx < numDiagonals; diagonalIdx++) {
            //half the diagons will be on an earlier diag, half on a leter diag
            const int32_t diagonalOffset = (int32_t)diagonalIdx - (numDiagonals / 2);
            int32_t modelPositionsRemaining;
            int64_t sequencePositionsRemaining;

            /*
            this is pretty complicated, but here's what's going on. We have to handle prior and post flanking regions,
            as well as reverse compliment search.  For standard search, the flanks look like the following

                    sequence                                sequence
             ---------------------------------        ---------------------------------
           m |           \                   |      m |                               |
           o |            \                  |      o |             /                 |
           d |             \                 |      d |            /                  |
           e |       \      \                |      e |           /        /          |
           l |        \   *                  |      l |          /        /           |
             |         \   \                 |        |                  /            |
             |              \                |        |             /   /             |
             |               \  \            |        |            /                  |
             |                \  \           |        |           /                   |
             |             \      \          |        |     /    *                    |
             |              \      \         |        |    /                          |
             |               \               |        |   /        /                  |
             |                \              |        |           /                   |
             |                               |        |          /                    |
             ---------------------------------        ---------------------------------
            Here, the star is the reported position. The paddles at the ends of the main diagonal are the regions we are going to explore.
            The left diagram shows the standard alignment, and the right shows search on the reverse compliment strand.
            */
            if (isPriorFlank) {
                if (isReverseCompliment) {
                    modelPositionsRemaining = phmm.header.modelLength - modelPosition - diagonalOffset;
                }
                else {
                    modelPositionsRemaining = modelPosition - diagonalOffset;
                }

                sequencePositionsRemaining = localSequencePosition + diagonalOffset;
            }
            else {
                if (isReverseCompliment) {
                    modelPositionsRemaining = modelPosition - searchParams.maximumHitLength + diagonalOffset;
                }
                else {
                    modelPositionsRemaining = phmm.header.modelLength - modelPosition -
                        searchParams.maximumHitLength + diagonalOffset;
                }
                sequencePositionsRemaining = sequenceLength - localSequencePosition -
                    searchParams.maximumHitLength - diagonalOffset;
            }

            modelPositionsRemaining = std::max(modelPositionsRemaining, (int32_t)0);
            modelPositionsRemaining = std::min(modelPositionsRemaining, (int32_t)searchParams.flankExtensionLength);
            sequencePositionsRemaining = std::max(sequencePositionsRemaining, (int64_t)0);
            sequencePositionsRemaining = std::min(sequencePositionsRemaining, (int64_t)searchParams.flankExtensionLength);

            const int64_t diagLength = std::min((int64_t)modelPositionsRemaining, sequencePositionsRemaining);
            for (auto diagCellIdx = 0; diagCellIdx < diagLength; diagCellIdx++) {
                float matchScore;
                uint8_t symbolIdx;
                uint32_t  cellModelPosition;

                if (isPriorFlank) {
                    const auto cellSequencePosition = (localSequencePosition - diagonalIdx) - diagCellIdx;
                    const char sequenceSymbol = sequencePtr[cellSequencePosition];
                    symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);

                    cellModelPosition = isReverseCompliment ?
                        modelPosition - diagonalIdx + diagCellIdx :
                        modelPosition - diagonalIdx - diagCellIdx;

                }
                else {
                    const auto cellSequencePosition = localSequencePosition + searchParams.flankExtensionLength + diagCellIdx + diagonalIdx;
                    const char sequenceSymbol = sequencePtr[cellSequencePosition];
                    symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);
                    cellModelPosition = isReverseCompliment ?
                        modelPosition - searchParams.flankExtensionLength + diagonalIdx - diagCellIdx :
                        modelPosition + searchParams.flankExtensionLength - diagonalIdx + diagCellIdx;
                }
                matchScore = phmm.model.matchEmissionScores[(cellModelPosition * alphabetSize) + symbolIdx];

                diagonalMaxes[diagCellIdx] += std::max(diagonalMaxes[diagCellIdx] + matchScore, 0.0f);
            }


            return *std::max(diagonalMaxes.begin(), diagonalMaxes.end());
        }

        const auto maxIter = std::max(diagonalMaxes.begin(), diagonalMaxes.end());
        return *maxIter;
    }

}