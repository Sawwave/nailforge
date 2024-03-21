#include "SeedExtension.hpp"
#include "../Alphabet/LetterConversion.hpp"
#include "iostream"
#include <array>
#include <algorithm>

namespace NailForge::SeedExtension {


    constexpr int32_t numDiagonals = 8;


    float findFlankingDiagMax(const char* sequencePtr, const uint64_t sequenceLength,
        const P7Hmm& phmm, const std::vector<float>& phmmMatchScores, const NailForge::SearchParams& searchParams,
        const uint64_t localSequencePosition, const uint32_t modelPosition, const bool isReverseCompliment, const bool isPriorFlank);

    bool verifySeedViaExtension(const uint64_t localSequencePosition, const uint32_t modelPosition,
        const uint64_t sequenceIdx, const NailForge::SearchParams& searchParams, const bool isReverseCompliment,
        const FastaVector& fastaVector, const P7Hmm& phmm, const std::vector<float>& phmmMatchScores) {
        char* sequencePtr;
        size_t sequenceLength;
        fastaVectorGetSequence(&fastaVector, sequenceIdx, &sequencePtr, &sequenceLength);
        if (sequencePtr == NULL) {
            std::cerr << "Error: could not read sequence from fasta at given sequence idx" << std::endl;
            return false;
        }

        float diagMax = findFlankingDiagMax(sequencePtr, sequenceLength,
            phmm, phmmMatchScores, searchParams, localSequencePosition, modelPosition, isReverseCompliment, true);

        if (__builtin_expect(diagMax >= searchParams.extensionThresholdScore, false)) {
            return true;
        }
        diagMax += findFlankingDiagMax(sequencePtr, sequenceLength,
            phmm, phmmMatchScores, searchParams, localSequencePosition, modelPosition, isReverseCompliment, false);
        return diagMax >= searchParams.extensionThresholdScore;
    }


    float findFlankingDiagMax(const char* sequencePtr, const uint64_t sequenceLength,
        const P7Hmm& phmm, const std::vector<float>& phmmMatchScores, const NailForge::SearchParams& searchParams,
        const uint64_t localSequencePosition, const uint32_t modelPosition, const bool isReverseCompliment, const bool isPriorFlank) {

        const uint8_t symbolEncodingComplimentBitmask = isReverseCompliment ? 0x03 : 0x00;
        const NailForge::Alphabet alphabet = phmm.header.alphabet == P7HmmReaderAlphabetAmino ?
            NailForge::Alphabet::Amino : NailForge::Alphabet::Dna;
        const auto alphabetSize = NailForge::getAlphabetSize(alphabet);

        std::array<float, numDiagonals> diagonalMaxes;
        std::fill(diagonalMaxes.begin(), diagonalMaxes.end(), 0.0f);

        for (int32_t diagonalIdx = 0; diagonalIdx < numDiagonals; diagonalIdx++) {
            float thisDiagScore = 0.0f;
            //half the diagons will be on an earlier diag, half on a leter diag
            const int32_t diagonalOffset = (int32_t)diagonalIdx - (numDiagonals / 2);
            int64_t modelPositionsRemaining;
            int64_t sequencePositionsRemaining;

            /*
            this is pretty complicated, but here's what's going on. We have to handle prior and post flanking regions,
            as well as reverse compliment search.  For standard search, the flanks look like the following

                    sequence                                sequence
             ---------------------------------        ---------------------------------
           m |           \                   |      m |                               |
           o |            \  n-1             |      o |          0  /                 |
           d |             \                 |      d |            /                  |
           e |       \      \                |      e |           /        /          |
           l |     0  \   *                  |      l |          /        /           |
             |         \   \                 |        |                  /  n-1       |
             |              \                |        |             /   /             |
             |               \  \            |        |            /                  |
             |                \  \           |        |      /    /                   |
             |             \      \          |        |     /    *                    |
             |              \      \ n-1     |        |  0 /                          |
             |               \               |        |   /        /                  |
             |           0    \              |        |           /  n-1              |
             |                               |        |          /                    |
             ---------------------------------        ---------------------------------
            Here, the star is the reported position. The paddles at the ends of the main diagonal are the regions we are going to explore.
            The left diagram shows the standard alignment, and the right shows search on the reverse compliment strand.
            */

            int64_t cellModelPosition;
            int64_t cellSequencePosition;
            if (isPriorFlank) {
                if (isReverseCompliment) {
                    cellModelPosition = modelPosition + diagonalOffset;
                    modelPositionsRemaining = phmm.header.modelLength - cellModelPosition;
                }
                else {
                    cellModelPosition = (int64_t)modelPosition - (int64_t)diagonalOffset;
                    modelPositionsRemaining = cellModelPosition;
                }
                cellSequencePosition = localSequencePosition + diagonalOffset;
                sequencePositionsRemaining = cellSequencePosition;
            }
            else {
                if (isReverseCompliment) {
                    cellModelPosition = (int64_t)modelPosition - (int64_t)searchParams.maximumHitLength + (int64_t)diagonalOffset;
                    modelPositionsRemaining = cellModelPosition;
                }
                else {
                    cellModelPosition = (int64_t)modelPosition + (int64_t)searchParams.maximumHitLength - (int64_t)diagonalOffset;
                    modelPositionsRemaining = phmm.header.modelLength - cellModelPosition - 1;
                }
                cellSequencePosition = localSequencePosition + searchParams.maximumHitLength + diagonalOffset;
                sequencePositionsRemaining = sequenceLength - cellSequencePosition - 1;
            }


            modelPositionsRemaining = std::max(modelPositionsRemaining, (int64_t)0);
            modelPositionsRemaining = std::min(modelPositionsRemaining, (int64_t)searchParams.flankExtensionLength);
            sequencePositionsRemaining = std::max(sequencePositionsRemaining, (int64_t)0);
            sequencePositionsRemaining = std::min(sequencePositionsRemaining, (int64_t)searchParams.flankExtensionLength);

            //compute the score along this diagonal
            const int64_t diagLength = std::min((int64_t)modelPositionsRemaining, sequencePositionsRemaining);
            if (isPriorFlank) {
                const auto cellModelPositionIncrement = isReverseCompliment ? 1 : -1;
                for (auto diagCellIdx = 0; diagCellIdx < diagLength; diagCellIdx++) {
                    //update the positions
                    cellSequencePosition--;
                    cellModelPosition += cellModelPositionIncrement;

                    const char sequenceSymbol = sequencePtr[cellSequencePosition] ^ symbolEncodingComplimentBitmask;
                    const uint8_t symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);
                    const float matchScore = phmmMatchScores[(cellModelPosition * alphabetSize) + symbolIdx];

                    thisDiagScore = std::max(thisDiagScore + matchScore, 0.0f);
                    diagonalMaxes[diagonalIdx] = std::max(diagonalMaxes[diagonalIdx], thisDiagScore);

                }
            }
            else {
                const auto cellModelPositionIncrement = isReverseCompliment ? -1 : 1;
                for (auto diagCellIdx = 0; diagCellIdx < diagLength; diagCellIdx++) {
                    //update the positions
                    cellSequencePosition++;
                    cellModelPosition += isReverseCompliment ? -1 : 1;

                    const char sequenceSymbol = sequencePtr[cellSequencePosition] ^ symbolEncodingComplimentBitmask;
                    const uint8_t symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);
                    const float matchScore = phmmMatchScores[(cellModelPosition * alphabetSize) + symbolIdx];

                    thisDiagScore = std::max(thisDiagScore + matchScore, 0.0f);
                    diagonalMaxes[diagonalIdx] = std::max(diagonalMaxes[diagonalIdx], thisDiagScore);

                }
            }

            // for (auto diagCellIdx = 0; diagCellIdx < diagLength; diagCellIdx++) {
            //     float matchScore;
            //     uint8_t symbolIdx;
            //     uint32_t  cellModelPosition;

            //     if (isPriorFlank) {
            //         const auto cellSequencePosition = (localSequencePosition - diagonalOffset) - diagCellIdx;
            //         const char sequenceSymbol = sequencePtr[cellSequencePosition] ^ symbolEncodingComplimentBitmask;
            //         symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);

            //         cellModelPosition = isReverseCompliment ?
            //             modelPosition - diagonalOffset + diagCellIdx :
            //             modelPosition - diagonalOffset - diagCellIdx;

            //     }
            //     else {
            //         const auto cellSequencePosition = localSequencePosition + searchParams.flankExtensionLength + diagCellIdx + diagonalOffset;
            //         const char sequenceSymbol = sequencePtr[cellSequencePosition] ^ symbolEncodingComplimentBitmask;
            //         symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);
            //         cellModelPosition = isReverseCompliment ?
            //             modelPosition - searchParams.flankExtensionLength + diagonalOffset - diagCellIdx :
            //             modelPosition + searchParams.flankExtensionLength - diagonalOffset + diagCellIdx;
            //     }
            //     matchScore = phmmMatchScores[(cellModelPosition * alphabetSize) + symbolIdx];

            //     thisDiagScore = std::max(thisDiagScore + matchScore, 0.0f);

            //     diagonalMaxes[diagonalIdx] = std::max(diagonalMaxes[diagonalIdx], thisDiagScore);
            // }
        }

        const auto maxIter = std::max_element(diagonalMaxes.begin(), diagonalMaxes.end());
        return *maxIter;
    }

}