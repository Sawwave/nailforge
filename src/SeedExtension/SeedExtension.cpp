#include "SeedExtension.hpp"
#include "../PhmmProcessor/PhmmProcessor.hpp"
#include "../Alphabet/LetterConversion.hpp"
#include "iostream"
#include <array>
#include <algorithm>
#include <immintrin.h>

namespace NailForge::SeedExtension {
    const __m256i aminoIdxLookup = _mm256_set_epi8(
        20, 0, 20, 1, 2, 3, 4, 5,
        6, 7, 20, 8, 9, 10, 11, 20,
        12, 13, 14, 15, 16, 20, 17, 18,
        20, 19, 20, 20, 20, 20, 20, 20);

    const __m256i dnaIdxLookup = _mm256_set_epi8(
        0, 0, 0, 1, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 3, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
    );
    // constexpr int32_t numDiagonals = 8;


    float findFlankingDiagMax(const char* sequencePtr, const uint64_t sequenceLength,
        const P7Hmm& phmm, const std::vector<float>& phmmMatchScores, const NailForge::SearchParams& searchParams,
        const uint64_t localSequencePosition, const uint32_t modelPosition, const bool isReverseCompliment, const bool isPriorFlank);

    inline bool canSimdExtendFlank(const uint64_t localSequencePosition, const uint32_t modelPosition,
        const uint64_t sequenceLength, const uint32_t modelLength, const uint8_t extensionLength,
        const bool isReverseCompliment, const bool isPriorFlank);


    bool verifySeedViaExtension(const uint64_t localSequencePosition, const uint32_t modelPosition,
        const float maxScoreAlongDiagonal, const uint64_t sequenceIdx, const NailForge::SearchParams& searchParams,
        const bool isReverseCompliment, const FastaVector& fastaVector, const P7Hmm& phmm,
        const std::vector<float>& phmmMatchScores) {

        char* sequencePtr;
        size_t sequenceLength;
        fastaVectorGetSequence(&fastaVector, sequenceIdx, &sequencePtr, &sequenceLength);
        if (sequencePtr == NULL) {
            std::cerr << "Error: could not read sequence from fasta at given sequence idx" << std::endl;
            return false;
        }

        const float extensionThreshold = PhmmProcessor::findThresholdForAmino(phmm, sequenceLength, searchParams.extensionPValue);

        const float priorDiagMax = findFlankingDiagMax(sequencePtr, sequenceLength,
            phmm, phmmMatchScores, searchParams, localSequencePosition, modelPosition, isReverseCompliment, true);
        const float postDiagMax = findFlankingDiagMax(sequencePtr, sequenceLength,
            phmm, phmmMatchScores, searchParams, localSequencePosition, modelPosition, isReverseCompliment, false);

        return (std::max(priorDiagMax, postDiagMax) + maxScoreAlongDiagonal) >= extensionThreshold;
    }


    float findFlankingDiagMax(const char* sequencePtr, const uint64_t sequenceLength,
        const P7Hmm& phmm, const std::vector<float>& phmmMatchScores, const NailForge::SearchParams& searchParams,
        const uint64_t localSequencePosition, const uint32_t modelPosition, const bool isReverseCompliment, const bool isPriorFlank) {

        const uint8_t symbolEncodingComplimentBitmask = isReverseCompliment ? 0x03 : 0x00;
        const NailForge::Alphabet alphabet = phmm.header.alphabet == P7HmmReaderAlphabetAmino ?
            NailForge::Alphabet::Amino : NailForge::Alphabet::Dna;
        const auto alphabetSize = NailForge::getAlphabetSize(alphabet);

        std::vector<float> diagonalMaxes;
        diagonalMaxes.resize(searchParams.extentionGroupWidth);
        std::fill(diagonalMaxes.begin(), diagonalMaxes.end(), 0.0f);

        for (int32_t diagonalIdx = 0; diagonalIdx < searchParams.extentionGroupWidth; diagonalIdx++) {
            float thisDiagScore = 0.0f;
            //half the diagons will be on an earlier diag, half on a leter diag
            const int32_t diagonalOffset = (int32_t)diagonalIdx - (searchParams.extentionGroupWidth / 2);
            int64_t modelPositionsRemaining;
            int64_t sequencePositionsRemaining;

            /*
            this is pretty complicated, but here's what's going on. We have to handle prior and post flanking regions,
            as well as reverse compliment search.

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
                    cellModelPosition = (int64_t)modelPosition + (int64_t)diagonalOffset;
                    modelPositionsRemaining = (int64_t)phmm.header.modelLength - (int64_t)cellModelPosition;
                }
                else {
                    cellModelPosition = (int64_t)modelPosition - (int64_t)diagonalOffset;
                    modelPositionsRemaining = cellModelPosition;
                }
                cellSequencePosition = (int64_t)localSequencePosition + (int64_t)diagonalOffset;
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
                cellSequencePosition = (int64_t)localSequencePosition + (int64_t)searchParams.maximumHitLength + (int64_t)diagonalOffset;
                sequencePositionsRemaining = (int64_t)sequenceLength - (int64_t)cellSequencePosition - 1;
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

                    const char sequenceSymbol = sequencePtr[cellSequencePosition];
                    const uint8_t symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);
                    const uint8_t phmmSymbolIdx = symbolIdx ^ symbolEncodingComplimentBitmask;
                    const float matchScore = phmmMatchScores[(cellModelPosition * alphabetSize) + phmmSymbolIdx];

                    thisDiagScore = std::max(thisDiagScore + matchScore, 0.0f);
                    diagonalMaxes[diagonalIdx] = std::max(diagonalMaxes[diagonalIdx], thisDiagScore);
                }
            }
            else {
                const auto cellModelPositionIncrement = isReverseCompliment ? -1 : 1;
                for (auto diagCellIdx = 0; diagCellIdx < diagLength; diagCellIdx++) {
                    //update the positions
                    cellSequencePosition++;
                    cellModelPosition += cellModelPositionIncrement;

                    const char sequenceSymbol = sequencePtr[cellSequencePosition];
                    const uint8_t symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);
                    const uint8_t phmmSymbolIdx = symbolIdx ^ symbolEncodingComplimentBitmask;
                    const float matchScore = phmmMatchScores[(cellModelPosition * alphabetSize) + phmmSymbolIdx];

                    thisDiagScore = std::max(thisDiagScore + matchScore, 0.0f);
                    diagonalMaxes[diagonalIdx] = std::max(diagonalMaxes[diagonalIdx], thisDiagScore);
                }
            }
        }

        const auto maxIter = std::max_element(diagonalMaxes.begin(), diagonalMaxes.end());
        return *maxIter;
    }

    // inline bool canSimdExtendFlank(const uint64_t localSequencePosition, const uint32_t modelPosition,
    //     const uint64_t sequenceLength, const uint32_t modelLength, const uint8_t extensionLength,
    //     const bool isReverseCompliment, const bool isPriorFlank) {
    //     const uint32_t requiredLength = extensionLength + numDiagonals / 2;
    //     uint32_t modelLengthAvailable;
    //     uint64_t sequenceLengthAvailable;
    //     //if 
    //     if (isPriorFlank) {
    //         sequenceLengthAvailable = localSequencePosition;
    //         modelLengthAvailable = isReverseCompliment ?
    //             modelLength - modelPosition : modelPosition;
    //     }
    //     else {
    //         sequenceLengthAvailable = sequenceLength - localSequencePosition;
    //         modelLengthAvailable = isReverseCompliment ?
    //             modelPosition : modelLength;
    //     }

    //     return modelLengthAvailable > requiredLength && sequenceLengthAvailable > (requiredLength + 4);
    // }

    // float findFlankingDiagMaxSimd(const char* sequencePtr, const uint64_t sequenceLength,
    //     const P7Hmm& phmm, const std::vector<float>& phmmMatchScores, const NailForge::SearchParams& searchParams,
    //     const uint64_t localSequencePosition, const uint32_t modelPosition, const bool isReverseCompliment, const bool isPriorFlank) {

    //     __m256 diagMaxes = _mm256_set1_ps(0.0f);
    //     __m256 diagonalScores = _mm256_set1_ps(0.0f);

    //     //find the lowest symbol to extract for this diagonal
    //     uint64_t lowestDiagSymbolIdx = localSequencePosition - 4;
    //     if (!isPriorFlank) {
    //         lowestDiagSymbolIdx += searchParams.flankExtensionLength;
    //     }

    //     const __m256i incrementVector = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);

    //     for (uint32_t diagCellIdx = 0; diagCellIdx < searchParams.flankExtensionLength; diagCellIdx++) {
    //         const __m256i symbolVector = _mm256_i32gather_epi32(sequencePtr + lowestDiagSymbolIdx, incrementVector, 1);
    //         const __m256i capitalBitmaskVector = _mm256_set1_epi32(0x1F);
    //         const __m256i capitalSymbolVector = _mm256_and_si256(symbolVector, capitalBitmaskVector);
    //         const __m256i symbolIdxVector = phmm.header.alphabet == P7HmmReaderAlphabetAmino ?
    //             _mm256_i32gather_epi32(aminoIdxLookup, incrementVector, 1) :
    //             _mm256_i32gather_epi32(aminoIdxLookup, incrementVector, 1);

    //         const __m256i modelOffsetVector = phmm.header.alphabet == P7HmmReaderAlphabetAmino ?
    //             _mm256_set_epi32(20 * 0, 20 * 1, 20 * 2, 20 * 3, 20 * 4, 20 * 5, 20 * 6, 20 * 7) :
    //             _mm256_set_epi32(4 * 0, 4 * 1, 4 * 2, 4 * 3, 4 * 4, 4 * 5, 4 * 6, 4 * 7);
    //     }
    // }
}
