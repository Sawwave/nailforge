#ifndef NAIL_FORGE_SEED_EXTENSION_HPP
#define NAIL_FORGE_SEED_EXTENSION_HPP

#include "../nailforge.hpp"
#include "../Alphabet/LetterConversion.hpp"
#include "../StringTree/StringTree.hpp"
#include <iostream>

namespace NailForge::SeedExtension {



    //lsp, mp, params, phmm, matchScores, hitLen
    bool verifySeedViaExtension(const StringTree::Context& context,
        const StringTree::HitPosition, const float maxScoreAlongDiagonal);

    template<bool isReverseCompliment>
    float findFlankingDiag(const char* sequencePtr, const uint64_t sequenceLength, const StringTree::Context& context,
        const StringTree::HitPosition& hitPosition, const bool isPriorFlank) noexcept;


    //hitLength should be currentDepth + 1
    template<bool isReverseCompliment>
    float findFlankingDiag(const char* sequencePtr, const uint64_t sequenceLength, const StringTree::Context& context,
        const StringTree::HitPosition& hitPosition, const bool isPriorFlank) noexcept {

        const uint8_t symbolEncodingComplimentBitmask = isReverseCompliment ? 0x03 : 0x00;
        const NailForge::Alphabet alphabet = context.phmm.header.alphabet == P7HmmReaderAlphabetAmino ?
            NailForge::Alphabet::Amino : NailForge::Alphabet::Dna;
        const auto alphabetSize = NailForge::getAlphabetSize(alphabet);

        float accumulatedScore = 0.0f;
        float maxAccumulatedScore = 0.0f;
        uint64_t flankLength;

        if constexpr (isReverseCompliment) {
            flankLength = isPriorFlank ?
                std::min(hitPosition.localSequencePosition,
                    static_cast<uint64_t>(context.phmm.header.modelLength - hitPosition.modelPosition - 1)) :
                std::min(sequenceLength - hitPosition.localSequencePosition - hitPosition.hitLength,
                    static_cast<uint64_t>(hitPosition.modelPosition - hitPosition.hitLength));
        }
        else {
            flankLength = isPriorFlank ?
                std::min(hitPosition.localSequencePosition, static_cast<uint64_t>(hitPosition.modelPosition)) :
                std::min(sequenceLength - hitPosition.localSequencePosition - hitPosition.hitLength,
                    static_cast<uint64_t>(context.phmm.header.modelLength - hitPosition.modelPosition - hitPosition.hitLength));

        }
        flankLength = std::min(flankLength, static_cast<uint64_t>(context.searchParams.flankExtensionLength));

        for (uint32_t flankPosition = 0; flankPosition < flankLength; flankPosition++) {
            uint32_t flankModelPosition;
            uint32_t flankSequencePosition;

            flankSequencePosition = isPriorFlank ?
                hitPosition.localSequencePosition - flankPosition - 1 :
                hitPosition.localSequencePosition + flankPosition + hitPosition.hitLength;

            if constexpr (isReverseCompliment) {
                flankModelPosition = isPriorFlank ?
                    hitPosition.modelPosition + flankPosition + 1 :
                    hitPosition.modelPosition - flankPosition - hitPosition.hitLength;
            }
            else {
                flankModelPosition = isPriorFlank ?
                    hitPosition.modelPosition - flankPosition - 1 :
                    hitPosition.modelPosition + flankPosition + hitPosition.hitLength;
            }

            const char sequenceSymbol = sequencePtr[flankSequencePosition];
            const uint8_t symbolIdx = NailForge::LetterConversion::asciiLetterToLetterIndex(sequenceSymbol, alphabet);
            const uint8_t phmmSymbolIdx = symbolIdx ^ symbolEncodingComplimentBitmask;
            if (flankModelPosition >= context.phmm.header.modelLength) {
                std::cout << "\nflank mod pos " << flankModelPosition << "greater than len " << context.phmm.header.modelLength << std::endl;
                std::cout << "iscomp " << (int)isReverseCompliment << ", isprior? " << (int)isPriorFlank << std::endl;
                std::cout << "modpos" << hitPosition.modelPosition << ", flankPosition " << flankPosition << std::endl;
                std::cout << "hitlen " << (int)hitPosition.hitLength << std::endl;
            }
            const float matchScore = context.matchScores[(flankModelPosition * alphabetSize) + phmmSymbolIdx];

            accumulatedScore += matchScore;
            maxAccumulatedScore = std::max(maxAccumulatedScore, accumulatedScore);
        }
        return maxAccumulatedScore;
    }

}

#endif