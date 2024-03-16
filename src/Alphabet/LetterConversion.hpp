#ifndef NAIL_FILTER_LETTER_CONVERSION_HPP
#define NAIL_FILTER_LETTER_CONVERSION_HPP
#include <cstdint>
#include <string>
#include <vector>
#include "Alphabet.hpp"
extern "C" {
#include "AwFmIndex.h"
}

namespace NailForge::LetterConversion {
    char letterIndexToAscii(const uint8_t letterIndex,
        const NailForge::Alphabet alphabet) noexcept;

    uint8_t asciiLetterToLetterIndex(const char asciiLetter,
        const NailForge::Alphabet alphabet) noexcept;
    char letterIndexToComplimentAscii(const uint8_t letterIndex) noexcept;

    char reverseComplimentAscii(const char letter) noexcept;

    void reverseComplimentCstr(char* cString, const uint32_t stringLength)noexcept;

    std::string letterIndexVectorToString(const std::vector<uint8_t> letterIndexVector, const NailForge::Alphabet alphabet)noexcept;

    template <uint8_t arrayLen>
    void arrayToCstr(char* to, std::array<uint8_t, arrayLen> fromArray)noexcept;

    template <uint8_t arrayLen>
    void arrayToReverseComplimentCstr(char* to,
        std::array<uint8_t, arrayLen> fromArray)noexcept;
}

#endif