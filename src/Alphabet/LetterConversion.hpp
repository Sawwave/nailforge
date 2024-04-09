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
    /// @brief converts a letter index, either 0-3 or 0-19 to ascii representation
    /// @param letterIndex index between 0 and alphabet cardinality -1
    /// @param alphabet alphabet that the letter index is from, AA, RNA, or DNA.
    /// @return lowercase ascii for the given symbol. 
    char letterIndexToAscii(const uint8_t letterIndex,
        const NailForge::Alphabet alphabet) noexcept;

    /// @brief converts the given ascii letter to the symbol index from the given alphabet
    /// @param asciiLetter ascii of the symbol, either capital or lowercase
    /// @param alphabet alphabe the symbol is from
    /// @return symbol index, 0-19 for AA, 0-3 for DNA/RNA
    uint8_t asciiLetterToLetterIndex(const char asciiLetter,
        const NailForge::Alphabet alphabet) noexcept;

    /// @brief finds the compliment of the given DNA/RNA letter, and converts it to ascii
    /// @param letterIndex letter to compliment, should be in range 0-3
    /// @return ascii for the letter's compliment.
    char letterIndexToComplimentAscii(const uint8_t letterIndex) noexcept;

    /// @brief finds the compliment for the ascii letter
    /// @param letter DNA/RNA character in ascii
    /// @return compliment of the given character, in ascii
    char reverseComplimentAscii(const char letter) noexcept;

    /// @brief performs an in-place conversion of the given C-style string to its reverse compliment
    /// @param cString start of the C-style string to reverse compliment
    /// @param stringLength length of the given string
    void reverseComplimentCstr(char* cString, const uint32_t stringLength)noexcept;

}

#endif