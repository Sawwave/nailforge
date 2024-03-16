#include "LetterConversion.hpp"
#include <ctype.h>


namespace NailForge::LetterConversion {

    char letterIndexToAscii(const uint8_t letterIndex, const NailForge::Alphabet alphabet) noexcept {
        if (alphabet == NailForge::Alphabet::Dna) {
            static constexpr char letters[4] = { 'a', 'c', 'g', 't' };
            return letters[letterIndex];

        }
        else if (alphabet == NailForge::Alphabet::Rna) {
            static constexpr char letters[4] = { 'a', 'c', 'g', 'u' };
            return letters[letterIndex];
        }
        else {
            static constexpr char letters[20] = { 'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l',
                                                'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y' };
            return letters[letterIndex];
        }
    }

    uint8_t asciiLetterToLetterIndex(const char asciiLetter, const NailForge::Alphabet alphabet) noexcept {
        const char lowercaseLetter = asciiLetter | 0x20;
        if (alphabet == NailForge::Alphabet::Amino) {
            switch (lowercaseLetter) {
            case 'a': return 0;
            case 'c': return 1;
            case 'd': return 2;
            case 'e': return 3;
            case 'f': return 4;
            case 'g': return 5;
            case 'h': return 6;
            case 'i': return 7;
            case 'k': return 8;
            case 'l': return 9;
            case 'm': return 10;
            case 'n': return 11;
            case 'p': return 12;
            case 'q': return 13;
            case 'r': return 14;
            case 's': return 15;
            case 't': return 16;
            case 'v': return 17;
            case 'w': return 18;
            case 'y': return 19;
            default: return 20;
            }
        }
        else {
            switch (lowercaseLetter) {
            case 'a': return 0;
            case 'c': return 1;
            case 'g': return 2;
            default: return 3;
            }
        }
    }

    char letterIndexToComplimentAscii(const uint8_t letterIndex) noexcept {
        static constexpr char letters[4] = { 't', 'g', 'c', 'a' };
        return letters[letterIndex];
    }


    char reverseComplimentAscii(const char letter) noexcept {
        switch (letter) {
        case 'a': case 'A': return 't';
        case 'c': case 'C': return 'g';
        case 'g': case 'G': return 'c';
        default: return 'a';
        }
    }


    void reverseComplimentCstr(char* cString, const uint32_t stringLength) noexcept {
        for (uint32_t letterIdx = 0; letterIdx <= (stringLength / 2); letterIdx++) {

            char* firstChar = &cString[letterIdx];
            char* secondChar = &cString[stringLength - 1 - letterIdx];
            switch (*firstChar & 0x7) {   // agct can be resolved using only the first 3 bits, no matter capitalization
            case 0x001: *firstChar = 't';   break;
            case 0x011: *firstChar = 'g';   break;
            case 0x111: *firstChar = 'c';   break;
            default: *firstChar = 'a';
            }
            switch (*secondChar & 0x7) {   // agct can be resolved using only the first 3 bits, no matter capitalization
            case 0x001: *secondChar = 't';  break;
            case 0x011: *secondChar = 'g';  break;
            case 0x111: *secondChar = 'c';  break;
            default: *secondChar = 'a';
            }
            std::swap(*firstChar, *secondChar);
        }
    }


    std::string letterIndexVectorToString(const std::vector<uint8_t> letterIndexVector, const NailForge::Alphabet alphabet) noexcept {
        std::string s;
        for (const auto& letterIndex : letterIndexVector) {
            s.push_back(letterIndexToAscii(letterIndex, alphabet));
        }
        return s;
    }

    template <uint8_t arrayLen>
    void rrayToCstr(char* to, std::array<uint8_t, arrayLen> fromArray) noexcept {
        for (uint8_t i = 0; i < arrayLen; i++) {
            to[i] = letterIndexToAscii(fromArray[i]);
        }
    }

    template <uint8_t arrayLen>
    void arrayToReverseComplimentCstr(char* to,
        std::array<uint8_t, arrayLen> fromArray) noexcept {
        for (uint8_t i = 0; i < arrayLen; i++) {
            to[i] = letterIndexToComplimentAscii(fromArray[arrayLen - i - 1]);
        }
    }

}
