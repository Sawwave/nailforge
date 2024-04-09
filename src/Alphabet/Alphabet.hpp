#ifndef NAIL_FILTER_ALPHABET_HPP
#define NAIL_FILTER_ALPHABET_HPP

extern "C" {
#include "AwFmIndex.h"
#include "p7HmmReader.h"
}

namespace NailForge {

    enum class Alphabet {
        Dna = 0,
        Rna = 1,
        Amino = 2
    };

    /// @brief determines if the alphabet in the phmm list matches the fm index. Only checks the first model of the phmm list.
    /// @param phmmList list of phmms
    /// @param fmIndex fm index
    /// @return true if the alphabets match 
    bool alphabetsMatch(const P7HmmList& phmmList, const AwFmIndex* fmIndex);

    /// @brief finds the cardinality of the given alphabet
    /// @param alphabet  alphabet to check
    /// @return 20 for amino acid, or 4 for DNA/RNA
    uint32_t getAlphabetSize(const NailForge::Alphabet alphabet);

}

#endif