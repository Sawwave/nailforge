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

    bool alphabetsMatch(const P7HmmList& phmmList, const AwFmIndex* fmIndex);
    uint32_t getAlphabetSize(const NailForge::Alphabet alphabet);

}

#endif