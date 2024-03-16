#include "Alphabet.hpp"

namespace NailForge {

    bool alphabetsMatch(const P7HmmList& phmmList, const AwFmIndex* fmIndex) {
        if (phmmList.count == 0) {
            return true;
        }
        else {
            //it's okay if one is rna and the other is dna, but you can't combine amino and either of the others together.
            P7Alphabet firstPhmmAlphabet = phmmList.phmms[0].header.alphabet;
            AwFmAlphabetType fmAlphabet = fmIndex->config.alphabetType;


            //check to make sure all models in the hmm file have the same alphabet
            for (uint32_t i = 1; i < phmmList.count; i++) {
                if (phmmList.phmms[i].header.alphabet != firstPhmmAlphabet) {
                    return false;
                }
            }
            switch (firstPhmmAlphabet) {
            case P7HmmReaderAlphabetAmino:
                return fmAlphabet == AwFmAlphabetAmino;
            case P7HmmReaderAlphabetDna:    //fall-through case
            case P7HmmReaderAlphabetRna:
                return fmAlphabet == AwFmAlphabetDna || fmAlphabet == AwFmAlphabetRna;
            default: return false;
            }
        }
    }

    uint32_t getAlphabetSize(const NailForge::Alphabet alphabet) {
        if (alphabet == Alphabet::Amino) {
            return 20;
        }
        else {
            return 4;
        }
    }
}