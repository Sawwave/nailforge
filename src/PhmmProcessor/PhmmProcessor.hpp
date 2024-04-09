#ifndef NAIL_Forge_PHMM_PROCESSOR_HPPvalues
#define NAIL_Forge_PHMM_PROCESSOR_HPP
#include <vector>
#include <cstdint>
#include "../nailforge.hpp"

namespace NailForge::PhmmProcessor {
    //Converts the phmm match scores to bits represented by floats
    std::vector<float> toFloatMatchScores(const P7Hmm& phmm);

    //finds a phmm's threshold value for a given p-value using the model's MAXL field
    float findThresholdForNucleotide(const P7Hmm& phmm, const float pValue);

    //find a phmm's threshold value for a given p-value using the target sequence length
    float findThresholdForAmino(const P7Hmm& phmm, const uint64_t sequenceLength, const float pValue);

    //The inverse gumbel survival function, aka, what score is needed to hit the given p-value?
    float gumbelInverseSurvival(const float lambda, const float mu, const float pValue);


    //finds the number of match scores per phmm position, aka, the cardinality of phmm's alphabet.
    uint64_t getNumMatchScores(const P7Hmm& phmm);

};

#endif