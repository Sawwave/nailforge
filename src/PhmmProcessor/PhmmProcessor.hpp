#ifndef NAIL_Forge_PHMM_PROCESSOR_HPPvalues
#define NAIL_Forge_PHMM_PROCESSOR_HPP
#include <vector>
#include <cstdint>
#include "../nailforge.hpp"

namespace NailForge::PhmmProcessor {
    //Converts the phmm match scores into bits represented by 8-bit values
    std::vector<int8_t> to8BitMatchScores(const P7Hmm& phmm);

    //Converts the phmm match scores into projected 8-bit values, where a score of 'thresholdHitValue' represents a threshold hit
    std::vector<int8_t> toProjected8BitMatchScores(const P7Hmm& phmm, const float pValue, const float thresholdHitValue);

    //Converts the phmm match scores to bits represented by floats
    std::vector<float> toFloatMatchScores(const P7Hmm& phmm);

    //currently unimplemented, might implement later
    //  std::vector<int16_t> to16BitMatchScores(const P7Hmm* phmm);
    //  std::vector<int16_t> toProjected16BitMatchScores(const P7Hmm* phmm, const float pValue);
    //  std::vector<int32_t> to32BitMatchScores(const P7Hmm* phmm);
    //  std::vector<int32_t> toProjected32BitMatchScores(const P7Hmm* phmm, const float pValue);

    //finds a phmm's threshold value for a given p-value
    float findThreshold(const P7Hmm& phmm, const float pValue);

    //The inverse gumbel survival function, aka, what score is needed to hit the given p-value?
    float gumbelInverseSurvival(const float lambda, const float mu, const float pValue);

    //finds the tau scaling factor, as documented in https://www.biorxiv.org/content/10.1101/2023.09.20.558701v1
    float findTauScalingFactor(const P7Hmm& phmm, const float pValue, const float projectedThreshold);

    //finds the number of match scores per phmm position, aka, the cardinality of phmm's alphabet.
    uint64_t getNumMatchScores(const P7Hmm& phmm);

};

#endif