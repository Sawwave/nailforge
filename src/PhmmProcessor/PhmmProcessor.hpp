#ifndef NAIL_Forge_PHMM_PROCESSOR_HPP
#define NAIL_Forge_PHMM_PROCESSOR_HPP
#include <vector>
#include <cstdint>
#include "../nailforge.hpp"

class PhmmProcessor {
public:
//Converts the phmm match scores into bits represented by 8-bit values
    static std::vector<int8_t> to8BitMatchScores(const P7Hmm& phmm);

    //Converts the phmm match scores into projected 8-bit values, where a score of 'thresholdHitValue' represents a threshold hit
    static std::vector<int8_t> toProjected8BitMatchScores(const P7Hmm& phmm, const float pValue, const float thresholdHitValue);

    //Converts the phmm match scores to bits represented by floats
    static std::vector<float> toFloatMatchScores(const P7Hmm& phmm);

    //currently unimplemented, might implement later
    // static std::vector<int16_t> to16BitMatchScores(const P7Hmm* phmm);
    // static std::vector<int16_t> toProjected16BitMatchScores(const P7Hmm* phmm, const float pValue);
    // static std::vector<int32_t> to32BitMatchScores(const P7Hmm* phmm);
    // static std::vector<int32_t> toProjected32BitMatchScores(const P7Hmm* phmm, const float pValue);

    //finds a phmm's threshold value for a given p-value
    static float findThreshold(const P7Hmm& phmm, const float pValue);

private:
    //The inverse gumbel survival function, aka, what score is needed to hit the given p-value?
    static float gumbelInverseSurvival(const float lambda, const float mu, const float pValue);

    //finds the tau scaling factor, as documented in https://www.biorxiv.org/content/10.1101/2023.09.20.558701v1
    static float findTauScalingFactor(const P7Hmm& phmm, const float pValue, const float projectedThreshold);

    //finds the number of match scores per phmm position, aka, the cardinality of phmm's alphabet.
    static uint64_t getNumMatchScores(const P7Hmm& phmm);

};

#endif