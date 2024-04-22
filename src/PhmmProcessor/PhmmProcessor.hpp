#ifndef NAIL_Forge_PHMM_PROCESSOR_HPPvalues
#define NAIL_Forge_PHMM_PROCESSOR_HPP
#include <vector>
#include <cstdint>
#include "../nailforge.hpp"

namespace NailForge::PhmmProcessor {

    /// @brief converts the negative nat log likelihood scores in the phmm into bits scores for use in alignment.
    /// @param phmm phmm to convert score
    /// @return flattened float vector of scores. Score for symbol s at position p is located at vector idx [p * alphabetCardinality + s]
    std::vector<float> toFloatMatchScores(const P7Hmm& phmm);

    /// @brief finds the threshold score necessary to meet a given p-value. Will use the hmm MAXL property if it exists.
    /// @param phmm phmm to find the threshold for
    /// @param sequenceLength length of the target sequence. this will be used if the MAXL property is missing.
    /// @param pValue p-value that will be hit at the returned score.
    /// @return score, in bits, necessary to hit the given p-value.
    float findThreshold(const P7Hmm& phmm, const uint64_t sequenceLength, const float pValue);

    /// @brief finds the threshold score necessary to hit to meet the given p-value for a nucleotide phmm.
    /// @param phmm phmm to find the threshold for
    /// @param pValue p-value that will be hit at the returned score.
    /// @return score, in bits, necessary to hit the given p-value.
    float findThresholdForNucleotide(const P7Hmm& phmm, const float pValue);

    /// @brief finds the threshold score necessary to hit to meet the given p-value for an amino phmm.
    /// @param phmm phmm to find the threshold for
    /// @param sequenceLength length of the target sequence. Nucleotide Phmms use its MAXL property, but since AA models don't have this property, use the target sequence length.
    /// @param pValue p-value that will be hit at the returned score.
    /// @return score, in bits, necessary to hit the given p-value.
    float findThresholdForAmino(const P7Hmm& phmm, const uint64_t sequenceLength, const float pValue);

    /// @brief finds the inverse survival of the gumbel function
    /// @param lambda lambda value from the STATS LOCAL MSV line of the model
    /// @param mu mu value from the STATS LOCAL MSV line of the model
    /// @param pValue p-value to hit
    /// @return  inverse survival value
    float gumbelInverseSurvival(const float lambda, const float mu, const float pValue);

    /// @brief finds the number of match scores in the phmm, defined as model length * alphabet cardinality
    /// @param phmm model to query
    /// @return length of the model multiplied by the size of the model's alphabet
    uint64_t getNumMatchScores(const P7Hmm& phmm);

};

#endif