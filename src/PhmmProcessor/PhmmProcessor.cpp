#include "PhmmProcessor.hpp"
#include <cmath>
#include <limits.h>
#include <stdexcept>
#include <iostream>

#define NAIL_NAT_LOG_2 0.69314718055994529
#define NAIL_NAT_LOG_ONE_FOURTH -1.3862943611198906
#define NAIL_NAT_LOG_2_E 1.4426950408889634073599246



float log2InverseAminoBackgroundDistribution[20] = {
    //these values represnt log2(1/b), where b is 
    //the background distribution for each letter
    3.6657612593,   //A
    6.0435864363,   //C
    4.2237187714,   //D
    3.9033646331,   //E
    4.6544918925,   //F
    3.8466958362,   //G
    5.4472617347,   //H
    4.0829162902,   //I
    4.0723686777,   //K
    3.3752301682,   //L
    5.3946050415,   //M
    4.5928809273,   //N
    4.3721197761,   //P
    4.6596715434,   //Q
    4.2082862647,   //R
    3.8712019420,   //S
    4.2090625204,   //T
    3.8923560482,   //V
    6.4531149215,   //W
    5.0391538251,   //Y
};

std::vector<int8_t> PhmmProcessor::to8BitMatchScores(const P7Hmm& phmm) {
    const float LOG2_E = 1.44269504089;
    const float constantAlpha = 2;
    const float constantBeta = LOG2_E;

    uint64_t numMatchScore = PhmmProcessor::getNumMatchScores(phmm);

    std::vector<int8_t> projectedScores;
    projectedScores.reserve(numMatchScore);
    for (uint32_t scoreIndex = 0; scoreIndex < numMatchScore; scoreIndex++) {
        float negLogScore = phmm.model.matchEmissionScores[scoreIndex];
        float projectedScore = constantAlpha - (negLogScore * constantBeta);

        if (__builtin_expect(projectedScore < SCHAR_MIN, 0)) {
            projectedScore = SCHAR_MIN;
        }
        else if (__builtin_expect(projectedScore > SCHAR_MAX, 0)) {
            projectedScore = SCHAR_MAX;
        }
        projectedScores.push_back((int8_t)projectedScore);
    }

    return projectedScores;
}


std::vector<int8_t> PhmmProcessor::toProjected8BitMatchScores(const P7Hmm& phmm, const float pValue, const float thresholdHitValue) {
    const float scoreMultiplier = PhmmProcessor::findTauScalingFactor(phmm, pValue, thresholdHitValue);
    const float LOG2_E = 1.44269504089;
    const float constantAlpha = 2 * scoreMultiplier;
    const float constantBeta = LOG2_E * scoreMultiplier;


    uint64_t numMatchScore = PhmmProcessor::getNumMatchScores(phmm);

    std::vector<int8_t> projectedScores;
    projectedScores.reserve(numMatchScore);
    for (uint32_t scoreIndex = 0; scoreIndex < numMatchScore; scoreIndex++) {
        float negLogScore = phmm.model.matchEmissionScores[scoreIndex];
        float projectedScore = constantAlpha - (negLogScore * constantBeta);

        if (__builtin_expect(projectedScore < SCHAR_MIN, 0)) {
            projectedScore = SCHAR_MIN;
        }
        else if (__builtin_expect(projectedScore > SCHAR_MAX, 0)) {
            projectedScore = SCHAR_MAX;
        }
        projectedScores.push_back((int8_t)projectedScore);

    }

    return projectedScores;
}



std::vector<float> PhmmProcessor::toFloatMatchScores(const P7Hmm& phmm) {
    const bool isAminoAlphabet = phmm.header.alphabet == P7HmmReaderAlphabetAmino;
    const float constantAlpha = isAminoAlphabet ? std::log2(20) : 2;
    const uint64_t numMatchScore = PhmmProcessor::getNumMatchScores(phmm);

    std::vector<float> projectedScores;
    projectedScores.reserve(numMatchScore);

    if (isAminoAlphabet) {
        for (uint32_t scoreIndex = 0; scoreIndex < numMatchScore; scoreIndex++) {
            const float negLogScore = phmm.model.matchEmissionScores[scoreIndex];
            const float background = log2InverseAminoBackgroundDistribution[scoreIndex % 20];
            const float projectedScore = background - (negLogScore * NAIL_NAT_LOG_2_E);
            projectedScores.push_back(projectedScore);
        }
    }
    else {
        for (uint32_t scoreIndex = 0; scoreIndex < numMatchScore; scoreIndex++) {
            const float negLogScore = phmm.model.matchEmissionScores[scoreIndex];
            const float background = 2.0f; //DNA or RNA, here log2(1/(1/4)) =log2(4) = 2
            const float projectedScore = background - (negLogScore * NAIL_NAT_LOG_2_E);
            projectedScores.push_back(projectedScore);
        }
    }

    return projectedScores;
}


float PhmmProcessor::findThreshold(const P7Hmm& phmm, const float pValue) {

    std::cout << "NOTICE: This function is currently Deprecated. It cannot function for Amino Acid searching. "
        << "ignoring threshold computation, using hard-coded value of 15" << std::endl;
    return 15.0f;

    const float mu = phmm.stats.msvGumbelMu;
    const float lambda = phmm.stats.msvGumbelLambda;
    const float maxLength = phmm.header.maxLength;
    const float modelLength = phmm.header.modelLength;

    //given the probability of survival, give me the score. what score gives the p value?
    //what score do we need to hit the required p value? this value is for the full model. 
    float scoreRequiredForFullModelPvalue = PhmmProcessor::gumbelInverseSurvival(lambda, mu, pValue); //inverse survival (gumbel distribution)
    //we're only working with the single-hit model, so we need to add some penalties to make the probabilities match.
    float nStateLoopPenalty = logf(maxLength / (maxLength + 3));  //probablility of staying in the n state (or c state), looping.
    float nStateLoopPenaltyTotal = nStateLoopPenalty * maxLength;  //total length penalty for the max sequence length
    float nStateEscapePenalty = logf(3.0f / (maxLength + 3));  //penalty for escaping the n or c states (increases by 1 every time seq len doubles)
    float bStateToAnyMStatePenalty = logf(2.0f / (modelLength * (modelLength + 1)));  //penalty for transition from b to kth 'm' option (msubk)
    float transitionEToC = logf(1.0f / 2);
    float coreModelAdjustment = (nStateEscapePenalty + nStateLoopPenaltyTotal +
        nStateEscapePenalty + bStateToAnyMStatePenalty + transitionEToC);

    float backgroundLoopProbability = maxLength / (maxLength + 1); //background loop probability
    float backgroundLoopPenaltyTotal = maxLength * log(backgroundLoopProbability);  //total length penalty for the max sequence length background
    float backgroundMovePenalty = log(1.0 - backgroundLoopProbability);
    float backgroundScore = backgroundLoopPenaltyTotal + backgroundMovePenalty;

    float thresholdScoreInNats =
        (scoreRequiredForFullModelPvalue * NAIL_NAT_LOG_2) + backgroundScore - coreModelAdjustment;
    // that's in nats; let's shift to bits
    float thresholdScoreInBits = thresholdScoreInNats / NAIL_NAT_LOG_2;

    return thresholdScoreInBits;
}

float PhmmProcessor::findTauScalingFactor(const P7Hmm& phmm, const float pValue, const float projectedThreshold) {
    float thresholdScoreInBits = PhmmProcessor::findThreshold(phmm, pValue);
    float scaleFactor = projectedThreshold / thresholdScoreInBits; //a hit should be triggered if a cell hits 256 (causes an 8-bit int overflow)

    return scaleFactor;
}

float PhmmProcessor::gumbelInverseSurvival(const float lambda, const float mu, const float pValue) {
    /* The real calculation is mu - ( log(-1. * log(1-p)) / lambda).
    *  But there's a problem with small p:
    *     for p<1e-15, 1-p will be viewed as 1, so
    *     log ( -log(1-p) ) == log (0) -> inf
    *  Instead, use two approximations;
    *    (1) log( 1-p) ~= -p   for small p (e.g. p<0.001)
    *      so log(-1. * log(1-p)) ~= log(p)
    *    (2) log (p) ~= (p^p - 1) / p
    *
    *    See notes Mar 1, 2010.
    */

    const float NAIL_GUMBEL_EPSILON = 5e-9;
    float log_part = (pValue < NAIL_GUMBEL_EPSILON) ?
        (std::pow(pValue, pValue) - 1) / pValue :
        log(-1. * log(1 - pValue));
    return mu - (log_part / lambda);
}


uint64_t PhmmProcessor::getNumMatchScores(const P7Hmm& phmm) {
    uint64_t alphabetCardinality;
    switch (phmm.header.alphabet) {
    case P7HmmReaderAlphabetDna:alphabetCardinality = 4;        break;
    case P7HmmReaderAlphabetRna:alphabetCardinality = 4;        break;
    case P7HmmReaderAlphabetAmino:alphabetCardinality = 20;     break;
    case P7HmmReaderAlphabetCoins: alphabetCardinality = 2;     break;
    case P7HmmReaderAlphabetDice:  alphabetCardinality = 6;     break;
    default:
        throw std::runtime_error("Unsupported phmm alphabet");
    }

    return alphabetCardinality * phmm.header.modelLength;
}

