#include "PhmmProcessor.hpp"
#include <cmath>
#include <limits.h>
#include <stdexcept>
#include <iostream>


namespace NailForge::PhmmProcessor {

#define NAIL_NAT_LOG_2 0.69314718055994529f
#define NAIL_NAT_LOG_ONE_HALF -0.69314718055994529f
#define NAIL_NAT_LOG_2_E 1.4426950408889634073599246f

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

    std::vector<float> toFloatMatchScores(const P7Hmm& phmm) {
        const bool isAminoAlphabet = phmm.header.alphabet == P7HmmReaderAlphabetAmino;
        const uint64_t numMatchScore = PhmmProcessor::getNumMatchScores(phmm);

        std::vector<float> projectedScores;
        projectedScores.resize(numMatchScore);
        if (isAminoAlphabet) {
            for (uint32_t scoreIndex = 0; scoreIndex < numMatchScore; scoreIndex++) {
                const float negLogScore = phmm.model.matchEmissionScores[scoreIndex];
                const float background = log2InverseAminoBackgroundDistribution[scoreIndex % 20];
                const float projectedScore = background - (negLogScore * NAIL_NAT_LOG_2_E);
                projectedScores[scoreIndex] = projectedScore;
            }
        }
        else {
            for (uint32_t scoreIndex = 0; scoreIndex < numMatchScore; scoreIndex++) {
                const float negLogScore = phmm.model.matchEmissionScores[scoreIndex];
                const float background = 2.0f; //DNA or RNA, here log2(1/(1/4)) =log2(4) = 2
                const float projectedScore = background - (negLogScore * NAIL_NAT_LOG_2_E);
                projectedScores[scoreIndex] = projectedScore;
            }
        }

        return projectedScores;
    }


    float findThreshold(const P7Hmm& phmm, const float pValue, const size_t sequenceLength) {
        const float mu = phmm.stats.msvGumbelMu;
        const float lambda = phmm.stats.msvGumbelLambda;
        const float modelLength = phmm.header.modelLength;
        const float maxAlignLength = phmm.header.maxLength != 0 ?
            phmm.header.maxLength : sequenceLength;

        //given the probability of survival, give me the score. what score gives the p value?
        //what score do we need to hit the required p value? this value is for the full model. 
        const float scoreRequiredForFullModelPvalue = PhmmProcessor::gumbelInverseSurvival(lambda, mu, pValue); //inverse survival (gumbel distribution)
        //we're only working with the single-hit model, so we need to add some penalties to make the probabilities match.
        //finds the probablility of staying in the n state (or c state), looping, summed for all positions in the target
        const float lengthWithAdditionalStates = maxAlignLength + 3.0f;
        const float nStateLoopPenaltyTotal = logf(maxAlignLength / lengthWithAdditionalStates) * maxAlignLength;
        const float nStateEscapePenalty = logf(3.0f / lengthWithAdditionalStates);  //penalty for escaping the n or c states (increases by 1 every time seq len doubles)
        const float bStateToAnyMStatePenalty = logf(2.0f / (modelLength * (modelLength + 1.0f)));  //penalty for transition from b to kth 'm' option (msubk)
        const float transitionEToC = logf(1.0f / 2.0f);
        const float coreModelAdjustment = (nStateEscapePenalty + nStateLoopPenaltyTotal +
            nStateEscapePenalty + bStateToAnyMStatePenalty + transitionEToC);

        const double backgroundLoopProbability = (double)maxAlignLength / ((double)maxAlignLength + 1.0L); //background loop probability
        const float backgroundLoopPenaltyTotal = maxAlignLength * log(backgroundLoopProbability);  //total length penalty for the max sequence length background
        //this is equivalent to log(1.0 - (floatSeqLen / (floatSeqLen+1))), but has better numerical stability
        const float backgroundMovePenalty = -logf(maxAlignLength + 1.0L);
        const float backgroundScore = backgroundLoopPenaltyTotal + backgroundMovePenalty;

        const float thresholdScoreInNats =
            (scoreRequiredForFullModelPvalue * NAIL_NAT_LOG_2) + backgroundScore - coreModelAdjustment;
        // that's in nats; let's shift to bits
        float thresholdScoreInBits = thresholdScoreInNats / NAIL_NAT_LOG_2;

        return thresholdScoreInBits;
    }

    float gumbelInverseSurvival(const float lambda, const float mu, const float pValue) {
        /* The real calculation is mu - ( log(-1. * log(1-p)) / lambda).
        *  But there's a problem with small p:
        *     for p<1e-15, 1-p will be viewed as 1, so
        *     log ( -log(1-p) ) == log (0) -> inf
        *  Instead, use two approximations;
        *    (1) log( 1-p) ~= -p   for small p (e.g. p<0.001)
        *      so log(-1. * log(1-p)) ~= log(p)
        *    (2) log (p) ~= (p^p - 1) / p
        */

        const float NAIL_GUMBEL_EPSILON = 5e-9;
        float log_part = (pValue < NAIL_GUMBEL_EPSILON) ?
            (std::pow(pValue, pValue) - 1) / pValue :
            log(-1. * log(1 - pValue));
        return mu - (log_part / lambda);
    }


    uint64_t getNumMatchScores(const P7Hmm& phmm) {
        uint64_t alphabetCardinality;
        switch (phmm.header.alphabet) {
        case P7HmmReaderAlphabetDna:alphabetCardinality = 4;        break;
        case P7HmmReaderAlphabetRna:alphabetCardinality = 4;        break;
        case P7HmmReaderAlphabetAmino:alphabetCardinality = 20;     break;
        default:
            throw std::runtime_error("Unsupported phmm alphabet. Alphabet must be Dna, Rna, or Amino");
        }

        return alphabetCardinality * phmm.header.modelLength;
    }

}