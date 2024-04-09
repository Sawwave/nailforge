#include "MaxExtensionTable.hpp"
#include <iostream>

namespace NailForge::StringTree::Table {

    MaxExtensionTable::MaxExtensionTable(const uint32_t modelLength, const uint8_t depth, const uint8_t alphabetSize)noexcept
        :modelLength(modelLength), depth(depth), alphabetSize(alphabetSize) {
        maxScoreTable.resize(modelLength * depth);
        //fill in the first row with zeros. This conceptually means that extending any position by 0 characters won't change the score.
        std::fill_n(maxScoreTable.begin(), modelLength, 0.0f);
    }


    float& MaxExtensionTable::scoreAt(const uint32_t modelPosition, const uint8_t depth) noexcept {
        return maxScoreTable[(depth * this->modelLength) + modelPosition];
    }

    float MaxExtensionTable::getMaxExtensionScore(const uint32_t modelPosition, const uint8_t extensionLength) noexcept {
        return scoreAt(modelPosition, extensionLength);
    }

    //this is in the header to make the compiler happier with the template argument.
    /**
     * Fills in the maximum extension table. At row i, col j of the table, the value contained is the maximum value attainable
     * by extending model position j by up to i+1 letters. This can be built one of two ways: for standard search, or reverse compliment.
     * Imagine a hmm model with the following scores
     * [ [5,-1,-1,-1], [-1,10,-1,-1], [-1,-1,15,-1], [-1,-1,-1,20], [25,-1,-1,-1], [-1,30,-1,-1], [-1,-1,35,-1], [-1,-1,-1-40]]
     *
     * The table for a depth of 4 would appear as follows for standard search:
     *  [0, 0,  0,  0,  0,  0,  0,  0]
     *  [0, 5, 10, 15, 20, 25, 30, 35]
     *  [0, 5, 15, 25, 35, 45, 55, 65]
     *  [0, 5, 15, 30, 45, 60, 75, 90]
     *
     * Note that the table for standard search extends backwards. this is because the FM-index algorithm works from the end of the string
     * The reverse compliment table would look like this:
     * [ 0,  0,  0,  0,  0,  0,  0, 0]
     * [10, 15, 20, 25, 30, 35, 40, 0]
     * [25, 35, 45, 55, 65, 75, 40, 0]
     * [30, 45, 60, 75, 90, 75, 40, 0]
     *
    */
    void MaxExtensionTable::populateTable(const std::vector<float>& scores, bool isReverseCompliment) noexcept {

        //the first row of the table is all zeros, because extending anything by 0 characters won't change the score
        //the second row of the table is what you get by extending the position by one character.
        //since FM-index works from the end of the string, this means that standard search will extend backwards,
        //and reverse compliment will extend forward.
        if (isReverseCompliment) {
            //find the maximum score possible by extending every position 1 character
            for (uint32_t modelPosition = 0; modelPosition < (modelLength - 1); modelPosition++) {
                const float* scoreVector = &scores[(modelPosition + 1) * alphabetSize];
                float& thisVectorMaxScore = scoreAt(modelPosition, 1);
                thisVectorMaxScore = scoreVector[0];
                for (uint32_t scoreIdx = 1; scoreIdx < alphabetSize; scoreIdx++) {
                    if (scoreVector[scoreIdx] > thisVectorMaxScore) {
                        thisVectorMaxScore = scoreVector[scoreIdx];
                    }
                }
            }
            //the final model position can't be extended forward, so we give it a score of 0
            scoreAt(modelLength - 1, 1) = 0;

            for (uint8_t tableDepth = 2; tableDepth < this->depth; tableDepth++) {
                for (uint32_t modelPosition = 0; modelPosition < (modelLength - tableDepth); modelPosition++) {
                    scoreAt(modelPosition, tableDepth) =                //this score equals
                        scoreAt(modelPosition, tableDepth - 1) +        //best score from extending it one less than we need
                        scoreAt(modelPosition + (tableDepth - 1), 1);     //plus the final letter extension
                    //this final extension is found by walking diagonally back down the table until the 1st extension row
                }
                //now, complete the row by finishing the cells that would have gone past the bounds of the model
                //As the cells represent the max score from extending UP TO that many positions, we just copy
                //the score from the cell extending one less than what we already have.
                for (uint32_t modelPosition = (modelLength - tableDepth); modelPosition < modelLength;modelPosition++) {
                    scoreAt(modelPosition, tableDepth) = scoreAt(modelPosition, tableDepth - 1);
                }
            }
        }
        else {
            //we're extending backwards for standard search.
            //the first model position can't be extended forward, so we give it a score of 0
            scoreAt(0, 1) = 0;

            //find the maximum score possible by extending every position 1 character
            for (uint32_t modelPosition = 1; modelPosition < modelLength; modelPosition++) {
                const float* scoreVector = &scores[(modelPosition - 1) * alphabetSize];
                float& thisVectorMaxScore = scoreAt(modelPosition, 1);
                thisVectorMaxScore = scoreVector[0];
                for (uint32_t scoreIdx = 1; scoreIdx < alphabetSize; scoreIdx++) {
                    if (scoreVector[scoreIdx] > thisVectorMaxScore) {
                        thisVectorMaxScore = scoreVector[scoreIdx];
                    }
                }
            }

            for (uint8_t tableDepth = 2; tableDepth < this->depth; tableDepth++) {
                for (uint32_t modelPosition = 0; modelPosition < modelLength; modelPosition++) {
                    if (tableDepth > modelPosition) {
                        scoreAt(modelPosition, tableDepth) = scoreAt(modelPosition, tableDepth - 1);
                    }
                    else {
                        scoreAt(modelPosition, tableDepth) =                //this score equals
                            scoreAt(modelPosition, tableDepth - 1) +        //best score from extending it one less than we need
                            scoreAt(modelPosition - (tableDepth - 1), 1);     //plus the final letter extension
                        //this final extension is found by walking diagonally back down the table until the 1st extension row
                    }
                }
            }
        }
    }
}
