#ifndef NAIL_FORGE_MAX_EXTENSION_TABLE_HPP
#define NAIL_FORGE_MAX_EXTENSION_TABLE_HPP
#include <cstdint>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>

namespace NailForge::StringTree::Table
{
    class MaxExtensionTable {
    public:
        MaxExtensionTable(const uint32_t modelLength, const uint8_t depth, const uint8_t alphabetSize)noexcept;
        void populateTable(const std::vector<float>& scores, bool isReverseComplement) noexcept;
        float getMaxExtensionScore(const uint32_t modelPosition, const uint8_t extensionLength) noexcept;

    protected:

        uint32_t modelLength;
        uint8_t depth;
        uint8_t alphabetSize;
        std::vector<float> maxScoreTable;

        float& scoreAt(const uint32_t modelPosition, const uint8_t depth) noexcept;
    };
}

#endif