#ifndef VCFX_DUPLICATE_REMOVER_H
#define VCFX_DUPLICATE_REMOVER_H

#include <iostream>
#include <string>
#include <unordered_set>

// Function to display help message
void printHelp();

// Structure to represent a unique variant key
struct VariantKey {
    std::string chrom;
    int pos;
    std::string ref;
    std::string alt;

    bool operator==(const VariantKey& other) const {
        return chrom == other.chrom && pos == other.pos && ref == other.ref && alt == other.alt;
    }
};

// Hash function for VariantKey
struct VariantKeyHash {
    std::size_t operator()(const VariantKey& k) const {
        return std::hash<std::string>()(k.chrom) ^ std::hash<int>()(k.pos) ^ std::hash<std::string>()(k.ref) ^ std::hash<std::string>()(k.alt);
    }
};

// Function to remove duplicate variants
bool removeDuplicates(std::istream& in, std::ostream& out);

#endif // VCFX_DUPLICATE_REMOVER_H
