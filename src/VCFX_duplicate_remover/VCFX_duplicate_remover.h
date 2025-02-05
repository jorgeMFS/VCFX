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
    std::string alt; // normalized: sorted, comma-separated alleles

    bool operator==(const VariantKey& other) const {
        return chrom == other.chrom &&
               pos == other.pos &&
               ref == other.ref &&
               alt == other.alt;
    }
};

// Custom hash function for VariantKey using a hash-combine approach
struct VariantKeyHash {
    std::size_t operator()(const VariantKey& k) const {
        std::size_t h1 = std::hash<std::string>()(k.chrom);
        std::size_t h2 = std::hash<int>()(k.pos);
        std::size_t h3 = std::hash<std::string>()(k.ref);
        std::size_t h4 = std::hash<std::string>()(k.alt);
        std::size_t seed = h1;
        seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h4 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

// Function to remove duplicate variants from a VCF file
bool removeDuplicates(std::istream& in, std::ostream& out);

#endif // VCFX_DUPLICATE_REMOVER_H
