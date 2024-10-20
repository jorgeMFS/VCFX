#ifndef VCFX_VARIANT_CLASSIFIER_H
#define VCFX_VARIANT_CLASSIFIER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Enumeration for variant types
enum class VariantType {
    SNP,
    INDEL,
    MNV,
    STRUCTURAL_VARIANT,
    UNKNOWN
};

// Function to convert VariantType enum to string
std::string variantTypeToString(VariantType type);

// Structure to represent a VCF variant
struct VCFVariant {
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::vector<std::string> alt;
    std::string qual;
    std::string filter;
    std::string info;
    VariantType type;
};

// Function to classify a single allele pair
VariantType classifyAllele(const std::string& ref, const std::string& alt);

// Function to classify a VCF variant
VariantType classifyVariant(const std::string& ref, const std::vector<std::string>& alt);

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant);

// Function to classify variants in VCF
bool classifyVariants(std::istream& in, std::ostream& out);

#endif // VCFX_VARIANT_CLASSIFIER_H
