#include "VCFX_variant_classifier.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <unordered_set>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_variant_classifier\n"
              << "Usage: VCFX_variant_classifier [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Classifies variants in a VCF file as SNPs, indels, MNVs, or structural variants based on the REF and ALT alleles.\n\n"
              << "Examples:\n"
              << "  ./VCFX_variant_classifier < input.vcf > classified_variants.tsv\n";
}

// Function to convert VariantType enum to string
std::string variantTypeToString(VariantType type) {
    switch (type) {
        case VariantType::SNP:
            return "SNP";
        case VariantType::INDEL:
            return "Indel";
        case VariantType::MNV:
            return "MNV";
        case VariantType::STRUCTURAL_VARIANT:
            return "Structural_Variant";
        default:
            return "Unknown";
    }
}

// Function to classify a single allele pair
VariantType classifyAllele(const std::string& ref, const std::string& alt) {
    // Check for symbolic alleles indicating structural variants
    if (!alt.empty() && alt[0] == '<' && alt.back() == '>') {
        return VariantType::STRUCTURAL_VARIANT;
    }

    // Check for simple SNP
    if (ref.length() == 1 && alt.length() == 1 &&
        std::isalpha(ref[0]) && std::isalpha(alt[0])) {
        return VariantType::SNP;
    }

    // Check for indel
    if (ref.length() != alt.length()) {
        // Length difference significant enough to consider as structural variant
        if (ref.length() > 50 || alt.length() > 50) { // Arbitrary threshold
            return VariantType::STRUCTURAL_VARIANT;
        }
        return VariantType::INDEL;
    }

    // Check for multi-nucleotide variant (MNV)
    if (ref.length() > 1 && alt.length() > 1) {
        return VariantType::MNV;
    }

    return VariantType::UNKNOWN;
}

// Function to classify a VCF variant
VariantType classifyVariant(const std::string& ref, const std::vector<std::string>& alt) {
    std::unordered_set<VariantType> types;

    for (const auto& allele : alt) {
        VariantType type = classifyAllele(ref, allele);
        types.insert(type);
    }

    // Prioritize variant types
    if (types.find(VariantType::STRUCTURAL_VARIANT) != types.end()) {
        return VariantType::STRUCTURAL_VARIANT;
    }
    if (types.find(VariantType::MNV) != types.end()) {
        return VariantType::MNV;
    }
    if (types.find(VariantType::INDEL) != types.end()) {
        return VariantType::INDEL;
    }
    if (types.find(VariantType::SNP) != types.end()) {
        return VariantType::SNP;
    }

    return VariantType::UNKNOWN;
}

// Function to split a string by a delimiter and trim whitespace
std::vector<std::string> splitAndTrim(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        // Trim whitespace
        size_t start = token.find_first_not_of(" \t");
        size_t end = token.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            tokens.emplace_back(token.substr(start, end - start + 1));
        } else if (start != std::string::npos) {
            tokens.emplace_back(token.substr(start));
        } else {
            tokens.emplace_back(""); // All whitespace
        }
    }
    return tokens;
}

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant) {
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    while (std::getline(ss, field, '\t')) {
        fields.emplace_back(field);
    }

    if (fields.size() < 8) {
        std::cerr << "Warning: Skipping invalid VCF line with fewer than 8 fields.\n";
        return false;
    }

    variant.chrom = fields[0];
    try {
        variant.pos = std::stoi(fields[1]);
    } catch (...) {
        std::cerr << "Warning: Invalid POS value. Skipping line.\n";
        return false;
    }
    variant.id = fields[2];
    variant.ref = fields[3];
    
    // ALT can have multiple alleles separated by ','
    variant.alt = splitAndTrim(fields[4], ',');
    variant.qual = fields[5];
    variant.filter = fields[6];
    variant.info = fields[7];

    // Classify variant
    variant.type = classifyVariant(variant.ref, variant.alt);

    return true;
}

// Function to classify variants in VCF
bool classifyVariants(std::istream& in, std::ostream& out) {
    std::string line;
    bool header_found = false;

    // Print header for classification
    out << "CHROM\tPOS\tID\tREF\tALT\tVARIANT_TYPE\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            continue; // Skip header lines
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        VCFVariant variant;
        if (!parseVCFLine(line, variant)) {
            continue; // Skip invalid lines
        }

        // Prepare ALT field as comma-separated string
        std::string alt_str = "";
        for (size_t i = 0; i < variant.alt.size(); ++i) {
            alt_str += variant.alt[i];
            if (i != variant.alt.size() - 1) {
                alt_str += ",";
            }
        }

        out << variant.chrom << "\t"
            << variant.pos << "\t"
            << variant.id << "\t"
            << variant.ref << "\t"
            << alt_str << "\t"
            << variantTypeToString(variant.type) << "\n";
    }

    return true;
}

int main(int argc, char* argv[]) {
    // Argument parsing for help
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Classify variants
    bool success = classifyVariants(std::cin, std::cout);
    return success ? 0 : 1;
}
