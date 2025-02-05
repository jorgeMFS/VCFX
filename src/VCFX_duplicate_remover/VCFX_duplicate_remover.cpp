#include "VCFX_duplicate_remover.h"
#include <sstream>
#include <vector>
#include <algorithm>

// Utility function: Splits a string by a given delimiter.
static std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to display help message
void printHelp() {
    std::cout << "VCFX_duplicate_remover\n"
              << "Usage: VCFX_duplicate_remover [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Removes duplicate variants from a VCF file based on the combination of\n"
              << "  chromosome, position, REF, and ALT alleles. For multi-allelic records, the\n"
              << "  ALT field is normalized by sorting the comma-separated alleles so that the\n"
              << "  ordering does not affect duplicate detection.\n\n"
              << "Example:\n"
              << "  ./VCFX_duplicate_remover < input.vcf > unique_variants.vcf\n";
}

// Generates a unique key for a variant based on chrom, pos, ref, and alt.
// For multi-allelic ALT fields, the ALT alleles are split, sorted, and rejoined.
std::string generateNormalizedVariantKey(const std::string& chrom,
                                          const std::string& pos,
                                          const std::string& ref,
                                          const std::string& alt) {
    // Split ALT field on commas.
    std::vector<std::string> alts = splitString(alt, ',');
    // Sort alleles lexicographically.
    std::sort(alts.begin(), alts.end());
    // Rejoin sorted alleles into a single string.
    std::ostringstream oss;
    for (size_t i = 0; i < alts.size(); ++i) {
        if (i > 0) {
            oss << ",";
        }
        oss << alts[i];
    }
    std::string normalizedAlt = oss.str();
    return chrom + ":" + pos + ":" + ref + ":" + normalizedAlt;
}

// Helper function to generate a VariantKey from parsed values.
static VariantKey generateVariantKey(const std::string& chrom,
                                     const std::string& pos,
                                     const std::string& ref,
                                     const std::string& alt) {
    VariantKey key;
    key.chrom = chrom;
    try {
        key.pos = std::stoi(pos);
    } catch (...) {
        key.pos = 0;
    }
    key.ref = ref;
    key.alt = "";  // Will be set to normalized ALT.
    // Normalize ALT: sort multi-allelic entries.
    key.alt = generateNormalizedVariantKey(chrom, pos, ref, alt).substr(chrom.size() + pos.size() + ref.size() + 3); // skip prefix "chrom:pos:ref:"
    // Alternatively, simply:
    key.alt = generateNormalizedVariantKey(chrom, pos, ref, alt);
    // However, since generateNormalizedVariantKey already concatenates chrom:pos:ref:normalizedAlt,
    // we extract the normalizedAlt portion if needed. For simplicity, we can just store the full key.
    // For our VariantKey, we want: chrom, pos, ref, normalizedAlt.
    // We'll do that by re-parsing:
    std::vector<std::string> parts = splitString(generateNormalizedVariantKey(chrom, pos, ref, alt), ':');
    if (parts.size() >= 4) {
        key.chrom = parts[0];
        try {
            key.pos = std::stoi(parts[1]);
        } catch (...) {
            key.pos = 0;
        }
        key.ref = parts[2];
        key.alt = parts[3];
    }
    return key;
}

// Function to remove duplicate variants from a VCF file.
bool removeDuplicates(std::istream& in, std::ostream& out) {
    std::string line;
    // Print header lines as-is.
    // Use an unordered_set with our custom hash function to track seen variants.
    std::unordered_set<VariantKey, VariantKeyHash> seen_variants;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // Print header lines unchanged.
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info;
        // Expect at least 8 columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        // Generate a variant key that normalizes ALT ordering.
        VariantKey key = generateVariantKey(chrom, pos, ref, alt);
        if (seen_variants.find(key) == seen_variants.end()) {
            // Variant is unique; record and output the line.
            seen_variants.insert(key);
            out << line << "\n";
        } else {
            // Duplicate variant; skip.
        }
    }
    return true;
}

// ----------------------------------------------------------------------
// main: Parse command-line arguments and call removeDuplicates.
// ----------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Simple argument parsing: if --help or -h is provided, print help.
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Remove duplicates from standard input and write to standard output.
    bool success = removeDuplicates(std::cin, std::cout);
    return success ? 0 : 1;
}
