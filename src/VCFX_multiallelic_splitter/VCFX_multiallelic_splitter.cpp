#include "VCFX_multiallelic_splitter.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_multiallelic_splitter\n"
              << "Usage: VCFX_multiallelic_splitter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Splits multi-allelic variants in a VCF file into multiple bi-allelic records.\n\n"
              << "Examples:\n"
              << "  ./VCFX_multiallelic_splitter < input.vcf > split_biallelic.vcf\n";
}

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant) {
    if (line.empty() || line[0] == '#') return false;

    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    // Split the line by tab
    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    // VCF must have at least 8 fields before FORMAT
    if (fields.size() < 8) return false;

    variant.chrom = fields[0];
    try {
        variant.pos = std::stoi(fields[1]);
    } catch (...) {
        return false;
    }
    variant.id = fields[2];
    variant.ref = fields[3];
    // Split ALT alleles by comma
    std::stringstream alt_ss(fields[4]);
    std::string alt_allele;
    while (std::getline(alt_ss, alt_allele, ',')) {
        variant.alt.push_back(alt_allele);
    }
    variant.qual = fields[5];
    variant.filter = fields[6];
    variant.info = fields[7];

    // Collect FORMAT and SAMPLE fields if present
    if (fields.size() > 8) {
        variant.other_fields.assign(fields.begin() + 8, fields.end());
    }

    return true;
}

// Function to reconstruct INFO field
std::string reconstructInfo(const std::string& chrom, int pos, const VCFVariant& variant) {
    return variant.info; // For simplicity, retain the original INFO field
}

// Function to split multi-allelic variants into bi-allelic
bool splitMultiAllelicVariants(std::istream& in, std::ostream& out) {
    std::string line;
    bool header_found = false;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Pass through header lines unchanged
            out << line << "\n";
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        VCFVariant variant;
        if (!parseVCFLine(line, variant)) {
            std::cerr << "Warning: Skipping malformed VCF line: " << line << "\n";
            continue;
        }

        // If only one ALT allele, write the record as is
        if (variant.alt.size() <= 1) {
            out << line << "\n";
            continue;
        }

        // More than one ALT allele, split into multiple records
        for (size_t i = 0; i < variant.alt.size(); ++i) {
            VCFVariant split_variant = variant;
            split_variant.alt = { variant.alt[i] };
            split_variant.info = reconstructInfo(variant.chrom, variant.pos, variant);

            // Reconstruct the VCF line
            std::stringstream split_ss;
            split_ss << split_variant.chrom << "\t"
                     << split_variant.pos << "\t"
                     << split_variant.id << "\t"
                     << split_variant.ref << "\t";

            // ALT field
            split_ss << split_variant.alt[0] << "\t"
                     << split_variant.qual << "\t"
                     << split_variant.filter << "\t"
                     << split_variant.info;

            // Append FORMAT and sample fields if present
            if (!split_variant.other_fields.empty()) {
                for (const auto& sample_field : split_variant.other_fields) {
                    split_ss << "\t" << sample_field;
                }
            }

            split_ss << "\n";
            out << split_ss.str();
        }
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

    // Split multi-allelic variants
    bool success = splitMultiAllelicVariants(std::cin, std::cout);
    return success ? 0 : 1;
}
