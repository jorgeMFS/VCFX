#include "vcfx_core.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cctype>

// ---------------------------------------------------------
// Print help message
// ---------------------------------------------------------
static void printHelp() {
    std::cout 
        << "VCFX_allele_freq_calc\n"
        << "Usage: VCFX_allele_freq_calc [OPTIONS]\n\n"
        << "Options:\n"
        << "  --help, -h   Display this help message and exit.\n\n"
        << "Description:\n"
        << "  Reads a VCF from stdin and outputs a TSV file:\n"
        << "    CHROM  POS  ID  REF  ALT  Allele_Frequency\n\n"
        << "  Allele frequency is computed as (#ALT alleles / total #alleles),\n"
        << "  counting any non-zero numeric allele (1,2,3,...) as ALT.\n\n"
        << "Example:\n"
        << "  ./VCFX_allele_freq_calc < input.vcf > allele_frequencies.tsv\n";
}

// ---------------------------------------------------------
// Utility: split a string by a delimiter
// ---------------------------------------------------------
static std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// ---------------------------------------------------------
// parseGenotype: counts how many are alt vs total among numeric alleles
// ---------------------------------------------------------
static void parseGenotype(const std::string& genotype, int &alt_count, int &total_count) {
    // e.g. genotype might be "0/1" or "2|3", etc.
    // First, unify the separators by replacing '|' with '/'
    std::string gt = genotype;
    for (char &c : gt) {
        if (c == '|') c = '/';
    }

    // Split on '/'
    auto alleles = split(gt, '/');

    for (const auto & allele : alleles) {
        if (allele.empty() || allele == ".") {
            // Skip missing
            continue;
        }
        // Check if it's numeric
        bool numeric = true;
        for (char c : allele) {
            if (!std::isdigit(c)) {
                numeric = false;
                break;
            }
        }
        if (!numeric) {
            // Not a numeric allele -> skip
            continue;
        }
        // If it's "0", it's reference; otherwise it's alt
        if (allele == "0") {
            // The reference allele still contributes to the total
            total_count += 1;
        } else {
            // Any non-zero numeric allele => alt
            alt_count += 1;
            total_count += 1;
        }
    }
}

// ---------------------------------------------------------
// Main calculation: read VCF from 'in', write results to 'out'
// ---------------------------------------------------------
static void calculateAlleleFrequency(std::istream& in, std::ostream& out) {
    std::string line;

    // We'll track how many samples are in the VCF, though we only need it
    // to confirm that columns 10..N contain sample data
    int sample_count = 0;
    bool foundChromHeader = false;

    // Print a single TSV header for our results
    out << "CHROM\tPOS\tID\tREF\tALT\tAllele_Frequency\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // If this is the #CHROM line, figure out how many samples
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                auto fields = split(line, '\t');
                if (fields.size() > 9) {
                    sample_count = static_cast<int>(fields.size() - 9);
                }
            }
            // Skip printing VCF headers in our output
            continue;
        }

        // We expect #CHROM line before actual data lines
        if (!foundChromHeader) {
            std::cerr << "Warning: Data line encountered before #CHROM header. Skipping line:\n"
                      << line << "\n";
            continue;
        }

        // Split the data line
        auto fields = split(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line (fewer than 9 fields):\n"
                      << line << "\n";
            continue;
        }

        // Standard VCF columns:
        //  0:CHROM, 1:POS, 2:ID, 3:REF, 4:ALT, 5:QUAL, 6:FILTER, 7:INFO, 8:FORMAT, 9+:samples
        const std::string &chrom = fields[0];
        const std::string &pos   = fields[1];
        const std::string &id    = fields[2];
        const std::string &ref   = fields[3];
        const std::string &alt   = fields[4];
        const std::string &fmt   = fields[8];

        // Look for the "GT" field index in the FORMAT column
        auto formatFields = split(fmt, ':');
        int gtIndex = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }
        if (gtIndex < 0) {
            // No GT field found => we cannot parse genotypes
            // skip this line
            continue;
        }

        int alt_count = 0;
        int total_count = 0;

        // For each sample column (starting at index 9)
        for (size_t i = 9; i < fields.size(); ++i) {
            // example sample string: "0/1:30,10:99:..." if GT is first
            auto sampleTokens = split(fields[i], ':');
            if (static_cast<size_t>(gtIndex) >= sampleTokens.size()) {
                // This sample has no GT data
                continue;
            }
            // parse its genotype
            parseGenotype(sampleTokens[gtIndex], alt_count, total_count);
        }

        // Compute allele frequency = alt_count / total_count
        double freq = 0.0;
        if (total_count > 0) {
            freq = static_cast<double>(alt_count) / static_cast<double>(total_count);
        }

        // Print the row with a fixed decimal precision
        out << chrom << "\t" << pos << "\t" << id << "\t"
            << ref   << "\t" << alt << "\t"
            << std::fixed << std::setprecision(4) << freq << "\n";
    }
}

// ---------------------------------------------------------
// main()
// ---------------------------------------------------------
int main(int argc, char* argv[]) {
    if (vcfx::handle_version_flag(argc, argv, "VCFX_allele_freq_calc")) return 0;
    // Parse arguments for help
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // If no arguments (and no piped input), just show help
    // (This is optional; if the user pipes data, we can proceed.)
    if (argc == 1) {
        // Not strictly required, but user-friendly
        // Check if there's something on stdin, or just show help
        if (std::cin.peek() == EOF) {
            printHelp();
            return 1;
        }
    }

    // Perform allele frequency calculation
    calculateAlleleFrequency(std::cin, std::cout);
    return 0;
}
