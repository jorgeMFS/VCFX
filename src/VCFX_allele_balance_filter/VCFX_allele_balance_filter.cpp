#include "vcfx_core.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <getopt.h>

// ------------------------------------------------------
// Class Declaration
// ------------------------------------------------------
class VCFXAlleleBalanceFilter {
public:
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on allele balance threshold
    void filterByAlleleBalance(std::istream& in, std::ostream& out, double threshold);

    // Parses the genotype to calculate allele balance
    double calculateAlleleBalance(const std::string& genotype);
};

// ------------------------------------------------------
// Implementation
// ------------------------------------------------------

// Entry point for the tool
int VCFXAlleleBalanceFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    double threshold = -1.0; // invalid until set

    static struct option long_options[] = {
        {"help",                 no_argument,       0, 'h'},
        {"filter-allele-balance", required_argument, 0, 'f'},
        {0,                      0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hf:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                try {
                    threshold = std::stod(optarg);
                } catch (const std::invalid_argument&) {
                    std::cerr << "Error: Invalid threshold value.\n";
                    displayHelp();
                    return 1;
                }
                break;
            default:
                showHelp = true;
        }
    }

    // Validate threshold and possibly show help
    if (showHelp || threshold < 0.0 || threshold > 1.0) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // Perform allele balance filtering from stdin to stdout
    filterByAlleleBalance(std::cin, std::cout, threshold);

    return 0;
}

void VCFXAlleleBalanceFilter::displayHelp() {
    std::cout 
        << "VCFX_allele_balance_filter: Filter VCF variants based on allele balance ratios.\n\n"
        << "Usage:\n"
        << "  VCFX_allele_balance_filter --filter-allele-balance <THRESHOLD> [options]\n\n"
        << "Options:\n"
        << "  -h, --help                       Display this help message and exit\n"
        << "  -f, --filter-allele-balance VAL  Specify the allele balance threshold (0.0 - 1.0)\n\n"
        << "Example:\n"
        << "  VCFX_allele_balance_filter --filter-allele-balance 0.3 < input.vcf > filtered.vcf\n\n"
        << "Note:\n"
        << "  This filter lumps all non-'0' alleles (1,2,3,...) as ALT when calculating the ratio.\n"
        << "  If any sample's allele balance is < THRESHOLD, the entire variant line is skipped.\n";
}

// The core filter function
void VCFXAlleleBalanceFilter::filterByAlleleBalance(std::istream& in, std::ostream& out, double threshold) {
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // Print header lines unchanged
        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        // Parse the fixed VCF fields up to FORMAT
        // VCF line: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [samples...]
        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        // The rest are genotype columns
        std::vector<std::string> genotypes;
        std::string genotype;
        while (ss >> genotype) {
            genotypes.push_back(genotype);
        }

        // "all or nothing" filter: if ANY genotype fails threshold => skip
        bool pass = true;
        for (const auto& gtField : genotypes) {
            double ab = calculateAlleleBalance(gtField);
            if (ab < threshold) {
                pass = false;
                break;
            }
        }

        if (pass) {
            out << line << "\n";
        }
    }
}

// Calculates allele balance as refCount / (refCount + altCount)
// *all* non-zero numeric alleles are counted as alt
double VCFXAlleleBalanceFilter::calculateAlleleBalance(const std::string& genotype) {
    // Genotype might be "0/1:...". We want the portion before the first ':'
    size_t colonPos = genotype.find(':');
    std::string gt = (colonPos == std::string::npos) ? genotype : genotype.substr(0, colonPos);

    // Replace '|' with '/' for simpler splitting
    for (auto &ch : gt) {
        if (ch == '|') ch = '/';
    }

    // Split on '/'
    std::stringstream ss(gt);
    std::string token;
    std::vector<std::string> alleles;
    while (std::getline(ss, token, '/')) {
        alleles.push_back(token);
    }

    int ref_count = 0;
    int alt_count = 0;

    // Count '0' as ref, any other number as alt
    for (auto & allele : alleles) {
        // If allele is empty or ".", skip it
        if (allele.empty() || allele == ".") {
            continue;
        }
        // If we can parse it as an integer:
        //   0 => ref, 1 or 2 or 3 => alt
        // If it's non-numeric, we skip it as missing.
        bool numeric = true;
        for (char c : allele) {
            if (!isdigit(c)) {
                numeric = false;
                break;
            }
        }
        if (!numeric) {
            // Non-numeric => skip?
            continue; 
        }
        // numeric => check if it's "0" or not
        if (allele == "0") {
            ref_count++;
        } else {
            alt_count++;
        }
    }

    int total = ref_count + alt_count;
    if (total == 0) {
        // No valid calls => ratio 0.0 or treat as "no data"
        return 0.0;
    }
    return static_cast<double>(ref_count) / total;
}

// ------------------------------------------------------
// main() linking to class
// ------------------------------------------------------
int main(int argc, char* argv[]) {
    if (vcfx::handle_version_flag(argc, argv, "VCFX_allele_balance_filter")) return 0;
    VCFXAlleleBalanceFilter alleleBalanceFilter;
    return alleleBalanceFilter.run(argc, argv);
}
