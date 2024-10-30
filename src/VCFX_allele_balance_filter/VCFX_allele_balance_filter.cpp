#include "VCFX_allele_balance_filter.h"
#include <getopt.h>
#include <algorithm>
#include <iomanip>

int VCFXAlleleBalanceFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    double threshold = 0.0;

    static struct option long_options[] = {
        {"help",                no_argument,       0, 'h'},
        {"filter-allele-balance", required_argument, 0, 'f'},
        {0,                     0,                 0,  0 }
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

    if (showHelp || threshold < 0.0 || threshold > 1.0) {
        displayHelp();
        return 1;
    }

    // Perform allele balance filtering
    filterByAlleleBalance(std::cin, std::cout, threshold);

    return 0;
}

void VCFXAlleleBalanceFilter::displayHelp() {
    std::cout << "VCFX_allele_balance_filter: Filter VCF variants based on allele balance ratios.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_allele_balance_filter --filter-allele-balance <THRESHOLD> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                          Display this help message and exit\n";
    std::cout << "  -f, --filter-allele-balance <THRESHOLD>  Specify the allele balance threshold (0.0 - 1.0)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_allele_balance_filter --filter-allele-balance 0.3 < input.vcf > filtered.vcf\n";
}

void VCFXAlleleBalanceFilter::filterByAlleleBalance(std::istream& in, std::ostream& out, double threshold) {
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        std::vector<std::string> genotypes;
        std::string genotype;
        while (ss >> genotype) {
            genotypes.push_back(genotype);
        }

        bool pass = true;
        for (const auto& gt : genotypes) {
            double ab = calculateAlleleBalance(gt);
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

double VCFXAlleleBalanceFilter::calculateAlleleBalance(const std::string& genotype) {
    // Assuming genotype format GT:... (e.g., 0/1:...)
    size_t gt_end = genotype.find(':');
    std::string gt = (gt_end != std::string::npos) ? genotype.substr(0, gt_end) : genotype;

    // Count alleles
    int ref_count = 0;
    int alt_count = 0;
    for (char allele : gt) {
        if (allele == '0') {
            ref_count++;
        } else if (allele == '1') {
            alt_count++;
        }
    }

    int total = ref_count + alt_count;
    if (total == 0) {
        return 0.0;
    }

    return static_cast<double>(ref_count) / total;
}

int main(int argc, char* argv[]) {
    VCFXAlleleBalanceFilter alleleBalanceFilter;
    return alleleBalanceFilter.run(argc, argv);
}