#include "VCFX_cross_sample_concordance.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <set>

// Implementation of VCFXCrossSampleConcordance
int VCFXCrossSampleConcordance::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help", no_argument,       0, 'h'},
        {0,      0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Calculate concordance
    calculateConcordance(std::cin, std::cout);

    return 0;
}

void VCFXCrossSampleConcordance::displayHelp() {
    std::cout << "VCFX_cross_sample_concordance: Check variant concordance across multiple samples.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_cross_sample_concordance [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_cross_sample_concordance < input.vcf > concordance_results.txt\n";
}

void VCFXCrossSampleConcordance::calculateConcordance(std::istream& in, std::ostream& out) {
    std::string line;
    std::vector<std::string> sampleNames;
    bool headerParsed = false;

    // Read header to get sample names
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line.substr(0, 6) == "#CHROM") {
            std::stringstream ss(line);
            std::string field;
            // Parse header fields
            std::vector<std::string> headers;
            while (std::getline(ss, field, '\t')) {
                headers.push_back(field);
            }
            // Sample names start from the 10th column
            for (size_t i = 9; i < headers.size(); ++i) {
                sampleNames.push_back(headers[i]);
            }
            headerParsed = true;
            break;
        }
    }

    if (!headerParsed) {
        std::cerr << "Error: VCF header line with #CHROM not found.\n";
        return;
    }

    // Initialize concordance counts
    size_t totalVariants = 0;
    size_t concordantVariants = 0;
    size_t discordantVariants = 0;

    // Read and process each variant
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
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
        while (std::getline(ss, genotype, '\t')) {
            genotypes.push_back(genotype);
        }

        if (genotypes.size() != sampleNames.size()) {
            std::cerr << "Warning: Number of genotypes does not match number of samples in line: " << line << "\n";
            continue;
        }

        // Check concordance across samples
        std::set<std::string> uniqueGenotypes(genotypes.begin(), genotypes.end());
        if (uniqueGenotypes.size() == 1) {
            concordantVariants++;
        } else {
            discordantVariants++;
        }

        totalVariants++;
    }

    // Output concordance results
    out << "Total Variants Processed: " << totalVariants << "\n";
    out << "Concordant Variants (All samples agree): " << concordantVariants << "\n";
    out << "Discordant Variants (Samples differ): " << discordantVariants << "\n";
}

int main(int argc, char* argv[]) {
    VCFXCrossSampleConcordance crossSampleConcordance;
    return crossSampleConcordance.run(argc, argv);
}