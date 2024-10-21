#include "VCFX_inbreeding_calculator.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

// Implementation

int VCFXInbreedingCalculator::run(int argc, char* argv[]) {
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

    // Calculate inbreeding coefficients
    calculateInbreedingCoefficients(std::cin, std::cout);

    return 0;
}

void VCFXInbreedingCalculator::displayHelp() {
    std::cout << "VCFX_inbreeding_calculator: Calculate inbreeding coefficients (F-statistics) for each individual in a population based on VCF genotypes.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_inbreeding_calculator [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help               Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_inbreeding_calculator < input.vcf > inbreeding_coefficients.txt\n";
}

bool VCFXInbreedingCalculator::parseGenotype(const std::string& genotype, int& a1, int& a2) {
    // Parses genotype in the form GT:... e.g., 0/1, 1/1, 0/0, ./.
    if (genotype.empty()) return false;

    size_t sep = genotype.find(':');
    std::string gt = (sep != std::string::npos) ? genotype.substr(0, sep) : genotype;

    // Handle missing data
    if (gt.find('.') != std::string::npos) {
        return false;
    }

    // Split alleles
    size_t slash = gt.find('/');
    size_t pipe = gt.find('|');
    if (slash != std::string::npos) {
        a1 = std::stoi(gt.substr(0, slash));
        a2 = std::stoi(gt.substr(slash + 1));
        return true;
    }
    else if (pipe != std::string::npos) {
        a1 = std::stoi(gt.substr(0, pipe));
        a2 = std::stoi(gt.substr(pipe + 1));
        return true;
    }

    return false;
}

double VCFXInbreedingCalculator::calculateExpectedHet(int totalAlleles, int homRef, int homAlt, int het) {
    int refAlleles = (homRef * 2) + het;
    int altAlleles = (homAlt * 2) + het;
    double p = static_cast<double>(refAlleles) / totalAlleles;
    double q = static_cast<double>(altAlleles) / totalAlleles;
    return 2.0 * p * q;
}

double VCFXInbreedingCalculator::calculateF(int homozygous, int heterozygous, double expectedHet) {
    if (expectedHet == 0.0) return 0.0;
    return 1.0 - static_cast<double>(heterozygous) / expectedHet;
}

void VCFXInbreedingCalculator::calculateInbreedingCoefficients(std::istream& in, std::ostream& out) {
    std::string line;
    std::vector<std::string> sampleNames;
    bool headerParsed = false;

    // Initialize counts
    std::unordered_map<std::string, int> homozygousCounts;
    std::unordered_map<std::string, int> heterozygousCounts;
    std::unordered_map<std::string, int> totalCounts;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

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
                homozygousCounts[headers[i]] = 0;
                heterozygousCounts[headers[i]] = 0;
                totalCounts[headers[i]] = 0;
            }
            headerParsed = true;
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        // Read genotype fields if present
        std::vector<std::string> genotypes;
        std::string genotype;
        while (ss >> genotype) {
            genotypes.push_back(genotype);
        }

        for (size_t i = 0; i < genotypes.size() && i < sampleNames.size(); ++i) {
            int a1, a2;
            if (parseGenotype(genotypes[i], a1, a2)) {
                if (a1 == a2) {
                    homozygousCounts[sampleNames[i]] += 1;
                }
                else {
                    heterozygousCounts[sampleNames[i]] += 1;
                }
                totalCounts[sampleNames[i]] += 1;
            }
        }
    }

    // Calculate F-statistics
    out << "Sample\tInbreeding_Coefficient(F)\n";
    for (const auto& sample : sampleNames) {
        int homo = homozygousCounts[sample];
        int het = heterozygousCounts[sample];
        int total = totalCounts[sample];
        if (total == 0) {
            out << sample << "\tNA\n";
            continue;
        }
        int totalAlleles = total * 2;
        int homRef = homo;
        int homAlt = homo; // Assuming homozygous can be either homRef or homAlt; lacking information
        double expectedHet = calculateExpectedHet(totalAlleles, homRef, homAlt, het);
        double F = calculateF(homo, het, expectedHet);
        out << sample << "\t" << F << "\n";
    }
}
