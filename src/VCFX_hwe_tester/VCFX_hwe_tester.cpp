#include "VCFX_hwe_tester.h"
#include <getopt.h>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>

// Implementation of VCFX_hwe_tester
int VCFXHWETester::run(int argc, char* argv[]) {
    // Parse command-line arguments (if any in future extensions)
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0,      0,           0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Perform HWE tests on stdin
    performHWE(std::cin);

    return 0;
}

void VCFXHWETester::displayHelp() {
    std::cout << "VCFX_hwe_tester: Perform Hardy-Weinberg Equilibrium (HWE) tests on a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_hwe_tester [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_hwe_tester < input.vcf\n";
}

void VCFXHWETester::performHWE(std::istream& in) {
    std::string line;
    // Print header for output
    std::cout << "CHROM\tPOS\tID\tREF\tALT\tHWE_p-value\n";

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue; // Skip header lines

        // Split the VCF line into fields
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 10) {
            std::cerr << "Invalid VCF line with fewer than 10 fields.\n";
            continue;
        }

        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string format = fields[8];

        // Find the index of the GT field
        std::vector<std::string> format_fields;
        std::stringstream fmt_ss(format);
        std::string fmt_field;
        while (std::getline(fmt_ss, fmt_field, ':')) {
            format_fields.push_back(fmt_field);
        }

        int gt_index = -1;
        for (size_t i = 0; i < format_fields.size(); ++i) {
            if (format_fields[i] == "GT") {
                gt_index = static_cast<int>(i);
                break;
            }
        }

        if (gt_index == -1) {
            std::cerr << "GT field not found in FORMAT column.\n";
            continue;
        }

        // Collect all genotype information
        std::vector<std::string> genotypes;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_fields;
            std::stringstream samp_ss(fields[i]);
            std::string samp_field;
            while (std::getline(samp_ss, samp_field, ':')) {
                sample_fields.push_back(samp_field);
            }

            if (gt_index >= static_cast<int>(sample_fields.size())) {
                std::cerr << "GT index out of range in sample fields.\n";
                continue;
            }

            std::string genotype = sample_fields[gt_index];
            genotypes.push_back(genotype);
        }

        // Parse genotypes to count homozygotes and heterozygotes
        int homRef = 0;
        int het = 0;
        int homAlt = 0;
        bool valid = parseGenotypes(genotypes, homRef, het, homAlt);
        if (!valid) {
            std::cerr << "Error parsing genotypes for variant at position " << chrom << ":" << pos << "\n";
            continue;
        }

        // Calculate HWE p-value
        double pValue = calculateHWE(homRef, het, homAlt);

        // Output the result
        std::cout << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" 
                  << std::fixed << std::setprecision(5) << pValue << "\n";
    }
}

bool VCFXHWETester::parseGenotypes(const std::vector<std::string>& genotypes, int& homRef, int& het, int& homAlt) {
    homRef = 0;
    het = 0;
    homAlt = 0;

    for (const auto& gt : genotypes) {
        if (gt.empty() || gt == "." || gt == "./." || gt == ".|.") continue;

        std::string gt_clean = gt;
        // Replace '|' with '/' for consistency
        std::replace(gt_clean.begin(), gt_clean.end(), '|', '/');

        // Split genotype into alleles
        std::vector<std::string> alleles;
        std::stringstream ss(gt_clean);
        std::string allele;
        while (std::getline(ss, allele, '/')) {
            alleles.push_back(allele);
        }

        if (alleles.size() != 2) {
            // Not a diploid genotype
            return false;
        }

        if (alleles[0] == "0" && alleles[1] == "0") {
            homRef++;
        } else if (alleles[0] == "1" && alleles[1] == "1") {
            homAlt++;
        } else if ((alleles[0] == "0" && alleles[1] == "1") || (alleles[0] == "1" && alleles[1] == "0")) {
            het++;
        } else {
            // Invalid allele
            return false;
        }
    }

    return true;
}

double VCFXHWETester::calculateHWE(int homRef, int het, int homAlt) {
    int N = homRef + het + homAlt;
    int x = 2 * homAlt + het; // Count of alternate alleles
    int y = 2 * homRef + het; // Count of reference alleles

    // Ensure the site is biallelic
    if (x + y != 2 * N) {
        return 1.0; // Return 1.0 for non-biallelic sites
    }

    return genotypeProbability(homRef, het, homAlt);
}

double VCFXHWETester::genotypeProbability(int homRef, int het, int homAlt) {
    int N = homRef + het + homAlt;
    int x = 2 * homAlt + het; // Count of alternate alleles
    int y = 2 * homRef + het; // Count of reference alleles

    // Initialize variables for exact test
    int n = y; // Number of ref alleles
    int N_total = 2 * N;

    // Initialize log factorials
    // Precompute log factorials up to 2N
    std::vector<double> logFac(2 * N + 1, 0.0);
    for (int i = 1; i <= 2 * N; ++i) {
        logFac[i] = logFac[i - 1] + std::log(static_cast<double>(i));
    }

    // Function to calculate log(n!), using precomputed values
    auto logFactorial = [&](int n) -> double {
        if (n < 0 || n > 2 * N) return 0.0; // Undefined
        return logFac[n];
    };

    // Function to calculate the log probability of a genotype configuration
    auto logProb = [&](int a) -> double {
        // a = number of heterozygotes
        // number of homozygous refs = (y - a) / 2
        // number of homozygous alts = (x - a) / 2
        if ((y - a) < 0 || (x - a) < 0) return -INFINITY;
        if ((y - a) % 2 != 0 || (x - a) % 2 != 0) return -INFINITY;

        int homRef = (y - a) / 2;
        int homAlt = (x - a) / 2;

        // Multinomial coefficient: (N)! / (n_AA! * n_Aa! * n_Aa!)
        double logCoef = logFactorial(N) - logFactorial(homRef) - logFactorial(homAlt) - logFactorial(a);

        // Probability under HWE: (p^2)^homRef * (2pq)^het * (q^2)^homAlt
        // Taking log:
        // 2*log(p) * homRef + log(2) + log(p) + log(q) * het + 2*log(q) * homAlt
        // Simplified as:
        // 2*homRef*log(p) + het*log(2*p*q) + 2*homAlt*log(q)
        // Where p = y / (y + x), q = x / (y + x)
        double p = static_cast<double>(y) / (y + x);
        double q = static_cast<double>(x) / (y + x);

        if (p == 0.0 || q == 0.0) return -INFINITY;

        double logP = 2.0 * homRef * std::log(p) + a * std::log(2.0 * p * q) + 2.0 * homAlt * std::log(q);

        return logCoef + logP;
    };

    // Calculate observed heterozygotes
    int observedHet = het;

    // Calculate log probability of observed genotype
    double logProbObs = logProb(observedHet);

    // Find minimum log probability to avoid underflow
    double minLogProb = logProbObs;
    for (int a = 0; a <= std::min(y, x); ++a) {
        double lp = logProb(a);
        if (lp < minLogProb && lp != -INFINITY) {
            minLogProb = lp;
        }
    }

    // Calculate the exponents relative to minLogProb for numerical stability
    double sum = 0.0;
    double probObs = 0.0;

    for (int a = 0; a <= std::min(y, x); ++a) {
        double lp = logProb(a);
        if (lp == -INFINITY) continue;

        double relative = lp - minLogProb;
        double prob = std::exp(relative);
        sum += prob;

        if (a == observedHet) {
            probObs = prob;
        }
    }

    // Normalize
    double probNormalized = probObs / sum;

    // Calculate p-value: sum of probabilities <= probObs
    double pValue = 0.0;
    for (int a = 0; a <= std::min(y, x); ++a) {
        double lp = logProb(a);
        if (lp == -INFINITY) continue;

        double relative = lp - minLogProb;
        double prob = std::exp(relative);

        if (prob <= probObs + 1e-8) { // Added a small epsilon for floating point precision
            pValue += prob;
        }
    }

    // Normalize p-value
    pValue /= sum;

    // Ensure p-value is between 0 and 1
    if (pValue > 1.0) pValue = 1.0;
    if (pValue < 0.0) pValue = 0.0;

    return pValue;
}

int main(int argc, char* argv[]) {
    VCFXHWETester tester;
    return tester.run(argc, argv);
}
