#include "VCFX_haplotype_phaser.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <numeric>

int VCFXHaplotypePhaser::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    double ldThreshold = 0.8; // Default LD threshold

    static struct option long_options[] = {
        {"help",          no_argument,       0, 'h'},
        {"ld-threshold", required_argument, 0, 'l'},
        {0,               0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hl:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'l':
                try {
                    ldThreshold = std::stod(optarg);
                } catch (const std::invalid_argument&) {
                    std::cerr << "Error: Invalid LD threshold value.\n";
                    displayHelp();
                    return 1;
                }
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || ldThreshold < 0.0 || ldThreshold > 1.0) {
        displayHelp();
        return 1;
    }

    // Perform haplotype phasing
    phaseHaplotypes(std::cin, std::cout, ldThreshold);

    return 0;
}

void VCFXHaplotypePhaser::displayHelp() {
    std::cout << "VCFX_haplotype_phaser: Group variants into phased haplotype blocks based on linkage disequilibrium.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_haplotype_phaser [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                    Display this help message and exit\n";
    std::cout << "  -l, --ld-threshold <VALUE>    Specify the LD threshold for grouping (0.0 - 1.0, default: 0.8)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_haplotype_phaser --ld-threshold 0.9 < input.vcf > phased_blocks.txt\n";
}

void VCFXHaplotypePhaser::phaseHaplotypes(std::istream& in, std::ostream& out, double ldThreshold) {
    std::string line;
    std::vector<std::string> sampleNames;
    bool headerParsed = false;
    std::vector<std::vector<int>> genotypeMatrix; // Rows: Variants, Columns: Samples

    // Read header to get sample names
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line.substr(0, 6) == "#CHROM") {
            std::stringstream ss(line);
            std::string field;
            std::vector<std::string> headers;
            while (std::getline(ss, field, '\t')) {
                headers.push_back(field);
            }
            // Sample names start from the 10th column
            for (size_t i = 9; i < headers.size(); ++i) {
                sampleNames.push_back(headers[i]);
            }
            headerParsed = true;
            // Output header
            out << line << "\n";
            break;
        }
    }

    if (!headerParsed) {
        std::cerr << "Error: VCF header line with #CHROM not found.\n";
        return;
    }

    // Read and store all variants
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

        std::vector<int> genotypes;
        std::string genotype;
        while (ss >> genotype) {
            // Parse genotype (assuming GT field first)
            size_t sep = genotype.find_first_of("/|");
            if (sep == std::string::npos) {
                genotypes.push_back(-1); // Missing or unknown genotype
                continue;
            }
            std::string gt1 = genotype.substr(0, sep);
            std::string gt2 = genotype.substr(sep + 1);
            int allele1 = (gt1.empty()) ? -1 : std::stoi(gt1);
            int allele2 = (gt2.empty()) ? -1 : std::stoi(gt2);
            // For simplicity, store count of reference alleles
            if (allele1 == -1 || allele2 == -1) {
                genotypes.push_back(-1); // Missing genotype
            } else {
                genotypes.push_back(allele1 + allele2);
            }
        }

        genotypeMatrix.push_back(genotypes);
    }

    if (genotypeMatrix.empty()) {
        std::cerr << "Error: No variant data found in VCF.\n";
        return;
    }

    // Group variants into haplotype blocks based on LD
    std::vector<std::vector<int>> haplotypeBlocks = groupVariants(genotypeMatrix, ldThreshold);

    // Output haplotype blocks
    out << "#HAPLOTYPE_BLOCKS_START\n";
    for (size_t i = 0; i < haplotypeBlocks.size(); ++i) {
        out << "Block " << (i + 1) << ": ";
        for (size_t j = 0; j < haplotypeBlocks[i].size(); ++j) {
            out << haplotypeBlocks[i][j];
            if (j != haplotypeBlocks[i].size() - 1) {
                out << ", ";
            }
        }
        out << "\n";
    }
    out << "#HAPLOTYPE_BLOCKS_END\n";
}

std::vector<std::vector<int>> VCFXHaplotypePhaser::groupVariants(const std::vector<std::vector<int>>& genotypeMatrix, double ldThreshold) {
    std::vector<std::vector<int>> blocks;
    std::vector<int> currentBlock;

    for (size_t i = 0; i < genotypeMatrix.size(); ++i) {
        if (currentBlock.empty()) {
            currentBlock.push_back(i);
        } else {
            // Calculate LD between the last variant in the block and the current variant
            int lastVariantIndex = currentBlock.back();
            double ld = calculateLD(genotypeMatrix[lastVariantIndex], genotypeMatrix[i]);

            if (ld >= ldThreshold) {
                currentBlock.push_back(i);
            } else {
                blocks.push_back(currentBlock);
                currentBlock.clear();
                currentBlock.push_back(i);
            }
        }
    }

    if (!currentBlock.empty()) {
        blocks.push_back(currentBlock);
    }

    return blocks;
}

double VCFXHaplotypePhaser::calculateLD(const std::vector<int>& variant1, const std::vector<int>& variant2) {
    // Calculate r² linkage disequilibrium between two variants
    int n = 0;
    int sum_X = 0, sum_Y = 0;
    int sum_XY = 0;
    int sum_X2 = 0, sum_Y2 = 0;

    for (size_t i = 0; i < variant1.size(); ++i) {
        int x = variant1[i];
        int y = variant2[i];
        if (x == -1 || y == -1) {
            continue; // Skip missing genotypes
        }

        n++;
        sum_X += x;
        sum_Y += y;
        sum_XY += x * y;
        sum_X2 += x * x;
        sum_Y2 += y * y;
    }

    if (n == 0) {
        return 0.0; // No data to calculate LD
    }

    double mean_X = static_cast<double>(sum_X) / n;
    double mean_Y = static_cast<double>(sum_Y) / n;

    double cov = (static_cast<double>(sum_XY) / n) - (mean_X * mean_Y);
    double var_X = (static_cast<double>(sum_X2) / n) - (mean_X * mean_X);
    double var_Y = (static_cast<double>(sum_Y2) / n) - (mean_Y * mean_Y);

    if (var_X == 0.0 || var_Y == 0.0) {
        return 0.0; // Avoid division by zero
    }

    double r = cov / (std::sqrt(var_X) * std::sqrt(var_Y));
    double r_squared = r * r;

    return r_squared;
}

int main(int argc, char* argv[]) {
    VCFXHaplotypePhaser haplotypePhaser;
    return haplotypePhaser.run(argc, argv);
}