#include "VCFX_ancestry_inferrer.h"
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>

// Implementation of VCFXAncestryInferer
int VCFXAncestryInferer::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string freqFilePath;

    static struct option long_options[] = {
        {"help",       no_argument,       0, 'h'},
        {"frequency",  required_argument, 0, 'f'},
        {0,            0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hf:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                freqFilePath = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || freqFilePath.empty()) {
        displayHelp();
        return 1;
    }

    // Load population allele frequencies
    if (!loadPopulationFrequencies(freqFilePath)) {
        std::cerr << "Error: Failed to load population frequencies from " << freqFilePath << "\n";
        return 1;
    }

    // Perform ancestry inference
    inferAncestry(std::cin, std::cout);

    return 0;
}

void VCFXAncestryInferer::displayHelp() {
    std::cout << "VCFX_ancestry_inferrer: Infer population ancestry based on allele frequencies.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_ancestry_inferrer --frequency <freq_file> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                   Display this help message and exit\n";
    std::cout << "  -f, --frequency <file>       Path to population allele frequencies file\n\n";
    std::cout << "Frequency File Format:\n";
    std::cout << "  Each line should contain: chromosome, position, reference allele, alternate allele, population, frequency\n";
    std::cout << "  Fields are tab-separated. Example:\n";
    std::cout << "  chr1\t12345\tA\tG\tEUR\t0.45\n";
    std::cout << "  chr1\t12345\tA\tG\tAFR\t0.30\n\n";
    std::cout << "Example Command:\n";
    std::cout << "  VCFX_ancestry_inferrer --frequency pop_frequencies.txt < input.vcf > ancestry_results.txt\n";
}

bool VCFXAncestryInferer::loadPopulationFrequencies(const std::string& freqFilePath) {
    std::ifstream freqFile(freqFilePath);
    if (!freqFile.is_open()) {
        std::cerr << "Error: Cannot open frequency file: " << freqFilePath << "\n";
        return false;
    }

    std::string line;
    while (std::getline(freqFile, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string chrom, pos, ref, alt, population, freqStr;
        if (!std::getline(ss, chrom, '\t') ||
            !std::getline(ss, pos, '\t') ||
            !std::getline(ss, ref, '\t') ||
            !std::getline(ss, alt, '\t') ||
            !std::getline(ss, population, '\t') ||
            !std::getline(ss, freqStr, '\t')) {
            std::cerr << "Warning: Invalid line in frequency file: " << line << "\n";
            continue;
        }

        double freq = 0.0;
        try {
            freq = std::stod(freqStr);
        } catch (...) {
            std::cerr << "Warning: Invalid frequency value in line: " << line << "\n";
            continue;
        }

        // Find if the population already exists
        auto it = std::find_if(populations.begin(), populations.end(),
            [&](const PopulationFrequencies& pf) { return pf.population == population; });

        if (it == populations.end()) {
            // Create new population entry
            PopulationFrequencies pf;
            pf.population = population;
            pf.variantFrequencies["" + chrom + ":" + pos + ":" + ref + ":" + alt] = freq;
            populations.push_back(pf);
        } else {
            // Add to existing population
            it->variantFrequencies["" + chrom + ":" + pos + ":" + ref + ":" + alt] = freq;
        }
    }

    freqFile.close();

    if (populations.empty()) {
        std::cerr << "Error: No valid population frequencies loaded.\n";
        return false;
    }

    return true;
}

void VCFXAncestryInferer::inferAncestry(std::istream& vcfInput, std::ostream& ancestryOutput) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    std::vector<std::string> sampleNames;

    // To store sample genotype counts per variant
    std::unordered_map<std::string, std::unordered_map<std::string, int>> sampleAlleleCounts; // sample -> population -> count
    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores; // sample -> population -> score

    while (std::getline(vcfInput, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.substr(0, 6) == "#CHROM") {
                std::stringstream ss(line);
                std::string field;
                // Extract header fields
                while (std::getline(ss, field, '\t')) {
                    headerFields.push_back(field);
                }

                // Extract sample names
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    sampleNames.push_back(headerFields[i]);
                    // Initialize scores
                    for (const auto& pf : populations) {
                        sampleScores[headerFields[i]][pf.population] = 0.0;
                    }
                }

                headerParsed = true;
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse VCF data lines
        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        // Get first 9 columns
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Invalid VCF line with fewer than 9 fields: " << line << "\n";
            continue;
        }

        // Split ALT alleles
        std::vector<std::string> altAlleles;
        std::stringstream alt_ss(alt);
        std::string allele;
        while (std::getline(alt_ss, allele, ',')) {
            altAlleles.push_back(allele);
        }

        // Split FORMAT fields
        std::vector<std::string> formatFields;
        std::stringstream fmt_ss(format);
        std::string fmt_field;
        while (std::getline(fmt_ss, fmt_field, ':')) {
            formatFields.push_back(fmt_field);
        }

        // Find the index of GT (Genotype)
        int gt_index = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gt_index = static_cast<int>(i);
                break;
            }
        }

        if (gt_index == -1) {
            std::cerr << "Warning: GT field not found in FORMAT column.\n";
            continue;
        }

        // For each sample, get the genotype
        for (size_t i = 0; i < sampleNames.size(); ++i) {
            std::string sample = "";
            std::getline(ss, sample, '\t');
            if (sample.empty()) {
                continue;
            }

            std::vector<std::string> sampleFields;
            std::stringstream samp_ss(sample);
            std::string samp_field;
            while (std::getline(samp_ss, samp_field, ':')) {
                sampleFields.push_back(samp_field);
            }

            if (gt_index >= static_cast<int>(sampleFields.size())) {
                std::cerr << "Warning: GT index out of range in sample fields for sample " << sampleNames[i] << ".\n";
                continue;
            }

            std::string genotype = sampleFields[gt_index];
            // Replace '|' with '/' for consistency
            std::replace(genotype.begin(), genotype.end(), '|', '/');

            if (genotype.empty() || genotype == ".") {
                // Missing genotype
                continue;
            }

            std::vector<std::string> allelesList;
            std::stringstream genotype_ss(genotype);
            std::string alleleNum;
            while (std::getline(genotype_ss, alleleNum, '/')) {
                allelesList.push_back(alleleNum);
            }

            if (allelesList.size() != 2) {
                // Not diploid
                continue;
            }

            // Construct variant key
            std::string variantKey = chrom + ":" + pos + ":" + ref + ":" + (altAlleles.empty() ? "NA" : altAlleles[0]);

            // Iterate over populations and accumulate scores based on allele frequencies
            for (const auto& pf : populations) {
                auto it = pf.variantFrequencies.find(variantKey);
                if (it != pf.variantFrequencies.end()) {
                    double freq = it->second;
                    // Assuming allele encoding as 0 for ref, 1 for alt
                    // Genotype: 0/0 -> 0, 0/1 -> 1, 1/1 -> 2
                    int genotypeScore = 0;
                    try {
                        int allele1 = std::stoi(allelesList[0]);
                        int allele2 = std::stoi(allelesList[1]);
                        genotypeScore = allele1 + allele2;
                    } catch (...) {
                        // Invalid genotype
                        continue;
                    }

                    // Update score (could be weighted by frequency)
                    sampleScores[sampleNames[i]][pf.population] += genotypeScore * freq;
                }
            }
        }
    }

    // Infer ancestry for each sample
    ancestryOutput << "Sample\tInferred_Population\n";
    for (const auto& sample : sampleNames) {
        if (sampleScores.find(sample) == sampleScores.end()) {
            ancestryOutput << sample << "\tUnknown\n";
            continue;
        }

        const auto& scores = sampleScores[sample];
        std::string bestPopulation = "Unknown";
        double maxScore = -1e9;

        for (const auto& [pop, score] : scores) {
            if (score > maxScore) {
                maxScore = score;
                bestPopulation = pop;
            }
        }

        ancestryOutput << sample << "\t" << bestPopulation << "\n";
    }
}
