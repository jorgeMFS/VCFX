#include "VCFX_ancestry_assigner.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

// Implementation

int VCFXAncestryAssigner::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string freqFile = "";

    static struct option long_options[] = {
        {"help",                no_argument,       0, 'h'},
        {"assign-ancestry",     required_argument, 0, 'a'},
        {0,                     0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                freqFile = std::string(optarg);
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || freqFile.empty()) {
        displayHelp();
        return 1;
    }

    // Open frequency file
    std::ifstream freqStream(freqFile);
    if (!freqStream.is_open()) {
        std::cerr << "Error: Unable to open allele frequency file: " << freqFile << "\n";
        return 1;
    }

    // Load ancestral frequencies
    if (!loadAncestralFrequencies(freqStream)) {
        std::cerr << "Error: Failed to load ancestral frequencies.\n";
        return 1;
    }

    // Assign ancestry based on VCF input
    assignAncestry(std::cin, std::cout);

    return 0;
}

void VCFXAncestryAssigner::displayHelp() {
    std::cout << "VCFX_ancestry_assigner: Assign samples to ancestral populations based on variant allele frequencies.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_ancestry_assigner --assign-ancestry <frequency_file> < input.vcf > ancestry_assignments.txt\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                   Display this help message and exit\n";
    std::cout << "  -a, --assign-ancestry <FILE> Specify the ancestral allele frequency file\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_ancestry_assigner --assign-ancestry ancestral_freq.txt < input.vcf > ancestry_assignments.txt\n";
}

bool VCFXAncestryAssigner::parseFrequencyLine(const std::string& line, std::string& chrom, int& pos, char& ref, char& alt, std::unordered_map<std::string, double>& freqMap, const std::vector<std::string>& populations) {
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    // Expecting at least 4 columns: CHROM, POS, REF, ALT, Population1, Population2, ...
    if (fields.size() < 4) {
        return false;
    }

    chrom = fields[0];
    try {
        pos = std::stoi(fields[1]);
    } catch (...) {
        return false;
    }
    if (fields[2].empty() || fields[3].empty()) {
        return false;
    }
    ref = fields[2].at(0);
    alt = fields[3].at(0);

    // Populate frequency map for each population
    for (size_t i = 0; i < populations.size(); ++i) {
        if (4 + i >= fields.size()) {
            freqMap[populations[i]] = 0.0; // Assign 0.0 if frequency is missing
            continue;
        }
        try {
            double freq = std::stod(fields[4 + i]);
            freqMap[populations[i]] = freq;
        } catch (...) {
            freqMap[populations[i]] = 0.0; // Assign 0.0 if conversion fails
        }
    }

    return true;
}

bool VCFXAncestryAssigner::loadAncestralFrequencies(std::istream& in) {
    std::string line;
    // Assume the first line is a header with population names
    if (!std::getline(in, line)) {
        std::cerr << "Error: Allele frequency file is empty.\n";
        return false;
    }

    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> headers;

    // Parse header to get population names
    while (std::getline(ss, field, '\t')) {
        headers.push_back(field);
    }

    // Check if there are at least 5 columns
    if (headers.size() < 5) {
        std::cerr << "Error: Allele frequency file header must contain at least 5 columns (CHROM, POS, REF, ALT, Populations...).\n";
        return false;
    }

    // Extract population names starting from the 5th column
    for (size_t i = 4; i < headers.size(); ++i) {
        populations.push_back(headers[i]);
    }

    // Parse frequency data
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ssLine(line);
        std::string chrom;
        int pos;
        char ref, alt;
        std::unordered_map<std::string, double> freqMap;

        if (!parseFrequencyLine(line, chrom, pos, ref, alt, freqMap, populations)) {
            std::cerr << "Warning: Skipping invalid frequency line: " << line << "\n";
            continue;
        }

        // Key for variant
        std::string key = chrom + ":" + std::to_string(pos) + ":" + std::string(1, ref) + ":" + std::string(1, alt);
        variantFrequencies[key] = freqMap;
    }

    return true;
}

void VCFXAncestryAssigner::assignAncestry(std::istream& vcfIn, std::ostream& out) {
    std::string line;
    bool headerParsed = false;
    int chrIndex = -1, posIndex = -1, refIndex = -1, altIndex = -1;
    std::vector<std::string> sampleNames;

    // Initialize sample ancestry scores: Sample -> Population -> Cumulative Frequency
    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores;
    // Initialize sample variant counts per population
    std::unordered_map<std::string, std::unordered_map<std::string, int>> sampleVariantCounts;

    while (std::getline(vcfIn, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.substr(0, 6) == "#CHROM") {
                // Parse header to get column indices and sample names
                std::stringstream ss(line);
                std::string field;
                std::vector<std::string> headers;
                while (std::getline(ss, field, '\t')) {
                    headers.push_back(field);
                }
                for (size_t i = 0; i < headers.size(); ++i) {
                    if (headers[i] == "CHROM") chrIndex = static_cast<int>(i);
                    else if (headers[i] == "POS") posIndex = static_cast<int>(i);
                    else if (headers[i] == "REF") refIndex = static_cast<int>(i);
                    else if (headers[i] == "ALT") altIndex = static_cast<int>(i);
                }
                if (chrIndex == -1 || posIndex == -1 || refIndex == -1 || altIndex == -1) {
                    std::cerr << "Error: VCF header does not contain required fields (CHROM, POS, REF, ALT).\n";
                    return;
                }
                // Sample names start from the 10th column (index 9)
                for (size_t i = 9; i < headers.size(); ++i) {
                    sampleNames.push_back(headers[i]);
                    // Initialize population scores and counts
                    for (const auto& pop : populations) {
                        sampleScores[headers[i]][pop] = 0.0;
                        sampleVariantCounts[headers[i]][pop] = 0;
                    }
                }
                headerParsed = true;
            }
            // Continue to next line
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;

        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields): " << line << "\n";
            continue;
        }

        std::string chrom = fields[chrIndex];
        int pos = 0;
        try {
            pos = std::stoi(fields[posIndex]);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value. Skipping line: " << line << "\n";
            continue;
        }

        std::string ref = fields[refIndex];
        std::string alt = fields[altIndex];

        // Handle multi-allelic ALT fields
        std::vector<std::string> alts;
        std::stringstream altSS(alt);
        while (std::getline(altSS, field, ',')) {
            alts.push_back(field);
        }

        for (const auto& allele : alts) {
            // Key for variant
            std::string key = chrom + ":" + std::to_string(pos) + ":" + ref + ":" + allele;
            if (variantFrequencies.find(key) == variantFrequencies.end()) {
                // Variant not found in frequency data
                continue;
            }

            // For each population, get frequency
            std::unordered_map<std::string, double> freqMap = variantFrequencies[key];

            // Find the population with the highest frequency for this variant
            std::string assignedPop = "";
            double maxFreq = -1.0;
            for (const auto& pop : populations) {
                if (freqMap.find(pop) != freqMap.end()) {
                    if (freqMap[pop] > maxFreq) {
                        maxFreq = freqMap[pop];
                        assignedPop = pop;
                    }
                }
            }

            if (assignedPop.empty()) {
                continue;
            }

            // Get the GT (genotype) field index
            std::string format = fields[8];
            std::vector<std::string> formatFields;
            std::stringstream formatSS(format);
            while (std::getline(formatSS, field, ':')) {
                formatFields.push_back(field);
            }

            int gtIndex = -1;
            for (size_t i = 0; i < formatFields.size(); ++i) {
                if (formatFields[i] == "GT") {
                    gtIndex = static_cast<int>(i);
                    break;
                }
            }

            if (gtIndex == -1) {
                // GT field not found
                continue;
            }

            // Iterate over each sample
            for (size_t i = 9; i < fields.size() && (i - 9) < sampleNames.size(); ++i) {
                std::vector<std::string> genotypeFields;
                std::stringstream gtSS(fields[i]);
                while (std::getline(gtSS, field, ':')) {
                    genotypeFields.push_back(field);
                }

                if (gtIndex >= static_cast<int>(genotypeFields.size())) {
                    continue;
                }

                std::string genotype = genotypeFields[gtIndex];
                int a1 = 0, a2 = 0;
                // Parse genotype
                size_t sep = genotype.find_first_of("/|");
                if (sep == std::string::npos) {
                    continue;
                }

                try {
                    a1 = std::stoi(genotype.substr(0, sep));
                    a2 = std::stoi(genotype.substr(sep + 1));
                } catch (...) {
                    continue;
                }

                // Update sample scores based on genotype
                if (a1 == 0 && a2 == 0) {
                    // Reference homozygous: No contribution to alternate allele frequency
                } else if ((a1 == 0 && a2 == 1) || (a1 == 1 && a2 == 0)) {
                    // Heterozygous
                    sampleScores[sampleNames[i - 9]][assignedPop] += maxFreq;
                    sampleVariantCounts[sampleNames[i - 9]][assignedPop] += 1;
                } else if (a1 == 1 && a2 == 1) {
                    // Alternate homozygous
                    sampleScores[sampleNames[i - 9]][assignedPop] += 2 * maxFreq;
                    sampleVariantCounts[sampleNames[i - 9]][assignedPop] += 1;
                }
            }
        }
    }

    // Output results
    out << "Sample\tAssigned_Ancestry\n";
    for (const auto& sample : sampleScores) {
        const std::string& sampleName = sample.first;
        const auto& popScores = sample.second;

        std::string assignedPop = "NA";
        double maxScore = -1.0;

        for (const auto& popScore : popScores) {
            if (popScore.second > maxScore) {
                maxScore = popScore.second;
                assignedPop = popScore.first;
            }
        }

        out << sampleName << "\t" << assignedPop << "\n";
    }
}
