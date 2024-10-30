#include "VCFX_population_filter.h"
#include <getopt.h>
#include <sstream>
#include <unordered_set>
#include <fstream>
#include <algorithm>

// Implementation of VCFXPopulationFilter
int VCFXPopulationFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string populationTag;
    std::string popMapFile;

    static struct option long_options[] = {
        {"help",      no_argument,       0, 'h'},
        {"population", required_argument, 0, 'p'},
        {"pop-map",   required_argument, 0, 'm'},
        {0,           0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hp:m:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'p':
                populationTag = optarg;
                break;
            case 'm':
                popMapFile = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || populationTag.empty() || popMapFile.empty()) {
        displayHelp();
        return 1;
    }

    // Perform population filtering on stdin and output to stdout
    filterPopulation(std::cin, std::cout, populationTag, popMapFile);

    return 0;
}

void VCFXPopulationFilter::displayHelp() {
    std::cout << "VCFX_population_filter: Filter VCF to include only samples from a specified population.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_population_filter --population \"<POP_TAG>\" --pop-map <pop_map_file> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                Display this help message and exit\n";
    std::cout << "  -p, --population <POP_TAG> Specify the population tag to filter (e.g., EUR, AFR)\n";
    std::cout << "  -m, --pop-map <file>      Path to population mapping file (format: sample\tpopulation)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_population_filter --population \"EUR\" --pop-map populations.txt < input.vcf > filtered.vcf\n";
}

void VCFXPopulationFilter::filterPopulation(std::istream& in, std::ostream& out, const std::string& populationTag, const std::string& popMapFile) {
    std::unordered_set<std::string> samplesToInclude;
    std::ifstream popMap(popMapFile);
    if (!popMap.is_open()) {
        std::cerr << "Error: Unable to open population mapping file: " << popMapFile << "\n";
        return;
    }

    std::string line;
    while (std::getline(popMap, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string sample, pop;
        if (!std::getline(ss, sample, '\t') || !std::getline(ss, pop, '\t')) {
            std::cerr << "Warning: Invalid line in population mapping file: " << line << "\n";
            continue;
        }
        if (pop == populationTag) {
            samplesToInclude.insert(sample);
        }
    }
    popMap.close();

    if (samplesToInclude.empty()) {
        std::cerr << "Warning: No samples found for population tag: " << populationTag << "\n";
    }

    // Process VCF
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    std::vector<int> sampleIndices;

    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
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

                // Identify samples to include
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    if (samplesToInclude.find(headerFields[i]) != samplesToInclude.end()) {
                        sampleIndices.push_back(static_cast<int>(i));
                    }
                }

                // Write new header with filtered samples
                std::stringstream newHeader;
                for (size_t i = 0; i < 9; ++i) { // First 9 columns
                    newHeader << headerFields[i] << "\t";
                }
                for (size_t i = 0; i < sampleIndices.size(); ++i) {
                    newHeader << headerFields[sampleIndices[i]];
                    if (i != sampleIndices.size() - 1) {
                        newHeader << "\t";
                    }
                }
                out << newHeader.str() << "\n";
                headerParsed = true;
            } else {
                // Write other header lines as-is
                out << line << "\n";
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse VCF data lines
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fieldsVec;
        while (std::getline(ss, field, '\t')) {
            fieldsVec.push_back(field);
        }

        if (fieldsVec.size() < 9) {
            std::cerr << "Warning: Invalid VCF line with fewer than 9 fields: " << line << "\n";
            continue;
        }

        // Prepare new VCF line with filtered samples
        std::stringstream newLine;
        for (size_t i = 0; i < 9; ++i) { // First 9 columns
            newLine << fieldsVec[i] << "\t";
        }
        for (size_t i = 0; i < sampleIndices.size(); ++i) {
            newLine << fieldsVec[sampleIndices[i]];
            if (i != sampleIndices.size() - 1) {
                newLine << "\t";
            }
        }
        out << newLine.str() << "\n";
    }

    // If no variants were output, ensure at least the header was written
    if (!headerParsed) {
        std::cerr << "Error: No header line found in VCF input.\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXPopulationFilter populationFilter;
    return populationFilter.run(argc, argv);
}
