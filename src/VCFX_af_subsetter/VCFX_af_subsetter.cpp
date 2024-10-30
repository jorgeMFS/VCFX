#include "VCFX_af_subsetter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>

int VCFXAfSubsetter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    double minAF = 0.0;
    double maxAF = 1.0;

    static struct option long_options[] = {
        {"help",      no_argument,       0, 'h'},
        {"af-filter", required_argument, 0, 'a'},
        {0,           0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                {
                    std::string range(optarg);
                    size_t dashPos = range.find('-');
                    if (dashPos == std::string::npos) {
                        std::cerr << "Error: Invalid AF range format. Use <minAF>-<maxAF>.\n";
                        displayHelp();
                        return 1;
                    }
                    try {
                        minAF = std::stod(range.substr(0, dashPos));
                        maxAF = std::stod(range.substr(dashPos + 1));
                        if (minAF < 0.0 || maxAF > 1.0 || minAF > maxAF) {
                            throw std::invalid_argument("AF values out of range or minAF > maxAF.");
                        }
                    } catch (const std::invalid_argument&) {
                        std::cerr << "Error: Invalid AF range values. Ensure they are numbers between 0.0 and 1.0 with minAF <= maxAF.\n";
                        displayHelp();
                        return 1;
                    }
                }
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Perform AF-based subsetting
    subsetByAlleleFrequency(std::cin, std::cout, minAF, maxAF);

    return 0;
}

void VCFXAfSubsetter::displayHelp() {
    std::cout << "VCFX_af_subsetter: Subset variants based on alternate allele frequency (AF) ranges.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_af_subsetter [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help            Display this help message and exit\n";
    std::cout << "  -a, --af-filter <minAF>-<maxAF>  Specify the AF range for filtering (e.g., 0.01-0.05)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_af_subsetter --af-filter 0.01-0.05 < input.vcf > subsetted.vcf\n";
}

bool VCFXAfSubsetter::parseAF(const std::string& infoField, double& af) {
    size_t pos = infoField.find("AF=");
    if (pos == std::string::npos) {
        return false;
    }
    size_t start = pos + 3; // Length of "AF="
    size_t end = infoField.find(';', start);
    std::string afStr = (end == std::string::npos) ? infoField.substr(start) : infoField.substr(start, end - start);
    try {
        af = std::stod(afStr);
    } catch (...) {
        return false;
    }
    return true;
}

void VCFXAfSubsetter::subsetByAlleleFrequency(std::istream& in, std::ostream& out, double minAF, double maxAF) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // Header lines are printed as-is
            out << line << "\n";
            continue;
        }

        // Split the line into fields
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields): " << line << "\n";
            continue;
        }

        std::string info = fields[7];
        double af = 0.0;
        if (!parseAF(info, af)) {
            std::cerr << "Warning: AF not found or invalid in INFO field. Skipping variant: " << line << "\n";
            continue;
        }

        if (af >= minAF && af <= maxAF) {
            // Variant falls within the specified AF range
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXAfSubsetter afSubsetter;
    return afSubsetter.run(argc, argv);
}