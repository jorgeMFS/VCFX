#include "VCFX_quality_adjuster.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <cmath>

// Implementation

int VCFXQualityAdjuster::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string transformationStr;

    static struct option long_options[] = {
        {"help",            no_argument,        0, 'h'},
        {"adjust-qual",     required_argument,  0, 'a'},
        {0,                 0,                  0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                transformationStr = std::string(optarg);
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || transformationStr.empty()) {
        displayHelp();
        return 1;
    }

    // Parse the transformation function
    std::function<double(double)> transFunc;
    if (!parseTransformationFunction(transformationStr, transFunc)) {
        std::cerr << "Error: Unsupported transformation function '" << transformationStr << "'.\n";
        displayHelp();
        return 1;
    }

    // Adjust quality scores
    adjustQualityScores(std::cin, std::cout, transFunc);

    return 0;
}

void VCFXQualityAdjuster::displayHelp() {
    std::cout << "VCFX_quality_adjuster: Adjust quality scores in a VCF file using a specified transformation function.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_quality_adjuster [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                 Display this help message and exit\n";
    std::cout << "  -a, --adjust-qual <FUNC>   Specify the transformation function for QUAL scores (e.g., log, sqrt, square, identity)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_quality_adjuster --adjust-qual log < input.vcf > adjusted_quality.vcf\n";
}

bool VCFXQualityAdjuster::parseTransformationFunction(const std::string& funcStr, std::function<double(double)> &transFunc) {
    // Check if the function string matches a supported function
    if (supportedFunctions.find(funcStr) != supportedFunctions.end()) {
        transFunc = supportedFunctions[funcStr];
        return true;
    }

    // Attempt to parse functions with parameters (e.g., custom scaling)
    // For simplicity, only predefined functions are supported
    return false;
}

void VCFXQualityAdjuster::adjustQualityScores(std::istream& in, std::ostream& out, std::function<double(double)> transFunc) {
    std::string line;
    while (std::getline(in, line)) {
        // Output header lines as-is
        if (line.empty() || line[0] == '#') {
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;

        // Split the line into fields
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields): " << line << "\n";
            continue;
        }

        // Adjust the QUAL field (6th column, index 5)
        try {
            double qual = std::stod(fields[5]);
            double adjustedQual = transFunc(qual);
            // Ensure QUAL is non-negative
            adjustedQual = std::max(adjustedQual, 0.0);
            fields[5] = std::to_string(adjustedQual);
        } catch (const std::invalid_argument&) {
            std::cerr << "Warning: Invalid QUAL value. Skipping line: " << line << "\n";
            continue;
        }

        // Reconstruct the line
        std::ostringstream oss;
        for (size_t i = 0; i < fields.size(); ++i) {
            oss << fields[i];
            if (i != fields.size() - 1) {
                oss << "\t";
            }
        }
        oss << "\n";
        out << oss.str();
    }
}
