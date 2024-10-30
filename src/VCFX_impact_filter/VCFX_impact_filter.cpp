#include "VCFX_impact_filter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <regex>

// Implementation of VCFXImpactFilter
int VCFXImpactFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string targetImpact;

    static struct option long_options[] = {
        {"help",        no_argument,       0, 'h'},
        {"filter-impact", required_argument, 0, 'i'},
        {0,             0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hi:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'i':
                targetImpact = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || targetImpact.empty()) {
        displayHelp();
        return 1;
    }

    // Perform impact filtering on stdin and output to stdout
    filterByImpact(std::cin, std::cout, targetImpact);

    return 0;
}

void VCFXImpactFilter::displayHelp() {
    std::cout << "VCFX_impact_filter: Filter VCF variants based on predicted impact from annotations.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_impact_filter --filter-impact \"<IMPACT_LEVEL>\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                 Display this help message and exit\n";
    std::cout << "  -i, --filter-impact <level> Specify the impact level to filter (e.g., HIGH, MODERATE)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_impact_filter --filter-impact \"HIGH\" < input.vcf > filtered.vcf\n";
}

void VCFXImpactFilter::filterByImpact(std::istream& in, std::ostream& out, const std::string& targetImpact) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    size_t impactIndex = std::string::npos;

    // Define possible impact levels in hierarchy
    enum ImpactLevel { UNKNOWN, MODIFIER, LOW, MODERATE, HIGH };
    ImpactLevel targetLevel;

    if (targetImpact == "HIGH") {
        targetLevel = HIGH;
    } else if (targetImpact == "MODERATE") {
        targetLevel = MODERATE;
    } else if (targetImpact == "LOW") {
        targetLevel = LOW;
    } else if (targetImpact == "MODIFIER") {
        targetLevel = MODIFIER;
    } else {
        std::cerr << "Error: Invalid impact level \"" << targetImpact << "\". Choose from HIGH, MODERATE, LOW, MODIFIER.\n";
        return;
    }

    // Regex to extract Impact=LEVEL from INFO field
    std::regex impactRegex(R"(Impact=([A-Z]+))");

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Handle header lines
            if (line.substr(0, 6) == "#CHROM") {
                // Parse header to find the index of 'Impact' in INFO field if necessary
                out << line << "\tREF_IMPACT\n"; // Add a new INFO field for filtered impact
                headerParsed = true;
            } else {
                // Other header lines
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

        if (fieldsVec.size() < 8) {
            std::cerr << "Warning: Invalid VCF line with fewer than 8 fields: " << line << "\n";
            continue;
        }

        std::string infoField = fieldsVec[7];
        std::smatch matches;
        std::string impactValue = "UNKNOWN";

        if (std::regex_search(infoField, matches, impactRegex)) {
            impactValue = matches[1];
        }

        // Determine the impact level
        ImpactLevel variantLevel;
        if (impactValue == "HIGH") {
            variantLevel = HIGH;
        } else if (impactValue == "MODERATE") {
            variantLevel = MODERATE;
        } else if (impactValue == "LOW") {
            variantLevel = LOW;
        } else if (impactValue == "MODIFIER") {
            variantLevel = MODIFIER;
        } else {
            variantLevel = UNKNOWN;
        }

        // Define hierarchy: HIGH > MODERATE > LOW > MODIFIER
        // Include variants >= targetImpact level
        bool includeVariant = false;
        switch (targetLevel) {
            case HIGH:
                if (variantLevel == HIGH) includeVariant = true;
                break;
            case MODERATE:
                if (variantLevel == HIGH || variantLevel == MODERATE) includeVariant = true;
                break;
            case LOW:
                if (variantLevel == HIGH || variantLevel == MODERATE || variantLevel == LOW) includeVariant = true;
                break;
            case MODIFIER:
                includeVariant = true; // Include all
                break;
            default:
                includeVariant = false;
        }

        if (includeVariant) {
            out << line << "\t" << impactValue << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXImpactFilter impactFilter;
    return impactFilter.run(argc, argv);
}   