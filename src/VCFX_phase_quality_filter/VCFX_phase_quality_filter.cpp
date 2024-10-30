#include "VCFX_phase_quality_filter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>

// Implementation of VCFXPhaseQualityFilter
int VCFXPhaseQualityFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string condition;

    static struct option long_options[] = {
        {"help",                no_argument,       0, 'h'},
        {"filter-pq",           required_argument, 0, 'f'},
        {0,                     0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hf:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                condition = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || condition.empty()) {
        displayHelp();
        return 1;
    }

    // Parse the condition (e.g., PQ>30)
    double threshold = 0.0;
    char op;
    if (sscanf(condition.c_str(), "PQ%c%lf", &op, &threshold) != 2) {
        std::cerr << "Error: Invalid condition format. Use format PQ>30\n";
        displayHelp();
        return 1;
    }

    // Perform PQ filtering
    filterByPQ(std::cin, std::cout, threshold);

    return 0;
}

void VCFXPhaseQualityFilter::displayHelp() {
    std::cout << "VCFX_phase_quality_filter: Filter VCF variants based on phasing quality scores.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_phase_quality_filter --filter-pq \"<CONDITION>\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                  Display this help message and exit\n";
    std::cout << "  -f, --filter-pq \"<CONDITION>\" Specify the PQ condition (e.g., PQ>30)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_phase_quality_filter --filter-pq \"PQ>30\" < input.vcf > filtered.vcf\n";
}

void VCFXPhaseQualityFilter::filterByPQ(std::istream& in, std::ostream& out, double threshold) {
    std::string line;
    bool headerPassed = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        double pqScore = parsePQScore(info);
        if (pqScore >= threshold) {
            out << line << "\n";
        }
    }
}

double VCFXPhaseQualityFilter::parsePQScore(const std::string& infoField) {
    // INFO field contains key=value pairs separated by semicolons
    std::stringstream ss(infoField);
    std::string token;
    while (std::getline(ss, token, ';')) {
        if (token.find("PQ=") == 0) {
            std::string valueStr = token.substr(3); // Remove "PQ="
            try {
                return std::stod(valueStr);
            } catch (const std::invalid_argument&) {
                std::cerr << "Warning: Invalid PQ score \"" << valueStr << "\". Treating as 0.\n";
                return 0.0;
            }
        }
    }
    // If PQ not found, treat as 0
    return 0.0;
}

int main(int argc, char* argv[]) {
    VCFXPhaseQualityFilter phaseQualityFilter;
    return phaseQualityFilter.run(argc, argv);
} 