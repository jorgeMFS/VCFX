
#include "VCFX_info_aggregator.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <vector>

// Implementation of VCFXInfoAggregator
int VCFXInfoAggregator::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string infoFieldsStr;

    static struct option long_options[] = {
        {"help",          no_argument,       0, 'h'},
        {"aggregate-info", required_argument, 0, 'a'},
        {0,               0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                infoFieldsStr = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || infoFieldsStr.empty()) {
        displayHelp();
        return 1;
    }

    // Split infoFieldsStr by comma to get individual INFO fields
    std::vector<std::string> infoFields;
    std::stringstream ss(infoFieldsStr);
    std::string field;
    while (std::getline(ss, field, ',')) {
        // Trim whitespace
        field.erase(field.find_last_not_of(" \n\r\t")+1);
        field.erase(0, field.find_first_not_of(" \n\r\t"));
        if (!field.empty()) {
            infoFields.push_back(field);
        }
    }

    if (infoFields.empty()) {
        std::cerr << "Error: No valid INFO fields specified for aggregation.\n";
        return 1;
    }

    // Perform INFO field aggregation on stdin and output to stdout
    aggregateInfo(std::cin, std::cout, infoFields);

    return 0;
}

void VCFXInfoAggregator::displayHelp() {
    std::cout << "VCFX_info_aggregator: Aggregate numeric values in the INFO field across samples.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_info_aggregator --aggregate-info \"<INFO_FIELDS>\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                       Display this help message and exit\n";
    std::cout << "  -a, --aggregate-info <fields>    Comma-separated list of INFO fields to aggregate (e.g., DP,AF)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_info_aggregator --aggregate-info \"DP,AF\" < input.vcf > aggregated_info.txt\n";
}

void VCFXInfoAggregator::aggregateInfo(std::istream& in, std::ostream& out, const std::vector<std::string>& infoFields) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    std::map<std::string, std::vector<double>> aggregates; // Key: INFO field, Value: list of numeric values

    // Initialize aggregates map
    for (const auto& field : infoFields) {
        aggregates[field] = std::vector<double>();
    }

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Handle header lines
            if (line.substr(0, 6) == "#CHROM") {
                // Output header with additional aggregated fields
                out << line;
                for (const auto& field : infoFields) {
                    out << "\t" << "AGG_" << field;
                }
                out << "\n";
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
        std::string fieldEntry;
        std::vector<std::string> fieldsVec;

        while (std::getline(ss, fieldEntry, '\t')) {
            fieldsVec.push_back(fieldEntry);
        }

        if (fieldsVec.size() < 8) {
            std::cerr << "Warning: Invalid VCF line with fewer than 8 fields: " << line << "\n";
            continue;
        }

        std::string infoField = fieldsVec[7];
        std::map<std::string, double> currentVariantAggregates;

        // Extract specified INFO fields
        for (const auto& aggField : infoFields) {
            size_t pos = infoField.find(aggField + "=");
            if (pos != std::string::npos) {
                size_t start = pos + aggField.length() + 1;
                size_t end = infoField.find(';', start);
                std::string valueStr = (end != std::string::npos) ? infoField.substr(start, end - start) : infoField.substr(start);
                try {
                    double value = std::stod(valueStr);
                    aggregates[aggField].push_back(value);
                    currentVariantAggregates[aggField] = value;
                } catch (...) {
                    std::cerr << "Warning: Non-numeric value for INFO field \"" << aggField << "\": " << valueStr << "\n";
                    currentVariantAggregates[aggField] = 0.0;
                }
            } else {
                // INFO field not present; assign default value
                aggregates[aggField].push_back(0.0);
                currentVariantAggregates[aggField] = 0.0;
            }
        }

        // Append aggregated values to the VCF line
        out << line;
        for (const auto& aggField : infoFields) {
            out << "\t" << currentVariantAggregates[aggField];
        }
        out << "\n";
    }

    std::cout << "Aggregated INFO Fields:\n";
    for (const auto& agg : aggregates) {
        double sum = 0.0;
        for (const auto& val : agg.second) {
            sum += val;
        }
        double average = (agg.second.empty()) ? 0.0 : sum / agg.second.size();
        std::cout << agg.first << ": Sum = " << sum << ", Average = " << average << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXInfoAggregator infoAggregator;
    return infoAggregator.run(argc, argv);
}