#include "VCFX_outlier_detector.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

// Implementation

int VCFXOutlierDetector::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string metric = "AF"; // Default metric
    double threshold = 0.0;
    bool isVariant = true; // true: variant outliers, false: sample outliers

    static struct option long_options[] = {
        {"help",           no_argument,       0, 'h'},
        {"metric",         required_argument, 0, 'm'},
        {"threshold",      required_argument, 0, 't'},
        {"variant",        no_argument,       0, 'v'},
        {"sample",         no_argument,       0, 's'},
        {0,                0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hm:t:vs", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'm':
                metric = std::string(optarg);
                break;
            case 't':
                try {
                    threshold = std::stod(optarg);
                } catch (const std::invalid_argument&) {
                    std::cerr << "Error: Invalid threshold value.\n";
                    displayHelp();
                    return 1;
                }
                break;
            case 'v':
                isVariant = true;
                break;
            case 's':
                isVariant = false;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || threshold <= 0.0) {
        displayHelp();
        return 1;
    }

    if (isVariant) {
        // Detect variant outliers
        detectOutliers(std::cin, std::cout, metric, threshold, true);
    } else {
        // Detect sample outliers
        detectOutliers(std::cin, std::cout, metric, threshold, false);
    }

    return 0;
}

void VCFXOutlierDetector::displayHelp() {
    std::cout << "VCFX_outlier_detector: Detect outlier variants or samples based on specified quality metrics or allele frequencies.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_outlier_detector [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                Display this help message and exit\n";
    std::cout << "  -m, --metric <METRIC>     Specify the metric to use for outlier detection (e.g., AF, DP, QUAL)\n";
    std::cout << "  -t, --threshold <VALUE>   Specify the threshold for outlier detection\n";
    std::cout << "  -v, --variant             Detect outlier variants based on the specified metric\n";
    std::cout << "  -s, --sample              Detect outlier samples based on the specified metric\n\n";
    std::cout << "Examples:\n";
    std::cout << "  VCFX_outlier_detector --metric AF --threshold 0.05 --variant < input.vcf > variant_outliers.txt\n";
    std::cout << "  VCFX_outlier_detector --metric DP --threshold 200 --sample < input.vcf > sample_outliers.txt\n";
}

bool VCFXOutlierDetector::parseMetricFromInfo(const std::string& infoField, const std::string& metric, double& value) {
    // INFO fields are semicolon-separated key=value pairs
    std::stringstream ss(infoField);
    std::string token;
    while (std::getline(ss, token, ';')) {
        size_t eqPos = token.find('=');
        if (eqPos != std::string::npos) {
            std::string key = token.substr(0, eqPos);
            std::string valStr = token.substr(eqPos + 1);
            if (key == metric) {
                try {
                    value = std::stod(valStr);
                    return true;
                } catch (...) {
                    return false;
                }
            }
        }
    }
    return false;
}

bool VCFXOutlierDetector::parseMetricFromGenotype(const std::string& genotypeField, const std::string& metric, double& value) {
    // FORMAT fields are colon-separated key=value pairs
    std::stringstream ss(genotypeField);
    std::string token;
    while (std::getline(ss, token, ':')) {
        size_t eqPos = token.find('=');
        if (eqPos != std::string::npos) {
            std::string key = token.substr(0, eqPos);
            std::string valStr = token.substr(eqPos + 1);
            if (key == metric) {
                try {
                    value = std::stod(valStr);
                    return true;
                } catch (...) {
                    return false;
                }
            }
        }
    }
    return false;
}

void VCFXOutlierDetector::detectOutliers(std::istream& in, std::ostream& out, const std::string& metric, double threshold, bool isVariant) {
    std::string line;
    if (isVariant) {
        // Detect variant outliers
        std::cout << "Detecting variant outliers based on metric '" << metric << "' with threshold '" << threshold << "'.\n";
        out << "Chromosome\tPosition\tID\t" << metric << "\n";

        while (std::getline(in, line)) {
            if (line.empty()) continue;

            if (line[0] == '#') continue;

            std::stringstream ss(line);
            std::vector<std::string> fields;
            std::string field;

            // Split fields by tab
            while (std::getline(ss, field, '\t')) {
                fields.push_back(field);
            }

            if (fields.size() < 8) {
                std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields): " << line << "\n";
                continue;
            }

            std::string chrom = fields[0];
            std::string pos = fields[1];
            std::string id = fields[2];
            std::string info = fields[7];

            double metricValue = 0.0;
            if (parseMetricFromInfo(info, metric, metricValue)) {
                if (metricValue > threshold) {
                    out << chrom << "\t" << pos << "\t" << id << "\t" << metricValue << "\n";
                }
            }
        }
    }
    else {
        // Detect sample outliers
        std::cout << "Detecting sample outliers based on metric '" << metric << "' with threshold '" << threshold << "'.\n";

        std::vector<std::string> sampleNames;
        bool headerParsed = false;

        // Initialize sample-wise metrics
        std::unordered_map<std::string, double> sampleMetrics;
        std::unordered_map<std::string, int> sampleCounts;

        while (std::getline(in, line)) {
            if (line.empty()) continue;

            if (line.substr(0, 6) == "#CHROM") {
                std::stringstream ss(line);
                std::string field;
                // Parse header fields
                std::vector<std::string> headers;
                while (std::getline(ss, field, '\t')) {
                    headers.push_back(field);
                }
                // Sample names start from the 10th column
                for (size_t i = 9; i < headers.size(); ++i) {
                    sampleNames.push_back(headers[i]);
                    sampleMetrics[headers[i]] = 0.0;
                    sampleCounts[headers[i]] = 0;
                }
                headerParsed = true;
                continue;
            }

            if (!headerParsed) {
                std::cerr << "Error: VCF header line with #CHROM not found.\n";
                return;
            }

            std::stringstream ss(line);
            std::string chrom, pos, id, ref, alt, qual, filter, info, format;
            if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
                std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
                continue;
            }

            // Read genotype fields if present
            std::vector<std::string> genotypes;
            std::string genotype;
            while (ss >> genotype) {
                genotypes.push_back(genotype);
            }

            // Split FORMAT field
            std::vector<std::string> formatFields;
            std::stringstream formatSS(format);
            std::string fmt;
            while (std::getline(formatSS, fmt, ':')) {
                formatFields.push_back(fmt);
            }

            // Determine the index of the specified metric
            int metricIndex = -1;
            for (size_t i = 0; i < formatFields.size(); ++i) {
                if (formatFields[i] == metric) {
                    metricIndex = static_cast<int>(i);
                    break;
                }
            }

            if (metricIndex == -1) continue; // Metric not found

            for (size_t i = 0; i < genotypes.size() && i < sampleNames.size(); ++i) {
                std::stringstream gtSS(genotypes[i]);
                std::string gtField;
                std::vector<std::string> gtValues;
                while (std::getline(gtSS, gtField, ':')) {
                    gtValues.push_back(gtField);
                }

                if (metricIndex >= 0 && metricIndex < static_cast<int>(gtValues.size())) {
                    try {
                        double value = std::stod(gtValues[metricIndex]);
                        sampleMetrics[sampleNames[i]] += value;
                        sampleCounts[sampleNames[i]] += 1;
                    } catch (...) {
                        // Ignore invalid values
                        continue;
                    }
                }
            }
        }

        // Calculate average metric per sample and identify outliers
        out << "Sample\tAverage_" << metric << "\n";
        for (const auto& sample : sampleNames) {
            if (sampleCounts[sample] == 0) {
                out << sample << "\tNA\n";
                continue;
            }
            double avg = sampleMetrics[sample] / sampleCounts[sample];
            if (avg > threshold) {
                out << sample << "\t" << avg << "\n";
            }
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXOutlierDetector outlierDetector;
    return outlierDetector.run(argc, argv);
}