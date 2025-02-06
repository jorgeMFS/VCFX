./VCFX_nonref_filter/VCFX_nonref_filter.cpp
#include "VCFX_nonref_filter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>

// Implementation of VCFXNonRefFilter
int VCFXNonRefFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0,      0,           0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Filter VCF input from stdin and output to stdout
    filterNonRef(std::cin, std::cout);

    return 0;
}

void VCFXNonRefFilter::displayHelp() {
    std::cout << "VCFX_nonref_filter: Filter out variants where all genotypes are homozygous reference.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_nonref_filter [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_nonref_filter < input.vcf > filtered.vcf\n";
}

void VCFXNonRefFilter::filterNonRef(std::istream& in, std::ostream& out) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerLines;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            out << line << "\n"; // Preserve header lines
            if (line.substr(0, 6) == "#CHROM") {
                headerParsed = true;
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse VCF data lines
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fieldsVec;
        while (std::getline(ss, field, '\t')) {
            fieldsVec.push_back(field);
        }

        if (fieldsVec.size() < 10) {
            std::cerr << "Invalid VCF line with fewer than 10 fields.\n";
            continue;
        }

        std::string format = fieldsVec[8];

        // Split FORMAT field to find the index of genotype (GT)
        std::vector<std::string> formatFields;
        std::stringstream fmt_ss(format);
        std::string fmt_field;
        while (std::getline(fmt_ss, fmt_field, ':')) {
            formatFields.push_back(fmt_field);
        }

        int gt_index = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gt_index = static_cast<int>(i);
                break;
            }
        }

        if (gt_index == -1) {
            std::cerr << "GT field not found in FORMAT column.\n";
            continue;
        }

        bool allHomRef = true;

        // Iterate over each sample to check genotypes
        for (size_t i = 9; i < fieldsVec.size(); ++i) {
            std::string sample = fieldsVec[i];
            std::vector<std::string> sampleFields;
            std::stringstream samp_ss(sample);
            std::string samp_field;
            while (std::getline(samp_ss, samp_field, ':')) {
                sampleFields.push_back(samp_field);
            }

            if (gt_index >= static_cast<int>(sampleFields.size())) {
                std::cerr << "GT index out of range in sample fields.\n";
                allHomRef = false;
                break;
            }

            std::string genotype = sampleFields[gt_index];
            // Replace '|' with '/' for consistency
            std::replace(genotype.begin(), genotype.end(), '|', '/');

            // Check if genotype is homozygous reference
            if (genotype == "0/0" || genotype == "0") {
                continue;
            } else {
                allHomRef = false;
                break;
            }
        }

        if (!allHomRef) {
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXNonRefFilter filter;
    return filter.run(argc, argv);
}


./VCFX_nonref_filter/VCFX_nonref_filter.h
#ifndef VCFX_NONREF_FILTER_H
#define VCFX_NONREF_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXNonRefFilter: Header file for Non-Reference Variant Filter tool
class VCFXNonRefFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input to exclude variants with all homozygous reference genotypes
    void filterNonRef(std::istream& in, std::ostream& out);
};

#endif // VCFX_NONREF_FILTER_H


./VCFX_outlier_detector/VCFX_outlier_detector.cpp
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

./VCFX_outlier_detector/VCFX_outlier_detector.h
#ifndef VCFX_OUTLIER_DETECTOR_H
#define VCFX_OUTLIER_DETECTOR_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>

// VCFXOutlierDetector: Header file for Outlier Detection Tool
class VCFXOutlierDetector {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Detects outliers based on user-specified criteria
    void detectOutliers(std::istream& in, std::ostream& out, const std::string& metric, double threshold, bool isVariant);

    // Parses the INFO field to extract the specified metric
    bool parseMetricFromInfo(const std::string& infoField, const std::string& metric, double& value);

    // Parses genotype fields to extract specified metrics (if needed)
    bool parseMetricFromGenotype(const std::string& genotypeField, const std::string& metric, double& value);
};

#endif // VCFX_OUTLIER_DETECTOR_H


./VCFX_phase_checker/VCFX_phase_checker.cpp
#include "VCFX_phase_checker.h"
#include <getopt.h>

// Implementation of VCFX_phase_checker
int VCFXPhaseChecker::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    while ((opt = getopt(argc, argv, "h")) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Process VCF input from stdin
    processVCF(std::cin);

    return 0;
}

void VCFXPhaseChecker::displayHelp() {
    std::cout << "VCFX_phase_checker: Check if variants are phased in a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_phase_checker [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_phase_checker < input.vcf\n";
}

void VCFXPhaseChecker::processVCF(std::istream& in) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            std::cout << line << "\n"; // Preserve header lines
            continue;
        }

        // Split the VCF line into fields
        std::vector<std::string> fields;
        std::string field;
        size_t pos = 0;
        while ((pos = line.find('\t')) != std::string::npos) {
            field = line.substr(0, pos);
            fields.push_back(field);
            line.erase(0, pos + 1);
        }
        fields.push_back(line); // Add the last field

        if (fields.size() < 10) {
            std::cerr << "Invalid VCF line with fewer than 10 fields.\n";
            continue;
        }

        // GT is typically in FORMAT field followed by sample fields
        std::string format = fields[8];
        std::vector<std::string> format_fields;
        pos = 0;
        while ((pos = format.find(':')) != std::string::npos) {
            field = format.substr(0, pos);
            format_fields.push_back(field);
            format.erase(0, pos + 1);
        }
        format_fields.push_back(format);

        // Find the index of the GT field
        int gt_index = -1;
        for (size_t i = 0; i < format_fields.size(); ++i) {
            if (format_fields[i] == "GT") {
                gt_index = static_cast<int>(i);
                break;
            }
        }

        if (gt_index == -1) {
            std::cerr << "GT field not found in FORMAT column.\n";
            continue;
        }

        // Iterate over sample columns
        bool all_phased = true;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sample_fields;
            std::string sample = fields[i];
            size_t s_pos = 0;
            while ((s_pos = sample.find(':')) != std::string::npos) {
                field = sample.substr(0, s_pos);
                sample_fields.push_back(field);
                sample.erase(0, s_pos + 1);
            }
            sample_fields.push_back(sample);

            if (gt_index >= static_cast<int>(sample_fields.size())) {
                std::cerr << "GT index out of range in sample fields.\n";
                continue;
            }

            std::string genotype = sample_fields[gt_index];
            if (!isPhased(genotype)) {
                all_phased = false;
                break;
            }
        }

        if (all_phased) {
            std::cout << line << "\n";
        } else {
            std::cerr << "Unphased genotype found at position " << fields[1] << "\n";
        }
    }
}

bool VCFXPhaseChecker::isPhased(const std::string& genotype) {
    return genotype.find('|') != std::string::npos;
}

int main(int argc, char* argv[]) {
    VCFXPhaseChecker checker;
    return checker.run(argc, argv);
}



./VCFX_phase_checker/VCFX_phase_checker.h
#ifndef VCFX_PHASE_CHECKER_H
#define VCFX_PHASE_CHECKER_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_phase_checker: Header file for phased variant checking tool
class VCFXPhaseChecker {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input from the given stream
    void processVCF(std::istream& in);

    // Checks if a genotype is phased
    bool isPhased(const std::string& genotype);
};

#endif // VCFX_PHASE_CHECKER_H


./VCFX_phase_quality_filter/VCFX_phase_quality_filter.cpp
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

./VCFX_phase_quality_filter/VCFX_phase_quality_filter.h
#ifndef VCFX_PHASE_QUALITY_FILTER_H
#define VCFX_PHASE_QUALITY_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXPhaseQualityFilter: Header file for Variant Phasing Quality Filter Tool
class VCFXPhaseQualityFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on phasing quality scores
    void filterByPQ(std::istream& in, std::ostream& out, double threshold);

    // Parses the PQ score from the INFO field
    double parsePQScore(const std::string& infoField);
};

#endif // VCFX_PHASE_QUALITY_FILTER_H


./VCFX_phred_filter/VCFX_phred_filter.cpp
#include "VCFX_phred_filter.h"
#include <getopt.h>
#include <sstream>

// Implementation of VCFX_phred_filter
int VCFXPhredFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    double threshold = 30.0; // Default threshold
    bool showHelp = false;

    static struct option long_options[] = {
        {"phred-filter", required_argument, 0, 'p'},
        {"help",         no_argument,       0, 'h'},
        {0,              0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "p:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'p':
                {
                    std::stringstream ss(optarg);
                    if (!(ss >> threshold)) {
                        std::cerr << "Invalid threshold value: " << optarg << "\n";
                        return 1;
                    }
                }
                break;
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Process VCF input from stdin
    processVCF(std::cin, threshold);

    return 0;
}

void VCFXPhredFilter::displayHelp() {
    std::cout << "VCFX_phred_filter: Filter variants based on Phred quality score.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_phred_filter --phred-filter <threshold> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -p, --phred-filter   Phred quality score threshold (e.g., 30)\n";
    std::cout << "  -h, --help           Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_phred_filter --phred-filter 30 < input.vcf > filtered.vcf\n";
}

void VCFXPhredFilter::processVCF(std::istream& in, double threshold) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            std::cout << line << "\n"; // Preserve header lines
            continue;
        }

        // Split the VCF line into fields
        std::vector<std::string> fields;
        std::string field;
        size_t pos = 0;
        while ((pos = line.find('\t')) != std::string::npos) {
            field = line.substr(0, pos);
            fields.push_back(field);
            line.erase(0, pos + 1);
        }
        fields.push_back(line); // Add the last field

        if (fields.size() < 5) {
            std::cerr << "Invalid VCF line with fewer than 5 fields.\n";
            continue;
        }

        std::string qualStr = fields[5];
        double qual = parseQUAL(qualStr);

        if (qual >= threshold) {
            std::cout << line << "\n";
        }
    }
}

double VCFXPhredFilter::parseQUAL(const std::string& qualStr) {
    if (qualStr == ".") {
        return 0.0; // Treat missing QUAL as 0
    }
    try {
        return std::stod(qualStr);
    } catch (const std::invalid_argument&) {
        std::cerr << "Invalid QUAL value: " << qualStr << "\n";
        return 0.0;
    }
}

int main(int argc, char* argv[]) {
    VCFXPhredFilter filter;
    return filter.run(argc, argv);
}


./VCFX_phred_filter/VCFX_phred_filter.h
#ifndef VCFX_PHRED_FILTER_H
#define VCFX_PHRED_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_phred_filter: Header file for Phred score filtering tool
class VCFXPhredFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input from the given stream
    void processVCF(std::istream& in, double threshold);

    // Parses the QUAL field
    double parseQUAL(const std::string& qualStr);
};

#endif // VCFX_PHRED_FILTER_H


./VCFX_population_filter/VCFX_population_filter.cpp
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


./VCFX_population_filter/VCFX_population_filter.h
#ifndef VCFX_POPULATION_FILTER_H
#define VCFX_POPULATION_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXPopulationFilter: Header file for Population Subset Filter tool
class VCFXPopulationFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input to include only samples from the specified population
    void filterPopulation(std::istream& in, std::ostream& out, const std::string& populationTag, const std::string& popMapFile);
};

#endif // VCFX_POPULATION_FILTER_H


./VCFX_position_subsetter/VCFX_position_subsetter.cpp
#include "VCFX_position_subsetter.h"
#include <sstream>
#include <vector>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_position_subsetter\n"
              << "Usage: VCFX_position_subsetter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --region, -r \"CHR:START-END\"   Specify the genomic region to subset (e.g., \"chr1:10000-20000\").\n"
              << "  --help, -h                      Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Subsets VCF records based on the specified genomic region.\n\n"
              << "Examples:\n"
              << "  ./VCFX_position_subsetter --region \"chr1:10000-20000\" < input.vcf > subset.vcf\n";
}

// Helper structure to represent a genomic region
struct GenomicRegion {
    std::string chrom;
    int start;
    int end;
};

// Function to parse the region string
bool parseRegion(const std::string& region_str, GenomicRegion& region) {
    size_t colon = region_str.find(':');
    size_t dash = region_str.find('-');

    if (colon == std::string::npos || dash == std::string::npos || dash < colon) {
        std::cerr << "Error: Invalid region format. Expected format \"chrX:start-end\".\n";
        return false;
    }

    region.chrom = region_str.substr(0, colon);
    try {
        region.start = std::stoi(region_str.substr(colon + 1, dash - colon - 1));
        region.end = std::stoi(region_str.substr(dash + 1));
    } catch (...) {
        std::cerr << "Error: Unable to parse start or end positions.\n";
        return false;
    }

    if (region.start > region.end) {
        std::cerr << "Error: Start position is greater than end position.\n";
        return false;
    }

    return true;
}

// Function to subset VCF records based on genomic range
bool subsetVCFByPosition(std::istream& in, std::ostream& out, const std::string& region_str) {
    GenomicRegion region;
    if (!parseRegion(region_str, region)) {
        return false;
    }

    std::string line;
    bool header_found = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n"; // Preserve header lines
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        std::stringstream ss(line);
        std::string chrom, pos_str;
        // Extract CHROM and POS fields
        if (!(ss >> chrom >> pos_str)) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        int pos = 0;
        try {
            pos = std::stoi(pos_str);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value. Skipping line.\n";
            continue;
        }

        if (chrom == region.chrom && pos >= region.start && pos <= region.end) {
            out << line << "\n";
        }
    }

    return true;
}

int main(int argc, char* argv[]) {
    std::string region_str;

    // Argument parsing for help and region
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--region" || arg == "-r") && i + 1 < argc) {
            region_str = argv[++i];
        } else if (arg.find("--region=") == 0) {
            region_str = arg.substr(9);
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    if (region_str.empty()) {
        std::cerr << "Error: Genomic region not specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    bool success = subsetVCFByPosition(std::cin, std::cout, region_str);
    return success ? 0 : 1;
}


./VCFX_position_subsetter/VCFX_position_subsetter.h
#ifndef VCFX_POSITION_SUBSETTER_H
#define VCFX_POSITION_SUBSETTER_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Function to subset VCF records based on genomic range
bool subsetVCFByPosition(std::istream& in, std::ostream& out, const std::string& region);

#endif // VCFX_POSITION_SUBSETTER_H


./VCFX_probability_filter/VCFX_probability_filter.cpp
#include "VCFX_probability_filter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <regex>

// Implementation of VCFXProbabilityFilter
int VCFXProbabilityFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string condition;

    static struct option long_options[] = {
        {"help",           no_argument,       0, 'h'},
        {"filter-probability", required_argument, 0, 'f'},
        {0,                0,                 0,  0 }
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

    // Perform probability filtering on stdin and output to stdout
    filterByProbability(std::cin, std::cout, condition);

    return 0;
}

void VCFXProbabilityFilter::displayHelp() {
    std::cout << "VCFX_probability_filter: Filter VCF based on genotype probability scores.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_probability_filter --filter-probability \"<CONDITION>\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                        Display this help message and exit\n";
    std::cout << "  -f, --filter-probability <cond>    Specify the genotype probability filter condition (e.g., GP>0.9)\n\n";
    std::cout << "Supported Operators: >, <, >=, <=, ==, !=\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_probability_filter --filter-probability \"GP>0.9\" < input.vcf > filtered.vcf\n";
}

void VCFXProbabilityFilter::filterByProbability(std::istream& in, std::ostream& out, const std::string& condition) {
    // Parse the filter condition using regex (e.g., "GP>0.9")
    std::regex conditionRegex(R"((\w+)\s*(>=|<=|>|<|==|!=)\s*([0-9]*\.?[0-9]+))");
    std::smatch matches;
    if (!std::regex_match(condition, matches, conditionRegex)) {
        std::cerr << "Error: Invalid filter condition format. Expected format like \"GP>0.9\".\n";
        return;
    }

    std::string field = matches[1];
    std::string op = matches[2];
    double threshold = std::stod(matches[3]);

    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    size_t formatIndex = std::string::npos;
    size_t fieldIndex = std::string::npos;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Handle header lines
            if (line.substr(0, 6) == "#CHROM") {
                // Parse header to identify the FORMAT column and sample columns
                out << line << "\n";
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

        if (fieldsVec.size() < 9) {
            std::cerr << "Warning: Invalid VCF line with fewer than 9 fields: " << line << "\n";
            continue;
        }

        std::string formatField = fieldsVec[8];
        std::vector<std::string> formatFields;
        std::stringstream fmt_ss(formatField);
        std::string fmt_field;
        while (std::getline(fmt_ss, fmt_field, ':')) {
            formatFields.push_back(fmt_field);
        }

        // Find the index of the specified field in the FORMAT column
        if (fieldIndex == std::string::npos) {
            for (size_t i = 0; i < formatFields.size(); ++i) {
                if (formatFields[i] == field) {
                    fieldIndex = i;
                    break;
                }
            }

            if (fieldIndex == std::string::npos) {
                std::cerr << "Error: Specified field \"" << field << "\" not found in FORMAT column.\n";
                return;
            }
        }

        bool pass = true;

        // Iterate over each sample to check the genotype probability
        for (size_t i = 9; i < fieldsVec.size(); ++i) {
            std::string sample = fieldsVec[i];
            std::stringstream samp_ss(sample);
            std::string samp_field;
            std::vector<std::string> sampleFields;

            while (std::getline(samp_ss, samp_field, ':')) {
                sampleFields.push_back(samp_field);
            }

            if (fieldIndex >= sampleFields.size()) {
                std::cerr << "Warning: Field index out of range in sample fields.\n";
                pass = false;
                break;
            }

            std::string valueStr = sampleFields[fieldIndex];
            if (valueStr.empty() || valueStr == ".") {
                pass = false;
                break;
            }

            double value;
            try {
                value = std::stod(valueStr);
            } catch (...) {
                std::cerr << "Warning: Unable to convert value \"" << valueStr << "\" to number.\n";
                pass = false;
                break;
            }

            // Apply the filter condition
            if (op == ">") {
                if (!(value > threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "<") {
                if (!(value < threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == ">=") {
                if (!(value >= threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "<=") {
                if (!(value <= threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "==") {
                if (!(value == threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "!=") {
                if (!(value != threshold)) {
                    pass = false;
                    break;
                }
            }
        }

        if (pass) {
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXProbabilityFilter probabilityFilter;
    return probabilityFilter.run(argc, argv);
}

./VCFX_probability_filter/VCFX_probability_filter.h
#ifndef VCFX_PROBABILITY_FILTER_H
#define VCFX_PROBABILITY_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXProbabilityFilter: Header file for Genotype Probability Filter tool
class VCFXProbabilityFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on the specified genotype probability condition
    void filterByProbability(std::istream& in, std::ostream& out, const std::string& condition);
};

#endif // VCFX_PROBABILITY_FILTER_H


./VCFX_quality_adjuster/VCFX_quality_adjuster.cpp
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

int main(int argc, char* argv[]) {
    VCFXQualityAdjuster qualityAdjuster;
    return qualityAdjuster.run(argc, argv);
}   

./VCFX_quality_adjuster/VCFX_quality_adjuster.h
#ifndef VCFX_QUALITY_ADJUSTER_H
#define VCFX_QUALITY_ADJUSTER_H

#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <unordered_map>

// VCFXQualityAdjuster: Header file for Quality Score Adjuster Tool
class VCFXQualityAdjuster {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Parses the transformation function string and returns a corresponding function
    bool parseTransformationFunction(const std::string& funcStr, std::function<double(double)> &transFunc);

    // Adjusts the QUAL scores based on the transformation function
    void adjustQualityScores(std::istream& in, std::ostream& out, std::function<double(double)> transFunc);

    // Supported transformation functions
    std::unordered_map<std::string, std::function<double(double)>> supportedFunctions = {
        {"log", [](double x) -> double { return std::log(x + 1e-10); }}, // Added epsilon to avoid log(0)
        {"sqrt", [](double x) -> double { return std::sqrt(x); }},
        {"square", [](double x) -> double { return x * x; }},
        {"identity", [](double x) -> double { return x; }}
    };
};

#endif // VCFX_QUALITY_ADJUSTER_H


./VCFX_record_filter/VCFX_record_filter.cpp
#include "VCFX_record_filter.h"
#include <sstream>
#include <vector>
#include <algorithm>

// Helper function to trim whitespace from both ends of a string
static inline std::string trim(const std::string& s) {
    auto start = s.begin();
    while (start != s.end() && isspace(*start)) {
        start++;
    }

    auto end = s.end();
    do {
        end--;
    } while (distance(start, end) > 0 && isspace(*end));

    return std::string(start, end + 1);
}

// Function to parse individual filter criterion
bool parseSingleCriterion(const std::string& token, FilterCriterion& criterion) {
    size_t pos = token.find_first_of("><=!");
    if (pos == std::string::npos) {
        return false;
    }

    // Determine the operator
    Operator op;
    std::string op_str;
    if (token[pos] == '>') {
        if (pos + 1 < token.size() && token[pos + 1] == '=') {
            op = Operator::GREATER_EQUAL;
            op_str = ">=";
            pos += 1;
        } else {
            op = Operator::GREATER_THAN;
            op_str = ">";
        }
    } else if (token[pos] == '<') {
        if (pos + 1 < token.size() && token[pos + 1] == '=') {
            op = Operator::LESS_EQUAL;
            op_str = "<=";
            pos += 1;
        } else {
            op = Operator::LESS_THAN;
            op_str = "<";
        }
    } else if (token[pos] == '=') {
        if (pos + 1 < token.size() && token[pos + 1] == '=') {
            op = Operator::EQUAL;
            op_str = "==";
            pos += 1;
        } else {
            // Single '=' treated as '=='
            op = Operator::EQUAL;
            op_str = "==";
        }
    } else {
        // Unsupported operator
        return false;
    }

    // Extract field and value
    std::string field = trim(token.substr(0, pos));
    std::string value_str = trim(token.substr(pos + 1));

    if (field.empty() || value_str.empty()) {
        return false;
    }

    try {
        double value = std::stod(value_str);
        criterion.field = field;
        criterion.op = op;
        criterion.value = value;
    } catch (const std::invalid_argument&) {
        return false;
    }

    return true;
}

// Function to parse filter criteria string into structured criteria
bool parseCriteria(const std::string& criteria_str, std::vector<FilterCriterion>& criteria) {
    std::stringstream ss(criteria_str);
    std::string token;

    while (std::getline(ss, token, ';')) {
        token = trim(token);
        if (token.empty()) {
            continue;
        }

        FilterCriterion criterion;
        if (!parseSingleCriterion(token, criterion)) {
            std::cerr << "Failed to parse filter criterion: '" << token << "'" << std::endl;
            return false;
        }

        criteria.push_back(criterion);
    }

    if (criteria.empty()) {
        std::cerr << "No valid filter criteria found." << std::endl;
        return false;
    }

    return true;
}

// Function to split a string by a delimiter
static inline std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::stringstream ss(s);
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to retrieve the value of a specific field from a VCF record
bool getFieldValue(const std::vector<std::string>& fields, const std::string& field_name, double& value) {
    static const std::vector<std::string> standard_fields = {
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
    };

    auto it = std::find(standard_fields.begin(), standard_fields.end(), field_name);
    if (it != standard_fields.end()) {
        size_t index = std::distance(standard_fields.begin(), it);
        if (index >= fields.size()) {
            return false;
        }

        if (field_name == "POS") {
            try {
                value = std::stod(fields[index]);
            } catch (...) {
                return false;
            }
            return true;
        } else if (field_name == "QUAL") {
            if (fields[index] == ".") {
                value = 0.0;
            } else {
                try {
                    value = std::stod(fields[index]);
                } catch (...) {
                    return false;
                }
            }
            return true;
        } else {
            // Currently, only POS and QUAL are supported as numeric fields
            return false;
        }
    } else {
        // Attempt to extract from INFO field
        if (standard_fields.size() < 8) {
            return false;
        }
        std::string info_field = fields[7];
        std::vector<std::string> info_entries = split(info_field, ';');
        for (const auto& entry : info_entries) {
            size_t eq = entry.find('=');
            if (eq != std::string::npos) {
                std::string key = trim(entry.substr(0, eq));
                std::string val_str = trim(entry.substr(eq + 1));
                if (key == field_name) {
                    try {
                        value = std::stod(val_str);
                        return true;
                    } catch (...) {
                        return false;
                    }
                }
            } else {
                if (entry == field_name) {
                    // Flag present, treat as true (1.0)
                    value = 1.0;
                    return true;
                }
            }
        }
    }

    return false;
}

// Function to apply all filters to a single record
bool applyFilters(const std::string& record, const std::vector<FilterCriterion>& criteria) {
    std::vector<std::string> fields = split(record, '\t');

    for (const auto& criterion : criteria) {
        double field_value = 0.0;
        if (!getFieldValue(fields, criterion.field, field_value)) {
            // If the field is not found or not numeric, skip this criterion
            return false;
        }

        bool condition_met = false;
        switch (criterion.op) {
            case Operator::GREATER_THAN:
                condition_met = (field_value > criterion.value);
                break;
            case Operator::LESS_THAN:
                condition_met = (field_value < criterion.value);
                break;
            case Operator::GREATER_EQUAL:
                condition_met = (field_value >= criterion.value);
                break;
            case Operator::LESS_EQUAL:
                condition_met = (field_value <= criterion.value);
                break;
            case Operator::EQUAL:
                condition_met = (field_value == criterion.value);
                break;
        }

        if (!condition_met) {
            return false; // If any criterion is not met, reject the record
        }
    }

    return true; // All criteria met
}

// Function to process and filter records
void processRecords(std::istream& in, std::ostream& out, const std::vector<FilterCriterion>& criteria) {
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            continue; // Skip header lines
        }

        if (applyFilters(line, criteria)) {
            out << line << "\n";
        }
    }
}

void printHelp() {
    std::cout << "VCFX_record_filter\n"
              << "Usage: VCFX_record_filter --filter \"CRITERIA\" [OPTIONS]\n\n"
              << "Options:\n"
              << "  --filter, -f          Specify filter criteria (e.g., \"QUAL>30;DP<100\").\n"
              << "  --help, -h            Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Filters VCF records based on specified criteria.\n\n"
              << "Example:\n"
              << "  ./VCFX_record_filter --filter \"QUAL>30;DP<100\" < input.vcf > filtered.vcf\n";
}

int main(int argc, char* argv[]) {
    // Argument parsing
    std::string criteria_str;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
        if (arg == "--filter" || arg == "-f") {
            if (i + 1 < argc) {
                criteria_str = argv[++i];
            } else {
                std::cerr << "Error: --filter option requires an argument.\n";
                return 1;
            }
        } else if (arg.find("--filter=") == 0) {
            criteria_str = arg.substr(9);
        }
    }

    if (criteria_str.empty()) {
        std::cerr << "No filter criteria provided.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    std::vector<FilterCriterion> criteria;
    if (!parseCriteria(criteria_str, criteria)) {
        std::cerr << "Failed to parse filter criteria.\n";
        return 1;
    }

    processRecords(std::cin, std::cout, criteria);
    return 0;
}


./VCFX_record_filter/VCFX_record_filter.h
#ifndef VCFX_RECORD_FILTER_H
#define VCFX_RECORD_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// Enumeration for comparison operators
enum class Operator {
    GREATER_THAN,    // >
    LESS_THAN,       // <
    GREATER_EQUAL,   // >=
    LESS_EQUAL,      // <=
    EQUAL            // ==
};

// Structure to hold individual filter criteria
struct FilterCriterion {
    std::string field;
    Operator op;
    double value;
};

// Function to parse and populate filter criteria from a string
bool parseCriteria(const std::string& criteria_str, std::vector<FilterCriterion>& criteria);

// Function to apply all filter criteria to a single VCF record
bool applyFilters(const std::string& record, const std::vector<FilterCriterion>& criteria);

// Function to process and filter VCF records based on the provided criteria
void processRecords(std::istream& in, std::ostream& out, const std::vector<FilterCriterion>& criteria);

#endif // VCFX_RECORD_FILTER_H


./VCFX_ref_comparator/VCFX_ref_comparator.cpp
#include "VCFX_ref_comparator.h"
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <algorithm>

// Implementation of VCFXRefComparator
int VCFXRefComparator::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string referencePath;

    static struct option long_options[] = {
        {"help",        no_argument,       0, 'h'},
        {"reference",   required_argument, 0, 'r'},
        {0,             0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hr:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'r':
                referencePath = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || referencePath.empty()) {
        displayHelp();
        return 1;
    }

    // Load the reference genome
    if (!loadReference(referencePath)) {
        std::cerr << "Error: Failed to load reference genome from " << referencePath << "\n";
        return 1;
    }

    // Compare VCF variants against the reference genome
    compareWithReference(std::cin, std::cout);

    return 0;
}

void VCFXRefComparator::displayHelp() {
    std::cout << "VCFX_ref_comparator: Compare VCF variants against a reference genome.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_ref_comparator --reference <reference.fasta> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                  Display this help message and exit\n";
    std::cout << "  -r, --reference <file>      Path to the reference genome FASTA file\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_ref_comparator --reference reference.fasta < input.vcf > comparison_output.vcf\n";
}

bool VCFXRefComparator::loadReference(const std::string& referencePath) {
    std::ifstream refFile(referencePath);
    if (!refFile.is_open()) {
        std::cerr << "Error: Cannot open reference genome file: " << referencePath << "\n";
        return false;
    }

    std::string line;
    std::string currentChrom;
    std::ostringstream sequenceStream;

    while (std::getline(refFile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Save the previous chromosome sequence
            if (!currentChrom.empty()) {
                referenceGenome[currentChrom] = sequenceStream.str();
                sequenceStream.str("");
                sequenceStream.clear();
            }
            // Extract chromosome name
            std::stringstream ss(line.substr(1));
            ss >> currentChrom;
        } else {
            // Append sequence lines, removing any whitespace and converting to uppercase
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            sequenceStream << line;
        }
    }

    // Save the last chromosome sequence
    if (!currentChrom.empty()) {
        referenceGenome[currentChrom] = sequenceStream.str();
    }

    refFile.close();
    return true;
}

void VCFXRefComparator::compareWithReference(std::istream& vcfInput, std::ostream& vcfOutput) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    int infoIndex = -1;

    while (std::getline(vcfInput, line)) {
        if (line.empty()) {
            vcfOutput << "\n";
            continue;
        }

        if (line[0] == '#') {
            if (line.substr(0, 6) == "#CHROM") {
                // Parse header line
                std::stringstream ss(line);
                std::string field;
                headerFields.clear();
                while (std::getline(ss, field, '\t')) {
                    headerFields.push_back(field);
                }

                // Add new INFO field for comparison result
                std::cout << "##INFO=<ID=REF_COMPARISON,Number=1,Type=String,Description=\"Comparison of variant alleles against the reference genome\">\n";
                
                // Write modified header
                std::stringstream newHeader;
                for (size_t i = 0; i < headerFields.size(); ++i) {
                    newHeader << headerFields[i];
                    if (i != headerFields.size() - 1) {
                        newHeader << "\t";
                    }
                }
                newHeader << "\tREF_COMPARISON\n";
                vcfOutput << newHeader.str() << "\n";
                headerParsed = true;
            } else {
                // Write other header lines as-is
                vcfOutput << line << "\n";
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse VCF data line
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

        std::string chrom = fieldsVec[0];
        int pos = std::stoi(fieldsVec[1]);
        std::string id = fieldsVec[2];
        std::string ref = fieldsVec[3];
        std::string alt = fieldsVec[4];
        std::string qual = fieldsVec[5];
        std::string filter = fieldsVec[6];
        std::string info = fieldsVec[7];

        // Retrieve reference allele from the reference genome
        if (referenceGenome.find(chrom) == referenceGenome.end()) {
            std::cerr << "Warning: Chromosome " << chrom << " not found in reference genome.\n";
            info += ";REF_COMPARISON=UNKNOWN_CHROM";
            vcfOutput << line << "\t" << "UNKNOWN_CHROM" << "\n";
            continue;
        }

        const std::string& refSeq = referenceGenome[chrom];
        if (pos < 1 || pos > static_cast<int>(refSeq.size())) {
            std::cerr << "Warning: Position " << pos << " out of bounds for chromosome " << chrom << ".\n";
            info += ";REF_COMPARISON=INVALID_POS";
            vcfOutput << line << "\t" << "INVALID_POS" << "\n";
            continue;
        }

        // Extract reference allele from the reference genome
        std::string refFromGenome = refSeq.substr(pos - 1, ref.size());

        // Compare VCF ref allele with reference genome
        bool refMatch = (ref == refFromGenome);

        // Compare alternate alleles with reference genome
        std::vector<std::string> altAlleles;
        std::stringstream alt_ss(alt);
        std::string alt_allele;
        while (std::getline(alt_ss, alt_allele, ',')) {
            altAlleles.push_back(alt_allele);
        }

        std::vector<std::string> comparisonResults;
        for (const auto& altA : altAlleles) {
            if (altA == refFromGenome) {
                comparisonResults.push_back("REF_MATCH");
            } else {
                comparisonResults.push_back("NOVEL");
            }
        }

        // Prepare comparison result string
        std::string comparisonStr;
        for (size_t i = 0; i < comparisonResults.size(); ++i) {
            comparisonStr += comparisonResults[i];
            if (i != comparisonResults.size() - 1) {
                comparisonStr += ",";
            }
        }

        // Append to INFO field
        if (info.back() != ';') {
            info += ";";
        }
        info += "REF_COMPARISON=" + comparisonStr;

        // Reconstruct VCF line with updated INFO field
        std::stringstream newLine;
        for (size_t i = 0; i < 8; ++i) {
            newLine << fieldsVec[i] << "\t";
        }
        newLine << info << "\t";

        // Append any additional fields (e.g., FORMAT and samples)
        for (size_t i = 8; i < fieldsVec.size(); ++i) {
            newLine << fieldsVec[i];
            if (i != fieldsVec.size() - 1) {
                newLine << "\t";
            }
        }

        // Append the comparison result
        newLine << "\t" << comparisonStr;

        vcfOutput << newLine.str() << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXRefComparator refComparator;
    return refComparator.run(argc, argv);
}

./VCFX_ref_comparator/VCFX_ref_comparator.h
#ifndef VCFX_REF_COMPARATOR_H
#define VCFX_REF_COMPARATOR_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// VCFXRefComparator: Header file for Reference Genome Comparator tool
class VCFXRefComparator {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads the reference genome from a FASTA file
    bool loadReference(const std::string& referencePath);

    // Compares VCF variants against the reference genome
    void compareWithReference(std::istream& vcfInput, std::ostream& vcfOutput);

    // Reference genome data: chromosome -> sequence
    std::unordered_map<std::string, std::string> referenceGenome;
};

#endif // VCFX_REF_COMPARATOR_H


./VCFX_reformatter/VCFX_reformatter.cpp
#include "VCFX_reformatter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cctype>

// Implementation

int VCFXReformatter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::vector<std::string> compressFields;
    std::vector<std::string> reorderInfoFields;
    std::vector<std::string> reorderFormatFields;

    static struct option long_options[] = {
        {"help",            no_argument,       0, 'h'},
        {"compress-info",   required_argument, 0, 'c'},
        {"compress-format", required_argument, 0, 'f'},
        {"reorder-info",    required_argument, 0, 'i'},
        {"reorder-format",  required_argument, 0, 'o'},
        {0,                 0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hc:f:i:o:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'c':
                {
                    // Compress INFO fields separated by commas
                    std::string fields(optarg);
                    std::stringstream ss(fields);
                    std::string field;
                    while (std::getline(ss, field, ',')) {
                        // Trim whitespace
                        field.erase(field.find_last_not_of(" \n\r\t")+1);
                        field.erase(0, field.find_first_not_of(" \n\r\t"));
                        if (!field.empty()) {
                            compressFields.push_back(field);
                        }
                    }
                }
                break;
            case 'f':
                {
                    // Compress FORMAT fields separated by commas
                    std::string fields(optarg);
                    std::stringstream ss(fields);
                    std::string field;
                    while (std::getline(ss, field, ',')) {
                        // Trim whitespace
                        field.erase(field.find_last_not_of(" \n\r\t")+1);
                        field.erase(0, field.find_first_not_of(" \n\r\t"));
                        if (!field.empty()) {
                            compressFields.push_back(field);
                        }
                    }
                }
                break;
            case 'i':
                {
                    // Reorder INFO fields, order separated by commas
                    std::string fields(optarg);
                    std::stringstream ss(fields);
                    std::string field;
                    while (std::getline(ss, field, ',')) {
                        // Trim whitespace
                        field.erase(field.find_last_not_of(" \n\r\t")+1);
                        field.erase(0, field.find_first_not_of(" \n\r\t"));
                        if (!field.empty()) {
                            reorderInfoFields.push_back(field);
                        }
                    }
                }
                break;
            case 'o':
                {
                    // Reorder FORMAT fields, order separated by commas
                    std::string fields(optarg);
                    std::stringstream ss(fields);
                    std::string field;
                    while (std::getline(ss, field, ',')) {
                        // Trim whitespace
                        field.erase(field.find_last_not_of(" \n\r\t")+1);
                        field.erase(0, field.find_first_not_of(" \n\r\t"));
                        if (!field.empty()) {
                            reorderFormatFields.push_back(field);
                        }
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

    // Perform VCF reformatting
    reformatVCF(std::cin, std::cout, compressFields, reorderInfoFields, reorderFormatFields);

    return 0;
}

void VCFXReformatter::displayHelp() {
    std::cout << "VCFX_reformatter: Reformat VCF fields (e.g., compressing or reordering INFO/FORMAT fields).\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_reformatter [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                      Display this help message and exit\n";
    std::cout << "  -c, --compress-info <FIELDS>    Compress specified INFO fields (comma-separated)\n";
    std::cout << "  -f, --compress-format <FIELDS>  Compress specified FORMAT fields (comma-separated)\n";
    std::cout << "  -i, --reorder-info <ORDER>      Reorder INFO fields as per specified order (comma-separated)\n";
    std::cout << "  -o, --reorder-format <ORDER>    Reorder FORMAT fields as per specified order (comma-separated)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_reformatter --compress-info AF,DP --reorder-info AF,DP,INFO < input.vcf > reformatted.vcf\n";
}

void VCFXReformatter::reformatVCF(std::istream& in, std::ostream& out,
                                 const std::vector<std::string>& compressFields,
                                 const std::vector<std::string>& reorderInfoFields,
                                 const std::vector<std::string>& reorderFormatFields) {
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

        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string qual = fields[5];
        std::string filter = fields[6];
        std::string info = fields[7];

        // If there are genotype fields, they start from index 8
        std::vector<std::string> genotypeFields;
        std::string formatField;
        if (fields.size() > 8) {
            formatField = fields[8];
            genotypeFields.assign(fields.begin() + 9, fields.end());
        }

        // Compress specified INFO fields
        if (!compressFields.empty()) {
            info = compressFieldsFunction(info, compressFields);
        }

        // Reorder INFO fields
        if (!reorderInfoFields.empty()) {
            info = reorderInfo(info, reorderInfoFields);
        }

        // Compress specified FORMAT fields
        if (!compressFields.empty() && !formatField.empty()) {
            formatField = compressFieldsFunction(formatField, compressFields);
        }

        // Reorder FORMAT fields
        if (!reorderFormatFields.empty() && !formatField.empty()) {
            formatField = reorderFormat(formatField, reorderFormatFields);
        }

        // Reconstruct the VCF line
        std::ostringstream oss;
        oss << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
            << "\t" << qual << "\t" << filter << "\t" << info;

        if (!formatField.empty()) {
            oss << "\t" << formatField;
            for (const auto& gt : genotypeFields) {
                oss << "\t" << gt;
            }
        }

        oss << "\n";
        out << oss.str();
    }
}

std::string VCFXReformatter::compressFieldsFunction(const std::string& fieldValue, const std::vector<std::string>& fieldsToCompress) {
    // Example compression: Remove specified key-value pairs from INFO or FORMAT fields
    std::stringstream ss(fieldValue);
    std::string token;
    std::vector<std::string> compressedTokens;

    while (std::getline(ss, token, ';')) {
        bool toCompress = false;
        for (const auto& key : fieldsToCompress) {
            if (token.find(key + "=") == 0) { // Starts with key=
                toCompress = true;
                break;
            }
        }
        if (!toCompress) {
            compressedTokens.push_back(token);
        }
    }

    // Reconstruct the field
    std::string compressedField;
    for (size_t i = 0; i < compressedTokens.size(); ++i) {
        compressedField += compressedTokens[i];
        if (i != compressedTokens.size() - 1) {
            compressedField += ";";
        }
    }

    return compressedField;
}

std::string VCFXReformatter::reorderInfo(const std::string& infoField, const std::vector<std::string>& reorderOrder) {
    // Split INFO field into key-value pairs
    std::stringstream ss(infoField);
    std::string token;
    std::unordered_map<std::string, std::string> infoMap;
    std::vector<std::string> keys;

    while (std::getline(ss, token, ';')) {
        size_t eqPos = token.find('=');
        if (eqPos != std::string::npos) {
            std::string key = token.substr(0, eqPos);
            std::string value = token.substr(eqPos + 1);
            infoMap[key] = value;
            keys.push_back(key);
        } else {
            // Flag without value
            infoMap[token] = "";
            keys.push_back(token);
        }
    }

    // Reorder based on reorderOrder
    std::ostringstream oss;
    for (size_t i = 0; i < reorderOrder.size(); ++i) {
        const std::string& key = reorderOrder[i];
        if (infoMap.find(key) != infoMap.end()) {
            if (!infoMap[key].empty()) {
                oss << key << "=" << infoMap[key];
            } else {
                oss << key;
            }
            if (i != reorderOrder.size() - 1) {
                oss << ";";
            }
            infoMap.erase(key);
        }
    }

    // Append remaining keys in original order
    for (const auto& key : keys) {
        if (infoMap.find(key) != infoMap.end()) {
            if (!oss.str().empty() && oss.str().back() != ';') {
                oss << ";";
            }
            if (!infoMap[key].empty()) {
                oss << key << "=" << infoMap[key];
            } else {
                oss << key;
            }
            infoMap.erase(key);
        }
    }

    return oss.str();
}

std::string VCFXReformatter::reorderFormat(const std::string& formatField, const std::vector<std::string>& reorderOrder) {
    // Split FORMAT field into keys
    std::stringstream ss(formatField);
    std::string token;
    std::vector<std::string> formatKeys;
    std::vector<std::string> remainingKeys;

    while (std::getline(ss, token, ':')) {
        formatKeys.push_back(token);
    }

    // Reorder based on reorderOrder
    std::vector<std::string> reorderedFormat;
    for (const auto& key : reorderOrder) {
        auto it = std::find(formatKeys.begin(), formatKeys.end(), key);
        if (it != formatKeys.end()) {
            reorderedFormat.push_back(key);
            formatKeys.erase(it);
        }
    }

    // Append remaining keys in original order
    for (const auto& key : formatKeys) {
        reorderedFormat.push_back(key);
    }

    // Reconstruct the FORMAT field
    std::string reorderedFormatStr;
    for (size_t i = 0; i < reorderedFormat.size(); ++i) {
        reorderedFormatStr += reorderedFormat[i];
        if (i != reorderedFormat.size() - 1) {
            reorderedFormatStr += ":";
        }
    }

    return reorderedFormatStr;
}

int main(int argc, char* argv[]) {
    VCFXReformatter reformatter;
    return reformatter.run(argc, argv);
}

./VCFX_reformatter/VCFX_reformatter.h
#ifndef VCFX_REFORMATTER_H
#define VCFX_REFORMATTER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXReformatter: Header file for VCF Reformatting Tool
class VCFXReformatter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Reformat VCF input based on user-specified rules
    void reformatVCF(std::istream& in, std::ostream& out,
                    const std::vector<std::string>& compressFields,
                    const std::vector<std::string>& reorderInfoFields,
                    const std::vector<std::string>& reorderFormatFields);

    // Compresses specified INFO or FORMAT fields
    std::string compressFieldsFunction(const std::string& fieldValue, const std::vector<std::string>& fieldsToCompress);

    // Reorders INFO fields based on user-specified order
    std::string reorderInfo(const std::string& infoField, const std::vector<std::string>& reorderOrder);

    // Reorders FORMAT fields based on user-specified order
    std::string reorderFormat(const std::string& formatField, const std::vector<std::string>& reorderOrder);
};

#endif // VCFX_REFORMATTER_H


./VCFX_region_subsampler/VCFX_region_subsampler.cpp
#include "VCFX_region_subsampler.h"
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <algorithm>

// Implementation of VCFXRegionSubsampler
int VCFXRegionSubsampler::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string bedFilePath;

    static struct option long_options[] = {
        {"help",         no_argument,       0, 'h'},
        {"region-bed",   required_argument, 0, 'b'},
        {0,              0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hb:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'b':
                bedFilePath = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || bedFilePath.empty()) {
        displayHelp();
        return 1;
    }

    // Load regions from BED file
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> regions;
    if (!loadRegions(bedFilePath, regions)) {
        std::cerr << "Error: Failed to load regions from " << bedFilePath << "\n";
        return 1;
    }

    // Subsample regions from VCF
    subsampleRegions(std::cin, std::cout, regions);

    return 0;
}

void VCFXRegionSubsampler::displayHelp() {
    std::cout << "VCFX_region_subsampler: Subsample variants from specific genomic regions defined in a BED file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_region_subsampler --region-bed <regions.bed> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                Display this help message and exit\n";
    std::cout << "  -b, --region-bed <bed>    Specify the BED file with genomic regions\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_region_subsampler --region-bed regions.bed < input.vcf > subsampled.vcf\n";
}

bool VCFXRegionSubsampler::loadRegions(const std::string& bedFilePath, std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions) {
    std::ifstream infile(bedFilePath);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open BED file " << bedFilePath << "\n";
        return false;
    }

    std::string line;
    size_t line_num = 0;
    while (std::getline(infile, line)) {
        line_num++;
        if (line.empty() || line[0] == '#') {
            continue; // Skip comments and empty lines
        }

        std::stringstream ss(line);
        std::string chrom;
        int start, end;

        if (!(ss >> chrom >> start >> end)) {
            std::cerr << "Warning: Skipping invalid BED line " << line_num << ": " << line << "\n";
            continue;
        }

        // BED format uses 0-based start and 1-based end
        regions[chrom].emplace_back(std::make_pair(start + 1, end));
    }

    infile.close();
    return true;
}

bool VCFXRegionSubsampler::isVariantInRegions(const std::string& chrom, int pos, const std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions) {
    auto it = regions.find(chrom);
    if (it == regions.end()) {
        return false;
    }

    const auto& regionList = it->second;
    for (const auto& region : regionList) {
        if (pos >= region.first && pos <= region.second) {
            return true;
        }
    }

    return false;
}

void VCFXRegionSubsampler::subsampleRegions(std::istream& in, std::ostream& out, const std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions) {
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, posStr, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> posStr >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        int pos;
        try {
            pos = std::stoi(posStr);
        } catch (const std::invalid_argument&) {
            std::cerr << "Warning: Invalid position \"" << posStr << "\" in line: " << line << "\n";
            continue;
        }

        if (isVariantInRegions(chrom, pos, regions)) {
            // Output the variant line
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXRegionSubsampler regionSubsampler;
    return regionSubsampler.run(argc, argv);
}

./VCFX_region_subsampler/VCFX_region_subsampler.h
#ifndef VCFX_REGION_SUBSAMPLER_H
#define VCFX_REGION_SUBSAMPLER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXRegionSubsampler: Header file for Region-based Subsampling Tool
class VCFXRegionSubsampler {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads regions from a BED file into a map
    bool loadRegions(const std::string& bedFilePath, std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions);

    // Checks if a variant falls within any of the specified regions
    bool isVariantInRegions(const std::string& chrom, int pos, const std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions);

    // Subsamples variants from the VCF input based on regions
    void subsampleRegions(std::istream& in, std::ostream& out, const std::unordered_map<std::string, std::vector<std::pair<int, int>>>& regions);
};

#endif // VCFX_REGION_SUBSAMPLER_H


./VCFX_sample_extractor/VCFX_sample_extractor.cpp
#include "VCFX_sample_extractor.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_sample_extractor\n"
              << "Usage: VCFX_sample_extractor [OPTIONS]\n\n"
              << "Options:\n"
              << "  --sample, -s \"SampleName\"  Specify the sample name to extract data for.\n"
              << "  --help, -h                  Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Extracts genotype and related data for a specified sample from a VCF file.\n\n"
              << "Example:\n"
              << "  ./VCFX_sample_extractor --sample \"Sample1\" < input.vcf > sample1_data.tsv\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::string& sample_name) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--sample" || arg == "-s") && i + 1 < argc) {
            sample_name = argv[++i];
            return true;
        } else if (arg.find("--sample=") == 0) {
            sample_name = arg.substr(9);
            return true;
        }
    }
    return false;
}

// Function to extract sample data from VCF
void extractSampleData(std::istream& in, std::ostream& out, const std::string& sample_name) {
    std::string line;
    int sample_index = -1;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                // Header line with sample names
                std::stringstream ss(line);
                std::string field;
                int index = 0;
                while (std::getline(ss, field, '\t')) {
                    if (field == sample_name) {
                        sample_index = index;
                        break;
                    }
                    index++;
                }

                if (sample_index == -1) {
                    std::cerr << "Sample name '" << sample_name << "' not found in VCF header." << std::endl;
                    return;
                }

                // Print header
                out << "CHROM\tPOS\tID\tREF\tALT\t" << sample_name << "\n";
            }
            continue; // Skip other header lines
        }

        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (sample_index == -1 || sample_index >= static_cast<int>(fields.size())) {
            // Sample index not set or out of range
            continue;
        }

        // Extract required fields
        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string sample_data = fields[sample_index];

        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << sample_data << "\n";
    }
}

int main(int argc, char* argv[]) {
    std::string sample_name;
    if (!parseArguments(argc, argv, sample_name)) {
        std::cerr << "No sample name specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    extractSampleData(std::cin, std::cout, sample_name);
    return 0;
}


./VCFX_sample_extractor/VCFX_sample_extractor.h
#ifndef VCFX_SAMPLE_EXTRACTOR_H
#define VCFX_SAMPLE_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::string& sample_name);

// Function to extract sample data from VCF
void extractSampleData(std::istream& in, std::ostream& out, const std::string& sample_name);

// Function to display help message
void printHelp();

#endif // VCFX_SAMPLE_EXTRACTOR_H


./VCFX_sorter/VCFX_sorter.cpp
#include "VCFX_sorter.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_sorter\n"
              << "Usage: VCFX_sorter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h            Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Sorts VCF records based on chromosome and position.\n\n"
              << "Example:\n"
              << "  ./VCFX_sorter < unsorted.vcf > sorted.vcf\n";
}

// Implement comparator based on chromosome and position
bool VCFRecord::operator<(const VCFRecord& other) const {
    if (chrom != other.chrom) {
        return chrom < other.chrom;
    }
    return pos < other.pos;
}

// Function to parse a VCF line into a VCFRecord
bool parseVCFLine(const std::string& line, VCFRecord& record) {
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    while (std::getline(ss, field, '\t')) {
        fields.push_back(field);
    }

    if (fields.size() < 8) {
        return false; // Invalid VCF line
    }

    record.chrom = fields[0];
    try {
        record.pos = std::stoi(fields[1]);
    } catch (...) {
        return false; // Invalid position
    }
    record.id = fields[2];
    record.ref = fields[3];
    record.alt = fields[4];
    record.qual = fields[5];
    record.filter = fields[6];
    record.info = fields[7];

    // Parse sample fields if present
    if (fields.size() > 8) {
        record.samples.assign(fields.begin() + 8, fields.end());
    }

    return true;
}

// Function to sort VCF records
void sortVCFRecords(std::vector<VCFRecord>& records) {
    std::sort(records.begin(), records.end());
}

// Function to print sorted VCF records
void printSortedVCF(const std::vector<VCFRecord>& records, const std::string& header) {
    std::cout << header << "\n";
    for (const auto& record : records) {
        std::cout << record.chrom << "\t"
                  << record.pos << "\t"
                  << record.id << "\t"
                  << record.ref << "\t"
                  << record.alt << "\t"
                  << record.qual << "\t"
                  << record.filter << "\t"
                  << record.info;

        // Print sample fields if present
        if (!record.samples.empty()) {
            for (const auto& sample : record.samples) {
                std::cout << "\t" << sample;
            }
        }
        std::cout << "\n";
    }
}

int main(int argc, char* argv[]) {
    // Argument parsing for help
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    std::string line;
    std::string header;
    std::vector<VCFRecord> records;

    while (std::getline(std::cin, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            header = line; // Save header
            continue;
        }

        VCFRecord record;
        if (parseVCFLine(line, record)) {
            records.push_back(record);
        } else {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
        }
    }

    sortVCFRecords(records);
    printSortedVCF(records, header);

    return 0;
}


./VCFX_sorter/VCFX_sorter.h
#ifndef VCFX_SORTER_H
#define VCFX_SORTER_H

#include <iostream>
#include <string>
#include <vector>

// Structure to represent a VCF record
struct VCFRecord {
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::string alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::vector<std::string> samples;

    // Comparator for sorting
    bool operator<(const VCFRecord& other) const;
};

// Function to parse a VCF line into a VCFRecord
bool parseVCFLine(const std::string& line, VCFRecord& record);

// Function to sort VCF records
void sortVCFRecords(std::vector<VCFRecord>& records);

// Function to print sorted VCF records
void printSortedVCF(const std::vector<VCFRecord>& records, const std::string& header);

// Function to display help message
void printHelp();

#endif // VCFX_SORTER_H


./VCFX_subsampler/VCFX_subsampler.cpp
#include "VCFX_subsampler.h"
#include <sstream>
#include <vector>
#include <random>
#include <ctime>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_subsampler\n"
              << "Usage: VCFX_subsampler [OPTIONS]\n\n"
              << "Options:\n"
              << "  --subsample, -s <number>  Specify the number of variants to sample.\n"
              << "  --help, -h                Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Performs reservoir sampling on a VCF file to extract a subset of variants.\n\n"
              << "Example:\n"
              << "  ./VCFX_subsampler --subsample 1000 < input.vcf > sampled.vcf\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], int& sample_size) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--subsample" || arg == "-s") && i + 1 < argc) {
            try {
                sample_size = std::stoi(argv[++i]);
                if (sample_size <= 0) {
                    throw std::invalid_argument("Sample size must be positive.");
                }
                return true;
            } catch (...) {
                std::cerr << "Error: Invalid sample size provided.\n";
                return false;
            }
        } else if (arg.find("--subsample=") == 0) {
            try {
                sample_size = std::stoi(arg.substr(12));
                if (sample_size <= 0) {
                    throw std::invalid_argument("Sample size must be positive.");
                }
                return true;
            } catch (...) {
                std::cerr << "Error: Invalid sample size provided.\n";
                return false;
            }
        }
    }
    return false;
}

// Function to perform reservoir sampling on VCF records
void subsampleVariants(std::istream& in, std::ostream& out, int sample_size) {
    std::vector<std::string> reservoir;
    std::string line;
    int count = 0;

    // Preserve header lines
    std::vector<std::string> headers;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            headers.push_back(line);
            continue;
        }
        break;
    }

    // Print headers
    for (const auto& header : headers) {
        out << header << "\n";
    }

    if (line.empty() || line[0] == '#') {
        // No variant records
        return;
    }

    // Initialize reservoir with first sample_size records
    while (count < sample_size && !line.empty() && line[0] != '#') {
        reservoir.push_back(line);
        count++;
        if (!std::getline(in, line)) {
            break;
        }
    }

    // Continue with remaining records
    std::default_random_engine rng(static_cast<unsigned int>(std::time(nullptr)));
    std::uniform_int_distribution<int> dist;

    while (!line.empty() && line[0] != '#') {
        count++;
        dist = std::uniform_int_distribution<int>(0, count - 1);
        int j = dist(rng);
        if (j < sample_size) {
            reservoir[j] = line;
        }
        if (!std::getline(in, line)) {
            break;
        }
    }

    // Output the sampled records
    for (const auto& record : reservoir) {
        out << record << "\n";
    }
}

int main(int argc, char* argv[]) {
    // Argument parsing for help and subsample size
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    int sample_size = 0;
    if (!parseArguments(argc, argv, sample_size)) {
        std::cerr << "Usage: " << argv[0] << " --subsample <number_of_variants> < input.vcf > output.vcf\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    subsampleVariants(std::cin, std::cout, sample_size);
    return 0;
}


./VCFX_subsampler/VCFX_subsampler.h
#ifndef VCFX_SUBSAMPLER_H
#define VCFX_SUBSAMPLER_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], int& sample_size);

// Function to perform reservoir sampling on VCF records
void subsampleVariants(std::istream& in, std::ostream& out, int sample_size);

#endif // VCFX_SUBSAMPLER_H


./VCFX_sv_handler/VCFX_sv_handler.cpp
#include "VCFX_sv_handler.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

// Implementations

int VCFXSvHandler::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    bool filterOnly = false;
    bool modifySV = false;

    static struct option long_options[] = {
        {"help",            no_argument,       0, 'h'},
        {"sv-filter-only",  no_argument,       0, 'f'},
        {"sv-modify",       no_argument,       0, 'm'},
        {0,                 0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hfm", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                filterOnly = true;
                break;
            case 'm':
                modifySV = true;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Handle structural variants
    handleStructuralVariants(std::cin, std::cout, filterOnly, modifySV);

    return 0;
}

void VCFXSvHandler::displayHelp() {
    std::cout << "VCFX_sv_handler: Parse and manipulate structural variants (SVs) in a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_sv_handler [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help               Display this help message and exit\n";
    std::cout << "  -f, --sv-filter-only     Filter and output only structural variants\n";
    std::cout << "  -m, --sv-modify          Modify structural variant INFO fields\n\n";
    std::cout << "Examples:\n";
    std::cout << "  VCFX_sv_handler < input.vcf > output.vcf\n";
    std::cout << "  VCFX_sv_handler --sv-filter-only < input.vcf > sv_only.vcf\n";
    std::cout << "  VCFX_sv_handler --sv-modify < input.vcf > sv_modified.vcf\n";
}

bool VCFXSvHandler::isStructuralVariant(const std::string& infoField) {
    // Check if INFO field contains SVTYPE
    return infoField.find("SVTYPE=") != std::string::npos;
}

std::string VCFXSvHandler::parseSVType(const std::string& infoField) {
    size_t pos = infoField.find("SVTYPE=");
    if (pos == std::string::npos) {
        return "";
    }
    size_t start = pos + 7; // Length of "SVTYPE="
    size_t end = infoField.find(';', start);
    if (end == std::string::npos) {
        end = infoField.length();
    }
    return infoField.substr(start, end - start);
}

int VCFXSvHandler::parseEndPosition(const std::string& infoField) {
    size_t pos = infoField.find("END=");
    if (pos == std::string::npos) {
        return -1;
    }
    size_t start = pos + 4; // Length of "END="
    size_t end = infoField.find(';', start);
    std::string endStr = (end == std::string::npos) ? infoField.substr(start) : infoField.substr(start, end - start);
    try {
        return std::stoi(endStr);
    } catch (...) {
        return -1;
    }
}

int VCFXSvHandler::parsePos(const std::string& posField) {
    try {
        return std::stoi(posField);
    } catch (...) {
        return -1;
    }
}

std::string VCFXSvHandler::manipulateSVInfo(const std::string& infoField, const std::string& svType, int pos, int endPos) {
    // Example manipulation: Add a new INFO field indicating validation status
    std::string modifiedInfo = infoField;
    if (!modifiedInfo.empty() && modifiedInfo.back() != ';') {
        modifiedInfo += ";";
    }
    modifiedInfo += "SV_VALIDATED=1";

    // Add SV_SIZE if END is present and POS is valid
    if (endPos != -1 && pos != -1 && (svType == "DEL" || svType == "DUP")) {
        int svSize = endPos - pos;
        modifiedInfo += ";SV_SIZE=" + std::to_string(svSize);
    }

    // Further manipulations based on SV type
    if (svType == "INV") {
        // Example: Add inversion specific INFO field
        modifiedInfo += ";INV_TYPE=PARALLEL";
    } else if (svType == "BND") {
        // Example: Add breakend specific INFO field
        modifiedInfo += ";BND_ORIENTATION=PAIR";
    }

    // Additional SV types and their manipulations can be added here

    return modifiedInfo;
}

void VCFXSvHandler::handleStructuralVariants(std::istream& in, std::ostream& out, bool filterOnly, bool modifySV) {
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

        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string qual = fields[5];
        std::string filter = fields[6];
        std::string info = fields[7];

        // If there are genotype fields, they start from index 8
        std::vector<std::string> genotypeFields;
        if (fields.size() > 8) {
            genotypeFields.assign(fields.begin() + 8, fields.end());
        }

        if (isStructuralVariant(info)) {
            if (filterOnly && !modifySV) {
                // Output only structural variants
                out << line << "\n";
                continue;
            }

            if (modifySV) {
                // Modify SV INFO fields
                std::string svType = parseSVType(info);
                if (svType.empty()) {
                    std::cerr << "Warning: SVTYPE not found. Skipping variant.\n";
                    continue;
                }

                int posInt = parsePos(pos);
                int endPos = parseEndPosition(info);
                if (posInt == -1) {
                    std::cerr << "Warning: Invalid POS value. Skipping variant.\n";
                    continue;
                }

                std::string modifiedInfo = manipulateSVInfo(info, svType, posInt, endPos);

                // Reconstruct the VCF line with modified INFO
                std::ostringstream oss;
                oss << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
                    << "\t" << qual << "\t" << filter << "\t" << modifiedInfo;

                // Append genotype fields if any
                for (const auto& gt : genotypeFields) {
                    oss << "\t" << gt;
                }
                oss << "\n";

                out << oss.str();
                continue;
            }

            // If neither filtering nor modifying, default behavior (output as-is)
            out << line << "\n";
        } else {
            // If not a structural variant
            if (!filterOnly) {
                // Output non-SV records if not filtering
                out << line << "\n";
            }
            // Else, skip non-SV records when filtering
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXSvHandler svHandler;
    return svHandler.run(argc, argv);
}

./VCFX_sv_handler/VCFX_sv_handler.h
#ifndef VCFX_SV_HANDLER_H
#define VCFX_SV_HANDLER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXSvHandler: Header file for Structural Variant Handler Tool
class VCFXSvHandler {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Checks if a variant is a structural variant based on the INFO field
    bool isStructuralVariant(const std::string& infoField);

    // Parses the SVTYPE from the INFO field
    std::string parseSVType(const std::string& infoField);

    // Parses the END position from the INFO field
    int parseEndPosition(const std::string& infoField);

    // Parses the POS field and converts to integer
    int parsePos(const std::string& posField);

    // Manipulates the INFO field for structural variants
    std::string manipulateSVInfo(const std::string& infoField, const std::string& svType, int pos, int endPos);

    // Handles structural variants: filtering and modification
    void handleStructuralVariants(std::istream& in, std::ostream& out, bool filterOnly, bool modifySV);
};

#endif // VCFX_SV_HANDLER_H

./VCFX_validator/VCFX_validator.cpp
#include "VCFX_validator.h"
#include <sstream>
#include <vector>
#include <cctype>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_validator\n"
              << "Usage: VCFX_validator [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h            Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Validates the integrity and format of a VCF file.\n\n"
              << "Example:\n"
              << "  ./VCFX_validator < input.vcf\n";
}

// Function to trim whitespace from both ends of a string
static inline std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    if (start == std::string::npos)
        return "";
    return s.substr(start, end - start + 1);
}

// Function to validate VCF meta-information headers
bool validateVCFHeader(const std::string& line) {
    // Basic validation: check if header starts with ##
    if (line.size() < 2 || line[0] != '#' || line[1] != '#') {
        return false;
    }

    // Further validations can be implemented as per VCF specifications
    // For simplicity, we'll assume headers starting with ## are valid
    return true;
}

// Function to validate a single VCF record
bool validateVCFRecord(const std::string& line, int line_number) {
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    // Split the line by tabs
    while (std::getline(ss, field, '\t')) {
        fields.push_back(trim(field));
    }

    // VCF records should have at least 8 fields
    if (fields.size() < 8) {
        std::cerr << "Error: Line " << line_number << " has fewer than 8 fields.\n";
        return false;
    }

    // Validate CHROM: non-empty
    if (fields[0].empty()) {
        std::cerr << "Error: Line " << line_number << " has empty CHROM field.\n";
        return false;
    }

    // Validate POS: positive integer
    try {
        int pos = std::stoi(fields[1]);
        if (pos <= 0) {
            std::cerr << "Error: Line " << line_number << " has invalid POS value.\n";
            return false;
        }
    } catch (...) {
        std::cerr << "Error: Line " << line_number << " has invalid POS value.\n";
        return false;
    }

    // REF: non-empty
    if (fields[3].empty()) {
        std::cerr << "Error: Line " << line_number << " has empty REF field.\n";
        return false;
    }

    // ALT: non-empty
    if (fields[4].empty()) {
        std::cerr << "Error: Line " << line_number << " has empty ALT field.\n";
        return false;
    }

    // QUAL: Either '.' or a non-negative floating-point number
    if (fields[5] != ".") {
        try {
            float qual = std::stof(fields[5]);
            if (qual < 0) {
                std::cerr << "Error: Line " << line_number << " has negative QUAL value.\n";
                return false;
            }
        } catch (...) {
            std::cerr << "Error: Line " << line_number << " has invalid QUAL value.\n";
            return false;
        }
    }

    // FILTER: non-empty (can be "PASS" or specific filter names)
    if (fields[6].empty()) {
        std::cerr << "Error: Line " << line_number << " has empty FILTER field.\n";
        return false;
    }

    // INFO: can be empty but should follow key=value format separated by semicolons
    // Basic check: if not empty, contains at least one '=' or contains flags (no '=')
    if (!fields[7].empty()) {
        bool has_key_value = false;
        bool has_flag = false;
        std::stringstream info_ss(fields[7]);
        std::string info_field;
        while (std::getline(info_ss, info_field, ';')) {
            if (info_field.find('=') != std::string::npos) {
                has_key_value = true;
            } else {
                has_flag = true;
            }
        }
        if (!has_key_value && !has_flag) {
            std::cerr << "Error: Line " << line_number << " has invalid INFO field.\n";
            return false;
        }
    }

    // FORMAT and sample fields can be further validated if necessary
    // This implementation focuses on basic validations

    return true;
}

// Function to validate the entire VCF file
bool validateVCF(std::istream& in, std::ostream& out) {
    std::string line;
    int line_number = 0;
    bool header_found = false;
   
    while (std::getline(in, line)) {
        line_number++;
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_found = true;
                continue; // Header line found, continue to records
            } else {
                // Validate meta-information headers
                if (!validateVCFHeader(line)) {
                    std::cerr << "Error: Invalid VCF meta-information header at line " << line_number << ".\n";
                    return false;
                }
            }
            continue; // Continue processing headers
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records at line " << line_number << ".\n";
            return false;
        }

        // Validate the VCF record
        if (!validateVCFRecord(line, line_number)) {
            // Detailed error already printed in validateVCFRecord
            return false;
        }
    }

    if (!header_found) {
        std::cerr << "Error: VCF header (#CHROM) not found in the file.\n";
        return false;
    }

    std::cout << "VCF file is valid.\n";
    return true;
}

int main(int argc, char* argv[]) {
    // Argument parsing for help
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Validate VCF
    bool is_valid = validateVCF(std::cin, std::cout);
    return is_valid ? 0 : 1;
}

./VCFX_validator/VCFX_validator.h
#ifndef VCFX_VALIDATOR_H
#define VCFX_VALIDATOR_H

#include <iostream>
#include <string>

// Function to display help message
void printHelp();

// Function to validate VCF header
bool validateVCFHeader(const std::string& line);

// Function to validate a single VCF record
bool validateVCFRecord(const std::string& line, int line_number);

// Function to validate the entire VCF file
bool validateVCF(std::istream& in, std::ostream& out);

#endif // VCFX_VALIDATOR_H


./VCFX_variant_classifier/VCFX_variant_classifier.cpp
#include "VCFX_variant_classifier.h"
#include <sstream>
#include <algorithm>
#include <cctype>
#include <unordered_set>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_variant_classifier\n"
              << "Usage: VCFX_variant_classifier [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Classifies variants in a VCF file as SNPs, indels, MNVs, or structural variants based on the REF and ALT alleles.\n\n"
              << "Examples:\n"
              << "  ./VCFX_variant_classifier < input.vcf > classified_variants.tsv\n";
}

// Function to convert VariantType enum to string
std::string variantTypeToString(VariantType type) {
    switch (type) {
        case VariantType::SNP:
            return "SNP";
        case VariantType::INDEL:
            return "Indel";
        case VariantType::MNV:
            return "MNV";
        case VariantType::STRUCTURAL_VARIANT:
            return "Structural_Variant";
        default:
            return "Unknown";
    }
}

// Function to classify a single allele pair
VariantType classifyAllele(const std::string& ref, const std::string& alt) {
    // Check for symbolic alleles indicating structural variants
    if (!alt.empty() && alt[0] == '<' && alt.back() == '>') {
        return VariantType::STRUCTURAL_VARIANT;
    }

    // Check for simple SNP
    if (ref.length() == 1 && alt.length() == 1 &&
        std::isalpha(ref[0]) && std::isalpha(alt[0])) {
        return VariantType::SNP;
    }

    // Check for indel
    if (ref.length() != alt.length()) {
        // Length difference significant enough to consider as structural variant
        if (ref.length() > 50 || alt.length() > 50) { // Arbitrary threshold
            return VariantType::STRUCTURAL_VARIANT;
        }
        return VariantType::INDEL;
    }

    // Check for multi-nucleotide variant (MNV)
    if (ref.length() > 1 && alt.length() > 1) {
        return VariantType::MNV;
    }

    return VariantType::UNKNOWN;
}

// Function to classify a VCF variant
VariantType classifyVariant(const std::string& ref, const std::vector<std::string>& alt) {
    std::unordered_set<VariantType> types;

    for (const auto& allele : alt) {
        VariantType type = classifyAllele(ref, allele);
        types.insert(type);
    }

    // Prioritize variant types
    if (types.find(VariantType::STRUCTURAL_VARIANT) != types.end()) {
        return VariantType::STRUCTURAL_VARIANT;
    }
    if (types.find(VariantType::MNV) != types.end()) {
        return VariantType::MNV;
    }
    if (types.find(VariantType::INDEL) != types.end()) {
        return VariantType::INDEL;
    }
    if (types.find(VariantType::SNP) != types.end()) {
        return VariantType::SNP;
    }

    return VariantType::UNKNOWN;
}

// Function to split a string by a delimiter and trim whitespace
std::vector<std::string> splitAndTrim(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        // Trim whitespace
        size_t start = token.find_first_not_of(" \t");
        size_t end = token.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            tokens.emplace_back(token.substr(start, end - start + 1));
        } else if (start != std::string::npos) {
            tokens.emplace_back(token.substr(start));
        } else {
            tokens.emplace_back(""); // All whitespace
        }
    }
    return tokens;
}

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant) {
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;

    while (std::getline(ss, field, '\t')) {
        fields.emplace_back(field);
    }

    if (fields.size() < 8) {
        std::cerr << "Warning: Skipping invalid VCF line with fewer than 8 fields.\n";
        return false;
    }

    variant.chrom = fields[0];
    try {
        variant.pos = std::stoi(fields[1]);
    } catch (...) {
        std::cerr << "Warning: Invalid POS value. Skipping line.\n";
        return false;
    }
    variant.id = fields[2];
    variant.ref = fields[3];
    
    // ALT can have multiple alleles separated by ','
    variant.alt = splitAndTrim(fields[4], ',');
    variant.qual = fields[5];
    variant.filter = fields[6];
    variant.info = fields[7];

    // Classify variant
    variant.type = classifyVariant(variant.ref, variant.alt);

    return true;
}

// Function to classify variants in VCF
bool classifyVariants(std::istream& in, std::ostream& out) {
    std::string line;
    bool header_found = false;

    // Print header for classification
    out << "CHROM\tPOS\tID\tREF\tALT\tVARIANT_TYPE\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                header_found = true;
            }
            continue; // Skip header lines
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        VCFVariant variant;
        if (!parseVCFLine(line, variant)) {
            continue; // Skip invalid lines
        }

        // Prepare ALT field as comma-separated string
        std::string alt_str = "";
        for (size_t i = 0; i < variant.alt.size(); ++i) {
            alt_str += variant.alt[i];
            if (i != variant.alt.size() - 1) {
                alt_str += ",";
            }
        }

        out << variant.chrom << "\t"
            << variant.pos << "\t"
            << variant.id << "\t"
            << variant.ref << "\t"
            << alt_str << "\t"
            << variantTypeToString(variant.type) << "\n";
    }

    return true;
}

int main(int argc, char* argv[]) {
    // Argument parsing for help
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Classify variants
    bool success = classifyVariants(std::cin, std::cout);
    return success ? 0 : 1;
}


./VCFX_variant_classifier/VCFX_variant_classifier.h
#ifndef VCFX_VARIANT_CLASSIFIER_H
#define VCFX_VARIANT_CLASSIFIER_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Enumeration for variant types
enum class VariantType {
    SNP,
    INDEL,
    MNV,
    STRUCTURAL_VARIANT,
    UNKNOWN
};

// Function to convert VariantType enum to string
std::string variantTypeToString(VariantType type);

// Structure to represent a VCF variant
struct VCFVariant {
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::vector<std::string> alt;
    std::string qual;
    std::string filter;
    std::string info;
    VariantType type;
};

// Function to classify a single allele pair
VariantType classifyAllele(const std::string& ref, const std::string& alt);

// Function to classify a VCF variant
VariantType classifyVariant(const std::string& ref, const std::vector<std::string>& alt);

// Function to parse a VCF line into VCFVariant structure
bool parseVCFLine(const std::string& line, VCFVariant& variant);

// Function to classify variants in VCF
bool classifyVariants(std::istream& in, std::ostream& out);

#endif // VCFX_VARIANT_CLASSIFIER_H


./VCFX_variant_counter/VCFX_variant_counter.cpp
#include "VCFX_variant_counter.h"
#include <string>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_variant_counter\n"
              << "Usage: VCFX_variant_counter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h            Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Counts the total number of variants in a VCF file.\n\n"
              << "Example:\n"
              << "  ./VCFX_variant_counter < input.vcf > variant_count.txt\n";
}

// Function to count variants in VCF
int countVariants(std::istream& in) {
    int count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip header lines
        }
        count++;
    }
    return count;
}

int main(int argc, char* argv[]) {
    // Argument parsing for help
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Count variants
    int total_variants = countVariants(std::cin);
    std::cout << "Total Variants: " << total_variants << std::endl;
    return 0;
}


./VCFX_variant_counter/VCFX_variant_counter.h
#ifndef VCFX_VARIANT_COUNTER_H
#define VCFX_VARIANT_COUNTER_H

#include <iostream>
#include <string>

// Function to count variants in VCF
int countVariants(std::istream& in);

// Function to display help message
void printHelp();

#endif // VCFX_VARIANT_COUNTER_H


