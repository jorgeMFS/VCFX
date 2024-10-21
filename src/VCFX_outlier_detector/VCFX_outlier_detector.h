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
