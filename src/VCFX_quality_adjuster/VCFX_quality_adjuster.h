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
