#ifndef VCFX_QUALITY_ADJUSTER_H
#define VCFX_QUALITY_ADJUSTER_H

#include <iostream>
#include <string>
#include <functional>
#include <unordered_map>

class VCFXQualityAdjuster {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();
    // parse the transformation name and build a transform function
    bool parseTransformationFunction(const std::string &funcStr,
                                     std::function<double(double)> &transFunc);

    // read lines from 'in', transform the QUAL field (6th col),
    // write lines to 'out'. If clamp is false, do not clamp negative or large.
    void adjustQualityScores(std::istream &in, std::ostream &out,
                             std::function<double(double)> transFunc,
                             bool clamp);

    // map of supported transformations
    std::unordered_map<std::string, std::function<double(double)>> supportedFunctions;

    // helper to initialize the map
    void initSupportedFunctions();
};

#endif // VCFX_QUALITY_ADJUSTER_H
