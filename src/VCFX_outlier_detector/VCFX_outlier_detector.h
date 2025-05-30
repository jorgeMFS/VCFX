#ifndef VCFX_OUTLIER_DETECTOR_H
#define VCFX_OUTLIER_DETECTOR_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class VCFXOutlierDetector {
  public:
    // Main entry point
    int run(int argc, char *argv[]);

  private:
    // Print usage
    void displayHelp();

    // The function that does the reading and analysis
    void detectOutliers(std::istream &in, std::ostream &out, const std::string &metric, double threshold,
                        bool isVariantMode);

    // Parse the user-specified metric from INFO
    bool parseMetricFromInfo(const std::string &info, const std::string &key, double &val) const;

    // Parse the user-specified metric from genotype subfields
    bool parseMetricFromGenotype(const std::string &genotypeField, const std::string &metric, double &value) const;
};

#endif
