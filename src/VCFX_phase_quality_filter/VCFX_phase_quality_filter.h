#ifndef VCFX_PHASE_QUALITY_FILTER_H
#define VCFX_PHASE_QUALITY_FILTER_H

#include <iostream>
#include <string>

class VCFXPhaseQualityFilter {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();
    // Takes the final operator (>, >=, <, <=, ==, !=) plus threshold, filters lines
    void filterByPQ(std::istream &in, std::ostream &out,
                    const std::string &op, double threshold);

    // Extracts the "PQ=" value from INFO or returns 0.0 if missing or invalid
    double parsePQScore(const std::string &info);

    // Internal helper to parse the condition string, e.g. "PQ>=30" => op=">=", thr=30
    bool parseCondition(const std::string &condition,
                        std::string &op, double &threshold);
};

#endif // VCFX_PHASE_QUALITY_FILTER_H
