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
