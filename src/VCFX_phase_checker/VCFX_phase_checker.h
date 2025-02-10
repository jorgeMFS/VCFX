#ifndef VCFX_PHASE_CHECKER_H
#define VCFX_PHASE_CHECKER_H

#include <iostream>
#include <string>
#include <vector>

/*
 * VCFXPhaseChecker:
 *   Reads a VCF, outputs only lines in which all samples' GT are fully phased.
 *   If any sample is unphased (or missing GT), logs a warning and discards line.
 */
class VCFXPhaseChecker {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    // Processes the input line-by-line. 
    // For each data line, if all samples are phased => write to stdout, else skip.
    void processVCF(std::istream &in, std::ostream &out);

    // Checks if a genotype string is "fully phased" (no '/' seen, multiple alleles separated by '|').
    // If GT is missing or partial => returns false (unphased).
    bool isFullyPhased(const std::string &gt) const;
};

#endif // VCFX_PHASE_CHECKER_H
