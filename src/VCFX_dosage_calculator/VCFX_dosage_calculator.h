#ifndef VCFX_DOSAGE_CALCULATOR_H
#define VCFX_DOSAGE_CALCULATOR_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_dosage_calculator: Header file for Genotype Dosage Calculation tool
class VCFXDosageCalculator {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Calculates genotype dosage from VCF input and writes output
    void calculateDosage(std::istream& in, std::ostream& out);
};

#endif // VCFX_DOSAGE_CALCULATOR_H
