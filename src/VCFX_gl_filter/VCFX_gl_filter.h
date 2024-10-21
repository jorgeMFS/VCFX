#ifndef VCFX_GL_FILTER_H
#define VCFX_GL_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXGLFilter: Header file for Genotype Likelihood Filter tool
class VCFXGLFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on genotype likelihood scores (e.g., GQ > threshold)
    void filterByGL(std::istream& in, std::ostream& out, const std::string& filterCondition);
};

#endif // VCFX_GL_FILTER_H
