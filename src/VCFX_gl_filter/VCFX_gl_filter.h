#ifndef VCFX_GL_FILTER_H
#define VCFX_GL_FILTER_H

#include <iostream>
#include <string>

// VCFXGLFilter: Filters VCF records by a genotype-likelihood field (e.g. "GQ>20").
class VCFXGLFilter {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on genotype-likelihood expression (stdin fallback)
    void filterByGL(std::istream &in, std::ostream &out, const std::string &filterCondition, bool anyMode);

    // Memory-mapped file processing (fast path)
    bool filterByGLMmap(const char* filepath, std::ostream& out,
                         const std::string& field, int opTypeInt,
                         double threshold, bool anyMode);
};

#endif // VCFX_GL_FILTER_H