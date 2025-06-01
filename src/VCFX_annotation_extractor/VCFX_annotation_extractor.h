#ifndef VCFX_ANNOTATION_EXTRACTOR_H
#define VCFX_ANNOTATION_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_annotation_extractor: Header file for variant annotation extraction tool
class VCFXAnnotationExtractor {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Processes VCF input from the given stream
    void processVCF(std::istream &in, const std::vector<std::string> &annotations);

    // Parses annotations from the INFO field
    std::vector<std::string> parseINFO(const std::string &info);

    // Extracts specified annotations
    std::vector<std::string> extractAnnotations(const std::vector<std::string> &info_fields,
                                                const std::vector<std::string> &annotations);
};

#endif // VCFX_ANNOTATION_EXTRACTOR_H