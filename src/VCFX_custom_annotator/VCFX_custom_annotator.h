#ifndef VCFX_CUSTOM_ANNOTATOR_H
#define VCFX_CUSTOM_ANNOTATOR_H

#include <iostream>
#include <string>
#include <unordered_map>

// VCFXCustomAnnotator: Header file for Custom Annotation Addition Tool
class VCFXCustomAnnotator {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Loads annotations from a file into a map
    bool loadAnnotations(const std::string& annotationFilePath, std::unordered_map<std::string, std::string>& annotations);

    // Adds annotations to the VCF input
    void addAnnotations(std::istream& in, std::ostream& out, const std::unordered_map<std::string, std::string>& annotations);

    // Generates a unique key for a variant based on chromosome, position, ref, and alt
    std::string generateVariantKey(const std::string& chrom, const std::string& pos, const std::string& ref, const std::string& alt);
};

#endif // VCFX_CUSTOM_ANNOTATOR_H
