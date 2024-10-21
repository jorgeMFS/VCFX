#ifndef VCFX_SV_HANDLER_H
#define VCFX_SV_HANDLER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXSvHandler: Header file for Structural Variant Handler Tool
class VCFXSvHandler {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Checks if a variant is a structural variant based on the INFO field
    bool isStructuralVariant(const std::string& infoField);

    // Parses the SVTYPE from the INFO field
    std::string parseSVType(const std::string& infoField);

    // Parses the END position from the INFO field
    int parseEndPosition(const std::string& infoField);

    // Parses the POS field and converts to integer
    int parsePos(const std::string& posField);

    // Manipulates the INFO field for structural variants
    std::string manipulateSVInfo(const std::string& infoField, const std::string& svType, int pos, int endPos);

    // Handles structural variants: filtering and modification
    void handleStructuralVariants(std::istream& in, std::ostream& out, bool filterOnly, bool modifySV);
};

#endif // VCFX_SV_HANDLER_H