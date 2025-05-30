#ifndef VCFX_SV_HANDLER_H
#define VCFX_SV_HANDLER_H

#include <iostream>
#include <string>

class VCFXSvHandler {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();

    // The main method to read lines from 'in' and apply filtering/modify logic, then write to 'out'.
    void handleStructuralVariants(std::istream &in, std::ostream &out, bool filterOnly, bool modifySV);

    // Checks if a line's INFO indicates an SV (i.e. has "SVTYPE=").
    bool isStructuralVariant(const std::string &infoField) const;

    // Extract the SVTYPE=... substring from INFO; returns "" if not found.
    std::string parseSVType(const std::string &infoField) const;

    // Extract 'END=' from INFO; returns -1 if not found or invalid.
    int parseEndPosition(const std::string &infoField) const;

    // Parse position from string. Return -1 on error.
    int parsePos(const std::string &posField) const;

    // If we are modifying, do the manipulations (like adding SV_SIZE=..., etc.)
    std::string manipulateSVInfo(const std::string &infoField, const std::string &svType, int pos, int endPos) const;
};

#endif
