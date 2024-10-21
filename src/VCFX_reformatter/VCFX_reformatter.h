#ifndef VCFX_REFORMATTER_H
#define VCFX_REFORMATTER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// VCFXReformatter: Header file for VCF Reformatting Tool
class VCFXReformatter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Reformat VCF input based on user-specified rules
    void reformatVCF(std::istream& in, std::ostream& out,
                    const std::vector<std::string>& compressFields,
                    const std::vector<std::string>& reorderInfoFields,
                    const std::vector<std::string>& reorderFormatFields);

    // Compresses specified INFO or FORMAT fields
    std::string compressFieldsFunction(const std::string& fieldValue, const std::vector<std::string>& fieldsToCompress);

    // Reorders INFO fields based on user-specified order
    std::string reorderInfo(const std::string& infoField, const std::vector<std::string>& reorderOrder);

    // Reorders FORMAT fields based on user-specified order
    std::string reorderFormat(const std::string& formatField, const std::vector<std::string>& reorderOrder);
};

#endif // VCFX_REFORMATTER_H
