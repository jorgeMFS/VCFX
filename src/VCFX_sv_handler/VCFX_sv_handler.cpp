#include "VCFX_sv_handler.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

// Implementations

int VCFXSvHandler::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    bool filterOnly = false;
    bool modifySV = false;

    static struct option long_options[] = {
        {"help",            no_argument,       0, 'h'},
        {"sv-filter-only",  no_argument,       0, 'f'},
        {"sv-modify",       no_argument,       0, 'm'},
        {0,                 0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hfm", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                filterOnly = true;
                break;
            case 'm':
                modifySV = true;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Handle structural variants
    handleStructuralVariants(std::cin, std::cout, filterOnly, modifySV);

    return 0;
}

void VCFXSvHandler::displayHelp() {
    std::cout << "VCFX_sv_handler: Parse and manipulate structural variants (SVs) in a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_sv_handler [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help               Display this help message and exit\n";
    std::cout << "  -f, --sv-filter-only     Filter and output only structural variants\n";
    std::cout << "  -m, --sv-modify          Modify structural variant INFO fields\n\n";
    std::cout << "Examples:\n";
    std::cout << "  VCFX_sv_handler < input.vcf > output.vcf\n";
    std::cout << "  VCFX_sv_handler --sv-filter-only < input.vcf > sv_only.vcf\n";
    std::cout << "  VCFX_sv_handler --sv-modify < input.vcf > sv_modified.vcf\n";
}

bool VCFXSvHandler::isStructuralVariant(const std::string& infoField) {
    // Check if INFO field contains SVTYPE
    return infoField.find("SVTYPE=") != std::string::npos;
}

std::string VCFXSvHandler::parseSVType(const std::string& infoField) {
    size_t pos = infoField.find("SVTYPE=");
    if (pos == std::string::npos) {
        return "";
    }
    size_t start = pos + 7; // Length of "SVTYPE="
    size_t end = infoField.find(';', start);
    if (end == std::string::npos) {
        end = infoField.length();
    }
    return infoField.substr(start, end - start);
}

int VCFXSvHandler::parseEndPosition(const std::string& infoField) {
    size_t pos = infoField.find("END=");
    if (pos == std::string::npos) {
        return -1;
    }
    size_t start = pos + 4; // Length of "END="
    size_t end = infoField.find(';', start);
    std::string endStr = (end == std::string::npos) ? infoField.substr(start) : infoField.substr(start, end - start);
    try {
        return std::stoi(endStr);
    } catch (...) {
        return -1;
    }
}

int VCFXSvHandler::parsePos(const std::string& posField) {
    try {
        return std::stoi(posField);
    } catch (...) {
        return -1;
    }
}

std::string VCFXSvHandler::manipulateSVInfo(const std::string& infoField, const std::string& svType, int pos, int endPos) {
    // Example manipulation: Add a new INFO field indicating validation status
    std::string modifiedInfo = infoField;
    if (!modifiedInfo.empty() && modifiedInfo.back() != ';') {
        modifiedInfo += ";";
    }
    modifiedInfo += "SV_VALIDATED=1";

    // Add SV_SIZE if END is present and POS is valid
    if (endPos != -1 && pos != -1 && (svType == "DEL" || svType == "DUP")) {
        int svSize = endPos - pos;
        modifiedInfo += ";SV_SIZE=" + std::to_string(svSize);
    }

    // Further manipulations based on SV type
    if (svType == "INV") {
        // Example: Add inversion specific INFO field
        modifiedInfo += ";INV_TYPE=PARALLEL";
    } else if (svType == "BND") {
        // Example: Add breakend specific INFO field
        modifiedInfo += ";BND_ORIENTATION=PAIR";
    }

    // Additional SV types and their manipulations can be added here

    return modifiedInfo;
}

void VCFXSvHandler::handleStructuralVariants(std::istream& in, std::ostream& out, bool filterOnly, bool modifySV) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // Header lines are printed as-is
            out << line << "\n";
            continue;
        }

        // Split the line into fields
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields): " << line << "\n";
            continue;
        }

        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string qual = fields[5];
        std::string filter = fields[6];
        std::string info = fields[7];

        // If there are genotype fields, they start from index 8
        std::vector<std::string> genotypeFields;
        if (fields.size() > 8) {
            genotypeFields.assign(fields.begin() + 8, fields.end());
        }

        if (isStructuralVariant(info)) {
            if (filterOnly && !modifySV) {
                // Output only structural variants
                out << line << "\n";
                continue;
            }

            if (modifySV) {
                // Modify SV INFO fields
                std::string svType = parseSVType(info);
                if (svType.empty()) {
                    std::cerr << "Warning: SVTYPE not found. Skipping variant.\n";
                    continue;
                }

                int posInt = parsePos(pos);
                int endPos = parseEndPosition(info);
                if (posInt == -1) {
                    std::cerr << "Warning: Invalid POS value. Skipping variant.\n";
                    continue;
                }

                std::string modifiedInfo = manipulateSVInfo(info, svType, posInt, endPos);

                // Reconstruct the VCF line with modified INFO
                std::ostringstream oss;
                oss << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
                    << "\t" << qual << "\t" << filter << "\t" << modifiedInfo;

                // Append genotype fields if any
                for (const auto& gt : genotypeFields) {
                    oss << "\t" << gt;
                }
                oss << "\n";

                out << oss.str();
                continue;
            }

            // If neither filtering nor modifying, default behavior (output as-is)
            out << line << "\n";
        } else {
            // If not a structural variant
            if (!filterOnly) {
                // Output non-SV records if not filtering
                out << line << "\n";
            }
            // Else, skip non-SV records when filtering
        }
    }
}