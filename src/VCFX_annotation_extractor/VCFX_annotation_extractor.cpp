#include "VCFX_annotation_extractor.h"
#include <getopt.h>
#include <sstream>

// Implementation of VCFX_annotation_extractor
int VCFXAnnotationExtractor::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    std::vector<std::string> annotations;
    bool showHelp = false;

    static struct option long_options[] = {
        {"annotation-extract", required_argument, 0, 'a'},
        {"help",               no_argument,       0, 'h'},
        {0,                    0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "a:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'a':
                {
                    std::stringstream ss(optarg);
                    std::string annot;
                    while (std::getline(ss, annot, ',')) {
                        annotations.emplace_back(annot);
                    }
                }
                break;
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp || annotations.empty()) {
        displayHelp();
        return 0;
    }

    // Process VCF input from stdin
    processVCF(std::cin, annotations);

    return 0;
}

void VCFXAnnotationExtractor::displayHelp() {
    std::cout << "VCFX_annotation_extractor: Extract variant annotations from a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_annotation_extractor --annotation-extract \"ANN,Gene\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -a, --annotation-extract   Comma-separated list of annotations to extract (e.g., ANN,Gene)\n";
    std::cout << "  -h, --help                 Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_annotation_extractor --annotation-extract \"ANN,Gene\" < input.vcf\n";
}

void VCFXAnnotationExtractor::processVCF(std::istream& in, const std::vector<std::string>& annotations) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            std::cout << line << "\n"; // Preserve header lines
            continue;
        }

        // Split the VCF line into fields
        std::vector<std::string> fields;
        std::string field;
        size_t pos = 0;
        while ((pos = line.find('\t')) != std::string::npos) {
            field = line.substr(0, pos);
            fields.push_back(field);
            line.erase(0, pos + 1);
        }
        fields.push_back(line); // Add the last field

        if (fields.size() < 8) {
            std::cerr << "Invalid VCF line with fewer than 8 fields.\n";
            continue;
        }

        std::string info = fields[7];
        std::vector<std::string> info_fields = parseINFO(info);
        std::vector<std::string> extracted = extractAnnotations(info_fields, annotations);

        // Output the extracted annotations
        for (const auto& annot : extracted) {
            std::cout << annot << "\t";
        }
        std::cout << "\n";
    }
}

std::vector<std::string> VCFXAnnotationExtractor::parseINFO(const std::string& info) {
    std::vector<std::string> info_fields;
    std::string field;
    size_t pos = 0;
    std::string temp = info;
    while ((pos = temp.find(';')) != std::string::npos) {
        field = temp.substr(0, pos);
        info_fields.push_back(field);
        temp.erase(0, pos + 1);
    }
    info_fields.push_back(temp);
    return info_fields;
}

std::vector<std::string> VCFXAnnotationExtractor::extractAnnotations(const std::vector<std::string>& info_fields, const std::vector<std::string>& annotations) {
    std::vector<std::string> extracted;
    for (const auto& annot : annotations) {
        bool found = false;
        for (const auto& field : info_fields) {
            if (field.find(annot + "=") == 0) {
                extracted.push_back(field.substr(annot.length() + 1));
                found = true;
                break;
            }
        }
        if (!found) {
            extracted.push_back("NA");
        }
    }
    return extracted;
}

int main(int argc, char* argv[]) {
    VCFXAnnotationExtractor extractor;
    return extractor.run(argc, argv);
}