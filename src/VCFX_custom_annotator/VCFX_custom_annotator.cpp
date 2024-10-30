#include "VCFX_custom_annotator.h"
#include <getopt.h>
#include <sstream>
#include <fstream>

// Implementation of VCFXCustomAnnotator
int VCFXCustomAnnotator::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string annotationFilePath;

    static struct option long_options[] = {
        {"help",          no_argument,       0, 'h'},
        {"add-annotation", required_argument, 0, 'a'},
        {0,               0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'a':
                annotationFilePath = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp || annotationFilePath.empty()) {
        displayHelp();
        return 1;
    }

    // Load annotations
    std::unordered_map<std::string, std::string> annotations;
    if (!loadAnnotations(annotationFilePath, annotations)) {
        std::cerr << "Error: Failed to load annotations from " << annotationFilePath << "\n";
        return 1;
    }

    // Add annotations to VCF
    addAnnotations(std::cin, std::cout, annotations);

    return 0;
}

void VCFXCustomAnnotator::displayHelp() {
    std::cout << "VCFX_custom_annotator: Add custom annotations to the INFO field in a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_custom_annotator --add-annotation <annotations.txt> [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                  Display this help message and exit\n";
    std::cout << "  -a, --add-annotation <file> Specify the annotation file\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_custom_annotator --add-annotation annotations.txt < input.vcf > annotated.vcf\n";
}

bool VCFXCustomAnnotator::loadAnnotations(const std::string& annotationFilePath, std::unordered_map<std::string, std::string>& annotations) {
    std::ifstream infile(annotationFilePath);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open annotation file " << annotationFilePath << "\n";
        return false;
    }

    std::string line;
    size_t line_num = 0;
    while (std::getline(infile, line)) {
        line_num++;
        if (line.empty() || line[0] == '#') {
            continue; // Skip comments and empty lines
        }

        std::stringstream ss(line);
        std::string chrom, pos, ref, alt, annotation;

        // Assuming the annotation file is tab-delimited: chrom, pos, ref, alt, annotation
        if (!(ss >> chrom >> pos >> ref >> alt)) {
            std::cerr << "Warning: Skipping invalid annotation line " << line_num << ": " << line << "\n";
            continue;
        }

        // The rest of the line is the annotation value
        if (!std::getline(ss, annotation, '\n')) {
            std::cerr << "Warning: Missing annotation value in line " << line_num << ": " << line << "\n";
            continue;
        }

        // Trim leading whitespace from annotation
        annotation.erase(0, annotation.find_first_not_of(" \t"));

        std::string key = generateVariantKey(chrom, pos, ref, alt);
        annotations[key] = annotation;
    }

    infile.close();
    return true;
}

std::string VCFXCustomAnnotator::generateVariantKey(const std::string& chrom, const std::string& pos, const std::string& ref, const std::string& alt) {
    return chrom + ":" + pos + ":" + ref + ":" + alt;
}

void VCFXCustomAnnotator::addAnnotations(std::istream& in, std::ostream& out, const std::unordered_map<std::string, std::string>& annotations) {
    std::string line;
    bool headerPassed = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // Add new INFO header line if not already present
            if (line.substr(0, 6) == "#CHROM") {
                out << "##INFO=<ID=CustomAnnotation,Number=.,Type=String,Description=\"Custom annotations added by VCFX_custom_annotator\">\n";
            }
            out << line << "\n";
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info, format;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        // Generate key and check for annotation
        std::string key = generateVariantKey(chrom, pos, ref, alt);
        auto it = annotations.find(key);
        if (it != annotations.end()) {
            // Append the custom annotation to the INFO field
            if (info != ".") {
                info += ";";
            }
            info += "CustomAnnotation=" + it->second;
        }

        // Reconstruct the VCF line with updated INFO field
        std::string rest_of_line;
        getline(ss, rest_of_line);
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t"
            << qual << "\t" << filter << "\t" << info << "\t" << format << rest_of_line << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXCustomAnnotator customAnnotator;
    return customAnnotator.run(argc, argv);
}