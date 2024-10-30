#include "VCFX_file_splitter.h"
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <algorithm>

// Implementation of VCFXFileSplitter
int VCFXFileSplitter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string outputPrefix = "split";

    static struct option long_options[] = {
        {"help",    no_argument,       0, 'h'},
        {"prefix",  required_argument, 0, 'p'},
        {0,         0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hp:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'p':
                outputPrefix = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Split VCF by chromosome
    splitVCFByChromosome(std::cin, outputPrefix);

    return 0;
}

void VCFXFileSplitter::displayHelp() {
    std::cout << "VCFX_file_splitter: Split a large VCF file into smaller files based on chromosome.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_file_splitter [options] [--prefix <output_prefix>]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help           Display this help message and exit\n";
    std::cout << "  -p, --prefix <prefix> Specify the output file prefix (default: split)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_file_splitter --prefix \"chr\" < input.vcf\n";
}

void VCFXFileSplitter::splitVCFByChromosome(std::istream& in, const std::string& outputPrefix) {
    std::unordered_map<std::string, std::ofstream*> chromosomeFiles;
    std::string line;
    bool headerParsed = false;
    std::string headerLine;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.substr(0, 6) == "#CHROM") {
                headerLine = line;
                headerParsed = true;
            } else {
                // Write other header lines to all chromosome files
                for (auto& [chrom, ofs] : chromosomeFiles) {
                    if (ofs && ofs->is_open()) {
                        *ofs << line << "\n";
                    }
                }
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse chromosome from VCF line
        std::stringstream ss(line);
        std::string chrom;
        if (!std::getline(ss, chrom, '\t')) {
            std::cerr << "Warning: Unable to parse chromosome from line: " << line << "\n";
            continue;
        }

        // Check if file for this chromosome already exists
        if (chromosomeFiles.find(chrom) == chromosomeFiles.end()) {
            std::string filename = outputPrefix + "_" + chrom + ".vcf";
            chromosomeFiles[chrom] = new std::ofstream(filename);
            if (!chromosomeFiles[chrom]->is_open()) {
                std::cerr << "Error: Unable to create file: " << filename << "\n";
                delete chromosomeFiles[chrom];
                chromosomeFiles.erase(chrom);
                continue;
            }
            // Write header to the new chromosome file
            *chromosomeFiles[chrom] << headerLine << "\n";
        }

        // Write the line to the corresponding chromosome file
        if (chromosomeFiles[chrom] && chromosomeFiles[chrom]->is_open()) {
            *chromosomeFiles[chrom] << line << "\n";
        } else {
            std::cerr << "Warning: File stream for chromosome " << chrom << " is not open.\n";
        }
    }

    // Close all chromosome files
    for (auto& [chrom, ofs] : chromosomeFiles) {
        if (ofs && ofs->is_open()) {
            ofs->close();
        }
        delete ofs;
    }
}

int main(int argc, char* argv[]) {
    VCFXFileSplitter fileSplitter;
    return fileSplitter.run(argc, argv);
}
