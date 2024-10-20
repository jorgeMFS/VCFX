#include "VCFX_sample_extractor.h"
#include <sstream>
#include <algorithm>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_sample_extractor\n"
              << "Usage: VCFX_sample_extractor [OPTIONS]\n\n"
              << "Options:\n"
              << "  --sample, -s \"SampleName\"  Specify the sample name to extract data for.\n"
              << "  --help, -h                  Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Extracts genotype and related data for a specified sample from a VCF file.\n\n"
              << "Example:\n"
              << "  ./VCFX_sample_extractor --sample \"Sample1\" < input.vcf > sample1_data.tsv\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], std::string& sample_name) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--sample" || arg == "-s") && i + 1 < argc) {
            sample_name = argv[++i];
            return true;
        } else if (arg.find("--sample=") == 0) {
            sample_name = arg.substr(9);
            return true;
        }
    }
    return false;
}

// Function to extract sample data from VCF
void extractSampleData(std::istream& in, std::ostream& out, const std::string& sample_name) {
    std::string line;
    int sample_index = -1;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.find("#CHROM") == 0) {
                // Header line with sample names
                std::stringstream ss(line);
                std::string field;
                int index = 0;
                while (std::getline(ss, field, '\t')) {
                    if (field == sample_name) {
                        sample_index = index;
                        break;
                    }
                    index++;
                }

                if (sample_index == -1) {
                    std::cerr << "Sample name '" << sample_name << "' not found in VCF header." << std::endl;
                    return;
                }

                // Print header
                out << "CHROM\tPOS\tID\tREF\tALT\t" << sample_name << "\n";
            }
            continue; // Skip other header lines
        }

        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        if (sample_index == -1 || sample_index >= static_cast<int>(fields.size())) {
            // Sample index not set or out of range
            continue;
        }

        // Extract required fields
        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string sample_data = fields[sample_index];

        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << sample_data << "\n";
    }
}

int main(int argc, char* argv[]) {
    std::string sample_name;
    if (!parseArguments(argc, argv, sample_name)) {
        std::cerr << "No sample name specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    extractSampleData(std::cin, std::cout, sample_name);
    return 0;
}
