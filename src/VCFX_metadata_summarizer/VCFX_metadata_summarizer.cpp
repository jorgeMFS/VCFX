#include "VCFX_metadata_summarizer.h"
#include <getopt.h>
#include <sstream>

// Implementation of VCFX_metadata_summarizer
int VCFXMetadataSummarizer::run(int argc, char* argv[]) {
    // Parse command-line arguments (if any in future extensions)
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {0,      0,           0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Summarize metadata from stdin
    summarizeMetadata(std::cin);

    return 0;
}

void VCFXMetadataSummarizer::displayHelp() {
    std::cout << "VCFX_metadata_summarizer: Summarize key metadata from a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_metadata_summarizer [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_metadata_summarizer < input.vcf\n";
}

void VCFXMetadataSummarizer::summarizeMetadata(std::istream& in) {
    std::string line;
    std::map<std::string, int> metadata;
    int numSamples = 0;
    int numVariants = 0;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line.find("##") == 0) {
                // Meta-information lines
                parseHeader(line, metadata);
            } else if (line.find("#CHROM") == 0) {
                // Header line with column names and sample information
                std::vector<std::string> fields;
                std::stringstream ss(line);
                std::string field;
                while (std::getline(ss, field, '\t')) {
                    fields.push_back(field);
                }
                numSamples = static_cast<int>(fields.size()) - 9; // Standard VCF has 9 fixed columns
                metadata["Number of Samples"] = numSamples;
            }
            continue;
        }

        // Count variants
        numVariants++;
    }

    metadata["Number of Variants"] = numVariants;

    // Additional summaries can be added here (e.g., contigs, variant types)
    
    printSummary(metadata);
}

void VCFXMetadataSummarizer::parseHeader(const std::string& line, std::map<std::string, int>& metadata) {
    // Example meta-information:
    // ##fileformat=VCFv4.2
    // ##source=...
    // ##contig=<ID=1,length=248956422>
    if (line.find("##contig=") != std::string::npos) {
        metadata["Number of Contigs"]++;
    }
    if (line.find("##INFO=") != std::string::npos) {
        metadata["Number of INFO Fields"]++;
    }
    if (line.find("##FILTER=") != std::string::npos) {
        metadata["Number of FILTER Fields"]++;
    }
    if (line.find("##FORMAT=") != std::string::npos) {
        metadata["Number of FORMAT Fields"]++;
    }
}

void VCFXMetadataSummarizer::printSummary(const std::map<std::string, int>& metadata) {
    std::cout << "VCF Metadata Summary:\n";
    std::cout << "---------------------\n";
    for (const auto& [key, value] : metadata) {
        std::cout << key << ": " << value << "\n";
    }
}

int main(int argc, char* argv[]) {
    VCFXMetadataSummarizer summarizer;
    return summarizer.run(argc, argv);
}
