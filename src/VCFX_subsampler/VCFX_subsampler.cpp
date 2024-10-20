#include "VCFX_subsampler.h"
#include <sstream>
#include <vector>
#include <random>
#include <ctime>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_subsampler\n"
              << "Usage: VCFX_subsampler [OPTIONS]\n\n"
              << "Options:\n"
              << "  --subsample, -s <number>  Specify the number of variants to sample.\n"
              << "  --help, -h                Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Performs reservoir sampling on a VCF file to extract a subset of variants.\n\n"
              << "Example:\n"
              << "  ./VCFX_subsampler --subsample 1000 < input.vcf > sampled.vcf\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char* argv[], int& sample_size) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--subsample" || arg == "-s") && i + 1 < argc) {
            try {
                sample_size = std::stoi(argv[++i]);
                if (sample_size <= 0) {
                    throw std::invalid_argument("Sample size must be positive.");
                }
                return true;
            } catch (...) {
                std::cerr << "Error: Invalid sample size provided.\n";
                return false;
            }
        } else if (arg.find("--subsample=") == 0) {
            try {
                sample_size = std::stoi(arg.substr(12));
                if (sample_size <= 0) {
                    throw std::invalid_argument("Sample size must be positive.");
                }
                return true;
            } catch (...) {
                std::cerr << "Error: Invalid sample size provided.\n";
                return false;
            }
        }
    }
    return false;
}

// Function to perform reservoir sampling on VCF records
void subsampleVariants(std::istream& in, std::ostream& out, int sample_size) {
    std::vector<std::string> reservoir;
    std::string line;
    int count = 0;

    // Preserve header lines
    std::vector<std::string> headers;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            headers.push_back(line);
            continue;
        }
        break;
    }

    // Print headers
    for (const auto& header : headers) {
        out << header << "\n";
    }

    if (line.empty() || line[0] == '#') {
        // No variant records
        return;
    }

    // Initialize reservoir with first sample_size records
    while (count < sample_size && !line.empty() && line[0] != '#') {
        reservoir.push_back(line);
        count++;
        if (!std::getline(in, line)) {
            break;
        }
    }

    // Continue with remaining records
    std::default_random_engine rng(static_cast<unsigned int>(std::time(nullptr)));
    std::uniform_int_distribution<int> dist;

    while (!line.empty() && line[0] != '#') {
        count++;
        dist = std::uniform_int_distribution<int>(0, count - 1);
        int j = dist(rng);
        if (j < sample_size) {
            reservoir[j] = line;
        }
        if (!std::getline(in, line)) {
            break;
        }
    }

    // Output the sampled records
    for (const auto& record : reservoir) {
        out << record << "\n";
    }
}

int main(int argc, char* argv[]) {
    // Argument parsing for help and subsample size
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    int sample_size = 0;
    if (!parseArguments(argc, argv, sample_size)) {
        std::cerr << "Usage: " << argv[0] << " --subsample <number_of_variants> < input.vcf > output.vcf\n";
        std::cerr << "Use --help for usage information.\n";
        return 1;
    }

    subsampleVariants(std::cin, std::cout, sample_size);
    return 0;
}
