#include "VCFX_variant_counter.h"
#include <string>

// Function to count variants in VCF
int countVariants(std::istream& in) {
    int count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip header lines
        }
        count++;
    }
    return count;
}

int main(int argc, char* argv[]) {
    int total_variants = countVariants(std::cin);
    std::cout << "Total Variants: " << total_variants << std::endl;
    return 0;
}
