#include "VCFX_header_parser.h"
#include <iostream>
#include <sstream>

void processHeader(std::istream& in, std::ostream& out) {
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') {
            out << line << std::endl;
        } else {
            break; // Stop reading after header
        }
    }
}

int main(int argc, char* argv[]) {
    processHeader(std::cin, std::cout);
    return 0;
}
