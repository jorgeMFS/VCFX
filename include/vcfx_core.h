#ifndef VCFX_CORE_H
#define VCFX_CORE_H

#include <iostream>
#include <string>
#include <vector>

namespace vcfx {

// Trim leading and trailing whitespace from a string
std::string trim(const std::string& str);

// Split a string on the given delimiter
std::vector<std::string> split(const std::string& str, char delimiter);

// Convenience helpers for printing common messages
void print_error(const std::string& msg, std::ostream& os = std::cerr);
void print_version(const std::string& tool, const std::string& version,
                   std::ostream& os = std::cout);

}  // namespace vcfx

#endif // VCFX_CORE_H
