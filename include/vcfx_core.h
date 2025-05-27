#ifndef VCFX_CORE_H
#define VCFX_CORE_H

#include <iostream>
#include <string>
#include <vector>
#include <cstring>

namespace vcfx {

// Trim leading and trailing whitespace from a string
std::string trim(const std::string& str);

// Split a string on the given delimiter
std::vector<std::string> split(const std::string& str, char delimiter);

// Convenience helpers for printing common messages
void print_error(const std::string& msg, std::ostream& os = std::cerr);
void print_version(const std::string& tool, const std::string& version,
                   std::ostream& os = std::cout);

inline std::string get_version() {
#ifdef VCFX_VERSION
    return VCFX_VERSION;
#else
    return "unknown";
#endif
}

inline bool handle_version_flag(int argc, char* argv[], const std::string& tool,
                                std::ostream& os = std::cout) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--version") == 0 ||
            std::strcmp(argv[i], "-v") == 0) {
            print_version(tool, get_version(), os);
            return true;
        }
    }
    return false;
}

// Handle common flags shared by all tools. Returns true if "--version" or
// "--help"/"-h" was processed and the program should exit.
bool handle_common_flags(int argc, char* argv[],
                         const std::string& tool,
                         void (*print_help)(),
                         std::ostream& os = std::cout);

// Read entire input stream, automatically decompressing if gzip/BGZF
// compressed. Returns true on success and stores the resulting text in
// 'out'.
bool read_maybe_compressed(std::istream& in, std::string& out);

// Convenience helper to read a file that may be gzip/BGZF compressed. The file
// is loaded completely into memory and stored in 'out'. Returns true on
// success.
bool read_file_maybe_compressed(const std::string& path, std::string& out);

}  // namespace vcfx

#endif // VCFX_CORE_H
