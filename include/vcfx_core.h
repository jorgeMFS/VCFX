#ifndef VCFX_CORE_H
#define VCFX_CORE_H

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace vcfx {

// Trim leading and trailing whitespace from a string
std::string trim(const std::string &str);

// Split a string on the given delimiter
std::vector<std::string> split(const std::string &str, char delimiter);

// Convenience helpers for printing common messages
void print_error(const std::string &msg, std::ostream &os = std::cerr);
void print_version(const std::string &tool, const std::string &version, std::ostream &os = std::cout);

inline std::string get_version() {
#ifdef VCFX_VERSION
    return VCFX_VERSION;
#else
    return "unknown";
#endif
}

inline bool handle_version_flag(int argc, char *argv[], const std::string &tool, std::ostream &os = std::cout) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--version") == 0 || std::strcmp(argv[i], "-v") == 0) {
            print_version(tool, get_version(), os);
            return true;
        }
    }
    return false;
}

// Check if a specific flag (long or short form) is present
bool flag_present(int argc, char *argv[], const char *long_flag, const char *short_flag = nullptr);

// Handle the --help flag using the provided callback. Returns true if the flag
// was found and handled.
inline bool handle_help_flag(int argc, char *argv[], void (*print_help)()) {
    if (flag_present(argc, argv, "--help", "-h")) {
        if (print_help)
            print_help();
        return true;
    }
    return false;
}

// Handle both --help and --version flags. Returns true if either flag was found
// and processed (in which case the caller should exit).
inline bool handle_common_flags(int argc, char *argv[], const std::string &tool, void (*print_help)(),
                                std::ostream &os = std::cout) {
    if (handle_help_flag(argc, argv, print_help))
        return true;
    return handle_version_flag(argc, argv, tool, os);
}

// Read entire input stream, automatically decompressing if gzip/BGZF
// compressed. Returns true on success and stores the resulting text in
// 'out'.
bool read_maybe_compressed(std::istream &in, std::string &out);

// Convenience helper to read a file that may be gzip/BGZF compressed. The file
// is loaded completely into memory and stored in 'out'. Returns true on
// success.
bool read_file_maybe_compressed(const std::string &path, std::string &out);

} // namespace vcfx

#endif // VCFX_CORE_H
