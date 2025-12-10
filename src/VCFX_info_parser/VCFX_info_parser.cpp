#include "VCFX_info_parser.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <vector>

#ifdef __unix__
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __AVX2__
#include <immintrin.h>
#endif

// Function to display help message
void printHelp() {
    std::cout
        << "VCFX_info_parser\n"
        << "Usage: VCFX_info_parser [OPTIONS]\n\n"
        << "Options:\n"
        << "  -i, --info \"FIELD1,FIELD2\"   Specify the INFO fields to display (e.g., \"DP,AF\").\n"
        << "  -I, --input FILE             Input VCF file (default: stdin).\n"
        << "  -q, --quiet                  Suppress warnings.\n"
        << "  -h, --help                   Display this help message and exit.\n\n"
        << "Description:\n"
        << "  Parses the INFO field of a VCF file and displays the selected INFO fields in a user-friendly format.\n\n"
        << "Examples:\n"
        << "  VCFX_info_parser --info \"DP,AF\" < input.vcf > output_info.tsv\n"
        << "  VCFX_info_parser -i \"DP,AF\" -I input.vcf > output_info.tsv\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char *argv[], std::vector<std::string> &info_fields) {
    bool foundAnyField = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if ((arg == "--info" || arg == "-i") && i + 1 < argc) {
            std::string fields_str = argv[++i];
            std::stringstream ss(fields_str);
            std::string field;
            while (std::getline(ss, field, ',')) {
                // Trim whitespace
                field.erase(0, field.find_first_not_of(" \t\n\r\f\v"));
                field.erase(field.find_last_not_of(" \t\n\r\f\v") + 1);
                if (!field.empty()) {
                    info_fields.push_back(field);
                    foundAnyField = true;
                }
            }
        } else if (arg.rfind("--info=", 0) == 0) {
            // e.g. --info=DP,AF
            std::string fields_str = arg.substr(7);
            std::stringstream ss(fields_str);
            std::string field;
            while (std::getline(ss, field, ',')) {
                // Trim
                field.erase(0, field.find_first_not_of(" \t\n\r\f\v"));
                field.erase(field.find_last_not_of(" \t\n\r\f\v") + 1);
                if (!field.empty()) {
                    info_fields.push_back(field);
                    foundAnyField = true;
                }
            }
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        }
        // ignore other unrecognized args...
    }

    return foundAnyField;
}

// Splits a string by a delimiter (optimized with find-based approach)
std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    tokens.reserve(8);
    size_t start = 0;
    size_t end;
    while ((end = s.find(delimiter, start)) != std::string::npos) {
        tokens.emplace_back(s, start, end - start);
        start = end + 1;
    }
    tokens.emplace_back(s, start);
    return tokens;
}

// Parses the INFO field and display selected fields
bool parseInfoFields(std::istream &in, std::ostream &out, const std::vector<std::string> &info_fields) {
    // If we have any fields, print header
    if (!info_fields.empty()) {
        out << "CHROM\tPOS\tID\tREF\tALT";
        for (const auto &field : info_fields) {
            out << "\t" << field;
        }
        out << "\n";
    }

    std::string line;

    // Performance: reuse containers across iterations
    std::vector<std::string> fields;
    fields.reserve(16);
    std::unordered_map<std::string, std::string> info_map;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // skip header lines
            continue;
        }

        // Performance: use find-based tab splitting instead of stringstream
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        // parse INFO => key-value or flag
        info_map.clear();
        auto info_entries = split(fields[7], ';');
        for (auto &entry : info_entries) {
            size_t eqPos = entry.find('=');
            if (eqPos != std::string::npos) {
                std::string key = entry.substr(0, eqPos);
                std::string value = entry.substr(eqPos + 1);
                info_map[key] = value;
            } else {
                // This is a flag (no '=' => empty string)
                info_map[entry] = "";
            }
        }

        // Print row
        out << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3] << "\t" << fields[4];
        for (const auto &field : info_fields) {
            auto it = info_map.find(field);
            if (it != info_map.end()) {
                // If it's a flag (empty string), print "."
                if (it->second.empty()) {
                    out << "\t.";
                } else {
                    out << "\t" << it->second;
                }
            } else {
                // Field not present
                out << "\t.";
            }
        }
        out << "\n";
    }
    return true;
}

// ============================================================================
// MMAP-based high-performance implementation
// ============================================================================

#ifdef __unix__

namespace {

// RAII wrapper for memory-mapped files
struct MappedFile {
    const char* data = nullptr;
    size_t size = 0;
    int fd = -1;

    MappedFile() = default;
    ~MappedFile() {
        if (data && data != MAP_FAILED) {
            munmap(const_cast<char*>(data), size);
        }
        if (fd >= 0) {
            close(fd);
        }
    }
    MappedFile(const MappedFile&) = delete;
    MappedFile& operator=(const MappedFile&) = delete;

    bool open(const char* filepath) {
        fd = ::open(filepath, O_RDONLY);
        if (fd < 0) return false;

        struct stat st;
        if (fstat(fd, &st) < 0) return false;
        size = static_cast<size_t>(st.st_size);
        if (size == 0) return true;

        data = static_cast<const char*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            data = nullptr;
            return false;
        }

        // Hint for sequential access
        madvise(const_cast<char*>(data), size, MADV_SEQUENTIAL | MADV_WILLNEED);
        return true;
    }
};

// Output buffer for efficient writes
class OutputBuffer {
public:
    explicit OutputBuffer(std::ostream& os, size_t bufSize = 1024 * 1024)
        : out_(os), buffer_(bufSize), pos_(0) {}

    ~OutputBuffer() { flush(); }

    void append(std::string_view sv) {
        if (pos_ + sv.size() > buffer_.size()) {
            flush();
        }
        if (sv.size() > buffer_.size()) {
            out_.write(sv.data(), sv.size());
        } else {
            memcpy(buffer_.data() + pos_, sv.data(), sv.size());
            pos_ += sv.size();
        }
    }

    void append(char c) {
        if (pos_ >= buffer_.size()) {
            flush();
        }
        buffer_[pos_++] = c;
    }

    void flush() {
        if (pos_ > 0) {
            out_.write(buffer_.data(), pos_);
            pos_ = 0;
        }
    }

private:
    std::ostream& out_;
    std::vector<char> buffer_;
    size_t pos_;
};

// SIMD-accelerated newline search
inline const char* findNewline(const char* start, const char* end) {
#ifdef __AVX2__
    const __m256i newline = _mm256_set1_epi8('\n');
    while (start + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(start));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, newline);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask != 0) {
            return start + __builtin_ctz(mask);
        }
        start += 32;
    }
#endif
#ifdef __SSE2__
    const __m128i newline = _mm_set1_epi8('\n');
    while (start + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(start));
        __m128i cmp = _mm_cmpeq_epi8(chunk, newline);
        int mask = _mm_movemask_epi8(cmp);
        if (mask != 0) {
            return start + __builtin_ctz(mask);
        }
        start += 16;
    }
#endif
    // Fallback
    const char* p = static_cast<const char*>(memchr(start, '\n', end - start));
    return p ? p : end;
}

// Find next tab or end
inline const char* findTab(const char* start, const char* end) {
    const char* p = static_cast<const char*>(memchr(start, '\t', end - start));
    return p ? p : end;
}

// Find a specific INFO key and return its value
// Returns empty string_view if not found
inline std::string_view findInfoValue(std::string_view info, std::string_view key) {
    size_t pos = 0;
    while (pos < info.size()) {
        // Find the end of this entry (';' or end of string)
        size_t semi = info.find(';', pos);
        if (semi == std::string_view::npos) semi = info.size();

        std::string_view entry = info.substr(pos, semi - pos);

        // Check if this entry matches the key
        size_t eq = entry.find('=');
        std::string_view entryKey = (eq != std::string_view::npos) ? entry.substr(0, eq) : entry;

        if (entryKey == key) {
            if (eq != std::string_view::npos) {
                return entry.substr(eq + 1);
            } else {
                // Flag field - return empty but indicate found
                return std::string_view("", 0);
            }
        }

        pos = semi + 1;
    }
    return std::string_view();  // Not found
}

} // anonymous namespace

bool parseInfoFieldsMmap(const char* filepath, std::ostream& out,
                         const std::vector<std::string>& info_fields, bool quiet) {
    MappedFile file;
    if (!file.open(filepath)) {
        if (!quiet) {
            std::cerr << "Error: Cannot open file: " << filepath << "\n";
        }
        return false;
    }

    if (file.size == 0) {
        return true;
    }

    const char* data = file.data;
    const char* end = data + file.size;
    const char* pos = data;

    OutputBuffer buf(out);

    // Write header
    if (!info_fields.empty()) {
        buf.append("CHROM\tPOS\tID\tREF\tALT");
        for (const auto& field : info_fields) {
            buf.append('\t');
            buf.append(field);
        }
        buf.append('\n');
    }

    while (pos < end) {
        const char* lineEnd = findNewline(pos, end);
        std::string_view line(pos, lineEnd - pos);
        pos = (lineEnd < end) ? lineEnd + 1 : end;

        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Parse first 8 columns
        const char* p = line.data();
        const char* lend = p + line.size();

        // CHROM (col 0)
        const char* tab = findTab(p, lend);
        std::string_view chrom(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // POS (col 1)
        tab = findTab(p, lend);
        std::string_view posField(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // ID (col 2)
        tab = findTab(p, lend);
        std::string_view id(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // REF (col 3)
        tab = findTab(p, lend);
        std::string_view ref(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // ALT (col 4)
        tab = findTab(p, lend);
        std::string_view alt(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // QUAL (col 5) - skip
        tab = findTab(p, lend);
        if (tab >= lend) continue;
        p = tab + 1;

        // FILTER (col 6) - skip
        tab = findTab(p, lend);
        if (tab >= lend) continue;
        p = tab + 1;

        // INFO (col 7)
        tab = findTab(p, lend);
        std::string_view info(p, tab - p);

        // Output: CHROM POS ID REF ALT + requested fields
        buf.append(chrom);
        buf.append('\t');
        buf.append(posField);
        buf.append('\t');
        buf.append(id);
        buf.append('\t');
        buf.append(ref);
        buf.append('\t');
        buf.append(alt);

        for (const auto& field : info_fields) {
            buf.append('\t');
            std::string_view value = findInfoValue(info, field);
            if (value.data() == nullptr) {
                // Not found
                buf.append('.');
            } else if (value.empty()) {
                // Flag (present but no value)
                buf.append('.');
            } else {
                buf.append(value);
            }
        }
        buf.append('\n');
    }

    return true;
}

#else
// Non-Unix fallback
bool parseInfoFieldsMmap(const char* filepath, std::ostream& out,
                         const std::vector<std::string>& info_fields, bool quiet) {
    (void)quiet;
    std::ifstream in(filepath);
    if (!in) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }
    return parseInfoFields(in, out, info_fields);
}
#endif

static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_info_parser", show_help))
        return 0;

    std::vector<std::string> info_fields;
    const char* inputFile = nullptr;
    bool quiet = false;

    static struct option long_options[] = {
        {"info", required_argument, nullptr, 'i'},
        {"input", required_argument, nullptr, 'I'},
        {"quiet", no_argument, nullptr, 'q'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    optind = 1;
    while ((opt = getopt_long(argc, argv, "i:I:qh", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'i': {
                std::string fields_str = optarg;
                std::stringstream ss(fields_str);
                std::string field;
                while (std::getline(ss, field, ',')) {
                    field.erase(0, field.find_first_not_of(" \t\n\r\f\v"));
                    field.erase(field.find_last_not_of(" \t\n\r\f\v") + 1);
                    if (!field.empty()) {
                        info_fields.push_back(field);
                    }
                }
                break;
            }
            case 'I':
                inputFile = optarg;
                break;
            case 'q':
                quiet = true;
                break;
            case 'h':
                printHelp();
                return 0;
            default:
                break;
        }
    }

    if (info_fields.empty()) {
        std::cerr << "Error: INFO fields not specified.\n"
                  << "Use --help for usage information.\n";
        return 1;
    }

    bool ok;
    if (inputFile) {
        ok = parseInfoFieldsMmap(inputFile, std::cout, info_fields, quiet);
    } else {
        ok = parseInfoFields(std::cin, std::cout, info_fields);
    }

    return ok ? 0 : 1;
}
