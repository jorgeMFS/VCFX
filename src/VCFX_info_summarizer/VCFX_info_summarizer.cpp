#include "VCFX_info_summarizer.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iomanip>
#include <map>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
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
    std::cout << "VCFX_info_summarizer\n"
              << "Usage: VCFX_info_summarizer [OPTIONS]\n\n"
              << "Options:\n"
              << "  -i, --info \"FIELD1,FIELD2\"   Specify the INFO fields to summarize (e.g., \"DP,AF\").\n"
              << "  -I, --input FILE             Input VCF file (default: stdin).\n"
              << "  -q, --quiet                  Suppress warnings.\n"
              << "  -h, --help                   Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Summarizes numeric fields in the INFO column of a VCF file by calculating\n"
              << "  statistics such as mean, median, and mode.\n\n"
              << "Examples:\n"
              << "  VCFX_info_summarizer --info \"DP,AF\" < input.vcf > summary_stats.tsv\n"
              << "  VCFX_info_summarizer -i \"DP,AF\" -I input.vcf > summary_stats.tsv\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char *argv[], std::vector<std::string> &info_fields) {
    bool foundAnyField = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // Check for --info or -i with the next argument
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
        }
        // Check for --info=<FIELDS>
        else if (arg.rfind("--info=", 0) == 0) {
            std::string fields_str = arg.substr(7);
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
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        }
        // else ignore unrecognized options
    }

    if (!foundAnyField) {
        std::cerr << "Error: INFO fields not specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return false;
    }
    return true;
}

// Function to calculate mean
double calculateMean(const std::vector<double> &data) {
    if (data.empty())
        return 0.0;
    double sum = 0.0;
    for (auto val : data) {
        sum += val;
    }
    return sum / static_cast<double>(data.size());
}

// Function to calculate median
double calculateMedian(std::vector<double> data) {
    if (data.empty())
        return 0.0;
    std::sort(data.begin(), data.end());
    size_t n = data.size();
    if (n % 2 == 0) {
        return (data[n / 2 - 1] + data[n / 2]) / 2.0;
    } else {
        return data[n / 2];
    }
}

// Function to calculate mode
double calculateMode(const std::vector<double> &data) {
    if (data.empty())
        return 0.0;
    std::unordered_map<double, int> frequency;
    int maxFreq = 0;
    double modeValue = data[0];

    for (auto val : data) {
        frequency[val]++;
        if (frequency[val] > maxFreq) {
            maxFreq = frequency[val];
            modeValue = val;
        }
    }
    return modeValue;
}

// Function to parse the INFO field and collect specified fields
bool summarizeInfoFields(std::istream &in, std::ostream &out, const std::vector<std::string> &info_fields) {
    std::string line;
    bool header_found = false;

    // Map to store vectors of values for each requested INFO field
    std::map<std::string, std::vector<double>> info_data;
    for (const auto &field : info_fields) {
        info_data[field]; // ensures key is created
    }

    // Performance: reuse containers across iterations
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty())
            continue;

        if (line[0] == '#') {
            // Check if it is #CHROM
            if (line.rfind("#CHROM", 0) == 0) {
                header_found = true;
            }
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        // parse columns
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping malformed VCF line: " << line << "\n";
            continue;
        }

        // Parse INFO field => store key->value
        std::unordered_map<std::string, std::string> info_map;
        {
            std::stringstream info_ss(fields[7]);
            std::string kv;
            while (std::getline(info_ss, kv, ';')) {
                if (kv.empty())
                    continue;
                size_t eq = kv.find('=');
                if (eq != std::string::npos) {
                    std::string key = kv.substr(0, eq);
                    std::string value = kv.substr(eq + 1);
                    info_map[key] = value;
                } else {
                    // A flag => store "1"
                    info_map[kv] = "1";
                }
            }
        }

        // For each requested field
        for (const auto &field : info_fields) {
            auto it = info_map.find(field);
            if (it == info_map.end()) {
                // not present
                continue;
            }
            // If present, possibly multiple comma-separated values
            std::stringstream valSS(it->second);
            std::string val;
            while (std::getline(valSS, val, ',')) {
                try {
                    double v = std::stod(val);
                    // Skip NaN / Inf
                    if (std::isnan(v) || std::isinf(v)) {
                        std::cerr << "Warning: Non-finite value for field " << field << " in line: " << line << "\n";
                        continue; // skip
                    }
                    info_data[field].push_back(v);
                } catch (...) {
                    std::cerr << "Warning: Non-numeric value for field " << field << " in line: " << line << "\n";
                }
            }
        }
    }

    // Print summary table
    out << "INFO_Field\tMean\tMedian\tMode\n";
    for (const auto &field : info_fields) {
        const auto &data = info_data.at(field);
        if (data.empty()) {
            // no numeric data => NA
            out << field << "\tNA\tNA\tNA\n";
            continue;
        }
        double mean = calculateMean(data);
        double median = calculateMedian(data);
        double mode = calculateMode(data);
        out << field << "\t" << std::fixed << std::setprecision(4) << mean << "\t" << median << "\t" << mode << "\n";
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
inline std::string_view findInfoValue(std::string_view info, std::string_view key) {
    size_t pos = 0;
    while (pos < info.size()) {
        size_t semi = info.find(';', pos);
        if (semi == std::string_view::npos) semi = info.size();

        std::string_view entry = info.substr(pos, semi - pos);
        size_t eq = entry.find('=');
        std::string_view entryKey = (eq != std::string_view::npos) ? entry.substr(0, eq) : entry;

        if (entryKey == key) {
            if (eq != std::string_view::npos) {
                return entry.substr(eq + 1);
            } else {
                // Flag field
                return std::string_view("1", 1);
            }
        }

        pos = semi + 1;
    }
    return std::string_view();  // Not found
}

// Fast double parsing from string_view
inline bool parseDouble(std::string_view sv, double& result) {
    if (sv.empty()) return false;

    // Quick validation and parsing
    const char* start = sv.data();
    const char* end = start + sv.size();
    char* parseEnd = nullptr;

    result = std::strtod(start, &parseEnd);

    // Check if we parsed the entire string (or close to it)
    if (parseEnd <= start) return false;

    // Check for NaN/Inf
    if (std::isnan(result) || std::isinf(result)) return false;

    return true;
}

} // anonymous namespace

bool summarizeInfoFieldsMmap(const char* filepath, std::ostream& out,
                              const std::vector<std::string>& info_fields, bool quiet) {
    MappedFile file;
    if (!file.open(filepath)) {
        if (!quiet) {
            std::cerr << "Error: Cannot open file: " << filepath << "\n";
        }
        return false;
    }

    if (file.size == 0) {
        // Empty file - still output header
        out << "INFO_Field\tMean\tMedian\tMode\n";
        for (const auto& field : info_fields) {
            out << field << "\tNA\tNA\tNA\n";
        }
        return true;
    }

    const char* data = file.data;
    const char* end = data + file.size;
    const char* pos = data;

    // Map to store vectors of values for each requested INFO field
    std::map<std::string, std::vector<double>> info_data;
    for (const auto& field : info_fields) {
        info_data[field].reserve(10000);  // Pre-allocate for large files
    }

    bool header_found = false;

    while (pos < end) {
        const char* lineEnd = findNewline(pos, end);
        std::string_view line(pos, lineEnd - pos);
        pos = (lineEnd < end) ? lineEnd + 1 : end;

        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                header_found = true;
            }
            continue;
        }

        if (!header_found) {
            if (!quiet) {
                std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            }
            return false;
        }

        // Skip to INFO field (column 7)
        const char* p = line.data();
        const char* lend = p + line.size();

        // Skip first 7 columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER)
        for (int col = 0; col < 7; ++col) {
            const char* tab = findTab(p, lend);
            if (tab >= lend) {
                p = lend;
                break;
            }
            p = tab + 1;
        }

        if (p >= lend) continue;

        // INFO field
        const char* tab = findTab(p, lend);
        std::string_view info(p, tab - p);

        // For each requested field
        for (const auto& field : info_fields) {
            std::string_view value = findInfoValue(info, field);
            if (value.data() == nullptr) continue;

            // Parse comma-separated values
            size_t vpos = 0;
            while (vpos < value.size()) {
                size_t comma = value.find(',', vpos);
                if (comma == std::string_view::npos) comma = value.size();

                std::string_view val = value.substr(vpos, comma - vpos);
                double d;
                if (parseDouble(val, d)) {
                    info_data[field].push_back(d);
                }

                vpos = comma + 1;
            }
        }
    }

    // Print summary table
    out << "INFO_Field\tMean\tMedian\tMode\n";
    for (const auto& field : info_fields) {
        const auto& fieldData = info_data.at(field);
        if (fieldData.empty()) {
            out << field << "\tNA\tNA\tNA\n";
            continue;
        }
        double mean = calculateMean(fieldData);
        double median = calculateMedian(fieldData);
        double mode = calculateMode(fieldData);
        out << field << "\t" << std::fixed << std::setprecision(4)
            << mean << "\t" << median << "\t" << mode << "\n";
    }

    return true;
}

#else
// Non-Unix fallback
bool summarizeInfoFieldsMmap(const char* filepath, std::ostream& out,
                              const std::vector<std::string>& info_fields, bool quiet) {
    (void)quiet;
    std::ifstream in(filepath);
    if (!in) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }
    return summarizeInfoFields(in, out, info_fields);
}
#endif

static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_info_summarizer", show_help))
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

    bool success;
    if (inputFile) {
        success = summarizeInfoFieldsMmap(inputFile, std::cout, info_fields, quiet);
    } else {
        success = summarizeInfoFields(std::cin, std::cout, info_fields);
    }

    return success ? 0 : 1;
}
