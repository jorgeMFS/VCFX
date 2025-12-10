#include "VCFX_info_aggregator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string_view>
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

// ----------------------------------------------------------------------
// Displays help
// ----------------------------------------------------------------------
void VCFXInfoAggregator::displayHelp() {
    std::cout << "VCFX_info_aggregator: Aggregate numeric INFO field values from a VCF.\n\n"
              << "Usage:\n"
              << "  VCFX_info_aggregator --aggregate-info \"DP,AF,...\" < input.vcf > output.vcf\n"
              << "  VCFX_info_aggregator -a \"DP,AF,...\" -i input.vcf > output.vcf\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin, prints it unmodified, and at the end, appends a\n"
              << "  summary section of the form:\n"
              << "    #AGGREGATION_SUMMARY\n"
              << "    DP: Sum=..., Average=...\n"
              << "    AF: Sum=..., Average=...\n"
              << "  The VCF portion remains fully valid. The final lines start with '#' so most\n"
              << "  VCF parsers will ignore them.\n\n"
              << "Options:\n"
              << "  -h, --help                     Print this help message.\n"
              << "  -a, --aggregate-info <fields>  Comma-separated list of INFO fields to aggregate.\n"
              << "  -i, --input FILE               Input VCF file (default: stdin).\n"
              << "  -q, --quiet                    Suppress warnings.\n\n"
              << "Example:\n"
              << "  VCFX_info_aggregator --aggregate-info \"DP,AF\" < input.vcf > aggregated.vcf\n"
              << "  VCFX_info_aggregator -a \"DP,AF\" -i input.vcf > aggregated.vcf\n";
}

// ----------------------------------------------------------------------
// parse command line, call aggregator
// ----------------------------------------------------------------------
int VCFXInfoAggregator::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    std::string infoFieldsStr;
    const char* inputFile = nullptr;
    bool quiet = false;

    static struct option longOptions[] = {
        {"help", no_argument, 0, 'h'},
        {"aggregate-info", required_argument, 0, 'a'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    optind = 1;
    while ((opt = getopt_long(argc, argv, "ha:i:q", longOptions, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'a':
            infoFieldsStr = optarg;
            break;
        case 'i':
            inputFile = optarg;
            break;
        case 'q':
            quiet = true;
            break;
        default:
            showHelp = true;
            break;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }
    if (infoFieldsStr.empty()) {
        std::cerr << "Error: Must specify --aggregate-info with at least one field.\n";
        return 1;
    }

    // parse comma separated fields
    std::vector<std::string> infoFields;
    {
        std::stringstream ss(infoFieldsStr);
        std::string f;
        while (std::getline(ss, f, ',')) {
            // trim
            while (!f.empty() && isspace((unsigned char)f.back()))
                f.pop_back();
            while (!f.empty() && isspace((unsigned char)f.front()))
                f.erase(f.begin());
            if (!f.empty()) {
                infoFields.push_back(f);
            }
        }
    }
    if (infoFields.empty()) {
        std::cerr << "Error: no valid fields in --aggregate-info\n";
        return 1;
    }

    if (inputFile) {
        return aggregateInfoMmap(inputFile, std::cout, infoFields, quiet) ? 0 : 1;
    } else {
        aggregateInfo(std::cin, std::cout, infoFields);
        return 0;
    }
}

// ----------------------------------------------------------------------
// aggregator function (stdin fallback)
// ----------------------------------------------------------------------
void VCFXInfoAggregator::aggregateInfo(std::istream &in, std::ostream &out,
                                       const std::vector<std::string> &infoFields) {
    // We'll store each field's collected numeric values in a vector
    std::map<std::string, std::vector<double>> collected;
    for (auto &fld : infoFields) {
        collected[fld] = {};
    }

    bool foundChromHeader = false;

    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);
    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            out << line << "\n";
            continue;
        }

        if (line[0] == '#') {
            // pass it through unmodified
            out << line << "\n";
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Error: encountered data line before #CHROM header.\n";
            break; // break or continue, but we'll break
        }

        // parse columns
        // minimal vcf => CHROM POS ID REF ALT QUAL FILTER INFO ...
        // We only need the 8th col => INFO
        vcfx::split_tabs(line, fields);
        // pass line unmodified
        out << line << "\n";

        if (fields.size() < 8) {
            // skip aggregator for this line
            continue;
        }
        std::string &info = fields[7];

        // parse info col => find each "key=val" or "key" flag
        // we only look for key=val with numeric val if the key is in infoFields
        // We'll do a naive parse: split by ';', then split each part at '='
        {
            std::stringstream infoSS(info);
            std::string item;
            while (std::getline(infoSS, item, ';')) {
                if (item.empty())
                    continue;
                // find '='
                size_t eqPos = item.find('=');
                if (eqPos == std::string::npos) {
                    // no '=' => skip
                    continue;
                }
                std::string key = item.substr(0, eqPos);
                std::string val = item.substr(eqPos + 1);
                if (key.empty() || val.empty())
                    continue;
                // trim spaces
                while (!key.empty() && isspace((unsigned char)key.back()))
                    key.pop_back();
                while (!key.empty() && isspace((unsigned char)key.front()))
                    key.erase(key.begin());
                while (!val.empty() && isspace((unsigned char)val.back()))
                    val.pop_back();
                while (!val.empty() && isspace((unsigned char)val.front()))
                    val.erase(val.begin());

                auto it = collected.find(key);
                if (it == collected.end()) {
                    // not one of requested fields
                    continue;
                }
                // parse numeric, skip if NaN or Inf
                try {
                    double d = std::stod(val);
                    // check for finite
                    if (std::isfinite(d)) {
                        it->second.push_back(d);
                    }
                    // else skip
                } catch (...) {
                    // skip
                }
            }
        }
    }

    // now print aggregator summary
    out << "#AGGREGATION_SUMMARY\n";
    for (auto &kv : collected) {
        const std::string &field = kv.first;
        const auto &vals = kv.second;
        double sum = 0.0;
        for (auto &v : vals) {
            sum += v;
        }
        double mean = 0.0;
        if (!vals.empty()) {
            mean = sum / static_cast<double>(vals.size());
        }
        out << field << ": Sum=" << sum << ", Average=" << mean << "\n";
    }
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
inline std::string_view findInfoValue(std::string_view info, std::string_view key) {
    size_t pos = 0;
    while (pos < info.size()) {
        size_t semi = info.find(';', pos);
        if (semi == std::string_view::npos) semi = info.size();

        std::string_view entry = info.substr(pos, semi - pos);
        size_t eq = entry.find('=');
        if (eq != std::string_view::npos) {
            std::string_view entryKey = entry.substr(0, eq);
            if (entryKey == key) {
                return entry.substr(eq + 1);
            }
        }

        pos = semi + 1;
    }
    return std::string_view();  // Not found
}

// Fast double parsing from string_view
inline bool parseDouble(std::string_view sv, double& result) {
    if (sv.empty()) return false;

    const char* start = sv.data();
    char* parseEnd = nullptr;

    result = std::strtod(start, &parseEnd);

    if (parseEnd <= start) return false;
    if (!std::isfinite(result)) return false;

    return true;
}

} // anonymous namespace

bool VCFXInfoAggregator::aggregateInfoMmap(const char* filepath, std::ostream& out,
                                            const std::vector<std::string>& infoFields, bool quiet) {
    MappedFile file;
    if (!file.open(filepath)) {
        if (!quiet) {
            std::cerr << "Error: Cannot open file: " << filepath << "\n";
        }
        return false;
    }

    if (file.size == 0) {
        out << "#AGGREGATION_SUMMARY\n";
        for (const auto& field : infoFields) {
            out << field << ": Sum=0, Average=0\n";
        }
        return true;
    }

    const char* data = file.data;
    const char* end = data + file.size;
    const char* pos = data;

    // We'll store each field's collected numeric values
    std::map<std::string, std::vector<double>> collected;
    for (const auto& fld : infoFields) {
        collected[fld].reserve(10000);
    }

    OutputBuffer buf(out);
    bool foundChromHeader = false;

    while (pos < end) {
        const char* lineEnd = findNewline(pos, end);
        std::string_view line(pos, lineEnd - pos);

        // Output line unchanged
        buf.append(line);
        buf.append('\n');

        pos = (lineEnd < end) ? lineEnd + 1 : end;

        if (line.empty()) continue;

        if (line[0] == '#') {
            if (!foundChromHeader && line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            if (!quiet) {
                std::cerr << "Error: encountered data line before #CHROM header.\n";
            }
            return false;
        }

        // Skip to INFO field (column 7)
        const char* p = line.data();
        const char* lend = p + line.size();

        for (int col = 0; col < 7; ++col) {
            const char* tab = findTab(p, lend);
            if (tab >= lend) {
                p = lend;
                break;
            }
            p = tab + 1;
        }

        if (p >= lend) continue;

        const char* tab = findTab(p, lend);
        std::string_view info(p, tab - p);

        // For each requested field
        for (const auto& field : infoFields) {
            std::string_view value = findInfoValue(info, field);
            if (value.data() == nullptr) continue;

            double d;
            if (parseDouble(value, d)) {
                collected[field].push_back(d);
            }
        }
    }

    // Flush the passthrough output
    buf.flush();

    // Now print aggregator summary
    out << "#AGGREGATION_SUMMARY\n";
    for (const auto& kv : collected) {
        const std::string& field = kv.first;
        const auto& vals = kv.second;
        double sum = 0.0;
        for (auto v : vals) {
            sum += v;
        }
        double mean = 0.0;
        if (!vals.empty()) {
            mean = sum / static_cast<double>(vals.size());
        }
        out << field << ": Sum=" << sum << ", Average=" << mean << "\n";
    }

    return true;
}

#else
// Non-Unix fallback
bool VCFXInfoAggregator::aggregateInfoMmap(const char* filepath, std::ostream& out,
                                            const std::vector<std::string>& infoFields, bool quiet) {
    (void)quiet;
    std::ifstream in(filepath);
    if (!in) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }
    aggregateInfo(in, out, infoFields);
    return true;
}
#endif

static void show_help() {
    VCFXInfoAggregator obj;
    char arg0[] = "VCFX_info_aggregator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_info_aggregator", show_help))
        return 0;
    VCFXInfoAggregator app;
    return app.run(argc, argv);
}
