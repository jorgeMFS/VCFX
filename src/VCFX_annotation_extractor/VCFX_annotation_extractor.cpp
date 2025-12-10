#include "VCFX_annotation_extractor.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
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

// --------------------------------------------------------------
// A small struct to store command-line options
// --------------------------------------------------------------
struct AnnotationOptions {
    std::vector<std::string> annotations; // e.g. ["ANN", "Gene"]
    const char* inputFile = nullptr;
    bool quiet = false;
};

// --------------------------------------------------------------
// Utility: split a string by a delimiter into a vector
// --------------------------------------------------------------
static std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

// --------------------------------------------------------------
// Utility: parse the INFO field into a map key->value
//   e.g. "ANN=xxx;Gene=YYY;DP=100" => {ANN:xxx, Gene:YYY, DP:100}
// --------------------------------------------------------------
static std::unordered_map<std::string, std::string> parseInfoToMap(const std::string &info) {
    std::unordered_map<std::string, std::string> infoMap;
    // split by ';'
    auto fields = split(info, ';');
    for (auto &f : fields) {
        if (f.empty()) {
            continue;
        }
        // e.g. f="ANN=..."
        auto eqPos = f.find('=');
        if (eqPos == std::string::npos) {
            // key without value? e.g. "SOMATIC"
            // You could store it as {SOMATIC: ""} if you want
            infoMap[f] = "";
            continue;
        }
        std::string key = f.substr(0, eqPos);
        std::string val = f.substr(eqPos + 1);
        infoMap[key] = val;
    }
    return infoMap;
}

// --------------------------------------------------------------
// Show usage/help
// --------------------------------------------------------------
static void printHelp() {
    std::cout << "VCFX_annotation_extractor: Extract variant annotations from a VCF file.\n\n"
              << "Usage:\n"
              << "  VCFX_annotation_extractor --annotation-extract \"ANN,Gene\" < input.vcf > out.tsv\n"
              << "  VCFX_annotation_extractor -a \"ANN,Gene\" -i input.vcf > out.tsv\n\n"
              << "Options:\n"
              << "  -a, --annotation-extract   Comma-separated list of annotations to extract (e.g., ANN,Gene)\n"
              << "  -i, --input FILE           Input VCF file (default: stdin)\n"
              << "  -q, --quiet                Suppress warnings\n"
              << "  -h, --help                 Display this help message and exit\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin and prints a tab-delimited output. For multi-ALT\n"
              << "  lines, each ALT allele is printed on its own line. If an annotation field (like\n"
              << "  'ANN=') has multiple comma-separated sub-entries, we attempt to align them with\n"
              << "  the ALT alleles in order.\n\n"
              << "Example:\n"
              << "  VCFX_annotation_extractor --annotation-extract \"ANN,Gene\" < input.vcf > out.tsv\n"
              << "  VCFX_annotation_extractor -a \"ANN,Gene\" -i input.vcf > out.tsv\n";
}

// --------------------------------------------------------------
// parseArguments: fill in AnnotationOptions
// --------------------------------------------------------------
static bool parseArguments(int argc, char *argv[], AnnotationOptions &opts) {
    bool showHelp = false;

    static struct option long_options[] = {
        {"annotation-extract", required_argument, 0, 'a'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    optind = 1;
    while (true) {
        int optIdx = 0;
        int c = getopt_long(argc, argv, "a:i:qh", long_options, &optIdx);
        if (c == -1)
            break;
        switch (c) {
        case 'a': {
            // comma-separated annotation names
            auto items = split(optarg, ',');
            for (auto &it : items) {
                // trim spaces
                while (!it.empty() && (it.front() == ' ' || it.front() == '\t')) {
                    it.erase(it.begin());
                }
                while (!it.empty() && (it.back() == ' ' || it.back() == '\t')) {
                    it.pop_back();
                }
                opts.annotations.push_back(it);
            }
        } break;
        case 'i':
            opts.inputFile = optarg;
            break;
        case 'q':
            opts.quiet = true;
            break;
        case 'h':
        default:
            showHelp = true;
            break;
        }
    }

    if (showHelp) {
        printHelp();
        // Return false to indicate we should exit.
        return false;
    }
    // If no annotations, also show help
    if (opts.annotations.empty()) {
        printHelp();
        return false;
    }
    return true;
}

// --------------------------------------------------------------
// main extraction logic (stdin fallback)
// --------------------------------------------------------------
static void processVCF(std::istream &in, const AnnotationOptions &opts) {
    std::string line;
    bool foundChromHeader = false;

    // Print a single header row for the TSV:
    std::cout << "CHROM\tPOS\tID\tREF\tALT";
    for (auto &annName : opts.annotations) {
        std::cout << "\t" << annName;
    }
    std::cout << "\n";

    // Declare fields vector outside loop for reuse
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // If header line
        if (line[0] == '#') {
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
            }
            continue;
        }

        // We expect #CHROM line before data
        if (!foundChromHeader) {
            if (!opts.quiet) {
                std::cerr << "Warning: Data encountered before #CHROM header: skipping\n";
            }
            continue;
        }

        // Split the line by tabs
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            if (!opts.quiet) {
                std::cerr << "Warning: Invalid VCF line (fewer than 8 fields): " << line << "\n";
            }
            continue;
        }

        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        const std::string &altStr = fields[4];
        const std::string &info = fields[7];

        // Split ALT on comma for multiple alt alleles
        std::vector<std::string> alts = split(altStr, ',');

        // Parse INFO into a map: KEY -> value
        auto infoMap = parseInfoToMap(info);

        // For each annotation requested, retrieve the value or "NA"
        std::unordered_map<std::string, std::string> rawAnnValues;
        for (auto &annName : opts.annotations) {
            auto it = infoMap.find(annName);
            if (it == infoMap.end()) {
                rawAnnValues[annName] = "NA";
            } else {
                rawAnnValues[annName] = it->second;
            }
        }

        // For ANN field, split on comma for multi-allele support
        std::unordered_map<std::string, std::vector<std::string>> splittedAnnValues;
        for (auto &annName : opts.annotations) {
            std::string val = rawAnnValues[annName];
            if (val == "NA") {
                splittedAnnValues[annName] = {};
                continue;
            }
            if (annName == "ANN") {
                auto subVals = split(val, ',');
                splittedAnnValues[annName] = subVals;
            }
        }

        // Now let's produce lines
        for (size_t altIndex = 0; altIndex < alts.size(); ++altIndex) {
            const std::string &thisAlt = alts[altIndex];

            std::ostringstream outLine;
            outLine << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << thisAlt;

            for (auto &annName : opts.annotations) {
                outLine << "\t";
                auto it = splittedAnnValues.find(annName);
                if (it != splittedAnnValues.end() && !it->second.empty()) {
                    if (altIndex < it->second.size()) {
                        outLine << it->second[altIndex];
                    } else {
                        outLine << "NA";
                    }
                } else {
                    auto rawIt = rawAnnValues.find(annName);
                    if (rawIt != rawAnnValues.end()) {
                        outLine << rawIt->second;
                    } else {
                        outLine << "NA";
                    }
                }
            }
            outLine << "\n";
            std::cout << outLine.str();
        }
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
    const char* p = static_cast<const char*>(memchr(start, '\n', end - start));
    return p ? p : end;
}

// Find next tab or end
inline const char* findTab(const char* start, const char* end) {
    const char* p = static_cast<const char*>(memchr(start, '\t', end - start));
    return p ? p : end;
}

// Find next comma or end
inline const char* findComma(const char* start, const char* end) {
    const char* p = static_cast<const char*>(memchr(start, ',', end - start));
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
                return std::string_view("", 0);
            }
        }

        pos = semi + 1;
    }
    return std::string_view();  // Not found (nullptr)
}

// Get nth comma-separated value
inline std::string_view getNthCommaSeparated(std::string_view sv, size_t n) {
    size_t pos = 0;
    size_t idx = 0;
    while (pos < sv.size()) {
        size_t comma = sv.find(',', pos);
        if (comma == std::string_view::npos) comma = sv.size();

        if (idx == n) {
            return sv.substr(pos, comma - pos);
        }

        pos = comma + 1;
        ++idx;
    }
    return std::string_view();
}

} // anonymous namespace

static bool processVCFMmap(const char* filepath, std::ostream& out,
                           const AnnotationOptions& opts) {
    MappedFile file;
    if (!file.open(filepath)) {
        if (!opts.quiet) {
            std::cerr << "Error: Cannot open file: " << filepath << "\n";
        }
        return false;
    }

    if (file.size == 0) {
        // Just output header
        out << "CHROM\tPOS\tID\tREF\tALT";
        for (const auto& annName : opts.annotations) {
            out << "\t" << annName;
        }
        out << "\n";
        return true;
    }

    const char* data = file.data;
    const char* end = data + file.size;
    const char* pos = data;

    OutputBuffer buf(out);

    // Print header
    buf.append("CHROM\tPOS\tID\tREF\tALT");
    for (const auto& annName : opts.annotations) {
        buf.append('\t');
        buf.append(annName);
    }
    buf.append('\n');

    bool foundChromHeader = false;

    while (pos < end) {
        const char* lineEnd = findNewline(pos, end);
        std::string_view line(pos, lineEnd - pos);
        pos = (lineEnd < end) ? lineEnd + 1 : end;

        if (line.empty()) continue;

        if (line[0] == '#') {
            if (!foundChromHeader && line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            if (!opts.quiet) {
                std::cerr << "Warning: Data encountered before #CHROM header: skipping\n";
            }
            continue;
        }

        // Parse first 8 columns
        const char* p = line.data();
        const char* lend = p + line.size();

        // CHROM
        const char* tab = findTab(p, lend);
        std::string_view chrom(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // POS
        tab = findTab(p, lend);
        std::string_view posField(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // ID
        tab = findTab(p, lend);
        std::string_view id(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // REF
        tab = findTab(p, lend);
        std::string_view ref(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // ALT
        tab = findTab(p, lend);
        std::string_view altStr(p, tab - p);
        if (tab >= lend) continue;
        p = tab + 1;

        // QUAL - skip
        tab = findTab(p, lend);
        if (tab >= lend) continue;
        p = tab + 1;

        // FILTER - skip
        tab = findTab(p, lend);
        if (tab >= lend) continue;
        p = tab + 1;

        // INFO
        tab = findTab(p, lend);
        std::string_view info(p, tab - p);

        // Count ALT alleles
        size_t numAlts = 1;
        for (char c : altStr) {
            if (c == ',') ++numAlts;
        }

        // For each ALT allele, output a line
        size_t altPos = 0;
        for (size_t altIdx = 0; altIdx < numAlts; ++altIdx) {
            size_t comma = altStr.find(',', altPos);
            if (comma == std::string_view::npos) comma = altStr.size();
            std::string_view thisAlt = altStr.substr(altPos, comma - altPos);
            altPos = comma + 1;

            buf.append(chrom);
            buf.append('\t');
            buf.append(posField);
            buf.append('\t');
            buf.append(id);
            buf.append('\t');
            buf.append(ref);
            buf.append('\t');
            buf.append(thisAlt);

            for (const auto& annName : opts.annotations) {
                buf.append('\t');
                std::string_view value = findInfoValue(info, annName);
                if (value.data() == nullptr) {
                    buf.append("NA");
                } else if (annName == "ANN" && numAlts > 1) {
                    // ANN is multi-allele, get the nth entry
                    std::string_view subVal = getNthCommaSeparated(value, altIdx);
                    if (subVal.data() == nullptr || subVal.empty()) {
                        buf.append("NA");
                    } else {
                        buf.append(subVal);
                    }
                } else {
                    if (value.empty()) {
                        buf.append("NA");
                    } else {
                        buf.append(value);
                    }
                }
            }
            buf.append('\n');
        }
    }

    return true;
}

#else
// Non-Unix fallback
static bool processVCFMmap(const char* filepath, std::ostream& out,
                           const AnnotationOptions& opts) {
    std::ifstream in(filepath);
    if (!in) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }
    processVCF(in, opts);
    return true;
}
#endif

// --------------------------------------------------------------
// main()
// --------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_annotation_extractor", show_help))
        return 0;
    AnnotationOptions opts;
    if (!parseArguments(argc, argv, opts)) {
        return 1;
    }

    if (opts.inputFile) {
        return processVCFMmap(opts.inputFile, std::cout, opts) ? 0 : 1;
    } else {
        processVCF(std::cin, opts);
        return 0;
    }
}
