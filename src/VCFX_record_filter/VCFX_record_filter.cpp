#include "VCFX_record_filter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// ============================================================================
// SIMD includes for fast newline scanning
// ============================================================================
#if defined(__x86_64__) && defined(__AVX2__)
#include <immintrin.h>
#define USE_AVX2 1
#elif defined(__x86_64__) && defined(__SSE2__)
#include <emmintrin.h>
#define USE_SSE2 1
#endif

// ============================================================================
// SIMD-optimized newline finder
// ============================================================================
#if defined(USE_AVX2)
static inline const char* findNewline(const char* p, const char* end) {
    const __m256i nl = _mm256_set1_epi8('\n');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, nl);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 32;
    }
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#elif defined(USE_SSE2)
static inline const char* findNewline(const char* p, const char* end) {
    const __m128i nl = _mm_set1_epi8('\n');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, nl);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 16;
    }
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#else
static inline const char* findNewline(const char* p, const char* end) {
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#endif

// ============================================================================
// Fast string_view trimming (no allocation)
// ============================================================================
static inline std::string_view trimView(std::string_view sv) {
    while (!sv.empty() && (sv.front() == ' ' || sv.front() == '\t'))
        sv.remove_prefix(1);
    while (!sv.empty() && (sv.back() == ' ' || sv.back() == '\t'))
        sv.remove_suffix(1);
    return sv;
}

// ============================================================================
// Parse operator from string
// ============================================================================
static bool parseOperator(std::string_view opStr, FilterOp& op) {
    if (opStr == ">") { op = FilterOp::GT; return true; }
    if (opStr == ">=") { op = FilterOp::GE; return true; }
    if (opStr == "<") { op = FilterOp::LT; return true; }
    if (opStr == "<=") { op = FilterOp::LE; return true; }
    if (opStr == "==") { op = FilterOp::EQ; return true; }
    if (opStr == "!=") { op = FilterOp::NE; return true; }
    return false;
}

// ============================================================================
// Parse a single criterion string like "POS>100" or "FILTER==PASS"
// ============================================================================
static bool parseSingleCriterion(std::string_view token, FilterCriterion& crit) {
    // Find operator - check 2-char operators first
    static const std::pair<std::string_view, FilterOp> ops[] = {
        {">=", FilterOp::GE}, {"<=", FilterOp::LE},
        {"==", FilterOp::EQ}, {"!=", FilterOp::NE},
        {">", FilterOp::GT}, {"<", FilterOp::LT}
    };

    size_t opPos = std::string_view::npos;
    size_t opLen = 0;
    FilterOp foundOp;

    for (const auto& [opStr, opVal] : ops) {
        auto p = token.find(opStr);
        if (p != std::string_view::npos) {
            opPos = p;
            opLen = opStr.size();
            foundOp = opVal;
            break;
        }
    }

    if (opPos == std::string_view::npos) {
        std::cerr << "Error: no operator found in '" << token << "'.\n";
        return false;
    }

    // Extract field name and value
    std::string_view fieldName = trimView(token.substr(0, opPos));
    std::string_view valPart = trimView(token.substr(opPos + opLen));

    if (fieldName.empty()) {
        std::cerr << "Error: empty field name in '" << token << "'.\n";
        return false;
    }
    if (valPart.empty()) {
        std::cerr << "Error: no value in '" << token << "'.\n";
        return false;
    }

    crit.fieldName = std::string(fieldName);
    crit.op = foundOp;

    // Compile target field for fast dispatch
    if (fieldName == "POS") {
        crit.target = TargetField::POS;
    } else if (fieldName == "QUAL") {
        crit.target = TargetField::QUAL;
    } else if (fieldName == "FILTER") {
        crit.target = TargetField::FILTER;
    } else {
        crit.target = TargetField::INFO_KEY;
    }

    // Try to parse as numeric using fast charconv
    double d;
    auto result = std::from_chars(valPart.data(), valPart.data() + valPart.size(), d);
    if (result.ec == std::errc{} && result.ptr == valPart.data() + valPart.size()) {
        crit.fieldType = FieldType::NUMERIC;
        crit.numericValue = d;
        crit.stringValue.clear();
    } else {
        crit.fieldType = FieldType::STRING;
        crit.stringValue = std::string(valPart);
        crit.numericValue = 0.0;
    }

    return true;
}

// ============================================================================
// VCFXRecordFilter::parseCriteria - Parse semicolon-separated criteria
// ============================================================================
bool VCFXRecordFilter::parseCriteria(const std::string& criteriaStr,
                                      std::vector<FilterCriterion>& criteria) {
    criteria.clear();
    std::string_view sv(criteriaStr);

    size_t start = 0;
    while (start < sv.size()) {
        size_t end = sv.find(';', start);
        if (end == std::string_view::npos) end = sv.size();

        std::string_view token = trimView(sv.substr(start, end - start));
        if (!token.empty()) {
            FilterCriterion c;
            if (!parseSingleCriterion(token, c)) {
                return false;
            }
            criteria.push_back(std::move(c));
        }
        start = end + 1;
    }

    if (criteria.empty()) {
        std::cerr << "Error: no valid criteria in '" << criteriaStr << "'.\n";
        return false;
    }
    return true;
}

// ============================================================================
// Fast field extraction - extract Nth tab-delimited field (zero-copy)
// ============================================================================
std::string_view VCFXRecordFilter::extractField(std::string_view line, int fieldIndex) {
    const char* ptr = line.data();
    const char* end = ptr + line.size();
    int currentField = 0;

    // Fast scan to find field start
    while (currentField < fieldIndex && ptr < end) {
        if (*ptr == '\t') currentField++;
        ptr++;
    }

    if (currentField < fieldIndex) {
        return {};  // Not enough fields
    }

    // Find field end
    const char* fieldStart = ptr;
    while (ptr < end && *ptr != '\t') {
        ptr++;
    }

    return std::string_view(fieldStart, ptr - fieldStart);
}

// ============================================================================
// Fast INFO value extraction - find key=value in INFO field (zero-copy)
// ============================================================================
bool VCFXRecordFilter::extractInfoValue(std::string_view info, std::string_view key,
                                         std::string_view& valueOut) {
    if (info.empty() || info == ".") return false;

    // Scan INFO field for key
    size_t pos = 0;
    while (pos < info.size()) {
        // Find end of this token (semicolon or end)
        size_t tokenEnd = info.find(';', pos);
        if (tokenEnd == std::string_view::npos) tokenEnd = info.size();

        std::string_view token = info.substr(pos, tokenEnd - pos);

        // Check for key=value or flag
        size_t eq = token.find('=');
        if (eq != std::string_view::npos) {
            std::string_view k = token.substr(0, eq);
            if (k == key) {
                valueOut = token.substr(eq + 1);
                return true;
            }
        } else {
            // Flag (no value) - return the flag itself as value
            if (token == key) {
                valueOut = token;
                return true;
            }
        }

        pos = tokenEnd + 1;
    }

    return false;
}

// ============================================================================
// Fast numeric parsing using std::from_chars (no exceptions, no allocations)
// ============================================================================
bool VCFXRecordFilter::parseDouble(std::string_view sv, double& out) {
    if (sv.empty()) return false;
    auto result = std::from_chars(sv.data(), sv.data() + sv.size(), out);
    return result.ec == std::errc{} && result.ptr == sv.data() + sv.size();
}

bool VCFXRecordFilter::parseInt(std::string_view sv, int64_t& out) {
    if (sv.empty()) return false;
    auto result = std::from_chars(sv.data(), sv.data() + sv.size(), out);
    return result.ec == std::errc{} && result.ptr == sv.data() + sv.size();
}

// ============================================================================
// Comparison functions (branchless-friendly)
// ============================================================================
bool VCFXRecordFilter::compareDouble(double x, FilterOp op, double y) {
    switch (op) {
        case FilterOp::GT: return x > y;
        case FilterOp::GE: return x >= y;
        case FilterOp::LT: return x < y;
        case FilterOp::LE: return x <= y;
        case FilterOp::EQ: return x == y;
        case FilterOp::NE: return x != y;
    }
    return false;
}

bool VCFXRecordFilter::compareString(std::string_view s, FilterOp op, std::string_view t) {
    switch (op) {
        case FilterOp::EQ: return s == t;
        case FilterOp::NE: return s != t;
        default: return false;  // < and > not supported for strings
    }
}

// ============================================================================
// Evaluate a single criterion against a line (hot path - must be fast)
// ============================================================================
bool VCFXRecordFilter::evaluateCriterion(std::string_view line,
                                          const FilterCriterion& c) const {
    // Use compiled target for fast dispatch
    switch (c.target) {
        case TargetField::POS: {
            std::string_view posStr = extractField(line, 1);
            if (posStr.empty()) return false;
            double pos;
            if (!parseDouble(posStr, pos)) return false;
            return compareDouble(pos, c.op, c.numericValue);
        }

        case TargetField::QUAL: {
            std::string_view qualStr = extractField(line, 5);
            if (qualStr.empty() || qualStr == ".") {
                return compareDouble(0.0, c.op, c.numericValue);
            }
            double qual;
            if (!parseDouble(qualStr, qual)) return false;
            return compareDouble(qual, c.op, c.numericValue);
        }

        case TargetField::FILTER: {
            std::string_view filterStr = extractField(line, 6);
            if (c.fieldType == FieldType::NUMERIC) return false;
            return compareString(filterStr, c.op, c.stringValue);
        }

        case TargetField::INFO_KEY: {
            std::string_view info = extractField(line, 7);
            std::string_view value;
            if (!extractInfoValue(info, c.fieldName, value)) {
                return false;  // Key not found
            }

            if (c.fieldType == FieldType::NUMERIC) {
                double num;
                if (!parseDouble(value, num)) return false;
                return compareDouble(num, c.op, c.numericValue);
            } else {
                return compareString(value, c.op, c.stringValue);
            }
        }
    }
    return false;
}

// ============================================================================
// Evaluate all criteria against a line
// ============================================================================
bool VCFXRecordFilter::evaluateLine(std::string_view line) const {
    if (useAndLogic_) {
        // AND logic: all must pass
        for (const auto& c : criteria_) {
            if (!evaluateCriterion(line, c)) {
                return false;
            }
        }
        return true;
    } else {
        // OR logic: any must pass
        for (const auto& c : criteria_) {
            if (evaluateCriterion(line, c)) {
                return true;
            }
        }
        return false;
    }
}

// ============================================================================
// Memory-mapped file processing (fast path)
// ============================================================================
bool VCFXRecordFilter::processFileMmap(const char* filepath) {
    int fd = open(filepath, O_RDONLY);
    if (fd < 0) {
        std::cerr << "Error: cannot open file '" << filepath << "'\n";
        return false;
    }

    struct stat st;
    if (fstat(fd, &st) < 0) {
        close(fd);
        std::cerr << "Error: cannot stat file\n";
        return false;
    }

    size_t fileSize = st.st_size;
    if (fileSize == 0) {
        close(fd);
        return true;  // Empty file
    }

    void* mapped = mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        std::cerr << "Error: cannot mmap file\n";
        return false;
    }

    // Advise kernel for sequential access
    madvise(mapped, fileSize, MADV_SEQUENTIAL);

    const char* data = static_cast<const char*>(mapped);
    const char* ptr = data;
    const char* end = data + fileSize;

    // Output buffer for efficient writes
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);  // 1MB buffer

    bool foundChrom = false;

    while (ptr < end) {
        // Find end of line using SIMD
        const char* lineEnd = findNewline(ptr, end);
        if (!lineEnd) lineEnd = end;

        std::string_view line(ptr, lineEnd - ptr);

        // Handle CR-LF
        if (!line.empty() && line.back() == '\r') {
            line.remove_suffix(1);
        }

        // Process line
        if (line.empty()) {
            outputBuffer.push_back('\n');
        } else if (line[0] == '#') {
            // Header line - pass through
            outputBuffer.append(line.data(), line.size());
            outputBuffer.push_back('\n');
            if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChrom = true;
            }
        } else if (foundChrom) {
            // Data line - apply filter
            if (evaluateLine(line)) {
                outputBuffer.append(line.data(), line.size());
                outputBuffer.push_back('\n');
            }
        }

        // Flush buffer periodically
        if (outputBuffer.size() > 512 * 1024) {
            std::cout.write(outputBuffer.data(), outputBuffer.size());
            outputBuffer.clear();
        }

        ptr = lineEnd + 1;
    }

    // Flush remaining
    if (!outputBuffer.empty()) {
        std::cout.write(outputBuffer.data(), outputBuffer.size());
    }

    munmap(mapped, fileSize);
    close(fd);
    return true;
}

// ============================================================================
// Stdin processing (fallback for pipes)
// ============================================================================
void VCFXRecordFilter::processStdin() {
    std::string line;
    line.reserve(65536);

    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);

    bool foundChrom = false;

    while (std::getline(std::cin, line)) {
        // Handle CR-LF
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (line.empty()) {
            outputBuffer.push_back('\n');
            continue;
        }

        if (line[0] == '#') {
            outputBuffer.append(line);
            outputBuffer.push_back('\n');
            if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChrom = true;
            }
            continue;
        }

        if (!foundChrom) {
            std::cerr << "Warning: data line before #CHROM => skipping.\n";
            continue;
        }

        // Apply filter
        if (evaluateLine(line)) {
            outputBuffer.append(line);
            outputBuffer.push_back('\n');
        }

        // Flush buffer periodically
        if (outputBuffer.size() > 512 * 1024) {
            std::cout.write(outputBuffer.data(), outputBuffer.size());
            outputBuffer.clear();
        }
    }

    // Flush remaining
    if (!outputBuffer.empty()) {
        std::cout.write(outputBuffer.data(), outputBuffer.size());
    }
}

// ============================================================================
// Help display
// ============================================================================
void VCFXRecordFilter::displayHelp() {
    std::cout
        << "VCFX_record_filter: Filter VCF data lines by multiple criteria.\n\n"
           "Usage:\n"
           "  VCFX_record_filter [options] --filter \"CRITERIA\" [input.vcf]\n"
           "  VCFX_record_filter [options] --filter \"CRITERIA\" < input.vcf > output.vcf\n\n"
           "Options:\n"
           "  -f, --filter \"...\"   One or more criteria separated by semicolons, e.g.\n"
           "                        \"POS>10000; QUAL>=30; AF<0.05; FILTER==PASS\"\n"
           "                        Each criterion must use an operator among >,>=,<,<=,==,!=\n\n"
           "  -l, --logic and|or    'and' => a line must pass all criteria (default)\n"
           "                        'or'  => pass if any criterion is satisfied.\n"
           "  -i <file>             Input file (uses memory-mapped I/O for speed)\n"
           "  -q, --quiet           Suppress warnings\n"
           "  -h, --help            Show this help.\n\n"
           "Fields:\n"
           "  POS => numeric, QUAL => numeric, FILTER => string.\n"
           "  Others => assumed to be an INFO key. We try numeric parse if the criterion is numeric, else string.\n\n"
           "Performance:\n"
           "  Pass file directly for memory-mapped I/O (fastest).\n"
           "  Uses SIMD-optimized parsing on x86_64.\n"
           "  Zero-copy string_view parsing eliminates allocations.\n\n"
           "Example:\n"
           "  VCFX_record_filter --filter \"POS>=1000;FILTER==PASS;DP>10\" --logic and input.vcf\n"
           "  VCFX_record_filter -f \"QUAL>=30\" < in.vcf > out.vcf\n";
}

// ============================================================================
// Main entry point
// ============================================================================
int VCFXRecordFilter::run(int argc, char* argv[]) {
    std::string criteriaStr;
    std::string logicStr = "and";
    bool showHelp = false;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"filter", required_argument, 0, 'f'},
        {"logic", required_argument, 0, 'l'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    // Reset getopt
    optind = 1;

    while (true) {
        int c = ::getopt_long(argc, argv, "hf:l:i:q", long_opts, nullptr);
        if (c == -1) break;
        switch (c) {
            case 'h': showHelp = true; break;
            case 'f': criteriaStr = optarg; break;
            case 'l': logicStr = optarg; break;
            case 'i': inputFile_ = optarg; break;
            case 'q': quietMode_ = true; break;
            default: showHelp = true;
        }
    }

    // Check for positional file argument
    if (optind < argc && inputFile_.empty()) {
        inputFile_ = argv[optind];
    }

    if (showHelp || argc == 1) {
        displayHelp();
        return 0;
    }

    if (criteriaStr.empty()) {
        std::cerr << "Error: must provide --filter \"CRITERIA\".\n";
        displayHelp();
        return 1;
    }

    // Parse logic
    if (logicStr == "and") {
        useAndLogic_ = true;
    } else if (logicStr == "or") {
        useAndLogic_ = false;
    } else {
        std::cerr << "Error: logic must be 'and' or 'or'.\n";
        return 1;
    }

    // Parse and compile criteria
    if (!parseCriteria(criteriaStr, criteria_)) {
        std::cerr << "Error: failed to parse criteria.\n";
        return 1;
    }

    // Reserve space for reusable buffers
    infoTokens_.reserve(32);

    // Process input
    if (!inputFile_.empty() && inputFile_ != "-") {
        if (!processFileMmap(inputFile_.c_str())) {
            return 1;
        }
    } else {
        processStdin();
    }

    return 0;
}

// ============================================================================
// Legacy API compatibility functions
// ============================================================================
bool parseCriteria(const std::string& criteriaStr, std::vector<FilterCriterion>& criteria) {
    VCFXRecordFilter filter;
    return filter.parseCriteria(criteriaStr, criteria);
}

bool recordPasses(const std::string& record, const std::vector<FilterCriterion>& criteria, bool useAndLogic) {
    // Create a temporary filter with the criteria
    // This is less efficient but maintains API compatibility
    std::string_view line(record);

    if (useAndLogic) {
        for (const auto& c : criteria) {
            std::string_view fieldValue;

            switch (c.target) {
                case TargetField::POS: {
                    fieldValue = VCFXRecordFilter::extractField(line, 1);
                    if (fieldValue.empty()) return false;
                    double pos;
                    if (!VCFXRecordFilter::parseDouble(fieldValue, pos)) return false;
                    if (!VCFXRecordFilter::compareDouble(pos, c.op, c.numericValue)) return false;
                    break;
                }
                case TargetField::QUAL: {
                    fieldValue = VCFXRecordFilter::extractField(line, 5);
                    double qual = 0.0;
                    if (!fieldValue.empty() && fieldValue != ".") {
                        if (!VCFXRecordFilter::parseDouble(fieldValue, qual)) return false;
                    }
                    if (!VCFXRecordFilter::compareDouble(qual, c.op, c.numericValue)) return false;
                    break;
                }
                case TargetField::FILTER: {
                    fieldValue = VCFXRecordFilter::extractField(line, 6);
                    if (!VCFXRecordFilter::compareString(fieldValue, c.op, c.stringValue)) return false;
                    break;
                }
                case TargetField::INFO_KEY: {
                    std::string_view info = VCFXRecordFilter::extractField(line, 7);
                    std::string_view value;
                    if (!VCFXRecordFilter::extractInfoValue(info, c.fieldName, value)) return false;
                    if (c.fieldType == FieldType::NUMERIC) {
                        double num;
                        if (!VCFXRecordFilter::parseDouble(value, num)) return false;
                        if (!VCFXRecordFilter::compareDouble(num, c.op, c.numericValue)) return false;
                    } else {
                        if (!VCFXRecordFilter::compareString(value, c.op, c.stringValue)) return false;
                    }
                    break;
                }
            }
        }
        return true;
    } else {
        // OR logic
        for (const auto& c : criteria) {
            std::string_view fieldValue;
            bool pass = false;

            switch (c.target) {
                case TargetField::POS: {
                    fieldValue = VCFXRecordFilter::extractField(line, 1);
                    if (!fieldValue.empty()) {
                        double pos;
                        if (VCFXRecordFilter::parseDouble(fieldValue, pos)) {
                            pass = VCFXRecordFilter::compareDouble(pos, c.op, c.numericValue);
                        }
                    }
                    break;
                }
                case TargetField::QUAL: {
                    fieldValue = VCFXRecordFilter::extractField(line, 5);
                    double qual = 0.0;
                    if (!fieldValue.empty() && fieldValue != ".") {
                        VCFXRecordFilter::parseDouble(fieldValue, qual);
                    }
                    pass = VCFXRecordFilter::compareDouble(qual, c.op, c.numericValue);
                    break;
                }
                case TargetField::FILTER: {
                    fieldValue = VCFXRecordFilter::extractField(line, 6);
                    pass = VCFXRecordFilter::compareString(fieldValue, c.op, c.stringValue);
                    break;
                }
                case TargetField::INFO_KEY: {
                    std::string_view info = VCFXRecordFilter::extractField(line, 7);
                    std::string_view value;
                    if (VCFXRecordFilter::extractInfoValue(info, c.fieldName, value)) {
                        if (c.fieldType == FieldType::NUMERIC) {
                            double num;
                            if (VCFXRecordFilter::parseDouble(value, num)) {
                                pass = VCFXRecordFilter::compareDouble(num, c.op, c.numericValue);
                            }
                        } else {
                            pass = VCFXRecordFilter::compareString(value, c.op, c.stringValue);
                        }
                    }
                    break;
                }
            }

            if (pass) return true;
        }
        return false;
    }
}

void processVCF(std::istream& in, std::ostream& out, const std::vector<FilterCriterion>& criteria, bool useAndLogic) {
    std::string line;
    line.reserve(65536);
    bool foundChrom = false;

    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);

    while (std::getline(in, line)) {
        if (line.empty()) {
            outputBuffer.push_back('\n');
            continue;
        }
        if (line[0] == '#') {
            outputBuffer.append(line);
            outputBuffer.push_back('\n');
            if (line.rfind("#CHROM", 0) == 0) foundChrom = true;
            continue;
        }
        if (!foundChrom) {
            std::cerr << "Warning: data line before #CHROM => skipping.\n";
            continue;
        }
        if (recordPasses(line, criteria, useAndLogic)) {
            outputBuffer.append(line);
            outputBuffer.push_back('\n');
        }

        if (outputBuffer.size() > 512 * 1024) {
            out.write(outputBuffer.data(), outputBuffer.size());
            outputBuffer.clear();
        }
    }

    if (!outputBuffer.empty()) {
        out.write(outputBuffer.data(), outputBuffer.size());
    }
}

void printHelp() {
    VCFXRecordFilter filter;
    filter.displayHelp();
}

// ============================================================================
// Main
// ============================================================================
static void show_help() { printHelp(); }

int main(int argc, char* argv[]) {
    vcfx::init_io();

    if (vcfx::handle_common_flags(argc, argv, "VCFX_record_filter", show_help))
        return 0;

    VCFXRecordFilter filter;
    return filter.run(argc, argv);
}
