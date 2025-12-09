#ifndef VCFX_RECORD_FILTER_H
#define VCFX_RECORD_FILTER_H

#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

// ============================================================================
// Filter operator types
// ============================================================================
enum class FilterOp : uint8_t {
    GT, // >
    GE, // >=
    LT, // <
    LE, // <=
    EQ, // ==
    NE  // !=
};

// Field types for comparison
enum class FieldType : uint8_t { NUMERIC, STRING };

// Target field type for optimized access
enum class TargetField : uint8_t {
    POS,      // Column 1 (numeric)
    QUAL,     // Column 5 (numeric)
    FILTER,   // Column 6 (string)
    INFO_KEY  // INFO field key lookup
};

// ============================================================================
// Compiled filter criterion (optimized for fast evaluation)
// ============================================================================
struct FilterCriterion {
    std::string fieldName;       // Original field name (for INFO keys)
    FilterOp op;
    double numericValue;         // Pre-parsed numeric threshold
    std::string stringValue;     // Pre-parsed string value
    FieldType fieldType;
    TargetField target;          // Compiled target for fast dispatch
};

// ============================================================================
// High-performance VCF Record Filter
// ============================================================================
class VCFXRecordFilter {
public:
    VCFXRecordFilter() = default;

    // Main entry point
    int run(int argc, char* argv[]);

    // Parse criteria from filter string
    bool parseCriteria(const std::string& criteriaStr, std::vector<FilterCriterion>& criteria);

    // Display help
    void displayHelp();

    // Fast field extraction (zero-copy) - public for legacy API compatibility
    static std::string_view extractField(std::string_view line, int fieldIndex);
    static bool extractInfoValue(std::string_view info, std::string_view key,
                                  std::string_view& valueOut);

    // Fast numeric parsing (no exceptions) - public for legacy API compatibility
    static bool parseDouble(std::string_view sv, double& out);
    static bool parseInt(std::string_view sv, int64_t& out);

    // Comparison functions - public for legacy API compatibility
    static bool compareDouble(double x, FilterOp op, double y);
    static bool compareString(std::string_view s, FilterOp op, std::string_view t);

private:
    // Configuration
    bool useAndLogic_ = true;
    bool quietMode_ = false;
    std::string inputFile_;

    // Pre-compiled criteria
    std::vector<FilterCriterion> criteria_;

    // Reusable buffers (avoid allocations in hot loop)
    std::vector<std::string_view> infoTokens_;

    // Processing methods
    bool processFileMmap(const char* filepath);
    void processStdin();

    // Core filtering (inlined for performance)
    bool evaluateLine(std::string_view line) const;
    bool evaluateCriterion(std::string_view line, const FilterCriterion& c) const;
};

// ============================================================================
// Legacy API compatibility (for existing code)
// ============================================================================
bool parseCriteria(const std::string& criteriaStr, std::vector<FilterCriterion>& criteria);
bool recordPasses(const std::string& record, const std::vector<FilterCriterion>& criteria, bool useAndLogic);
void processVCF(std::istream& in, std::ostream& out, const std::vector<FilterCriterion>& criteria, bool useAndLogic);
void printHelp();

#endif // VCFX_RECORD_FILTER_H
