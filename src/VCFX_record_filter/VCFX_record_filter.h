#ifndef VCFX_RECORD_FILTER_H
#define VCFX_RECORD_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// Operator kinds
enum class FilterOp {
    GT,  // >
    GE,  // >=
    LT,  // <
    LE,  // <=
    EQ,  // ==
    NE   // !=
};

// We store whether the field is numeric or string, as we do different compares
enum class FieldType {
    NUMERIC,
    STRING
};

// A single criterion: e.g. "POS > 100", or "FILTER == PASS", or "AF >= 0.1".
struct FilterCriterion {
    std::string fieldName;
    FilterOp op;
    double numericValue;    // used if fieldType==NUMERIC
    std::string stringValue;// used if fieldType==STRING
    FieldType fieldType;
};

// parse multiple criteria from a single string (separated by semicolon)
bool parseCriteria(const std::string &criteriaStr,
                   std::vector<FilterCriterion> &criteria);

// apply the (AND/OR) logic to a single VCF line
bool recordPasses(const std::string &record,
                  const std::vector<FilterCriterion> &criteria,
                  bool useAndLogic);

// read lines from 'in', filter, write pass lines to 'out'
void processVCF(std::istream &in,
                std::ostream &out,
                const std::vector<FilterCriterion> &criteria,
                bool useAndLogic);

// display usage
void printHelp();

#endif
