#ifndef VCFX_MISSING_DATA_HANDLER_H
#define VCFX_MISSING_DATA_HANDLER_H

#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Structure to hold command-line arguments for the missing data handler tool.
 */
struct Arguments {
    bool fill_missing = false;                 ///< Flag indicating whether to impute missing genotypes.
    std::string default_genotype = "./.";      ///< Default genotype to use for imputation.
    std::vector<std::string> input_files;      ///< List of input VCF files (if any).
};

/**
 * @brief Displays the help message for the missing data handler tool.
 */
void printHelp();

/**
 * @brief Parses command-line arguments.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @param args Reference to Arguments structure to populate.
 * @return true if parsing is successful, false otherwise.
 */
bool parseArguments(int argc, char* argv[], Arguments& args);

/**
 * @brief Splits a string by a given delimiter.
 *
 * @param str The input string to split.
 * @param delimiter The character delimiter.
 * @return A vector of split substrings.
 */
std::vector<std::string> splitString(const std::string& str, char delimiter);

/**
 * @brief Processes the VCF file to handle missing genotype data.
 *
 * @param in Input stream (VCF file).
 * @param out Output stream (Modified VCF).
 * @param args Command-line arguments specifying behavior.
 * @return true if processing is successful, false otherwise.
 */
bool handleMissingData(std::istream& in, std::ostream& out, const Arguments& args);

#endif // VCFX_MISSING_DATA_HANDLER_H
