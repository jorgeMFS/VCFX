#include "VCFX_genotype_query.h"
#include "vcfx_core.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

// ------------------------------------------------------------------
// printHelp
// ------------------------------------------------------------------
void printHelp() {
    std::cout << "VCFX_genotype_query\n"
              << "Usage: VCFX_genotype_query [OPTIONS]\n\n"
              << "Options:\n"
              << "  --genotype-query, -g \"GENOTYPE\"  Specify the genotype to query (e.g., \"0/1\", \"1/1\").\n"
              << "  --strict                        Use strict string compare (no phasing unify or allele sorting).\n"
              << "  --help, -h                      Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin, outputs only the lines (plus all header lines) where\n"
              << "  at least one sample has the specified genotype in the 'GT' subfield.\n\n"
              << "Examples:\n"
              << "  # Flexible matching 0/1 or 0|1 => both become 0/1\n"
              << "  ./VCFX_genotype_query --genotype-query \"0/1\" < input.vcf > out.vcf\n\n"
              << "  # Strict matching => \"0|1\" won't match \"0/1\"\n"
              << "  ./VCFX_genotype_query --genotype-query \"0|1\" --strict < input.vcf > out.vcf\n";
}

// ------------------------------------------------------------------
// parseArguments
// ------------------------------------------------------------------
bool parseArguments(int argc, char *argv[], std::string &genotype_query, bool &strictCompare) {
    genotype_query.clear();
    strictCompare = false;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "--genotype-query" || arg == "-g") && i + 1 < argc) {
            genotype_query = argv[++i];
        } else if (arg.rfind("--genotype-query=", 0) == 0) {
            // e.g. --genotype-query=0/1
            genotype_query = arg.substr(17);
        } else if (arg == "--strict") {
            strictCompare = true;
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        }
    }
    return !genotype_query.empty();
}

// ------------------------------------------------------------------
// unifyGenotype(gt, strict):
//   - If strict==true, return gt unchanged.
//   - Else unify phasing by replacing '|' with '/', and sort alleles numerically
//   - Return empty if malformed
// ------------------------------------------------------------------
static std::string unifyGenotype(const std::string &gt, bool strict) {
    if (gt.empty())
        return "";

    // If user wants strict compare, do no changes
    if (strict) {
        return gt;
    }

    // unify phasing
    std::string g = gt;
    for (char &c : g) {
        if (c == '|')
            c = '/';
    }

    // Split on '/'
    std::vector<std::string> tokens;
    {
        std::stringstream ss(g);
        std::string t;
        while (std::getline(ss, t, '/')) {
            tokens.push_back(t);
        }
    }
    // If not diploid or missing, just return g
    if (tokens.size() != 2) {
        return g;
    }

    // Check numeric; if not numeric or ".", leave as is
    std::vector<int> vals;
    for (auto &tk : tokens) {
        if (tk == "." || tk.empty()) {
            // missing
            return g;
        }
        for (char c : tk) {
            if (!std::isdigit(static_cast<unsigned char>(c))) {
                return g;
            }
        }
        vals.push_back(std::stoi(tk));
    }
    // sort numeric
    std::sort(vals.begin(), vals.end());

    // reassemble
    std::stringstream out;
    out << vals[0] << "/" << vals[1];
    return out.str();
}

// ------------------------------------------------------------------
// genotypeQuery: stream-based approach
// ------------------------------------------------------------------
void genotypeQuery(std::istream &in, std::ostream &out, const std::string &genotype_query, bool strictCompare) {
    bool foundChrom = false;
    bool printedHeader = false;
    std::vector<std::string> headerLines;

    // We'll unify the user-specified genotype if !strict
    std::string unifiedQuery = unifyGenotype(genotype_query, strictCompare);
    if (unifiedQuery.empty()) {
        // fallback to whatever user typed
        unifiedQuery = genotype_query;
    }

    // Read the file line by line
    std::string line;
    while (true) {
        if (!std::getline(in, line)) {
            // EOF
            break;
        }
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // This is a header line
            headerLines.push_back(line);
            // Check if it's the #CHROM line
            if (line.rfind("#CHROM", 0) == 0) {
                foundChrom = true;
            }
        } else {
            // Data line
            if (!foundChrom) {
                std::cerr << "Error: No #CHROM header found before data lines.\n";
                return;
            }
            // Print the header lines once, right before the first data line
            if (!printedHeader) {
                for (const auto &h : headerLines) {
                    out << h << "\n";
                }
                printedHeader = true;
            }

            // Parse columns
            std::stringstream ss(line);
            std::vector<std::string> fields;
            {
                std::string token;
                while (std::getline(ss, token, '\t')) {
                    fields.push_back(token);
                }
            }
            if (fields.size() < 10) {
                std::cerr << "Warning: skipping line with <10 fields: " << line << "\n";
                continue;
            }

            // FORMAT column is fields[8]. We look for "GT" inside it
            std::string formatStr = fields[8];
            std::vector<std::string> fmts;
            {
                std::stringstream fmtsSS(formatStr);
                std::string f;
                while (std::getline(fmtsSS, f, ':')) {
                    fmts.push_back(f);
                }
            }
            int gtIndex = -1;
            for (int i = 0; i < (int)fmts.size(); i++) {
                if (fmts[i] == "GT") {
                    gtIndex = i;
                    break;
                }
            }
            if (gtIndex < 0) {
                // no GT => skip
                continue;
            }

            // Check each sample's GT
            bool match = false;
            for (size_t c = 9; c < fields.size(); c++) {
                std::stringstream sampSS(fields[c]);
                std::vector<std::string> sampleTokens;
                {
                    std::string x;
                    while (std::getline(sampSS, x, ':')) {
                        sampleTokens.push_back(x);
                    }
                }
                if (gtIndex >= (int)sampleTokens.size()) {
                    continue;
                }
                // unify genotype
                std::string gtVal = unifyGenotype(sampleTokens[gtIndex], strictCompare);
                if (gtVal == unifiedQuery) {
                    match = true;
                    break;
                }
            }

            if (match) {
                out << line << "\n";
            }
        }
    }

    // If we never found #CHROM, error out
    if (!foundChrom) {
        std::cerr << "Error: No #CHROM line found in VCF.\n";
        return;
    }
    // If we found #CHROM but no data lines at all,
    // we still need to print the header lines if not printed
    if (!printedHeader) {
        for (const auto &h : headerLines) {
            out << h << "\n";
        }
    }
}

// ------------------------------------------------------------------
// main
// ------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    if (vcfx::handle_common_flags(argc, argv, "VCFX_genotype_query", show_help))
        return 0;
    std::string genotypeQueryStr;
    bool strictCompare = false;
    if (!parseArguments(argc, argv, genotypeQueryStr, strictCompare)) {
        std::cerr << "Usage: " << argv[0] << " --genotype-query \"0/1\" [--strict] < input.vcf > output.vcf\n";
        std::cerr << "Use --help for usage.\n";
        return 1;
    }
    genotypeQuery(std::cin, std::cout, genotypeQueryStr, strictCompare);
    return 0;
}
