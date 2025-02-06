#include "VCFX_genotype_query.h"
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>

// ------------------------------------------------------------------
// printHelp
// ------------------------------------------------------------------
void printHelp() {
    std::cout 
        << "VCFX_genotype_query\n"
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
bool parseArguments(int argc, char* argv[], std::string& genotype_query, bool &strictCompare) {
    genotype_query.clear();
    strictCompare = false;
    for (int i=1; i<argc; i++) {
        std::string arg = argv[i];
        if ((arg == "--genotype-query" || arg == "-g") && i+1<argc) {
            genotype_query = argv[++i];
        } else if (arg.rfind("--genotype-query=",0)==0) {
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
// Helper: unify genotype by replacing '|' with '/', then possibly sorting alleles if not strict
// e.g. "0|1" => "0/1", "2|1" => "1/2" if not strict
// returns "" if malformed
// ------------------------------------------------------------------
static std::string unifyGenotype(const std::string &gt, bool strict) {
    if (gt.empty()) return "";

    // If user wants strict compare, do no changes
    if (strict) {
        return gt;
    }
    // else unify phasing, sort alleles
    std::string g = gt;
    // unify separators
    for (char &c : g) {
        if (c=='|') c='/';
    }
    // split
    std::vector<std::string> tokens;
    {
        std::stringstream ss(g);
        std::string t;
        while (std::getline(ss,t,'/')) {
            tokens.push_back(t);
        }
    }
    if (tokens.size()<2) {
        // not diploid => can just return g or empty
        return g; 
    }
    // check if both are numeric
    std::vector<int> vals;
    vals.reserve(tokens.size());
    for (auto &tk : tokens) {
        if (tk=="." || tk.empty()) {
            // missing
            return g; // unify as is or return empty
        }
        // check numeric
        for (char c : tk) {
            if (!std::isdigit(c)) {
                // not numeric => just return g
                return g;
            }
        }
        vals.push_back(std::stoi(tk));
    }
    // sort numeric
    std::sort(vals.begin(), vals.end());
    // reassemble
    std::stringstream out;
    out << vals[0];
    for (size_t i=1; i<vals.size(); i++) {
        out << "/" << vals[i];
    }
    return out.str();
}

// ------------------------------------------------------------------
// genotypeQuery
// ------------------------------------------------------------------
void genotypeQuery(std::istream& in, std::ostream& out,
                   const std::string& genotype_query,
                   bool strictCompare) {
    std::string line;
    bool headerFound = false;

    // We'll unify the user-specified genotype if !strict
    std::string unifiedQuery = unifyGenotype(genotype_query, strictCompare);
    if (unifiedQuery.empty()) {
        // fallback
        unifiedQuery = genotype_query;
    }

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // Print all header lines
            if (!headerFound && line.rfind("#CHROM",0)==0) {
                headerFound = true;
            }
            out << line << "\n";
            continue;
        }
        if (!headerFound) {
            std::cerr << "Error: No #CHROM header found before data lines.\n";
            return;
        }
        // parse columns
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string token;
            while (std::getline(ss, token, '\t')) {
                fields.push_back(token);
            }
        }
        if (fields.size()<10) {
            // not enough columns for at least 1 sample
            std::cerr << "Warning: skipping line with <10 fields: " << line << "\n";
            continue;
        }
        // The 9th column is the format: e.g. GT:DP:...
        std::string formatStr = fields[8];
        // find GT index
        std::vector<std::string> fmts;
        {
            std::stringstream fmtsSS(formatStr);
            std::string f;
            while (std::getline(fmtsSS, f, ':')) {
                fmts.push_back(f);
            }
        }
        int gtIndex = -1;
        for (int i=0; i<(int)fmts.size(); i++) {
            if (fmts[i] == "GT") {
                gtIndex = i;
                break;
            }
        }
        if (gtIndex<0) {
            // no GT => skip
            continue;
        }
        bool match = false;
        // from column 10 onward are sample columns
        for (size_t c=9; c<fields.size(); c++) {
            // parse sample data by ':'
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

int main(int argc, char* argv[]) {
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
