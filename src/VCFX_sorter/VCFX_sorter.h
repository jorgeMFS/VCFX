#ifndef VCFX_SORTER_H
#define VCFX_SORTER_H

#include <string>
#include <vector>

struct VCFRecord {
    std::string chrom;
    int pos;
    std::vector<std::string> fields; // store entire splitted line so we can rebuild

    // We'll do a custom comparator that depends on whether we do "natural" or "lex"
    static bool lexCompare(const VCFRecord &a, const VCFRecord &b);
    static bool naturalCompare(const VCFRecord &a, const VCFRecord &b);
};

// A class to hold main logic
class VCFXSorter {
public:
    int run(int argc, char* argv[]);

private:
    // parse arguments, then read lines
    // store header lines, store data lines in memory
    // sort data lines
    // output header + sorted lines
    void displayHelp();
    void loadVCF(std::istream &in);
    void sortRecords();
    void outputVCF(std::ostream &out);

    // read an environment var or a CLI option to decide lexicographic or natural
    bool naturalChromOrder= false;

    // store all lines that begin with '#'
    std::vector<std::string> headerLines;
    // store data lines in memory
    std::vector<VCFRecord> records;
};

#endif
