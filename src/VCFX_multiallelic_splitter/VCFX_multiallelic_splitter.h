#ifndef VCFX_MULTIALLELIC_SPLITTER_H
#define VCFX_MULTIALLELIC_SPLITTER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

/* Simple structure describing whether an ID is an INFO or FORMAT field, plus its Number field */
struct SubfieldMeta {
    bool isInfo = false, isFormat = false;
    std::string id, number;
};

/* We store an ID->SubfieldMeta map for known fields from the header. */
struct VCFHeaderInfo {
    std::unordered_map<std::string, SubfieldMeta> meta;
};

/* A single variant record with possibly multiple ALT alleles. */
struct VCFVariant {
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::vector<std::string> alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::string formatStr;
    std::vector<std::string> samples; // one element per sample
};

/* Print minimal help/usage. */
void printHelp();

/* Reads a VCF from 'in', writes lines to 'out' with multi-allelic sites split,
 * rewriting subfields (GT, AD, PL) for each splitted line.
 * Header lines are passed unmodified. */
bool splitMultiAllelicVariants(std::istream &in, std::ostream &out);

#endif
