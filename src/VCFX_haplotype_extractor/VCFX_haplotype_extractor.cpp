#include "vcfx_core.h"
#include "VCFX_haplotype_extractor.h"
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cstdlib>

// Constructor not needed here, we do no special initialization
// HaplotypeExtractor::HaplotypeExtractor() {}

// ---------------------------------------------------------------------
// printHelp
// ---------------------------------------------------------------------
void printHelp() {
    std::cout << "VCFX_haplotype_extractor\n"
              << "Usage: VCFX_haplotype_extractor [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h                 Display this help message and exit.\n"
             << "  --block-size <int>         Maximum distance for grouping consecutive variants (default 100000).\n"
             << "  --check-phase-consistency  If set, try a minimal check across variants.\n"
             << "  --debug                   Output verbose debug information.\n\n"
             << "Description:\n"
              << "  Extracts phased haplotype blocks from genotype data in a VCF file. "
              << "It reconstructs haplotypes for each sample by analyzing phased genotype fields.\n\n"
              << "Examples:\n"
              << "  ./VCFX_haplotype_extractor --block-size 50000 < phased.vcf > haplotypes.tsv\n"
              << "  ./VCFX_haplotype_extractor --check-phase-consistency < phased.vcf > haplotypes.tsv\n";
}

// ---------------------------------------------------------------------
// parseHeader: parse #CHROM line => sample columns
// ---------------------------------------------------------------------
bool HaplotypeExtractor::parseHeader(const std::string& headerLine) {
    auto fields = splitString(headerLine, '\t');
    if (fields.size() <= 8) {
        std::cerr << "Error: VCF header does not contain sample columns.\n";
        return false;
    }
    // from col=9 onward => sample names
    sampleNames.assign(fields.begin() + 9, fields.end());
    numSamples = sampleNames.size();
    return true;
}

// ---------------------------------------------------------------------
// splitString
// ---------------------------------------------------------------------
std::vector<std::string> HaplotypeExtractor::splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string tmp;
    while (std::getline(ss, tmp, delimiter)) {
        tokens.push_back(tmp);
    }
    return tokens;
}

// ---------------------------------------------------------------------
// areAllSamplesPhased: checks each genotype for '|'
// ---------------------------------------------------------------------
bool HaplotypeExtractor::areAllSamplesPhased(const std::vector<std::string>& genotypes) {
    for (auto &g : genotypes) {
        if (g.find('|') == std::string::npos) {
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------
// phaseIsConsistent: minimal check across variants
//   e.g. if new variant has "0|1" but existing block's last variant for that sample was "1|0"
//   we might call that inconsistent if we want a simplistic approach. 
//   This logic is trivially commented out or replaced with a more advanced approach.
// ---------------------------------------------------------------------
bool HaplotypeExtractor::phaseIsConsistent(const HaplotypeBlock& block,
                                           const std::vector<std::string>& newGenotypes)
{
    // We'll do a naive check: if the new genotype is e.g. "0|1" 
    // but the last appended genotype portion for the sample is "1|0", we call it inconsistent.
    // Actually, we store haplotypes with each variant appended by "|" + newGT?
    // We'll parse the last variant's GT for each sample from the block's haplotypes.
    // This is simplistic.

    if (block.haplotypes.size() != newGenotypes.size()) {
        // mismatch in #samples => inconsistent
        return false;
    }
    
    // Optional debugging output
    if (debugMode) {
        std::cerr << "Checking phase consistency\n";
    }

    for (size_t s=0; s<block.haplotypes.size(); s++) {
        // block's haplotypes[s] is a big string with variants separated by '|', e.g. "0|1|0|1"
        // let's get the last chunk after the last '|'
        const std::string &allVar = block.haplotypes[s];
        
        // The haplotype string is stored like: "0|1|1|0|0|1" where each pair is one genotype
        // We need to extract the last genotype, which is the last 3 characters
        std::string lastGT;
        if (allVar.size() < 3) {
            lastGT = allVar;
        } else {
            // Get the last genotype (3 characters: allele|allele)
            size_t lastPipePos = allVar.rfind('|', allVar.size() - 2);
            if (lastPipePos == std::string::npos) {
                // No pipe found in position except the last one, so this is the first genotype
                lastGT = allVar;
            } else {
                // Extract the last genotype (e.g., "0|1" from "0|1|1|0")
                lastGT = allVar.substr(lastPipePos - 1);
            }
        }
        
        if (debugMode) {
            std::cerr << "Sample " << s << " last GT: " << lastGT << " new GT: " << newGenotypes[s] << "\n";
        }

        // compare lastGT with newGenotypes[s]
        // if they differ in a 2-allele reversed manner, we might call it inconsistent
        // e.g. lastGT="0|1", new="1|0"
        // We'll interpret them as numeric pairs, ignoring '.' or weirdness
        // If can't parse => skip checking
        if (lastGT.size()<3 || newGenotypes[s].size()<3) {
            continue;
        }
        
        // Get the first and second alleles for comparison
        char lastAllele1 = (lastGT.size() >= 1) ? lastGT[lastGT.size() - 3] : '.';
        char lastAllele2 = (lastGT.size() >= 1) ? lastGT[lastGT.size() - 1] : '.';
        char newAllele1 = newGenotypes[s][0];
        char newAllele2 = newGenotypes[s][2];
        
        if (debugMode) {
            std::cerr << "Comparing alleles: " << lastAllele1 << "|" << lastAllele2
                      << " vs " << newAllele1 << "|" << newAllele2 << "\n";
        }
        
        // If e.g. lastGT=="0|1", newGenotypes=="1|0" => inconsistent
        // Check for phase flips - when both alleles flip positions
        if (lastAllele1 != newAllele1 && lastAllele2 != newAllele2 && 
            lastAllele1 == newAllele2 && lastAllele2 == newAllele1) {
            if (debugMode) {
                std::cerr << "Phase flip detected in sample " << s << "\n";
            }
            return false;
        }
    }

    if (debugMode) {
        std::cerr << "All phases consistent\n";
    }
    return true;
}

// ---------------------------------------------------------------------
// updateBlocks: merges new variant into last block or starts new block
//   we rely on the blockDistanceThreshold and (optionally) checkPhaseConsistency
// ---------------------------------------------------------------------
void HaplotypeExtractor::updateBlocks(std::vector<HaplotypeBlock>& haplotypeBlocks,
                                      const std::string& chrom, int pos,
                                      const std::vector<std::string>& genotypes)
{
    if (haplotypeBlocks.empty()) {
        // start first block
        HaplotypeBlock b;
        b.chrom = chrom;
        b.start = pos;
        b.end   = pos;
        // haplotypes => each sample's genotype
        b.haplotypes = genotypes;
        haplotypeBlocks.push_back(std::move(b));
        return;
    }

    // examine last block
    HaplotypeBlock &lastB = haplotypeBlocks.back();
    // if same chrom, and pos-lastB.end <= blockDistanceThreshold
    // and if checkPhaseConsistency => phaseIsConsistent
    bool canExtend = (chrom == lastB.chrom) &&
                     (pos - lastB.end <= blockDistanceThreshold);

    if (canExtend && checkPhaseConsistency) {
        canExtend = phaseIsConsistent(lastB, genotypes);
    }

    if (canExtend) {
        // extend
        lastB.end = pos;
        // for each sample, append "|" + new genotype
        // if they are guaranteed phased, we can do that
        for (size_t s=0; s<lastB.haplotypes.size(); s++) {
            lastB.haplotypes[s] += "|" + genotypes[s];
        }
    } else {
        // start new block
        HaplotypeBlock nb;
        nb.chrom = chrom;
        nb.start = pos;
        nb.end   = pos;
        nb.haplotypes = genotypes;
        haplotypeBlocks.push_back(std::move(nb));
    }
}

// ---------------------------------------------------------------------
// processVariant: parse data line, check phased, update blocks
// ---------------------------------------------------------------------
bool HaplotypeExtractor::processVariant(const std::vector<std::string>& fields,
                                        std::vector<HaplotypeBlock>& haplotypeBlocks)
{
    if (fields.size()<9) {
        std::cerr << "Warning: skipping invalid VCF line (<9 fields)\n";
        return false;
    }
    const std::string &chrom = fields[0];
    int pos=0;
    try {
        pos = std::stoi(fields[1]);
    } catch(...) {
        std::cerr << "Warning: invalid POS => skip variant\n";
        return false;
    }
    // find GT index
    const std::string &formatStr = fields[8];
    auto formatToks = splitString(formatStr, ':');
    int gtIndex = -1;
    for (int i=0; i<(int)formatToks.size(); i++) {
        if (formatToks[i]=="GT") {
            gtIndex = i;
            break;
        }
    }
    if (gtIndex<0) {
        // no GT => skip
        return false;
    }
    // gather genotypes
    // from col=9 onward => sample columns
    std::vector<std::string> genotypeFields(numSamples);
    bool allPhased = true;
    for (size_t s=0; s<numSamples; s++) {
        if (9+s >= fields.size()) {
            // missing sample => set .|.
            genotypeFields[s] = ".|.";
            allPhased = false;
            continue;
        }
        auto sampleToks = splitString(fields[9+s], ':');
        if (gtIndex >= (int)sampleToks.size()) {
            genotypeFields[s] = ".|.";
            allPhased = false;
            continue;
        }
        // e.g. "0|1"
        std::string gt = sampleToks[gtIndex];
        if (gt.find('|')==std::string::npos) {
            allPhased = false;
        }
        genotypeFields[s] = gt;
    }

    if (!allPhased) {
        // for demonstration, we skip if not fully phased
        // or we can allow partial. We'll skip
        std::cerr << "Warning: Not all samples phased at " << chrom << ":" << pos << ".\n";
        return false;
    }

    // if we get here => we have a fully phased variant => update blocks
    updateBlocks(haplotypeBlocks, chrom, pos, genotypeFields);
    return true;
}

// ---------------------------------------------------------------------
// extractHaplotypes: main function to read from 'in', produce blocks => 'out'
// ---------------------------------------------------------------------
bool HaplotypeExtractor::extractHaplotypes(std::istream& in, std::ostream& out) {
    bool foundHeader = false;
    std::vector<HaplotypeBlock> haplotypeBlocks;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0]=='#') {
            // if #CHROM => parse
            if (!foundHeader && line.rfind("#CHROM",0)==0) {
                foundHeader = true;
                if (!parseHeader(line)) {
                    return false;
                }
            }
            // skip output
            continue;
        }
        if (!foundHeader) {
            std::cerr << "Error: no #CHROM header found before data.\n";
            return false;
        }
        // parse columns
        auto fields = splitString(line, '\t');
        processVariant(fields, haplotypeBlocks);
    }

    // output
    //  CHROM   START   END   Sample1   Sample2 ...
    out << "CHROM\tSTART\tEND";
    for (auto &s : sampleNames) {
        out << "\t" << s;
    }
    out << "\n";

    for (auto &block : haplotypeBlocks) {
        out << block.chrom << "\t" << block.start << "\t" << block.end;
        for (auto &hap : block.haplotypes) {
            out << "\t" << hap;
        }
        out << "\n";
    }

    return true;
}

// ---------------------------------------------------------------------
// main
// ---------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (vcfx::handle_version_flag(argc, argv, "VCFX_haplotype_extractor")) return 0;
    int blockSize = 100000;
    bool doCheck = false;
    bool debug = false;

    // simple arg parse
    for (int i=1; i<argc; i++) {
        std::string a = argv[i];
        if (a=="--help" || a=="-h") {
            printHelp();
            return 0;
        } else if (a=="--block-size" && i+1<argc) {
            blockSize = std::stoi(argv[++i]);
        } else if (a=="--check-phase-consistency") {
            doCheck = true;
        } else if (a=="--debug") {
            debug = true;
        }
    }

    HaplotypeExtractor extractor;
    extractor.setBlockDistanceThreshold(blockSize);
    extractor.setCheckPhaseConsistency(doCheck);
    extractor.setDebug(debug);

    bool ok = extractor.extractHaplotypes(std::cin, std::cout);
    return (ok ? 0 : 1);
}
