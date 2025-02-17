#include "VCFX_sample_extractor.h"
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cctype>

// Helper to trim whitespace
static std::string trim(const std::string &s){
    size_t start=0; 
    while(start<s.size() && isspace((unsigned char)s[start])) start++;
    if(start== s.size()) return "";
    size_t end= s.size()-1;
    while(end> start && isspace((unsigned char)s[end])) end--;
    return s.substr(start, end-start+1);
}

int VCFXSampleExtractor::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    std::vector<std::string> samples;

    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"samples", required_argument, 0, 's'},
        {0,0,0,0}
    };
    while(true){
        int c=::getopt_long(argc, argv, "hs:", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 's': {
                // parse comma or space delimited
                std::stringstream ss(optarg);
                // we allow e.g. --samples "SampleA,SampleB"
                // or --samples "SampleA SampleB"
                std::string token;
                while(ss>> token){
                    // also check for comma splitting
                    std::stringstream s2(token);
                    std::string sub;
                    while(std::getline(s2, sub, ',')){
                        sub= trim(sub);
                        if(!sub.empty()) samples.push_back(sub);
                    }
                }
            }break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    if(samples.empty()){
        std::cerr<<"Error: must specify at least one sample with --samples.\n";
        return 1;
    }
    // do it
    extractSamples(std::cin, std::cout, samples);
    return 0;
}

void VCFXSampleExtractor::displayHelp(){
    std::cout <<
"VCFX_sample_extractor: Subset a VCF to a chosen set of samples.\n\n"
"Usage:\n"
"  VCFX_sample_extractor --samples \"Sample1,Sample2\" < input.vcf > output.vcf\n\n"
"Options:\n"
"  -h, --help              Print this help.\n"
"  -s, --samples <LIST>    Comma or space separated list of sample names.\n\n"
"Description:\n"
"  Reads #CHROM line to identify sample columns. Keeps only user-specified samples.\n"
"  Rewrites #CHROM line with that subset. For each variant data line, we keep only the\n"
"  chosen sample columns. If a sample isn't found in the header, logs a warning.\n\n"
"Example:\n"
"  VCFX_sample_extractor --samples \"IndivA IndivB\" < input.vcf > subset.vcf\n";
}

void VCFXSampleExtractor::extractSamples(std::istream &in, std::ostream &out,
                                         const std::vector<std::string> &samples)
{
    // store the samples in a set for quick membership check
    std::unordered_set<std::string> sampleSet(samples.begin(), samples.end());

    bool foundChromLine=false;
    std::string line;
    // We'll keep track of the sample name -> index
    // We'll build a vector of indices we want to keep in data lines
    // Also track how many of user’s samples were found
    std::vector<int> keepSampleIndices;
    keepSampleIndices.reserve(samples.size());
    // We'll also store them in the final #CHROM order
    std::vector<std::string> finalSampleNames;
    bool headerProcessed=false;

    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            // pass header lines
            if(line.rfind("#CHROM",0)==0){
                foundChromLine= true;
                // parse
                std::stringstream ss(line);
                std::vector<std::string> headerFields;
                {
                    std::string col;
                    while(std::getline(ss, col, '\t')){
                        headerFields.push_back(col);
                    }
                }
                // sample columns => from index=9..end
                // find which ones are in sampleSet
                keepSampleIndices.clear();
                finalSampleNames.clear();
                if(headerFields.size()<9){
                    // means no samples in this VCF => skip
                    out<< line<<"\n";
                    continue;
                }
                // standard columns 0..8 => CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT
                for(size_t i=9;i< headerFields.size(); i++){
                    // if it's in sampleSet => keep
                    if(sampleSet.count(headerFields[i])>0){
                        keepSampleIndices.push_back((int)i);
                        finalSampleNames.push_back(headerFields[i]);
                    }
                }
                // for any sample in sampleSet that wasn't found => warn
                for(auto &s: samples){
                    if(std::find(finalSampleNames.begin(), finalSampleNames.end(), s) == finalSampleNames.end()){
                        std::cerr<<"Warning: sample '"<< s <<"' not found in header.\n";
                    }
                }
                // now rewrite the #CHROM line
                // keep the first 9 columns, then only the chosen sample columns
                std::ostringstream newHeader;
                // 0..8
                for(int i=0;i<9;i++){
                    if(i>0) newHeader<<"\t";
                    newHeader<< headerFields[i];
                }
                // then each sample
                for(size_t i=0;i< finalSampleNames.size(); i++){
                    newHeader<<"\t"<< finalSampleNames[i];
                }
                out<< newHeader.str()<<"\n";
                headerProcessed=true;
            } else {
                out<< line<<"\n";
            }
            continue;
        }
        if(!foundChromLine){
            std::cerr<<"Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }
        // parse data line
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string col;
            while(std::getline(ss,col,'\t')){
                fields.push_back(col);
            }
        }
        if(fields.size()<8){
            std::cerr<<"Warning: line has <8 columns => skipping.\n";
            continue;
        }
        // if there's no sample columns => we can just keep line, but that wouldn't make sense
        // actually if there's < 10 columns => no sample data. Then there's no subset to do
        // We'll keep it if there's no sample columns? It's up to us. We'll do the standard approach:
        // we need at least 9 => CHROM..FORMAT plus at least one sample => if not => just print as is
        if(fields.size()<9){
            // no sample columns
            // keep or skip? We'll keep it to keep it consistent. Or we skip because it's an invalid VCF?
            // We'll just skip
            std::cerr<<"Warning: data line with no sample columns => skipping.\n";
            continue;
        }
        // now we keep 0..8 plus only the chosen samples
        std::ostringstream newLine;
        // 0..8
        for(int i=0;i<9;i++){
            if(i>0) newLine<<"\t";
            newLine<< fields[i];
        }
        // now each chosen sample
        for(auto idx: keepSampleIndices){
            if((size_t)idx < fields.size()){
                newLine<<"\t"<< fields[idx];
            } else {
                // out of range => sample column not present => put "."
                newLine<<"\t.";
            }
        }
        out<< newLine.str()<<"\n";
    }
}
