#include "VCFX_region_subsampler.h"
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cctype>
#include <cstdlib>

// A helper to merge intervals
static void mergeIntervals(std::vector<Region> &ivs){
    if(ivs.empty()) return;
    std::sort(ivs.begin(), ivs.end(),
              [](const Region &a, const Region &b){
                  return a.start < b.start;
              });
    std::vector<Region> merged;
    merged.push_back(ivs[0]);
    for(size_t i=1; i< ivs.size(); i++){
        Region &last= merged.back();
        const Region &curr= ivs[i];
        if(curr.start <= last.end+1){
            // overlap or contiguous
            if(curr.end> last.end) last.end= curr.end;
        } else {
            merged.push_back(curr);
        }
    }
    ivs= merged;
}

// For a sorted list of intervals, do a binary search to see if pos is in any
// typical approach: find interval with start <= pos, then check if pos<= that interval's end
static bool inRegions(const std::vector<Region> &ivs, int pos){
    // binary search approach
    int left=0, right= (int)ivs.size()-1;
    while(left<= right){
        int mid= (left+right)/2;
        if(pos < ivs[mid].start){
            right= mid-1;
        } else if(pos > ivs[mid].end){
            left= mid+1;
        } else {
            // pos in [ivs[mid].start, ivs[mid].end]
            return true;
        }
    }
    return false;
}

int VCFXRegionSubsampler::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    std::string bedFile;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"region-bed", required_argument, 0, 'b'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv,"hb:", long_opts,nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'b':
                bedFile= optarg;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    if(bedFile.empty()){
        std::cerr<<"Error: Must specify --region-bed <FILE>.\n";
        displayHelp();
        return 1;
    }
    // load
    if(!loadRegions(bedFile, regions)){
        std::cerr<<"Error: failed to load regions from "<<bedFile<<"\n";
        return 1;
    }
    // sort & merge intervals
    sortAndMergeIntervals(regions);

    // now process vcf from stdin
    processVCF(std::cin, std::cout);
    return 0;
}

void VCFXRegionSubsampler::displayHelp(){
    std::cout <<
"VCFX_region_subsampler: Keep only variants whose (CHROM,POS) is in a set of regions.\n\n"
"Usage:\n"
"  VCFX_region_subsampler --region-bed FILE < input.vcf > out.vcf\n\n"
"Options:\n"
"  -h, --help             Show help.\n"
"  -b, --region-bed FILE  BED file listing multiple regions.\n\n"
"Description:\n"
"  Reads the BED, which is <chrom> <start> <end> in 0-based. This tool converts\n"
"  them to 1-based [start+1 .. end]. Then merges intervals per chrom.\n"
"  Then only lines in the VCF that fall in those intervals for that CHROM are printed.\n\n"
"Example:\n"
"  VCFX_region_subsampler --region-bed myregions.bed < input.vcf > out.vcf\n";
}

bool VCFXRegionSubsampler::loadRegions(const std::string &bedFilePath,
                                       std::unordered_map<std::string, std::vector<Region>> &chromRegions)
{
    std::ifstream in(bedFilePath);
    if(!in.is_open()){
        std::cerr<<"Error: cannot open BED "<<bedFilePath<<"\n";
        return false;
    }
    std::string line;
    int lineCount=0;
    while(true){
        if(!std::getline(in,line)) break;
        lineCount++;
        if(line.empty()|| line[0]=='#') continue;
        std::stringstream ss(line);
        std::string chrom; 
        int start=0, end=0;
        if(!(ss>>chrom>>start>>end)){
            std::cerr<<"Warning: skipping invalid bed line "<<lineCount<<": "<< line<<"\n";
            continue;
        }
        // store => we do start+1 => 1-based
        if(start<0) start=0; 
        Region r;
        r.start= start+1;
        r.end= end;
        if(r.end< r.start){
            // ignore negative intervals
            continue;
        }
        chromRegions[chrom].push_back(r);
    }
    return true;
}

void VCFXRegionSubsampler::sortAndMergeIntervals(std::unordered_map<std::string, std::vector<Region>> &chromRegions){
    for(auto &kv: chromRegions){
        std::vector<Region> &ivs= kv.second;
        // sort & merge
        std::sort(ivs.begin(), ivs.end(), [](const Region &a, const Region &b){
            return a.start< b.start;
        });
        // merge
        std::vector<Region> merged;
        merged.reserve(ivs.size());
        merged.push_back(ivs[0]);
        for(size_t i=1; i< ivs.size(); i++){
            Region &last= merged.back();
            const Region &curr= ivs[i];
            if(curr.start <= last.end+1){
                if(curr.end> last.end) last.end= curr.end;
            } else {
                merged.push_back(curr);
            }
        }
        ivs= merged;
    }
}

// check if chrom is known => do binary search for pos in intervals
bool VCFXRegionSubsampler::isInAnyRegion(const std::string &chrom, int pos) const{
    auto it= regions.find(chrom);
    if(it== regions.end()) return false;
    // do binary search among intervals
    const auto &ivs= it->second;
    // typical approach:
    int left=0, right=(int)ivs.size()-1;
    while(left<= right){
        int mid= (left+right)/2;
        if(pos< ivs[mid].start) {
            right= mid-1;
        } else if(pos> ivs[mid].end){
            left= mid+1;
        } else {
            // pos in [ivs[mid].start..ivs[mid].end]
            return true;
        }
    }
    return false;
}

void VCFXRegionSubsampler::processVCF(std::istream &in, std::ostream &out){
    bool foundChromHeader=false;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line<<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line<<"\n";
            if(line.rfind("#CHROM",0)==0) foundChromHeader=true;
            continue;
        }
        if(!foundChromHeader){
            std::cerr<<"Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }
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
        std::string &chrom= fields[0];
        std::string &posStr= fields[1];
        int pos=0;
        try {
            pos= std::stoi(posStr);
        } catch(...){
            std::cerr<<"Warning: invalid POS => skipping.\n";
            continue;
        }
        // uppercase chrom if you want
        // check membership
        if(isInAnyRegion(chrom, pos)){
            out<< line<<"\n";
        }
    }
}


int main(int argc, char* argv[]){
    VCFXRegionSubsampler app;
    return app.run(argc, argv);
}
