#include "VCFX_outlier_detector.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdlib>

static void split(const std::string &s, char d, std::vector<std::string> &tokens) {
    tokens.clear();
    std::stringstream ss(s);
    std::string t;
    while(std::getline(ss,t,d)) {
        tokens.push_back(t);
    }
}

int VCFXOutlierDetector::run(int argc, char* argv[]){
    bool showHelp=false;
    std::string metric="AF";
    double threshold= 0.0;
    bool isVariantMode= true;

    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"metric", required_argument, 0, 'm'},
        {"threshold", required_argument, 0, 't'},
        {"variant", no_argument, 0, 'v'},
        {"sample", no_argument, 0, 's'},
        {0,0,0,0}
    };
    while(true){
        int c= getopt_long(argc, argv, "hm:t:vs", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true; 
                break;
            case 'm':
                metric= optarg; 
                break;
            case 't':
                try {
                    threshold= std::stod(optarg);
                } catch(...){
                    std::cerr<< "Error: invalid threshold.\n";
                    return 1;
                }
                break;
            case 'v':
                isVariantMode= true;
                break;
            case 's':
                isVariantMode= false;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp|| threshold<=0.0){
        displayHelp();
        return 0;
    }
    detectOutliers(std::cin, std::cout, metric, threshold, isVariantMode);
    return 0;
}

void VCFXOutlierDetector::displayHelp(){
    std::cout <<
"VCFX_outlier_detector: Identify outliers among variants or samples based on a numeric metric.\n\n"
"Usage:\n"
"  VCFX_outlier_detector --metric <KEY> --threshold <VAL> [--variant|--sample]\n"
"  < input.vcf > out\n\n"
"Options:\n"
"  --help, -h           Print this help.\n"
"  --metric, -m <KEY>   Name of the metric to use (e.g. AF, DP, GQ...).\n"
"  --threshold, -t <VAL> Numeric threshold.\n"
"  --variant, -v        Evaluate each variant's <KEY> in INFO>threshold => print.\n"
"  --sample, -s         Evaluate sample averages of <KEY> in genotype subfield => print outliers.\n\n"
"Examples:\n"
"  1) Outlier variants with AF>0.05:\n"
"     VCFX_outlier_detector --metric AF --threshold 0.05 --variant < in.vcf > out.txt\n"
"  2) Outlier samples if average GQ>30:\n"
"     VCFX_outlier_detector --metric GQ --threshold 30 --sample < in.vcf > sample_outliers.txt\n";
}

bool VCFXOutlierDetector::parseMetricFromInfo(const std::string &info,
                                              const std::string &key,
                                              double &val) const
{
    if(info.empty()||info=="." ) return false;
    std::stringstream ss(info);
    std::string kv;
    while(std::getline(ss, kv, ';')) {
        auto eq= kv.find('=');
        if(eq==std::string::npos) continue;
        auto k= kv.substr(0,eq);
        auto v= kv.substr(eq+1);
        if(k==key) {
            try {
                val= std::stod(v);
                return true;
            } catch(...){
                return false;
            }
        }
    }
    return false;
}

bool VCFXOutlierDetector::parseMetricFromGenotype(const std::string &genotypeField,
                                                  const std::string &metric,
                                                  double &value) const
{
    std::stringstream ss(genotypeField);
    std::string token;
    while(std::getline(ss, token, ':')) {
        auto eq= token.find('=');
        if(eq!=std::string::npos) {
            auto left= token.substr(0,eq);
            auto right= token.substr(eq+1);
            if(left== metric) {
                try {
                    value= std::stod(right);
                    return true;
                } catch(...) {
                    return false;
                }
            }
        }
    }
    return false;
}

void VCFXOutlierDetector::detectOutliers(std::istream &in,
                                         std::ostream &out,
                                         const std::string &metric,
                                         double threshold,
                                         bool isVariantMode)
{
    if(isVariantMode) {
        out<<"#CHROM\tPOS\tID\t"<<metric<<"\n";
        bool headerFound=false;
        std::string line;
        bool anyMetricFound= false;
        while(true){
            if(!std::getline(in,line)) break;
            if(line.empty())continue;
            if(line[0]=='#'){
                if(line.rfind("#CHROM",0)==0) headerFound=true;
                continue;
            }
            if(!headerFound){
                // skip or pass, but we can't parse columns
                continue;
            }
            std::vector<std::string> fields; {
                std::stringstream ss(line);
                std::string f;
                while(std::getline(ss,f,'\t')) fields.push_back(f);
            }
            if(fields.size()<8) continue;
            std::string &chrom=fields[0];
            std::string &pos= fields[1];
            std::string &id= fields[2];
            std::string &info= fields[7];
            double val= 0.0;
            if(parseMetricFromInfo(info, metric, val)) {
                anyMetricFound= true;
                if(val> threshold) {
                    out<< chrom<<"\t"<< pos<<"\t"<< id<<"\t"<< val<<"\n";
                }
            }
        }
        if(!anyMetricFound) {
            std::cerr<<"Warning: metric '"<<metric<<"' not found in any INFO field.\n";
        }
    } else {
        // sample approach: we read entire file, gather sample names => sum metric => count => compute average => compare threshold
        bool headerFound=false;
        std::vector<std::string> sampleNames;
        std::unordered_map<std::string, double> sums;
        std::unordered_map<std::string, int> counts;
        std::string line;
        bool anyMetricFound= false;
        while(true){
            if(!std::getline(in,line)) break;
            if(line.empty())continue;
            if(line[0]=='#'){
                if(line.rfind("#CHROM",0)==0){
                    headerFound=true;
                    std::vector<std::string> f; {
                        std::stringstream ss(line);
                        std::string x;
                        while(std::getline(ss,x,'\t')) f.push_back(x);
                    }
                    for(size_t i=9;i<f.size();i++){
                        sampleNames.push_back(f[i]);
                        sums[f[i]]=0.0;
                        counts[f[i]]=0;
                    }
                }
                continue;
            }
            if(!headerFound) {
                continue;
            }
            std::vector<std::string> fields; {
                std::stringstream ss(line);
                std::string fx;
                while(std::getline(ss, fx, '\t')) fields.push_back(fx);
            }
            if(fields.size()<9) continue;
            std::string &format= fields[8];
            // parse format subfields? We'll do a quick parse
            // If the user is looking for "DP=###" or "GQ=###" => parse parseMetricFromGenotype
            // We do that for each sample
            // We do not confirm the 'metric' is in the format header => naive approach
            // We'll attempt parse
            std::vector<std::string> sampleColumns;
            for(size_t i=9;i<fields.size();i++){
                sampleColumns.push_back(fields[i]);
            }
            for(size_t s=0;s<sampleNames.size() && s<sampleColumns.size();s++){
                double val=0.0;
                if(parseMetricFromGenotype(sampleColumns[s], metric, val)) {
                    sums[sampleNames[s]]+= val;
                    counts[sampleNames[s]]++;
                    anyMetricFound= true;
                }
            }
        }
        out<<"#Sample\tAverage_"<< metric <<"\n";
        for(auto &nm: sampleNames) {
            if(counts[nm]==0) {
                out<< nm<<"\tNA\n";
            } else {
                double avg= sums[nm]/ (double)counts[nm];
                if(avg> threshold) {
                    out<< nm<<"\t"<< avg<<"\n";
                }
            }
        }
        if(!anyMetricFound) {
            std::cerr<<"Warning: metric '"<<metric<<"' was not found in any sample genotype.\n";
        }
    }
}
