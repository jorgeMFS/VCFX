./VCFX_nonref_filter/VCFX_nonref_filter.cpp
#include "VCFX_nonref_filter.h"
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <vector>
#include <string>

int VCFXNonRefFilter::run(int argc, char* argv[]){
    bool showHelp=false;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };
    while(true){
        int c= getopt_long(argc, argv, "h", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
            default:
                showHelp= true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    filterNonRef(std::cin, std::cout);
    return 0;
}

void VCFXNonRefFilter::displayHelp(){
    std::cout <<
"VCFX_nonref_filter: Exclude variants if all samples are homozygous reference.\n\n"
"Usage:\n"
"  VCFX_nonref_filter [options] < input.vcf > output.vcf\n\n"
"Description:\n"
"  Reads VCF lines. For each variant, we check each sample's genotype. If a\n"
"  genotype is polyploid, all alleles must be '0'. If a genotype is missing\n"
"  or partial, we consider it not guaranteed hom-ref => keep variant.\n"
"  If we find at least one sample not hom-ref, we print the variant. Otherwise,\n"
"  we skip it.\n\n"
"Example:\n"
"  VCFX_nonref_filter < input.vcf > filtered.vcf\n\n";
}

bool VCFXNonRefFilter::isDefinitelyHomRef(const std::string &genotypeField) const {
    if(genotypeField.empty() || genotypeField=="." || genotypeField=="./." || genotypeField==".|.") return false;
    std::string g= genotypeField;
    for(char &c:g) if(c=='|') c='/';
    // split by '/'
    std::vector<std::string> alleles;
    {
        std::stringstream ss(g);
        std::string tok;
        while(std::getline(ss, tok, '/')) alleles.push_back(tok);
    }
    if(alleles.empty()) return false;
    // if any allele is not "0", => not homRef
    // if allele is "." => missing => not guaranteed homRef => false
    for(auto &al : alleles){
        if(al!="0") return false;
    }
    return true;
}

void VCFXNonRefFilter::filterNonRef(std::istream& in, std::ostream& out){
    bool headerFound=false;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line << "\n";
            if(line.rfind("#CHROM",0)==0) headerFound=true;
            continue;
        }
        if(!headerFound){
            std::cerr<<"Warning: VCF data line encountered before #CHROM. Passing line.\n";
            out<< line <<"\n";
            continue;
        }
        std::stringstream ss(line);
        std::vector<std::string> fields; {
            std::string f;
            while(std::getline(ss,f,'\t')) fields.push_back(f);
        }
        if(fields.size()<10){
            out<< line <<"\n";
            continue;
        }
        std::string formatStr= fields[8];
        std::vector<std::string> fmtParts;
        {
            std::stringstream fs(formatStr);
            std::string ff;
            while(std::getline(fs,ff,':')) fmtParts.push_back(ff);
        }
        int gtIndex=-1;
        for(size_t i=0;i<fmtParts.size();i++){
            if(fmtParts[i]=="GT"){ gtIndex=(int)i; break;}
        }
        if(gtIndex<0){
            // no genotype => cannot confirm all hom-ref => we keep
            out<< line <<"\n";
            continue;
        }
        bool allHomRef=true;
        for(size_t s=9; s< fields.size(); s++){
            std::string &sampleCol= fields[s];
            std::vector<std::string> subf;
            {
                std::stringstream sampleSS(sampleCol);
                std::string token;
                while(std::getline(sampleSS, token, ':')) subf.push_back(token);
            }
            if(gtIndex>=(int)subf.size()){
                allHomRef= false; 
                break;
            }
            if(!isDefinitelyHomRef(subf[gtIndex])){
                allHomRef= false;
                break;
            }
        }
        if(!allHomRef) out<< line <<"\n";
    }
}

int main(int argc, char* argv[]){
    VCFXNonRefFilter app;
    return app.run(argc, argv);
}


./VCFX_nonref_filter/VCFX_nonref_filter.h
#ifndef VCFX_NONREF_FILTER_H
#define VCFX_NONREF_FILTER_H

#include <iostream>
#include <string>
#include <vector>

/*
 * VCFX_nonref_filter:
 *   Reads a VCF, discards any variant for which EVERY sample is hom-ref
 *   (all alleles = '0'). If a genotype is missing or partial, treat that
 *   as “not guaranteed hom-ref,” so we keep that variant.
 */
class VCFXNonRefFilter {
public:
    // main entry point
    int run(int argc, char* argv[]);

private:
    // prints usage
    void displayHelp();

    // checks the input line-by-line, printing lines that pass the filter
    void filterNonRef(std::istream& in, std::ostream& out);

    // parse a sample’s genotype from a field => returns true if definitely hom-ref
    bool isDefinitelyHomRef(const std::string &genotypeField) const;
};

#endif // VCFX_NONREF_FILTER_H


./VCFX_outlier_detector/VCFX_outlier_detector.cpp
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
    if (genotypeField.empty() || genotypeField == ".") {
        return false;
    }

    // First, try the equals-sign based approach for non-standard formats
    std::stringstream ss1(genotypeField);
    std::string token;
    while(std::getline(ss1, token, ':')) {
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

    // If that didn't work, check if we're looking for a standard field in the FORMAT column
    // Need to find the corresponding FORMAT field and extract the value at that position
    std::vector<std::string> fieldValues;
    std::stringstream ss2(genotypeField);
    while(std::getline(ss2, token, ':')) {
        fieldValues.push_back(token);
    }

    // We need to check the FORMAT column, which we don't have here
    // Return false and handle FORMAT fields at the variant level in detectOutliers
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
            
            // Parse the FORMAT field to find the index of our metric
            int metricIndex = -1;
            std::vector<std::string> formatFields;
            {
                std::stringstream fmt(format);
                std::string token;
                int idx = 0;
                while(std::getline(fmt, token, ':')) {
                    if(token == metric) {
                        metricIndex = idx;
                    }
                    formatFields.push_back(token);
                    idx++;
                }
            }
            
            // Get sample columns
            std::vector<std::string> sampleColumns;
            for(size_t i=9;i<fields.size();i++){
                sampleColumns.push_back(fields[i]);
            }
            
            // Process each sample
            for(size_t s=0;s<sampleNames.size() && s<sampleColumns.size();s++){
                double val=0.0;
                // First, try parsing with the custom equals-sign based approach
                if(parseMetricFromGenotype(sampleColumns[s], metric, val)) {
                    sums[sampleNames[s]] += val;
                    counts[sampleNames[s]]++;
                    anyMetricFound = true;
                }
                // If that didn't work and we found the metric in FORMAT, try standard approach
                else if(metricIndex >= 0) {
                    std::vector<std::string> sampleValues;
                    std::stringstream sampleStream(sampleColumns[s]);
                    std::string sampleToken;
                    while(std::getline(sampleStream, sampleToken, ':')) {
                        sampleValues.push_back(sampleToken);
                    }
                    
                    if(metricIndex < (int)sampleValues.size() && sampleValues[metricIndex] != "." && !sampleValues[metricIndex].empty()) {
                        try {
                            val = std::stod(sampleValues[metricIndex]);
                            sums[sampleNames[s]] += val;
                            counts[sampleNames[s]]++;
                            anyMetricFound = true;
                        } catch(...) {
                            // Invalid numeric value, skip
                        }
                    }
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
                } else {
                    out<< nm<<"\tNA\n";
                }
            }
        }
        if(!anyMetricFound) {
            std::cerr<<"Warning: metric '"<<metric<<"' was not found in any sample genotype.\n";
        }
    }
}

int main(int argc, char* argv[]){
    VCFXOutlierDetector app;
    return app.run(argc, argv);
}

./VCFX_outlier_detector/VCFX_outlier_detector.h
#ifndef VCFX_OUTLIER_DETECTOR_H
#define VCFX_OUTLIER_DETECTOR_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

class VCFXOutlierDetector {
public:
    // Main entry point
    int run(int argc, char* argv[]);

private:
    // Print usage
    void displayHelp();

    // The function that does the reading and analysis
    void detectOutliers(std::istream &in, std::ostream &out,
                        const std::string &metric,
                        double threshold,
                        bool isVariantMode);

    // Parse the user-specified metric from INFO
    bool parseMetricFromInfo(const std::string &info,
                             const std::string &key,
                             double &val) const;

    // Parse the user-specified metric from genotype subfields
    bool parseMetricFromGenotype(const std::string &genotypeField,
                                 const std::string &metric,
                                 double &value) const;
};

#endif


./VCFX_phase_checker/VCFX_phase_checker.cpp
#include "VCFX_phase_checker.h"
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <sys/select.h>
#include <unistd.h>

int VCFXPhaseChecker::run(int argc, char* argv[]) {
    bool showHelp = false;
    
    // Check if stdin has data available
    struct timeval tv;
    fd_set fds;
    tv.tv_sec = 0;
    tv.tv_usec = 0;
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds);
    bool hasStdinInput = select(STDIN_FILENO+1, &fds, NULL, NULL, &tv) > 0;
    
    // If no arguments are provided and no stdin input, display help.
    if (argc == 1 && !hasStdinInput) {
        displayHelp();
        return 0;
    }
    
    int opt;
    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    while ((opt = getopt(argc, argv, "h")) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            default:
                showHelp = true;
                break;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }
    processVCF(std::cin, std::cout);
    return 0;
}

void VCFXPhaseChecker::displayHelp() {
    std::cout <<
"VCFX_phase_checker: Output only VCF variant lines in which every sample genotype is fully phased.\n\n"
"Usage:\n"
"  ./VCFX_phase_checker [options] < input.vcf > phased_output.vcf\n\n"
"Options:\n"
"  -h, --help   Display this help message and exit\n\n"
"Description:\n"
"  The tool reads a VCF from standard input and checks the GT field (genotype) for each sample.\n"
"  A genotype is considered fully phased if it uses the '|' separator (e.g., \"0|1\") and contains\n"
"  no missing alleles. If every sample in a variant line is fully phased, the line is printed to\n"
"  standard output; otherwise, it is skipped and a warning is written to standard error.\n\n"
"Examples:\n"
"  To extract only fully phased variants:\n"
"    ./VCFX_phase_checker < input.vcf > phased_output.vcf\n";
}

bool VCFXPhaseChecker::isFullyPhased(const std::string &gt) const {
    if(gt.empty() || gt == "." || gt == "./." || gt == ".|.") return false;
    
    // Haploid genotypes (like "0" or "1") are not considered phased
    if(gt.find('|') == std::string::npos) return false;
    
    // A fully phased genotype should use '|' exclusively.
    if(gt.find('/') != std::string::npos) return false;
    
    std::vector<std::string> alleles;
    std::istringstream ss(gt);
    std::string token;
    while(std::getline(ss, token, '|')) {
        alleles.push_back(token);
    }
    if(alleles.empty()) return false;
    for(const auto &al : alleles) {
        if(al.empty() || al == ".") return false;
    }
    return true;
}

void VCFXPhaseChecker::processVCF(std::istream &in, std::ostream &out) {
    bool headerFound = false;
    std::string line;
    while (std::getline(in, line)) {
        if(line.empty()) {
            out << line << "\n";
            continue;
        }
        if(line[0] == '#') {
            out << line << "\n";
            if(line.rfind("#CHROM", 0) == 0) {
                headerFound = true;
            }
            continue;
        }
        if(!headerFound) {
            std::cerr << "Warning: Data line encountered before #CHROM header; skipping line.\n";
            continue;
        }
        std::vector<std::string> fields;
        std::istringstream ss(line);
        std::string f;
        while(std::getline(ss, f, '\t')) {
            fields.push_back(f);
        }
        if(fields.size() < 10) {
            std::cerr << "Warning: Invalid VCF line with fewer than 10 columns; skipping line.\n";
            continue;
        }
        std::vector<std::string> fmt;
        {
            std::istringstream fs(fields[8]);
            std::string token;
            while(std::getline(fs, token, ':')) {
                fmt.push_back(token);
            }
        }
        int gtIndex = -1;
        for (size_t i = 0; i < fmt.size(); ++i) {
            if (fmt[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }
        if (gtIndex == -1) {
            std::cerr << "Warning: GT field not found; skipping line.\n";
            continue;
        }
        bool allPhased = true;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::vector<std::string> sampleFields;
            std::istringstream ssSample(fields[i]);
            std::string sub;
            while(std::getline(ssSample, sub, ':')) {
                sampleFields.push_back(sub);
            }
            if (gtIndex >= static_cast<int>(sampleFields.size())) {
                std::cerr << "Warning: GT index out of range in sample; skipping line.\n";
                allPhased = false;
                break;
            }
            if (!isFullyPhased(sampleFields[gtIndex])) {
                allPhased = false;
                break;
            }
        }
        if (allPhased) {
            out << line << "\n";
        } else {
            std::cerr << "Unphased genotype found at CHROM=" << fields[0] 
                      << ", POS=" << fields[1] << "; line skipped.\n";
        }
    }
}

int main(int argc, char* argv[]){
    VCFXPhaseChecker checker;
    return checker.run(argc, argv);
}


./VCFX_phase_checker/VCFX_phase_checker.h
#ifndef VCFX_PHASE_CHECKER_H
#define VCFX_PHASE_CHECKER_H

#include <iostream>
#include <string>
#include <vector>

/*
 * VCFXPhaseChecker:
 *   Reads a VCF, outputs only lines in which all samples' GT are fully phased.
 *   If any sample is unphased (or missing GT), logs a warning and discards line.
 */
class VCFXPhaseChecker {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    // Processes the input line-by-line. 
    // For each data line, if all samples are phased => write to stdout, else skip.
    void processVCF(std::istream &in, std::ostream &out);

    // Checks if a genotype string is "fully phased" (no '/' seen, multiple alleles separated by '|').
    // If GT is missing or partial => returns false (unphased).
    bool isFullyPhased(const std::string &gt) const;
};

#endif // VCFX_PHASE_CHECKER_H


./VCFX_phase_quality_filter/VCFX_phase_quality_filter.cpp
#include "VCFX_phase_quality_filter.h"
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cctype>

static void split(const std::string &s, char delim, std::vector<std::string> &tokens) {
    tokens.clear();
    std::stringstream ss(s);
    std::string t;
    while(std::getline(ss,t,delim)) {
        tokens.push_back(t);
    }
}

int VCFXPhaseQualityFilter::run(int argc, char* argv[]) {
    if(argc==1) {
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    std::string condition;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"filter-pq", required_argument, 0, 'f'},
        {0,0,0,0}
    };
    while(true) {
        int c= ::getopt_long(argc, argv, "hf:", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'f':
                condition= optarg;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp) {
        displayHelp();
        return 0;
    }
    if(condition.empty()) {
        std::cerr<<"Error: Must specify condition with --filter-pq\n";
        displayHelp();
        return 1;
    }
    std::string op;
    double threshold=0.0;
    if(!parseCondition(condition, op, threshold)) {
        std::cerr<<"Error: Unable to parse condition '"<< condition<<"'. e.g. PQ>=30\n";
        return 1;
    }
    filterByPQ(std::cin, std::cout, op, threshold);
    return 0;
}

void VCFXPhaseQualityFilter::displayHelp() {
    std::cout <<
"VCFX_phase_quality_filter: Filter variants by phasing quality (PQ) in the INFO field.\n\n"
"Usage:\n"
"  VCFX_phase_quality_filter --filter-pq \"PQ<OP><THRESHOLD>\" < input.vcf > output.vcf\n\n"
"Options:\n"
"  -h, --help             Print this help message.\n"
"  -f, --filter-pq <COND> Condition like 'PQ>30', 'PQ>=20', 'PQ!=10', etc.\n\n"
"Description:\n"
"  Reads each variant line, extracts 'PQ=' from INFO. If missing or invalid, PQ=0.\n"
"  Keeps lines if 'PQ <OP> THRESHOLD' is true. Otherwise, discards.\n\n"
"Supported operators: >, >=, <, <=, ==, !=\n\n"
"Examples:\n"
"  1) Keep variants with PQ>30:\n"
"     VCFX_phase_quality_filter --filter-pq \"PQ>30\" < in.vcf > out.vcf\n"
"  2) Keep PQ<=15:\n"
"     VCFX_phase_quality_filter --filter-pq \"PQ<=15\" < in.vcf > out.vcf\n";
}

bool VCFXPhaseQualityFilter::parseCondition(const std::string &condition,
                                            std::string &op,
                                            double &threshold) {
    if(condition.size()<3) return false;
    if(condition.rfind("PQ",0)!=0) {
        return false;
    }
    auto sub= condition.substr(2); 
    std::vector<std::string> ops={">=","<=","==","!=",">","<"};
    std::string matched;
    size_t foundPos=std::string::npos;
    for(auto &o : ops){
        auto p= sub.find(o);
        if(p==0) { matched=o; break; }
    }
    if(matched.empty()) {
        return false;
    }
    op= matched;
    auto valStr= sub.substr(op.size());
    if(valStr.empty()) return false;
    try {
        threshold= std::stod(valStr);
    } catch(...){
        return false;
    }
    return true;
}

void VCFXPhaseQualityFilter::filterByPQ(std::istream &in, std::ostream &out,
                                        const std::string &op,
                                        double threshold) {
    bool headerFound=false;
    std::string line;
    while(true) {
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line <<"\n";
            if(line.rfind("#CHROM",0)==0) headerFound=true;
            continue;
        }
        if(!headerFound){
            std::cerr<<"Warning: data line before #CHROM => skipping\n";
            continue;
        }
        std::vector<std::string> fields;
        {
            std::stringstream ss(line);
            std::string f;
            while(std::getline(ss,f,'\t')) fields.push_back(f);
        }
        if(fields.size()<8){
            std::cerr<<"Warning: VCF line with fewer than 8 columns => skipping.\n";
            continue;
        }
        // parse info
        double pq= parsePQScore(fields[7]);
        bool keep=false;
        if(op==">") {
            if(pq> threshold) keep=true;
        } else if(op==">=") {
            if(pq>=threshold) keep=true;
        } else if(op=="<") {
            if(pq< threshold) keep=true;
        } else if(op=="<=") {
            if(pq<=threshold) keep=true;
        } else if(op=="==") {
            if(pq==threshold) keep=true;
        } else if(op=="!=") {
            if(pq!=threshold) keep=true;
        }
        if(keep) {
            out<< line <<"\n";
        }
    }
}

double VCFXPhaseQualityFilter::parsePQScore(const std::string &info) {
    if(info.empty() || info=="." ) return 0.0;
    std::stringstream ss(info);
    std::string kv;
    while(std::getline(ss, kv, ';')) {
        if(kv.rfind("PQ=",0)==0) {
            auto valStr= kv.substr(3);
            try {
                return std::stod(valStr);
            } catch(...) {
                std::cerr<<"Warning: invalid PQ= '"<<valStr<<"'. Using 0.0.\n";
                return 0.0;
            }
        }
    }
    return 0.0;
}

int main(int argc, char* argv[]){
    VCFXPhaseQualityFilter f;
    return f.run(argc, argv);
}


./VCFX_phase_quality_filter/VCFX_phase_quality_filter.h
#ifndef VCFX_PHASE_QUALITY_FILTER_H
#define VCFX_PHASE_QUALITY_FILTER_H

#include <iostream>
#include <string>

class VCFXPhaseQualityFilter {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();
    // Takes the final operator (>, >=, <, <=, ==, !=) plus threshold, filters lines
    void filterByPQ(std::istream &in, std::ostream &out,
                    const std::string &op, double threshold);

    // Extracts the "PQ=" value from INFO or returns 0.0 if missing or invalid
    double parsePQScore(const std::string &info);

    // Internal helper to parse the condition string, e.g. "PQ>=30" => op=">=", thr=30
    bool parseCondition(const std::string &condition,
                        std::string &op, double &threshold);
};

#endif // VCFX_PHASE_QUALITY_FILTER_H


./VCFX_phred_filter/VCFX_phred_filter.cpp
#include "VCFX_phred_filter.h"
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>

int VCFXPhredFilter::run(int argc, char* argv[]){
    if(argc==1) {
        displayHelp();
        return 0;
    }
    double threshold= 30.0;
    bool showHelp=false;
    bool keepMissingAsPass= false;
    static struct option long_opts[]={
        {"phred-filter", required_argument, 0, 'p'},
        {"keep-missing-qual", no_argument, 0, 'k'},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "p:kh", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'p': {
                try { threshold= std::stod(optarg); }
                catch(...){
                    std::cerr<<"Error: Invalid threshold '"<< optarg<<"'.\n";
                    return 1;
                }
            }break;
            case 'k':
                keepMissingAsPass= true;
                break;
            case 'h':
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    processVCF(std::cin, threshold, keepMissingAsPass);
    return 0;
}

void VCFXPhredFilter::displayHelp(){
    std::cout <<
"VCFX_phred_filter: Filter VCF lines by their QUAL field.\n\n"
"Usage:\n"
"  VCFX_phred_filter [options] < input.vcf > output.vcf\n\n"
"Options:\n"
"  -p, --phred-filter <VAL>      Phred QUAL threshold (default=30)\n"
"  -k, --keep-missing-qual       Treat '.' (missing QUAL) as pass\n"
"  -h, --help                    Display this help and exit\n\n"
"Description:\n"
"  Reads VCF lines from stdin. For each data line, parse the QUAL field.\n"
"  If QUAL >= threshold => print line. Otherwise, skip. By default, missing\n"
"  QUAL ('.') is treated as 0. Use --keep-missing-qual to treat '.' as pass.\n\n"
"Examples:\n"
"  1) Keep variants with QUAL>=30:\n"
"     VCFX_phred_filter -p 30 < in.vcf > out.vcf\n"
"  2) Keep missing QUAL lines:\n"
"     VCFX_phred_filter -p 30 --keep-missing-qual < in.vcf > out.vcf\n";
}

void VCFXPhredFilter::processVCF(std::istream &in, double threshold, bool keepMissingAsPass){
    std::string line;
    bool foundChrom=false;
    while(true){
        if(!std::getline(in, line)) break;
        if(line.empty()){
            std::cout<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            std::cout<< line <<"\n";
            if(line.rfind("#CHROM",0)==0) foundChrom=true;
            continue;
        }
        if(!foundChrom) {
            std::cerr<<"Warning: data line before #CHROM => skipping line.\n";
            continue;
        }
        std::vector<std::string> fields;
        {
            std::stringstream ss(line);
            std::string f;
            while(std::getline(ss,f,'\t')) {
                fields.push_back(f);
            }
        }
        // we need at least CHROM,POS,ID,REF,ALT,QUAL => 6 columns
        if(fields.size()<6) {
            std::cerr<<"Warning: line has <6 columns => skipping.\n";
            continue;
        }
        double q= parseQUAL(fields[5], keepMissingAsPass);
        if(q>= threshold) {
            std::cout<< line <<"\n";
        }
    }
}

double VCFXPhredFilter::parseQUAL(const std::string &qualStr, bool keepMissingAsPass){
    if(qualStr=="." || qualStr.empty()){
        if(keepMissingAsPass) return 1e9; 
        else return 0.0;
    }
    try {
        return std::stod(qualStr);
    } catch(...) {
        std::cerr<<"Warning: Invalid QUAL '"<<qualStr<<"'. Using 0.\n";
        return 0.0;
    }
}

int main(int argc, char* argv[]){
    VCFXPhredFilter pf;
    return pf.run(argc,argv);
}


./VCFX_phred_filter/VCFX_phred_filter.h
#ifndef VCFX_PHRED_FILTER_H
#define VCFX_PHRED_FILTER_H

#include <iostream>
#include <string>

class VCFXPhredFilter {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();
    void processVCF(std::istream &in, double threshold, bool keepMissingAsPass);
    double parseQUAL(const std::string &qualStr, bool keepMissingAsPass);
};

#endif


./VCFX_population_filter/VCFX_population_filter.cpp
#include "VCFX_population_filter.h"
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <cstdlib>

int VCFXPopulationFilter::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    std::string populationTag;
    std::string popMapFile;

    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"population", required_argument, 0, 'p'},
        {"pop-map", required_argument, 0, 'm'},
        {0,0,0,0}
    };

    while(true){
        int c= ::getopt_long(argc, argv, "hp:m:", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'p':
                populationTag= optarg;
                break;
            case 'm':
                popMapFile= optarg;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    if(populationTag.empty()|| popMapFile.empty()){
        std::cerr<<"Error: Must specify --population <TAG> and --pop-map <file>.\n";
        displayHelp();
        return 1;
    }

    std::unordered_set<std::string> samplesToInclude;
    if(!loadPopulationMap(popMapFile, populationTag, samplesToInclude)){
        std::cerr<<"Error: Unable to load or parse pop map.\n";
        return 1;
    }
    if(samplesToInclude.empty()){
        std::cerr<<"Warning: No samples found for population tag: "<< populationTag<<"\n";
    }
    filterPopulation(std::cin, std::cout, samplesToInclude, populationTag);
    return 0;
}

void VCFXPopulationFilter::displayHelp(){
    std::cout <<
"VCFX_population_filter: Subset VCF to samples in specified population.\n\n"
"Usage:\n"
"  VCFX_population_filter [options] < input.vcf > output.vcf\n\n"
"Options:\n"
"  --help, -h               Print this help.\n"
"  --population, -p <TAG>   Population tag to keep (e.g. 'EUR','AFR', etc.)\n"
"  --pop-map, -m <FILE>     Tab-delimited file: 'SampleName <tab> Population'\n\n"
"Description:\n"
"  Reads the pop map, finds samples that match the chosen population.\n"
"  Then reads the VCF from stdin and prints lines with only those sample columns.\n"
"  If a sample is not in that population, it's dropped from the #CHROM header and data columns.\n\n"
"Example:\n"
"  VCFX_population_filter --population AFR --pop-map pops.txt < input.vcf > out.vcf\n";
}

bool VCFXPopulationFilter::loadPopulationMap(const std::string &popMapFile,
                                             const std::string &popTag,
                                             std::unordered_set<std::string> &samplesToInclude)
{
    std::ifstream f(popMapFile);
    if(!f.is_open()){
        std::cerr<<"Error: cannot open "<< popMapFile<<"\n";
        return false;
    }
    std::string line;
    while(true){
        if(!std::getline(f,line)) break;
        if(line.empty())continue;
        std::stringstream ss(line);
        std::string samp, pop;
        if(!(ss>>samp>>pop)){
            std::cerr<<"Warning: popmap line invalid: "<< line<<"\n";
            continue;
        }
        if(pop== popTag){
            samplesToInclude.insert(samp);
        }
    }
    return true;
}

void VCFXPopulationFilter::filterPopulation(std::istream &in,
                                            std::ostream &out,
                                            const std::unordered_set<std::string> &samplesToInclude,
                                            const std::string &popTag)
{
    bool foundChromLine=false;
    std::string line;
    std::vector<std::string> finalHeader;
    std::vector<int> sampleIndices;

    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            if(line.rfind("#CHROM",0)==0){
                foundChromLine=true;
                std::stringstream ss(line);
                std::vector<std::string> fields;
                {
                    std::string col;
                    while(std::getline(ss,col,'\t')){
                        fields.push_back(col);
                    }
                }
                // fields[0..8] => fixed, fields[9..] => samples
                sampleIndices.clear();
                for(size_t i=9; i<fields.size(); i++){
                    if(samplesToInclude.count(fields[i])>0){
                        sampleIndices.push_back((int)i);
                    }
                }
                // build new #CHROM line
                std::ostringstream newChrom;
                for(int i=0; i<9; i++){
                    newChrom<< fields[i];
                    if(i<8 || !sampleIndices.empty()) newChrom<<"\t";
                }
                for(size_t i=0;i< sampleIndices.size(); i++){
                    newChrom << fields[sampleIndices[i]];
                    if(i+1< sampleIndices.size()) newChrom<<"\t";
                }
                out<< newChrom.str()<<"\n";
            } else {
                out<< line <<"\n";
            }
            continue;
        }
        if(!foundChromLine) {
            std::cerr<<"Warning: data line before #CHROM => skipping.\n";
            continue;
        }
        std::vector<std::string> columns;
        {
            std::stringstream ss(line);
            std::string x;
            while(std::getline(ss,x,'\t')) {
                columns.push_back(x);
            }
        }
        if(columns.size()<9) {
            std::cerr<<"Warning: line with fewer than 9 columns => skipping.\n";
            continue;
        }
        // build new line
        std::ostringstream newLine;
        for(int i=0; i<9; i++){
            newLine<< columns[i];
            if(i<8 || !sampleIndices.empty()) newLine<<"\t";
        }
        for(size_t i=0;i< sampleIndices.size(); i++){
            newLine<< columns[sampleIndices[i]];
            if(i+1< sampleIndices.size()) newLine<<"\t";
        }
        out<< newLine.str()<<"\n";
    }
    if(!foundChromLine) {
        std::cerr<<"Error: No #CHROM header found in VCF.\n";
    }
}

int main(int argc, char* argv[]){
    VCFXPopulationFilter pf;
    return pf.run(argc, argv);
}


./VCFX_population_filter/VCFX_population_filter.h
#ifndef VCFX_POPULATION_FILTER_H
#define VCFX_POPULATION_FILTER_H

#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

class VCFXPopulationFilter {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    bool loadPopulationMap(const std::string &popMapFile,
                           const std::string &popTag,
                           std::unordered_set<std::string> &samplesToInclude);
    void filterPopulation(std::istream &in,
                          std::ostream &out,
                          const std::unordered_set<std::string> &samplesToInclude,
                          const std::string &popTag);
};

#endif


./VCFX_position_subsetter/VCFX_position_subsetter.cpp
#include "VCFX_position_subsetter.h"
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>

int VCFXPositionSubsetter::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    std::string regionStr;

    static struct option long_opts[] = {
        {"region", required_argument, 0, 'r'},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };

    while(true){
        int c= ::getopt_long(argc, argv, "r:h", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'r':
                regionStr= optarg;
                break;
            case 'h':
            default:
                showHelp=true;
        }
    }

    if(showHelp){
        displayHelp();
        return 0;
    }
    if(regionStr.empty()){
        std::cerr<<"Error: --region <chrX:start-end> is required.\n";
        displayHelp();
        return 1;
    }
    std::string chrom;
    int start=0, end=0;
    if(!parseRegion(regionStr, chrom, start, end)){
        return 1;
    }
    bool ok= subsetVCFByPosition(std::cin, std::cout, chrom, start, end);
    return ok? 0: 1;
}

void VCFXPositionSubsetter::displayHelp(){
    std::cout <<
"VCFX_position_subsetter: Subset VCF by a single genomic region.\n\n"
"Usage:\n"
"  VCFX_position_subsetter --region \"chr1:10000-20000\" < in.vcf > out.vcf\n\n"
"Options:\n"
"  -r, --region \"CHR:START-END\"   The region to keep.\n"
"  -h, --help                     Print this help.\n\n"
"Description:\n"
"  Reads lines from VCF input, and only prints data lines where:\n"
"    1) CHROM matches 'CHR' exactly, and\n"
"    2) POS is in [START,END].\n"
"  All header lines (#...) are passed unmodified.\n\n"
"Example:\n"
"  VCFX_position_subsetter --region \"chr2:500-1000\" < input.vcf > subset.vcf\n";
}

bool VCFXPositionSubsetter::parseRegion(const std::string &regionStr,
                                        std::string &chrom,
                                        int &start,
                                        int &end){
    // find colon, dash
    auto cpos= regionStr.find(':');
    auto dpos= regionStr.find('-');
    if(cpos== std::string::npos || dpos==std::string::npos || dpos<= cpos){
        std::cerr<<"Error: invalid region: "<< regionStr<<". Expected e.g. chr1:10000-20000.\n";
        return false;
    }
    chrom= regionStr.substr(0, cpos);
    std::string startStr= regionStr.substr(cpos+1, dpos-(cpos+1));
    std::string endStr= regionStr.substr(dpos+1);
    try {
        start= std::stoi(startStr);
        end= std::stoi(endStr);
    } catch(...){
        std::cerr<<"Error: cannot parse region start/end.\n";
        return false;
    }
    if(start> end){
        std::cerr<<"Error: region start> end.\n";
        return false;
    }
    return true;
}

bool VCFXPositionSubsetter::subsetVCFByPosition(std::istream &in,
                                               std::ostream &out,
                                               const std::string &regionChrom,
                                               int regionStart,
                                               int regionEnd)
{
    bool headerFound=false;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<<line<<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line <<"\n";
            if(line.rfind("#CHROM",0)==0) headerFound=true;
            continue;
        }
        if(!headerFound){
            std::cerr<<"Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }
        std::vector<std::string> fields;
        {
            std::stringstream ss(line);
            std::string f;
            while(std::getline(ss,f,'\t')) fields.push_back(f);
        }
        if(fields.size()<2){
            std::cerr<<"Warning: line has <2 columns => skipping.\n";
            continue;
        }
        const std::string &chrom= fields[0];
        const std::string &posStr= fields[1];
        int pos=0;
        try {
            pos= std::stoi(posStr);
        } catch(...){
            std::cerr<<"Warning: invalid POS '"<< posStr<<"'. Skipping.\n";
            continue;
        }
        if(chrom== regionChrom && pos>= regionStart && pos<= regionEnd){
            out<< line <<"\n";
        }
    }
    return true;
}

int main(int argc, char* argv[]){
    VCFXPositionSubsetter subsetter;
    return subsetter.run(argc, argv);
}


./VCFX_position_subsetter/VCFX_position_subsetter.h
#ifndef VCFX_POSITION_SUBSETTER_H
#define VCFX_POSITION_SUBSETTER_H

#include <iostream>
#include <string>

class VCFXPositionSubsetter {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    // parse "chr1:10000-20000", store in regionChrom, regionStart, regionEnd
    bool parseRegion(const std::string &regionStr,
                     std::string &chrom,
                     int &start, int &end);

    // do the actual subsetting from in->out
    bool subsetVCFByPosition(std::istream &in, std::ostream &out,
                             const std::string &regionChrom,
                             int regionStart,
                             int regionEnd);
};

#endif // VCFX_POSITION_SUBSETTER_H


./VCFX_probability_filter/VCFX_probability_filter.cpp
#include "VCFX_probability_filter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <regex>

// Implementation of VCFXProbabilityFilter
int VCFXProbabilityFilter::run(int argc, char* argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string condition;

    static struct option long_options[] = {
        {"help",           no_argument,       0, 'h'},
        {"filter-probability", required_argument, 0, 'f'},
        {0,                0,                 0,  0 }
    };

    while ((opt = getopt_long(argc, argv, "hf:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'h':
                showHelp = true;
                break;
            case 'f':
                condition = optarg;
                break;
            default:
                showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }
    
    if (condition.empty()) {
        displayHelp();
        return 1;
    }

    // Perform probability filtering on stdin and output to stdout
    filterByProbability(std::cin, std::cout, condition);

    return 0;
}

void VCFXProbabilityFilter::displayHelp() {
    std::cout << "VCFX_probability_filter: Filter VCF based on genotype probability scores.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_probability_filter --filter-probability \"<CONDITION>\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                        Display this help message and exit\n";
    std::cout << "  -f, --filter-probability <cond>    Specify the genotype probability filter condition (e.g., GP>0.9)\n\n";
    std::cout << "Supported Operators: >, <, >=, <=, ==, !=\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_probability_filter --filter-probability \"GP>0.9\" < input.vcf > filtered.vcf\n";
}

void VCFXProbabilityFilter::filterByProbability(std::istream& in, std::ostream& out, const std::string& condition) {
    // Parse the filter condition using regex (e.g., "GP>0.9")
    std::regex conditionRegex(R"((\w+)\s*(>=|<=|>|<|==|!=)\s*([0-9]*\.?[0-9]+))");
    std::smatch matches;
    if (!std::regex_match(condition, matches, conditionRegex)) {
        std::cerr << "Error: Invalid filter condition format. Expected format like \"GP>0.9\".\n";
        return;
    }

    std::string field = matches[1];
    std::string op = matches[2];
    double threshold = std::stod(matches[3]);

    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    size_t formatIndex = std::string::npos;
    size_t fieldIndex = std::string::npos;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Handle header lines
            if (line.substr(0, 6) == "#CHROM") {
                // Parse header to identify the FORMAT column and sample columns
                out << line << "\n";
                headerParsed = true;
            } else {
                // Other header lines
                out << line << "\n";
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse VCF data lines
        std::stringstream ss(line);
        std::string fieldEntry;
        std::vector<std::string> fieldsVec;

        while (std::getline(ss, fieldEntry, '\t')) {
            fieldsVec.push_back(fieldEntry);
        }

        if (fieldsVec.size() < 9) {
            std::cerr << "Warning: Invalid VCF line with fewer than 9 fields: " << line << "\n";
            continue;
        }

        std::string formatField = fieldsVec[8];
        std::vector<std::string> formatFields;
        std::stringstream fmt_ss(formatField);
        std::string fmt_field;
        while (std::getline(fmt_ss, fmt_field, ':')) {
            formatFields.push_back(fmt_field);
        }

        // Find the index of the specified field in the FORMAT column
        if (fieldIndex == std::string::npos) {
            for (size_t i = 0; i < formatFields.size(); ++i) {
                if (formatFields[i] == field) {
                    fieldIndex = i;
                    break;
                }
            }

            if (fieldIndex == std::string::npos) {
                std::cerr << "Error: Specified field \"" << field << "\" not found in FORMAT column.\n";
                return;
            }
        }

        bool pass = true;

        // Iterate over each sample to check the genotype probability
        for (size_t i = 9; i < fieldsVec.size(); ++i) {
            std::string sample = fieldsVec[i];
            std::stringstream samp_ss(sample);
            std::string samp_field;
            std::vector<std::string> sampleFields;

            while (std::getline(samp_ss, samp_field, ':')) {
                sampleFields.push_back(samp_field);
            }

            if (fieldIndex >= sampleFields.size()) {
                std::cerr << "Warning: Field index out of range in sample fields.\n";
                pass = false;
                break;
            }

            std::string valueStr = sampleFields[fieldIndex];
            if (valueStr.empty() || valueStr == ".") {
                pass = false;
                break;
            }

            double value;
            try {
                value = std::stod(valueStr);
            } catch (...) {
                std::cerr << "Warning: Unable to convert value \"" << valueStr << "\" to number.\n";
                pass = false;
                break;
            }

            // Apply the filter condition
            if (op == ">") {
                if (!(value > threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "<") {
                if (!(value < threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == ">=") {
                if (!(value >= threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "<=") {
                if (!(value <= threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "==") {
                if (!(value == threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "!=") {
                if (!(value != threshold)) {
                    pass = false;
                    break;
                }
            }
        }

        if (pass) {
            out << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    VCFXProbabilityFilter probabilityFilter;
    return probabilityFilter.run(argc, argv);
}

./VCFX_probability_filter/VCFX_probability_filter.h
#ifndef VCFX_PROBABILITY_FILTER_H
#define VCFX_PROBABILITY_FILTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXProbabilityFilter: Header file for Genotype Probability Filter tool
class VCFXProbabilityFilter {
public:
    // Entry point for the tool
    int run(int argc, char* argv[]);

private:
    // Displays the help message
    void displayHelp();

    // Filters VCF input based on the specified genotype probability condition
    void filterByProbability(std::istream& in, std::ostream& out, const std::string& condition);
};

#endif // VCFX_PROBABILITY_FILTER_H


./VCFX_quality_adjuster/VCFX_quality_adjuster.cpp
#include "VCFX_quality_adjuster.h"
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdlib>

int VCFXQualityAdjuster::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    bool clamp=true;
    std::string transformStr;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"adjust-qual", required_argument, 0, 'a'},
        {"no-clamp", no_argument, 0, 'n'},
        {0,0,0,0}
    };
    while(true){
        int c=::getopt_long(argc, argv, "ha:n", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'a':
                transformStr= optarg;
                break;
            case 'n':
                clamp=false;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    if(transformStr.empty()){
        std::cerr<<"Error: Must specify a transformation with --adjust-qual <FUNC>.\n";
        displayHelp();
        return 1;
    }
    initSupportedFunctions();
    std::function<double(double)> transFunc;
    if(!parseTransformationFunction(transformStr, transFunc)){
        std::cerr<<"Error: unsupported transformation '"<< transformStr<<"'.\n";
        return 1;
    }
    adjustQualityScores(std::cin, std::cout, transFunc, clamp);
    return 0;
}

void VCFXQualityAdjuster::displayHelp(){
    std::cout <<
"VCFX_quality_adjuster: Apply a transformation to the QUAL field of a VCF.\n\n"
"Usage:\n"
"  VCFX_quality_adjuster [options] < input.vcf > output.vcf\n\n"
"Options:\n"
"  -h, --help               Show this help.\n"
"  -a, --adjust-qual <FUNC> Required. One of: log, sqrt, square, identity.\n"
"  -n, --no-clamp           Do not clamp negative or large values.\n\n"
"Description:\n"
"  Reads each line from VCF. If it's a data line with >=8 columns, we parse\n"
"  the QUAL field (6th col). We transform it with <FUNC>, e.g.:\n"
"    log => log(QUAL + 1e-10)\n"
"    sqrt=> sqrt(QUAL)\n"
"    square=> (QUAL * QUAL)\n"
"    identity=> no change\n"
"  By default, negative results from e.g. log are clamped to 0, and large\n"
"  results are capped at 1e12. If you do not want clamping, use --no-clamp.\n\n"
"Examples:\n"
"  1) Log-transform:\n"
"     VCFX_quality_adjuster --adjust-qual log < in.vcf > out.vcf\n"
"  2) Square, keep negative or big values as is:\n"
"     VCFX_quality_adjuster --adjust-qual square --no-clamp < in.vcf > out.vcf\n";
}

void VCFXQualityAdjuster::initSupportedFunctions(){
    supportedFunctions.clear();
    supportedFunctions["log"] = [](double x){
        // protect near zero with epsilon
        return std::log(x + 1e-10);
    };
    supportedFunctions["sqrt"]= [](double x){
        return std::sqrt(std::max(0.0, x)); // can't sqrt negative
    };
    supportedFunctions["square"]= [](double x){
        return x*x;
    };
    supportedFunctions["identity"]= [](double x){
        return x;
    };
}

bool VCFXQualityAdjuster::parseTransformationFunction(const std::string &funcStr,
                                                      std::function<double(double)> &transFunc)
{
    auto it= supportedFunctions.find(funcStr);
    if(it== supportedFunctions.end()){
        return false;
    }
    transFunc= it->second;
    return true;
}

void VCFXQualityAdjuster::adjustQualityScores(std::istream &in, std::ostream &out,
                                              std::function<double(double)> transFunc,
                                              bool clamp)
{
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line <<"\n";
            continue;
        }
        // parse fields
        std::vector<std::string> fields;
        {
            std::stringstream ss(line);
            std::string f;
            while(std::getline(ss,f,'\t')){
                fields.push_back(f);
            }
        }
        if(fields.size()<8){
            std::cerr<<"Warning: line with <8 fields => skipping.\n";
            continue;
        }
        double oldQual=0.0;
        bool valid=true;
        if(fields[5]=="." || fields[5].empty()){
            // treat missing as 0?
            oldQual=0.0;
        } else {
            try {
                oldQual= std::stod(fields[5]);
            } catch(...){
                std::cerr<<"Warning: invalid QUAL '"<<fields[5]<<"'. Skipping.\n";
                valid=false;
            }
        }
        if(!valid) continue;
        double newQual= transFunc(oldQual);
        if(clamp){
            if(newQual< 0.0) newQual= 0.0;
            // clamp large values
            if(newQual>1e12) newQual= 1e12;
        }
        fields[5]= std::to_string(newQual);
        std::ostringstream oss;
        for(size_t i=0; i<fields.size(); i++){
            if(i>0) oss<<"\t";
            oss<< fields[i];
        }
        out<< oss.str()<<"\n";
    }
}

int main(int argc, char* argv[]){
    VCFXQualityAdjuster app;
    return app.run(argc, argv);
}

./VCFX_quality_adjuster/VCFX_quality_adjuster.h
#ifndef VCFX_QUALITY_ADJUSTER_H
#define VCFX_QUALITY_ADJUSTER_H

#include <iostream>
#include <string>
#include <functional>
#include <unordered_map>

class VCFXQualityAdjuster {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();
    // parse the transformation name and build a transform function
    bool parseTransformationFunction(const std::string &funcStr,
                                     std::function<double(double)> &transFunc);

    // read lines from 'in', transform the QUAL field (6th col),
    // write lines to 'out'. If clamp is false, do not clamp negative or large.
    void adjustQualityScores(std::istream &in, std::ostream &out,
                             std::function<double(double)> transFunc,
                             bool clamp);

    // map of supported transformations
    std::unordered_map<std::string, std::function<double(double)>> supportedFunctions;

    // helper to initialize the map
    void initSupportedFunctions();
};

#endif // VCFX_QUALITY_ADJUSTER_H


./VCFX_record_filter/VCFX_record_filter.cpp
#include "VCFX_record_filter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cctype>
#include <cstdlib>

// Trim leading/trailing whitespace
static std::string trim(const std::string &s){
    size_t start=0; while(start<s.size() && isspace((unsigned char)s[start]))start++;
    if(start==s.size()) return "";
    size_t end=s.size()-1; 
    while(end>start && isspace((unsigned char)s[end]))end--;
    return s.substr(start, end-start+1);
}

static void split(const std::string &s, char delim, std::vector<std::string> &out){
    out.clear();
    std::stringstream ss(s);
    std::string t;
    while(std::getline(ss,t,delim)){
        out.push_back(t);
    }
}

// parse an operator token: "==", "!=", ">=", "<=", ">", "<"
static bool parseOperator(const std::string &opStr, FilterOp &op){
    if(opStr==">") {op= FilterOp::GT;return true;}
    if(opStr==">="){op= FilterOp::GE;return true;}
    if(opStr=="<") {op= FilterOp::LT;return true;}
    if(opStr=="<="){op= FilterOp::LE;return true;}
    if(opStr=="=="){op= FilterOp::EQ;return true;}
    if(opStr=="!="){op= FilterOp::NE;return true;}
    return false;
}

// parse single token e.g. "POS>100" or "FILTER!=PASS" 
// returns false if invalid
static bool parseSingleCriterion(const std::string &token, FilterCriterion &crit){
    // find the operator among [>=,<=,==,!=,>,<]
    // a robust approach is to search for them in descending length order
    static std::vector<std::string> ops={">=", "<=", "==", "!=", ">", "<"};
    size_t bestPos= std::string::npos;
    std::string bestOp;
    for(auto &o: ops){
        auto p= token.find(o);
        if(p!=std::string::npos){
            bestPos= p; 
            bestOp= o;
            break; // first match
        }
    }
    if(bestPos==std::string::npos){
        std::cerr<<"Error: no operator found in '"<< token<<"'.\n";
        return false;
    }
    // parse fieldName= token.substr(0,bestPos)
    // parse value= token.substr(bestPos+ bestOp.size())
    crit.fieldName= trim(token.substr(0,bestPos));
    if(crit.fieldName.empty()){
        std::cerr<<"Error: empty field name in '"<< token<<"'.\n";
        return false;
    }
    if(!parseOperator(bestOp, crit.op)){
        std::cerr<<"Error: unknown operator '"<< bestOp<<"' in '"<< token<<"'.\n";
        return false;
    }
    std::string valPart= trim(token.substr(bestPos+ bestOp.size()));
    if(valPart.empty()){
        std::cerr<<"Error: no value in '"<< token<<"'.\n";
        return false;
    }
    // try parse numeric
    // if parse fails => treat as string
    try{
        double d= std::stod(valPart);
        // numeric
        crit.fieldType= FieldType::NUMERIC;
        crit.numericValue= d;
        crit.stringValue.clear();
    } catch(...){
        // treat as string
        crit.fieldType= FieldType::STRING;
        crit.stringValue= valPart;
        crit.numericValue= 0.0;
    }
    return true;
}

bool parseCriteria(const std::string &criteriaStr,
                   std::vector<FilterCriterion> &criteria)
{
    criteria.clear();
    std::vector<std::string> tokens;
    split(criteriaStr, ';', tokens);
    for(auto &t: tokens){
        auto trimmed= trim(t);
        if(trimmed.empty()) continue;
        FilterCriterion c;
        if(!parseSingleCriterion(trimmed, c)) {
            return false;
        }
        criteria.push_back(c);
    }
    if(criteria.empty()){
        std::cerr<<"Error: no valid criteria in '"<<criteriaStr<<"'.\n";
        return false;
    }
    return true;
}

// split by tab
static std::vector<std::string> tabSplit(const std::string &line){
    std::vector<std::string> f;
    std::stringstream ss(line);
    std::string x;
    while(std::getline(ss,x,'\t')) f.push_back(x);
    return f;
}

// parse info => either return key= or a flag
// if numeric, parse double
// if string => keep string
static bool getInfoValue(const std::string &infoField,
                         const std::string &key,
                         bool wantNumeric,
                         double &numericOut,
                         std::string &stringOut)
{
    if(infoField.empty()|| infoField=="." ) return false;
    std::vector<std::string> tokens;
    split(infoField, ';', tokens);
    for(auto &kv: tokens){
        auto eq= kv.find('=');
        if(eq!=std::string::npos){
            std::string k= trim(kv.substr(0,eq));
            std::string v= trim(kv.substr(eq+1));
            if(k== key){
                if(wantNumeric){
                    try{
                        double d= std::stod(v);
                        numericOut= d;
                        return true;
                    }catch(...){
                        return false;
                    }
                } else {
                    stringOut= v;
                    return true;
                }
            }
        } else {
            // it's a flag
            std::string fl= trim(kv);
            if(fl== key){
                if(wantNumeric){
                    numericOut= 1.0;
                } else {
                    stringOut= fl; 
                }
                return true;
            }
        }
    }
    return false;
}

static bool compareDouble(double x, FilterOp op, double y){
    switch(op){
        case FilterOp::GT: return (x> y);
        case FilterOp::GE: return (x>=y);
        case FilterOp::LT: return (x< y);
        case FilterOp::LE: return (x<=y);
        case FilterOp::EQ: return (x==y);
        case FilterOp::NE: return (x!=y);
    }
    return false; 
}

static bool compareString(const std::string &s, FilterOp op, const std::string &t){
    switch(op){
        case FilterOp::EQ: return (s==t);
        case FilterOp::NE: return (s!=t);
        // if user tries < or > => invalid for strings in this code => fail
        default: 
            // we do no partial compare. We fail this criterion
            return false;
    }
}

// We handle standard fields: POS(numeric), QUAL(numeric), FILTER(string?), plus info
// If fieldName== "POS", parse fields[1] => numeric
// If fieldName== "QUAL", parse fields[5] => numeric
// If fieldName== "FILTER", parse fields[6] => string
// else => treat as an INFO key
static bool evaluateCriterion(const std::vector<std::string> &fields,
                              const FilterCriterion &c){
    if(fields.size()<8) return false;
    if(c.fieldName=="POS"){
        if(fields[1].empty()) return false;
        try{
            double p= std::stod(fields[1]);
            return compareDouble(p, c.op, c.numericValue);
        } catch(...){
            return false;
        }
    } else if(c.fieldName=="QUAL"){
        // index 5
        if(fields[5].empty()|| fields[5]=="." ) {
            // treat missing as 0? or fail?
            double q=0.0;
            return compareDouble(q, c.op, c.numericValue);
        }
        try{
            double q= std::stod(fields[5]);
            return compareDouble(q, c.op, c.numericValue);
        }catch(...){
            return false;
        }
    } else if(c.fieldName=="FILTER"){
        // index=6 => string compare
        if(c.fieldType== FieldType::NUMERIC) {
            // invalid usage => fail
            return false;
        }
        return compareString(fields[6], c.op, c.stringValue);
    } else {
        // treat as info key
        // c.fieldType => numeric/string
        double num=0.0; 
        std::string st;
        bool got= getInfoValue(fields[7], c.fieldName, (c.fieldType==FieldType::NUMERIC), num, st);
        if(!got) return false; // can't find key => fail
        if(c.fieldType==FieldType::NUMERIC){
            return compareDouble(num, c.op, c.numericValue);
        } else {
            return compareString(st, c.op, c.stringValue);
        }
    }
}

bool recordPasses(const std::string &record,
                  const std::vector<FilterCriterion> &criteria,
                  bool useAndLogic)
{
    // parse fields
    auto fields= tabSplit(record);
    if(fields.size()<8) return false;
    bool anyPass= false;
    // check each criterion
    for(auto &crit: criteria){
        bool pass= evaluateCriterion(fields, crit);
        if(useAndLogic){
            if(!pass) return false; 
        } else {
            // or logic
            if(pass) anyPass= true;
        }
    }
    if(useAndLogic) {
        return true; // all pass
    } else {
        return anyPass; 
    }
}

void processVCF(std::istream &in,
                std::ostream &out,
                const std::vector<FilterCriterion> &criteria,
                bool useAndLogic)
{
    bool foundChrom=false;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line <<"\n"; 
            if(line.rfind("#CHROM",0)==0) foundChrom=true;
            continue;
        }
        if(!foundChrom){
            // data line before #CHROM => skip with a warning?
            std::cerr<<"Warning: data line before #CHROM => skipping.\n";
            continue;
        }
        // apply
        if(recordPasses(line, criteria, useAndLogic)){
            out<< line <<"\n";
        }
    }
}

void printHelp(){
    std::cout <<
"VCFX_record_filter: Filter VCF data lines by multiple criteria.\n\n"
"Usage:\n"
"  VCFX_record_filter [options] --filter \"CRITERIA\"\n"
"  < input.vcf > output.vcf\n\n"
"Options:\n"
"  --filter, -f \"...\"   One or more criteria separated by semicolons, e.g.\n"
"                        \"POS>10000; QUAL>=30; AF<0.05; FILTER==PASS\"\n"
"                        Each criterion must use an operator among >,>=,<,<=,==,!=\n\n"
"  --logic and|or        'and' => a line must pass all criteria (default)\n"
"                        'or'  => pass if any criterion is satisfied.\n"
"  --help, -h            Show this help.\n\n"
"Fields:\n"
"  POS => numeric, QUAL => numeric, FILTER => string.\n"
"  Others => assumed to be an INFO key. We try numeric parse if the criterion is numeric, else string.\n\n"
"Example:\n"
"  VCFX_record_filter --filter \"POS>=1000;FILTER==PASS;DP>10\" --logic and < in.vcf > out.vcf\n";
}

// main with typical argument parse
int main(int argc, char* argv[]){
    if(argc==1){
        printHelp();
        return 0;
    }
    bool showHelp=false;
    std::string criteriaStr;
    std::string logicStr="and";
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"filter", required_argument, 0, 'f'},
        {"logic", required_argument, 0, 'l'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hf:l:", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'f':
                criteriaStr= optarg;
                break;
            case 'l':
                logicStr= optarg;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        printHelp();
        return 0;
    }
    if(criteriaStr.empty()){
        std::cerr<<"Error: must provide --filter \"CRITERIA\".\n";
        printHelp();
        return 1;
    }
    bool useAndLogic=true;
    if(logicStr=="and"){
        useAndLogic=true;
    } else if(logicStr=="or"){
        useAndLogic=false;
    } else {
        std::cerr<<"Error: logic must be 'and' or 'or'.\n";
        return 1;
    }
    std::vector<FilterCriterion> criteria;
    if(!parseCriteria(criteriaStr, criteria)){
        std::cerr<<"Error: failed to parse criteria.\n";
        return 1;
    }
    processVCF(std::cin, std::cout, criteria, useAndLogic);
    return 0;
}


./VCFX_record_filter/VCFX_record_filter.h
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


./VCFX_ref_comparator/VCFX_ref_comparator.cpp
#include "VCFX_ref_comparator.h"
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cctype>

// A small helper to remove whitespace
static inline void stripSpaces(std::string &s){
    s.erase(std::remove_if(s.begin(), s.end(),
        [](unsigned char c){return std::isspace(c);}), s.end());
}

int VCFXRefComparator::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    std::string referencePath;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"reference", required_argument, 0, 'r'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv,"hr:", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'r':
                referencePath= optarg;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    if(referencePath.empty()){
        std::cerr<<"Error: must specify --reference <FASTA>.\n";
        displayHelp();
        return 1;
    }
    if(!loadReference(referencePath)){
        std::cerr<<"Error: failed to load reference from "<< referencePath<<"\n";
        return 1;
    }
    compareVCF(std::cin, std::cout);
    return 0;
}

void VCFXRefComparator::displayHelp(){
    std::cout <<
"VCFX_ref_comparator: Compare VCF REF/ALT with a reference genome.\n\n"
"Usage:\n"
"  VCFX_ref_comparator --reference ref.fasta < input.vcf > output.vcf\n\n"
"Description:\n"
"  Reads a reference FASTA into memory. Then reads each variant line:\n"
"   - If chromosome or position is invalid, logs a warning and sets REF_COMPARISON=UNKNOWN_CHROM or INVALID_POS.\n"
"   - Otherwise, compares the VCF's REF vs the reference substring. Then for each ALT, indicates 'REF_MATCH' if ALT= reference substring or 'NOVEL'.\n"
"  The result is appended to the 'INFO' field as REF_COMPARISON=...\n\n"
"Example:\n"
"  VCFX_ref_comparator --reference genome.fa < in.vcf > out.vcf\n";
}

bool VCFXRefComparator::loadReference(const std::string &referenceFastaPath){
    std::ifstream in(referenceFastaPath);
    if(!in.is_open()){
        std::cerr<<"Error: cannot open reference "<< referenceFastaPath<<"\n";
        return false;
    }
    referenceGenome.clear();
    std::string line, currentChrom;
    std::ostringstream seq;
    while(true){
        if(!std::getline(in, line)) break;
        if(line.empty()) continue;
        if(line[0]=='>'){
            // store old chrom
            if(!currentChrom.empty()){
                referenceGenome[currentChrom] = seq.str();
            }
            seq.str("");
            seq.clear();
            // parse new chrom
            currentChrom= line.substr(1);
            // if there's whitespace, strip it after first token
            {
                std::istringstream iss(currentChrom);
                iss >> currentChrom;
            }
            // uppercase
            std::transform(currentChrom.begin(), currentChrom.end(), currentChrom.begin(), ::toupper);
        } else {
            stripSpaces(line);
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq<< line;
        }
    }
    if(!currentChrom.empty()){
        referenceGenome[currentChrom]= seq.str();
    }
    return true;
}

void VCFXRefComparator::compareVCF(std::istream &vcfIn, std::ostream &vcfOut){
    bool foundChromHeader=false;
    infoHeaderInserted= false;
    std::string line;
    while(true){
        if(!std::getline(vcfIn, line)) break;
        if(line.empty()){
            vcfOut<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            // check if #CHROM
            if(line.rfind("#CHROM",0)==0){
                foundChromHeader= true;
                // Insert new INFO line BEFORE we write #CHROM
                if(!infoHeaderInserted){
                    vcfOut<<"##INFO=<ID=REF_COMPARISON,Number=1,Type=String,Description=\"Comparison of REF/ALT vs reference genome substring\">\n";
                    infoHeaderInserted=true;
                }
                vcfOut<< line<<"\n";
            } else {
                vcfOut<< line<<"\n";
            }
            continue;
        }
        if(!foundChromHeader){
            std::cerr<<"Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }
        // parse fields
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string col;
            while(std::getline(ss,col,'\t')){
                fields.push_back(col);
            }
        }
        if(fields.size()<8){
            std::cerr<<"Warning: VCF line has <8 columns => skipping.\n";
            continue;
        }
        // fields: 0=CHROM,1=POS,2=ID,3=REF,4=ALT,5=QUAL,6=FILTER,7=INFO,...
        std::string &chrom= fields[0];
        std::string &posStr= fields[1];
        std::string &ref= fields[3];
        std::string &alt= fields[4];
        std::string &info= fields[7];

        // uppercase chrom
        std::transform(chrom.begin(), chrom.end(), chrom.begin(), ::toupper);

        int pos=0;
        try{
            pos= std::stoi(posStr);
        }catch(...){
            // out of range
            if(!info.empty() && info.back()!=';') info+=';';
            info+= "REF_COMPARISON=INVALID_POS";
            // rewrite
            std::ostringstream newLine;
            for(int i=0; i<7; i++){
                newLine<< fields[i];
                if(i<7) newLine<<"\t";
            }
            newLine<< info;
            for(size_t i=8; i<fields.size(); i++){
                newLine<<"\t"<< fields[i];
            }
            vcfOut<< newLine.str()<<"\n";
            continue;
        }
        // find reference
        auto it= referenceGenome.find(chrom);
        if(it== referenceGenome.end()){
            // unknown chrom
            if(!info.empty() && info.back()!=';') info+=';';
            info+= "REF_COMPARISON=UNKNOWN_CHROM";
            // rewrite
            std::ostringstream newLine;
            for(int i=0; i<7; i++){
                newLine<< fields[i];
                if(i<7) newLine<<"\t";
            }
            newLine<< info;
            for(size_t i=8; i<fields.size(); i++){
                newLine<<"\t"<< fields[i];
            }
            vcfOut<< newLine.str()<<"\n";
            continue;
        }
        const std::string &seq= it->second;
        if(pos<1 || pos>(int)seq.size()){
            // invalid pos
            if(!info.empty() && info.back()!=';') info+=';';
            info+= "REF_COMPARISON=INVALID_POS";
            // rewrite
            std::ostringstream newLine;
            for(int i=0; i<7; i++){
                newLine<< fields[i];
                if(i<7) newLine<<"\t";
            }
            newLine<< info;
            for(size_t i=8; i<fields.size(); i++){
                newLine<<"\t"<< fields[i];
            }
            vcfOut<< newLine.str()<<"\n";
            continue;
        }
        // ref from genome
        std::string genomeRef= seq.substr(pos-1, ref.size()); // 1-based
        // uppercase
        std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);

        // compare ref vs genomeRef
        bool refMatch= (ref== genomeRef);

        // split alt by comma
        std::vector<std::string> altAlleles;
        {
            std::stringstream as(alt);
            std::string a;
            while(std::getline(as,a,',')){
                std::transform(a.begin(), a.end(), a.begin(), ::toupper);
                altAlleles.push_back(a);
            }
        }
        // for each alt, see if alt== genomeRef => "REF_MATCH" else "NOVEL"
        std::vector<std::string> comparisons;
        for(auto &a: altAlleles){
            if(a== genomeRef) comparisons.push_back("REF_MATCH");
            else comparisons.push_back("NOVEL");
        }
        // build comparisonStr
        std::string comparisonStr;
        for(size_t i=0; i< comparisons.size(); i++){
            if(i>0) comparisonStr+=",";
            comparisonStr+= comparisons[i];
        }
        if(!info.empty() && info.back()!=';') info+=';';
        // e.g. "REF_COMPARISON=REF_MATCH,REF_MATCH" or "NOVEL" etc.
        if(!refMatch){
            // optionally label mismatch => we won't specifically do that
            // the alt status is in the comparisons
        }
        info+= "REF_COMPARISON=" + comparisonStr;

        // rebuild line
        std::ostringstream outLine;
        for(int i=0;i<7;i++){
            outLine<< fields[i];
            if(i<6) outLine<<"\t";  // Changed from i<7 to i<6 to avoid adding an extra tab after FILTER
        }
        outLine<<"\t"<< info;  // Add tab after FILTER, then add INFO without spaces
        // any other columns
        for(size_t i=8; i< fields.size(); i++){
            outLine<<"\t"<< fields[i];
        }
        vcfOut<< outLine.str()<<"\n";
    }
}

int main(int argc, char* argv[]){
    VCFXRefComparator refComp;
    return refComp.run(argc, argv);
}


./VCFX_ref_comparator/VCFX_ref_comparator.h
#ifndef VCFX_REF_COMPARATOR_H
#define VCFX_REF_COMPARATOR_H

#include <iostream>
#include <string>
#include <unordered_map>

class VCFXRefComparator {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();
    bool loadReference(const std::string &referenceFastaPath);
    void compareVCF(std::istream &vcfIn, std::ostream &vcfOut);

    // chromosome -> uppercase sequence
    std::unordered_map<std::string, std::string> referenceGenome;

    // store whether we have already inserted the INFO line
    bool infoHeaderInserted = false;
};

#endif


./VCFX_reformatter/VCFX_reformatter.cpp
#include "VCFX_reformatter.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cctype>
#include <unordered_set>

// Helper to trim whitespace
static std::string trim(const std::string &s){
    size_t start=0; while(start<s.size() && isspace((unsigned char)s[start])) start++;
    if(start== s.size()) return "";
    size_t end= s.size()-1; 
    while(end>start && isspace((unsigned char)s[end])) end--;
    return s.substr(start, end-start+1);
}

// Splits a string by a delimiter into tokens
static void split(const std::string &s, char d, std::vector<std::string> &out){
    out.clear();
    std::stringstream ss(s);
    std::string t;
    while(std::getline(ss,t,d)){
        out.push_back(t);
    }
}

int VCFXReformatter::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    std::vector<std::string> compressInfoFields;
    std::vector<std::string> compressFormatFields;
    std::vector<std::string> reorderInfoFields;
    std::vector<std::string> reorderFormatFields;

    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"compress-info", required_argument, 0, 'c'},
        {"compress-format", required_argument, 0, 'f'},
        {"reorder-info", required_argument, 0, 'i'},
        {"reorder-format", required_argument, 0, 'o'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hc:f:i:o:", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'c': {
                // parse comma-separated
                std::stringstream ss(optarg);
                std::string val;
                while(std::getline(ss,val,',')){
                    val= trim(val);
                    if(!val.empty()) compressInfoFields.push_back(val);
                }
            }break;
            case 'f': {
                std::stringstream ss(optarg);
                std::string val;
                while(std::getline(ss,val,',')){
                    val=trim(val);
                    if(!val.empty()) compressFormatFields.push_back(val);
                }
            }break;
            case 'i': {
                std::stringstream ss(optarg);
                std::string val;
                while(std::getline(ss,val,',')){
                    val=trim(val);
                    if(!val.empty()) reorderInfoFields.push_back(val);
                }
            }break;
            case 'o': {
                std::stringstream ss(optarg);
                std::string val;
                while(std::getline(ss,val,',')){
                    val= trim(val);
                    if(!val.empty()) reorderFormatFields.push_back(val);
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
    reformatVCF(std::cin, std::cout, 
                compressInfoFields, 
                compressFormatFields,
                reorderInfoFields, 
                reorderFormatFields);
    return 0;
}

void VCFXReformatter::displayHelp(){
    std::cout <<
"VCFX_reformatter: Reformat INFO/FORMAT fields in a VCF.\n\n"
"Usage:\n"
"  VCFX_reformatter [options] < input.vcf > output.vcf\n\n"
"Options:\n"
"  -h, --help                     Show this help.\n"
"  -c, --compress-info <keys>     Remove these INFO keys, comma-separated.\n"
"  -f, --compress-format <keys>   Remove these FORMAT keys, comma-separated.\n"
"  -i, --reorder-info <keys>      Reorder these INFO keys at the front, leftover appended.\n"
"  -o, --reorder-format <keys>    Reorder these FORMAT keys at the front, leftover appended.\n\n"
"Example:\n"
"  VCFX_reformatter --compress-info AF,DP --reorder-info AF,DP < in.vcf > out.vcf\n"
"Description:\n"
"  This tool modifies data lines:\n"
"   * 'compress-info': remove specified keys from the semicolon INFO field.\n"
"   * 'compress-format': remove specified keys from the colon FORMAT field,\n"
"      and also remove them from each sample's subfield.\n"
"   * 'reorder-info': place specified keys in that order at the front, then\n"
"      append leftover keys in the order encountered.\n"
"   * 'reorder-format': reorder the FORMAT colon-delimited keys in #8 col,\n"
"      then reorder each sample's subfields accordingly.\n"
"  Lines with <8 columns are skipped with a warning. Header lines (#) are\n"
"  passed unmodified.\n";
}

void VCFXReformatter::reformatVCF(std::istream &in, std::ostream &out,
                                  const std::vector<std::string> &compressInfoFields,
                                  const std::vector<std::string> &compressFormatFields,
                                  const std::vector<std::string> &reorderInfoFields,
                                  const std::vector<std::string> &reorderFormatFields)
{
    // We'll store them in sets for faster removal checks
    std::unordered_set<std::string> infoToRemove(compressInfoFields.begin(), compressInfoFields.end());
    std::unordered_set<std::string> formatToRemove(compressFormatFields.begin(), compressFormatFields.end());

    bool foundChrom=false;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line <<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line <<"\n";
            if(line.rfind("#CHROM",0)==0) foundChrom=true;
            continue;
        }
        if(!foundChrom){
            std::cerr<<"Warning: data line before #CHROM => skipping.\n";
            continue;
        }
        std::vector<std::string> fields;
        {
            std::stringstream ss(line);
            std::string col;
            while(std::getline(ss,col,'\t')){
                fields.push_back(col);
            }
        }
        if(fields.size()<8){
            std::cerr<<"Warning: line with <8 columns => skipping.\n";
            continue;
        }
        std::string &infoField= fields[7];
        // compress info
        if(!infoToRemove.empty() && !infoField.empty() && infoField!="."){
            infoField= compressInfo(infoField, infoToRemove);
        }
        // reorder info
        if(!reorderInfoFields.empty() && !infoField.empty() && infoField!="."){
            infoField= reorderInfo(infoField, reorderInfoFields);
        }

        // if there's format => fields.size() might be >=9
        if(fields.size()>8){
            // compress format
            std::string &formatField= fields[8];
            if(!formatField.empty() && formatField!="."){
                // remove keys
                std::vector<int> keepIndices;
                if(!formatToRemove.empty()){
                    formatField= compressFormat(formatField, formatToRemove, keepIndices);
                } else {
                    // if not removing anything, keep all
                    std::vector<std::string> fsplit;
                    split(formatField, ':', fsplit);
                    keepIndices.resize(fsplit.size());
                    for(size_t i=0;i<fsplit.size();i++){
                        keepIndices[i]=(int)i;
                    }
                }
                // reorder format
                if(!reorderFormatFields.empty() && !formatField.empty() && formatField!="."){
                    std::vector<int> newIndices;
                    formatField= reorderFormat(formatField, reorderFormatFields, newIndices);
                    // apply new reorder to keepIndices
                    // we now want to chain keepIndices -> newIndices
                    // but if keepIndices had removed some, we must unify
                    std::vector<int> chained;
                    // newIndices has size = final # of keys
                    // each newIndices[i] => oldIndex
                    // but that oldIndex references the 'kept' keys from the original?
                    // Actually let's do a simpler approach: we do reorder first, then compress second. 
                    // But we do the user-specified approach in the order we coded => let's keep it.
                    // We'll do: if we reorder after compress, we must reorder the 'formatField' that came from compress, 
                    // but that's fine. We'll just do the reorder ignoring keepIndices for now.
                    keepIndices.clear();
                    // We must re-split the formatField to see how many remain
                    std::vector<std::string> finalFmtSplit;
                    split(formatField,':', finalFmtSplit);
                    keepIndices.resize(finalFmtSplit.size());
                    for(size_t i=0;i< finalFmtSplit.size(); i++){
                        keepIndices[i]= (int)i; 
                    }
                }

                // Now we must apply the final format structure to each sample
                for(size_t sampleI=9; sampleI< fields.size(); sampleI++){
                    fields[sampleI]= applyFormatReorderToSample(fields[sampleI], /*currently we do not reorder again*/ keepIndices);
                }
            }
        }
        // rebuild line
        std::ostringstream oss;
        for(size_t i=0; i<fields.size(); i++){
            if(i>0) oss<<"\t";
            oss<< fields[i];
        }
        out<< oss.str()<<"\n";
    }
}

std::string VCFXReformatter::compressInfo(const std::string &infoStr,
                                         const std::unordered_set<std::string> &keysToRemove)
{
    if(infoStr=="."|| infoStr.empty()) return infoStr;
    std::vector<std::string> tokens;
    split(infoStr, ';', tokens);
    std::vector<std::string> keep;
    for(auto &kv : tokens){
        if(kv.empty()) continue;
        // parse key=...
        auto eq= kv.find('=');
        if(eq== std::string::npos){
            // flag
            if(keysToRemove.count(kv)==0){
                keep.push_back(kv);
            }
        } else {
            std::string k= kv.substr(0,eq);
            if(keysToRemove.count(k)==0){
                keep.push_back(kv);
            }
        }
    }
    if(keep.empty()) return ".";
    std::ostringstream oss;
    for(size_t i=0; i< keep.size(); i++){
        if(i>0) oss<<";";
        oss<< keep[i];
    }
    return oss.str();
}

// remove keys from format column => rebuild format => store the indices of the keys we keep in keepIndices
std::string VCFXReformatter::compressFormat(const std::string &formatStr,
                                           const std::unordered_set<std::string> &keysToRemove,
                                           std::vector<int> &keepIndices)
{
    if(formatStr=="."||formatStr.empty()){
        keepIndices.clear(); 
        return formatStr;
    }
    std::vector<std::string> keys;
    split(formatStr,':', keys);
    keepIndices.clear();
    for(size_t i=0;i<keys.size();i++){
        if(keysToRemove.count(keys[i])==0){
            keepIndices.push_back((int)i);
        }
    }
    if(keepIndices.empty()) return "."; 
    // rebuild
    std::ostringstream oss;
    bool first=true;
    for(size_t i=0;i< keepIndices.size(); i++){
        if(!first) oss<<":";
        else first=false;
        oss<< keys[ keepIndices[i] ];
    }
    return oss.str();
}

std::string VCFXReformatter::reorderInfo(const std::string &infoStr, 
                                         const std::vector<std::string> &order)
{
    if(infoStr=="."|| infoStr.empty()) return infoStr;
    std::vector<std::string> tokens;
    split(infoStr,';', tokens);
    // parse into map + remember original order
    std::unordered_map<std::string,std::string> kvMap;
    std::vector<std::string> originalKeys;
    for(auto &item: tokens){
        if(item.empty()) continue;
        auto eq= item.find('=');
        if(eq==std::string::npos){
            // flag
            kvMap[item]="";
            originalKeys.push_back(item);
        } else {
            std::string k= item.substr(0,eq);
            std::string v= item.substr(eq+1);
            kvMap[k]= v;
            originalKeys.push_back(k);
        }
    }
    // build new list
    std::vector<std::string> result;
    // 1) add items from 'order' if exist in kvMap
    for(auto &k : order){
        auto it= kvMap.find(k);
        if(it!= kvMap.end()){
            if(it->second==""){
                result.push_back(k); // flag
            } else {
                result.push_back(k+"="+ it->second);
            }
            kvMap.erase(it);
        }
    }
    // 2) append leftover in original order
    for(auto &k: originalKeys){
        auto it= kvMap.find(k);
        if(it!= kvMap.end()){
            if(it->second==""){
                result.push_back(k);
            } else {
                result.push_back(k+"="+ it->second);
            }
            kvMap.erase(it);
        }
    }
    if(result.empty()) return ".";
    // join
    std::ostringstream oss;
    for(size_t i=0;i< result.size(); i++){
        if(i>0) oss<<";";
        oss<< result[i];
    }
    return oss.str();
}

std::string VCFXReformatter::reorderFormat(const std::string &fmtStr,
                                          const std::vector<std::string> &order,
                                          std::vector<int> &oldToNew)
{
    if(fmtStr=="."|| fmtStr.empty()){
        oldToNew.clear();
        return fmtStr;
    }
    std::vector<std::string> keys;
    split(fmtStr,':', keys);
    // build result
    std::vector<std::string> newOrder;
    newOrder.reserve(keys.size());
    oldToNew.assign(keys.size(), -1); // -1 => removed

    // first place the requested keys in order
    int usedCount=0;
    // keep track of which old indices are used
    std::vector<bool> used(keys.size(), false);

    // For each key in 'order', see if we find it in keys
    for(auto &k : order){
        auto it= std::find(keys.begin(), keys.end(), k);
        if(it!= keys.end()){
            int oldI= (int)std::distance(keys.begin(), it);
            newOrder.push_back(k);
            oldToNew[oldI]= usedCount;
            usedCount++;
            used[oldI]=true;
        }
    }
    // then append leftover in original order
    for(int i=0; i<(int)keys.size(); i++){
        if(!used[i]){
            newOrder.push_back(keys[i]);
            oldToNew[i]= usedCount;
            usedCount++;
        }
    }
    // build final string
    if(newOrder.empty()) {
        oldToNew.clear();
        return ".";
    }
    std::ostringstream oss;
    for(size_t i=0;i< newOrder.size(); i++){
        if(i>0) oss<<":";
        oss<< newOrder[i];
    }
    return oss.str();
}

// We reorder or remove subfields in each sample column
std::string VCFXReformatter::applyFormatReorderToSample(const std::string &sampleStr,
                                                        const std::vector<int> &oldToNew)
{
    // if oldToNew is empty, means format is "." => no sample data
    if(oldToNew.empty()){
        return sampleStr; 
    }
    if(sampleStr=="."|| sampleStr.empty()) return sampleStr;
    std::vector<std::string> subs;
    split(sampleStr, ':', subs);
    // build a new list of same length as oldToNew, initially "."
    // if oldToNew[i] < 0 => that subfield is removed
    // else => newIndex= oldToNew[i]
    // we want final vector with size = maxIndex+1
    int maxIndex= -1;
    for(auto idx: oldToNew){
        if(idx> maxIndex) maxIndex= idx;
    }
    if(maxIndex<0) {
        // means everything removed
        return "."; 
    }
    std::vector<std::string> newSubs(maxIndex+1, ".");
    for(int oldI=0; oldI<(int)subs.size() && oldI<(int)oldToNew.size(); oldI++){
        int ni= oldToNew[oldI];
        if(ni>=0){
            if((size_t)oldI< subs.size()){
                newSubs[ni]= subs[oldI];
            } else {
                newSubs[ni]=".";
            }
        }
    }
    // if all are "." => return "."
    bool allDot=true;
    for(auto &x: newSubs){
        if(x!=".") { allDot=false; break;}
    }
    if(allDot) return ".";
    // join with colon
    std::ostringstream oss;
    bool first=true;
    for(size_t i=0;i< newSubs.size(); i++){
        if(!first) oss<<":";
        else first=false;
        oss<< newSubs[i];
    }
    return oss.str();
}

int main(int argc, char* argv[]){
    VCFXReformatter reformatter;
    return reformatter.run(argc, argv);
}


./VCFX_reformatter/VCFX_reformatter.h
#ifndef VCFX_REFORMATTER_H
#define VCFX_REFORMATTER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>

class VCFXReformatter {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    // Reformat the VCF from 'in' to 'out' using the user-specified lists
    void reformatVCF(std::istream &in, std::ostream &out,
                     const std::vector<std::string> &compressInfoFields,
                     const std::vector<std::string> &compressFormatFields,
                     const std::vector<std::string> &reorderInfo,
                     const std::vector<std::string> &reorderFormat);

    // Remove user-specified fields from the semicolon-based INFO
    // (Used for compressing info keys)
    std::string compressInfo(const std::string &infoStr,
                             const std::unordered_set<std::string> &keysToRemove);

    // Remove user-specified keys from colon-based FORMAT
    // and from each genotype subfield in that position
    // Returns the new format string, plus an index mapping
    std::string compressFormat(const std::string &formatStr,
                               const std::unordered_set<std::string> &keysToRemove,
                               std::vector<int> &keepIndices);

    // Reorder a semicolon-based INFO
    std::string reorderInfo(const std::string &infoStr,
                            const std::vector<std::string> &order);

    // Reorder a colon-based FORMAT; returns new format, plus old->new index mapping
    std::string reorderFormat(const std::string &fmtStr,
                              const std::vector<std::string> &order,
                              std::vector<int> &oldToNew);

    // reorder genotype subfields for each sample based on oldToNew
    // or remove subfields if the new index is negative
    std::string applyFormatReorderToSample(const std::string &sampleStr,
                                           const std::vector<int> &oldToNew);
};

#endif


./VCFX_region_subsampler/VCFX_region_subsampler.cpp
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


./VCFX_region_subsampler/VCFX_region_subsampler.h
#ifndef VCFX_REGION_SUBSAMPLER_H
#define VCFX_REGION_SUBSAMPLER_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Each region is (start, end), inclusive, 1-based
struct Region {
    int start;
    int end;
};

// VCFXRegionSubsampler: 
// Reads a BED file with multiple lines => chromosome -> sorted intervals
// Then reads a VCF and keeps lines whose POS is within any interval for that CHROM.
class VCFXRegionSubsampler {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();
    bool loadRegions(const std::string &bedFile,
                     std::unordered_map<std::string, std::vector<Region>> &chromRegions);
    // after loading, we sort the intervals for each chromosome, possibly merge them
    void sortAndMergeIntervals(std::unordered_map<std::string, std::vector<Region>> &chromRegions);

    bool isInAnyRegion(const std::string &chrom, int pos) const;

    // The map from chrom => sorted intervals
    std::unordered_map<std::string, std::vector<Region>> regions;

    void processVCF(std::istream &in, std::ostream &out);

};

#endif


./VCFX_sample_extractor/VCFX_sample_extractor.cpp
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

int main(int argc, char* argv[]){
    VCFXSampleExtractor app;
    return app.run(argc, argv);
}


./VCFX_sample_extractor/VCFX_sample_extractor.h
#ifndef VCFX_SAMPLE_EXTRACTOR_H
#define VCFX_SAMPLE_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>

class VCFXSampleExtractor {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    // read the user’s --samples list => store them
    // read VCF from in => output a valid VCF with only those samples
    void extractSamples(std::istream &in, std::ostream &out,
                        const std::vector<std::string> &samples);

};

#endif


./VCFX_sorter/VCFX_sorter.cpp
#include "VCFX_sorter.h"
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdlib>

void VCFXSorter::displayHelp(){
    std::cout <<
"VCFX_sorter: Sort a VCF by chromosome and position.\n\n"
"Usage:\n"
"  VCFX_sorter [options] < input.vcf > output.vcf\n\n"
"Options:\n"
"  -h, --help          Show help.\n"
"  -n, --natural-chr   Use a natural chromosome sort (chr1 < chr2 < chr10) instead of lexicographic.\n\n"
"Description:\n"
"  Reads all data lines into memory, sorts by (CHROM,POS). Preserves all header lines\n"
"  in original order, and outputs them first, then prints sorted data lines.\n\n"
"Examples:\n"
"  1) Lexicographic:\n"
"     VCFX_sorter < unsorted.vcf > sorted.vcf\n"
"  2) Natural order:\n"
"     VCFX_sorter --natural-chr < unsorted.vcf > sorted.vcf\n";
}

// A function to parse chromosome name in a natural manner if possible
// e.g. "chr10" => (chr,10) so that chr2 < chr10 in numeric sense
// We'll remove "chr" prefix if present, then parse leading digits
static bool parseChromNat(const std::string &chrom,
                          std::string &prefix,
                          long &num,
                          std::string &suffix)
{
    // e.g. "chr10_gl000" => prefix="chr", then 10, then suffix="_gl000"
    // We'll do a naive approach: remove leading "chr" or "Chr" or "CHR"
    std::string c= chrom;
    // uppercase to check
    std::string up;
    up.reserve(c.size());
    for(char ch: c) up.push_back(std::toupper(ch));
    if(up.rfind("CHR",0)==0) {
        // remove "chr" prefix
        prefix= c.substr(0,3); 
        c= c.substr(3);
    } else {
        prefix="";
    }
    // now parse leading digits
    // e.g. c="10_gl000"
    // we parse digits until non-digit
    size_t idx=0;
    while(idx< c.size() && std::isdigit((unsigned char)c[idx])) idx++;
    if(idx==0){
        // means no leading digit
        num=-1;
        suffix= c;
        return true;
    }
    // parse c[0..idx-1] as a number
    try {
        num= std::stol(c.substr(0,idx));
    } catch(...){
        return false;
    }
    suffix= c.substr(idx);
    return true;
}

// We'll define the lexCompare
bool VCFRecord::lexCompare(const VCFRecord &a, const VCFRecord &b){
    if(a.chrom != b.chrom) return a.chrom < b.chrom;
    return a.pos < b.pos;
}

// We'll define the naturalCompare
bool VCFRecord::naturalCompare(const VCFRecord &a, const VCFRecord &b){
    // parse each
    std::string apfx, asuf, bpfx, bsuf;
    long anum, bnum;
    if(!parseChromNat(a.chrom, apfx, anum, asuf) || 
       !parseChromNat(b.chrom, bpfx, bnum, bsuf)) {
        // fallback to lex
        if(a.chrom != b.chrom) return a.chrom < b.chrom;
        return a.pos< b.pos;
    }
    // if prefix differs, lex compare them
    if(apfx != bpfx) return (apfx< bpfx);
    // if both have num>=0
    if(anum>=0 && bnum>=0){
        if(anum!= bnum) return (anum< bnum);
        // if suffix differs => lex compare
        if(asuf!= bsuf) return asuf< bsuf;
        // else compare pos
        if(a.pos!= b.pos) return a.pos< b.pos;
        // all tie
        return false;
    } else if(anum>=0 && bnum<0){
        // numeric < no numeric
        return true;
    } else if(anum<0 && bnum>=0){
        return false;
    } else {
        // both have no numeric => fallback to lex compare
        if(a.chrom!= b.chrom) return a.chrom< b.chrom;
        return a.pos< b.pos;
    }
}

int VCFXSorter::run(int argc, char* argv[]){
    bool showHelp=false;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"natural-chr", no_argument, 0, 'n'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hn", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'n':
                naturalChromOrder= true;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    
    // Check if stdin has data
    if (std::cin.peek() == EOF) {
        displayHelp();
        return 0;
    }
    
    loadVCF(std::cin);
    sortRecords();
    outputVCF(std::cout);
    return 0;
}

// loads lines into memory
void VCFXSorter::loadVCF(std::istream &in){
    bool foundChrom=false;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            // we keep it in the header or skip? We'll keep it in header
            headerLines.push_back(line);
            continue;
        }
        if(line[0]=='#'){
            headerLines.push_back(line);
            if(line.rfind("#CHROM",0)==0) foundChrom=true;
            continue;
        }
        // parse data
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string col;
            while(std::getline(ss,col,'\t')){
                fields.push_back(col);
            }
        }
        if(fields.size()<8){
            std::cerr<<"Warning: skipping line with <8 columns.\n";
            continue;
        }
        VCFRecord rec;
        rec.chrom= fields[0];
        try {
            rec.pos= std::stoi(fields[1]);
        } catch(...){
            std::cerr<<"Warning: invalid POS => skipping line.\n";
            continue;
        }
        rec.fields= std::move(fields);
        records.push_back(rec);
    }
    if(!foundChrom){
        std::cerr<<"Warning: no #CHROM line found in input.\n";
    }
}

// sort them in memory
void VCFXSorter::sortRecords(){
    if(naturalChromOrder){
        std::sort(records.begin(), records.end(), VCFRecord::naturalCompare);
    } else {
        std::sort(records.begin(), records.end(), VCFRecord::lexCompare);
    }
}

void VCFXSorter::outputVCF(std::ostream &out){
    // print all header lines
    for(auto &l : headerLines){
        out<< l <<"\n";
    }
    // then print sorted data
    for(auto &rec: records){
        // rebuild line from rec.fields
        bool first=true;
        for(size_t i=0;i< rec.fields.size(); i++){
            if(!first) out<<"\t";
            else first=false;
            out<< rec.fields[i];
        }
        out<<"\n";
    }
}

int main(int argc, char* argv[]){
    VCFXSorter app;
    return app.run(argc, argv);
}

./VCFX_sorter/VCFX_sorter.h
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


./VCFX_subsampler/VCFX_subsampler.cpp
#include "VCFX_subsampler.h"
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <cstdlib>

void VCFXSubsampler::displayHelp(){
    std::cout <<
"VCFX_subsampler: Randomly pick N lines from a VCF data section.\n\n"
"Usage:\n"
"  VCFX_subsampler [options] < input.vcf > output.vcf\n\n"
"Options:\n"
"  -s, --subsample <N>   Required: number of data lines (variants) to keep.\n"
"  --seed <INT>          Use a reproducible random seed.\n"
"  -h, --help            Show this help.\n\n"
"Description:\n"
"  We read all header lines (#...) first and output them as-is. Then we do\n"
"  reservoir sampling on subsequent lines (the data lines). If the file has\n"
"  fewer than N lines, we keep them all. We skip lines with <8 columns.\n\n"
"Example:\n"
"  VCFX_subsampler --subsample 1000 < big.vcf > subset.vcf\n"
"  VCFX_subsampler --subsample 1000 --seed 1234 < big.vcf > subset2.vcf\n";
}

int VCFXSubsampler::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    int sampleSize=0;
    unsigned int seed= (unsigned int)std::time(nullptr);
    bool seedSpecified= false;

    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"subsample", required_argument, 0, 's'},
        {"seed", required_argument, 0, 1000},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hs:", long_opts,nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 's': {
                try {
                    sampleSize= std::stoi(optarg);
                    if(sampleSize<=0){
                        throw std::invalid_argument("must be >0");
                    }
                } catch(...){
                    std::cerr<<"Error: invalid subsample size.\n";
                    return 1;
                }
            } break;
            case 1000: { // --seed
                seedSpecified= true;
                try {
                    long long val= std::stoll(optarg);
                    seed= (unsigned int)val;
                } catch(...){
                    std::cerr<<"Error: invalid seed.\n";
                    return 1;
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
    if(sampleSize<=0){
        std::cerr<<"Error: must specify --subsample <N> with N>0.\n";
        return 1;
    }
    subsampleLines(std::cin, std::cout, sampleSize, seed);
    return 0;
}

void VCFXSubsampler::subsampleLines(std::istream &in, std::ostream &out,
                                    int sampleSize, unsigned int seed)
{
    std::string line;
    // store all # lines first
    while(true){
        std::streampos p= in.tellg();
        if(!std::getline(in,line)) break;
        if(line.empty()){
            out<< line<<"\n";
            continue;
        }
        if(line[0]=='#'){
            out<< line<<"\n";
            continue;
        } else {
            // not a header
            in.seekg(p);
            break;
        }
    }
    // now do reservoir sampling on data lines
    std::vector<std::string> reservoir;
    reservoir.reserve(sampleSize);

    int count=0;
    // random gen
    std::default_random_engine rng(seed);
    // read data lines
    while(true){
        if(!std::getline(in, line)) break;
        if(line.empty()) continue;
        // skip lines with <8 columns
        {
            std::stringstream s2(line);
            std::vector<std::string> fields;
            std::string tmp;
            while(std::getline(s2,tmp,'\t')){
                fields.push_back(tmp);
            }
            if(fields.size()<8){
                std::cerr<<"Warning: skipping line with <8 columns.\n";
                continue;
            }
        }

        if(count< sampleSize){
            reservoir.push_back(line);
        } else {
            // pick random int in [0..count]
            std::uniform_int_distribution<int> dist(0, count);
            int j= dist(rng);
            if(j< sampleSize){
                reservoir[j]= line;
            }
        }
        count++;
    }
    // output reservoir
    for(auto &r: reservoir){
        out<< r<<"\n";
    }
}


int main(int argc, char* argv[]){
    VCFXSubsampler app;
    return app.run(argc, argv);
}


./VCFX_subsampler/VCFX_subsampler.h
#ifndef VCFX_SUBSAMPLER_H
#define VCFX_SUBSAMPLER_H

#include <iostream>
#include <string>
#include <vector>

class VCFXSubsampler {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    // Reservoir sampling method
    void subsampleLines(std::istream &in, std::ostream &out,
                        int sampleSize, unsigned int seed);

};

#endif


./VCFX_sv_handler/VCFX_sv_handler.cpp
#include "VCFX_sv_handler.h"
#include <getopt.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

int VCFXSvHandler::run(int argc, char* argv[]){
    if(argc==1){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    bool filterOnly=false;
    bool modifySV=false;

    static struct option longOpts[]={
        {"help", no_argument, 0,'h'},
        {"sv-filter-only", no_argument, 0,'f'},
        {"sv-modify", no_argument, 0,'m'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hfm", longOpts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h': showHelp=true; break;
            case 'f': filterOnly=true; break;
            case 'm': modifySV=true; break;
            default: showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }

    handleStructuralVariants(std::cin, std::cout, filterOnly, modifySV);
    return 0;
}

void VCFXSvHandler::displayHelp(){
    std::cout <<
"VCFX_sv_handler: Filter or modify structural variants in a VCF.\n\n"
"Usage:\n"
"  VCFX_sv_handler [options] < input.vcf > output.vcf\n\n"
"Options:\n"
"  -h, --help           Show this help.\n"
"  -f, --sv-filter-only Keep only lines that have 'SVTYPE=' in their INFO.\n"
"  -m, --sv-modify      Modify the INFO field of structural variants.\n\n"
"Description:\n"
"  * If --sv-filter-only is set, we skip lines without structural variant.\n"
"  * If --sv-modify is set, we add 'SV_VALIDATED=1', 'SV_SIZE=...' for DEL/DUP.\n"
"    Also 'INV_TYPE=PARALLEL' for INV, 'BND_ORIENTATION=PAIR' for BND. etc.\n"
"  * If both are set, we do both filtering and modification.\n"
"  * Non-SV lines are only included if !filterOnly.\n\n"
"Example:\n"
"  1) Keep only structural variants:\n"
"     VCFX_sv_handler --sv-filter-only < in.vcf > out.vcf\n"
"  2) Modify structural variants:\n"
"     VCFX_sv_handler --sv-modify < in.vcf > out.vcf\n"
"  3) Do both:\n"
"     VCFX_sv_handler --sv-filter-only --sv-modify < in.vcf > out.vcf\n";
}

bool VCFXSvHandler::isStructuralVariant(const std::string &infoField) const{
    return infoField.find("SVTYPE=")!= std::string::npos;
}

std::string VCFXSvHandler::parseSVType(const std::string &infoField) const{
    // find "SVTYPE="
    auto pos= infoField.find("SVTYPE=");
    if(pos== std::string::npos) return "";
    size_t start= pos+7; 
    // read until next ; or end
    auto semicol= infoField.find(';', start);
    if(semicol== std::string::npos){
        return infoField.substr(start);
    } else {
        return infoField.substr(start, semicol - start);
    }
}

int VCFXSvHandler::parseEndPosition(const std::string &infoField) const{
    auto pos= infoField.find("END=");
    if(pos== std::string::npos) return -1;
    size_t start= pos+4;
    auto semicol= infoField.find(';', start);
    std::string endStr= (semicol== std::string::npos) ? infoField.substr(start)
                                                      : infoField.substr(start, semicol - start);
    try {
        return std::stoi(endStr);
    } catch(...){
        return -1;
    }
}

int VCFXSvHandler::parsePos(const std::string &posField) const{
    try{
        return std::stoi(posField);
    }catch(...){
        return -1;
    }
}

std::string VCFXSvHandler::manipulateSVInfo(const std::string &infoField,
                                            const std::string &svType,
                                            int pos,
                                            int endPos) const
{
    std::string modified= infoField;
    if(!modified.empty() && modified.back()!=';'){
        modified+=";";
    }
    modified+= "SV_VALIDATED=1"; // a sample annotation

    // if we have endPos>0 and pos>0 => maybe we compute size
    if(endPos>0 && pos>0 && (svType=="DEL"||svType=="DUP")){
        int svSize= endPos - pos;
        // possibly handle negative
        if(svSize<0) svSize= -svSize; 
        modified += ";SV_SIZE=" + std::to_string(svSize);
    }
    if(svType=="INV"){
        modified += ";INV_TYPE=PARALLEL";
    } else if(svType=="BND"){
        modified += ";BND_ORIENTATION=PAIR";
    }
    // other types?
    return modified;
}

void VCFXSvHandler::handleStructuralVariants(std::istream &in, std::ostream &out,
                                            bool filterOnly, bool modifySV)
{
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        if(line.empty()){
            continue;
        }
        if(line[0]=='#'){
            out<< line<<"\n";
            continue;
        }
        // parse
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string col;
            while(std::getline(ss,col,'\t')){
                fields.push_back(col);
            }
        }
        if(fields.size()<8){
            std::cerr<<"Warning: skipping line with <8 columns.\n";
            continue;
        }
        bool isSV= isStructuralVariant(fields[7]);
        if(isSV){
            // if filterOnly && !modify => just output
            // if modify => we parse and rewrite info
            if(filterOnly && !modifySV){
                out<< line<<"\n";
                continue;
            }
            if(modifySV){
                // parse svtype
                std::string svType= parseSVType(fields[7]);
                if(svType.empty()){
                    std::cerr<<"Warning: no SVTYPE => skipping line.\n";
                    continue;
                }
                int pos= parsePos(fields[1]);
                int endPos= parseEndPosition(fields[7]);
                if(pos<0){
                    std::cerr<<"Warning: invalid POS => skipping.\n";
                    continue;
                }
                std::string newInfo= manipulateSVInfo(fields[7], svType, pos, endPos);
                fields[7]= newInfo;
                // rebuild
                bool first=true;
                for(size_t i=0;i<fields.size(); i++){
                    if(!first) out<<"\t";
                    else first=false;
                    out<< fields[i];
                }
                out<<"\n";
                continue;
            }
            // if we get here => isSV, but not filterOnly or modify => just output
            out<< line<<"\n";
        } else {
            // not an SV
            if(!filterOnly){
                // keep it
                out<< line<<"\n";
            }
            // else skip
        }
    }
}

int main(int argc, char* argv[]){
    VCFXSvHandler app;
    return app.run(argc, argv);
}


./VCFX_sv_handler/VCFX_sv_handler.h
#ifndef VCFX_SV_HANDLER_H
#define VCFX_SV_HANDLER_H

#include <iostream>
#include <string>

class VCFXSvHandler {
public:
    int run(int argc, char* argv[]);

private:
    void displayHelp();

    // The main method to read lines from 'in' and apply filtering/modify logic, then write to 'out'.
    void handleStructuralVariants(std::istream &in, std::ostream &out, bool filterOnly, bool modifySV);

    // Checks if a line's INFO indicates an SV (i.e. has "SVTYPE=").
    bool isStructuralVariant(const std::string &infoField) const;

    // Extract the SVTYPE=... substring from INFO; returns "" if not found.
    std::string parseSVType(const std::string &infoField) const;

    // Extract 'END=' from INFO; returns -1 if not found or invalid.
    int parseEndPosition(const std::string &infoField) const;

    // Parse position from string. Return -1 on error.
    int parsePos(const std::string &posField) const;

    // If we are modifying, do the manipulations (like adding SV_SIZE=..., etc.)
    std::string manipulateSVInfo(const std::string &infoField,
                                 const std::string &svType,
                                 int pos, int endPos) const;
};

#endif


./VCFX_validator/VCFX_validator.cpp
#include "VCFX_validator.h"
#include <getopt.h>
#include <sstream>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <unistd.h>

static std::string trim(const std::string &s){
    size_t start=0; 
    while(start<s.size() && std::isspace((unsigned char)s[start])) start++;
    if(start== s.size()) return "";
    size_t end= s.size()-1;
    while(end> start && std::isspace((unsigned char)s[end])) end--;
    return s.substr(start, end-start+1);
}

int VCFXValidator::run(int argc, char* argv[]){
    bool hasStdin = !isatty(fileno(stdin));
    if(argc==1 && !hasStdin){
        displayHelp();
        return 0;
    }
    bool showHelp=false;
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"strict", no_argument, 0, 's'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv, "hs", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h': showHelp=true; break;
            case 's': strictMode= true; break;
            default: showHelp= true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    bool ok= validateVCF(std::cin);
    return (ok ? 0 : 1);
}

void VCFXValidator::displayHelp(){
    std::cout <<
"VCFX_validator: Checks basic validity of a VCF.\n\n"
"Usage:\n"
"  VCFX_validator [options] < input.vcf\n\n"
"Options:\n"
"  -h, --help     Show this help.\n"
"  -s, --strict   Enable stricter checks (not fully implemented, but reserved).\n\n"
"Description:\n"
"  Validates:\n"
"   * All '##' lines are recognized as meta lines.\n"
"   * #CHROM line is present, has at least 8 columns.\n"
"   * Each data line has >=8 columns, checks CHROM non-empty, POS>0,\n"
"     REF/ALT non-empty, QUAL is '.' or non-negative float, FILTER non-empty,\n"
"     INFO is minimal check. Logs errors/warnings.\n"
"  Exits 0 if pass, 1 if fail.\n";
}

bool VCFXValidator::validateMetaLine(const std::string &line, int lineNumber){
    // minimal check: must start with "##"
    if(line.size()<2) return false;
    if(line[0]=='#' && line[1]=='#'){
        return true;
    }
    std::cerr<<"Error: line "<< lineNumber <<" is a header line but doesn't start with '##'.\n";
    return false;
}

bool VCFXValidator::validateChromHeader(const std::string &line, int lineNumber){
    // parse by tab
    std::vector<std::string> f;
    {
        std::stringstream ss(line);
        std::string col;
        while(std::getline(ss,col,'\t')){
            f.push_back(col);
        }
    }
    if(f.size()<8){
        std::cerr<<"Error: #CHROM line at "<< lineNumber <<" has <8 columns.\n";
        return false;
    }
    // typically #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, (FORMAT?), ...
    if(f[0]!="#CHROM"){
        std::cerr<<"Error: #CHROM line doesn't start with '#CHROM' at line "<< lineNumber <<".\n";
        return false;
    }
    return true;
}

bool VCFXValidator::validateDataLine(const std::string &line, int lineNumber){
    // parse
    std::vector<std::string> f;
    {
        std::stringstream ss(line);
        std::string col;
        while(std::getline(ss,col,'\t')){
            f.push_back(trim(col));
        }
    }
    if(f.size()<8){
        std::cerr<<"Error: line "<< lineNumber <<" has <8 columns.\n";
        return false;
    }
    // CHROM
    if(f[0].empty()){
        std::cerr<<"Error: line "<<lineNumber<<" CHROM is empty.\n";
        return false;
    }
    // POS
    {
        int pos=0;
        try{
            pos= std::stoi(f[1]);
            if(pos<=0){
                std::cerr<<"Error: line "<<lineNumber<<" POS must be >0.\n";
                return false;
            }
        } catch(...){
            std::cerr<<"Error: line "<<lineNumber<<" POS not parseable.\n";
            return false;
        }
    }
    // ID can be empty => that's fine
    // REF
    if(f[3].empty()){
        std::cerr<<"Error: line "<< lineNumber<<" REF is empty.\n";
        return false;
    }
    // ALT
    if(f[4].empty()){
        std::cerr<<"Error: line "<< lineNumber<<" ALT is empty.\n";
        return false;
    }
    // QUAL => '.' or non-neg float
    if(f[5]!="."){
        try{
            double q= std::stod(f[5]);
            if(q<0){
                std::cerr<<"Error: line "<<lineNumber<<" negative QUAL.\n";
                return false;
            }
        } catch(...){
            std::cerr<<"Error: line "<<lineNumber<<" invalid QUAL.\n";
            return false;
        }
    }
    // FILTER => must not be empty
    if(f[6].empty()){
        std::cerr<<"Error: line "<<lineNumber<<" FILTER is empty.\n";
        return false;
    }
    // INFO => minimal check
    // could be '.' => ok
    // else splitted by ';' => each either 'key=val' or 'flag'
    if(f[7]!="."){
        std::stringstream infoSS(f[7]);
        std::string token;
        bool anyValid= false;
        while(std::getline(infoSS, token, ';')){
            token= trim(token);
            if(token.empty()) continue;
            auto eq= token.find('=');
            if(eq== std::string::npos){
                // treat as flag
                anyValid= true;
            } else {
                // key=val
                std::string k= token.substr(0, eq);
                // std::string v= token.substr(eq+1);
                if(k.empty()){
                    std::cerr<<"Error: line "<<lineNumber<<" has INFO with empty key.\n";
                    return false;
                }
                anyValid= true;
            }
        }
        // if not anyValid => error
        if(!anyValid){
            std::cerr<<"Error: line "<<lineNumber<<" INFO not valid.\n";
            return false;
        }
    }
    return true;
}

bool VCFXValidator::validateVCF(std::istream &in){
    std::string line;
    int lineNum=0;
    bool foundChromLine= false;

    while(true){
        if(!std::getline(in, line)) break;
        lineNum++;
        if(line.empty()){
            // skip
            continue;
        }
        if(line[0]=='#'){
            // could be ## or #CHROM
            if(line.rfind("##",0)==0){
                // validate meta
                if(!validateMetaLine(line, lineNum)) return false;
            } else if(line.rfind("#CHROM",0)==0){
                // validate #CHROM
                if(!validateChromHeader(line, lineNum)) return false;
                foundChromLine= true;
            } else {
                // Any line starting with '#' but not "##" or "#CHROM" is invalid
                std::cerr<<"Error: line "<<lineNum<<" is a header line but neither starts with '##' nor is a #CHROM header line.\n";
                return false;
            }
        } else {
            // data line
            if(!foundChromLine){
                std::cerr<<"Error: data line encountered before #CHROM at line "<<lineNum<<".\n";
                return false;
            }
            if(!validateDataLine(line, lineNum)) return false;
        }
    }
    if(!foundChromLine){
        std::cerr<<"Error: no #CHROM line found in file.\n";
        return false;
    }
    std::cout<<"VCF file is valid.\n";
    return true;
}

int main(int argc, char* argv[]){
    VCFXValidator validator;
    return validator.run(argc, argv);
}


./VCFX_validator/VCFX_validator.h
#ifndef VCFX_VALIDATOR_H
#define VCFX_VALIDATOR_H

#include <iostream>
#include <string>

class VCFXValidator {
public:
    int run(int argc, char* argv[]);

private:
    // If we add advanced checks for e.g. "strict" mode, we store a bool
    bool strictMode = false;

    // Show usage
    void displayHelp();

    // Main function that reads lines from in, does validation
    bool validateVCF(std::istream &in);

    // Validate a meta line "##" and returns true if it’s correct
    bool validateMetaLine(const std::string &line, int lineNumber);

    // Check #CHROM line
    bool validateChromHeader(const std::string &line, int lineNumber);

    // Validate a data line with at least 8 columns
    bool validateDataLine(const std::string &line, int lineNumber);
};

#endif


./VCFX_variant_classifier/VCFX_variant_classifier.cpp
#include "VCFX_variant_classifier.h"
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>
#include <unordered_set>
#include <poll.h>  // Add this for poll() function

// A helper to trim spaces
static std::string trim(const std::string &s){
    size_t start=0; 
    while(start<s.size() && std::isspace((unsigned char)s[start])) start++;
    if(start== s.size()) return "";
    size_t end= s.size()-1;
    while(end> start && std::isspace((unsigned char)s[end])) end--;
    return s.substr(start, end-start+1);
}

int VCFXVariantClassifier::run(int argc, char* argv[]){
    bool showHelp = false;
    
    // Check if stdin has data available when no arguments
    if(argc == 1) {
        // Check if stdin has data using poll
        struct pollfd fds;
        fds.fd = 0; // stdin
        fds.events = POLLIN;
        
        // Poll stdin with 0 timeout to check if data is available
        int ret = poll(&fds, 1, 0);
        
        // Show help if no data is available on stdin
        if (ret <= 0 || !(fds.revents & POLLIN)) {
            displayHelp();
            return 0;
        }
    }
    
    static struct option long_opts[]={
        {"help", no_argument, 0, 'h'},
        {"append-info", no_argument, 0, 'a'},
        {0,0,0,0}
    };
    while(true){
        int c= ::getopt_long(argc, argv,"ha", long_opts, nullptr);
        if(c==-1) break;
        switch(c){
            case 'h':
                showHelp=true;
                break;
            case 'a':
                appendInfo=true;
                break;
            default:
                showHelp=true;
        }
    }
    if(showHelp){
        displayHelp();
        return 0;
    }
    classifyStream(std::cin, std::cout);
    return 0;
}

void VCFXVariantClassifier::displayHelp(){
    std::cout <<
"VCFX_variant_classifier: Classify variants in a VCF as SNP, INDEL, MNV, or STRUCTURAL.\n\n"
"Usage:\n"
"  VCFX_variant_classifier [options] < input.vcf > output.vcf_or_tsv\n\n"
"Options:\n"
"  -h, --help         Show help.\n"
"  -a, --append-info  Instead of producing a TSV, output a valid VCF\n"
"                     with a new 'VCF_CLASS' subfield in the INFO.\n\n"
"Description:\n"
"  Reads each variant line, determines if it is:\n"
"    SNP: single base ref & alt,\n"
"    INDEL: length mismatch (less than 50 bp difference) in ref vs alt,\n"
"    MNV: same length >1,\n"
"    STRUCTURAL: alt is symbolic (<DEL>, <INS>, <DUP>), or breakend ([chr etc.)\n"
"                or length difference >=50.\n"
"  If --append-info, prints original columns + updated INFO. Otherwise prints\n"
"  'CHROM POS ID REF ALT Classification' as TSV.\n\n"
"Examples:\n"
"  1) TSV classification:\n"
"     VCFX_variant_classifier < input.vcf > classified.tsv\n"
"  2) Modify INFO in output VCF:\n"
"     VCFX_variant_classifier --append-info < input.vcf > annotated.vcf\n";
}

std::vector<std::string> VCFXVariantClassifier::split(const std::string &s, char delim) const{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string t;
    while(std::getline(ss,t,delim)){
        out.push_back(t);
    }
    return out;
}

bool VCFXVariantClassifier::isStructuralAllele(const std::string &alt) const{
    // If alt is symbolic: <DEL>, <INS>, ...
    if(!alt.empty() && alt.front()=='<' && alt.back()=='>') return true;
    // If alt has breakend notation: [chr..., ]chr...
    if(alt.find('[')!=std::string::npos || alt.find(']')!= std::string::npos){
        return true;
    }
    // If alt length or ref length > 50 => structural
    // We'll handle that in classifyAllele by comparing lengths
    return false;
}

VariantType VCFXVariantClassifier::classifyAllele(const std::string &ref, const std::string &alt) const{
    // check structural
    if(isStructuralAllele(alt)){
        return VariantType::STRUCTURAL;
    }
    
    // if length difference >=50 => structural
    if(std::abs((int)ref.size() - (int)alt.size()) >= 50){
        return VariantType::STRUCTURAL;
    }
    
    // If REF and ALT are identical, this should be UNKNOWN
    if(ref == alt) {
        return VariantType::UNKNOWN;
    }
    
    // if both single base
    if(ref.size()==1 && alt.size()==1 &&
       std::isalpha((unsigned char)ref[0]) && std::isalpha((unsigned char)alt[0])){
        return VariantType::SNP;
    }
    
    // if length differs => small INDEL
    if(ref.size()!= alt.size()){
        // Additional check for large variants
        if(ref.size() >= 40 || alt.size() >= 40){
            return VariantType::STRUCTURAL;
        }
        return VariantType::INDEL;
    }
    
    // same length >1 => MNV
    if(ref.size()>1){
        return VariantType::MNV;
    }
    return VariantType::UNKNOWN;
}

VariantType VCFXVariantClassifier::classifyVariant(const std::string &ref,
                                                  const std::vector<std::string> &alts) const
{
    bool anyStruct= false, anyMNV=false, anyIndel=false, anySNP=false;
    for(auto &a: alts){
        VariantType vt= classifyAllele(ref,a);
        if(vt== VariantType::STRUCTURAL) anyStruct= true;
        else if(vt== VariantType::MNV) anyMNV= true;
        else if(vt== VariantType::INDEL) anyIndel= true;
        else if(vt== VariantType::SNP) anySNP= true;
    }
    // priority: STRUCTURAL> MNV> INDEL> SNP> UNKNOWN
    if(anyStruct) return VariantType::STRUCTURAL;
    if(anyMNV) return VariantType::MNV;
    if(anyIndel) return VariantType::INDEL;
    if(anySNP) return VariantType::SNP;
    return VariantType::UNKNOWN;
}

std::string VCFXVariantClassifier::typeToStr(VariantType t) const{
    switch(t){
        case VariantType::SNP: return "SNP";
        case VariantType::INDEL: return "INDEL";
        case VariantType::MNV: return "MNV";
        case VariantType::STRUCTURAL: return "STRUCTURAL";
        default: return "UNKNOWN";
    }
}

// If we want to keep the line as a valid VCF, we parse the line, classify, then append to INFO "VCF_CLASS=XYZ"
std::string VCFXVariantClassifier::appendClassification(const std::string &line){
    // split by tab
    auto fields= split(line, '\t');
    if(fields.size()<8){
        // skip
        return line;
    }
    // parse alt
    auto alts= split(fields[4], ',');
    VariantType vt= classifyVariant(fields[3], alts);
    std::string cstr= typeToStr(vt);

    // append to INFO
    if(fields[7]=="."|| fields[7].empty()){
        fields[7]= "VCF_CLASS="+ cstr;
    } else {
        // add semicolon if needed
        if(!fields[7].empty() && fields[7].back()!=';'){
            fields[7]+= ";";
        }
        fields[7]+= "VCF_CLASS="+ cstr;
    }
    // reconstruct line
    std::ostringstream oss;
    for(size_t i=0; i<fields.size(); i++){
        if(i>0) oss<<"\t";
        oss<< fields[i];
    }
    return oss.str();
}

void VCFXVariantClassifier::classifyStream(std::istream &in, std::ostream &out){
    bool foundChromHeader= false;
    std::string line;
    if(appendInfo){
        // we produce a valid VCF. We'll pass all # lines unmodified
        while(true){
            if(!std::getline(in,line)) break;
            if(line.empty()){
                out<< line<<"\n";
                continue;
            }
            if(line[0]=='#'){
                out<< line<<"\n";
                if(line.rfind("#CHROM",0)==0) foundChromHeader= true;
                continue;
            }
            // data line
            if(!foundChromHeader){
                std::cerr<<"Warning: data line encountered before #CHROM => skipping.\n";
                continue;
            }
            // parse, classify, re-output
            auto fields= split(line,'\t');
            if(fields.size()<8){
                std::cerr<<"Warning: skipping line <8 columns.\n";
                continue;
            }
            // append classification in INFO
            std::string newLine= appendClassification(line);
            out<< newLine<<"\n";
        }
    } else {
        // produce a TSV: CHROM POS ID REF ALT CLASS
        // Output the TSV header first
        out << "CHROM\tPOS\tID\tREF\tALT\tClassification\n";
        
        // skip all # lines
        while(true){
            if(!std::getline(in,line)) break;
            if(line.empty()) continue;
            if(line[0]=='#'){
                if(line.rfind("#CHROM",0)==0) foundChromHeader= true;
                continue;
            }
            if(!foundChromHeader){
                std::cerr<<"Warning: data line before #CHROM => skipping.\n";
                continue;
            }
            // parse
            auto fields= split(line,'\t');
            if(fields.size()<8){
                std::cerr<<"Warning: skipping line <8 columns.\n";
                continue;
            }
            
            // Additional validation for malformed input
            // 1. Validate CHROM format (should start with "chr")
            if (fields[0].substr(0, 3) != "chr") {
                std::cerr<<"Warning: invalid chromosome format => skipping.\n";
                continue;
            }
            
            // 2. Validate POS is numeric
            bool validPos = true;
            for (char c : fields[1]) {
                if (!std::isdigit(c)) {
                    validPos = false;
                    break;
                }
            }
            if (!validPos) {
                std::cerr<<"Warning: position is not numeric => skipping.\n";
                continue;
            }
            
            // 3. Validate REF and ALT are not empty
            if (fields[3].empty() || fields[4].empty()) {
                std::cerr<<"Warning: REF or ALT is empty => skipping.\n";
                continue;
            }
            
            // 4. Validate REF contains only valid bases
            bool validRef = true;
            for (char c : fields[3]) {
                if (!std::isalpha(c)) {
                    validRef = false;
                    break;
                }
            }
            if (!validRef) {
                std::cerr<<"Warning: REF contains non-alphabetic characters => skipping.\n";
                continue;
            }
            
            // 5. Validate ALT format
            // Allow special cases like <DEL>, but disallow trailing commas
            if (fields[4].back() == ',') {
                std::cerr<<"Warning: ALT ends with a comma => skipping.\n";
                continue;
            }
            
            // alt => split by comma
            auto altList= split(fields[4], ',');
            VariantType vt= classifyVariant(fields[3], altList);
            // build alt string
            std::string altJoined= fields[4];
            // Output the classification line
            out<< fields[0]<<"\t"<< fields[1]<<"\t"<< fields[2]
               <<"\t"<< fields[3]<<"\t"<< altJoined<<"\t"<< typeToStr(vt)<<"\n";
        }
    }
}

int main(int argc, char* argv[]){
    VCFXVariantClassifier app;
    return app.run(argc, argv);
}


./VCFX_variant_classifier/VCFX_variant_classifier.h
#ifndef VCFX_VARIANT_CLASSIFIER_H
#define VCFX_VARIANT_CLASSIFIER_H

#include <string>
#include <vector>

enum class VariantType {
    SNP,
    INDEL,
    MNV,
    STRUCTURAL,
    UNKNOWN
};

class VCFXVariantClassifier {
public:
    int run(int argc, char* argv[]);

private:
    // Show usage
    void displayHelp();

    // Identify the user’s requested output mode
    bool appendInfo = false; // if true, output a fully valid VCF with classification appended to INFO
    // otherwise produce a TSV with columns: CHROM POS ID REF ALT Classification

    // The main method that reads lines from input, classifies, and writes output
    void classifyStream(std::istream &in, std::ostream &out);

    // Detect alt allele as structural if it is symbolic (<DEL>) or breakend notation ([chr or ]chr) or length difference >=50
    bool isStructuralAllele(const std::string &alt) const;

    // Classify a single (ref, alt)
    VariantType classifyAllele(const std::string &ref, const std::string &alt) const;

    // From the set of alt alleles, find the final classification
    VariantType classifyVariant(const std::string &ref, const std::vector<std::string> &alts) const;

    // Stringify the variant type
    std::string typeToStr(VariantType t) const;

    // If we are in append‐info mode, we parse line into columns, parse alt as multiple, classify, then append e.g. VCF_CLASS=xxx
    std::string appendClassification(const std::string &line);

    // Splitting helper
    std::vector<std::string> split(const std::string &s, char delim) const;
};

#endif


./VCFX_variant_counter/VCFX_variant_counter.cpp
#include "VCFX_variant_counter.h"
#include <getopt.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

void VCFXVariantCounter::displayHelp(){
    std::cout <<
"VCFX_variant_counter: Counts the total number of valid variants in a VCF.\n\n"
"Usage:\n"
"  VCFX_variant_counter [options] < input.vcf\n\n"
"Options:\n"
"  -h, --help        Show this help.\n"
"  -s, --strict      Fail on any data line with <8 columns.\n\n"
"Description:\n"
"  Reads a VCF from stdin, ignores all header lines (#). For each data line,\n"
"  we check if it has >=8 columns; if it does, we count it; if fewer columns:\n"
"   * if --strict => we exit with error,\n"
"   * otherwise => we skip with a warning.\n"
"  Finally, we print 'Total Variants: X'.\n\n"
"Example:\n"
"  VCFX_variant_counter < input.vcf\n"
"  VCFX_variant_counter --strict < input.vcf\n";
}

int VCFXVariantCounter::run(int argc, char* argv[]){
    bool showHelp=false;
    static struct option long_opts[]={
        {"help", no_argument, 0,'h'},
        {"strict", no_argument, 0,'s'},
        {0,0,0,0}
    };
    
    // Check for options/flags
    if (argc > 1) {
        while(true){
            int c= ::getopt_long(argc, argv,"hs", long_opts,nullptr);
            if(c==-1) break;
            switch(c){
                case 'h':
                    showHelp=true;
                    break;
                case 's':
                    strictMode= true;
                    break;
                default:
                    showHelp=true;
            }
        }
    }
    
    if(showHelp){
        displayHelp();
        return 0;
    }
    
    int total= countVariants(std::cin);
    if(total<0){
        // indicates an error if strict
        return 1;
    }
    std::cout<<"Total Variants: "<< total <<"\n";
    return 0;
}

int VCFXVariantCounter::countVariants(std::istream &in){
    int count=0;
    int lineNumber=0;
    std::string line;
    while(true){
        if(!std::getline(in,line)) break;
        lineNumber++;
        if(line.empty()) continue;
        if(line[0]=='#') continue; 
        // parse columns
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string col;
            while(std::getline(ss,col,'\t')){
                fields.push_back(col);
            }
        }
        if(fields.size()<8){
            if(strictMode){
                std::cerr<<"Error: line "<< lineNumber <<" has <8 columns.\n";
                return -1; // indicates error
            } else {
                std::cerr<<"Warning: skipping line "<<lineNumber<<" with <8 columns.\n";
                continue;
            }
        }
        // if we get here => count it
        count++;
    }
    return count;
}

int main(int argc, char* argv[]) {
    VCFXVariantCounter app;
    return app.run(argc, argv);
}


./VCFX_variant_counter/VCFX_variant_counter.h
#ifndef VCFX_VARIANT_COUNTER_H
#define VCFX_VARIANT_COUNTER_H

#include <iostream>
#include <string>

class VCFXVariantCounter {
public:
    int run(int argc, char* argv[]);

private:
    // If true, any line with <8 columns is a fatal error
    bool strictMode = false;

    // Show usage
    void displayHelp();

    // The actual counting function
    int countVariants(std::istream &in);

};

#endif


