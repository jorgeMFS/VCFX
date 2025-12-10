#!/bin/bash
# Comprehensive VCFX Performance Benchmark
# Tests all tools on 1000 Genomes Phase 3 data
#
# Usage: ./run_comprehensive_benchmark.sh [--tiny|--small|--full]
#
# Options:
#   --tiny     Run on tiny dataset (~100 variants, ~5KB) - quick test
#   --small    Run on small dataset (5K variants, ~50MB) - medium test
#   --full     Run on full dataset (427K variants, 4.3GB) - comprehensive (default)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="$ROOT_DIR/build/src"
DATA_DIR="$ROOT_DIR/benchmarks/data"
RESULTS_DIR="$ROOT_DIR/benchmarks/results"

# Dataset configurations
TINY_FILE="$DATA_DIR/chr20_tiny.vcf"
SMALL_FILE="$DATA_DIR/chr21_small.vcf"
FULL_FILE="$DATA_DIR/chr21.1kg.phase3.v5a.vcf"

# Parse arguments - default to full
DATASET_SIZE="full"
case "${1:-full}" in
    --tiny)
        DATASET_SIZE="tiny"
        DATA_FILE="$TINY_FILE"
        DEFAULT_TIMEOUT=30
        LONG_TIMEOUT=60
        ;;
    --small)
        DATASET_SIZE="small"
        DATA_FILE="$SMALL_FILE"
        DEFAULT_TIMEOUT=60
        LONG_TIMEOUT=300
        ;;
    --full|*)
        DATASET_SIZE="full"
        DATA_FILE="$FULL_FILE"
        DEFAULT_TIMEOUT=300
        LONG_TIMEOUT=3600
        ;;
esac

REF_FASTA="$DATA_DIR/chr21.fa"
POP_FREQ_FILE="$DATA_DIR/pop_frequencies.tsv"

# Create results directory
mkdir -p "$RESULTS_DIR"

# =============================================================================
# Generate helper files for tools that need additional data
# =============================================================================
echo "=== Preparing helper files ==="

# Download chr21 reference FASTA if not present
if [ ! -f "$REF_FASTA" ]; then
    echo "Downloading chr21 reference FASTA..."
    curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz" | gunzip > "$REF_FASTA" 2>/dev/null || {
        echo "Warning: Could not download reference FASTA, creating minimal version"
        echo ">21" > "$REF_FASTA"
        # Generate minimal reference from VCF positions (for benchmark purposes)
        awk -F'\t' '!/^#/ {print $2, $4}' "$DATA_FILE" | head -1000 | while read pos ref; do
            printf "%s" "$ref"
        done >> "$REF_FASTA"
        echo "" >> "$REF_FASTA"
    }
fi

# Generate population frequencies file from the VCF (extract AF from first 10K variants)
if [ ! -f "$POP_FREQ_FILE" ]; then
    echo "Generating population frequencies file..."
    awk -F'\t' 'BEGIN {OFS="\t"}
    !/^#/ && NR <= 10000 {
        # Parse AF from INFO field
        split($8, info, ";")
        af = "0.5"
        for (i in info) {
            if (info[i] ~ /^AF=/) {
                split(info[i], afval, "=")
                af = afval[2]
                gsub(/,.*/, "", af)  # Take first AF value for multiallelic
            }
        }
        # Assign population based on AF range (simulated)
        if (af < 0.1) pop = "EUR"
        else if (af < 0.2) pop = "AFR"
        else if (af < 0.3) pop = "EAS"
        else if (af < 0.4) pop = "SAS"
        else pop = "AMR"
        print $1, $2, $4, $5, pop, af
    }' "$DATA_FILE" > "$POP_FREQ_FILE"
fi
echo ""

# Results file - include dataset size in filename
RESULTS_CSV="$RESULTS_DIR/comprehensive_benchmark_${DATASET_SIZE}_$(date +%Y%m%d_%H%M%S).csv"
echo "tool,category,time_seconds,status,description" > "$RESULTS_CSV"

echo "=============================================="
echo "   VCFX Comprehensive Performance Benchmark"
echo "=============================================="
echo ""
echo "Dataset: $DATASET_SIZE"
echo "Test file: $DATA_FILE"
if [ -f "$DATA_FILE" ]; then
    echo "Size: $(ls -lh "$DATA_FILE" | awk '{print $5}')"
else
    echo "WARNING: Test file not found!"
    echo "Run 'make datasets' in benchmarks/ to download test data."
    exit 1
fi
echo "Default timeout: ${DEFAULT_TIMEOUT}s, Long timeout: ${LONG_TIMEOUT}s"
echo "Date: $(date)"
echo "Results: $RESULTS_CSV"
echo ""

# Function to benchmark a tool
benchmark() {
    local name="$1"
    local category="$2"
    local cmd="$3"
    local desc="$4"
    local timeout_secs="${5:-$DEFAULT_TIMEOUT}"

    echo -n "Testing $name... "

    # Run with timeout
    local start=$(python3 -c "import time; print(time.time())")
    if timeout $timeout_secs bash -c "$cmd" > /dev/null 2>&1; then
        local end=$(python3 -c "import time; print(time.time())")
        local elapsed=$(python3 -c "print(f'{$end - $start:.3f}')")
        echo "${elapsed}s"
        echo "$name,$category,$elapsed,success,$desc" >> "$RESULTS_CSV"
    else
        echo "FAILED/TIMEOUT"
        echo "$name,$category,0,failed,$desc" >> "$RESULTS_CSV"
    fi
}

# Warm-up run
echo "=== Warm-up (priming file cache) ==="
cat "$DATA_FILE" > /dev/null
echo ""

echo "=== Category 1: Basic I/O & Validation ==="
# variant_counter: mmap via positional arg
benchmark "variant_counter" "basic_io" \
    "$BUILD_DIR/VCFX_variant_counter/VCFX_variant_counter $DATA_FILE" \
    "Count variants (mmap)"

# validator: mmap via -i flag
benchmark "validator" "basic_io" \
    "$BUILD_DIR/VCFX_validator/VCFX_validator -i $DATA_FILE" \
    "Validate VCF (mmap)"

benchmark "header_parser" "basic_io" \
    "$BUILD_DIR/VCFX_header_parser/VCFX_header_parser < $DATA_FILE" \
    "Parse header"

benchmark "reformatter" "basic_io" \
    "$BUILD_DIR/VCFX_reformatter/VCFX_reformatter < $DATA_FILE" \
    "Reformat VCF"

echo ""
echo "=== Category 2: Filtering Tools ==="
# phred_filter: mmap via positional arg
benchmark "phred_filter" "filtering" \
    "$BUILD_DIR/VCFX_phred_filter/VCFX_phred_filter -p 30 $DATA_FILE" \
    "Filter QUAL>=30 (mmap)"

# record_filter: mmap via positional file argument
benchmark "record_filter" "filtering" \
    "$BUILD_DIR/VCFX_record_filter/VCFX_record_filter --filter 'FILTER==PASS' $DATA_FILE" \
    "Filter PASS only (mmap)"

benchmark "nonref_filter" "filtering" \
    "$BUILD_DIR/VCFX_nonref_filter/VCFX_nonref_filter -i $DATA_FILE" \
    "Non-reference filter (mmap)"

benchmark "gl_filter" "filtering" \
    "$BUILD_DIR/VCFX_gl_filter/VCFX_gl_filter --filter 'GQ>20' -i $DATA_FILE" \
    "GQ filter (GQ>20) (mmap)"

benchmark "allele_balance_filter" "filtering" \
    "$BUILD_DIR/VCFX_allele_balance_filter/VCFX_allele_balance_filter --min-ab 0.2 --max-ab 0.8 < $DATA_FILE" \
    "Allele balance filter"

benchmark "probability_filter" "filtering" \
    "$BUILD_DIR/VCFX_probability_filter/VCFX_probability_filter --min-gp 0.9 < $DATA_FILE" \
    "GP filter"

benchmark "phase_quality_filter" "filtering" \
    "$BUILD_DIR/VCFX_phase_quality_filter/VCFX_phase_quality_filter --min-pq 30 < $DATA_FILE" \
    "Phase quality filter"

benchmark "impact_filter" "filtering" \
    "$BUILD_DIR/VCFX_impact_filter/VCFX_impact_filter --filter-impact HIGH -I $DATA_FILE" \
    "Impact filter (HIGH) (mmap)"

echo ""
echo "=== Category 3: Analysis & Calculation Tools ==="
# allele_freq_calc: mmap via -i flag
benchmark "allele_freq_calc" "analysis" \
    "$BUILD_DIR/VCFX_allele_freq_calc/VCFX_allele_freq_calc -i $DATA_FILE" \
    "Allele frequencies (mmap)"

# variant_classifier: mmap via -i flag
benchmark "variant_classifier" "analysis" \
    "$BUILD_DIR/VCFX_variant_classifier/VCFX_variant_classifier -i $DATA_FILE" \
    "Classify variants (mmap)"

# NOTE: These O(variants × samples) tools - use mmap where available
# hwe_tester: mmap via -i flag
benchmark "hwe_tester" "analysis" \
    "$BUILD_DIR/VCFX_hwe_tester/VCFX_hwe_tester -i $DATA_FILE" \
    "HWE test (mmap)" $LONG_TIMEOUT

# allele_counter: mmap via -i flag, aggregate mode (-a) for efficient output
# Aggregate mode reduces output from O(variants × samples) to O(variants)
# Note: Use -l 100 to limit samples for benchmarking (full dataset = ~1B genotypes = ~1 hour)
benchmark "allele_counter" "analysis" \
    "$BUILD_DIR/VCFX_allele_counter/VCFX_allele_counter -a -l 100 -i $DATA_FILE" \
    "Count alleles (mmap, aggregate, 100 samples)" $LONG_TIMEOUT

# allele_balance_calc: mmap via -i flag
benchmark "allele_balance_calc" "analysis" \
    "$BUILD_DIR/VCFX_allele_balance_calc/VCFX_allele_balance_calc -i $DATA_FILE" \
    "Allele balance calc (mmap)" $LONG_TIMEOUT

# inbreeding_calculator: mmap via -i flag
benchmark "inbreeding_calculator" "analysis" \
    "$BUILD_DIR/VCFX_inbreeding_calculator/VCFX_inbreeding_calculator -i $DATA_FILE" \
    "Inbreeding coef (mmap)" $LONG_TIMEOUT

# dosage_calculator: mmap via -i flag
benchmark "dosage_calculator" "analysis" \
    "$BUILD_DIR/VCFX_dosage_calculator/VCFX_dosage_calculator -i $DATA_FILE" \
    "Dosage calc (mmap)" $LONG_TIMEOUT

# distance_calculator: mmap via -i flag
benchmark "distance_calculator" "analysis" \
    "$BUILD_DIR/VCFX_distance_calculator/VCFX_distance_calculator -i $DATA_FILE" \
    "Genetic distance (mmap)"

echo ""
echo "=== Category 4: Quality Control Tools ==="
# missing_detector: mmap via -i flag
benchmark "missing_detector" "qc" \
    "$BUILD_DIR/VCFX_missing_detector/VCFX_missing_detector -i $DATA_FILE" \
    "Detect missing (mmap)"

# Phase checker: mmap via -i flag
benchmark "phase_checker" "qc" \
    "$BUILD_DIR/VCFX_phase_checker/VCFX_phase_checker -i $DATA_FILE" \
    "Check phasing (mmap)" $LONG_TIMEOUT

benchmark "outlier_detector" "qc" \
    "$BUILD_DIR/VCFX_outlier_detector/VCFX_outlier_detector < $DATA_FILE" \
    "Detect outliers"

# cross_sample_concordance: mmap via -i flag
benchmark "cross_sample_concordance" "qc" \
    "$BUILD_DIR/VCFX_cross_sample_concordance/VCFX_cross_sample_concordance -i $DATA_FILE" \
    "Cross-sample concordance (mmap)"

# Alignment checker needs reference FASTA
benchmark "alignment_checker" "qc" \
    "$BUILD_DIR/VCFX_alignment_checker/VCFX_alignment_checker --alignment-discrepancy $DATA_FILE $REF_FASTA" \
    "Alignment check"

echo ""
echo "=== Category 5: Transformation Tools ==="
# multiallelic_splitter: mmap via -i flag
benchmark "multiallelic_splitter" "transformation" \
    "$BUILD_DIR/VCFX_multiallelic_splitter/VCFX_multiallelic_splitter -i $DATA_FILE" \
    "Split multiallelic (mmap)"

# indel_normalizer: mmap via -i flag
benchmark "indel_normalizer" "transformation" \
    "$BUILD_DIR/VCFX_indel_normalizer/VCFX_indel_normalizer -i $DATA_FILE" \
    "Normalize indels (mmap)" $LONG_TIMEOUT

# duplicate_remover: mmap via -i flag
benchmark "duplicate_remover" "transformation" \
    "$BUILD_DIR/VCFX_duplicate_remover/VCFX_duplicate_remover -i $DATA_FILE" \
    "Remove duplicates (mmap)"

# sorter: mmap via -i flag
benchmark "sorter" "transformation" \
    "$BUILD_DIR/VCFX_sorter/VCFX_sorter -i $DATA_FILE" \
    "Sort variants (mmap)" $LONG_TIMEOUT

benchmark "quality_adjuster" "transformation" \
    "$BUILD_DIR/VCFX_quality_adjuster/VCFX_quality_adjuster < $DATA_FILE" \
    "Adjust quality"

# missing_data_handler: mmap via -i flag
benchmark "missing_data_handler" "transformation" \
    "$BUILD_DIR/VCFX_missing_data_handler/VCFX_missing_data_handler --fill-missing -i $DATA_FILE" \
    "Handle missing (mmap)"

benchmark "sv_handler" "transformation" \
    "$BUILD_DIR/VCFX_sv_handler/VCFX_sv_handler < $DATA_FILE" \
    "Handle SVs"

echo ""
echo "=== Category 6: Extraction & Subsetting Tools ==="
benchmark "sample_extractor" "extraction" \
    "$BUILD_DIR/VCFX_sample_extractor/VCFX_sample_extractor --samples HG00096,HG00097,HG00099 -i $DATA_FILE" \
    "Extract samples (mmap)"

benchmark "field_extractor" "extraction" \
    "$BUILD_DIR/VCFX_field_extractor/VCFX_field_extractor --fields CHROM,POS,REF,ALT -i $DATA_FILE" \
    "Extract fields (mmap)"

benchmark "genotype_query" "extraction" \
    "$BUILD_DIR/VCFX_genotype_query/VCFX_genotype_query --genotype-query '0/1' -i $DATA_FILE" \
    "Query genotypes (mmap)"

benchmark "position_subsetter" "extraction" \
    "$BUILD_DIR/VCFX_position_subsetter/VCFX_position_subsetter --region 21:1-10000000 -i $DATA_FILE" \
    "Subset by region (mmap)"

# af_subsetter: mmap via -i flag
benchmark "af_subsetter" "extraction" \
    "$BUILD_DIR/VCFX_af_subsetter/VCFX_af_subsetter --af-filter 0.01-0.99 -i $DATA_FILE" \
    "Subset by AF (mmap)"

benchmark "subsampler" "extraction" \
    "$BUILD_DIR/VCFX_subsampler/VCFX_subsampler --fraction 0.1 < $DATA_FILE" \
    "Random subsample"

echo ""
echo "=== Category 7: Annotation & INFO Tools ==="
benchmark "info_parser" "annotation" \
    "$BUILD_DIR/VCFX_info_parser/VCFX_info_parser --info 'AC,AN,AF' -I $DATA_FILE" \
    "Parse INFO (mmap)"

benchmark "info_summarizer" "annotation" \
    "$BUILD_DIR/VCFX_info_summarizer/VCFX_info_summarizer --info 'AC,AN,AF' -I $DATA_FILE" \
    "Summarize INFO (mmap)"

# Info aggregator needs --aggregate-info argument
benchmark "info_aggregator" "annotation" \
    "$BUILD_DIR/VCFX_info_aggregator/VCFX_info_aggregator --aggregate-info 'AC,AN,AF' -i $DATA_FILE" \
    "Aggregate INFO (mmap)"

# metadata_summarizer: mmap via -i flag
benchmark "metadata_summarizer" "annotation" \
    "$BUILD_DIR/VCFX_metadata_summarizer/VCFX_metadata_summarizer -i $DATA_FILE" \
    "Summarize metadata (mmap)"

benchmark "annotation_extractor" "annotation" \
    "$BUILD_DIR/VCFX_annotation_extractor/VCFX_annotation_extractor --annotation-extract 'AF' -i $DATA_FILE" \
    "Extract annotations (mmap)"

echo ""
echo "=== Category 8: Haplotype Tools ==="
# haplotype_phaser: mmap via -i flag
benchmark "haplotype_phaser" "haplotype" \
    "$BUILD_DIR/VCFX_haplotype_phaser/VCFX_haplotype_phaser -i $DATA_FILE" \
    "Phase haplotypes (mmap)"

# haplotype_extractor: mmap via -i flag
benchmark "haplotype_extractor" "haplotype" \
    "$BUILD_DIR/VCFX_haplotype_extractor/VCFX_haplotype_extractor -i $DATA_FILE" \
    "Extract haplotypes (mmap)"

echo ""
echo "=== Category 9: Ancestry Tools ==="
# Ancestry inferrer needs population frequencies file
benchmark "ancestry_inferrer" "ancestry" \
    "$BUILD_DIR/VCFX_ancestry_inferrer/VCFX_ancestry_inferrer --frequency $POP_FREQ_FILE < $DATA_FILE" \
    "Infer ancestry"

# Ancestry assigner needs population frequencies file
benchmark "ancestry_assigner" "ancestry" \
    "$BUILD_DIR/VCFX_ancestry_assigner/VCFX_ancestry_assigner --pop-freqs $DATA_DIR/pop_frequencies.txt < $DATA_FILE" \
    "Assign ancestry"

echo ""
echo "=== Category 9b: Population Tools ==="
# Population filter needs a population map file
POP_MAP_FILE="$DATA_DIR/pop_map.tsv"
if [ ! -f "$POP_MAP_FILE" ]; then
    echo "Generating population map file from VCF samples..."
    # Extract sample names and assign populations based on sample ID patterns
    head -1000 "$DATA_FILE" | grep "^#CHROM" | cut -f10- | tr '\t' '\n' | awk '{
        # Assign populations based on sample patterns (simulated for 1000G samples)
        n = NR % 5
        if (n == 0) pop = "EUR"
        else if (n == 1) pop = "AFR"
        else if (n == 2) pop = "EAS"
        else if (n == 3) pop = "SAS"
        else pop = "AMR"
        print $0 "\t" pop
    }' > "$POP_MAP_FILE"
fi
benchmark "population_filter" "population" \
    "$BUILD_DIR/VCFX_population_filter/VCFX_population_filter --population EUR --pop-map $POP_MAP_FILE -i $DATA_FILE" \
    "Filter by population (mmap)"

echo ""
echo "=== Category 10: File Management Tools ==="
# indexer: uses file argument directly (mmap)
benchmark "indexer" "file_management" \
    "$BUILD_DIR/VCFX_indexer/VCFX_indexer $DATA_FILE" \
    "Create index (mmap)"

# Compressor needs --compress flag
benchmark "compressor" "file_management" \
    "$BUILD_DIR/VCFX_compressor/VCFX_compressor --compress < $DATA_FILE" \
    "Compress VCF"

# File splitter
benchmark "file_splitter" "file_management" \
    "$BUILD_DIR/VCFX_file_splitter/VCFX_file_splitter --chunks 10 < $DATA_FILE" \
    "Split into chunks"

# Merger (merge with itself as a test)
benchmark "merger" "file_management" \
    "$BUILD_DIR/VCFX_merger/VCFX_merger $DATA_FILE $DATA_FILE" \
    "Merge VCF files"

# Region subsampler - needs BED file with regions
REGION_BED="$DATA_DIR/regions.bed"
if [ ! -f "$REGION_BED" ]; then
    echo "Creating regions BED file..."
    echo -e "21\t0\t10000000\n21\t20000000\t30000000" > "$REGION_BED"
fi
benchmark "region_subsampler" "file_management" \
    "$BUILD_DIR/VCFX_region_subsampler/VCFX_region_subsampler --region-bed $REGION_BED -i $DATA_FILE" \
    "Subsample by region (mmap)"

echo ""
echo "=== Category 11: Format Conversion ==="
benchmark "format_converter_bed" "conversion" \
    "$BUILD_DIR/VCFX_format_converter/VCFX_format_converter --to-bed -i $DATA_FILE" \
    "Convert to BED (mmap)"

# fasta_converter: mmap via -i flag - O(variants × samples)
benchmark "fasta_converter" "conversion" \
    "$BUILD_DIR/VCFX_fasta_converter/VCFX_fasta_converter -i $DATA_FILE" \
    "Convert to FASTA (mmap)" $LONG_TIMEOUT

echo ""
echo "=== Category 12: Comparison & Diff Tools ==="
# diff_tool: uses --file1 and --file2 flags (mmap internally)
benchmark "diff_tool" "comparison" \
    "$BUILD_DIR/VCFX_diff_tool/VCFX_diff_tool --file1 $DATA_FILE --file2 $DATA_FILE --assume-sorted" \
    "Diff VCF files (mmap)"

# concordance_checker: mmap via -i flag
# Use first two samples from VCF - for 1000 Genomes data, these are HG00096 and HG00097
benchmark "concordance_checker" "comparison" \
    "$BUILD_DIR/VCFX_concordance_checker/VCFX_concordance_checker -i $DATA_FILE -s \"HG00096 HG00097\"" \
    "Check concordance (mmap)"

# Reference comparator needs reference FASTA
benchmark "ref_comparator" "comparison" \
    "$BUILD_DIR/VCFX_ref_comparator/VCFX_ref_comparator --reference $REF_FASTA < $DATA_FILE" \
    "Compare to reference"

echo ""
echo "=== Category 13: Advanced Analysis ==="
# ld_calculator: mmap via -i flag - computationally intensive
benchmark "ld_calculator" "advanced_analysis" \
    "$BUILD_DIR/VCFX_ld_calculator/VCFX_ld_calculator -i $DATA_FILE --window 1000" \
    "Calculate LD (mmap)" $LONG_TIMEOUT

# Custom annotator
benchmark "custom_annotator" "advanced_analysis" \
    "$BUILD_DIR/VCFX_custom_annotator/VCFX_custom_annotator --annotation 'TEST=1' < $DATA_FILE" \
    "Custom annotation"

echo ""
echo "=== Baseline Comparison: bcftools ==="
if command -v bcftools &> /dev/null; then
    benchmark "bcftools_view" "baseline" \
        "bcftools view -H $DATA_FILE | wc -l" \
        "bcftools count"

    benchmark "bcftools_query" "baseline" \
        "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $DATA_FILE" \
        "bcftools query"

    benchmark "bcftools_filter" "baseline" \
        "bcftools view -i 'QUAL>=30' $DATA_FILE" \
        "bcftools filter"
else
    echo "bcftools not installed - skipping baseline"
fi

echo ""
echo "=== Baseline Comparison: vcftools ==="
if command -v vcftools &> /dev/null; then
    benchmark "vcftools_freq" "baseline" \
        "vcftools --vcf $DATA_FILE --freq --stdout 2>/dev/null" \
        "vcftools freq"

    benchmark "vcftools_missing" "baseline" \
        "vcftools --vcf $DATA_FILE --missing-indv --stdout 2>/dev/null" \
        "vcftools missing"
else
    echo "vcftools not installed - skipping baseline"
fi

echo ""
echo "=============================================="
echo "   Benchmark Complete!"
echo "=============================================="
echo ""
echo "Results saved to: $RESULTS_CSV"
echo ""

# Print summary
echo "=== Summary by Category ==="
echo ""
python3 << PYSCRIPT
import csv
from collections import defaultdict

results = defaultdict(list)
with open("$RESULTS_CSV") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['status'] == 'success':
            results[row['category']].append(float(row['time_seconds']))

print(f"{'Category':<25} {'Tools':<8} {'Avg (s)':<10} {'Min (s)':<10} {'Max (s)':<10}")
print("-" * 63)
for cat in sorted(results.keys()):
    times = results[cat]
    print(f"{cat:<25} {len(times):<8} {sum(times)/len(times):<10.3f} {min(times):<10.3f} {max(times):<10.3f}")

total_time = sum(t for times in results.values() for t in times)
total_tools = sum(len(times) for times in results.values())
print("-" * 63)
print(f"{'TOTAL':<25} {total_tools:<8} {total_time:<10.3f}")
PYSCRIPT
