# VCFX Toolkit Documentation

VCFX is a comprehensive toolkit for VCF (Variant Call Format) file manipulation, analysis, and transformation. Each tool in the VCFX suite is designed to perform a specific operation on VCF files, following the Unix philosophy of doing one thing well.

## Key Features

- Process VCF files from standard input and output to standard output
- Easily compose operations by piping tools together
- Efficient processing, suitable for large-scale genomic data
- Consistent interface across all tools

## Installation

VCFX tools are built using CMake. To build the entire toolkit:

```bash
mkdir -p build
cd build
cmake ..
make
```

To build a specific tool:

```bash
make VCFX_tool_name
```

## Tool Categories

The VCFX toolkit includes tools in the following categories:

### Data Analysis

- [VCFX_allele_freq_calc](VCFX_allele_freq_calc.md) - Calculate allele frequencies
- [VCFX_variant_classifier](VCFX_variant_classifier.md) - Classify variants into SNP, INDEL, MNV, or STRUCTURAL
- [VCFX_inbreeding_calculator](VCFX_inbreeding_calculator.md) - Calculate inbreeding coefficients
- [VCFX_dosage_calculator](VCFX_dosage_calculator.md) - Calculate allele dosage from genotypes
- [VCFX_hwe_tester](VCFX_hwe_tester.md) - Perform Hardy-Weinberg Equilibrium tests
- [VCFX_distance_calculator](VCFX_distance_calculator.md) - Calculate genetic distances between samples
- [VCFX_allele_counter](VCFX_allele_counter.md) - Count alleles across samples
- [VCFX_allele_balance_calc](VCFX_allele_balance_calc.md) - Calculate allele balance metrics
- [VCFX_variant_counter](VCFX_variant_counter.md) - Count variants with specific properties
- [VCFX_ancestry_inferrer](VCFX_ancestry_inferrer.md) - Infer ancestry from genotypes using population frequencies
- [VCFX_ancestry_assigner](VCFX_ancestry_assigner.md) - Assign samples to ancestral populations

### Data Filtering

- [VCFX_phase_checker](VCFX_phase_checker.md) - Filter variants to keep only fully phased genotypes
- [VCFX_phred_filter](VCFX_phred_filter.md) - Filter variants based on Phred-scaled quality scores
- [VCFX_record_filter](VCFX_record_filter.md) - Filter variants based on various VCF fields
- [VCFX_gl_filter](VCFX_gl_filter.md) - Filter variants based on genotype likelihood scores
- [VCFX_allele_balance_filter](VCFX_allele_balance_filter.md) - Filter variants based on allele balance
- [VCFX_population_filter](VCFX_population_filter.md) - Filter variants based on population-level metrics
- [VCFX_probability_filter](VCFX_probability_filter.md) - Filter variants based on probabilistic metrics
- [VCFX_nonref_filter](VCFX_nonref_filter.md) - Filter variants to keep only those with non-reference alleles
- [VCFX_impact_filter](VCFX_impact_filter.md) - Filter variants based on predicted impact

### Data Transformation

- [VCFX_multiallelic_splitter](VCFX_multiallelic_splitter.md) - Split multiallelic variants into biallelic records
- [VCFX_sample_extractor](VCFX_sample_extractor.md) - Extract specific samples from a VCF file
- [VCFX_position_subsetter](VCFX_position_subsetter.md) - Subset variants by position
- [VCFX_format_converter](VCFX_format_converter.md) - Convert between VCF and other formats
- [VCFX_genotype_query](VCFX_genotype_query.md) - Filter variants based on genotype patterns
- [VCFX_indel_normalizer](VCFX_indel_normalizer.md) - Normalize indel representations
- [VCFX_sv_handler](VCFX_sv_handler.md) - Process structural variant records
- [VCFX_fasta_converter](VCFX_fasta_converter.md) - Convert VCF files to FASTA format
- [VCFX_sorter](VCFX_sorter.md) - Sort VCF records by position
- [VCFX_af_subsetter](VCFX_af_subsetter.md) - Subset variants by allele frequency

### Quality Control

- [VCFX_concordance_checker](VCFX_concordance_checker.md) - Check concordance between VCF files
- [VCFX_missing_detector](VCFX_missing_detector.md) - Detect and report missing data
- [VCFX_outlier_detector](VCFX_outlier_detector.md) - Detect outlier samples
- [VCFX_alignment_checker](VCFX_alignment_checker.md) - Check alignment quality
- [VCFX_cross_sample_concordance](VCFX_cross_sample_concordance.md) - Check concordance between samples
- [VCFX_validator](VCFX_validator.md) - Validate VCF format compliance

### File Management

- [VCFX_indexer](VCFX_indexer.md) - Create an index file for random access
- [VCFX_file_splitter](VCFX_file_splitter.md) - Split VCF files into smaller chunks
- [VCFX_compressor](VCFX_compressor.md) - Compress VCF files efficiently
- [VCFX_diff_tool](VCFX_diff_tool.md) - Compare two VCF files

### Annotation and Reporting

- [VCFX_custom_annotator](VCFX_custom_annotator.md) - Add custom annotations to VCF files
- [VCFX_info_summarizer](VCFX_info_summarizer.md) - Summarize INFO fields
- [VCFX_header_parser](VCFX_header_parser.md) - Extract information from VCF headers
- [VCFX_annotation_extractor](VCFX_annotation_extractor.md) - Extract specific annotations
- [VCFX_ref_comparator](VCFX_ref_comparator.md) - Compare variants against a reference genome
- [VCFX_field_extractor](VCFX_field_extractor.md) - Extract specific fields from VCF records
- [VCFX_info_aggregator](VCFX_info_aggregator.md) - Aggregate numeric INFO field values

### Data Processing

- [VCFX_missing_data_handler](VCFX_missing_data_handler.md) - Handle missing genotype data
- [VCFX_quality_adjuster](VCFX_quality_adjuster.md) - Apply transformations to quality scores
- [VCFX_haplotype_phaser](VCFX_haplotype_phaser.md) - Phase genotypes into haplotypes
- [VCFX_haplotype_extractor](VCFX_haplotype_extractor.md) - Extract haplotypes from phased genotypes
- [VCFX_duplicate_remover](VCFX_duplicate_remover.md) - Remove duplicate variants

## Common Usage Patterns

VCFX tools are designed to be used in pipelines. Here are some common usage patterns:

### Basic Filtering and Analysis

```bash
# Extract phased variants, filter by quality, and calculate allele frequencies
cat input.vcf | \
  VCFX_phase_checker | \
  VCFX_phred_filter --min-qual 30 | \
  VCFX_allele_freq_calc > result.tsv
```

### Sample Selection and Comparison

```bash
# Extract samples and check concordance
cat input.vcf | \
  VCFX_sample_extractor --samples SAMPLE1,SAMPLE2 > samples.vcf

cat samples.vcf reference.vcf | \
  VCFX_concordance_checker > concordance_report.tsv
```

### Variant Classification and Filtering

```bash
# Classify variants and filter for specific types
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep -P "VCF_CLASS=SNP|VCF_CLASS=INDEL" > filtered_variants.vcf
```

### Quality Control Pipeline

```bash
# Complete QC pipeline
cat input.vcf | \
  VCFX_validator | \
  VCFX_missing_detector | \
  VCFX_hwe_tester | \
  VCFX_inbreeding_calculator > qc_results.txt
```

### Population Genetics Analysis

```bash
# Perform population structure analysis
cat input.vcf | \
  VCFX_ancestry_inferrer --frequency pop_freqs.txt > sample_ancestry.txt

# Compare with known assignments
cat input.vcf | \
  VCFX_ancestry_assigner --assign-ancestry reference_freqs.tsv > ancestry_check.txt
```

## Contributing

VCFX is an open-source project. Contributions are welcome. Please follow these guidelines:

1. Follow the existing code style
2. Add comprehensive tests for any new functionality
3. Update documentation to reflect changes
4. Submit a pull request with a clear description of changes

## License

VCFX is licensed under [LICENSE TYPE]. See the LICENSE file for details. 