# VCFX Tools Overview

VCFX is a collection of C/C++ tools for processing and analyzing VCF (Variant Call Format) files, with optional WebAssembly compatibility. Each tool is an independent command-line executable that can parse input from `stdin` and write to `stdout`, enabling flexible piping and integration into bioinformatics pipelines.

The suite also includes a convenience wrapper `vcfx` so you can run commands as `vcfx <subcommand>`. For example, `vcfx variant_counter` is equivalent to running `VCFX_variant_counter`. Use `vcfx --list` to see available subcommands. All individual `VCFX_*` binaries remain available if you prefer calling them directly.

## Tool Categories

### Data Analysis

These tools help extract statistical information and insights from variant data:

- [VCFX_allele_freq_calc](VCFX_allele_freq_calc.md) - Calculate allele frequencies
- [VCFX_variant_classifier](VCFX_variant_classifier.md) - Classify variants into SNP, INDEL, MNV, or STRUCTURAL
- [VCFX_inbreeding_calculator](VCFX_inbreeding_calculator.md) - Calculate inbreeding coefficients
- [VCFX_dosage_calculator](VCFX_dosage_calculator.md) - Calculate allele dosage from genotypes
- [VCFX_hwe_tester](VCFX_hwe_tester.md) - Test for Hardy-Weinberg equilibrium
- [VCFX_distance_calculator](VCFX_distance_calculator.md) - Calculate genetic distances between samples
- [VCFX_allele_counter](VCFX_allele_counter.md) - Count alleles in VCF files
- [VCFX_allele_balance_calc](VCFX_allele_balance_calc.md) - Calculate allele balance metrics
- [VCFX_variant_counter](VCFX_variant_counter.md) - Count variants in VCF files
- [VCFX_ancestry_inferrer](VCFX_ancestry_inferrer.md) - Infer ancestry from genetic data
- [VCFX_ancestry_assigner](VCFX_ancestry_assigner.md) - Assign ancestry to samples
- [VCFX_ld_calculator](VCFX_ld_calculator.md) - Calculate pairwise linkage disequilibrium (r²) between variants

### Data Filtering

Tools for selecting variants based on specific criteria:

- [VCFX_phase_checker](VCFX_phase_checker.md) - Filter variants to keep only fully phased genotypes
- [VCFX_phred_filter](VCFX_phred_filter.md) - Filter variants based on Phred-scaled quality scores
- [VCFX_record_filter](VCFX_record_filter.md) - Filter variants based on various VCF fields
- [VCFX_gl_filter](VCFX_gl_filter.md) - Filter variants based on genotype likelihoods
- [VCFX_allele_balance_filter](VCFX_allele_balance_filter.md) - Filter variants based on allele balance
- [VCFX_population_filter](VCFX_population_filter.md) - Filter variants based on population statistics
- [VCFX_probability_filter](VCFX_probability_filter.md) - Filter variants based on probability scores
- [VCFX_nonref_filter](VCFX_nonref_filter.md) - Filter to keep only non-reference variants
- [VCFX_impact_filter](VCFX_impact_filter.md) - Filter variants based on predicted impact
- [VCFX_phase_quality_filter](VCFX_phase_quality_filter.md) - Filter variants based on phasing quality scores
- [VCFX_region_subsampler](VCFX_region_subsampler.md) - Filter variants based on genomic regions

### Data Transformation

Tools for converting or reformatting VCF data:

- [VCFX_multiallelic_splitter](VCFX_multiallelic_splitter.md) - Split multiallelic variants into biallelic records
- [VCFX_sample_extractor](VCFX_sample_extractor.md) - Extract specific samples from a VCF file
- [VCFX_position_subsetter](VCFX_position_subsetter.md) - Extract variants at specific positions
- [VCFX_format_converter](VCFX_format_converter.md) - Convert VCF files to other formats
- [VCFX_genotype_query](VCFX_genotype_query.md) - Query specific genotype patterns
- [VCFX_indel_normalizer](VCFX_indel_normalizer.md) - Normalize indel representations
- [VCFX_sv_handler](VCFX_sv_handler.md) - Handle structural variants in VCF files
- [VCFX_fasta_converter](VCFX_fasta_converter.md) - Convert VCF files to FASTA format
- [VCFX_sorter](VCFX_sorter.md) - Sort VCF files by position
- [VCFX_af_subsetter](VCFX_af_subsetter.md) - Extract variants based on allele frequency
- [VCFX_reformatter](VCFX_reformatter.md) - Reformat VCF files for better readability

### Quality Control

Tools for validating and checking data quality:

- [VCFX_concordance_checker](VCFX_concordance_checker.md) - Check concordance between VCF files
- [VCFX_missing_detector](VCFX_missing_detector.md) - Detect and report missing data
- [VCFX_outlier_detector](VCFX_outlier_detector.md) - Detect outlier samples or variants
- [VCFX_alignment_checker](VCFX_alignment_checker.md) - Check alignment of variants
- [VCFX_cross_sample_concordance](VCFX_cross_sample_concordance.md) - Check concordance between samples
- [VCFX_validator](VCFX_validator.md) - Validate VCF format compliance

### File Management

Tools for handling VCF files:

- [VCFX_indexer](VCFX_indexer.md) - Create an index file for random access
- [VCFX_file_splitter](VCFX_file_splitter.md) - Split VCF files into smaller chunks
- [VCFX_compressor](VCFX_compressor.md) - Compress VCF files efficiently
- [VCFX_diff_tool](VCFX_diff_tool.md) - Find differences between VCF files
- [VCFX_subsampler](VCFX_subsampler.md) - Subsample variants from a VCF file
- [VCFX_duplicate_remover](VCFX_duplicate_remover.md) - Remove duplicate variants
- [VCFX_merger](VCFX_merger.md) - Merge multiple VCF files by position

### Annotation and Reporting

Tools for annotating and extracting information from VCF files:

- [VCFX_custom_annotator](VCFX_custom_annotator.md) - Add custom annotations to VCF files
- [VCFX_info_summarizer](VCFX_info_summarizer.md) - Summarize INFO fields in VCF files
- [VCFX_header_parser](VCFX_header_parser.md) - Parse and extract information from VCF headers
- [VCFX_annotation_extractor](VCFX_annotation_extractor.md) - Extract annotations from VCF files
- [VCFX_ref_comparator](VCFX_ref_comparator.md) - Compare variants against a reference genome
- [VCFX_field_extractor](VCFX_field_extractor.md) - Extract specific fields from VCF files
- [VCFX_info_aggregator](VCFX_info_aggregator.md) - Aggregate INFO fields across variants
- [VCFX_info_parser](VCFX_info_parser.md) - Parse INFO fields in VCF files
- [VCFX_metadata_summarizer](VCFX_metadata_summarizer.md) - Summarize key metadata from VCF files

### Data Processing

Tools for processing variants and samples:

- [VCFX_missing_data_handler](VCFX_missing_data_handler.md) - Handle missing data in VCF files
- [VCFX_quality_adjuster](VCFX_quality_adjuster.md) - Adjust quality scores in VCF files
- [VCFX_haplotype_phaser](VCFX_haplotype_phaser.md) - Phase haplotypes in VCF files
- [VCFX_haplotype_extractor](VCFX_haplotype_extractor.md) - Extract haplotype information

## Common Usage Patterns

VCFX tools are designed to be combined in pipelines. Here are some common usage patterns:

### Basic Filtering and Analysis

```bash
# Extract phased variants, filter by quality, and calculate allele frequencies
cat input.vcf | \
  VCFX_phase_checker | \
  VCFX_phred_filter --min-qual 30 | \
  VCFX_allele_freq_calc > result.tsv
```

### Variant Classification and Filtering

```bash
# Classify variants and filter for SNPs with high quality
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep 'VCF_CLASS=SNP' | \
  VCFX_phred_filter --min-qual 30 > high_quality_snps.vcf
```

### Sample Extraction and Comparison

```bash
# Extract samples and check concordance
cat input.vcf | VCFX_sample_extractor --samples SAMPLE1,SAMPLE2 > samples.vcf
cat samples.vcf reference.vcf | VCFX_concordance_checker > concordance_report.tsv
```

### Linkage Disequilibrium Analysis

```bash
# Calculate LD in a specific region after filtering for common variants
cat input.vcf | \
  VCFX_af_subsetter --min-af 0.05 | \
  VCFX_ld_calculator --region chr1:10000-20000 > ld_matrix.txt
```

### Normalization and Splitting

```bash
# Normalize indels and split multiallelic variants
cat input.vcf | \
  VCFX_indel_normalizer | \
  VCFX_multiallelic_splitter > normalized_biallelic.vcf
```

### Population Analysis

```bash
# Extract population-specific VCFs and calculate allele frequencies
cat input.vcf | VCFX_population_filter --population EUR --pop-map pop_map.txt > eur.vcf
cat eur.vcf | VCFX_allele_freq_calc > eur_afs.tsv
```

### Quality Control Pipeline

```bash
# Validate, classify, detect missing data, and filter by quality
cat input.vcf | \
  VCFX_validator | \
  VCFX_variant_classifier --append-info | \
  VCFX_missing_detector --max-missing 0.1 | \
  VCFX_phred_filter --min-qual 20 > qc_passed.vcf
``` 