# Quick Start Guide

This guide provides a quick introduction to using the VCFX toolkit with practical examples.

## Overview

VCFX consists of multiple small command-line tools that are designed to be combined in pipelines. Each tool follows these principles:

- Reads from standard input and writes to standard output
- Performs a single, focused task
- Can be combined with other tools via pipes
- Has a consistent command-line interface

## Basic Usage Pattern

Most VCFX tools follow this basic pattern:

```bash
VCFX_tool_name [options] < input.vcf > output.vcf
```

For example:

```bash
VCFX_variant_classifier < input.vcf > classified.vcf
```

## Common Examples

Here are some common use cases for VCFX tools:

### Example 1: Basic Filtering

Filter for high-quality SNPs:

```bash
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep 'VCF_CLASS=SNP' | \
  VCFX_phred_filter --phred-filter 30 > high_quality_snps.vcf
```

### Example 2: Population Analysis

Extract European samples and calculate allele frequencies:

```bash
cat input.vcf | \
  VCFX_population_filter --population EUR --pop-map populations.txt | \
  VCFX_allele_freq_calc > eur_frequencies.tsv
```

### Example 3: Data Transformation

Normalize indels and split multiallelic variants:

```bash
cat input.vcf | \
  VCFX_indel_normalizer | \
  VCFX_multiallelic_splitter > normalized_biallelic.vcf
```

### Example 4: Quality Control

Check concordance between two VCF files:

```bash
VCFX_concordance_checker --vcf1 sample1.vcf --vcf2 sample2.vcf > concordance_report.tsv
```

## Working with Compressed Files

Most VCFX tools don't directly support compressed input/output. Use standard Unix tools:

```bash
# Reading from compressed files
zcat input.vcf.gz | VCFX_tool_name > output.vcf

# Writing to compressed files
VCFX_tool_name < input.vcf | gzip > output.vcf.gz
```

## Getting Help

All VCFX tools provide help information:

```bash
VCFX_tool_name --help
```

This will display the tool's purpose, options, and usage examples.

## Tool Categories

VCFX tools are categorized by their function:

### Data Analysis
Tools for extracting information from VCF files (e.g., `VCFX_allele_freq_calc`)

### Data Filtering
Tools for selecting variants based on criteria (e.g., `VCFX_phred_filter`)

### Data Transformation
Tools for converting or reformatting VCF data (e.g., `VCFX_indel_normalizer`)

### Quality Control
Tools for validating and checking data quality (e.g., `VCFX_validator`)

### File Management
Tools for handling VCF files (e.g., `VCFX_indexer`)

## Common Workflows

Here are some common workflows that combine multiple VCFX tools:

### Variant QC Pipeline

```bash
cat input.vcf | \
  VCFX_validator | \
  VCFX_variant_classifier --append-info | \
  VCFX_missing_detector | \
  grep -v 'MISSING_GENOTYPES=1' | \
  VCFX_phred_filter --phred-filter 20 > qc_passed.vcf
```

### Sample Comparison

```bash
# Extract common samples
VCFX_sample_extractor --samples SAMPLE1,SAMPLE2 < input1.vcf > samples1.vcf
VCFX_sample_extractor --samples SAMPLE1,SAMPLE2 < input2.vcf > samples2.vcf

# Check concordance
VCFX_concordance_checker --vcf1 samples1.vcf --vcf2 samples2.vcf > concordance.tsv
```

### Population Structure Analysis

```bash
# Extract population-specific VCFs
cat input.vcf | VCFX_population_filter --population EUR --pop-map pop_map.txt > eur.vcf
cat input.vcf | VCFX_population_filter --population AFR --pop-map pop_map.txt > afr.vcf

# Calculate allele frequencies for each population
cat eur.vcf | VCFX_allele_freq_calc > eur_afs.tsv
cat afr.vcf | VCFX_allele_freq_calc > afr_afs.tsv
```

## Next Steps

After becoming familiar with the basic usage of VCFX tools, you can:

1. Explore the [complete tool documentation](tools_overview.md) for details on each tool
2. Check the [installation guide](installation.md) if you need to install additional tools
3. Browse through the example VCF files in the repository to practice

For more complex workflows and advanced examples, refer to the individual tool documentation pages. 