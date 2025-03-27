# VCFX: A Comprehensive VCF Manipulation Toolkit

VCFX is a collection of specialized command-line tools designed for efficient manipulation, analysis, and transformation of VCF (Variant Call Format) files used in genomic research and bioinformatics.

## What is VCFX?

VCFX follows the Unix philosophy of creating small, focused tools that do one thing well and can be combined together to form powerful workflows. Each tool in the VCFX suite is optimized for a specific VCF-related task, enabling researchers and bioinformaticians to:

- Extract specific information from VCF files
- Filter variants based on various criteria
- Transform VCF data into different formats
- Analyze genotypes and compute statistics
- Validate and check VCF file integrity
- Manipulate structural variants and complex records

## Getting Started

To begin using VCFX, first follow the [installation instructions](#installation) below, then explore the [tool categories](#tool-categories) to find the right components for your workflow.

### Installation

VCFX tools are built using CMake. To build the entire toolkit:

```bash
git clone https://github.com/jorgeMFS/VCFX.git
cd VCFX
mkdir -p build
cd build
cmake ..
make
```

To build a specific tool:

```bash
make VCFX_tool_name
```

### Basic Example

Here's a simple example of using VCFX to analyze variants:

```bash
# Calculate allele frequencies for SNPs only
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep 'VCF_CLASS=SNP' | \
  VCFX_allele_freq_calc > snp_frequencies.tsv
```

## Tool Categories

The VCFX toolkit includes tools in the following categories:

### Data Analysis

Tools for extracting statistical information and insights from variant data:

- [VCFX_allele_freq_calc](VCFX_allele_freq_calc.md) - Calculate allele frequencies
- [VCFX_variant_classifier](VCFX_variant_classifier.md) - Classify variants into SNP, INDEL, MNV, or STRUCTURAL
- [VCFX_inbreeding_calculator](VCFX_inbreeding_calculator.md) - Calculate inbreeding coefficients
- [VCFX_dosage_calculator](VCFX_dosage_calculator.md) - Calculate allele dosage from genotypes
- [View all analysis tools...](tools_overview.md#data-analysis)

### Data Filtering

Tools for selecting variants based on specific criteria:

- [VCFX_phase_checker](VCFX_phase_checker.md) - Filter variants to keep only fully phased genotypes
- [VCFX_phred_filter](VCFX_phred_filter.md) - Filter variants based on Phred-scaled quality scores
- [VCFX_record_filter](VCFX_record_filter.md) - Filter variants based on various VCF fields
- [View all filtering tools...](tools_overview.md#data-filtering)

### Data Transformation

Tools for converting or reformatting VCF data:

- [VCFX_multiallelic_splitter](VCFX_multiallelic_splitter.md) - Split multiallelic variants into biallelic records
- [VCFX_sample_extractor](VCFX_sample_extractor.md) - Extract specific samples from a VCF file
- [VCFX_indel_normalizer](VCFX_indel_normalizer.md) - Normalize indel representations
- [View all transformation tools...](tools_overview.md#data-transformation)

### Quality Control

Tools for validating and checking data quality:

- [VCFX_concordance_checker](VCFX_concordance_checker.md) - Check concordance between VCF files
- [VCFX_missing_detector](VCFX_missing_detector.md) - Detect and report missing data
- [VCFX_validator](VCFX_validator.md) - Validate VCF format compliance
- [View all quality control tools...](tools_overview.md#quality-control)

### File Management

Tools for handling VCF files:

- [VCFX_indexer](VCFX_indexer.md) - Create an index file for random access
- [VCFX_file_splitter](VCFX_file_splitter.md) - Split VCF files into smaller chunks
- [VCFX_compressor](VCFX_compressor.md) - Compress VCF files efficiently
- [View all file management tools...](tools_overview.md#file-management)

For a complete list of all tools and detailed usage examples, see the [full documentation](tools_overview.md).

## Who Should Use VCFX?

VCFX is designed for:

- **Bioinformaticians** working with genomic variant data
- **Researchers** analyzing VCF files from sequencing projects
- **Pipeline developers** creating reproducible genomic workflows
- **Data scientists** extracting information from genetic variants

## Key Features

- **Composability**: All tools work with standard input/output for easy pipeline integration
- **Efficiency**: Optimized for performance with large genomic datasets
- **Robustness**: Careful error handling and validation of VCF formatting
- **Flexibility**: Works with various VCF versions and extensions
- **Simplicity**: Clear, focused tools with consistent interfaces

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

See the [Tools Overview](tools_overview.md#common-usage-patterns) for more usage examples.

## Community and Support

- [GitHub Repository](https://github.com/jorgeMFS/VCFX)
- [Issue Tracker](https://github.com/jorgeMFS/VCFX/issues)
- [Contributing Guidelines](CONTRIBUTING.md)

## License

VCFX is available under MIT License. See the [LICENSE](LICENSE.md) file for details. 