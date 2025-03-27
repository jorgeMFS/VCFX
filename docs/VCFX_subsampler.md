# VCFX_subsampler

## Overview
`VCFX_subsampler` is a utility tool for randomly selecting a subset of variants from a VCF file while preserving all header information. It implements reservoir sampling to efficiently select a random sample of specified size from the input data.

## Usage
```bash
VCFX_subsampler [OPTIONS] < input.vcf > output.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit |
| `-s`, `--subsample <N>` | **Required**: Number of data lines (variants) to keep |
| `--seed <INT>` | Use a specific random seed for reproducible sampling |

## Description
`VCFX_subsampler` processes a VCF file to create a random subset by:

1. Reading the VCF file from standard input
2. Preserving all header lines (starting with '#') without modification
3. Using reservoir sampling to randomly select N data lines
4. Writing the header lines and the selected subset to standard output

This tool is particularly useful for:
- Creating smaller, representative datasets for testing and development
- Reducing the size of large VCF files for easier handling
- Creating balanced datasets with a controlled number of variants
- Generating multiple random subsets for statistical analyses

## Sampling Method

### Reservoir Sampling
The tool uses reservoir sampling, an algorithm that efficiently selects a random sample of k items from a population of unknown size n in a single pass. The algorithm:

1. Stores the first k items in a "reservoir"
2. For each subsequent item (k+1 to n):
   - Generates a random number j between 0 and the current item number
   - If j < k, replaces the jth item in the reservoir with the current item

This ensures that each variant has an equal probability of being selected in the final sample, regardless of its position in the file.

## Examples

### Basic Usage
Select 1000 random variants from a VCF file:
```bash
VCFX_subsampler --subsample 1000 < large.vcf > sample.vcf
```

### Reproducible Sampling
Select 1000 random variants with a fixed seed for reproducibility:
```bash
VCFX_subsampler --subsample 1000 --seed 12345 < large.vcf > reproducible_sample.vcf
```

### Creating Multiple Samples
Generate multiple random samples of the same size:
```bash
VCFX_subsampler --subsample 500 --seed 1 < input.vcf > sample1.vcf
VCFX_subsampler --subsample 500 --seed 2 < input.vcf > sample2.vcf
VCFX_subsampler --subsample 500 --seed 3 < input.vcf > sample3.vcf
```

## Special Case Handling

### Small Input Files
- If the input file contains fewer variants than the requested sample size, all variants are included in the output

### Malformed Lines
- Lines with fewer than 8 columns are skipped with a warning
- Empty lines in the header section are preserved

### Missing Sample Size
- The `--subsample` option is required; the program will exit with an error if not specified
- Sample size must be a positive integer; specifying zero or a negative number will result in an error

### Random Seed
- If no seed is specified, the current time is used as the seed, resulting in different samples each time
- When a seed is specified, the sampling is deterministic and reproducible

## Performance Considerations
- Reservoir sampling requires only a single pass through the data
- Memory usage is proportional to the requested sample size, not the input file size
- Very large requested sample sizes may require significant memory
- Processing time scales linearly with the size of the input file

## Limitations
- Cannot sample based on specific criteria or filters
- No preservation of related variants (e.g., those in linkage disequilibrium)
- Does not maintain the order of variants from the original file
- No special handling for multi-sample VCF files (sampling is done at the variant level, not the sample level)
- Cannot handle compressed (gzipped) VCF files directly 