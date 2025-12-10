# VCFX_allele_counter

## Overview
`VCFX_allele_counter` counts the number of reference and alternate alleles in each sample for each variant in a VCF file. This tool provides a simple way to quantify allele occurrences across samples.

## Usage
```bash
VCFX_allele_counter [OPTIONS] [FILE]
```

## Options
| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapping for best performance) |
| `-t`, `--threads N` | Number of threads for parallel processing (default: auto-detect CPU cores) |
| `-s`, `--samples "Sample1 Sample2..."` | Optional. Specify sample names to calculate allele counts for (space-separated). If omitted, all samples are processed. |
| `-l`, `--limit-samples N` | Limit output to first N samples (useful for large cohorts) |
| `-a`, `--aggregate` | Output per-variant aggregate statistics instead of per-sample counts |
| `-z`, `--gzip` | Compress output with gzip |
| `-b`, `--binary` | Output in compact binary format (VCAC header) |
| `-q`, `--quiet` | Suppress informational messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description
`VCFX_allele_counter` processes a VCF file and counts reference and alternate alleles for each variant in each specified sample. The tool:

1. Reads a VCF file from standard input or file
2. Identifies sample columns from the VCF header
3. For each variant and each sample:
   - Extracts the genotype information
   - Counts reference alleles (0) and alternate alleles (non-0)
   - Outputs both counts in a tabular format
4. Outputs a tab-separated file with allele counts for each variant-sample combination

This tool is particularly useful for:
- Analyzing allele distribution across samples
- Quantifying the presence of specific alleles
- Preparing data for population genetics analyses
- Validating genotype calls across samples

## Output Formats

### Default Mode (Per-Sample Output)
The default output is a tab-separated values (TSV) file with the following columns:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome of the variant |
| POS | Position of the variant |
| ID | Variant identifier |
| REF | Reference allele |
| ALT | Alternate allele(s) |
| Sample | Sample name |
| Ref_Count | Number of reference alleles (0) in the sample's genotype |
| Alt_Count | Number of alternate alleles (non-0) in the sample's genotype |

### Aggregate Mode (`-a/--aggregate`)
Outputs per-variant summaries instead of per-sample counts, reducing output by ~2500x for large cohorts:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome of the variant |
| POS | Position of the variant |
| ID | Variant identifier |
| REF | Reference allele |
| ALT | Alternate allele(s) |
| Total_Ref | Sum of reference allele counts across all samples |
| Total_Alt | Sum of alternate allele counts across all samples |
| Sample_Count | Number of samples with valid genotypes |

### Binary Mode (`-b/--binary`)
Compact binary format for machine consumption (~25x smaller than text):
- 16-byte header: magic "VCAC", version (uint32), sample count (uint32), variant count (uint64)
- Per variant: CHROM, POS, ID, REF, ALT as null-terminated strings
- Per sample: 2 bytes (ref_count, alt_count as int8)

### Gzip Output (`-z/--gzip`)
Any output format can be compressed with gzip for ~10x size reduction.

## Examples

### File Input Mode (Recommended)
Use memory-mapped file I/O for best performance with multi-threaded processing:
```bash
VCFX_allele_counter -i input.vcf > allele_counts.tsv
VCFX_allele_counter -t 8 -i input.vcf > allele_counts.tsv  # Use 8 threads
VCFX_allele_counter input.vcf > allele_counts.tsv
```

### Basic Usage (Stdin)
Count alleles for all samples in a VCF file:
```bash
VCFX_allele_counter < input.vcf > allele_counts_all.tsv
```

### Specific Samples
Count alleles for specific samples:
```bash
VCFX_allele_counter --samples "SAMPLE1 SAMPLE2" -i input.vcf > allele_counts_subset.tsv
```

### Aggregate Mode (Recommended for Large Cohorts)
Get per-variant summaries instead of per-sample output:
```bash
VCFX_allele_counter --aggregate -i input.vcf > allele_summary.tsv
VCFX_allele_counter -a -i input.vcf > allele_summary.tsv  # Short form
```

### Limit Samples
Process only the first N samples:
```bash
VCFX_allele_counter --limit-samples 100 -i input.vcf > first_100_samples.tsv
VCFX_allele_counter -l 50 -a -i input.vcf > aggregate_50_samples.tsv
```

### Compressed Output
Gzip compress the output:
```bash
VCFX_allele_counter --gzip -i input.vcf > allele_counts.tsv.gz
VCFX_allele_counter -z -a -i input.vcf > aggregate.tsv.gz
```

### Binary Output
Compact binary format for downstream processing:
```bash
VCFX_allele_counter --binary -i input.vcf > allele_counts.bin
```

### Using with Other Tools
Process the output for further analysis:
```bash
VCFX_allele_counter -i input.vcf | awk -F'\t' '$8 > 0' > samples_with_alt_alleles.tsv
```

## Allele Counting Method

### Reference Alleles
The tool counts an allele as a reference allele when it has the value "0" in the genotype field. For example:
- In genotype "0/0", there are 2 reference alleles
- In genotype "0/1", there is 1 reference allele
- In genotype "1/2", there are 0 reference alleles

### Alternate Alleles
The tool counts an allele as an alternate allele when it has any non-zero numeric value in the genotype field. For example:
- In genotype "0/0", there are 0 alternate alleles
- In genotype "0/1", there is 1 alternate allele
- In genotype "1/2", there are 2 alternate alleles
- In genotype "1/1", there are 2 alternate alleles

### Handling Special Cases
- Missing genotypes (e.g., "./.", ".|."): No counts are recorded for these samples
- Partial missing (e.g., "0/."): Only the valid allele is counted
- Non-numeric alleles: These are skipped and not counted

## Handling Special Cases

### Missing Data
- Genotypes with missing values (`./.`, `.`) are skipped
- Partial missing genotypes only count the valid alleles present

### Multi-allelic Sites
- All non-reference alleles are counted as "alternate" regardless of their specific number
- For example, in a genotype "1/2", both alleles count as alternate alleles
- The tool does not differentiate between different alternate alleles

### Phased Genotypes
- Phasing information is ignored for allele counting
- Phased genotypes (e.g., "0|1") are treated the same as unphased (e.g., "0/1")

### Invalid Genotypes
- Non-numeric allele values are skipped
- Empty genotype fields are skipped

## Performance Characteristics

The tool uses several optimizations for high performance:

- **Multi-threading**: Parallel processing using multiple CPU cores
- **Memory-mapped I/O**: Uses `mmap()` for file input to minimize syscall overhead
- **SIMD acceleration**: Uses NEON/SSE2/AVX2 instructions for fast line/tab scanning
- **Batch processing**: Pre-computes all sample positions per line for O(1) access
- **Zero-copy parsing**: Parses VCF fields without creating intermediate strings
- **Buffered output**: Uses 16MB per-thread output buffers with direct `write()` syscalls
- **Optimized integer formatting**: Fast integer-to-string conversion

### Benchmark Results (4GB VCF, chr21, 427K variants)

| Samples | Output Lines | Time |
|---------|--------------|------|
| 1 | 427K | 1.8s |
| 10 | 4.27M | 1.7s |
| 100 | 42.7M | 3.3s |
| 500 | 213.5M | 13s |

### Output Volume

This tool generates one output line per variant × sample combination. For large population-scale VCFs, this can result in massive output:

| Input Size | Samples | Variants | Output Lines |
|------------|---------|----------|--------------|
| 50 MB | 10 | 100K | 1 million |
| 1 GB | 100 | 500K | 50 million |
| 4 GB | 2504 | 427K | 1.07 billion |

For very large datasets, use one of these output reduction strategies:
- **`--aggregate`**: Reduce output by ~2500x (O(V×S) → O(V)) - one line per variant
- **`--limit-samples N`**: Process only first N samples
- **`--binary`**: ~25x smaller output (2 bytes vs ~50 bytes per sample)
- **`--gzip`**: ~10x compression on any output format
- **`--samples`**: Select specific samples of interest

### Throughput

On typical hardware with multi-threading:
- **VCF parsing**: ~1.5 GB/s with mmap
- **Output generation**: ~10-15M lines/sec (with 100+ samples)
- **Scaling**: Near-linear scaling up to ~4-6 cores

## Limitations
- Does not distinguish between different alternate alleles (e.g., "1" vs "2")
- No options for filtering by allele count thresholds
- Cannot account for genotype quality or read depth
- Limited to processing standard VCF genotype fields
- No direct integration with population genetics metrics
