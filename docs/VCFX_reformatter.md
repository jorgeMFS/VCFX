# VCFX_reformatter

## Overview
`VCFX_reformatter` is a utility tool for modifying the structure of VCF files by selectively removing or reordering fields in the INFO and FORMAT columns. This allows for more compact VCF files with standardized field ordering, which can improve readability and processing efficiency.

## Usage
```bash
VCFX_reformatter [OPTIONS] < input.vcf > output.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit |
| `-c`, `--compress-info <keys>` | Remove specified INFO keys (comma-separated) |
| `-f`, `--compress-format <keys>` | Remove specified FORMAT keys (comma-separated) |
| `-i`, `--reorder-info <keys>` | Reorder these INFO keys to appear first (comma-separated) |
| `-o`, `--reorder-format <keys>` | Reorder these FORMAT keys to appear first (comma-separated) |

## Description
`VCFX_reformatter` processes a VCF file to modify its structure while preserving its content by:

1. Reading the VCF file from standard input
2. Preserving all header lines without modification
3. For each data line, depending on the options specified:
   - Removing specified INFO keys from the semicolon-delimited INFO field
   - Removing specified FORMAT keys from the colon-delimited FORMAT field and corresponding sample fields
   - Reordering INFO keys to place specified keys first, followed by remaining keys
   - Reordering FORMAT keys and correspondingly reordering each sample's subfields
4. Writing the reformatted VCF to standard output

This tool is particularly useful for:
- Removing unnecessary INFO or FORMAT fields to reduce file size
- Standardizing the order of fields across multiple VCF files
- Preparing VCF files for easier parsing or visual inspection
- Ensuring the most important fields appear first in each record

## Field Handling

### INFO Field Manipulation
The INFO field (column 8) consists of semicolon-separated key-value pairs:
- **Compression**: Specified keys are completely removed from the INFO field
- **Reordering**: Specified keys are placed at the beginning of the INFO field in the order provided, with remaining keys appended in their original order

### FORMAT Field Manipulation
The FORMAT field (column 9) and sample genotype fields (columns 10+) consist of colon-separated values:
- **Compression**: Specified keys are removed from the FORMAT field and the corresponding values are removed from each sample's field
- **Reordering**: Specified keys are placed at the beginning of the FORMAT field in the order provided, with sample fields reordered accordingly

## Examples

### Basic Compression
Remove AF and DP from INFO field:
```bash
VCFX_reformatter --compress-info AF,DP < input.vcf > output.vcf
```

### Format Field Compression
Remove GQ and DP from FORMAT field and sample genotypes:
```bash
VCFX_reformatter --compress-format GQ,DP < input.vcf > output.vcf
```

### Reordering Fields
Place AF and DP at the beginning of INFO field:
```bash
VCFX_reformatter --reorder-info AF,DP < input.vcf > output.vcf
```

### Combined Operations
Remove some fields and reorder others:
```bash
VCFX_reformatter --compress-info AC,AN --reorder-info AF,DP --compress-format PL --reorder-format GT,AD < input.vcf > output.vcf
```

## Example Transformations

### INFO Compression
```
Before: AF=0.5;DP=30;AC=1;AN=2
After:  AC=1;AN=2
```
(When using `--compress-info AF,DP`)

### FORMAT Compression
```
Before: FORMAT=GT:AD:DP:GQ:PL  SAMPLE=0/1:15,15:30:99:200,0,200
After:  FORMAT=GT:AD:PL        SAMPLE=0/1:15,15:200,0,200
```
(When using `--compress-format DP,GQ`)

### INFO Reordering
```
Before: AN=2;AC=1;DP=30;AF=0.5
After:  AF=0.5;DP=30;AN=2;AC=1
```
(When using `--reorder-info AF,DP`)

### FORMAT Reordering
```
Before: FORMAT=AD:DP:GT:GQ:PL  SAMPLE=15,15:30:0/1:99:200,0,200
After:  FORMAT=GT:GQ:AD:DP:PL  SAMPLE=0/1:99:15,15:30:200,0,200
```
(When using `--reorder-format GT,GQ`)

## Handling Special Cases

### Empty Fields
- If all INFO or FORMAT fields are removed, a `.` is output in their place
- Empty input fields are preserved as empty in the output

### Malformed Lines
- Lines with fewer than 8 columns are skipped with a warning
- Data lines encountered before the header line (`#CHROM`) are skipped with a warning

### Missing Keys
- When reordering, if specified keys are not found in the original data, they are simply ignored
- When compressing, if specified keys are not found, no action is taken for those keys

## Performance Considerations
- The tool processes VCF files line by line, with minimal memory requirements
- Performance scales linearly with the size of the input file
- Large VCF files with many samples may experience slower processing when FORMAT field manipulation is used

## Limitations
- No validation of VCF format or field content
- Cannot modify the fixed columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER)
- Cannot add new fields or modify field values
- No handling of compressed (gzipped) VCF files directly
- No complex field transformations or conditional modifications 