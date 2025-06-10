# VCFX_validator

## Overview
`VCFX_validator` is a utility tool for checking the validity of VCF files according to the basic VCF format specifications. It performs various checks on the file structure, header format, and data lines to ensure the file is properly formatted and contains valid data.

## Usage

```bash
VCFX_validator [OPTIONS] < input.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-s`, `--strict` | Enable stricter validation checks |
| `-d`, `--report-dups` | Report duplicate records to stderr |

## Description
`VCFX_validator` processes a VCF file to verify its structural validity by:

1. Reading the VCF file from standard input (plain or gzip/BGZF compressed)
2. Checking that all meta-information lines (starting with '##') are properly formatted
3. Validating that the #CHROM header line is present and has at least 8 required columns
4. For each data line:
   - Ensuring it has at least 8 columns
   - Verifying that CHROM is not empty
   - Confirming POS is a positive integer
   - Checking that REF and ALT contain only valid bases
   - Validating that QUAL is either '.' or a non-negative float
   - Ensuring FILTER is not empty
   - Checking INFO and FORMAT fields against header definitions
   - Validating genotype syntax
   - Detecting duplicate records when `--report-dups` is used
5. Reporting errors for any validation failures
6. Returning exit code 0 if the file is valid, or 1 if it contains errors

This tool is useful for validating VCF files before processing them with other tools, ensuring they meet the basic requirements of the VCF format specification.

## Validation Details

### Meta-Information Lines
- Must start with '##'
- No specific content validation beyond the prefix

### #CHROM Header Line
- Must be present in the file
- Must start with '#CHROM'
- Must have at least 8 columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
- Must appear before any data lines

### Data Lines
- Must have at least 8 columns
- CHROM: Must not be empty
- POS: Must be a positive integer
- ID: Can be empty or '.' (not validated)
- REF: Must contain only A,C,G,T,N
- ALT: Must contain only A,C,G,T,N
- QUAL: Must be '.' or a non-negative float
- FILTER: Must not be empty
- INFO: Keys must be defined in the header. Numeric counts are validated when numeric. Flags are allowed.

### Strict Mode
When `--strict` is used, additional checks are applied:
- The number of columns in every data line must exactly match the `#CHROM` header.
- If FORMAT/sample columns are present, each sample field must contain the same
  number of sub-fields as specified in the FORMAT column.
- Any warning that would normally be emitted is treated as an error and causes
  the validator to exit with a non-zero status.

## Examples

### Basic Validation
Check if a VCF file is valid:
```bash
VCFX_validator < input.vcf > validated.vcf
```

### Using Strict Mode
Enable stricter validation with additional checks:
```bash
VCFX_validator --strict < input.vcf > validated.vcf
```

When the input is valid, the original VCF is written unchanged to standard output,
allowing `VCFX_validator` to be used as a filter in processing pipelines. Informational
messages such as `VCF file is valid.` are printed to standard error.

### Redirecting Error Messages
Save validation errors to a file:
```bash
VCFX_validator < input.vcf 2> validation_errors.txt
```

## Example Output

### For Valid Files
```
VCF file is valid.
```

### For Invalid Files
```
Error: line 15 has <8 columns.
```
```
Error: line 42 POS must be >0.
```
```
Error: no #CHROM line found in file.
```

## Special Case Handling

### Empty Lines
- Empty lines are ignored during validation

### Malformed Header Lines
- Lines starting with '#' that are neither '##' meta-information lines nor the '#CHROM' header line are considered errors

### Missing #CHROM Line
- If no #CHROM line is found in the file, an error is reported
- Data lines encountered before a #CHROM line will cause validation to fail

### Whitespace
- Leading and trailing whitespace is trimmed from field values before validation

## Performance Considerations
- The tool processes the VCF file line by line, with minimal memory requirements
- Performance scales linearly with the size of the input file
- No external dependencies or reference files are required

## Limitations
- Does not validate VCF version compatibility
- No validation of the content of meta-information lines beyond the '##' prefix

