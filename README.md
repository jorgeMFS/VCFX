# VCFX

VCFX is a collection of C/C++ tools for processing and analyzing VCF (Variant Call Format) files, with optional WebAssembly compatibility. Each tool is an independent command-line executable that can parse input from `stdin` and write to `stdout`, enabling flexible piping and integration into bioinformatics pipelines.

## Tools

Below is the list of tools provided by VCFX:

1. **VCFX_header_parser**
2. **VCFX_record_filter**
3. **VCFX_field_extractor**
4. **VCFX_format_converter**
5. **VCFX_variant_counter**
6. **VCFX_sample_extractor**
7. **VCFX_sorter**
8. **VCFX_validator**
9. **VCFX_subsampler**
10. **VCFX_genotype_query**
11. **VCFX_allele_freq_calc**
12. **VCFX_indexer**
13. **VCFX_compressor**
14. **VCFX_position_subsetter**
15. **VCFX_haplotype_extractor**
16. **VCFX_info_parser**
17. **VCFX_variant_classifier**
18. **VCFX_duplicate_remover**
19. **VCFX_info_summarizer**
20. **VCFX_distance_calculator**
21. **VCFX_multiallelic_splitter**
22. **VCFX_missing_data_handler**
23. **VCFX_concordance_checker**
24. **VCFX_allele_balance_calc**
25. **VCFX_allele_counter**
26. **VCFX_phase_checker**
27. **VCFX_annotation_extractor**
28. **VCFX_phred_filter**
29. **VCFX_merger**
30. **VCFX_metadata_summarizer**
31. **VCFX_hwe_tester**
32. **VCFX_fasta_converter**
33. **VCFX_nonref_filter**
34. **VCFX_dosage_calculator**
35. **VCFX_population_filter**
36. **VCFX_file_splitter**
37. **VCFX_gl_filter**
38. **VCFX_ref_comparator**
39. **VCFX_ancestry_inferrer**
40. **VCFX_impact_filter**
41. **VCFX_info_aggregator**
42. **VCFX_probability_filter**
43. **VCFX_diff_tool**
44. **VCFX_cross_sample_concordance**
45. **VCFX_phase_quality_filter**
46. **VCFX_indel_normalizer**
47. **VCFX_custom_annotator**
48. **VCFX_region_subsampler**
49. **VCFX_allele_balance_filter**
50. **VCFX_missing_detector**
51. **VCFX_haplotype_phaser**
52. **VCFX_af_subsetter**
53. **VCFX_sv_handler**
54. **VCFX_reformatter**
55. **VCFX_quality_adjuster**
56. **VCFX_inbreeding_calculator**
57. **VCFX_outlier_detector**
58. **VCFX_alignment_checker**
59. **VCFX_ancestry_assigner**
60. **VCFX_ld_calculator**

Each tool has its own minimal `CMakeLists.txt` in `src/VCFX_*`.

---

## Building (Native)

1. **Clone** or obtain this repository:

   ```bash
   git clone https://github.com/youruser/VCFX.git
   cd VCFX
   ```

2. **Create** a build directory and run CMake:

   ```bash
   mkdir build && cd build
   cmake ..  # default native build
   make      # or cmake --build .
   ```

3. **All** tools** will be built as individual executables in the `build/src/VCFX_*` directories (or similar, depending on your generator).

---

## Building (WebAssembly)

If you have [Emscripten](https://emscripten.org/) installed, you can build all tools to **WASM**. For instance:

1. **Install** or activate your Emscripten environment so the Emscripten compiler is on your path.
2. **Create** a separate build directory:

   ```bash
   mkdir build_wasm && cd build_wasm
   ```

3. **Configure** with the `-DBUILD_WASM=ON` option:

   ```bash
   cmake -DBUILD_WASM=ON ..
   cmake --build .
   ```

4. Each tool is now compiled to WebAssembly. Depending on your Emscripten configuration, you might get `.wasm` + `.js`, or `.html` output files. You can run them in a browser or Node environment as appropriate.

---

## Running the Tools

Each tool is a command-line program that reads from `stdin` and writes to `stdout`. For example:

```bash
cd build
./src/VCFX_variant_counter/VCFX_variant_counter < ../test_data/sample.vcf
```

*(Paths vary by your local structure. In general: `./src/<ToolName>/<ToolName> ...`.)*

**Similarly** for WASM** output, you might have `VCFX_variant_counter.html` or `.js` in your build directory, which you can run in a browser or Node environment.

---

## Running Tests

We use [CTest](https://cmake.org/cmake/help/latest/manual/ctest.1.html). After building:

```bash
cd build
ctest --verbose
```

This runs a suite of basic unit or integration tests from the `tests/` directory, verifying minimal input-output correctness or exit codes for each tool.

---

## License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.