# Analyze a Batch of CIF Files

A high-level wrapper that analyzes CIF files in batch, supporting
parallel processing and batch-to-disk operations.

## Usage

``` r
analyze_cif_files(
  file_paths,
  workers = 1,
  output_dir = NULL,
  batch_size = 1000,
  ...
)
```

## Arguments

- file_paths:

  Character vector of paths or list of data.tables.

- workers:

  Integer. Number of parallel workers.

- output_dir:

  Path to output directory (optional).

- batch_size:

  Integer.

- ...:

  Args passed to `analyze_single_cif`.

## Value

Data.table or output directory path.

## Examples

``` r
# \donttest{
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  out_dir <- file.path(tempdir(), "cif_analysis")
  res <- analyze_cif_files(cif_path, output_dir = out_dir, workers = 1)

  # Clean up
  unlink(out_dir, recursive = TRUE)
}
#> Starting analysis of 1 files in 1 batches.
#> 
#> --- Processing Batch 1 of 1 (Files 1 to 1) ---
#> Batch 1 complete. Saved 1 results to '/tmp/RtmpvmV5iv/cif_analysis/batch_1.rds'.
#> ----------------------------------
#> Analysis Complete!
#> Batch results have been saved in the '/tmp/RtmpvmV5iv/cif_analysis/' directory.
# }
```
