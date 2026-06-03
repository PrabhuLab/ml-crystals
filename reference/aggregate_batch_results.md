# Aggregate Batched Analysis Results

Reads all `batch_*.rds` files from a directory and combines them.

## Usage

``` r
aggregate_batch_results(input_dir, cols_to_keep = NULL)
```

## Arguments

- input_dir:

  The directory containing RDS files from `analyze_cif_files`.

- cols_to_keep:

  Optional character vector of column names to keep.

## Value

A single `data.table` containing the aggregated results.

## Examples

``` r
# \donttest{
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  old_threads <- data.table::setDTthreads(2)
  out_dir <- file.path(tempdir(), "cif_batches")
  analyze_cif_files(cif_path, output_dir = out_dir)

  aggr <- aggregate_batch_results(out_dir)

  unlink(out_dir, recursive = TRUE)
  data.table::setDTthreads(old_threads)
}
#> Starting analysis of 1 files in 1 batches.
#> 
#> --- Processing Batch 1 of 1 (Files 1 to 1) ---
#> Batch 1 complete. Saved 1 results to '/tmp/Rtmp5yW9R0/cif_batches/batch_1.rds'.
#> ----------------------------------
#> Analysis Complete!
#> Batch results have been saved in the '/tmp/Rtmp5yW9R0/cif_batches/' directory.
#> Found 1 batch files to aggregate.
#> Aggregation complete. Total rows: 1.
# }
```
