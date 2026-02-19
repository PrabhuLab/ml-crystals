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
