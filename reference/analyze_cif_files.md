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
