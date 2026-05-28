# Export Analysis Results to a Directory of CSVs

Exports analysis results to CSV structure.

## Usage

``` r
export_analysis_to_csv(analysis_results, output_dir, overwrite = FALSE)
```

## Arguments

- analysis_results:

  A `data.table` object.

- output_dir:

  Path to main output directory.

- overwrite:

  Logical.

## Value

Invisibly returns output_dir.

## Examples

``` r
# \donttest{
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  res <- analyze_single_cif(cif_path)

  out_dir <- file.path(tempdir(), "cif_csvs")
  export_analysis_to_csv(res, output_dir = out_dir)

  unlink(out_dir, recursive = TRUE)
}
#> Analysis successfully exported to: /tmp/Rtmp9bus3q/cif_csvs
# }
```
