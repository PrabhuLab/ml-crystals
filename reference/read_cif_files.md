# Read CIF Files into Memory

Reads one or more CIF files from disk and loads each into a
`data.table`. It includes encoding sanitization to handle legacy
characters.

## Usage

``` r
read_cif_files(file_paths)
```

## Arguments

- file_paths:

  A character vector of paths to the CIF files.

## Value

A list of `data.table` objects.
