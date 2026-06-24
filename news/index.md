# Changelog

## crystract 1.0.2

### Performance Improvements

- **Massive Speedup in Symmetry Operations:** Symmetry equations are now
  compiled into vectorized R functions once per structure, eliminating
  the severe CPU bottleneck caused by repeated string parsing
  (`eval(parse())`).
- **“Fast Fail” CIF Validation:** Implemented high-speed pre-flight
  checks using raw C-level string matching. Invalid or severely
  malformed CIF files are now rejected in milliseconds, skipping
  unnecessary and expensive full-table regex scans.
- **Optimized Extraction:** Replaced multiple heavy regular expression
  searches with highly optimized exact string matching during the
  mandatory data extraction phase.

### Bug Fixes

- **Batch Processing Memory Exhaustion Resolved:** Fixed a critical
  issue in
  [`analyze_cif_files()`](https://prabhulab.github.io/ml-crystals/reference/analyze_cif_files.md)
  where large datasets would crash R by attempting to load all
  structures into system RAM simultaneously. CIF files are now
  dynamically read from disk strictly within their specified chunk/batch
  size, ensuring stable memory footprints when processing tens of
  thousands of files.

## crystract 1.0.1

CRAN release: 2026-05-28

- Refactored internal CIF extraction logic to robustly check for
  multiple alternative CIF tags when parsing chemical formulas,
  structure types, and space group names/numbers.
- Updated `README` installation instructions to include options for both
  stable (CRAN) and development (GitHub) releases.
- Improved the package vignette by updating Wyckoff filtering examples
  (from `4a` to `4d`) and adding a new section demonstrating how
  experimental errors are extracted from CIF loops and propagated
  through calculations.

## crystract 1.0.0

CRAN release: 2026-03-20

- Initial CRAN submission.
