# Changelog

## crystract 1.0.1

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
