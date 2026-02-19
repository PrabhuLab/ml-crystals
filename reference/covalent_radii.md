# Atomic Radii Data for Bond-Length Estimation

A data.table containing atomic radii for elements. This data is used for
estimating plausible bond lengths to filter out non-physical "ghost"
distances that can occur in disordered structures. Radii are in
Angstroms (Å).

## Usage

``` r
covalent_radii
```

## Format

A data table with three columns:

- Symbol:

  The chemical symbol of the element.

- Radius:

  The atomic radius in Angstroms.

- Type:

  The type of radius, e.g., "covalent". The default table only contains
  covalent radii.

## Source

J. Emsley. The Elements. Third edition 1998, Oxford University Press. As
provided by Julia-Maria Huebner.
