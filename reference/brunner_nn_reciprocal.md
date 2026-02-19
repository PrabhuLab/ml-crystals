# Identify Atomic Bonds using Brunner's Method (Reciprocal)

An alternative bonding detection method using the largest reciprocal
gap. Matches `pymatgen.analysis.local_env.BrunnerNNReciprocal`.

## Usage

``` r
brunner_nn_reciprocal(distances, tol = 0)
```

## Arguments

- distances:

  A `data.table` of interatomic distances.

- tol:

  A small distance offset (default 0.0) to add to the cutoff.

## Value

A `data.table` of bonded pairs.

## See also

Other bonding algorithms:
[`crystal_nn()`](https://prabhulab.github.io/ml-crystals/reference/crystal_nn.md),
[`econ_nn()`](https://prabhulab.github.io/ml-crystals/reference/econ_nn.md),
[`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md),
[`voronoi_nn()`](https://prabhulab.github.io/ml-crystals/reference/voronoi_nn.md)
