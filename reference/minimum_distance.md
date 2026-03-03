# Identify Atomic Bonds using the Minimum Distance Method

Identifies bonded atoms by finding the nearest neighbor distance (d_min)
for each central atom and defining a cutoff distance (d_cut) as d_cut =
(1 + delta) \* d_min.

## Usage

``` r
minimum_distance(distances, delta = 0.1)
```

## Arguments

- distances:

  A `data.table` of interatomic distances from `calculate_distances`.

- delta:

  The relative tolerance parameter (default 0.1).

## Value

A `data.table` of bonded pairs.

## See also

Other bonding algorithms:
[`brunner_nn_reciprocal()`](https://prabhulab.github.io/ml-crystals/reference/brunner_nn_reciprocal.md),
[`crystal_nn()`](https://prabhulab.github.io/ml-crystals/reference/crystal_nn.md),
[`econ_nn()`](https://prabhulab.github.io/ml-crystals/reference/econ_nn.md),
[`voronoi_nn()`](https://prabhulab.github.io/ml-crystals/reference/voronoi_nn.md)

## Examples

``` r
dists <- data.table::data.table(Atom1 = c("A", "A"), Atom2 = c("B", "C"),
                                Distance = c(1.5, 2.5),
                                DeltaX = c(1, 0), DeltaY = c(0, 1), DeltaZ = c(0, 0))
minimum_distance(dists, delta = 0.1)
#>     Atom1  Atom2 Distance DeltaX DeltaY DeltaZ Weight
#>    <char> <char>    <num>  <num>  <num>  <num>  <num>
#> 1:      A      B      1.5      1      0      0      1
```
