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

## References

Brunner, G. O. (1977). "A Definition of Coordination and Its Relevance
in the Structure Types AlB2 and NiAs." *Acta Crystallographica Section
A*, 33(1), 226–227.
[doi:10.1107/S0567739477000461](https://doi.org/10.1107/S0567739477000461)

## See also

Other bonding algorithms:
[`crystal_nn()`](https://prabhulab.github.io/ml-crystals/reference/crystal_nn.md),
[`econ_nn()`](https://prabhulab.github.io/ml-crystals/reference/econ_nn.md),
[`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md),
[`voronoi_nn()`](https://prabhulab.github.io/ml-crystals/reference/voronoi_nn.md)

## Examples

``` r
dists <- data.table::data.table(Atom1 = c("A", "A"), Atom2 = c("B", "C"),
                                Distance = c(1.5, 2.5),
                                DeltaX = c(1, 0), DeltaY = c(0, 1), DeltaZ = c(0, 0))
brunner_nn_reciprocal(dists)
#>     Atom1  Atom2 Distance DeltaX DeltaY DeltaZ Weight
#>    <char> <char>    <num>  <num>  <num>  <num>  <num>
#> 1:      A      B      1.5      1      0      0      1
```
