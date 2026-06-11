# Using the crystract Package for Crystallographic Analysis

## Introduction

This vignette demonstrates the core workflow for using the `crystract`
package. The package is designed for the **batch processing** of
Crystallographic Information Files (.cif) to extract structural and
chemical information from the repeating atomic arrangement that defines
a crystal. It calculates key geometric properties like bond lengths and
angles and compiles the results into a structured format suitable for
large-scale data analysis.

We will first demonstrate the main batch-processing function,
`analyze_cif_files`. Then, to illustrate how it works, we will walk
through the analysis of a single file step-by-step, explaining the
underlying functions, the crystallographic principles, and the
mathematical formulas they employ.

### Setup: Loading the Package

First, we load the `crystract` package.

``` r

library(crystract)
library(data.table)

# Limit data.table threads to 2 to comply with CRAN policies on CPU time
old_dt_threads <- data.table::setDTthreads(2)
```

## 1. The Core `crystract` Workflow

The primary goal of `crystract` is to automate the analysis of many CIF
files at once. While the package provides granular functions for each
step of the crystallographic analysis, the most common use case is to
leverage the main wrapper function for an end-to-end pipeline.

### 1.1 The Full Pipeline for Batch Processing

The
[`analyze_cif_files()`](https://prabhulab.github.io/ml-crystals/reference/analyze_cif_files.md)
function (and its single-file counterpart
[`analyze_single_cif()`](https://prabhulab.github.io/ml-crystals/reference/analyze_single_cif.md))
is the cornerstone of the package. It performs the entire sequence of
operations—from reading files to calculating bond angles with error
propagation. You can easily specify which bonding algorithms to run
simultaneously.

``` r

# IMPORTANT: Update this path to point to your own downloaded CIF file.
# 1. Try to find the file in the installed package
cif_path <- system.file("extdata", "2405219.cif", package = "crystract")

# 2. If that fails, look in the source directory
if (cif_path == "") {
  cif_path <- "../inst/extdata/2405219.cif"
}

# 3. Final check to provide a clear error if both fail
if (!file.exists(cif_path)) {
  stop("Could not find 2405219.cif in installed package or inst/extdata/")
}

# Run the pipeline on our single example file, extracting multiple bonding types at once
analysis_results <- analyze_single_cif(
  cif_path,
  bonding_algorithms = c("minimum_distance", "brunner", "econ", "voronoi", "crystal_nn")
)

# Let's inspect the structure of the output table.
# It's a single row containing all our results in nested data.tables.
str(analysis_results, max.level = 2)
#> Classes 'data.table' and 'data.frame':   1 obs. of  27 variables:
#>  $ file_name              : chr "2405219.cif"
#>  $ database_code          : chr "depnum_ccdc_archive CCDC 2405219"
#>  $ chemical_formula       : chr "Ge Yb"
#>  $ structure_type         : logi NA
#>  $ space_group_name       : chr "P n m a"
#>  $ space_group_number     : chr "62"
#>  $ unit_cell_metrics      :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    1 obs. of  12 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ atomic_coordinates     :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    2 obs. of  14 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ symmetry_operations    :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    8 obs. of  3 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ transformed_coords     :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    8 obs. of  5 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ expanded_coords        :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    216 obs. of  8 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ distances              :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    430 obs. of  6 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ bonds_minimum_distance :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    9 obs. of  8 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>   .. ..- attr(*, "sorted")= chr [1:3] "Atom1" "Atom2" "Distance"
#>  $ cn_minimum_distance    :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    2 obs. of  3 variables:
#>   .. ..- attr(*, "sorted")= chr "Atom"
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ angles_minimum_distance:List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    22 obs. of  5 variables:
#>   .. ..- attr(*, "sorted")= chr [1:3] "CentralAtom" "Neighbor1" "Neighbor2"
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ bonds_brunner          :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    9 obs. of  8 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>   .. ..- attr(*, "sorted")= chr [1:3] "Atom1" "Atom2" "Distance"
#>  $ cn_brunner             :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    2 obs. of  3 variables:
#>   .. ..- attr(*, "sorted")= chr "Atom"
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ angles_brunner         :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    22 obs. of  5 variables:
#>   .. ..- attr(*, "sorted")= chr [1:3] "CentralAtom" "Neighbor1" "Neighbor2"
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ bonds_econ             :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    18 obs. of  8 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>   .. ..- attr(*, "sorted")= chr [1:3] "Atom1" "Atom2" "Distance"
#>  $ cn_econ                :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    2 obs. of  3 variables:
#>   .. ..- attr(*, "sorted")= chr "Atom"
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ angles_econ            :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    72 obs. of  5 variables:
#>   .. ..- attr(*, "sorted")= chr [1:3] "CentralAtom" "Neighbor1" "Neighbor2"
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ bonds_voronoi          :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    30 obs. of  11 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>   .. ..- attr(*, "sorted")= chr [1:3] "Atom1" "Atom2" "Distance"
#>  $ cn_voronoi             :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    2 obs. of  3 variables:
#>   .. ..- attr(*, "sorted")= chr "Atom"
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ angles_voronoi         :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    0 obs. of  4 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>   .. ..- attr(*, "sorted")= chr "CentralAtom"
#>  $ bonds_crystal_nn       :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    26 obs. of  8 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>   .. ..- attr(*, "sorted")= chr [1:3] "Atom1" "Atom2" "Distance"
#>  $ cn_crystal_nn          :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    2 obs. of  3 variables:
#>   .. ..- attr(*, "sorted")= chr "Atom"
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>  $ angles_crystal_nn      :List of 1
#>   ..$ :Classes 'data.table' and 'data.frame':    0 obs. of  4 variables:
#>   .. ..- attr(*, ".internal.selfref")=<pointer: 0x560af40b0110> 
#>   .. ..- attr(*, "sorted")= chr "CentralAtom"
#>  - attr(*, ".internal.selfref")=<pointer: 0x560af40b0110>
```

As shown in the structure output, the result is a tidy `data.table`
containing all extracted and calculated information. We can easily
access the nested data frames for further analysis. Note that the output
includes separate results for each bonding algorithm requested, such as
`bonds_minimum_distance`, `bonds_brunner`, `bonds_voronoi`, etc.

To get the final bonded pairs table from the Minimum Distance method
(which includes propagated errors):

``` r

# The result is a list-column, so we access the element with [[1]]
final_bonds <- analysis_results$bonds_minimum_distance[[1]]

print(head(final_bonds))
```

### 1.2 A Step-by-Step Walkthrough

To understand what happens under the hood, this section breaks down the
process using individual functions. We will use the same CIF file to
demonstrate each step.

#### 1.2.1 Loading CIF Data

We use the package’s
[`read_cif_files()`](https://prabhulab.github.io/ml-crystals/reference/read_cif_files.md)
function to load data into memory.

``` r

# The path was defined in the previous section:
cif_data_list <- read_cif_files(cif_path)

# We'll work with the content of the first file
cif_content <- cif_data_list[[1]]

# Let's look at the first few lines of the raw data
knitr::kable(
  head(cif_content),
  caption = "First 6 lines of the raw CIF data."
)
```

| V1                                                                       |
|:-------------------------------------------------------------------------|
| \####################################################################### |
| \#                                                                       |
| \# This file contains crystal structure data downloaded from the         |
| \# Cambridge Structural Database (CSD) hosted by the Cambridge           |
| \# Crystallographic Data Centre (CCDC).                                  |
| \#                                                                       |

First 6 lines of the raw CIF data. {.table}

#### 1.2.2 Extracting Metadata and Unit Cell Parameters

A crystal structure is defined by its **unit cell** and the arrangement
of atoms within it. We first extract this high-level information.

``` r

database_code <- extract_database_code(cif_content)
chemical_formula <- extract_chemical_formula(cif_content)
space_group_name <- extract_space_group_name(cif_content)
space_group_number <- extract_space_group_number(cif_content)

cat("Database Code:", database_code, "\n")
#> Database Code: depnum_ccdc_archive CCDC 2405219
cat("Chemical Formula:", chemical_formula, "\n")
#> Chemical Formula: Ge Yb
cat("Space Group:", space_group_name, "(No.", space_group_number, ")\n")
#> Space Group: P n m a (No. 62 )
```

Next,
[`extract_unit_cell_metrics()`](https://prabhulab.github.io/ml-crystals/reference/extract_unit_cell_metrics.md)
extracts the six parameters that define the shape and size of the unit
cell: the lengths of its three axes ($`a`$, $`b`$, $`c`$) and the angles
between them ($`\alpha`$, $`\beta`$, $`\gamma`$). Their experimental
uncertainties are also extracted.

``` r

unit_cell_metrics <- extract_unit_cell_metrics(cif_content)
print(unit_cell_metrics)
#>    _cell_length_a _cell_length_b _cell_length_c _cell_angle_alpha
#>             <num>          <num>          <num>             <num>
#> 1:         7.9006         3.8981         5.8729                90
#>    _cell_angle_beta _cell_angle_gamma _cell_length_a_error
#>               <num>             <num>                <num>
#> 1:               90                90               0.0019
#>    _cell_length_b_error _cell_length_c_error _cell_angle_alpha_error
#>                   <num>                <num>                   <num>
#> 1:                9e-04               0.0015                      NA
#>    _cell_angle_beta_error _cell_angle_gamma_error
#>                     <num>                   <num>
#> 1:                     NA                      NA
```

#### 1.2.3 Extracting Atomic and Symmetry Data

Instead of listing every atom in the unit cell, a CIF file efficiently
describes the structure using only the **asymmetric unit**: the smallest
set of unique atoms. All other atoms in the unit cell can be generated
by applying the crystal’s **symmetry operations** to this unique set.

``` r

# Extract the coordinates of the unique atoms in the asymmetric unit
atomic_coordinates <- extract_atomic_coordinates(cif_content)
print("Asymmetric Atomic Coordinates:")
#> [1] "Asymmetric Atomic Coordinates:"
print(atomic_coordinates)
#>     Label TypeSymbol WyckoffSymbol WyckoffMultiplicity OxidationState
#>    <char>     <char>        <char>               <num>          <num>
#> 1:    Yb1         Yb             d                   4             NA
#> 2:    Ge1         Ge             d                   4             NA
#>    ThermalParam Occupancy OccupancyError     x_a   y_b    z_c x_error
#>           <num>     <num>          <num>   <num> <num>  <num>   <num>
#> 1:       0.0215         1              0 0.17675  0.25 0.6121 0.00019
#> 2:       0.0238         1              0 0.03550  0.25 0.1268 0.00040
#>    y_error z_error
#>      <num>   <num>
#> 1:       0   2e-04
#> 2:       0   5e-04

# Extract the symmetry operations
symmetry_operations <- extract_symmetry_operations(cif_content)
print("Symmetry Operations (first 6 of 8):")
#> [1] "Symmetry Operations (first 6 of 8):"
print(head(symmetry_operations))
#>         x      y      z
#>    <char> <char> <char>
#> 1:      x      y      z
#> 2: -x+1/2     -y  z+1/2
#> 3:     -x  y+1/2     -z
#> 4:  x+1/2 -y+1/2 -z+1/2
#> 5:     -x     -y     -z
#> 6:  x+1/2      y -z+1/2
```

#### 1.2.4 Generating the Full Crystal Structure

Now we use the asymmetric atoms and symmetry operations to
computationally build the full crystal structure. This is a two-step
process.

**Formula: Symmetry Operations**

A symmetry operation is an affine transformation that maps an initial
fractional coordinate $`(x, y, z)`$ to a new coordinate
$`(x', y', z')`$. It consists of a rotation/reflection matrix
$`\mathbf{W}`$ and a translation vector $`\mathbf{w}`$.

``` math
\begin{pmatrix} x' \\ y' \\ z' \end{pmatrix} =
\begin{pmatrix} W_{11} & W_{12} & W_{13} \\ W_{21} & W_{22} & W_{23} \\ W_{31} & W_{32} & W_{33} \end{pmatrix}
\begin{pmatrix} x \\ y \\ z \end{pmatrix} +
\begin{pmatrix} w_1 \\ w_2 \\ w_3 \end{pmatrix}
\quad \text{or} \quad \mathbf{x'} = \mathbf{W}\mathbf{x} + \mathbf{w}
```

The `apply_symmetry_operations` function applies all of the crystal’s
symmetry operations to each atom in the asymmetric unit, generating the
complete set of atoms within the primary unit cell. Note that it
requires the unit cell metrics to resolve distance tolerances across
boundaries.

**Formula: Supercell Expansion**

To find all nearest neighbors of an atom, we must also consider atoms in
adjacent unit cells. The `expand_transformed_coords` function generates
a **supercell** (in this case, 3x3x3) by translating the primary unit
cell atoms by integer vectors $`(i, j, k)`$, where $`i, j, k`$ each
range from -1 to 1. An atom at fractional coordinate $`(x, y, z)`$
generates new coordinates:

``` math
(x_{exp}, y_{exp}, z_{exp}) = (x + i, y + j, z + k)
```

``` r

# Apply symmetry to generate all atoms in the primary unit cell
transformed_coords <- apply_symmetry_operations(atomic_coordinates, symmetry_operations, unit_cell_metrics)
print("Unique atoms in full unit cell (first 6):")
#> [1] "Unique atoms in full unit cell (first 6):"
print(head(transformed_coords))
#>     Label SymOp     x_a   y_b    z_c
#>    <char> <int>   <num> <num>  <num>
#> 1:  Yb1_1     1 0.17675  0.25 0.6121
#> 2:  Yb1_2     2 0.32325  0.75 0.1121
#> 3:  Yb1_3     3 0.82325  0.75 0.3879
#> 4:  Yb1_4     4 0.67675  0.25 0.8879
#> 5:  Ge1_1     1 0.03550  0.25 0.1268
#> 6:  Ge1_2     2 0.46450  0.75 0.6268

# Expand into a supercell for neighbor calculations
expanded_coords <- expand_transformed_coords(transformed_coords)
print("Atoms in supercell (first 6):")
#> [1] "Atoms in supercell (first 6):"
print(head(expanded_coords))
#>             Label SymOp      x_a   y_b     z_c    Tx    Ty    Tz
#>            <char> <int>    <num> <num>   <num> <int> <int> <int>
#> 1: Yb1_1_-1_-1_-1     1 -0.82325 -0.75 -0.3879    -1    -1    -1
#> 2: Yb1_2_-1_-1_-1     2 -0.67675 -0.25 -0.8879    -1    -1    -1
#> 3: Yb1_3_-1_-1_-1     3 -0.17675 -0.25 -0.6121    -1    -1    -1
#> 4: Yb1_4_-1_-1_-1     4 -0.32325 -0.75 -0.1121    -1    -1    -1
#> 5: Ge1_1_-1_-1_-1     1 -0.96450 -0.75 -0.8732    -1    -1    -1
#> 6: Ge1_2_-1_-1_-1     2 -0.53550 -0.25 -0.3732    -1    -1    -1
```

#### 1.2.5 Calculating Interatomic Distances

**Formula: Interatomic Distance in a Triclinic System**

Because unit cell axes are not always orthogonal, the simple Pythagorean
theorem is insufficient. The distance $`d`$ between two atoms at
fractional coordinates $`(x_{f1}, y_{f1}, z_{f1})`$ and
$`(x_{f2}, y_{f2}, z_{f2})`$ requires the full crystallographic distance
formula:

``` math
d = \left(
\begin{aligned}
& a^2(x_{f1}-x_{f2})^2 + b^2(y_{f1}-y_{f2})^2 + c^2(z_{f1}-z_{f2})^2 \\
& + 2ab(x_{f1}-x_{f2})(y_{f1}-y_{f2})\cos\gamma \\
& + 2ac(x_{f1}-x_{f2})(z_{f1}-z_{f2})\cos\beta \\
& + 2bc(y_{f1}-y_{f2})(z_{f1}-z_{f2})\cos\alpha
\end{aligned}
\right)^{1/2}
```

The `calculate_distances` function computes the distances from each
central atom (from the asymmetric unit) to all other atoms in the
generated supercell.

``` r

distances <- calculate_distances(atomic_coordinates, expanded_coords, unit_cell_metrics)
print("Calculated Distances (shortest 6):")
#> [1] "Calculated Distances (shortest 6):"
print(head(distances[order(Distance)]))
#>     Atom1          Atom2 Distance  DeltaX DeltaY  DeltaZ
#>    <char>         <char>    <num>   <num>  <num>   <num>
#> 1:    Ge1 Ge1_3_-1_-1_-1 2.516281 0.07100    0.5  0.2536
#> 2:    Ge1  Ge1_3_-1_0_-1 2.516281 0.07100   -0.5  0.2536
#> 3:    Ge1  Yb1_3_-1_-1_0 2.993686 0.21225    0.5 -0.2611
#> 4:    Yb1  Ge1_3_-1_-1_0 2.993686 0.21225    0.5 -0.2611
#> 5:    Ge1   Yb1_3_-1_0_0 2.993686 0.21225   -0.5 -0.2611
#> 6:    Yb1   Ge1_3_-1_0_0 2.993686 0.21225   -0.5 -0.2611
```

#### 1.2.6 Identifying Bonds and Neighbors

**Formula: Bonding and Coordination Number**

Identifying which atoms are “bonded” is key to determining the
**coordination number (CN)**. `crystract` implements several robust,
geometric, and mathematical approaches:

- **`minimum_distance`**: Defines a custom distance cutoff for each
  central atom $`i`$:
  ``` math
   
  d_i^\text{cut} = (1+\delta) d_i^\text{min} 
  ```
  Here, $`d_i^\text{min}`$ is the shortest distance from atom $`i`$ to
  any other atom, and $`\delta`$ is a tolerance parameter (default is
  0.1).
- **`brunner`**: Identifies bonds by finding the largest gap in the
  reciprocal distances. Matches
  `pymatgen.analysis.local_env.BrunnerNNReciprocal`.
- **`econ`**: Uses Hoppe’s iterative method to calculate Effective
  Coordination Numbers, weighting neighbor contributions exponentially
  by distance.
- **`voronoi`**: Performs a 3D Voronoi tessellation to find neighbors
  sharing a face. *Note: `crystract` calculates these faces using the
  Delaunay dual. Results may slightly differ from exact Voronoi
  tessellation implementations (like Voro++) due to geometric handling
  of degenerate, co-planar vertices.*
- **`crystal_nn`**: A sophisticated approach combining Voronoi solid
  angles, ionic radii distance cutoffs, and electronegativity weights.

Let’s run each algorithm to demonstrate:

``` r

# 1. Minimum Distance
bonds_min <- minimum_distance(distances, delta = 0.1)

# 2. Brunner's Method
bonds_brunner <- brunner_nn_reciprocal(distances)

# 3. Hoppe's ECoN
bonds_econ <- econ_nn(distances, atomic_coordinates)

# 4. Voronoi Tessellation
bonds_voronoi <- voronoi_nn(atomic_coordinates, expanded_coords, unit_cell_metrics)

# 5. CrystalNN
bonds_crystal_nn <- crystal_nn(distances, atomic_coordinates, expanded_coords, unit_cell_metrics)

cat("Minimum Distance found", nrow(bonds_min), "bonds.\n")
#> Minimum Distance found 9 bonds.
cat("CrystalNN found", nrow(bonds_crystal_nn), "bonds.\n")
#> CrystalNN found 26 bonds.

# Calculate integer neighbor counts based on the bonded pairs (e.g., using CrystalNN)
neighbor_counts <- calculate_neighbor_counts(bonds_crystal_nn)
print("CrystalNN Neighbor Counts:")
#> [1] "CrystalNN Neighbor Counts:"
print(neighbor_counts)
#>      Atom CoordinationNumber
#>    <char>              <int>
#> 1:  Yb1_1                 17
#> 2:  Ge1_1                  9
```

#### 1.2.7 Calculating Bond Angles

**Formula: Bond Angle Calculation**

To calculate a bond angle A-B-C (with B at the vertex), the fractional
coordinates must first be converted to an orthogonal Cartesian system.

1.  **Fractional to Cartesian Conversion:** A fractional coordinate
    $`(\mathbf{x}_f)`$ is converted to a Cartesian coordinate
    $`(\mathbf{x}_c)`$ via a transformation matrix $`\mathbf{M}`$:

``` math
    \begin{pmatrix} x_c \\ y_c \\ z_c \end{pmatrix} = \mathbf{M} \cdot \begin{pmatrix} x_f \\ y_f \\ z_f \end{pmatrix} \quad \text{where} \quad
    \mathbf{M} = \begin{pmatrix}
    a & b \cos\gamma & c \cos\beta \\
    0 & b \sin\gamma & \frac{c(\cos\alpha - \cos\beta \cos\gamma)}{\sin\gamma} \\
    0 & 0 & \frac{V}{ab\sin\gamma}
    \end{pmatrix}
```
where $`V`$ is the unit cell volume.

2.  **Angle via Dot Product:** With atoms in Cartesian space, the angle
    $`\theta`$ between vectors $`\vec{u}`$ (from B to A) and $`\vec{v}`$
    (from B to C) is:

``` math
\theta = \arccos\left( \frac{\vec{u} \cdot \vec{v}}{|\vec{u}| |\vec{v}|} \right)
```

The `calculate_angles` function implements this for all possible bond
angles around each central atom. We’ll pass it our `bonds_min` results:

``` r

bond_angles <- calculate_angles(
  bonds_min,
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics
)
print("Calculated Bond Angles (first 6):")
#> [1] "Calculated Bond Angles (first 6):"
print(head(bond_angles))
#>    CentralAtom      Neighbor1     Neighbor2    Angle
#>         <char>         <char>        <char>    <num>
#> 1:         Ge1 Ge1_3_-1_-1_-1 Ge1_3_-1_0_-1 101.5332
#> 2:         Yb1    Ge1_1_0_0_0   Ge1_1_0_0_1 138.3540
#> 3:         Yb1    Ge1_1_0_0_0  Ge1_2_0_-1_0 107.6689
#> 4:         Yb1    Ge1_1_0_0_0   Ge1_2_0_0_0 107.6689
#> 5:         Yb1    Ge1_1_0_0_0 Ge1_3_-1_-1_0 105.8268
#> 6:         Yb1    Ge1_1_0_0_0  Ge1_3_-1_0_0 105.8268
```

#### 1.2.8 Error Propagation

Finally, `crystract` propagates the experimental uncertainties (i.e.,
the estimated standard deviations of the unit cell parameters and
fractional atomic coordinates) from cell parameters and atomic
coordinates to the calculated distances and angles.

Errors in .cif files can be found in the respective tables for unit cell
parameters and atomic coordinates in parenthesis following the
positional values. Below, the errors for the atomic coordinates of the
test file are shown:

``` md
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_adp_type
_atom_site_U_iso_or_equiv
_atom_site_site_symmetry_multiplicity
_atom_site_occupancy
_atom_site_calc_flag
_atom_site_refinement_flags
_atom_site_disorder_assembly
_atom_site_disorder_group
Yb1 Yb 0.17675(19) 0.25 0.6121(2) Uani 0.0215(4) 4 1 d . . .
Ge1 Ge 0.0355(4) 0.25 0.1268(5) Uani 0.0238(10) 4 1 d . . .
```

**Formula: Error Propagation**

The uncertainty, $`\sigma_f`$, in a calculated value $`f`$ that depends
on several uncorrelated variables $`p_i`$ (each with uncertainty
$`\sigma_{p_i}`$) is found using the sum of squares of the partial
derivatives:

``` math
\sigma_f^2 = \sum_i \left( \frac{\partial f}{\partial p_i} \sigma_{p_i} \right)^2
```

##### Uncertainty in Interatomic Distance ($`\sigma_d`$)

The uncertainty in the distance, $`\sigma_d`$, is found by applying the
general formula to the distance $`d`$. The input parameters $`p_i`$ are
the 12 variables that define the distance:
$`\{a, b, c, \alpha, \beta, \gamma, x_{f1}, y_{f1}, z_{f1}, x_{f2}, y_{f2}, z_{f2}\}`$.
Letting $`\Delta x = x_{f1} - x_{f2}`$, etc., the partial derivatives
are, for example:

``` math
\frac{\partial d}{\partial a} = \frac{1}{2d} \left( 2a(\Delta x)^2 + 2b(\Delta x)(\Delta y)\cos\gamma + 2c(\Delta x)(\Delta z)\cos\beta \right)
```

##### Uncertainty in Bond Angle ($`\sigma_\theta`$)

Propagating error to the bond angle $`\theta`$ is a more complex,
multi-step process involving propagating initial errors first to
Cartesian coordinates, and then from those Cartesian uncertainties to
the final angle via the chain rule.

``` r

# Propagate errors for interatomic distances
bonded_pairs_with_error <- propagate_distance_error(
  bonds_min,
  atomic_coordinates,
  unit_cell_metrics
)
print("Bonded Pairs with Distance Error (first 6):")
#> [1] "Bonded Pairs with Distance Error (first 6):"
print(head(bonded_pairs_with_error))
#> Key: <Atom1, Atom2, Distance>
#>     Atom1          Atom2 Distance   DeltaX DeltaY  DeltaZ    Weight
#>    <char>         <char>    <num>    <num>  <num>   <num>     <num>
#> 1:    Ge1 Ge1_3_-1_-1_-1 2.516281  0.07100    0.5  0.2536 1.0000000
#> 2:    Ge1  Ge1_3_-1_0_-1 2.516281  0.07100   -0.5  0.2536 1.0000000
#> 3:    Yb1    Ge1_1_0_0_0 3.060807  0.14125    0.0  0.4853 0.9780708
#> 4:    Yb1    Ge1_1_0_0_1 3.222200  0.14125    0.0 -0.5147 0.9290813
#> 5:    Yb1   Ge1_2_0_-1_0 2.995761 -0.28775    0.5 -0.0147 0.9993073
#> 6:    Yb1    Ge1_2_0_0_0 2.995761 -0.28775   -0.5 -0.0147 0.9993073
#>    DistanceError
#>            <num>
#> 1:   0.002684669
#> 2:   0.002684669
#> 3:   0.003281605
#> 4:   0.003286949
#> 5:   0.002704675
#> 6:   0.002704675

# Propagate errors for bond angles
bond_angles_with_error <- propagate_angle_error(
  bond_angles,
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics
)
print("Bond Angles with Angle Error (first 6):")
#> [1] "Bond Angles with Angle Error (first 6):"
print(head(bond_angles_with_error))
#> Key: <CentralAtom, Neighbor1, Neighbor2>
#>    CentralAtom      Neighbor1     Neighbor2    Angle AngleError
#>         <char>         <char>        <char>    <num>      <num>
#> 1:         Ge1 Ge1_3_-1_-1_-1 Ge1_3_-1_0_-1 101.5332 0.12841756
#> 2:         Yb1    Ge1_1_0_0_0   Ge1_1_0_0_1 138.3540 0.09644671
#> 3:         Yb1    Ge1_1_0_0_0  Ge1_2_0_-1_0 107.6689 0.07978724
#> 4:         Yb1    Ge1_1_0_0_0   Ge1_2_0_0_0 107.6689 0.07978724
#> 5:         Yb1    Ge1_1_0_0_0 Ge1_3_-1_-1_0 105.8268 0.08131349
#> 6:         Yb1    Ge1_1_0_0_0  Ge1_3_-1_0_0 105.8268 0.08131349
```

## 2. Tools for Post-Processing and Analysis

After running the main analysis pipeline, `crystract` provides several
helper functions to filter and refine the results, allowing you to
isolate data based on chemical identity or key crystallographic
properties.

### 2.1 Filtering by Chemical Identity

The
[`filter_atoms_by_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_atoms_by_symbol.md)
function provides an interactive way to filter results (like bond or
angle tables) to focus on the coordination environment around a specific
type of element.

``` r

# In an interactive R session, you would run this:
filtered_bonds <- filter_atoms_by_symbol(
  data_table = bonded_pairs_with_error,
  atom_col = "Atom1" # Filter based on the central atom
)
```

If you were to type `Si` and press Enter, the function would return a
new `data.table` containing only the rows where the central atom
(`Atom1`) is Silicon.

### 2.2 Filtering by Wyckoff Position

For many crystallographic studies, it is crucial to analyze atoms based
on their specific site symmetry, described by their **Wyckoff position**
(e.g., “2a”, “6c”). The
[`filter_by_wyckoff_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_wyckoff_symbol.md)
function is designed for this purpose.

``` r

# 1. Let's look at the Wyckoff information already extracted from the CIF.
# All atoms in this structure are on the "a" site with a multiplicity of 4.
print("Atomic coordinates showing Wyckoff information:")
#> [1] "Atomic coordinates showing Wyckoff information:"
print(atomic_coordinates[, .(Label, WyckoffSymbol, WyckoffMultiplicity)])
#>     Label WyckoffSymbol WyckoffMultiplicity
#>    <char>        <char>               <num>
#> 1:    Yb1             d                   4
#> 2:    Ge1             d                   4
cat("\n")

# 2. Filter bonds where the central atom is on the "4d" Wyckoff site.
bonds_from_4d_site <- filter_by_wyckoff_symbol(
  data_table = bonded_pairs_with_error,
  atomic_coordinates = atomic_coordinates,
  atom_col = "Atom1",
  wyckoff_symbols = "4d"
)

print(paste("Number of rows in original bond table:", nrow(bonded_pairs_with_error)))
#> [1] "Number of rows in original bond table: 9"
print(paste("Number of rows after filtering for site '4d':", nrow(bonds_from_4d_site)))
#> [1] "Number of rows after filtering for site '4d': 9"
```

### 2.3 Filtering Ghost Distances Using Atomic Radii

A common issue in crystallography is **site-occupancy disorder**, where
a single crystallographic site is statistically co-occupied by different
atoms. Standard calculations can produce physically meaningless, short
“ghost” distances.

The
[`filter_ghost_distances()`](https://prabhulab.github.io/ml-crystals/reference/filter_ghost_distances.md)
function cleans the distance table by removing these artifacts. It uses
a built-in table of covalent radii to establish a plausible bond length
range. Any calculated distance falling outside this physical range
(defined by a `margin`) is filtered out.

*Note: Since this function evaluates the entire supercell distance
matrix, pairs of atoms that are simply far away will be logged as “TOO
LONG”.*

``` r

# A distance d is kept if: (r1+r2)*(1-margin) <= d <= (r1+r2)*(1+margin)
filtered_result <- filter_ghost_distances(
    distances = distances,
    atomic_coordinates = atomic_coordinates,
    margin = 0.1 # Default margin is 10%
)

kept_distances <- filtered_result$kept
removed_distances <- filtered_result$removed

cat("Total distances calculated:", nrow(distances), "\n")
#> Total distances calculated: 430
cat("Distances kept after filtering:", nrow(kept_distances), "\n")
#> Distances kept after filtering: 20
cat("Unlikely / non-bonded distances removed:", nrow(removed_distances), "\n\n")
#> Unlikely / non-bonded distances removed: 410

print("Subset of removed distances (showing physically impossible / too long connections):")
#> [1] "Subset of removed distances (showing physically impossible / too long connections):"
print(head(removed_distances))
#>     Atom1          Atom2  Distance expected_dist lower_bound upper_bound
#>    <char>         <char>     <num>         <num>       <num>       <num>
#> 1:    Ge1 Ge1_1_-1_-1_-1 10.587994          2.44       2.196       2.684
#> 2:    Ge1 Ge1_2_-1_-1_-1  5.724757          2.44       2.196       2.684
#> 3:    Ge1 Ge1_4_-1_-1_-1  7.098444          2.44       2.196       2.684
#> 4:    Ge1  Ge1_1_0_-1_-1  7.048839          2.44       2.196       2.684
#> 5:    Ge1  Ge1_2_0_-1_-1  4.889711          2.44       2.196       2.684
#> 6:    Ge1  Ge1_3_0_-1_-1  7.738707          2.44       2.196       2.684
#>                  Reason
#>                  <char>
#> 1: Distance is TOO LONG
#> 2: Distance is TOO LONG
#> 3: Distance is TOO LONG
#> 4: Distance is TOO LONG
#> 5: Distance is TOO LONG
#> 6: Distance is TOO LONG
```

### 2.4 Filtering by Element Exclusion

To focus an analysis on a specific part of a structure (e.g., a host
framework without guest atoms),
[`filter_by_elements()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_elements.md)
allows for the complete removal of bonds involving certain elements.

``` r

# Let's filter our bond table to exclude any bonds involving Tin ("Sn").
# The resulting table will only contain Cu-Se bonds.
bonds_without_sn <- filter_by_elements(
    distances = bonded_pairs_with_error,
    atomic_coordinates = atomic_coordinates,
    elements_to_exclude = "Ge"
)

cat("Number of bonds in original table:", nrow(bonded_pairs_with_error), "\n")
#> Number of bonds in original table: 9
cat("Number of bonds after excluding 'Ge':", nrow(bonds_without_sn), "\n")
#> Number of bonds after excluding 'Ge': 0
```

### 2.5 Calculating Weighted Average Network Bond Distance

A common goal is to compute a single summary statistic that represents a
structure.

The
[`calculate_weighted_average_network_distance()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_average_network_distance.md)
function calculates a representative **bond length** for a specified
atomic network (defined by a set of Wyckoff sites).

This calculation should always be performed on a table of identified
**bonds** to ensure the result is physically meaningful.

#### Formula: Weighted Average Network Distance

To obtain a single, representative bond length for an atomic network,
the function calculates a true weighted average over all individual
bonds in the unit cell. This is equivalent to calculating the **total
sum of all bond lengths** within the network and dividing by the **total
number of bonds**. This method correctly accounts for site multiplicity,
site occupancy, and varying coordination numbers.

The formula is defined as:

\$\$ \large \bar{d}\_\text{network} = \frac{\sum_j (m_j \cdot
\text{occ}\_j \cdot S_j)}{\sum_j (m_j \cdot \text{occ}\_j \cdot n_j)}
\$\$

Where, for each unique atomic site *j* in the asymmetric unit: - $`m_j`$
is the **Wyckoff multiplicity** (how many equivalent atoms are generated
by symmetry). - $`\text{occ}_j`$ is the **site occupancy** (accounts for
disorder; 1.0 for a fully occupied site). - $`n_j`$ is the
**coordination number** (the number of bonds for an atom at site *j*). -
$`S_j`$ is the **sum of all bond distances** for a single atom at site
*j*. It is calculated as:

\$\$ \large S_j = \sum\_{i=1}^{n_j} d\_{ij} \$\$

The numerator, $`\sum_j (m_j \cdot \text{occ}_j \cdot S_j)`$, represents
the grand total sum of the lengths of every bond in the unit cell. The
denominator, $`\sum_j (m_j \cdot \text{occ}_j \cdot n_j)`$, represents
the total number of bonds in the unit cell.

#### Best Practices for Accurate Calculation

For the most physically meaningful and accurate results, it is crucial
to perform this calculation as the final step in a filtering pipeline,
especially when dealing with complex or disordered structures. The
recommended workflow is:

1.  Calculate all potential interatomic distances using
    [`calculate_distances()`](https://prabhulab.github.io/ml-crystals/reference/calculate_distances.md).
2.  Clean this raw data by removing non-physical “ghost” distances with
    [`filter_ghost_distances()`](https://prabhulab.github.io/ml-crystals/reference/filter_ghost_distances.md).
3.  Identify the true bonded pairs from the cleaned data using an
    algorithm like
    [`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md).
4.  **Only then**, pass the resulting table of confirmed bonds to
    [`calculate_weighted_average_network_distance()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_average_network_distance.md).

This “Ghosts -\> Bonds -\> Average” sequence ensures the final metric is
calculated only over actual chemical bonds, preventing artifacts from
skewing the result.

``` r

# Calculate the weighted average BOND distance for the network.
# First, identify the bonds in the structure. We use the `bonds_min` table
# created in section 1.2.6.

# Then, define the Wyckoff sites belonging to the network. 
# For this file, the framework atoms are on the "4d" sites.
network_wyckoff_sites <- "4d"

# Apply the function to the table of identified bonds.
weighted_avg_bond_dist <- calculate_weighted_average_network_distance(
    distances = bonds_min, # Use the bond table as input
    atomic_coordinates = atomic_coordinates,
    wyckoff_symbols = network_wyckoff_sites
)

cat("Weighted average network bond distance for the '4d' sites:", weighted_avg_bond_dist, "Å\n")
#> Weighted average network bond distance for the '4d' sites: 2.939673 Å
```

## 3. End-to-End Example: High-Throughput Batch Processing & CSV Extraction

While the previous examples demonstrated `crystract` on a single file to
explain the underlying mechanics, the package is heavily optimized for
processing thousands of structures simultaneously using parallel
workers.

Below is a complete, real-world example demonstrating how you might
override the internal atomic radii tables with a custom dictionary, run
an exhaustive parallel analysis using **all five bonding algorithms**,
and then unnest the resulting tables to save the extracted Coordination
Numbers (CNs) out to CSV files.

``` r

# --- 1. Define and set the custom radii dictionary ---
radii_dict <- c(
  Ac=2.15, Ag=1.45, Al=1.21, Am=1.8, Ar=1.06, As=1.19, At=1.5, Au=1.36, B=0.84, 
  Ba=2.15, Be=0.96, Bi=1.48, Br=1.2, C=0.73, Ca=1.76, Cd=1.44, Ce=2.04, Cl=1.02, 
  Cm=1.69, Co=1.38, Cr=1.39, Cs=2.44, Cu=1.32, Dy=1.92, Er=1.89, Eu=1.98, F=0.57, 
  Fe=1.42, Fr=2.6, Ga=1.22, Gd=1.96, Ge=1.2, H=0.31, He=0.28, Hf=1.75, Hg=1.32, 
  Ho=1.92, I=1.39, In=1.42, Ir=1.41, K=2.03, Kr=1.16, La=2.07, Li=1.28, Lu=1.87, 
  Mg=1.41, Mn=1.5, Mo=1.54, N=0.71, Na=1.66, Nb=1.64, Nd=2.01, Ne=0.58, Ni=1.24, 
  Np=1.9, O=0.66, Os=1.44, P=1.07, Pa=2, Pb=1.46, Pd=1.39, Pm=1.99, Po=1.4, Pr=2.03, 
  Pt=1.36, Pu=1.87, Ra=2.21, Rb=2.2, Re=1.51, Rh=1.42, Rn=1.5, Ru=1.46, S=1.05, 
  Sb=1.39, Sc=1.7, Se=1.2, Si=1.11, Sm=1.98, Sn=1.39, Sr=1.95, Ta=1.7, Tb=1.94, 
  Tc=1.47, Te=1.38, Th=2.06, Ti=1.6, Tl=1.45, Tm=1.9, U=1.96, V=1.53, W=1.62, 
  Xe=1.4, Y=1.9, Yb=1.87, Zn=1.22, Zr=1.75
)

custom_radii_table <- data.table(
  Symbol = names(radii_dict),
  Radius = as.numeric(radii_dict),
  Type   = "covalent"
)

# Inject the custom table into the current crystract session
set_radii_data(custom_radii_table)

# --- 2. Run the batch analysis ---
cat("Starting batch analysis on CIF files...\n")
results <- analyze_cif_files(
  file_paths = "path/to/cif_directory",
  tolerance = 1e-4,
  bonding_algorithms = c("minimum_distance", "brunner", "econ", "voronoi", "crystal_nn"),
  calculate_bond_angles = FALSE,       # Skip angles to speed up extraction
  perform_error_propagation = FALSE,   # Skip uncertainties
  workers = 4                          # Adjust based on your available CPU cores
)

# --- 3. Extract and Save Coordination Numbers to CSV ---
# Create a dedicated folder for the outputs to avoid clutter
output_dir <- "path/to/output_directory"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Identify all coordination number columns generated by the pipeline
cn_columns <- grep("^cn_", names(results), value = TRUE)

cat("\nExtracting coordination numbers and saving to CSV...\n")

for (col in cn_columns) {
  
  # Extract the list of data.tables for this specific algorithm
  list_of_dts <- results[[col]]
  
  # Associate each table with its CIF file name
  names(list_of_dts) <- results$file_name
  
  # Unnest and bind all tables together into one giant table
  combined_dt <- rbindlist(list_of_dts, idcol = "file_name", fill = TRUE)
  
  if (nrow(combined_dt) > 0) {
    # Format the file path (e.g., "cn_crystal_nn_crystract.csv")
    file_path <- file.path(output_dir, paste0(col, "_crystract.csv"))
    
    # Save to CSV
    fwrite(combined_dt, file_path)
    
    cat(sprintf("Saved %d rows across %d files to: %s\n", 
                nrow(combined_dt), 
                uniqueN(combined_dt$file_name), 
                basename(file_path)))
  } else {
    cat(sprintf("No valid coordination data found for %s, skipping.\n", col))
  }
}

cat("\nAll operations finished successfully! Files are in:", output_dir, "\n")
```

This workflow ensures that you can adapt `crystract`’s core metrics to
your own chemical definitions while harnessing its parallel speed to
crunch massive datasets efficiently.

## Conclusion

This vignette has demonstrated the core workflow of the `crystract`
package, from the high-level, automated processing of multiple CIF files
with `analyze_cif_files` to a detailed, step-by-step breakdown of the
underlying calculations. We have seen how the package extracts
fundamental data, builds a complete crystal model, calculates key
geometric properties across five different mathematical definitions of
bonding, and rigorously propagates experimental uncertainties into these
final results.

Furthermore, with a powerful suite of post-processing tools such as
`filter_atoms_by_symbol`, `filter_by_wyckoff_symbol`,
`filter_ghost_distances`, `calculate_weighted_average_network_distance`,
and user-defined radii tables, you can easily refine, clean, analyze,
and share the generated datasets to focus on the specific chemical or
crystallographic features of interest.

The goal of `crystract` is to provide a robust and accessible platform
for crystallographic data mining in R.
