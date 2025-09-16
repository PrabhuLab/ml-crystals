## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

options(width = 80)

options(datatable.print.width = 80)

## ----setup--------------------------------------------------------------------
library(crystract)

## ----full-pipeline-with-comment-----------------------------------------------
# Find the path to the single example CIF file included in the package
cif_path <- system.file("extdata", "ICSD422.cif", package = "crystract")

# --- Example for batch processing a directory ---
# To load all CIF files from the package's example directory, you would 
# uncomment and run this:
# extdata_dir <- system.file("extdata", package = "crystract")
# all_cif_paths <- list.files(path = extdata_dir, pattern = "\\.cif$", full.names = TRUE)
# analysis_results_batch <- analyze_cif_files(all_cif_paths)

# For this vignette, we will run the pipeline on just our single example file
analysis_results <- analyze_cif_files(cif_path)

# Let's inspect the structure of the output table.
# It's a single row containing all our results in nested data.tables.
str(analysis_results, max.level = 2)

## ----access-nested-data-------------------------------------------------------
# The result is a list-column, so we access the element with [[1]]
final_bonds <- analysis_results$bonded_pairs[[1]]

print(head(final_bonds))

## ----load-data----------------------------------------------------------------
# The path was defined in the previous section:
cif_data_list <- read_cif_files(cif_path)

# We'll work with the content of the first file
cif_content <- cif_data_list[[1]]

# Let's look at the first few lines of the raw data
knitr::kable(
  head(cif_content),
  caption = "First 6 lines of the raw CIF data."
)

## ----extract-metadata---------------------------------------------------------
database_code <- extract_database_code(cif_content)
chemical_formula <- extract_chemical_formula(cif_content)
space_group_name <- extract_space_group_name(cif_content)
space_group_number <- extract_space_group_number(cif_content)

cat("Database Code:", database_code, "\n")
cat("Chemical Formula:", chemical_formula, "\n")
cat("Space Group:", space_group_name, "(No.", space_group_number, ")\n")

## ----extract-metrics----------------------------------------------------------
unit_cell_metrics <- extract_unit_cell_metrics(cif_content)
print(unit_cell_metrics)

## ----extract-coords-symm------------------------------------------------------
# Extract the coordinates of the unique atoms in the asymmetric unit
atomic_coordinates <- extract_atomic_coordinates(cif_content)
print("Asymmetric Atomic Coordinates:")
print(atomic_coordinates)

# Extract the symmetry operations
symmetry_operations <- extract_symmetry_operations(cif_content)
print("Symmetry Operations (first 6 of 8):")
print(head(symmetry_operations))

## ----generate-structure-------------------------------------------------------
# Apply symmetry to generate all atoms in the primary unit cell
transformed_coords <- apply_symmetry_operations(atomic_coordinates, symmetry_operations)
print("Unique atoms in full unit cell (first 6 of 12):")
print(head(transformed_coords))

# Expand into a 3x3x3 supercell for neighbor calculations
expanded_coords <- expand_transformed_coords(transformed_coords)
print("Atoms in supercell (first 6 of 324):")
print(head(expanded_coords))

## ----calc-distances-----------------------------------------------------------
distances <- calculate_distances(atomic_coordinates, expanded_coords, unit_cell_metrics)
print("Calculated Distances (shortest 6):")
print(head(distances[order(Distance)]))

## ----calc-bonding-neighbors---------------------------------------------------
# Identify bonded pairs using a tolerance of 10%
bonded_pairs <- minimum_distance(distances, delta = 0.1)
print("Bonded Pairs (first 6):")
print(head(bonded_pairs))

# Calculate neighbor counts based on the bonded pairs
neighbor_counts <- calculate_neighbor_counts(bonded_pairs)
print("Neighbor Counts:")
print(neighbor_counts)

## ----calc-angles--------------------------------------------------------------
bond_angles <- calculate_angles(
  bonded_pairs,
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics
)
print("Calculated Bond Angles (first 6):")
print(head(bond_angles))

## ----propagate-errors---------------------------------------------------------
# Propagate errors for interatomic distances
bonded_pairs_with_error <- propagate_distance_error(
  bonded_pairs,
  atomic_coordinates,
  unit_cell_metrics
)
print("Bonded Pairs with Distance Error (first 6):")
print(head(bonded_pairs_with_error))

# Propagate errors for bond angles
bond_angles_with_error <- propagate_angle_error(
  bond_angles,
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics
)
print("Bond Angles with Angle Error (first 6):")
print(head(bond_angles_with_error))

## ----eval=FALSE---------------------------------------------------------------
# # In an interactive R session, you would run this:
# filtered_bonds <- filter_atoms_by_symbol(
#   data_table = bonded_pairs_with_error,
#   atom_col = "Atom1" # Filter based on the central atom
# )

## ----filter-simulation--------------------------------------------------------
# Manually simulate the user typing "Si" at the prompt
user_input <- "Si"
symbols_to_keep <- trimws(strsplit(user_input, ",")[[1]])
patterns <- sapply(symbols_to_keep, function(sym) paste0("(^", sym, "$)|(^", sym, "[^A-Za-z])"))
full_pattern <- paste(patterns, collapse = "|")
filtered_bonds_by_symbol <- bonded_pairs_with_error[grepl(full_pattern, get("Atom1"))]

print("Original table contains bonds centered on both Sr and Si atoms:")
print(bonded_pairs_with_error$Atom1)
cat("\n")
print("Filtered table now only contains bonds centered on 'Si1':")
print(filtered_bonds_by_symbol$Atom1)

## ----filter-by-wyckoff--------------------------------------------------------
# 1. In our example, all asymmetric atoms occupy the Wyckoff site 'c' with multiplicity 4 ("4c").
print("Atomic coordinates showing Wyckoff information:")
print(atomic_coordinates[, .(Label, WyckoffSymbol, WyckoffMultiplicity)])
cat("\n")

# 2. Filter bonds where the central atom is on the "4c" Wyckoff site.
bonds_from_4c_site <- filter_by_wyckoff_symbol(
  data_table = bonded_pairs_with_error,
  atomic_coordinates = atomic_coordinates,
  atom_col = "Atom1",
  wyckoff_symbols = "4c"
)

print(paste("Number of rows in original bond table:", nrow(bonded_pairs_with_error)))
print(paste("Number of rows after filtering for site '4c':", nrow(bonds_from_4c_site)))

## ----filter-ghost-distances---------------------------------------------------
# We use the full, unfiltered distance table for this demonstration.
# The function returns a list with two tables: 'kept' and 'removed'.
# A distance is kept if: Distance > (Radius1 + Radius2) * tolerance
filtered_result <- filter_ghost_distances(
    distances = distances,
    atomic_coordinates = atomic_coordinates,
    tolerance = 0.4 # Default tolerance
)

kept_distances <- filtered_result$kept
removed_distances <- filtered_result$removed

cat("Total distances calculated:", nrow(distances), "\n")
cat("Distances kept after filtering:", nrow(kept_distances), "\n")
cat("Ghost distances removed:", nrow(removed_distances), "\n\n")

# For a well-ordered structure like Sr2Si, no ghost bonds are expected.
# The 'removed' table should be empty, confirming the filter's behavior.
print("Removed ghost distances (should be empty for this well-ordered example):")
print(removed_distances)

## ----filter-by-elements-------------------------------------------------------
# Let's filter our bond table to exclude any bonds involving Strontium ("Sr").
# Since all bonds in this structure are Si-Sr, the result should be an empty table.
bonds_without_sr <- filter_by_elements(
    distances = bonded_pairs_with_error,
    atomic_coordinates = atomic_coordinates,
    elements_to_exclude = "Sr"
)

cat("Number of bonds in original table:", nrow(bonded_pairs_with_error), "\n")
cat("Number of bonds after excluding 'Sr':", nrow(bonds_without_sr), "\n")

## ----calculate-weighted-average-distance--------------------------------------
# Calculate the weighted average distance for the entire Sr2Si network.
# First, we define the Wyckoff sites belonging to the network. Here, it's just "4c".
network_wyckoff_sites <- "4c"

# The function uses the full distance table and coordinate info.
# Note: For this ordered structure, occ_j = 1.0, so the result is the same
# as a multiplicity-only weighted average. The difference would appear in a
# disordered or vacancy-containing structure.
weighted_avg_dist <- calculate_weighted_average_network_distance(
    distances = distances,
    atomic_coordinates = atomic_coordinates,
    wyckoff_symbols = network_wyckoff_sites
)

cat("Weighted average network distance for the '4c' sites:", weighted_avg_dist, "Å\n")

## ----export-results, message=TRUE---------------------------------------------
# We will use the 'analysis_results' table from the previous sections.

# 1. Define a path for our export. We use a temporary directory for this example.
export_path <- file.path(tempdir(), "crysmal_export_demo")

# 2. Run the export function. 
# We use overwrite = TRUE to ensure the example runs cleanly every time.
export_analysis_to_csv(
  analysis_results = analysis_results,
  output_dir = export_path,
  overwrite = TRUE
)

# 3. We can list the files to see the created directory structure.
cat("\nGenerated file structure:\n")
list.files(export_path, recursive = TRUE)

# Clean up the temporary directory after the demonstration
unlink(export_path, recursive = TRUE)

