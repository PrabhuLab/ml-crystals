## ----include = FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# Save current options and apply new ones for the vignette rendering
old_options <- options(width = 75, datatable.print.width = 80)

## ----setup, message=FALSE, warning=FALSE---------------------------------
library(crystract)
library(data.table)

## ----full-pipeline-with-comments, eval=TRUE------------------------------
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

## ----access-nested-data, eval=FALSE--------------------------------------
# # The result is a list-column, so we access the element with [[1]]
# final_bonds <- analysis_results$bonds_minimum_distance[[1]]
# 
# print(head(final_bonds))

## ----load-data, eval=TRUE------------------------------------------------
# The path was defined in the previous section:
cif_data_list <- read_cif_files(cif_path)

# We'll work with the content of the first file
cif_content <- cif_data_list[[1]]

# Let's look at the first few lines of the raw data
knitr::kable(
  head(cif_content),
  caption = "First 6 lines of the raw CIF data."
)

## ----extract-metadata, eval=TRUE-----------------------------------------
database_code <- extract_database_code(cif_content)
chemical_formula <- extract_chemical_formula(cif_content)
space_group_name <- extract_space_group_name(cif_content)
space_group_number <- extract_space_group_number(cif_content)

cat("Database Code:", database_code, "\n")
cat("Chemical Formula:", chemical_formula, "\n")
cat("Space Group:", space_group_name, "(No.", space_group_number, ")\n")

## ----extract-metrics, eval=TRUE------------------------------------------
unit_cell_metrics <- extract_unit_cell_metrics(cif_content)
print(unit_cell_metrics)

## ----extract-coords-symm, eval=TRUE--------------------------------------
# Extract the coordinates of the unique atoms in the asymmetric unit
atomic_coordinates <- extract_atomic_coordinates(cif_content)
print("Asymmetric Atomic Coordinates:")
print(atomic_coordinates)

# Extract the symmetry operations
symmetry_operations <- extract_symmetry_operations(cif_content)
print("Symmetry Operations (first 6 of 8):")
print(head(symmetry_operations))

## ----generate-structure, eval=TRUE---------------------------------------
# Apply symmetry to generate all atoms in the primary unit cell
transformed_coords <- apply_symmetry_operations(atomic_coordinates, symmetry_operations, unit_cell_metrics)
print("Unique atoms in full unit cell (first 6):")
print(head(transformed_coords))

# Expand into a supercell for neighbor calculations
expanded_coords <- expand_transformed_coords(transformed_coords)
print("Atoms in supercell (first 6):")
print(head(expanded_coords))

## ----calc-distances, eval=TRUE-------------------------------------------
distances <- calculate_distances(atomic_coordinates, expanded_coords, unit_cell_metrics)
print("Calculated Distances (shortest 6):")
print(head(distances[order(Distance)]))

## ----calc-bonding-neighbors, eval=TRUE-----------------------------------
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
cat("CrystalNN found", nrow(bonds_crystal_nn), "bonds.\n")

# Calculate integer neighbor counts based on the bonded pairs (e.g., using CrystalNN)
neighbor_counts <- calculate_neighbor_counts(bonds_crystal_nn)
print("CrystalNN Neighbor Counts:")
print(neighbor_counts)

## ----calc-angles, eval=TRUE----------------------------------------------
bond_angles <- calculate_angles(
  bonds_min,
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics
)
print("Calculated Bond Angles (first 6):")
print(head(bond_angles))

## ----propagate-errors, eval=TRUE-----------------------------------------
# Propagate errors for interatomic distances
bonded_pairs_with_error <- propagate_distance_error(
  bonds_min,
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

## ----eval=FALSE----------------------------------------------------------
# # In an interactive R session, you would run this:
# filtered_bonds <- filter_atoms_by_symbol(
#   data_table = bonded_pairs_with_error,
#   atom_col = "Atom1" # Filter based on the central atom
# )

## ----filter-by-wyckoff, eval=TRUE----------------------------------------
# 1. Let's look at the Wyckoff information already extracted from the CIF.
# All atoms in this structure are on the "a" site with a multiplicity of 4.
print("Atomic coordinates showing Wyckoff information:")
print(atomic_coordinates[, .(Label, WyckoffSymbol, WyckoffMultiplicity)])
cat("\n")

# 2. Filter bonds where the central atom is on the "4d" Wyckoff site.
bonds_from_4d_site <- filter_by_wyckoff_symbol(
  data_table = bonded_pairs_with_error,
  atomic_coordinates = atomic_coordinates,
  atom_col = "Atom1",
  wyckoff_symbols = "4d"
)

print(paste("Number of rows in original bond table:", nrow(bonded_pairs_with_error)))
print(paste("Number of rows after filtering for site '4d':", nrow(bonds_from_4d_site)))

## ----filter-ghost-distances, eval=TRUE-----------------------------------
# A distance d is kept if: (r1+r2)*(1-margin) <= d <= (r1+r2)*(1+margin)
filtered_result <- filter_ghost_distances(
    distances = distances,
    atomic_coordinates = atomic_coordinates,
    margin = 0.1 # Default margin is 10%
)

kept_distances <- filtered_result$kept
removed_distances <- filtered_result$removed

cat("Total distances calculated:", nrow(distances), "\n")
cat("Distances kept after filtering:", nrow(kept_distances), "\n")
cat("Unlikely / non-bonded distances removed:", nrow(removed_distances), "\n\n")

print("Subset of removed distances (showing physically impossible / too long connections):")
print(head(removed_distances))

## ----filter-by-elements, eval=TRUE---------------------------------------
# Let's filter our bond table to exclude any bonds involving Tin ("Sn").
# The resulting table will only contain Cu-Se bonds.
bonds_without_sn <- filter_by_elements(
    distances = bonded_pairs_with_error,
    atomic_coordinates = atomic_coordinates,
    elements_to_exclude = "Ge"
)

cat("Number of bonds in original table:", nrow(bonded_pairs_with_error), "\n")
cat("Number of bonds after excluding 'Ge':", nrow(bonds_without_sn), "\n")

## ----calculate-weighted-average-distance, eval=TRUE----------------------
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

## ----custom-radii-and-batch, eval=FALSE----------------------------------
# # --- 1. Define and set the custom radii dictionary ---
# radii_dict <- c(
#   Ac=2.15, Ag=1.45, Al=1.21, Am=1.8, Ar=1.06, As=1.19, At=1.5, Au=1.36, B=0.84,
#   Ba=2.15, Be=0.96, Bi=1.48, Br=1.2, C=0.73, Ca=1.76, Cd=1.44, Ce=2.04, Cl=1.02,
#   Cm=1.69, Co=1.38, Cr=1.39, Cs=2.44, Cu=1.32, Dy=1.92, Er=1.89, Eu=1.98, F=0.57,
#   Fe=1.42, Fr=2.6, Ga=1.22, Gd=1.96, Ge=1.2, H=0.31, He=0.28, Hf=1.75, Hg=1.32,
#   Ho=1.92, I=1.39, In=1.42, Ir=1.41, K=2.03, Kr=1.16, La=2.07, Li=1.28, Lu=1.87,
#   Mg=1.41, Mn=1.5, Mo=1.54, N=0.71, Na=1.66, Nb=1.64, Nd=2.01, Ne=0.58, Ni=1.24,
#   Np=1.9, O=0.66, Os=1.44, P=1.07, Pa=2, Pb=1.46, Pd=1.39, Pm=1.99, Po=1.4, Pr=2.03,
#   Pt=1.36, Pu=1.87, Ra=2.21, Rb=2.2, Re=1.51, Rh=1.42, Rn=1.5, Ru=1.46, S=1.05,
#   Sb=1.39, Sc=1.7, Se=1.2, Si=1.11, Sm=1.98, Sn=1.39, Sr=1.95, Ta=1.7, Tb=1.94,
#   Tc=1.47, Te=1.38, Th=2.06, Ti=1.6, Tl=1.45, Tm=1.9, U=1.96, V=1.53, W=1.62,
#   Xe=1.4, Y=1.9, Yb=1.87, Zn=1.22, Zr=1.75
# )
# 
# custom_radii_table <- data.table(
#   Symbol = names(radii_dict),
#   Radius = as.numeric(radii_dict),
#   Type   = "covalent"
# )
# 
# # Inject the custom table into the current crystract session
# set_radii_data(custom_radii_table)
# 
# # --- 2. Run the batch analysis ---
# cat("Starting batch analysis on CIF files...\n")
# results <- analyze_cif_files(
#   file_paths = "path/to/cif_directory",
#   tolerance = 1e-4,
#   bonding_algorithms = c("minimum_distance", "brunner", "econ", "voronoi", "crystal_nn"),
#   calculate_bond_angles = FALSE,       # Skip angles to speed up extraction
#   perform_error_propagation = FALSE,   # Skip uncertainties
#   workers = 4                          # Adjust based on your available CPU cores
# )
# 
# # --- 3. Extract and Save Coordination Numbers to CSV ---
# # Create a dedicated folder for the outputs to avoid clutter
# output_dir <- "path/to/output_directory"
# if (!dir.exists(output_dir)) dir.create(output_dir)
# 
# # Identify all coordination number columns generated by the pipeline
# cn_columns <- grep("^cn_", names(results), value = TRUE)
# 
# cat("\nExtracting coordination numbers and saving to CSV...\n")
# 
# for (col in cn_columns) {
# 
#   # Extract the list of data.tables for this specific algorithm
#   list_of_dts <- results[[col]]
# 
#   # Associate each table with its CIF file name
#   names(list_of_dts) <- results$file_name
# 
#   # Unnest and bind all tables together into one giant table
#   combined_dt <- rbindlist(list_of_dts, idcol = "file_name", fill = TRUE)
# 
#   if (nrow(combined_dt) > 0) {
#     # Format the file path (e.g., "cn_crystal_nn_crystract.csv")
#     file_path <- file.path(output_dir, paste0(col, "_crystract.csv"))
# 
#     # Save to CSV
#     fwrite(combined_dt, file_path)
# 
#     cat(sprintf("Saved %d rows across %d files to: %s\n",
#                 nrow(combined_dt),
#                 uniqueN(combined_dt$file_name),
#                 basename(file_path)))
#   } else {
#     cat(sprintf("No valid coordination data found for %s, skipping.\n", col))
#   }
# }
# 
# cat("\nAll operations finished successfully! Files are in:", output_dir, "\n")

## ----include = FALSE----------------------------------------------------------
# Restore original options as required by CRAN
options(old_options)

