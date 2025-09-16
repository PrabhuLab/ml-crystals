# --------------------------------------------------------------------------
# Direct CIF Metadata Extraction using crystract
# --------------------------------------------------------------------------
# This script uses specific, exported functions from the crystract package
# to extract metadata .cif files and saves
# the result to a single summary CSV file.

# --- 1. SETUP: Load Libraries and Define Paths ---

# Ensure the crystract package is installed and loaded
if (!require("crystract")) install.packages("crystract")
if (!require("data.table")) install.packages("data.table")
library(crystract)
library(data.table)

# --- USER-DEFINED VARIABLES ---
root_data_dir <- "/Users/donngo/repos/ml-crystals/packages/crystract/inst/extdata"

output_csv_path <- "/Users/donngo/repos/ml-crystals/required_ICSD_files.csv"
# -----------------------------


# --- 2. FILE DISCOVERY ---

cat("Step 1: Searching for .cif files in:", root_data_dir, "\n")

if (!dir.exists(root_data_dir)) {
  stop("Error: The specified root data directory does not exist: ", root_data_dir)
}

# Recursively find all CIF files (case-insensitive)
all_cif_paths <- list.files(
  path = root_data_dir,
  pattern = "\\.cif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(all_cif_paths) == 0) {
  stop("Execution stopped: No .cif files were found in the specified directory.")
}
cat("... Found", length(all_cif_paths), "CIF files. Starting extraction...\n")


# --- 3. PROCESSING LOOP ---

# Use lapply to process each file and store the resulting data row in a list
results_list <- lapply(all_cif_paths, function(file_path) {
  tryCatch({
    # Read the content of a single CIF file
    cif_content <- data.table::fread(file_path, sep = "\n", header = FALSE, strip.white = FALSE)

    # --- Call the exported functions from crystract ---

    # Extract single-value metadata
    single_values <- data.table(
      file_name = basename(file_path),
      category = basename(dirname(file_path)),
      database_code = extract_database_code(cif_content),
      chemical_formula = extract_chemical_formula(cif_content),
      structure_type = extract_structure_type(cif_content),
      space_group_name = extract_space_group_name(cif_content),
      space_group_number = extract_space_group_number(cif_content)
    )

  }, error = function(e) {
    # If a file fails, print a warning and return NULL so the script doesn't stop
    warning("Failed to process file '", file_path, "': ", e$message)
    return(NULL)
  })
})

# Combine all the successfully processed data rows into a single data.table
final_results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
cat("... Extraction complete.", nrow(final_results), "files were processed successfully.\n\n")


# --- 4. EXPORT TO CSV ---

cat("Step 2: Exporting metadata summary to:", output_csv_path, "\n")

# Create the output directory if it doesn't already exist
output_dir <- dirname(output_csv_path)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Write the final table to a single CSV file
data.table::fwrite(final_results, file = output_csv_path)

cat("\n--- Workflow Finished ---\n")
