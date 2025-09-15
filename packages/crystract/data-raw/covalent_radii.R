# Script to create and save the covalent_radii data object for the package.

# --- 0. Setup ---
# Load the necessary library
library(data.table)

# --- 1. Load Raw Data ---
file_path <- "covalent_radii_emsley.csv"
if (!file.exists(file_path)) {
  stop(
    "Error: 'covalent_radii_emsley.csv' not found. Please place it in the package's root directory."
  )
}
radii_dt <- fread(file_path)

# --- 2. Clean and Prepare Data ---
# Use data.table::setnames to be explicit and avoid errors
data.table::setnames(radii_dt,
                     old = names(radii_dt)[1:2],
                     new = c("Symbol", "Radius"))

# Select only the needed columns
radii_dt <- radii_dt[, .(Symbol, Radius)]

# Trim whitespace and convert 0 or NA radii to NA
radii_dt[, Symbol := trimws(Symbol)]
radii_dt[Radius == 0 | is.na(Radius), Radius := NA_real_]

# Remove rows with NA radii to create the final object
covalent_radii <- na.omit(radii_dt)

# --- 3. Add 'Type' Column ---
# This is essential for the new set_radii_data() functionality to work.
covalent_radii[, Type := "covalent"]

# --- 4. Save Final Package Object using Best Practices ---
# This function takes the 'covalent_radii' object currently in memory
# and saves it correctly inside the 'data/' folder for the package.
usethis::use_data(covalent_radii, overwrite = TRUE)

# --- 5. Confirmation Message ---
cat(
  "Success! The 'covalent_radii' object has been saved to 'data/covalent_radii.rda' with",
  nrow(covalent_radii),
  "entries.\n"
)
print("First few rows of the final data object:")
print(head(covalent_radii))
