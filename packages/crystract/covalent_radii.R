# Script to create the final data file from Julia's recommended table (Emsley, 1998)
library(data.table)

# --- 1. Load Data ---
file_path <- "covalent_radii_emsley.csv"
if (!file.exists(file_path)) {
  stop("Error: 'covalent_radii_emsley.csv' not found. Please create it from Julia's Excel file and save it in the package's root directory.")
}
radii_dt <- fread(file_path)

# --- 2. Clean Data ---
setnames(radii_dt, old = names(radii_dt)[1:2], new = c("Symbol", "Radius"))
radii_dt <- radii_dt[, .(Symbol, Radius)]
radii_dt[, Symbol := trimws(Symbol)]
radii_dt[Radius == 0 | is.na(Radius), Radius := NA_real_]
covalent_radii <- na.omit(radii_dt)

# --- 3. Save Final Package Object ---
dir.create("data", showWarnings = FALSE)
save(covalent_radii, file = "data/covalent_radii.rda")

cat("Success! 'data/covalent_radii.rda' has been updated with", nrow(covalent_radii), "entries.\n")
