# Script to create and save the internal data objects for the package.
# This compiles the data into R/sysdata.rda so it is instantly available internally.

library(data.table)

# --- 1. Process Pauling Electronegativity ---
file_path_en <- "~/Documents/GitHub/ml-crystals/packages/crystract/data-raw/pauling_electronegativity_stable_elements.csv"
if (!file.exists(file_path_en)) stop(paste("Error:", file_path_en, "not found."))

en_dt <- fread(file_path_en)
# Select necessary columns and fix classes
pauling_en <- en_dt[, .(Symbol, PaulingElectronegativity)]
pauling_en[PaulingElectronegativity == "", PaulingElectronegativity := NA]
pauling_en[, PaulingElectronegativity := as.numeric(PaulingElectronegativity)]


# --- 2. Process Shannon Radii ---
file_path_shannon <- "~/Documents/GitHub/ml-crystals/packages/crystract/data-raw/Shannon_Radii.csv"
if (!file.exists(file_path_shannon)) stop(paste("Error:", file_path_shannon, "not found."))

shannon_dt <- fread(file_path_shannon)

# Rename columns to standard internal names
setnames(shannon_dt,
         old = c("Element", "Charge", "Ionic Radius"),
         new = c("Symbol", "OxidationState", "Radius"),
         skip_absent = TRUE)

# Clean Coordination (Convert Roman Numerals to Integers)
roman_map <- c(I=1, II=2, III=3, IV=4, V=5, VI=6, VII=7, VIII=8, IX=9, X=10, XI=11, XII=12, XIII=13, XIV=14)

clean_coord <- function(x) {
  x_clean <- gsub("[^A-Za-z0-9]", "", as.character(x))
  if (grepl("^[0-9]+$", x_clean)) return(as.numeric(x_clean))
  if (x_clean %in% names(roman_map)) return(roman_map[[x_clean]])
  num_part <- gsub("[^0-9]", "", as.character(x))
  if (nchar(num_part) > 0) return(as.numeric(num_part))
  return(NA_real_)
}

shannon_dt[, Coordination := sapply(Coordination, clean_coord)]
shannon_dt[, Radius := as.numeric(Radius)]

# Select final columns
shannon_radii <- shannon_dt[!is.na(Radius), .(Symbol, OxidationState, Coordination, Radius)]

# --- 3. Save as Internal Package Data ---
# Using internal = TRUE creates R/sysdata.rda.
# These variables will be available exclusively to internal package functions.
usethis::use_data(pauling_en, shannon_radii, overwrite = TRUE, internal = TRUE)

cat("Success! 'pauling_en' and 'shannon_radii' saved to R/sysdata.rda\n")
