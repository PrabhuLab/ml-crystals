# Define global variables to satisfy R CMD check
# These are typically column names used in data.table operations
# or other symbols that appear as unbound global variables.

utils::globalVariables(c(
  "V1", "x_a", "y_b", "z_c", "Label", "Distance", "Atom1", "Atom2",
  "dmin", "dcut", ".", ":=", ".N", "x", "y", "z", "BondStrength",
  "NeighborCount", "CentralAtom", "Neighbor1", "Neighbor2", "Angle",
  "DeltaX", "DeltaY", "DeltaZ", "CosAlpha", "CosBeta", "CosGamma",
  "_cell_length_a", "_cell_length_b", "_cell_length_c",
  "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"
))
