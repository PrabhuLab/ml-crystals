context("Geometric Property Calculation Functions")

# Use the pre-computed 'results_422' object from the helper file
# to test the output of each calculation step.

test_that("Bonding analysis (minimum_distance) is correct for ICSD422", {
  bonds <- results_422$bonded_pairs[[1]]

  # Ground truth from notebook is 14 bonds
  expect_equal(nrow(bonds), 14)

  # Check a specific bond: Si1 to Sr1_1_0_0_0 from notebook
  si1_sr1_bond <- bonds[Atom1 == "Si1" & Atom2 == "Sr1_1_0_0_0"]
  expect_equal(si1_sr1_bond$Distance, 3.163544, tolerance = 1e-6)
  expect_equal(si1_sr1_bond$DistanceError, 0.014160750, tolerance = 1e-6)
  expect_equal(si1_sr1_bond$dmin, 3.163544, tolerance = 1e-6)
})

test_that("Neighbor counting is correct for ICSD422", {
  counts <- results_422$neighbor_counts[[1]]

  expect_equal(nrow(counts), 3)
  expect_equal(counts[Atom == "Sr1", NeighborCount], 4)
  expect_equal(counts[Atom == "Sr2", NeighborCount], 3)
  expect_equal(counts[Atom == "Si1", NeighborCount], 7)
})

test_that("Bond angle calculations are correct for ICSD422", {
  angles <- results_422$bond_angles[[1]]

  # Ground truth from notebook is 30 unique angles
  expect_equal(nrow(angles), 30)

  # Check a specific angle from the notebook's sorted list
  angle_row <- angles[CentralAtom == "Sr1" &
                        ((Neighbor1 == "Si1_1_0_0_0" & Neighbor2 == "Si1_2_0_0_0") |
                           (Neighbor2 == "Si1_1_0_0_0" & Neighbor1 == "Si1_2_0_0_0"))]

  expect_equal(angle_row$Angle, 100.63957, tolerance = 1e-5)
  expect_equal(angle_row$AngleError, 0.3584927, tolerance = 1e-5)
})
