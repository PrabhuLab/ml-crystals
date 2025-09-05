# This function loads data from a package, which is necessary for testing.
utils::data("covalent_radii", package = "crystract", envir = environment())

test_that("filter_ghost_distances identifies plausible and implausible distances",
          {
            distances <- results_422$distances[[1]]

            # --- Test Case 1: Known plausible bond ---
            # Si-Sr bond distance is ~3.16 Å. Covalent radii are Si=1.11, Sr=1.92. Sum = 3.03 Å.
            # With a 20% margin, the plausible range is [2.42, 3.64]. The 3.16 Å bond should be kept.
            plausible_bond <- distances[Atom1 == "Si1" &
                                          Atom2 == "Sr1_1_0_0_0"] # Distance is 3.163544
            filtered_plausible <- filter_ghost_distances(plausible_bond, atoms_422, margin = 0.2)
            expect_equal(nrow(filtered_plausible$kept), 1)
            expect_equal(nrow(filtered_plausible$removed), 0)

            # --- Test Case 2: Known long, non-bonding distance ---
            # A long distance between two Sr atoms should be removed. Let's find one.
            long_distance <- distances[Atom1 == "Sr1" &
                                         Distance > 4.5][1, ] # e.g., Sr1 to Sr1_3_0_0_0 is 4.77 Å
            filtered_long <- filter_ghost_distances(long_distance, atoms_422, margin = 0.2)
            expect_equal(nrow(filtered_long$removed), 1)
            expect_equal(filtered_long$removed$Reason, "Distance is TOO LONG")

            # --- Test Case 3: Known short, implausible distance ---
            fake_ghost <- data.table(Atom1 = "Si1",
                                     Atom2 = "Sr1_1_0_0_0",
                                     Distance = 0.5)
            filtered_short <- filter_ghost_distances(fake_ghost, atoms_422, margin = 0.2)
            expect_equal(nrow(filtered_short$removed), 1)
            expect_equal(filtered_short$removed$Reason, "Distance is TOO SHORT")
          })

test_that("calculate_weighted_average_network_distance is correct", {
  bonds_md <- results_422$bonded_pairs_minimum_distance[[1]]
  avg_dist <- calculate_weighted_average_network_distance(bonds_md, atoms_422, "4c")

  # The manual run confirms the correct value is ~3.281 Å.
  # The test should check against this known correct value.
  expect_equal(avg_dist, 3.281382, tolerance = 1e-6)

  # The dynamic calculation confirms the function's internal logic is self-consistent.
  atom_sums <- bonds_md[, .(sum_d = sum(Distance), n = .N), by = Atom1]
  atom_sums[, ParentLabel := sub("_.*", "", Atom1)]
  merged <- merge(atoms_422,
                  unique(atom_sums[, .(ParentLabel, sum_d, n)]),
                  by.x = "Label",
                  by.y = "ParentLabel")
  num <- sum(merged$WyckoffMultiplicity * merged$Occupancy * merged$sum_d)
  den <- sum(merged$WyckoffMultiplicity * merged$Occupancy * merged$n)
  expected_avg <- num / den

  expect_equal(avg_dist, expected_avg, tolerance = 1e-6)
})
