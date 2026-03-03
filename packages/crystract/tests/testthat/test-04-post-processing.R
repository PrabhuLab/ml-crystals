# Ensure default data is available
utils::data("covalent_radii", package = "crystract", envir = environment())

test_that("filter_ghost_distances identifies plausible and implausible distances",
          {
            # Inject dummy radii for Si and Sr to ensure test is robust regardless of default table
            my_radii <- data.table::data.table(
              Symbol = c("Si", "Sr"),
              Radius = c(1.11, 1.92),
              Type = c("covalent", "covalent")
            )
            set_radii_data(my_radii)
            on.exit(set_radii_data(NULL), add = TRUE)

            distances <- distances_1590946

            # Extract a real distance from the structure to test the filter
            real_bond <- distances[Atom1 == "Si1" &
                                     grepl("Sr", Atom2)][order(Distance)][1, ]

            filtered_plausible <- filter_ghost_distances(real_bond, atoms_1590946, margin = 0.5)
            expect_equal(nrow(filtered_plausible$kept), 1)
            expect_equal(nrow(filtered_plausible$removed), 0)

            # --- Test Case: Known short, implausible distance ---
            fake_ghost <- data.table::data.table(Atom1 = "Si1",
                                                 Atom2 = "Si1_1_0_0_0",
                                                 Distance = 0.5)
            filtered_short <- filter_ghost_distances(fake_ghost, atoms_1590946, margin = 0.2)
            expect_equal(nrow(filtered_short$removed), 1)
            expect_equal(filtered_short$removed$Reason, "Distance is TOO SHORT")
          })

test_that("calculate_weighted_average_network_distance is correct", {
  bonds_md <- results_1590946$bonds_minimum_distance[[1]]

  # The Wyckoff positions for this file are missing in the CIF, so mock them
  test_atoms <- data.table::copy(atoms_1590946)
  test_atoms[, WyckoffSymbol := "c"]
  test_atoms[, WyckoffMultiplicity := 4]

  # Calculate distance based on the mocked 4c site
  avg_dist <- calculate_weighted_average_network_distance(bonds_md, test_atoms, "4c")

  expect_type(avg_dist, "double")
  expect_true(!is.na(avg_dist))
  expect_true(avg_dist > 1.0 &&
                avg_dist < 5.0) # Reasonable bond length sanity check
})

test_that("set_radii_data custom atomic radii management works", {
  # Create a fake user radii table
  my_radii <- data.table::data.table(
    Symbol = c("Si", "Sr"),
    Radius = c(1.0, 2.0),
    Type = c("covalent", "covalent")
  )

  # Set the custom radii
  set_radii_data(my_radii)

  # Fetch internal active radii
  active_radii <- crystract:::get_radii_data()
  expect_equal(nrow(active_radii), 2)
  expect_equal(active_radii[Symbol == "Si", Radius], 1.0)

  # Reset to default
  set_radii_data(NULL)

  # Verify reset
  default_radii <- crystract:::get_radii_data()
  expect_true(nrow(default_radii) > 2)
})
