test_that("Chemical bonding algorithms execute and return correct structure",
          {
            dists <- distances_1590946

            # Minimum Distance
            bonds_md <- minimum_distance(dists, delta = 0.1)
            expect_s3_class(bonds_md, "data.table")
            expect_true(nrow(bonds_md) > 0)
            expect_true("Weight" %in% names(bonds_md))

            # Brunner's Reciprocal
            bonds_br <- brunner_nn_reciprocal(dists)
            expect_s3_class(bonds_br, "data.table")
            expect_true(nrow(bonds_br) > 0)

            # Hoppe's EconNN
            bonds_ec <- econ_nn(dists, atoms_1590946, use_fictive_radius = TRUE)
            expect_s3_class(bonds_ec, "data.table")
            expect_true(nrow(bonds_ec) > 0)
          })

test_that("Geometric bonding algorithms (Voronoi/CrystalNN) execute successfully",
          {
            skip_if_not_installed("geometry")

            # Prepare collapsed coordinate data needed for geometric algorithms
            col_atoms <- crystract:::merge_partial_occupancy_sites(atoms_1590946)
            col_full <- apply_symmetry_operations(col_atoms, sym_ops_1590946, metrics_1590946)

            exp_factors <- calculate_expansion_factors(metrics_1590946, 13.0)
            col_super <- expand_transformed_coords(col_full, exp_factors)

            # Voronoi
            bonds_vo <- voronoi_nn(col_atoms, col_super, metrics_1590946)
            expect_s3_class(bonds_vo, "data.table")
            expect_true("SolidAngle" %in% names(bonds_vo))

            # CrystalNN
            bonds_cr <- crystal_nn(distances_1590946, col_atoms, col_super, metrics_1590946)
            expect_s3_class(bonds_cr, "data.table")
            expect_true("Weight" %in% names(bonds_cr))
          })

test_that("Neighbor counting is correct", {
  counts <- results_1590946$cn_minimum_distance[[1]]
  expect_equal(nrow(counts), 3)
  expect_true("WeightedCoordinationNumber" %in% names(counts))
})

test_that("Bond angle calculations run", {
  angles <- results_1590946$angles_minimum_distance[[1]]
  expect_s3_class(angles, "data.table")
  expect_true("Angle" %in% names(angles))
  expect_true(nrow(angles) > 0)
})
