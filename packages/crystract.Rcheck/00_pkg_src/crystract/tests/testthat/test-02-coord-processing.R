test_that("apply_symmetry_operations generates correct full cell for ICSD422",
          {
            # Ground truth is 12 atoms in the full unit cell
            expect_equal(nrow(full_cell_422), 12)

            # Check a specific transformed coordinate against the ground truth data.
            # The atom 'Sr1' transformed by the 2nd sym op ('x, -y+1/2, z')
            # becomes 'Sr1_2' with coordinates (0.6529, 0.25, 0.0769).
            sr1_2 <- full_cell_422[Label == "Sr1_2"]
            expect_equal(sr1_2$x_a, 0.6529, tolerance = 1e-6)
            expect_equal(sr1_2$y_b, 0.2500, tolerance = 1e-6)
            expect_equal(sr1_2$z_c, 0.0769, tolerance = 1e-6)
          })

test_that("expand_transformed_coords generates correct supercell for ICSD422",
          {
            # Ground truth is 12 atoms/cell * (3^3 cells) = 324 atoms
            expect_equal(nrow(super_cell_422), 324)
          })
