test_that("apply_symmetry_operations generates correct full cell for 1590946",
          {
            # Ground truth is 12 atoms in the full unit cell.
            # The space group is P n m a (general multiplicity 8), but the atoms in this
            # structure are on the 4c Wyckoff site (y=0.25). So 3 atoms * 4 = 12 atoms.
            expect_equal(nrow(full_cell_1590946), 12)

            # Check a specific transformed coordinate against the ground truth data.
            # The atom 'Sr1' (0.6529, 0.25, 0.0769) transformed by the 1st sym op ('x+1/2, y, -z+1/2')
            # becomes (1.1529, 0.25, 0.4231) -> wrapped to (0.1529, 0.25, 0.4231)
            sr1_transformed <- full_cell_1590946[abs(x_a - 0.1529) < 1e-4 &
                                                   abs(z_c - 0.4231) < 1e-4]
            expect_true(nrow(sr1_transformed) >= 1)
          })

test_that("expand_transformed_coords generates correct supercell for 1590946",
          {
            # Ground truth is 12 atoms/cell * (3x3x3 cells = 27) = 324 atoms
            expect_equal(nrow(super_cell_1590946), 324)
          })

test_that(
  "calculate_expansion_factors correctly determines padding for geometric algorithms",
  {
    # For a 13.0 Angstrom Voronoi cutoff on an 8.11x5.15x9.54 unit cell,
    # expansion factors should safely cover it (e.g., cell b is 5.15, 13 / 5.15 > 2.5 => ceiling + buffer => ~3)
    factors <- calculate_expansion_factors(metrics_1590946, radius = 13.0)
    expect_length(factors, 3)
    expect_true(all(factors >= 1))
  }
)
