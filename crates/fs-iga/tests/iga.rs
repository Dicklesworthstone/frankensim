//! Battery for isogeometric analysis (fs-iga). Covers the B-spline partition
//! of unity, exact reproduction of a polynomial Poisson solution (the Galerkin
//! consistency property), high-order (k-refinement) convergence on a smooth
//! solution, DOF counting, and error paths.

use std::f64::consts::PI;

use fs_iga::{BsplineSpace, IgaError, solve_poisson};

#[test]
fn the_bspline_basis_is_a_partition_of_unity() {
    // Σ Nᵢ(x) = 1 is a GLOBAL invariant (CONTRACT), so sweep the whole family
    // and — critically — sample the ENDPOINTS and the INTERIOR KNOTS, not just
    // span interiors as the old 5-point check did: a knot-span off-by-one in
    // Cox–de Boor is invisible away from knots but breaks the sum exactly at
    // them, and the clamped ends are the sharpest single-basis PoU points.
    // Worst |Σ-1| over degrees 1..=8 × elements 1..=6 × {0, 1, every interior
    // knot, a 41-pt grid} measured 7.77e-16 (~3.5 ulp) — BIT-IDENTICAL on
    // aarch64 (M4 Pro) and x86-64 (Threadripper 5975WX). Gate at 1e-13 (~130×
    // that roundoff floor), also tighter than the old 1e-12.
    for degree in 1..=8usize {
        for elements in 1..=6usize {
            let space = BsplineSpace::clamped_uniform(degree, elements);
            let mut xs: Vec<f64> = vec![0.0, 1.0];
            for e in 1..elements {
                xs.push(e as f64 / elements as f64); // interior knots
            }
            for k in 0..=40usize {
                xs.push(k as f64 / 40.0);
            }
            for x in xs {
                let sum: f64 = (0..space.num_basis()).map(|i| space.basis(i, x)).sum();
                assert!(
                    (sum - 1.0).abs() < 1e-13,
                    "partition of unity deg {degree} elem {elements} at {x}: {sum}"
                );
            }
        }
    }
}

#[test]
fn dof_count_is_degree_plus_elements() {
    assert_eq!(BsplineSpace::clamped_uniform(2, 4).num_basis(), 6);
    assert_eq!(BsplineSpace::clamped_uniform(3, 8).num_basis(), 11);
}

#[test]
fn a_polynomial_solution_is_reproduced_exactly() {
    // -u'' = 2, u(0)=u(1)=0  =>  u(x) = x(1-x) is in the degree-2 space.
    let space = BsplineSpace::clamped_uniform(2, 4);
    let sol = solve_poisson(&space, |_x| 2.0).unwrap();
    let exact = |x: f64| x * (1.0 - x);
    // Galerkin reproduces an exactly-representable solution TO ROUNDOFF (the
    // CONTRACT's consistency invariant), not merely "closely": measured L2 =
    // 7.0e-17 and pointwise ~1.1e-16, BIT-IDENTICAL on aarch64 (M4 Pro) and
    // x86-64 (Threadripper 5975WX). Gate at 1e-13 — 4 orders tighter than the
    // old 1e-9 smoke bound, which would have passed a real ~1e-10 systematic
    // consistency error (under-integration, a BC-lifting bug) while still
    // claiming "exact" — yet ~1000× above the roundoff floor, robust to benign
    // reordering of the dense solve.
    assert!(
        sol.l2_error(exact) < 1e-13,
        "L2 error {}",
        sol.l2_error(exact)
    );
    assert!((sol.eval(0.5) - 0.25).abs() < 1e-13);
    assert!((sol.eval(0.25) - 0.1875).abs() < 1e-13);
}

#[test]
fn k_refinement_converges_fast_on_a_smooth_solution() {
    // -u'' = π² sin(πx), u(0)=u(1)=0 => u(x) = sin(πx) (not polynomial).
    let g = |x: f64| PI * PI * (PI * x).sin();
    let exact = |x: f64| (PI * x).sin();
    let low = solve_poisson(&BsplineSpace::clamped_uniform(2, 4), g).unwrap();
    let high = solve_poisson(&BsplineSpace::clamped_uniform(4, 4), g).unwrap();
    // higher degree (the IGA k-refinement superpower) is dramatically more accurate.
    assert!(high.l2_error(exact) < low.l2_error(exact));
    assert!(
        high.l2_error(exact) < 1e-3,
        "deg-4 L2 error {}",
        high.l2_error(exact)
    );
}

#[test]
fn too_small_a_space_is_rejected() {
    // degree 1, 1 element -> 2 DOFs, both on the Dirichlet boundary.
    let tiny = BsplineSpace::clamped_uniform(1, 1);
    assert_eq!(tiny.num_basis(), 2);
    assert_eq!(solve_poisson(&tiny, |_| 1.0), Err(IgaError::TooFewDofs));
}

#[test]
fn the_solution_satisfies_the_boundary_conditions() {
    let space = BsplineSpace::clamped_uniform(3, 6);
    let sol = solve_poisson(&space, |x: f64| (PI * x).sin()).unwrap();
    assert!(sol.eval(0.0).abs() < 1e-12);
    assert!(sol.eval(1.0).abs() < 1e-9); // clamped-end interpolation
}

#[test]
fn solving_is_deterministic() {
    let space = BsplineSpace::clamped_uniform(3, 5);
    let a = solve_poisson(&space, |x: f64| (PI * x).sin()).unwrap();
    let b = solve_poisson(&space, |x: f64| (PI * x).sin()).unwrap();
    assert_eq!(a.coeffs(), b.coeffs());
}

#[test]
fn one_element_high_degree_spaces_are_not_underintegrated() {
    let exact = |x: f64| x * (1.0 - x);
    for degree in 2..=8 {
        let space = BsplineSpace::clamped_uniform(degree, 1);
        let solution = solve_poisson(&space, |_| 2.0)
            .unwrap_or_else(|error| panic!("degree {degree} solve failed: {error:?}"));
        assert!(
            solution.coeffs().iter().all(|value| value.is_finite()),
            "degree {degree} produced non-finite coefficients"
        );
        let error = solution.l2_error(exact);
        // reproduction stays at ROUNDOFF for every degree through 8 — the point
        // of the (p+1)-pt rule is that it keeps the one-element high-degree
        // space fully ranked (no under-integration). Worst measured L2 =
        // 1.57e-16 (degree 8), BIT-IDENTICAL on aarch64 and x86-64; gate at
        // 1e-13 (was 1e-9) so a real under-integration, which would cost
        // several orders here, is caught.
        assert!(
            error < 1e-13,
            "degree {degree} polynomial reproduction error {error}"
        );
    }
}
