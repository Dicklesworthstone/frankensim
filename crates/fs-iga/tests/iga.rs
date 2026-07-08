//! Battery for isogeometric analysis (fs-iga). Covers the B-spline partition
//! of unity, exact reproduction of a polynomial Poisson solution (the Galerkin
//! consistency property), high-order (k-refinement) convergence on a smooth
//! solution, DOF counting, and error paths.

use std::f64::consts::PI;

use fs_iga::{BsplineSpace, IgaError, solve_poisson};

#[test]
fn the_bspline_basis_is_a_partition_of_unity() {
    let space = BsplineSpace::clamped_uniform(3, 5);
    for &x in &[0.05, 0.3, 0.5, 0.73, 0.95] {
        let sum: f64 = (0..space.num_basis()).map(|i| space.basis(i, x)).sum();
        assert!(
            (sum - 1.0).abs() < 1e-12,
            "partition of unity at {x}: {sum}"
        );
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
    // reproduced to roundoff (Galerkin reproduces exactly-representable solutions).
    assert!(
        sol.l2_error(exact) < 1e-9,
        "L2 error {}",
        sol.l2_error(exact)
    );
    assert!((sol.eval(0.5) - 0.25).abs() < 1e-9);
    assert!((sol.eval(0.25) - 0.1875).abs() < 1e-9);
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
