//! Battery for symmetry harvesting (addendum Proposal 13). Covers the cyclic
//! character table, the circulant matvec/solve round-trip (block-diagonal DFT
//! solve, including a non-symmetric operator to pin the eigenvalue sign),
//! singular detection, the certified asymmetry residual, isotypic
//! symmetrization, and the certified perturbation bound for approximate
//! symmetry, plus determinism.

use fs_symmetry::{
    CyclicGroup, SymmetryError, circulant_matvec, cyclic_residual, solve_circulant, symmetrize,
    symmetrized_solve,
};

fn dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn the_cyclic_character_table_is_roots_of_unity() {
    let g = CyclicGroup::new(4);
    // trivial irrep is 1 everywhere.
    for e in 0..4 {
        let (re, im) = g.character(0, e);
        assert!((re - 1.0).abs() < 1e-12 && im.abs() < 1e-12);
    }
    // irrep 1 at element 1 = i.
    let (re, im) = g.character(1, 1);
    assert!(re.abs() < 1e-12 && (im - 1.0).abs() < 1e-12);
    // irrep 2 at element 1 = -1.
    let (re, im) = g.character(2, 1);
    assert!((re + 1.0).abs() < 1e-12 && im.abs() < 1e-12);
}

#[test]
fn the_block_diagonal_solve_inverts_the_circulant() {
    // a NON-symmetric, diagonally-dominant circulant (pins the eigenvalue sign).
    let c = [4.0, 1.0, 2.0, 0.0];
    let b = [1.0, -2.0, 3.0, 0.5];
    let x = solve_circulant(&c, &b).unwrap();
    // C x must reproduce b (the symmetry-adapted solve equals the full solve).
    let cx = circulant_matvec(&c, &x).unwrap();
    assert!(dist(&cx, &b) < 1e-9, "residual {}", dist(&cx, &b));
}

#[test]
fn the_solve_scales_to_a_larger_ring() {
    // a 12-node ring (3-fold friendly): symmetric, diagonally dominant.
    let mut c = vec![0.0; 12];
    c[0] = 6.0;
    c[1] = -1.0;
    c[11] = -1.0;
    let b: Vec<f64> = (0..12).map(|i| (f64::from(i) - 5.0).sin()).collect();
    let x = solve_circulant(&c, &b).unwrap();
    assert!(dist(&circulant_matvec(&c, &x).unwrap(), &b) < 1e-9);
}

#[test]
fn a_singular_circulant_is_rejected() {
    // all-ones first row -> a zero eigenvalue.
    assert_eq!(
        solve_circulant(&[1.0, 1.0, 1.0, 1.0], &[1.0, 0.0, 0.0, 0.0]),
        Err(SymmetryError::Singular)
    );
    // shape errors.
    assert!(matches!(
        solve_circulant(&[1.0, 2.0], &[1.0]),
        Err(SymmetryError::LengthMismatch { .. })
    ));
    assert_eq!(solve_circulant(&[], &[]), Err(SymmetryError::EmptyInput));
}

#[test]
fn the_asymmetry_residual_is_certified() {
    // exactly 2-fold symmetric: [1,2,1,2] is invariant under a shift by 2.
    let r = cyclic_residual(&[1.0, 2.0, 1.0, 2.0], 2).unwrap();
    assert!(r.is_exact && r.residual < 1e-12);
    // not symmetric -> a positive residual.
    let a = cyclic_residual(&[1.0, 2.0, 3.0, 4.0], 2).unwrap();
    assert!(!a.is_exact && a.residual > 0.0);
    // error paths.
    assert_eq!(cyclic_residual(&[], 2), Err(SymmetryError::EmptyInput));
    assert_eq!(
        cyclic_residual(&[1.0, 2.0, 3.0], 2),
        Err(SymmetryError::NotDivisible { len: 3, k_fold: 2 })
    );
    assert_eq!(
        cyclic_residual(&[1.0, 2.0], 0),
        Err(SymmetryError::ZeroFold)
    );
}

#[test]
fn symmetrize_projects_onto_the_symmetric_subspace() {
    let (sym, asym) = symmetrize(&[1.0, 2.0, 3.0, 4.0], 2).unwrap();
    assert_eq!(sym, vec![2.0, 3.0, 2.0, 3.0]);
    assert_eq!(asym, vec![-1.0, -1.0, 1.0, 1.0]);
    // the symmetric part is now EXACTLY 2-fold symmetric.
    assert!(cyclic_residual(&sym, 2).unwrap().is_exact);
    // sym + asym reconstructs the original.
    for i in 0..4 {
        assert!((sym[i] + asym[i] - [1.0, 2.0, 3.0, 4.0][i]).abs() < 1e-12);
    }
}

#[test]
fn approximate_symmetry_yields_a_bound_that_contains_the_true_correction() {
    let c = [4.0, 1.0, 2.0, 1.0]; // symmetric, well-conditioned (lambda_min = 2)
    let rhs = [1.0, 2.0, 3.0, 8.0]; // NOT 2-fold symmetric
    let pb = symmetrized_solve(&c, &rhs, 2).unwrap();
    // the true full solution and the symmetric-part solution.
    let x_full = solve_circulant(&c, &rhs).unwrap();
    let true_correction = dist(&x_full, &pb.symmetric_solution);
    assert!(pb.asymmetry_residual > 0.0);
    // the certificate must CONTAIN the true correction.
    assert!(
        true_correction <= pb.correction_bound + 1e-9,
        "correction {true_correction} exceeded bound {}",
        pb.correction_bound
    );
}

#[test]
fn solving_is_deterministic() {
    let c = [4.0, 1.0, 2.0, 0.0];
    let b = [1.0, -2.0, 3.0, 0.5];
    assert_eq!(solve_circulant(&c, &b), solve_circulant(&c, &b));
}
