use super::*;

fn apply(lower: &[f64], diagonal: &[f64], upper: &[f64], x: &[f64]) -> Vec<f64> {
    (0..x.len())
        .map(|i| {
            let mut value = diagonal[i] * x[i];
            if i > 0 {
                value += lower[i - 1] * x[i - 1];
            }
            if i + 1 < x.len() {
                value += upper[i] * x[i + 1];
            }
            value
        })
        .collect()
}

fn assert_solution(lower: &[f64], diagonal: &[f64], upper: &[f64], expected: &[f64]) {
    let b = apply(lower, diagonal, upper, expected);
    let mut solved = vec![f64::NAN; expected.len()];
    let mut work = TridiagonalWorkspace::new(expected.len()).unwrap();
    work.solve_into(lower, diagonal, upper, &b, &mut solved, || false)
        .unwrap();
    let actual_b = apply(lower, diagonal, upper, &solved);
    for i in 0..expected.len() {
        assert!((solved[i] - expected[i]).abs() <= 2e-12 * (1.0 + expected[i].abs()));
        assert!((actual_b[i] - b[i]).abs() <= 2e-13 * (1.0 + b[i].abs()));
    }
}

#[test]
fn adjacent_pivots_retain_second_superdiagonal_and_solve_original_nonsymmetric_system() {
    assert_solution(&[3.0], &[0.0, 4.0], &[2.0], &[-2.0, 3.0]);
    assert_solution(
        &[2.0, 4.0, 6.0],
        &[0.0; 4],
        &[1.0, 3.0, 5.0],
        &[1.0, -2.0, 3.0, -4.0],
    );
    assert_solution(
        &[4.0, -5.0, 6.0, 3.0],
        &[0.0, 1.0, -1.0, 0.0, 2.0],
        &[2.0, 7.0, 8.0, -4.0],
        &[1.25, -2.5, 0.75, 9.0, -0.125],
    );
    // Pivot ties retain the current row, matching the deterministic contract.
    let mut work = TridiagonalWorkspace::new(2).unwrap();
    work.factor(&[-2.0], &[2.0, 3.0], &[1.0], || false).unwrap();
    assert!(!work.swapped[0]);
}

#[test]
fn one_prepared_factor_serves_multiple_rhs_and_refactoring_reuses_linear_storage() {
    let n = 513;
    // A positive scalar identity plus a real skew-symmetric matrix stays
    // well conditioned while exercising many adjacent pivots without diagonal
    // dominance. Unequal off-diagonals can instead hide exponential scaling.
    let lower = vec![2.0; n - 1];
    let diagonal = vec![0.5; n];
    let upper = vec![-2.0; n - 1];
    let mut work = TridiagonalWorkspace::new(n).unwrap();
    let pointers = (
        work.diagonal.as_ptr(),
        work.lower.as_ptr(),
        work.upper.as_ptr(),
        work.upper_two.as_ptr(),
        work.rhs.as_ptr(),
    );
    work.factor(&lower, &diagonal, &upper, || false).unwrap();
    assert!(work.swapped.iter().any(|&swapped| swapped));
    for scale in [0.0, 0.25, -7.0] {
        let expected: Vec<_> = (0..n).map(|i| scale * (i % 11) as f64).collect();
        let rhs = apply(&lower, &diagonal, &upper, &expected);
        let mut in_place = rhs.clone();
        work.solve_factored_in_place(&mut in_place, || false)
            .unwrap();
        let restored = apply(&lower, &diagonal, &upper, &in_place);
        for (actual, b) in restored.iter().zip(&rhs) {
            assert!((actual - b).abs() <= 2e-12 * (1.0 + b.abs()));
        }
        let mut separate = vec![999.0; n];
        work.solve_factored_into(&rhs, &mut separate, || false)
            .unwrap();
        assert_eq!(in_place, separate);
    }
    let fresh = vec![4.0; n];
    work.factor(&lower, &fresh, &upper, || false).unwrap();
    assert_eq!(
        pointers,
        (
            work.diagonal.as_ptr(),
            work.lower.as_ptr(),
            work.upper.as_ptr(),
            work.upper_two.as_ptr(),
            work.rhs.as_ptr(),
        )
    );
    assert_eq!(work.dimension(), n);
}

#[test]
fn empty_scalar_scaled_and_singular_systems_have_explicit_semantics() {
    let mut empty = TridiagonalWorkspace::new(0).unwrap();
    empty
        .solve_into(&[], &[], &[], &[], &mut [], || false)
        .unwrap();
    assert_solution(&[], &[4.0], &[], &[-2.5]);
    // Small finite units do not imply singularity or justify an arbitrary floor.
    assert_solution(&[], &[1e-300], &[], &[0.5]);
    assert_solution(&[1e-200], &[4e-200, 3e-200], &[-2e-200], &[2.0, -3.0]);
    let mut work = TridiagonalWorkspace::new(2).unwrap();
    work.factor(&[1.0], &[3.0, 4.0], &[1.0], || false).unwrap();
    assert_eq!(
        work.factor(&[2.0], &[1.0, 4.0], &[2.0], || false),
        Err(TridiagonalError::Singular { index: 1 })
    );
    let mut out = [5.0, 6.0];
    assert_eq!(
        work.solve_factored_into(&[1.0, 2.0], &mut out, || false),
        Err(TridiagonalError::Unfactored)
    );
    assert_eq!(out, [5.0, 6.0]);
    assert_eq!(
        work.factor(&[0.0], &[0.0, 1.0], &[1.0], || false),
        Err(TridiagonalError::Singular { index: 0 })
    );
    assert!(matches!(
        TridiagonalWorkspace::new(usize::MAX),
        Err(TridiagonalError::Capacity)
    ));
}

#[test]
fn invalid_data_and_overflow_never_publish_output_or_reuse_rejected_factors() {
    let mut work = TridiagonalWorkspace::new(2).unwrap();
    let mut out = [4.0, 5.0];
    assert_eq!(
        work.solve_into(&[], &[1.0, 1.0], &[0.0], &[1.0, 2.0], &mut out, || false),
        Err(TridiagonalError::Dimension)
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            work.factor(&[0.0], &[1.0, bad], &[0.0], || false),
            Err(TridiagonalError::NonFinite)
        );
        assert!(!work.factored);
    }
    assert_eq!(
        work.factor(&[f64::MAX], &[f64::MAX, f64::MAX], &[-f64::MAX], || false),
        Err(TridiagonalError::NonFinite)
    );
    work.factor(&[0.0], &[1e-300, 1.0], &[0.0], || false)
        .unwrap();
    assert_eq!(
        work.solve_factored_into(&[1e300, 0.0], &mut out, || false),
        Err(TridiagonalError::NonFinite)
    );
    assert_eq!(out, [4.0, 5.0]);
    // A refused RHS does not invalidate a previously successful factor.
    work.solve_factored_into(&[1e-300, 2.0], &mut out, || false)
        .unwrap();
    assert_eq!(out, [1.0, 2.0]);
    let mut rhs = [f64::NAN, -0.0];
    let before = rhs.map(f64::to_bits);
    assert_eq!(
        work.solve_factored_in_place(&mut rhs, || false),
        Err(TridiagonalError::NonFinite)
    );
    assert_eq!(rhs.map(f64::to_bits), before);
}

#[test]
fn bounded_cancellation_preserves_output_at_each_factor_and_substitution_stage() {
    let n = 259;
    let lower = vec![-1.0; n - 1];
    let diagonal = vec![4.0; n];
    let upper = vec![-2.0; n - 1];
    let b = vec![1.0; n];
    let mut work = TridiagonalWorkspace::new(n).unwrap();
    let mut factor_polls = 0;
    work.factor(&lower, &diagonal, &upper, || {
        factor_polls += 1;
        false
    })
    .unwrap();
    assert!(factor_polls >= 8);
    for stop in 1..=factor_polls {
        let mut polls = 0;
        assert_eq!(
            work.factor(&lower, &diagonal, &upper, || {
                polls += 1;
                polls == stop
            }),
            Err(TridiagonalError::Cancelled)
        );
        assert!(!work.factored);
    }
    work.factor(&lower, &diagonal, &upper, || false).unwrap();
    let mut solve_polls = 0;
    let mut solved = b.clone();
    work.solve_factored_in_place(&mut solved, || {
        solve_polls += 1;
        false
    })
    .unwrap();
    assert!(solve_polls >= 11);
    for stop in 1..=solve_polls {
        let mut polls = 0;
        let mut rhs = b.clone();
        assert_eq!(
            work.solve_factored_in_place(&mut rhs, || {
                polls += 1;
                polls == stop
            }),
            Err(TridiagonalError::Cancelled)
        );
        assert_eq!(rhs, b);
        assert!(work.factored);
    }
    work.solve_factored_in_place(&mut solved, || false).unwrap();
}
