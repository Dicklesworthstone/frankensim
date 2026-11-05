//! G0/G1/G3/G4: bounded large SPD admission and original quadratic KKT gates.

use fs_solver::orthant::active_set::{
    MAX_SPD_VARIABLES, SpdOrthantConfig, SpdOrthantError, admit_spd_matrix, prepare_spd_matrix,
    solve_spd_orthant,
};

fn matrix(n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        a[i * n + i] = 2.5;
        if i > 0 {
            a[i * n + i - 1] = -0.5;
            a[(i - 1) * n + i] = -0.5;
        }
    }
    a
}

fn config(n: usize) -> SpdOrthantConfig {
    SpdOrthantConfig {
        max_variables: n.max(1),
        max_iterations: 16,
        ..SpdOrthantConfig::default()
    }
}

#[test]
fn g3_large_spd_admission_preserves_diagonal_equilibration_and_original_matrix() {
    let n = 192;
    let a = matrix(n);
    let before = a.clone();
    let admitted = admit_spd_matrix(&a, n, config(n), || false).unwrap();
    assert_eq!(admitted.dimension, n);
    assert!(admitted.min_normalized_pivot > 0.9);
    assert_eq!(a, before);
    let scales: Vec<_> = (0..n).map(|i| 2.0_f64.powi((i % 41) as i32 - 20)).collect();
    let mut transformed = a.clone();
    for i in 0..n {
        for j in 0..n {
            transformed[i * n + j] *= scales[i] * scales[j];
        }
    }
    let equivalent = admit_spd_matrix(&transformed, n, config(n), || false).unwrap();
    assert_eq!(admitted, equivalent);
}

#[test]
fn g1_large_orthant_problem_retains_manufactured_primal_and_dual_solution() {
    let n = 96;
    let a = matrix(n);
    let mut expected = vec![0.0; n];
    expected[1] = 1.0;
    expected[n / 2] = 2.0;
    expected[n - 2] = 0.5;
    let gradient: Vec<_> = expected
        .iter()
        .map(|&x| if x == 0.0 { 1.0 } else { 0.0 })
        .collect();
    let linear: Vec<_> = a
        .chunks_exact(n)
        .zip(&gradient)
        .map(|(row, &g)| g - row.iter().zip(&expected).map(|(a, x)| a * x).sum::<f64>())
        .collect();
    let solution = solve_spd_orthant(&a, &linear, config(n), || false).unwrap();
    for i in 0..n {
        assert!((solution.point[i] - expected[i]).abs() < 2e-13);
        assert!((solution.gradient[i] - gradient[i]).abs() < 2e-13);
    }
    assert!(solution.residual < 2e-13);
    assert!(solution.complementarity_residual < 2e-13);
}

#[test]
fn g0_large_admission_keeps_capacity_symmetry_and_positive_pivot_gates() {
    let empty = admit_spd_matrix(&[], 0, config(MAX_SPD_VARIABLES), || false).unwrap();
    assert_eq!(empty.dimension, 0);
    assert_eq!(empty.min_normalized_pivot, 1.0);
    assert!(matches!(
        admit_spd_matrix(
            &[],
            MAX_SPD_VARIABLES + 1,
            config(MAX_SPD_VARIABLES),
            || false
        ),
        Err(SpdOrthantError::Invalid(_))
    ));
    let n = 65;
    let a = matrix(n);
    assert!(matches!(
        admit_spd_matrix(&a, n, SpdOrthantConfig::default(), || false),
        Err(SpdOrthantError::Invalid(_))
    ));
    let mut bad = a.clone();
    bad[(n - 1) * n + n - 1] = 0.0;
    assert!(matches!(
        admit_spd_matrix(&bad, n, config(n), || false),
        Err(SpdOrthantError::NotPositiveDefinite { .. })
    ));
    bad = a.clone();
    bad[1] = -0.4;
    assert!(matches!(
        admit_spd_matrix(&bad, n, config(n), || false),
        Err(SpdOrthantError::Invalid(_))
    ));
    bad = a;
    bad[0] = 0.01;
    assert!(matches!(
        admit_spd_matrix(&bad, n, config(n), || false),
        Err(SpdOrthantError::NotPositiveDefinite { .. })
    ));
}

#[test]
fn g4_large_admission_cancels_during_normalization_factorization_and_publication() {
    let n = 97;
    let a = matrix(n);
    let before = a.clone();
    let mut calls = 0;
    admit_spd_matrix(&a, n, config(n), || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, n + 3, calls / 2, calls] {
        let mut count = 0;
        assert_eq!(
            admit_spd_matrix(&a, n, config(n), || {
                count += 1;
                count == stop
            }),
            Err(SpdOrthantError::Cancelled)
        );
        assert_eq!(a, before);
    }
}

#[test]
fn g1_prepared_spd_factor_solves_repeated_sources_in_original_scaled_units() {
    let n = 96;
    let scales: Vec<_> = (0..n).map(|i| 2.0_f64.powi((i % 13) as i32 - 6)).collect();
    let mut a = matrix(n);
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] *= scales[i] * scales[j];
        }
    }
    let prepared = prepare_spd_matrix(&a, n, config(n), || false).unwrap();
    assert_eq!(
        prepared.admission(),
        admit_spd_matrix(&a, n, config(n), || false).unwrap()
    );
    for offset in [0, 4, 11] {
        let expected: Vec<_> = (0..n)
            .map(|i| ((i * 7 + offset) % 17) as f64 / 8.0 - 1.0)
            .collect();
        let rhs: Vec<_> = a
            .chunks_exact(n)
            .map(|row| row.iter().zip(&expected).map(|(a, x)| a * x).sum())
            .collect();
        let before = rhs.clone();
        let solved = prepared.solve(&rhs, || false).unwrap();
        assert_eq!(rhs, before);
        for (actual, expected) in solved.iter().zip(expected) {
            assert!((actual - expected).abs() < 2e-11, "{actual} != {expected}");
        }
        assert_eq!(solved, prepared.solve(&rhs, || false).unwrap());
    }
}

#[test]
fn g4_prepared_spd_solve_refuses_invalid_sources_and_retries_after_cancellation() {
    let n = 65;
    let prepared = prepare_spd_matrix(&matrix(n), n, config(n), || false).unwrap();
    let rhs: Vec<_> = (0..n).map(|i| i as f64 / 8.0 - 4.0).collect();
    let before = rhs.clone();
    let mut calls = 0;
    let expected = prepared
        .solve(&rhs, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert_eq!(
            prepared.solve(&rhs, || {
                count += 1;
                count == stop
            }),
            Err(SpdOrthantError::Cancelled)
        );
        assert_eq!(rhs, before);
        assert_eq!(prepared.solve(&rhs, || false).unwrap(), expected);
    }
    assert!(matches!(
        prepared.solve(&rhs[..n - 1], || false),
        Err(SpdOrthantError::Invalid(_))
    ));
    assert!(matches!(
        prepared.solve(&vec![f64::NAN; n], || false),
        Err(SpdOrthantError::Invalid(_))
    ));
    let tiny = prepare_spd_matrix(&[1e-200], 1, config(1), || false).unwrap();
    assert_eq!(
        tiny.solve(&[1e300], || false),
        Err(SpdOrthantError::NonFinite)
    );
    let empty = prepare_spd_matrix(&[], 0, config(1), || false).unwrap();
    assert_eq!(empty.solve(&[], || false).unwrap(), Vec::<f64>::new());
}
