//! G0/G1/G3/G4: bounded large SPD admission and original quadratic KKT gates.

use fs_solver::orthant::active_set::{
    MAX_SPD_VARIABLES, SpdOrthantConfig, SpdOrthantError, admit_spd_matrix, solve_spd_orthant,
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
