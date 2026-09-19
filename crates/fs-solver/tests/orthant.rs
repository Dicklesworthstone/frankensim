use fs_solver::orthant::{solve_gram_orthant, GramOrthantConfig, GramOrthantError};
fn config() -> GramOrthantConfig {
    GramOrthantConfig { max_variables: 8, max_columns: 8, max_setup_products: 512,
        max_sweeps: 512, gradient_tolerance: 1e-12 }
}
fn near(a: f64, b: f64) { assert!((a-b).abs() < 1e-10, "{a} != {b}"); }

#[test]
fn shared_chain_matches_the_full_analytic_solution() {
    let rows = vec![vec![-1.0, 1.0, 0.0], vec![0.0, -1.0, 1.0]];
    let r = solve_gram_orthant(&rows, &[-1.0, -2.0], config(), || false).unwrap();
    near(r.point[0], 4.0/3.0); near(r.point[1], 5.0/3.0);
    assert!(r.sweeps > 1 && r.residual <= config().gradient_tolerance);
    near(2.0*r.point[0]-r.point[1]-1.0, r.gradient[0]);
    near(2.0*r.point[1]-r.point[0]-2.0, r.gradient[1]);
}
#[test]
fn inactive_variables_and_initial_optimum_are_not_forced_active() {
    let rows = vec![vec![1.0, 0.0], vec![0.0, 2.0]];
    let r = solve_gram_orthant(&rows, &[1.0, -8.0], config(), || false).unwrap();
    assert_eq!(r.point, [0.0, 2.0]); assert_eq!(r.gradient, [1.0, 0.0]);
    assert_eq!(solve_gram_orthant(&rows, &[1.0, 2.0], config(), || false).unwrap().sweeps, 0);
    assert!(solve_gram_orthant(&[], &[], config(), || false).unwrap().point.is_empty());
}
#[test]
fn redundant_rows_solve_without_a_false_uniqueness_claim() {
    let rows = vec![vec![1.0, -1.0], vec![1.0, -1.0]];
    let r = solve_gram_orthant(&rows, &[-2.0, -2.0], config(), || false).unwrap();
    near(r.point.iter().sum(), 1.0); assert!(r.residual <= 1e-12);
}
#[test]
fn cancellation_at_each_observed_boundary_and_retry_are_deterministic() {
    let rows = vec![vec![-1.0, 1.0, 0.0], vec![0.0, -1.0, 1.0]];
    let mut calls = 0;
    let expected = solve_gram_orthant(&rows, &[-1.0, -2.0], config(), || { calls += 1; false }).unwrap();
    for stop in 1..=calls {
        let mut at = 0;
        assert_eq!(solve_gram_orthant(&rows, &[-1.0, -2.0], config(), || { at += 1; at == stop }),
            Err(GramOrthantError::Cancelled));
    }
    assert_eq!(solve_gram_orthant(&rows, &[-1.0, -2.0], config(), || false).unwrap(), expected);
}
#[test]
fn dimensions_exact_setup_caps_and_nonconvergence_refuse() {
    let rows = vec![vec![-1.0, 1.0, 0.0], vec![0.0, -1.0, 1.0]];
    for cap in [11, 12] {
        assert_eq!(solve_gram_orthant(&rows, &[-1.0, -2.0],
            GramOrthantConfig { max_setup_products: cap, ..config() }, || false).is_ok(), cap == 12);
    }
    assert!(matches!(solve_gram_orthant(&rows, &[-1.0, -2.0],
        GramOrthantConfig { max_sweeps: 1, ..config() }, || false), Err(GramOrthantError::NotConverged { .. })));
    for bad in [vec![vec![]], vec![vec![0.0]], vec![vec![f64::NAN]], vec![vec![1.0], vec![1.0, 2.0]]] {
        assert!(solve_gram_orthant(&bad, &vec![-1.0; bad.len()], config(), || false).is_err());
    }
    assert!(solve_gram_orthant(&rows, &[1.0], config(), || false).is_err());
}
#[test]
fn positive_scaling_and_row_permutation_preserve_the_solution() {
    let a = vec![vec![-1.0, 1.0, 0.0], vec![0.0, -1.0, 1.0]];
    let x = solve_gram_orthant(&a, &[-1.0, -2.0], config(), || false).unwrap();
    let y = solve_gram_orthant(&[a[1].clone(), a[0].clone()], &[-2.0, -1.0], config(), || false).unwrap();
    near(x.point[0], y.point[1]); near(x.point[1], y.point[0]);
    let scaled: Vec<Vec<f64>> = a.iter().map(|r| r.iter().map(|v| 2.0*v).collect()).collect();
    let z = solve_gram_orthant(&scaled, &[-4.0, -8.0], config(), || false).unwrap();
    near(x.point[0], z.point[0]); near(x.point[1], z.point[1]);
}
