use super::*;
use crate::{Kernel, Matern};
use crate::noisy::{joint_normal_bank, q_noisy_expected_improvement};
use std::panic::{AssertUnwindSafe, catch_unwind};

fn model() -> Gp {
    Gp::fit(&[vec![10.0]], &[0.0], Kernel {
        family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.5],
    }, 0.1)
}

// A deterministic quadrature fixture: choose normal coordinates by solving
// L z = desired - mean. Expected utilities below are hand-computed, not
// recomputed by the production acquisition's feasibility/reduction path.
fn bank_for_draws(gps: &[&Gp], points: &[Vec<f64>], draws: &[Vec<Vec<f64>>]) -> Vec<f64> {
    let n = points.len();
    let m = gps.len();
    let posteriors: Vec<_> = gps.iter().map(|gp| gp.predict_joint(points)).collect();
    let mut bank = vec![0.0; draws.len() * n * m];
    for (s, draw) in draws.iter().enumerate() {
        for (output, (mean, lower)) in posteriors.iter().enumerate() {
            for i in 0..n {
                let mut rhs = draw[i][output] - mean[i];
                for j in 0..i { rhs -= lower[i * n + j] * bank[s * n * m + j * m + output]; }
                bank[s * n * m + i * m + output] = rhs / lower[i * n + i];
            }
        }
    }
    bank
}

#[test]
fn g0_infeasible_raw_minimum_is_not_an_incumbent() {
    let f = model(); let g = model();
    let points = vec![vec![0.0], vec![1.0]];
    let bank = bank_for_draws(&[&f, &g], &points, &[vec![vec![-100.0, 1.0], vec![3.0, -1.0]]]);
    let constraints = [ConstraintModel { gp: &g, upper_bound: 0.0 }];
    let gain = q_feasible_noisy_improvement(&f, &constraints, &points[..1], &points[1..], 10.0, &bank);
    assert!((gain - 7.0).abs() < 1e-12, "no feasible baseline uses reference, not -100: {gain}");
}

#[test]
fn g0_joint_draws_intersect_constraints_and_resample_baseline_feasibility() {
    let f = model(); let g = model(); let h = model();
    let points = vec![vec![0.0], vec![1.0], vec![2.0]];
    let draws = vec![
        // Baseline feasible: only the last candidate meets BOTH constraints.
        vec![vec![5.0, -1.0, -1.0], vec![-100.0, -1.0, 1.0], vec![3.0, -1.0, -1.0]],
        // No feasible baseline: first candidate gives utility 10 - 4.
        vec![vec![-100.0, 1.0, -1.0], vec![4.0, -1.0, -1.0], vec![2.0, 1.0, -1.0]],
        // No feasible point anywhere: zero.
        vec![vec![-100.0, 1.0, -1.0], vec![-20.0, -1.0, 1.0], vec![-30.0, 1.0, 1.0]],
        // Feasible candidate is worse than the reference: zero.
        vec![vec![-100.0, 1.0, -1.0], vec![12.0, -1.0, -1.0], vec![13.0, -1.0, -1.0]],
    ];
    let bank = bank_for_draws(&[&f, &g, &h], &points, &draws);
    let constraints = [ConstraintModel { gp: &g, upper_bound: 0.0 }, ConstraintModel { gp: &h, upper_bound: 0.0 }];
    let value = q_feasible_noisy_improvement(&f, &constraints, &points[..1], &points[1..], 10.0, &bank);
    assert!((value - 2.0).abs() < 1e-12, "(2 + 6 + 0 + 0)/4: {value}");
}

#[test]
fn g0_uncertain_feasibility_is_not_a_posterior_mean_filter() {
    let f = model(); let g = model();
    let points = vec![vec![0.0], vec![1.0]];
    let draws = vec![vec![vec![0.0, 1.0], vec![3.0, -1.0]],
        vec![vec![0.0, 1.0], vec![3.0, 1.0]]];
    let bank = bank_for_draws(&[&f, &g], &points, &draws);
    // The constraint GP has mean zero at both points, hence a mean-only filter
    // would mark both feasible. The specified draws give only one feasible gain.
    let constraints = [ConstraintModel { gp: &g, upper_bound: 0.1 }];
    assert_eq!(g.predict(&points[1]).0, 0.0);
    let value = q_feasible_noisy_improvement(&f, &constraints, &points[..1], &points[1..], 10.0, &bank);
    assert!((value - 3.5).abs() < 1e-12);
}

#[test]
fn g3_unconstrained_limit_matches_existing_noisy_acquisition_below_reference() {
    let f = model();
    let baseline = vec![vec![0.0]]; let candidates = vec![vec![1.0]];
    let bank = vec![-1.0, -2.0, 1.0, 0.0, 0.5, -0.5];
    let new = q_feasible_noisy_improvement(&f, &[], &baseline, &candidates, 10.0, &bank);
    let old = q_noisy_expected_improvement(&f, &baseline, &candidates, &bank);
    assert_eq!(new.to_bits(), old.to_bits());
    assert!(new > 0.1);
}

#[test]
fn g3_duplicate_coordinates_share_objective_and_constraint_latents() {
    let f = model(); let g = model();
    let constraints = [ConstraintModel { gp: &g, upper_bound: 0.0 }];
    assert_eq!(q_feasible_noisy_improvement(&f, &constraints, &[vec![0.0]], &[vec![-0.0]],
        10.0, &[100.0, 1.0, -100.0, -1.0]), 0.0);
    let bank = joint_normal_bank(128, 4, 77);
    let mut repeated = Vec::new();
    for row in bank.chunks_exact(4) { repeated.extend_from_slice(row); repeated.extend_from_slice(&[-1e10, -1e10]); }
    let single = q_feasible_noisy_improvement(&f, &constraints, &[vec![0.0]], &[vec![1.0]], 10.0, &bank);
    let double = q_feasible_noisy_improvement(&f, &constraints, &[vec![0.0]], &[vec![1.0], vec![1.0]],
        10.0, &repeated);
    assert_eq!(single.to_bits(), double.to_bits());
    assert!(single > 0.0);
}

#[test]
fn g0_upper_bounds_are_inclusive_and_reference_caps_utility() {
    let f = model(); let g = model();
    // Zero normal coordinates give exactly zero constraint draws, at the bound.
    let constraints = [ConstraintModel { gp: &g, upper_bound: 0.0 }];
    let bank = [0.0, 0.0, -1.0, 0.0];
    assert!(q_feasible_noisy_improvement(&f, &constraints, &[vec![0.0]], &[vec![1.0]], 10.0, &bank) > 0.0);
    assert_eq!(q_feasible_noisy_improvement(&f, &constraints, &[vec![0.0]], &[vec![1.0]], -10.0, &bank), 0.0);
}

#[test]
fn g3_nearby_candidates_keep_within_output_cross_covariance() {
    let f = model(); let g = model();
    let constraints = [ConstraintModel { gp: &g, upper_bound: 10.0 }];
    let bank = joint_normal_bank(4096, 4, 87);
    let value = q_feasible_noisy_improvement(&f, &constraints, &[vec![0.0]], &[vec![0.001]], 10.0, &bank);
    assert!(value < 0.01, "independent baseline/candidate draws would create false gain: {value}");
    assert_eq!(value.to_bits(), q_feasible_noisy_improvement(&f, &constraints, &[vec![0.0]],
        &[vec![0.001]], 10.0, &bank).to_bits());
}

#[test]
fn g0_invalid_constraint_banks_and_bounds_refuse() {
    let f = model(); let g = model();
    let b = vec![vec![0.0]]; let c = vec![vec![1.0]];
    let constraints = [ConstraintModel { gp: &g, upper_bound: 0.0 }];
    for bank in [vec![], vec![0.0; 3], vec![0.0, 0.0, 0.0, f64::NAN]] {
        assert!(catch_unwind(AssertUnwindSafe(|| q_feasible_noisy_improvement(&f, &constraints, &b, &c, 10.0, &bank))).is_err());
    }
    for point in [vec![], vec![0.0, 1.0], vec![f64::INFINITY]] {
        assert!(catch_unwind(AssertUnwindSafe(|| q_feasible_noisy_improvement(&f, &constraints, &b, &[point], 10.0, &[0.0; 4]))).is_err());
    }
    assert!(catch_unwind(AssertUnwindSafe(|| q_feasible_noisy_improvement(&f, &constraints, &b, &c, f64::NAN, &[0.0; 4]))).is_err());
    let bad = [ConstraintModel { gp: &g, upper_bound: f64::INFINITY }];
    assert!(catch_unwind(AssertUnwindSafe(|| q_feasible_noisy_improvement(&f, &bad, &b, &c, 10.0, &[0.0; 4]))).is_err());
}
