use fs_ascent::quadratic_ball::{
    QuadraticBallBranch, QuadraticBallError, QuadraticBallMaximum, QuadraticBallOptions,
    maximize_quadratic_ball,
};
use std::ops::ControlFlow;

fn solve(g: &[f64], h: &[f64]) -> QuadraticBallMaximum {
    maximize_quadratic_ball(g, h, Default::default(), || ControlFlow::Continue(())).unwrap()
}
fn value(g: &[f64], h: &[f64], a: &[f64]) -> f64 {
    let n = h.len();
    (0..n)
        .map(|i| a[i] * (0..n).map(|j| g[i * n + j] * a[j]).sum::<f64>() + 2.0 * h[i] * a[i])
        .sum()
}
fn evidence(report: &QuadraticBallMaximum) {
    assert!(report.coefficients.iter().map(|a| a * a).sum::<f64>() <= 1.0);
    assert!(report.stationarity_residual < 1e-10);
    assert!(report.feasibility_residual < 1e-10);
    assert!(report.complementarity_residual < 1e-10);
    assert!(report.relative_duality_gap < 1e-10);
    assert!(report.eigenpair_residual < 1e-10);
    assert!(report.orthogonality_defect < 1e-10);
}

#[test]
fn analytic_one_dimensional_regular_and_extreme_uniform_scales() {
    for scale in [1e-320, 1e-250, 1.0, 1e250] {
        let g = [3.0 * scale];
        let h = [-0.7 * scale];
        let report = solve(&g, &h);
        assert_eq!(report.branch, QuadraticBallBranch::Regular);
        assert!(report.gradient_admissible());
        assert_eq!(report.coefficients, [-1.0]);
        let expected = g[0] + 2.0 * h[0].abs();
        let allowance = (1e-13 * expected.abs()).max(4.0 * f64::from_bits(1));
        assert!((report.value - expected).abs() <= allowance);
        assert!((report.multiplier - (g[0] + h[0].abs())).abs() <= allowance);
        evidence(&report);
    }
}

#[test]
fn hard_case_and_repeated_top_eigenspace_have_deterministic_witnesses() {
    let g = [4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.0];
    let h = [0.0, 0.0, 1.0];
    let report = solve(&g, &h);
    assert_eq!(report.branch, QuadraticBallBranch::HardCase);
    assert!(!report.gradient_admissible());
    assert!((report.value - 13.0 / 3.0).abs() < 1e-12);
    assert_eq!(report.multiplier, 4.0);
    assert!(report.coefficients[0] > 0.0);
    assert_eq!(report.coefficients[1], 0.0);
    assert!((report.coefficients[2] - 1.0 / 3.0).abs() < 1e-14);
    assert_eq!(report, solve(&g, &h));
    evidence(&report);

    let zero_linear = solve(&[1.0, 0.0, 0.0, 4.0], &[0.0, 0.0]);
    assert_eq!(zero_linear.coefficients, [0.0, 1.0]);
    assert_eq!(zero_linear.value, 4.0);
    assert!(!zero_linear.gradient_admissible());
    let flat = solve(&[0.0; 4], &[0.0; 2]);
    assert_eq!(flat.branch, QuadraticBallBranch::Flat);
    assert_eq!(flat.coefficients, [1.0, 0.0]);
    assert_eq!(flat.value, 0.0);
    assert!(!flat.gradient_admissible());
}

#[test]
fn singular_psd_and_linear_only_inputs_use_the_same_global_owner() {
    let rank_one = solve(&[1.0; 4], &[0.3, 0.3]);
    assert!((rank_one.value - (2.0 + 0.6 * fs_math::det::sqrt(2.0))).abs() < 1e-11);
    assert!(rank_one.gradient_admissible());
    evidence(&rank_one);
    let linear = solve(&[0.0; 4], &[3.0, -4.0]);
    assert!((linear.value - 10.0).abs() < 1e-12);
    assert!((linear.coefficients[0] - 0.6).abs() < 1e-14);
    assert!((linear.coefficients[1] + 0.8).abs() < 1e-14);
    assert!(linear.gradient_admissible());
    evidence(&linear);
}

#[test]
fn near_hard_shift_resolves_the_witness_without_claiming_a_smooth_gradient() {
    // A tiny nonzero top-space forcing must not be rounded into an arbitrary
    // hard-case sign. Geometric bracketing spans the full exponent range.
    let report = solve(&[4.0, 0.0, 0.0, 1.0], &[-1e-300, 1.0]);
    assert_eq!(report.branch, QuadraticBallBranch::Regular);
    assert!(report.coefficients[0] < 0.0);
    assert!((report.coefficients[1] - 1.0 / 3.0).abs() < 1e-11);
    assert!(!report.gradient_admissible());
    assert!(report.iterations < QuadraticBallOptions::default().max_iterations);
    evidence(&report);
}

#[test]
fn rotated_regular_solution_dominates_circle_and_has_envelope_derivatives() {
    let g = [3.0, 1.0, 1.0, 2.0];
    let h = [0.35, -0.8];
    let report = solve(&g, &h);
    assert!(report.gradient_admissible());
    evidence(&report);
    assert!((report.value - value(&g, &h, &report.coefficients)).abs() < 1e-12);
    for i in 0..4096 {
        let angle = std::f64::consts::TAU * (i as f64) / 4096.0;
        let a = [fs_math::det::cos(angle), fs_math::det::sin(angle)];
        assert!(value(&g, &h, &a) <= report.value + 1e-10);
    }
    let step = 1e-5;
    for i in 0..2 {
        let (mut plus, mut minus) = (h, h);
        plus[i] += step;
        minus[i] -= step;
        let actual = (solve(&g, &plus).value - solve(&g, &minus).value) / (2.0 * step);
        assert!((actual - 2.0 * report.coefficients[i]).abs() < 1e-7);
        let (mut plus, mut minus) = (g, g);
        plus[i * 2 + i] += step;
        minus[i * 2 + i] -= step;
        let actual = (solve(&plus, &h).value - solve(&minus, &h).value) / (2.0 * step);
        assert!((actual - report.coefficients[i] * report.coefficients[i]).abs() < 1e-7);
    }
    let (mut plus, mut minus) = (g, g);
    for index in [1, 2] {
        plus[index] += step;
        minus[index] -= step;
    }
    let actual = (solve(&plus, &h).value - solve(&minus, &h).value) / (2.0 * step);
    assert!((actual - 2.0 * report.coefficients[0] * report.coefficients[1]).abs() < 1e-7);
}

#[test]
fn admission_does_not_hide_indefiniteness_behind_a_large_linear_term() {
    let result = maximize_quadratic_ball(
        &[1.0, 0.0, 0.0, -0.001],
        &[1e200, 0.0],
        Default::default(),
        || ControlFlow::Continue(()),
    );
    assert!(matches!(result, Err(QuadraticBallError::Spectrum(_))));
    assert!(matches!(
        maximize_quadratic_ball(&[1.0, 0.1, 0.0, 1.0], &[0.0; 2], Default::default(), || {
            ControlFlow::Continue(())
        }),
        Err(QuadraticBallError::Spectrum(_))
    ));
    // Numerical PSD admission retains this small negative value; it does not
    // silently replace the diagonal or the eigenvalue with zero.
    let tiny_negative = solve(&[1.0, 0.0, 0.0, -1e-12], &[0.2, 0.4]);
    assert_eq!(tiny_negative.minimum_eigenvalue, -1e-12);
    assert!(
        (tiny_negative.value
            - value(
                &[1.0, 0.0, 0.0, -1e-12],
                &[0.2, 0.4],
                &tiny_negative.coefficients
            ))
        .abs()
            < 1e-12
    );
    evidence(&tiny_negative);
}

#[test]
fn bounds_cancellation_and_unrepresentable_results_return_no_witness() {
    assert!(matches!(
        maximize_quadratic_ball(&[], &[], Default::default(), || ControlFlow::Continue(())),
        Err(QuadraticBallError::Invalid(_))
    ));
    assert!(matches!(
        maximize_quadratic_ball(&[0.0; 289], &[0.0; 17], Default::default(), || {
            ControlFlow::Continue(())
        }),
        Err(QuadraticBallError::Invalid(_))
    ));
    assert!(matches!(
        maximize_quadratic_ball(&[f64::NAN], &[1.0], Default::default(), || {
            ControlFlow::Continue(())
        }),
        Err(QuadraticBallError::Invalid(_))
    ));
    assert!(matches!(
        maximize_quadratic_ball(&[f64::MAX], &[f64::MAX], Default::default(), || {
            ControlFlow::Continue(())
        }),
        Err(QuadraticBallError::Unrepresentable(_))
    ));
    assert!(matches!(
        maximize_quadratic_ball(
            &[f64::MAX],
            &[f64::from_bits(1)],
            Default::default(),
            || ControlFlow::Continue(())
        ),
        Err(QuadraticBallError::Unrepresentable(_))
    ));
    let mut calls = 0;
    let result = maximize_quadratic_ball(
        &[4.0, 0.0, 0.0, 1.0],
        &[0.6, 1.0],
        Default::default(),
        || {
            calls += 1;
            if calls >= 5 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        },
    );
    assert_eq!(result, Err(QuadraticBallError::Cancelled));
    let options = QuadraticBallOptions {
        max_iterations: 1,
        ..Default::default()
    };
    assert_eq!(
        maximize_quadratic_ball(&[4.0, 0.0, 0.0, 1.0], &[0.6, 1.0], options, || {
            ControlFlow::Continue(())
        }),
        Err(QuadraticBallError::NoConvergence { iterations: 1 })
    );
}
