//! G1/G3: known-frequency harmonic amplitudes, phase, and sampled residuals.

use fs_fft::C64;
use fs_fft::harmonic::{HarmonicFitBudget, HarmonicFitError, HarmonicFitPlan};

fn budget() -> HarmonicFitBudget {
    HarmonicFitBudget {
        max_workspace_bytes: 8 * 1024 * 1024,
        min_reciprocal_condition: 1.0e-10,
    }
}

fn signal(
    n: usize,
    dt: f64,
    f0: f64,
    orders: &[usize],
    phasors: &[C64],
    dc: f64,
    t0: f64,
) -> Vec<f64> {
    (0..n)
        .map(|j| {
            let t = t0 + j as f64 * dt;
            dc + orders
                .iter()
                .zip(phasors)
                .map(|(&h, p)| {
                    let (s, c) = (core::f64::consts::TAU * h as f64 * f0 * t).sin_cos();
                    p.re * c + p.im * s
                })
                .sum::<f64>()
        })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:e}"
    );
}

#[test]
fn joint_fit_resolves_weak_even_harmonics_on_fractional_period_records() {
    let dt = 1.0 / 8192.0;
    let f0 = 432.125;
    let orders = [1, 2, 4, 7];
    let expected = [
        C64::new(1.2, -0.7),
        C64::new(2.0e-7, -5.0e-7),
        C64::new(-0.031, 0.044),
        C64::new(0.017, 0.009),
    ];
    for n in [127, 257, 513] {
        let plan = HarmonicFitPlan::new(n, dt, f0, &orders, budget()).unwrap();
        let x = signal(n, dt, f0, &orders, &expected, 4.75, 0.0);
        let fit = plan.fit(&x).unwrap();
        close(fit.dc, 4.75, 2.0e-12);
        for (actual, expected) in fit.phasors.iter().zip(expected) {
            close(actual.re, expected.re, 2.0e-12);
            close(actual.im, expected.im, 2.0e-12);
        }
        assert!(fit.residual_rms < 2.0e-12);
        assert!(plan.reciprocal_condition() > budget().min_reciprocal_condition);
    }
}

#[test]
fn fractional_period_fundamental_does_not_create_even_harmonics_or_bias_dc() {
    let plan = HarmonicFitPlan::new(61, 0.001, 17.3, &[1, 2, 3], budget()).unwrap();
    let x = signal(61, 0.001, 17.3, &[1], &[C64::new(3.1, -1.7)], 7.0, 0.0);
    // This record has a nonzero sampled sinusoid mean. Subtracting that mean
    // and fitting each nonorthogonal sinusoid independently would be biased.
    assert!((x.iter().sum::<f64>() / x.len() as f64 - 7.0).abs() > 0.01);
    let fit = plan.fit(&x).unwrap();
    close(fit.dc, 7.0, 1.0e-12);
    close(fit.phasors[0].re, 3.1, 1.0e-12);
    close(fit.phasors[0].im, -1.7, 1.0e-12);
    for p in &fit.phasors[1..] {
        assert!(p.re.abs() + p.im.abs() < 1.0e-12);
    }
}

#[test]
fn negative_time_phasors_follow_a_shift_of_the_record_origin() {
    let n = 203;
    let dt = 0.0004;
    let f0 = 63.7;
    let orders = [1, 3, 5];
    let phasors = [C64::new(0.7, -0.2), C64::new(-0.5, 0.4), C64::new(0.1, 0.3)];
    let t0 = 0.0131;
    let plan = HarmonicFitPlan::new(n, dt, f0, &orders, budget()).unwrap();
    let fit = plan
        .fit(&signal(n, dt, f0, &orders, &phasors, 0.2, t0))
        .unwrap();
    for ((actual, expected), h) in fit.phasors.iter().zip(phasors).zip(orders) {
        let (s, c) = (core::f64::consts::TAU * h as f64 * f0 * t0).sin_cos();
        close(actual.re, expected.re * c + expected.im * s, 2.0e-12);
        close(actual.im, expected.im * c - expected.re * s, 2.0e-12);
    }
}

#[test]
fn omitted_harmonic_remains_in_the_reported_sampled_residual() {
    let n = 256;
    let plan = HarmonicFitPlan::new(n, 1.0 / 1024.0, 8.0, &[1, 2], budget()).unwrap();
    let x = signal(
        n,
        1.0 / 1024.0,
        8.0,
        &[1, 2, 7],
        &[
            C64::new(1.0, 0.3),
            C64::new(-0.2, 0.1),
            C64::new(0.0, 0.031),
        ],
        0.6,
        0.0,
    );
    let fit = plan.fit(&x).unwrap();
    close(fit.residual_rms, 0.031 / 2.0_f64.sqrt(), 1.0e-12);
    close(fit.phasors[0].re, 1.0, 1.0e-12);
    close(fit.phasors[1].im, 0.1, 1.0e-12);
}

#[test]
fn reusable_scratch_replays_and_scales_across_physical_units() {
    let n = 93;
    let plan = HarmonicFitPlan::new(n, 0.001, 23.4, &[1, 4], budget()).unwrap();
    let x = signal(
        n,
        0.001,
        23.4,
        &[1, 4],
        &[C64::new(0.7, -0.2), C64::new(0.05, 0.03)],
        2.1,
        0.0,
    );
    let base = plan.fit(&x).unwrap();
    for scale in [1.0e-200, 1.0, 1.0e200] {
        let scaled: Vec<_> = x.iter().map(|v| v * scale).collect();
        let fit = plan.fit(&scaled).unwrap();
        close(fit.dc / scale, base.dc, 1.0e-12);
        let mut scratch = vec![f64::NAN; n];
        let mut phasors = vec![C64::default(); 2];
        let summary = plan.fit_into(&scaled, &mut scratch, &mut phasors).unwrap();
        assert_eq!(fit.phasors, phasors);
        assert_eq!(fit.dc.to_bits(), summary.dc.to_bits());
        assert_eq!(fit.residual_rms.to_bits(), summary.residual_rms.to_bits());
        for (p, b) in fit.phasors.iter().zip(&base.phasors) {
            close(p.re / scale, b.re, 1.0e-12);
            close(p.im / scale, b.im, 1.0e-12);
        }
    }
    let constant = plan.fit(&vec![f64::MAX; n]).unwrap();
    assert_eq!(constant.dc, f64::MAX);
    assert_eq!(constant.residual_rms, 0.0);
    assert!(constant.phasors.iter().all(|p| p.re == 0.0 && p.im == 0.0));
}

#[test]
fn large_preload_does_not_set_the_scale_of_small_harmonic_variation() {
    let plan = HarmonicFitPlan::new(205, 0.001, 37.1, &[1, 2], budget()).unwrap();
    let x = signal(
        205,
        0.001,
        37.1,
        &[1, 2],
        &[C64::new(1.0e-4, -3.0e-4), C64::new(-2.0e-5, 7.0e-5)],
        1.0e6,
        0.0,
    );
    let fit = plan.fit(&x).unwrap();
    close(fit.dc, 1.0e6, 2.0e-10);
    close(fit.phasors[0].re, 1.0e-4, 2.0e-11);
    close(fit.phasors[0].im, -3.0e-4, 2.0e-11);
    close(fit.phasors[1].re, -2.0e-5, 2.0e-11);
    close(fit.phasors[1].im, 7.0e-5, 2.0e-11);
}

#[test]
fn conditioning_nyquist_shape_and_memory_limits_refuse_before_fitting() {
    assert!(matches!(
        HarmonicFitPlan::new(32, 1.0e-8, 1.0, &[1, 2], budget()),
        Err(HarmonicFitError::IllConditioned { .. })
    ));
    for orders in [vec![], vec![0], vec![1, 1], vec![2, 1], vec![5]] {
        assert!(matches!(
            HarmonicFitPlan::new(32, 0.01, 10.0, &orders, budget()),
            Err(HarmonicFitError::InvalidParameter { .. })
        ));
    }
    assert!(matches!(
        HarmonicFitPlan::new(4, 0.01, 10.0, &[1, 2], budget()),
        Err(HarmonicFitError::InvalidParameter { .. })
    ));
    assert!(matches!(
        HarmonicFitPlan::new(32, f64::NAN, 10.0, &[1], budget()),
        Err(HarmonicFitError::InvalidParameter { .. })
    ));
    let mut small = budget();
    small.max_workspace_bytes = 1;
    assert!(matches!(
        HarmonicFitPlan::new(32, 0.01, 10.0, &[1], small),
        Err(HarmonicFitError::BudgetExceeded { .. })
    ));
    assert!(matches!(
        HarmonicFitPlan::new(usize::MAX, 0.01, 10.0, &[1], budget()),
        Err(HarmonicFitError::InvalidParameter { .. })
    ));
}

#[test]
fn refused_input_preserves_output_and_plan_remains_reusable() {
    let plan = HarmonicFitPlan::new(32, 0.001, 50.0, &[1], budget()).unwrap();
    let sentinel = C64::new(12.0, 34.0);
    let mut output = [sentinel];
    let mut scratch = vec![0.0; 32];
    let mut samples = vec![0.0; 32];
    samples[17] = f64::NAN;
    assert!(matches!(
        plan.fit_into(&samples, &mut scratch, &mut output),
        Err(HarmonicFitError::NonFinite { .. })
    ));
    assert_eq!(output, [sentinel]);
    assert!(
        plan.fit_into(&samples[..31], &mut scratch, &mut output)
            .is_err()
    );
    assert_eq!(output, [sentinel]);
    samples.fill(2.0);
    let result = plan.fit_into(&samples, &mut scratch, &mut output).unwrap();
    assert_eq!(result.dc, 2.0);
    assert_eq!(output, [C64::default()]);
}
