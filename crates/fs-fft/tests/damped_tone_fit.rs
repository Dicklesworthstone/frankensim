//! Physical pole envelopes, first-sample phasors, and joint double-decay fits.

use fs_fft::C64;
use fs_fft::harmonic::{
    DampedToneFitPlan, DampedTonePole, HarmonicFitBudget, HarmonicFitError, ToneFitPlan,
};

fn budget() -> HarmonicFitBudget {
    HarmonicFitBudget {
        max_workspace_bytes: 8 * 1024 * 1024,
        min_reciprocal_condition: 1.0e-10,
    }
}

fn pole(frequency_hz: f64, decay_rate_per_s: f64) -> DampedTonePole {
    DampedTonePole {
        frequency_hz,
        decay_rate_per_s,
    }
}

fn signal(
    n: usize,
    dt: f64,
    poles: &[DampedTonePole],
    phasors: &[C64],
    dc: f64,
    t0: f64,
) -> Vec<f64> {
    (0..n)
        .map(|j| {
            let t = t0 + j as f64 * dt;
            dc + poles
                .iter()
                .zip(phasors)
                .map(|(pole, p)| {
                    let (s, c) = (core::f64::consts::TAU * pole.frequency_hz * t).sin_cos();
                    (-pole.decay_rate_per_s * t).exp() * (p.re * c + p.im * s)
                })
                .sum::<f64>()
        })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:e}"
    );
}

#[test]
fn underdamped_oscillator_recovers_initial_state_quadratures() {
    // Independent solution of q'' + 2*zeta*wn*q' + wn^2*q = 0.
    // The oscillation frequency differs from the undamped natural frequency.
    let natural_omega = core::f64::consts::TAU * 123.4;
    let zeta = 0.013_f64;
    let decay = zeta * natural_omega;
    let damped_omega = natural_omega * (1.0 - zeta * zeta).sqrt();
    let (q0, v0) = (0.7, -30.0);
    let quadrature = (v0 + decay * q0) / damped_omega;
    let (n, dt) = (1009, 1.0 / 4096.0);
    let samples: Vec<_> = (0..n)
        .map(|j| {
            let t = j as f64 * dt;
            let (s, c) = (damped_omega * t).sin_cos();
            2.7 + (-decay * t).exp() * (q0 * c + quadrature * s)
        })
        .collect();
    let frequency = damped_omega / core::f64::consts::TAU;
    let fit = DampedToneFitPlan::new(n, dt, &[pole(frequency, decay)], budget())
        .unwrap()
        .fit(&samples)
        .unwrap();
    close(fit.dc, 2.7, 1e-12);
    close(fit.phasors[0].re, q0, 1e-12);
    close(fit.phasors[0].im, quadrature, 1e-12);
    close(
        damped_omega * fit.phasors[0].im - decay * fit.phasors[0].re,
        v0,
        1e-9,
    );
    assert!(fit.residual_rms < 1e-12);

    let stationary = ToneFitPlan::new(n, dt, &[frequency], budget())
        .unwrap()
        .fit(&samples)
        .unwrap();
    assert!(stationary.residual_rms > 0.05);
    assert!((stationary.phasors[0].re - q0).abs() > 0.2);
}

#[test]
fn same_frequency_components_retain_double_decay_and_caller_order() {
    let poles = [pole(243.17, 11.0), pole(243.17, 0.7), pole(883.1, 0.0)];
    let phasors = [
        C64::new(0.6, -0.4),
        C64::new(-0.25, 0.9),
        C64::new(1e-4, -3e-4),
    ];
    let (n, dt) = (2041, 1.0 / 4096.0);
    let samples = signal(n, dt, &poles, &phasors, 1.5, 0.0);
    let plan = DampedToneFitPlan::new(n, dt, &poles, budget()).unwrap();
    assert_eq!(plan.poles(), poles);
    assert_eq!(plan.sample_count(), n);
    assert_eq!(plan.timestep_s(), dt);
    let fit = plan.fit(&samples).unwrap();
    close(fit.dc, 1.5, 2e-12);
    for (actual, expected) in fit.phasors.iter().zip(phasors) {
        close(actual.re, expected.re, 2e-11);
        close(actual.im, expected.im, 2e-11);
    }
    assert!(fit.residual_rms < 2e-12);
}

#[test]
fn close_damped_tones_retain_beating_and_shifted_time_origin() {
    let poles = [pole(432.1, 1.2), pole(433.4, 4.1), pole(1739.7, 12.0)];
    let phasors = [
        C64::new(1.0, -0.4),
        C64::new(-0.7, 0.9),
        C64::new(0.01, 0.02),
    ];
    let (n, dt) = (1021, 1.0 / 8192.0);
    assert!(poles[1].frequency_hz - poles[0].frequency_hz < 1.0 / (n as f64 * dt));
    let plan = DampedToneFitPlan::new(n, dt, &poles, budget()).unwrap();
    for t0 in [0.0, 0.173] {
        let samples = signal(n, dt, &poles, &phasors, -3.0, t0);
        let fit = plan.fit(&samples).unwrap();
        close(fit.dc, -3.0, 2e-12);
        assert!(fit.residual_rms < 3e-12);
        for ((actual, p), pole) in fit.phasors.iter().zip(phasors).zip(poles) {
            let (s, c) = (core::f64::consts::TAU * pole.frequency_hz * t0).sin_cos();
            let attenuation = (-pole.decay_rate_per_s * t0).exp();
            close(actual.re, attenuation * (p.re * c + p.im * s), 2e-10);
            close(actual.im, attenuation * (p.im * c - p.re * s), 2e-10);
        }
    }
}

#[test]
fn zero_decay_uses_the_stationary_fit_bits_and_residual() {
    let frequencies = [37.19, 415.5, 909.27];
    let poles = frequencies.map(|f| pole(f, 0.0));
    let phasors = [
        C64::new(0.7, -0.3),
        C64::new(0.1, 0.25),
        C64::new(0.02, 0.04),
    ];
    let (n, dt) = (317, 1.0 / 8192.0);
    let mut samples = signal(n, dt, &poles, &phasors, 2.0, 0.0);
    for (j, x) in samples.iter_mut().enumerate() {
        *x += 0.03 * (j as f64 * 0.719).sin(); // intentionally unmodeled
    }
    let stationary = ToneFitPlan::new(n, dt, &frequencies, budget()).unwrap();
    let damped = DampedToneFitPlan::new(n, dt, &poles, budget()).unwrap();
    assert_eq!(
        stationary.reciprocal_condition().to_bits(),
        damped.reciprocal_condition().to_bits()
    );
    let expected = stationary.fit(&samples).unwrap();
    let actual = damped.fit(&samples).unwrap();
    assert_eq!(actual.dc.to_bits(), expected.dc.to_bits());
    assert_eq!(
        actual.residual_rms.to_bits(),
        expected.residual_rms.to_bits()
    );
    assert!(actual.residual_rms > 0.01);
    for (a, e) in actual.phasors.iter().zip(expected.phasors) {
        assert_eq!(a.re.to_bits(), e.re.to_bits());
        assert_eq!(a.im.to_bits(), e.im.to_bits());
    }
}

#[test]
fn shared_fit_preserves_physical_units_preloads_and_reusable_outputs() {
    let poles = [pole(37.1, 0.2), pole(73.4, 6.0)];
    let phasors = [C64::new(0.7, -0.2), C64::new(0.05, 0.03)];
    let (n, dt) = (411, 0.001);
    let plan = DampedToneFitPlan::new(n, dt, &poles, budget()).unwrap();
    let samples = signal(n, dt, &poles, &phasors, 2.1, 0.0);
    let mut scratch = vec![f64::NAN; n];
    let mut output = vec![C64::default(); poles.len()];
    for scale in [1e-200, 1.0, 1e200] {
        let scaled: Vec<_> = samples.iter().map(|x| x * scale).collect();
        let summary = plan.fit_into(&scaled, &mut scratch, &mut output).unwrap();
        let allocated = plan.fit(&scaled).unwrap();
        assert_eq!(output, allocated.phasors);
        assert_eq!(summary.dc.to_bits(), allocated.dc.to_bits());
        assert_eq!(
            summary.residual_rms.to_bits(),
            allocated.residual_rms.to_bits()
        );
        close(summary.dc / scale, 2.1, 2e-12);
        for (actual, expected) in output.iter().zip(phasors) {
            close(actual.re / scale, expected.re, 2e-12);
            close(actual.im / scale, expected.im, 2e-12);
        }
    }
    let weak = [C64::new(1e-4, -3e-4), C64::new(-2e-5, 7e-5)];
    let preload = plan.fit(&signal(n, dt, &poles, &weak, 1e6, 0.0)).unwrap();
    close(preload.dc, 1e6, 3e-10);
    for (actual, expected) in preload.phasors.iter().zip(weak) {
        close(actual.re, expected.re, 5e-11);
        close(actual.im, expected.im, 5e-11);
    }

    let saved = output.clone();
    let mut bad = samples;
    bad[17] = f64::NAN;
    assert!(matches!(
        plan.fit_into(&bad, &mut scratch, &mut output),
        Err(HarmonicFitError::NonFinite { .. })
    ));
    assert_eq!(output, saved);
}

#[test]
fn omitted_damped_components_remain_in_the_measured_residual() {
    let poles = [pole(83.7, 1.0), pole(191.2, 7.0), pole(327.8, 3.0)];
    let phasors = [
        C64::new(0.8, -0.4),
        C64::new(-0.2, 0.3),
        C64::new(0.03, 0.07),
    ];
    let (n, dt) = (1023, 0.001);
    let samples = signal(n, dt, &poles, &phasors, 0.4, 0.0);
    let fit = DampedToneFitPlan::new(n, dt, &poles[..2], budget())
        .unwrap()
        .fit(&samples)
        .unwrap();
    let reconstruction = signal(n, dt, &poles[..2], &fit.phasors, fit.dc, 0.0);
    let independently_measured_rms = (samples
        .iter()
        .zip(reconstruction)
        .map(|(&x, y)| (x - y).powi(2))
        .sum::<f64>()
        / n as f64)
        .sqrt();
    assert!(fit.residual_rms > 0.01);
    close(fit.residual_rms, independently_measured_rms, 2e-12);
}

#[test]
fn unresolved_duplicate_invalid_and_under_sampled_poles_refuse() {
    for poles in [
        vec![],
        vec![pole(100.0, 1.0), pole(100.0, 1.0)],
        vec![pole(0.0, 1.0)],
        vec![pole(500.0, 1.0)], // Nyquist
        vec![pole(100.0, -1.0)],
        vec![pole(f64::NAN, 1.0)],
        vec![pole(100.0, f64::INFINITY)],
    ] {
        assert!(matches!(
            DampedToneFitPlan::new(65, 0.001, &poles, budget()),
            Err(HarmonicFitError::InvalidParameter { .. })
        ));
    }
    assert!(matches!(
        DampedToneFitPlan::new(4, 0.001, &[pole(100.0, 1.0), pole(200.0, 2.0)], budget()),
        Err(HarmonicFitError::InvalidParameter { .. })
    ));
    for poles in [
        vec![pole(100.0, 1.0), pole(100.0, 1.0 + 1e-12)],
        vec![pole(100.0, 1e6)],      // sine dies before the next sample
        vec![pole(100.0, f64::MAX)], // finite declaration, unresolvable envelope
    ] {
        assert!(matches!(
            DampedToneFitPlan::new(65, 0.001, &poles, budget()),
            Err(HarmonicFitError::IllConditioned { .. })
        ));
    }
}

#[test]
fn admission_includes_both_physical_parameters_of_every_pole() {
    let frequencies = [37.19, 415.5];
    let poles = frequencies.map(|f| pole(f, 0.0));
    let (n, dt) = (65, 0.001);
    let stationary = ToneFitPlan::new(n, dt, &frequencies, budget()).unwrap();
    let required = DampedToneFitPlan::new(n, dt, &poles, budget())
        .unwrap()
        .workspace_bytes();
    assert_eq!(required, stationary.workspace_bytes() + poles.len() * 8);
    let exact = HarmonicFitBudget {
        max_workspace_bytes: required,
        ..budget()
    };
    DampedToneFitPlan::new(n, dt, &poles, exact).unwrap();
    assert!(matches!(
        DampedToneFitPlan::new(n, dt, &poles, HarmonicFitBudget {
            max_workspace_bytes: required - 1,
            ..exact
        }),
        Err(HarmonicFitError::BudgetExceeded { required_bytes, .. }) if required_bytes == required
    ));
}
