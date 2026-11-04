//! G1/G3: physical ringdown poles, double decay and shared-channel covariance.

use fs_fft::C64;
use fs_fft::harmonic::{
    DampedToneBounds, DampedToneFitPlan, DampedTonePole, DampedToneRefinementBudget,
    DampedToneRefinementError as Error, HarmonicFitBudget, HarmonicFitError, ToneFitPlan,
    refine_damped_tone_poles, refine_damped_tone_poles_channels,
};

fn budget() -> DampedToneRefinementBudget {
    DampedToneRefinementBudget {
        absolute_frequency_tolerance_hz: 1e-8,
        absolute_decay_tolerance_per_s: 1e-7,
        max_evaluations: 60,
        min_pole_reciprocal_condition: 1e-8,
        fit: HarmonicFitBudget {
            max_workspace_bytes: 8 * 1024 * 1024,
            min_reciprocal_condition: 1e-10,
        },
    }
}

fn pole(frequency_hz: f64, decay_rate_per_s: f64) -> DampedTonePole {
    DampedTonePole {
        frequency_hz,
        decay_rate_per_s,
    }
}

fn bound(frequency_hz: [f64; 2], decay_rate_per_s: [f64; 2]) -> DampedToneBounds {
    DampedToneBounds {
        frequency_hz,
        decay_rate_per_s,
    }
}

fn signal(n: usize, dt: f64, poles: &[DampedTonePole], amplitudes: &[C64], dc: f64) -> Vec<f64> {
    (0..n)
        .map(|j| {
            let t = j as f64 * dt;
            dc + poles
                .iter()
                .zip(amplitudes)
                .map(|(p, a)| {
                    let (s, c) = (core::f64::consts::TAU * p.frequency_hz * t).sin_cos();
                    (-p.decay_rate_per_s * t).exp() * (a.re * c + a.im * s)
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
fn oscillator_initial_value_solution_identifies_damped_frequency_and_amplitude_decay() {
    // Independent solution of q'' + 2*zeta*wn*q' + wn^2*q = 0. Neither
    // the undamped frequency nor the energy-decay rate is the fitted pole.
    let wn = core::f64::consts::TAU * 123.4;
    let zeta = 0.013_f64;
    let decay = zeta * wn;
    let wd = wn * (1.0 - zeta * zeta).sqrt();
    let (q0, v0) = (0.7, -30.0);
    let b = (v0 + decay * q0) / wd;
    let (n, dt) = (1009, 1.0 / 4096.0);
    let samples: Vec<_> = (0..n)
        .map(|j| {
            let t = j as f64 * dt;
            let (s, c) = (wd * t).sin_cos();
            2.7 + (-decay * t).exp() * (q0 * c + b * s)
        })
        .collect();
    let seeds = [pole(wd / core::f64::consts::TAU + 0.23, decay * 0.85)];
    let bounds = [bound([122.5, 124.5], [5.0, 15.0])];
    let result = refine_damped_tone_poles(&samples, dt, &seeds, &bounds, budget()).unwrap();
    close(result.poles[0].frequency_hz, wd / core::f64::consts::TAU, 1e-7);
    close(result.poles[0].decay_rate_per_s, decay, 1e-6);
    close(result.fit.phasors[0].re, q0, 1e-8);
    close(result.fit.phasors[0].im, b, 1e-8);
    close(
        wd * result.fit.phasors[0].im - decay * result.fit.phasors[0].re,
        v0,
        1e-6,
    );
    close(result.fit.dc, 2.7, 1e-10);
    assert!(result.fit.residual_rms < 1e-9);
    assert!(
        result.last_corrections[0].frequency_hz.abs()
            <= budget().absolute_frequency_tolerance_hz
    );
    assert!(
        result.last_corrections[0].decay_rate_per_s.abs()
            <= budget().absolute_decay_tolerance_per_s
    );
    assert_eq!(
        result,
        refine_damped_tone_poles(&samples, dt, &seeds, &bounds, budget()).unwrap()
    );
    let common =
        refine_damped_tone_poles_channels(&[&samples], dt, &seeds, &bounds, budget()).unwrap();
    assert_eq!(common.poles, result.poles);
    assert_eq!(common.fits[0], result.fit);
    assert_eq!(common.last_corrections, result.last_corrections);
}

#[test]
fn decaying_beating_unisons_and_a_weak_inharmonic_partial_are_refined_together() {
    let truth = [pole(233.17, 2.7), pole(235.04, 5.1), pole(731.238, 9.2)];
    let a = [
        C64::new(1.2, -0.7),
        C64::new(0.51, 0.83),
        C64::new(-0.014, 0.023),
    ];
    let seeds = [pole(233.05, 3.0), pole(235.19, 4.6), pole(731.418, 8.7)];
    let bounds = [
        bound([232.7, 233.7], [1.0, 4.0]),
        bound([234.6, 235.6], [3.0, 7.0]),
        bound([730.7, 731.8], [6.0, 12.0]),
    ];
    let dt = 1.0 / 8192.0;
    let samples = signal(3001, dt, &truth, &a, 2.6);
    assert!(truth[1].frequency_hz - truth[0].frequency_hz < 1.0 / (samples.len() as f64 * dt));
    let initial = DampedToneFitPlan::new(samples.len(), dt, &seeds, budget().fit)
        .unwrap()
        .fit(&samples)
        .unwrap();
    assert!(initial.residual_rms > 0.001);
    let result = refine_damped_tone_poles(&samples, dt, &seeds, &bounds, budget()).unwrap();
    for (k, expected) in truth.iter().enumerate() {
        close(result.poles[k].frequency_hz, expected.frequency_hz, 1e-7);
        close(result.poles[k].decay_rate_per_s, expected.decay_rate_per_s, 1e-6);
        close(result.fit.phasors[k].re, a[k].re, 1e-7);
        close(result.fit.phasors[k].im, a[k].im, 1e-7);
    }
    assert!(result.fit.residual_rms < 1e-9);
    assert!(result.evaluations < budget().max_evaluations);
    assert!(result.pole_reciprocal_condition >= budget().min_pole_reciprocal_condition);
}

#[test]
fn same_frequency_double_decay_keeps_the_original_pole_order() {
    let truth = [pole(418.3, 6.0), pole(418.3, 1.0)];
    let a = [C64::new(1.0, -0.3), C64::new(0.45, 0.8)];
    let seeds = [pole(418.36, 5.5), pole(418.25, 1.3)];
    let bounds = [
        bound([417.8, 418.8], [4.0, 8.0]),
        bound([417.8, 418.8], [0.3, 2.0]),
    ];
    let dt = 1.0 / 4096.0;
    let samples = signal(2001, dt, &truth, &a, 0.3);
    let result = refine_damped_tone_poles(&samples, dt, &seeds, &bounds, budget()).unwrap();
    for k in 0..2 {
        close(result.poles[k].frequency_hz, truth[k].frequency_hz, 1e-7);
        close(result.poles[k].decay_rate_per_s, truth[k].decay_rate_per_s, 1e-6);
        close(result.fit.phasors[k].re, a[k].re, 1e-7);
        close(result.fit.phasors[k].im, a[k].im, 1e-7);
    }
    let reversed = refine_damped_tone_poles(
        &samples,
        dt,
        &[seeds[1], seeds[0]],
        &[bounds[1], bounds[0]],
        budget(),
    )
    .unwrap();
    close(reversed.poles[0].decay_rate_per_s, result.poles[1].decay_rate_per_s, 1e-6);
    close(reversed.fit.phasors[0].re, result.fit.phasors[1].re, 1e-7);
}

#[test]
fn common_channels_preserve_units_clock_origin_and_vector_rotation() {
    let truth = [pole(147.2, 1.3), pole(148.9, 3.6)];
    let seeds = [pole(147.08, 1.1), pole(149.03, 3.9)];
    let bounds = [
        bound([146.7, 147.7], [0.5, 2.0]),
        bound([148.4, 149.4], [2.0, 5.0]),
    ];
    let shift = |a: C64, p: DampedTonePole| {
        let t0 = 0.01736;
        let (s, c) = (core::f64::consts::TAU * p.frequency_hz * t0).sin_cos();
        let e = (-p.decay_rate_per_s * t0).exp();
        C64::new(e * (a.re * c + a.im * s), e * (-a.re * s + a.im * c))
    };
    let a = [shift(C64::new(1.3, -0.4), truth[0]), C64::default()];
    let b = [
        shift(C64::new(0.1, 0.2), truth[0]),
        shift(C64::new(-1.4, 0.6), truth[1]),
    ];
    let dt = 1.0 / 4096.0;
    let x = signal(1601, dt, &truth, &a, 0.3);
    let y = signal(1601, dt, &truth, &b, -1.2);
    let z = vec![7.0; x.len()];
    for scale in [1e-150, 1.0, 1e150] {
        let sx: Vec<_> = x.iter().map(|x| x * scale).collect();
        let sy: Vec<_> = y.iter().map(|x| x * scale).collect();
        let sz: Vec<_> = z.iter().map(|x| x * scale).collect();
        let result =
            refine_damped_tone_poles_channels(&[&sx, &sy, &sz], dt, &seeds, &bounds, budget())
                .unwrap();
        for k in 0..2 {
            close(result.poles[k].frequency_hz, truth[k].frequency_hz, 1e-7);
            close(result.poles[k].decay_rate_per_s, truth[k].decay_rate_per_s, 1e-6);
            close(result.fits[0].phasors[k].re / scale, a[k].re, 1e-8);
            close(result.fits[1].phasors[k].im / scale, b[k].im, 1e-8);
            assert_eq!(result.fits[2].phasors[k], C64::default());
        }
        close(result.fits[0].dc / scale, 0.3, 1e-9);
        close(result.fits[2].dc / scale, 7.0, 1e-14);
    }
    let (s, c) = 0.67_f64.sin_cos();
    let rx: Vec<_> = x.iter().zip(&y).map(|(&x, &y)| c * x - s * y).collect();
    let ry: Vec<_> = x.iter().zip(&y).map(|(&x, &y)| s * x + c * y).collect();
    let rotated =
        refine_damped_tone_poles_channels(&[&rx, &ry, &z], dt, &seeds, &bounds, budget()).unwrap();
    for k in 0..2 {
        close(rotated.poles[k].frequency_hz, truth[k].frequency_hz, 1e-7);
        close(rotated.poles[k].decay_rate_per_s, truth[k].decay_rate_per_s, 1e-6);
        close(rotated.fits[0].phasors[k].re, c * a[k].re - s * b[k].re, 1e-8);
        close(rotated.fits[1].phasors[k].im, s * a[k].im + c * b[k].im, 1e-8);
    }
}

#[test]
fn unmodeled_motion_remains_in_the_residual_at_the_local_stationary_point() {
    let truth = [pole(181.37, 2.2), pole(376.81, 5.6), pole(713.17, 3.0)];
    let a = [
        C64::new(0.8, -0.6),
        C64::new(0.24, 0.13),
        C64::new(0.0, 0.01),
    ];
    let dt = 1.0 / 8000.0;
    let samples = signal(4093, dt, &truth, &a, 3.0);
    let bounds = [
        bound([180.9, 181.9], [1.0, 4.0]),
        bound([376.3, 377.3], [4.0, 8.0]),
    ];
    let mut limits = budget();
    limits.absolute_frequency_tolerance_hz = 1e-6;
    limits.absolute_decay_tolerance_per_s = 1e-5;
    let result = refine_damped_tone_poles(
        &samples,
        dt,
        &[pole(181.6, 2.0), pole(376.7, 6.0)],
        &bounds,
        limits,
    )
    .unwrap();
    assert!(result.fit.residual_rms > 0.001);
    let reconstructed =
        signal(samples.len(), dt, &result.poles, &result.fit.phasors, result.fit.dc);
    let rms = (samples
        .iter()
        .zip(&reconstructed)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        / samples.len() as f64)
        .sqrt();
    close(rms, result.fit.residual_rms, 1e-12);
    for k in 0..2 {
        for (df, dd) in [(1e-4, 0.0), (-1e-4, 0.0), (0.0, 1e-3), (0.0, -1e-3)] {
            let mut perturbed = result.poles.clone();
            perturbed[k].frequency_hz += df;
            perturbed[k].decay_rate_per_s += dd;
            let fit = DampedToneFitPlan::new(samples.len(), dt, &perturbed, limits.fit)
                .unwrap()
                .fit(&samples)
                .unwrap();
            assert!(fit.residual_rms > result.fit.residual_rms);
        }
    }
}

#[test]
fn zero_decay_retains_the_stationary_fit_bits() {
    let truth = [pole(180.3, 0.0), pole(371.2, 0.0)];
    let samples = signal(
        1023,
        0.001,
        &truth,
        &[C64::new(1.2, -0.7), C64::new(0.2, 0.1)],
        3.0,
    );
    let bounds = [
        bound([180.0, 180.6], [0.0, 1.0]),
        bound([370.8, 371.6], [0.0, 1.0]),
    ];
    let mut limits = budget();
    limits.max_evaluations = 1;
    let result = refine_damped_tone_poles(&samples, 0.001, &truth, &bounds, limits).unwrap();
    let stationary = ToneFitPlan::new(samples.len(), 0.001, &[180.3, 371.2], limits.fit)
        .unwrap()
        .fit(&samples)
        .unwrap();
    assert_eq!(result.fit, stationary);
    assert_eq!(result.poles, truth);
}

#[test]
fn absent_unresolved_and_constant_poles_are_not_identified() {
    let dt = 1.0 / 8000.0;
    let samples = signal(1023, dt, &[pole(180.3, 2.0)], &[C64::new(1.2, -0.7)], 3.0);
    let seeds = [pole(180.3, 2.0), pole(371.2, 4.0)];
    let bounds = [
        bound([180.0, 180.6], [1.0, 3.0]),
        bound([370.8, 371.6], [3.0, 5.0]),
    ];
    assert!(matches!(
        refine_damped_tone_poles(&samples, dt, &seeds, &bounds, budget()),
        Err(Error::Unidentifiable { .. })
    ));
    assert!(matches!(
        refine_damped_tone_poles(&[3.0; 31], dt, &seeds, &bounds, budget()),
        Err(Error::Unidentifiable {
            reciprocal_condition: 0.0,
            ..
        })
    ));
    let seeds = [pole(180.3, 1.0 + 1e-10), pole(180.3, 1.0 + 1e-9)];
    let bounds = [
        bound([180.0, 180.6], [1.0, 1.0 + 2e-10]),
        bound([180.0, 180.6], [1.0 + 8e-10, 1.0 + 1.2e-9]),
    ];
    assert!(matches!(
        refine_damped_tone_poles(&samples, dt, &seeds, &bounds, budget()),
        Err(Error::Fit(HarmonicFitError::IllConditioned { .. }))
    ));
}

#[test]
fn residual_curvature_does_not_identify_an_unexcited_declared_pole() {
    // Complete cycles make the 400 Hz amplitude exactly absent in real
    // arithmetic. The omitted 404 Hz tone gives D^T*r a nonzero contribution
    // at that pole, so full-Jacobian rank alone would invent observability.
    let dt = 1.0 / 4096.0;
    let samples = signal(
        1024,
        dt,
        &[pole(200.0, 0.0), pole(404.0, 0.0)],
        &[C64::new(1.0, 0.2), C64::new(0.3, -0.1)],
        0.5,
    );
    let seeds = [pole(200.0, 0.0), pole(400.0, 0.0)];
    let bounds = [
        bound([199.0, 201.0], [0.0, 1.0]),
        bound([399.0, 401.0], [0.0, 1.0]),
    ];
    let fixed = DampedToneFitPlan::new(samples.len(), dt, &seeds, budget().fit)
        .unwrap()
        .fit(&samples)
        .unwrap();
    assert!(fixed.phasors[1].re.abs() < 1e-13 && fixed.phasors[1].im.abs() < 1e-13);
    assert!(fixed.residual_rms > 0.1);
    assert!(matches!(
        refine_damped_tone_poles(&samples, dt, &seeds, &bounds, budget()),
        Err(Error::Unidentifiable { .. })
    ));
}

#[test]
fn a_nonconstant_omitted_signal_cannot_identify_a_roundoff_amplitude_pole() {
    // With only one declared pole, its uniformly tiny projected Jacobian can
    // be well-conditioned. A relative rank gate alone is insufficient.
    let dt = 1.0 / 4096.0;
    let seeds = [pole(400.0, 0.0)];
    let bounds = [bound([399.0, 401.0], [0.0, 1.0])];
    let omitted = signal(
        1024,
        dt,
        &[pole(404.0, 0.0)],
        &[C64::new(0.3, -0.1)],
        0.5,
    );
    let resolved = signal(1024, dt, &seeds, &[C64::new(0.3, -0.1)], 0.5);
    for scale in [1e-150, 1.0, 1e150] {
        let samples: Vec<_> = omitted.iter().map(|x| x * scale).collect();
        let fixed = DampedToneFitPlan::new(samples.len(), dt, &seeds, budget().fit)
            .unwrap()
            .fit(&samples)
            .unwrap();
        assert!(fixed.phasors[0].re.abs() / scale < 1e-13);
        assert!(fixed.phasors[0].im.abs() / scale < 1e-13);
        assert!(fixed.residual_rms / scale > 0.1);
        assert!(matches!(
            refine_damped_tone_poles(&samples, dt, &seeds, &bounds, budget()),
            Err(Error::Unidentifiable {
                reciprocal_condition: 0.0,
                ..
            })
        ));

        // The same physical amplitude, now at the declared pole, is resolved
        // regardless of its units. The floor is not an absolute input cutoff.
        let samples: Vec<_> = resolved.iter().map(|x| x * scale).collect();
        let result = refine_damped_tone_poles(&samples, dt, &seeds, &bounds, budget()).unwrap();
        assert_eq!(result.poles, seeds);
        close(result.fit.phasors[0].re / scale, 0.3, 1e-12);
        close(result.fit.phasors[0].im / scale, -0.1, 1e-12);
    }
}

#[test]
fn storage_work_bounds_and_input_admission_do_not_publish_partial_fits() {
    let dt = 1.0 / 8000.0;
    let truth = [pole(200.0, 4.0)];
    let samples = signal(127, dt, &truth, &[C64::new(1.2, -0.7)], 3.0);
    let bounds = [bound([199.5, 200.5], [2.0, 6.0])];
    let mut limits = budget();
    limits.fit.max_workspace_bytes = 1;
    let Err(Error::Fit(HarmonicFitError::BudgetExceeded { required_bytes, .. })) =
        refine_damped_tone_poles(&samples, dt, &truth, &bounds, limits)
    else {
        panic!("complete refinement storage was not admitted")
    };
    limits.fit.max_workspace_bytes = required_bytes;
    refine_damped_tone_poles(&samples, dt, &truth, &bounds, limits).unwrap();
    assert!(matches!(
        refine_damped_tone_poles_channels(&[&samples, &samples], dt, &truth, &bounds, limits),
        Err(Error::Fit(HarmonicFitError::BudgetExceeded { .. }))
    ));
    limits.fit.max_workspace_bytes -= 1;
    assert!(matches!(
        refine_damped_tone_poles(&samples, dt, &truth, &bounds, limits),
        Err(Error::Fit(HarmonicFitError::BudgetExceeded { .. }))
    ));
    let mut limits = budget();
    limits.max_evaluations = 1;
    assert_eq!(
        refine_damped_tone_poles(&samples, dt, &truth, &bounds, limits)
            .unwrap()
            .evaluations,
        1
    );
    assert!(matches!(
        refine_damped_tone_poles(&samples, dt, &[pole(200.2, 4.2)], &bounds, limits),
        Err(Error::NotConverged { evaluations: 1, .. })
    ));
    for b in [
        bound([200.1, 200.4], [2.0, 6.0]),
        bound([199.5, 200.5], [4.1, 4.4]),
    ] {
        assert!(matches!(
            refine_damped_tone_poles(&samples, dt, &[pole(200.2, 4.2)], &[b], budget()),
            Err(Error::NotConverged { .. })
        ));
    }
    for b in [
        bound([199.0, 4000.0], [2.0, 6.0]),
        bound([199.0, 201.0], [-1.0, 6.0]),
        bound([199.0, 201.0], [5.0, 6.0]),
    ] {
        assert!(matches!(
            refine_damped_tone_poles(&samples, dt, &truth, &[b], budget()),
            Err(Error::InvalidParameter { .. })
        ));
    }
    assert!(matches!(
        refine_damped_tone_poles(
            &samples,
            dt,
            &[truth[0], pole(200.1, 4.1)],
            &[bounds[0], bounds[0]],
            budget()
        ),
        Err(Error::InvalidParameter { .. })
    ));
    for bad_dt in [0.0, -1.0, f64::NAN] {
        assert!(matches!(
            refine_damped_tone_poles(&samples, bad_dt, &truth, &bounds, budget()),
            Err(Error::InvalidParameter { .. })
        ));
    }
    assert!(matches!(
        refine_damped_tone_poles_channels(
            &[&samples, &samples[1..]],
            dt,
            &truth,
            &bounds,
            budget()
        ),
        Err(Error::InvalidParameter { .. })
    ));
    assert!(matches!(
        refine_damped_tone_poles(&samples[..4], dt, &truth, &bounds, budget()),
        Err(Error::InvalidParameter { .. })
    ));
    let mut poisoned = samples.clone();
    poisoned[7] = f64::NAN;
    assert!(matches!(
        refine_damped_tone_poles(&poisoned, dt, &truth, &bounds, budget()),
        Err(Error::Fit(HarmonicFitError::NonFinite { .. }))
    ));
}
