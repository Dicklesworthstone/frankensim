//! Frequency drift must remain a physical phase law, not spectral broadening.

use fs_fft::C64;
use fs_fft::harmonic::{
    ChirpedTone as Tone, ChirpedToneBounds as Bounds, ChirpedToneFitPlan as Plan,
    ChirpedToneRefinementBudget as Budget, ChirpedToneRefinementError as Error, HarmonicFitBudget,
    HarmonicFitError, ToneFitPlan, refine_chirped_tones, refine_chirped_tones_channels,
};

fn budget() -> Budget {
    Budget {
        absolute_frequency_tolerance_hz: 1e-8,
        max_evaluations: 32,
        min_chirp_reciprocal_condition: 1e-8,
        fit: HarmonicFitBudget {
            max_workspace_bytes: 8 * 1024 * 1024,
            min_reciprocal_condition: 1e-11,
        },
    }
}
fn tone(start: f64, end: f64) -> Tone {
    Tone {
        start_frequency_hz: start,
        end_frequency_hz: end,
    }
}
fn bounds(tones: &[Tone]) -> Vec<Bounds> {
    tones
        .iter()
        .map(|p| Bounds {
            start_frequency_hz: [p.start_frequency_hz - 1.0, p.start_frequency_hz + 1.0],
            end_frequency_hz: [p.end_frequency_hz - 1.0, p.end_frequency_hz + 1.0],
        })
        .collect()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} (tolerance {tolerance:e})"
    );
}
fn close_tones(actual: &[Tone], expected: &[Tone], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (a, b) in actual.iter().zip(expected) {
        close(a.start_frequency_hz, b.start_frequency_hz, tolerance);
        close(a.end_frequency_hz, b.end_frequency_hz, tolerance);
    }
}

// Independent quadratic-phase law in seconds. Production instead integrates
// endpoint frequencies in DD and reduces cycles before evaluating trig.
fn signal(n: usize, dt: f64, tones: &[Tone], phasors: &[C64], dc: f64) -> Vec<f64> {
    let duration = (n - 1) as f64 * dt;
    (0..n)
        .map(|j| {
            let t = j as f64 * dt;
            dc + tones
                .iter()
                .zip(phasors)
                .map(|(p, a)| {
                    let rate = (p.end_frequency_hz - p.start_frequency_hz) / duration;
                    let phase =
                        core::f64::consts::TAU * (p.start_frequency_hz * t + 0.5 * rate * t * t);
                    let (s, c) = phase.sin_cos();
                    a.re * c + a.im * s
                })
                .sum::<f64>()
        })
        .collect()
}

#[test]
fn supplied_crossing_chirps_recover_phase_dc_and_unweakened_amplitudes() {
    let (n, dt) = (501, 1.0 / 2048.0);
    let truth = [tone(137.2, 159.4), tone(161.7, 132.6)];
    let phasors = [C64::new(0.8, -0.4), C64::new(-0.2, 0.65)];
    let samples = signal(n, dt, &truth, &phasors, 2.7);
    let plan = Plan::new(n, dt, &truth, budget().fit).unwrap();
    assert_eq!(plan.tones(), truth);
    assert_eq!(plan.sample_count(), n);
    assert_eq!(plan.timestep_s(), dt);
    let fit = plan.fit(&samples).unwrap();
    close(fit.dc, 2.7, 2e-12);
    for (a, b) in fit.phasors.iter().zip(phasors) {
        close(a.re, b.re, 2e-12);
        close(a.im, b.im, 2e-12);
    }
    assert!(fit.residual_rms < 2e-12);
    let stationary = ToneFitPlan::new(n, dt, &[147.15, 148.3], budget().fit)
        .unwrap()
        .fit(&samples)
        .unwrap();
    assert!(
        stationary.residual_rms > 0.1,
        "the stationary approximation hides real drift"
    );
    let mut scratch = vec![0.0; n];
    let mut output = [C64::default(); 2];
    let summary = plan.fit_into(&samples, &mut scratch, &mut output).unwrap();
    assert_eq!(output.as_slice(), fit.phasors.as_slice());
    assert_eq!(summary.residual_rms.to_bits(), fit.residual_rms.to_bits());
}

#[test]
fn joint_refinement_recovers_independently_crossing_frequency_trajectories() {
    let (n, dt) = (501, 1.0 / 2048.0);
    let truth = [tone(137.2, 159.4), tone(161.7, 132.6)];
    let seeds = [tone(137.57, 159.16), tone(161.41, 132.92)];
    let phasors = [C64::new(0.8, -0.4), C64::new(-0.2, 0.65)];
    let samples = signal(n, dt, &truth, &phasors, 0.7);
    let result = refine_chirped_tones(&samples, dt, &seeds, &bounds(&truth), budget()).unwrap();
    close_tones(&result.tones, &truth, 2e-8);
    close(result.fit.dc, 0.7, 1e-10);
    for (a, b) in result.fit.phasors.iter().zip(phasors) {
        close(a.re, b.re, 1e-9);
        close(a.im, b.im, 1e-9);
    }
    assert!(result.fit.residual_rms < 1e-10);
    assert!(result.evaluations <= 12);
    assert!(result.chirp_reciprocal_condition > budget().min_chirp_reciprocal_condition);
    assert!(result.last_corrections.iter().all(|d| {
        d.start_frequency_hz.abs() <= budget().absolute_frequency_tolerance_hz
            && d.end_frequency_hz.abs() <= budget().absolute_frequency_tolerance_hz
    }));
    assert_eq!(
        result,
        refine_chirped_tones(&samples, dt, &seeds, &bounds(&truth), budget()).unwrap()
    );
}

#[test]
fn simultaneous_sensors_retain_nodal_tones_frame_changes_and_physical_units() {
    let (n, dt) = (501, 1.0 / 2048.0);
    let truth = [tone(137.2, 159.4), tone(161.7, 132.6)];
    let seeds = [tone(137.57, 159.16), tone(161.41, 132.92)];
    let a = C64::new(0.8, -0.4);
    let b = C64::new(-0.2, 0.65);
    let x = signal(n, dt, &truth, &[a, C64::default()], 2.1);
    let y = signal(n, dt, &truth, &[C64::default(), b], -0.7);
    for scale in [1e-140, 1.0, 1e140] {
        let first: Vec<_> = x
            .iter()
            .zip(&y)
            .map(|(x, y)| scale * (0.6 * x + 0.8 * y))
            .collect();
        let second: Vec<_> = x
            .iter()
            .zip(&y)
            .map(|(x, y)| scale * (-0.8 * x + 0.6 * y))
            .collect();
        let constant = vec![3.6 * scale; n];
        let result = refine_chirped_tones_channels(
            &[&first, &second, &constant],
            dt,
            &seeds,
            &bounds(&truth),
            budget(),
        )
        .unwrap();
        close_tones(&result.tones, &truth, 2e-8);
        let expected = [
            [
                C64::new(a.re * 0.6, a.im * 0.6),
                C64::new(b.re * 0.8, b.im * 0.8),
            ],
            [
                C64::new(a.re * -0.8, a.im * -0.8),
                C64::new(b.re * 0.6, b.im * 0.6),
            ],
        ];
        for (fit, expected) in result.fits[..2].iter().zip(expected) {
            for (actual, p) in fit.phasors.iter().zip(expected) {
                close(actual.re / scale, p.re, 1e-9);
                close(actual.im / scale, p.im, 1e-9);
            }
            assert!(fit.residual_rms / scale < 1e-10);
        }
        assert_eq!(result.fits[2].dc, constant[0]);
        assert_eq!(result.fits[2].phasors, vec![C64::default(); 2]);
    }
}

#[test]
fn record_reslicing_changes_endpoint_frequencies_and_rotates_the_initial_phasor() {
    let (n, dt, offset, kept) = (701, 1.0 / 2048.0, 73, 401);
    let p = tone(120.0, 167.0);
    let a = C64::new(0.9, -0.4);
    let all = signal(n, dt, &[p], &[a], 1.3);
    let shift = offset as f64 * dt;
    let rate = (p.end_frequency_hz - p.start_frequency_hz) / ((n - 1) as f64 * dt);
    let local = tone(
        p.start_frequency_hz + rate * shift,
        p.start_frequency_hz + rate * (shift + (kept - 1) as f64 * dt),
    );
    let phase =
        core::f64::consts::TAU * (p.start_frequency_hz * shift + 0.5 * rate * shift * shift);
    let (s, c) = phase.sin_cos();
    let rotated = C64::new(a.re * c + a.im * s, a.im * c - a.re * s);
    let seed = tone(local.start_frequency_hz + 0.3, local.end_frequency_hz - 0.2);
    let result = refine_chirped_tones(
        &all[offset..offset + kept],
        dt,
        &[seed],
        &bounds(&[local]),
        budget(),
    )
    .unwrap();
    close_tones(&result.tones, &[local], 2e-8);
    close(result.fit.phasors[0].re, rotated.re, 1e-9);
    close(result.fit.phasors[0].im, rotated.im, 1e-9);
    close(result.fit.dc, 1.3, 1e-10);
}

#[test]
fn zero_drift_matches_stationary_bits_and_endpoint_units_avoid_rate_underflow() {
    let (n, dt) = (257, 1.0 / 2048.0);
    let frequencies = [137.2187, 300.125];
    let trajectories = frequencies.map(|f| tone(f, f));
    let samples = signal(
        n,
        dt,
        &trajectories,
        &[C64::new(0.7, 0.3), C64::new(-0.2, 0.9)],
        1.2,
    );
    let stationary = ToneFitPlan::new(n, dt, &frequencies, budget().fit).unwrap();
    let chirped = Plan::new(n, dt, &trajectories, budget().fit).unwrap();
    assert_eq!(
        stationary.reciprocal_condition().to_bits(),
        chirped.reciprocal_condition().to_bits()
    );
    assert_eq!(
        stationary.fit(&samples).unwrap(),
        chirped.fit(&samples).unwrap()
    );
    assert_eq!(
        chirped.workspace_bytes(),
        stationary.workspace_bytes() + 2 * core::mem::size_of::<f64>()
    );

    let p = tone(1.0, 3.0);
    let samples = signal(n, 0.125, &[p], &[C64::new(0.8, -0.4)], 0.7);
    // The same dimensionless sampled phase, although Hz/s would underflow.
    let slow = tone(1e-200, 3e-200);
    let slow_dt = 0.125e200;
    assert_eq!(
        (slow.end_frequency_hz - slow.start_frequency_hz) / ((n - 1) as f64 * slow_dt),
        0.0
    );
    let fit = Plan::new(n, slow_dt, &[slow], budget().fit)
        .unwrap()
        .fit(&samples)
        .unwrap();
    close(fit.dc, 0.7, 2e-11);
    close(fit.phasors[0].re, 0.8, 2e-11);
    close(fit.phasors[0].im, -0.4, 2e-11);
    assert!(fit.residual_rms < 2e-11);
}

#[test]
fn omitted_energy_stays_in_the_residual_of_the_joint_frequency_law() {
    let (n, dt) = (501, 1.0 / 2048.0);
    let truth = [tone(137.2, 159.4), tone(161.7, 132.6)];
    let mut samples = signal(
        n,
        dt,
        &truth,
        &[C64::new(0.8, -0.4), C64::new(-0.2, 0.65)],
        0.7,
    );
    for (j, x) in samples.iter_mut().enumerate() {
        *x += 0.015 * (core::f64::consts::TAU * 257.3 * j as f64 * dt).cos();
    }
    let mut b = budget();
    b.absolute_frequency_tolerance_hz = 1e-7;
    let result = refine_chirped_tones(&samples, dt, &truth, &bounds(&truth), b).unwrap();
    assert!(result.fit.residual_rms > 0.009 && result.fit.residual_rms < 0.012);
}

#[test]
fn bounds_conditioning_and_work_refuse_unresolved_or_unconverged_frequency_laws() {
    let (n, dt) = (501, 1.0 / 2048.0);
    let truth = [tone(137.2, 159.4), tone(161.7, 132.6)];
    let samples = signal(
        n,
        dt,
        &truth,
        &[C64::new(0.8, -0.4), C64::new(-0.2, 0.65)],
        0.7,
    );
    assert!(matches!(
        refine_chirped_tones(&vec![3.0; n], dt, &truth, &bounds(&truth), budget()),
        Err(Error::Unidentifiable { .. })
    ));
    let absent = signal(n, dt, &truth[..1], &[C64::new(0.8, -0.4)], 0.7);
    assert!(matches!(
        refine_chirped_tones(&absent, dt, &truth, &bounds(&truth), budget()),
        Err(Error::Unidentifiable { .. })
    ));
    assert!(matches!(
        refine_chirped_tones(&samples[..8], dt, &truth, &bounds(&truth), budget()),
        Err(Error::InvalidParameter { .. })
    ));

    let mut b = budget();
    b.max_evaluations = 1;
    assert_eq!(
        refine_chirped_tones(&samples, dt, &truth, &bounds(&truth), b)
            .unwrap()
            .evaluations,
        1
    );
    let seeds = [tone(137.57, 159.16), tone(161.41, 132.92)];
    assert!(matches!(
        refine_chirped_tones(&samples, dt, &seeds, &bounds(&truth), b),
        Err(Error::NotConverged { evaluations: 1, .. })
    ));
    b = budget();
    b.fit.max_workspace_bytes = 1;
    assert!(matches!(
        refine_chirped_tones(&samples, dt, &truth, &bounds(&truth), b),
        Err(Error::Fit(HarmonicFitError::BudgetExceeded { .. }))
    ));

    let single = signal(n, dt, &truth[..1], &[C64::new(0.8, -0.4)], 0.7);
    let wrong_box = [Bounds {
        start_frequency_hz: [137.3, 138.0],
        end_frequency_hz: [159.5, 160.0],
    }];
    assert!(matches!(
        refine_chirped_tones(&single, dt, &[tone(137.6, 159.8)], &wrong_box, budget()),
        Err(Error::NotConverged { .. })
    ));
}

#[test]
fn both_instantaneous_endpoints_and_tone_identity_are_admitted_before_fitting() {
    let (n, dt) = (31, 1.0 / 1024.0);
    for tones in [
        vec![],
        vec![tone(20.0, 512.0)],
        vec![tone(20.0, 0.0)],
        vec![tone(20.0, f64::NAN)],
        vec![tone(20.0, 30.0); 2],
    ] {
        assert!(Plan::new(n, dt, &tones, budget().fit).is_err());
    }
    let truth = [tone(120.0, 150.0), tone(180.0, 140.0)];
    let samples = signal(
        n,
        dt,
        &truth,
        &[C64::new(0.8, -0.4), C64::new(-0.2, 0.65)],
        0.7,
    );
    let overlap = vec![
        Bounds {
            start_frequency_hz: [100.0, 200.0],
            end_frequency_hz: [100.0, 200.0]
        };
        2
    ];
    assert!(matches!(
        refine_chirped_tones(&samples, dt, &truth, &overlap, budget()),
        Err(Error::InvalidParameter { .. })
    ));
    let mut invalid_end = bounds(&truth);
    invalid_end[0].end_frequency_hz[1] = 512.0;
    assert!(matches!(
        refine_chirped_tones(&samples, dt, &truth, &invalid_end, budget()),
        Err(Error::InvalidParameter { .. })
    ));
    assert!(matches!(
        refine_chirped_tones_channels(
            &[&samples, &samples[..30]],
            dt,
            &truth,
            &bounds(&truth),
            budget()
        ),
        Err(Error::InvalidParameter { .. })
    ));
    let plan = Plan::new(n, dt, &truth, budget().fit).unwrap();
    let mut bad = samples;
    bad[7] = f64::NAN;
    let mut output = [C64::new(3.0, -4.0); 2];
    assert!(plan.fit_into(&bad, &mut vec![0.0; n], &mut output).is_err());
    assert_eq!(output, [C64::new(3.0, -4.0); 2]);
}
