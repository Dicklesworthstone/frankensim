//! G1/G3: amplitude decay and instantaneous-frequency drift must be fitted together.

use fs_fft::C64;
use fs_fft::harmonic::{
    ChirpedTone, ChirpedToneFitPlan, DampedChirpedTone as Tone, DampedChirpedToneBounds as Bounds,
    DampedChirpedToneFitPlan as Plan, DampedChirpedToneRefinementBudget as Budget,
    DampedChirpedToneRefinementError as Error, DampedToneFitPlan, DampedTonePole,
    HarmonicFitBudget, HarmonicFitError, refine_damped_chirped_tones,
    refine_damped_chirped_tones_channels,
};

fn budget() -> Budget {
    Budget {
        absolute_frequency_tolerance_hz: 1e-8,
        absolute_decay_tolerance_per_s: 1e-7,
        max_evaluations: 48,
        min_trajectory_reciprocal_condition: 1e-8,
        fit: HarmonicFitBudget {
            max_workspace_bytes: 8 * 1024 * 1024,
            min_reciprocal_condition: 1e-11,
        },
    }
}

fn tone(start: f64, end: f64, decay: f64) -> Tone {
    Tone {
        start_frequency_hz: start,
        end_frequency_hz: end,
        decay_rate_per_s: decay,
    }
}

fn bounds(tones: &[Tone]) -> Vec<Bounds> {
    tones
        .iter()
        .map(|p| Bounds {
            start_frequency_hz: [p.start_frequency_hz - 1.0, p.start_frequency_hz + 1.0],
            end_frequency_hz: [p.end_frequency_hz - 1.0, p.end_frequency_hz + 1.0],
            decay_rate_per_s: [
                (p.decay_rate_per_s - 1.0).max(0.0),
                p.decay_rate_per_s + 1.0,
            ],
        })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (tolerance {tolerance:e})"
    );
}

fn close_tones(actual: &[Tone], expected: &[Tone]) {
    assert_eq!(actual.len(), expected.len());
    for (a, b) in actual.iter().zip(expected) {
        close(a.start_frequency_hz, b.start_frequency_hz, 2e-7);
        close(a.end_frequency_hz, b.end_frequency_hz, 2e-7);
        close(a.decay_rate_per_s, b.decay_rate_per_s, 2e-6);
    }
}

// Independent quadratic phase in seconds; production forms endpoint phase in
// DD and reduces cycles, and uses deterministic exponential/trigonometry owners.
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
                    (-p.decay_rate_per_s * t).exp() * (a.re * c + a.im * s)
                })
                .sum::<f64>()
        })
        .collect()
}

#[test]
fn decaying_crossing_trajectories_recover_both_endpoint_frequencies_and_initial_amplitudes() {
    let (n, dt) = (601, 1.0 / 2048.0);
    let truth = [tone(137.2, 159.4, 2.7), tone(161.7, 132.6, 5.1)];
    let seeds = [tone(137.57, 159.16, 2.9), tone(161.41, 132.92, 4.8)];
    let phasors = [C64::new(0.8, -0.4), C64::new(-0.2, 0.65)];
    let samples = signal(n, dt, &truth, &phasors, 0.7);
    let plan = Plan::new(n, dt, &truth, budget().fit).unwrap();
    assert_eq!(plan.tones(), truth);
    assert_eq!(plan.sample_count(), n);
    assert_eq!(plan.timestep_s(), dt);
    let supplied = plan.fit(&samples).unwrap();
    assert!(supplied.residual_rms < 2e-12);
    close(supplied.dc, 0.7, 2e-12);
    let mut output = [C64::default(); 2];
    let summary = plan
        .fit_into(&samples, &mut vec![0.0; n], &mut output)
        .unwrap();
    assert_eq!(output.as_slice(), supplied.phasors.as_slice());
    assert_eq!(
        summary.residual_rms.to_bits(),
        supplied.residual_rms.to_bits()
    );

    let undamped = truth.map(|p| ChirpedTone {
        start_frequency_hz: p.start_frequency_hz,
        end_frequency_hz: p.end_frequency_hz,
    });
    let stationary = truth.map(|p| DampedTonePole {
        frequency_hz: p.start_frequency_hz.midpoint(p.end_frequency_hz),
        decay_rate_per_s: p.decay_rate_per_s,
    });
    assert!(
        ChirpedToneFitPlan::new(n, dt, &undamped, budget().fit)
            .unwrap()
            .fit(&samples)
            .unwrap()
            .residual_rms
            > 0.03
    );
    assert!(
        DampedToneFitPlan::new(n, dt, &stationary, budget().fit)
            .unwrap()
            .fit(&samples)
            .unwrap()
            .residual_rms
            > 0.03
    );

    let result =
        refine_damped_chirped_tones(&samples, dt, &seeds, &bounds(&truth), budget()).unwrap();
    close_tones(&result.tones, &truth);
    close(result.fit.dc, 0.7, 1e-10);
    for (a, b) in result.fit.phasors.iter().zip(phasors) {
        close(a.re, b.re, 1e-9);
        close(a.im, b.im, 1e-9);
    }
    assert!(result.fit.residual_rms < 1e-10);
    assert!(result.evaluations <= 12);
    assert!(result.trajectory_reciprocal_condition >= budget().min_trajectory_reciprocal_condition);
    assert!(result.last_corrections.iter().all(|d| {
        d.start_frequency_hz.abs() <= budget().absolute_frequency_tolerance_hz
            && d.end_frequency_hz.abs() <= budget().absolute_frequency_tolerance_hz
            && d.decay_rate_per_s.abs() <= budget().absolute_decay_tolerance_per_s
    }));
    assert_eq!(
        result,
        refine_damped_chirped_tones(&samples, dt, &seeds, &bounds(&truth), budget()).unwrap()
    );
}

#[test]
fn identical_frequency_laws_with_distinct_decay_remain_separate_and_ordered() {
    let (n, dt) = (2001, 1.0 / 4096.0);
    let truth = [tone(418.3, 423.1, 6.0), tone(418.3, 423.1, 1.0)];
    let seeds = [tone(418.36, 423.06, 5.5), tone(418.25, 423.13, 1.3)];
    let domains = [
        Bounds {
            start_frequency_hz: [417.8, 418.8],
            end_frequency_hz: [422.6, 423.6],
            decay_rate_per_s: [4.0, 8.0],
        },
        Bounds {
            start_frequency_hz: [417.8, 418.8],
            end_frequency_hz: [422.6, 423.6],
            decay_rate_per_s: [0.3, 2.0],
        },
    ];
    let phasors = [C64::new(1.0, -0.3), C64::new(0.45, 0.8)];
    let samples = signal(n, dt, &truth, &phasors, 0.3);
    let result = refine_damped_chirped_tones(&samples, dt, &seeds, &domains, budget()).unwrap();
    close_tones(&result.tones, &truth);
    for (a, b) in result.fit.phasors.iter().zip(phasors) {
        close(a.re, b.re, 1e-7);
        close(a.im, b.im, 1e-7);
    }
    assert!(result.fit.residual_rms < 1e-9);
    let reversed = refine_damped_chirped_tones(
        &samples,
        dt,
        &[seeds[1], seeds[0]],
        &[domains[1], domains[0]],
        budget(),
    )
    .unwrap();
    close_tones(&reversed.tones, &[truth[1], truth[0]]);
}

#[test]
fn simultaneous_nodal_sensors_keep_common_unit_weights_under_frame_changes() {
    let (n, dt) = (601, 1.0 / 2048.0);
    let truth = [tone(137.2, 159.4, 2.7), tone(161.7, 132.6, 5.1)];
    let seeds = [tone(137.57, 159.16, 2.9), tone(161.41, 132.92, 4.8)];
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
        let result = refine_damped_chirped_tones_channels(
            &[&first, &second, &constant],
            dt,
            &seeds,
            &bounds(&truth),
            budget(),
        )
        .unwrap();
        close_tones(&result.tones, &truth);
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
fn reslicing_translates_frequency_endpoints_and_attenuates_and_rotates_the_phasor() {
    let (n, dt, offset, kept) = (901, 1.0 / 2048.0, 117, 601);
    let p = tone(231.17, 246.8, 3.7);
    let a = C64::new(0.73, -0.41);
    let all = signal(n, dt, &[p], &[a], 1.3);
    let shift = offset as f64 * dt;
    let rate = (p.end_frequency_hz - p.start_frequency_hz) / ((n - 1) as f64 * dt);
    let local = tone(
        p.start_frequency_hz + rate * shift,
        p.start_frequency_hz + rate * (shift + (kept - 1) as f64 * dt),
        p.decay_rate_per_s,
    );
    let phase =
        core::f64::consts::TAU * (p.start_frequency_hz * shift + 0.5 * rate * shift * shift);
    let (s, c) = phase.sin_cos();
    let attenuation = (-p.decay_rate_per_s * shift).exp();
    let expected = C64::new(
        attenuation * (a.re * c + a.im * s),
        attenuation * (a.im * c - a.re * s),
    );
    let seed = tone(
        local.start_frequency_hz + 0.3,
        local.end_frequency_hz - 0.2,
        3.5,
    );
    let result = refine_damped_chirped_tones(
        &all[offset..offset + kept],
        dt,
        &[seed],
        &bounds(&[local]),
        budget(),
    )
    .unwrap();
    close_tones(&result.tones, &[local]);
    close(result.fit.phasors[0].re, expected.re, 1e-9);
    close(result.fit.phasors[0].im, expected.im, 1e-9);
    close(result.fit.dc, 1.3, 1e-10);
    assert!(result.fit.residual_rms < 1e-10);
}

#[test]
fn zero_decay_and_zero_drift_preserve_the_existing_supplied_law_fit_bits() {
    let (n, dt) = (257, 1.0 / 2048.0);
    let trajectories = [tone(137.2187, 144.6, 0.0), tone(300.125, 297.0, 0.0)];
    let samples = signal(
        n,
        dt,
        &trajectories,
        &[C64::new(0.7, 0.3), C64::new(-0.2, 0.9)],
        1.2,
    );
    let chirps = trajectories.map(|p| ChirpedTone {
        start_frequency_hz: p.start_frequency_hz,
        end_frequency_hz: p.end_frequency_hz,
    });
    let old = ChirpedToneFitPlan::new(n, dt, &chirps, budget().fit).unwrap();
    let new = Plan::new(n, dt, &trajectories, budget().fit).unwrap();
    assert_eq!(
        old.reciprocal_condition().to_bits(),
        new.reciprocal_condition().to_bits()
    );
    assert_eq!(old.fit(&samples).unwrap(), new.fit(&samples).unwrap());
    assert_eq!(
        new.workspace_bytes(),
        old.workspace_bytes() + 2 * core::mem::size_of::<f64>()
    );

    let trajectories = [tone(137.2187, 137.2187, 1.2), tone(300.125, 300.125, 4.1)];
    let samples = signal(
        n,
        dt,
        &trajectories,
        &[C64::new(0.7, 0.3), C64::new(-0.2, 0.9)],
        1.2,
    );
    let poles = trajectories.map(|p| DampedTonePole {
        frequency_hz: p.start_frequency_hz,
        decay_rate_per_s: p.decay_rate_per_s,
    });
    let old = DampedToneFitPlan::new(n, dt, &poles, budget().fit).unwrap();
    let new = Plan::new(n, dt, &trajectories, budget().fit).unwrap();
    assert_eq!(
        old.reciprocal_condition().to_bits(),
        new.reciprocal_condition().to_bits()
    );
    assert_eq!(old.fit(&samples).unwrap(), new.fit(&samples).unwrap());
    assert_eq!(
        new.workspace_bytes(),
        old.workspace_bytes() + 2 * core::mem::size_of::<f64>()
    );
}

#[test]
fn physical_clock_rescaling_retains_decay_and_frequency_when_hz_per_second_underflows() {
    let (n, dt) = (257, 0.125);
    let p = tone(1.0, 3.0, 0.04);
    let a = C64::new(0.8, -0.4);
    let samples = signal(n, dt, &[p], &[a], 0.7);
    let time_scale = 1e200;
    let scaled = tone(
        p.start_frequency_hz / time_scale,
        p.end_frequency_hz / time_scale,
        p.decay_rate_per_s / time_scale,
    );
    let scaled_dt = dt * time_scale;
    assert_eq!(
        (scaled.end_frequency_hz - scaled.start_frequency_hz) / ((n - 1) as f64 * scaled_dt),
        0.0
    );
    let plan = Plan::new(n, scaled_dt, &[scaled], budget().fit).unwrap();
    let fit = plan.fit(&samples).unwrap();
    close(fit.phasors[0].re, a.re, 3e-11);
    close(fit.phasors[0].im, a.im, 3e-11);
    assert!(fit.residual_rms < 3e-11);
    let domain = [Bounds {
        start_frequency_hz: [0.9 / time_scale, 1.1 / time_scale],
        end_frequency_hz: [2.9 / time_scale, 3.1 / time_scale],
        decay_rate_per_s: [0.01 / time_scale, 0.08 / time_scale],
    }];
    let mut b = budget();
    b.absolute_frequency_tolerance_hz /= time_scale;
    b.absolute_decay_tolerance_per_s /= time_scale;
    let seed = tone(1.003 / time_scale, 2.997 / time_scale, 0.043 / time_scale);
    let result = refine_damped_chirped_tones(&samples, scaled_dt, &[seed], &domain, b).unwrap();
    close(result.tones[0].start_frequency_hz * time_scale, 1.0, 2e-8);
    close(result.tones[0].end_frequency_hz * time_scale, 3.0, 2e-8);
    close(result.tones[0].decay_rate_per_s * time_scale, 0.04, 2e-7);
}

#[test]
fn three_coordinate_rank_work_storage_and_active_bounds_refuse_without_inventing_fits() {
    let (n, dt) = (601, 1.0 / 2048.0);
    let truth = [tone(137.2, 159.4, 2.7), tone(161.7, 132.6, 5.1)];
    let a = [C64::new(0.8, -0.4), C64::new(-0.2, 0.65)];
    let samples = signal(n, dt, &truth, &a, 0.7);
    assert!(matches!(
        refine_damped_chirped_tones(&vec![3.0; n], dt, &truth, &bounds(&truth), budget()),
        Err(Error::Unidentifiable { .. })
    ));
    let absent = signal(n, dt, &truth[..1], &a[..1], 0.7);
    assert!(matches!(
        refine_damped_chirped_tones(&absent, dt, &truth, &bounds(&truth), budget()),
        Err(Error::Unidentifiable { .. })
    ));
    // Two trajectories need SIX residual dimensions, beyond their five linear
    // amplitude/DC columns. The previous two-coordinate minimum would admit N=10.
    assert!(matches!(
        refine_damped_chirped_tones(&samples[..10], dt, &truth, &bounds(&truth), budget()),
        Err(Error::InvalidParameter { .. })
    ));

    let mut b = budget();
    b.max_evaluations = 1;
    assert_eq!(
        refine_damped_chirped_tones(&samples, dt, &truth, &bounds(&truth), b)
            .unwrap()
            .evaluations,
        1
    );
    let seeds = [tone(137.57, 159.16, 2.9), tone(161.41, 132.92, 4.8)];
    assert!(matches!(
        refine_damped_chirped_tones(&samples, dt, &seeds, &bounds(&truth), b),
        Err(Error::NotConverged { evaluations: 1, .. })
    ));
    b.fit.max_workspace_bytes = 0;
    let Error::Fit(HarmonicFitError::BudgetExceeded { required_bytes, .. }) =
        refine_damped_chirped_tones(&samples, dt, &truth, &bounds(&truth), b).unwrap_err()
    else {
        panic!("numeric storage must be admitted before fitting");
    };
    b.fit.max_workspace_bytes = required_bytes;
    assert_eq!(
        refine_damped_chirped_tones(&samples, dt, &truth, &bounds(&truth), b)
            .unwrap()
            .evaluations,
        1
    );
    b.fit.max_workspace_bytes -= 1;
    assert!(matches!(
        refine_damped_chirped_tones(&samples, dt, &truth, &bounds(&truth), b),
        Err(Error::Fit(HarmonicFitError::BudgetExceeded { .. }))
    ));

    // Frequencies alone can converge while the true decay remains outside its
    // allowed interval; that active decay bound must never be called success.
    let mut wrong = bounds(&truth[..1]);
    wrong[0].decay_rate_per_s = [3.1, 4.0];
    assert!(matches!(
        refine_damped_chirped_tones(&absent, dt, &[tone(137.2, 159.4, 3.5)], &wrong, budget()),
        Err(Error::NotConverged { .. })
    ));
}

#[test]
fn endpoint_decay_identity_and_transactional_output_admission_remain_explicit() {
    let (n, dt) = (31, 1.0 / 1024.0);
    for tones in [
        vec![],
        vec![tone(20.0, 512.0, 1.0)],
        vec![tone(0.0, 30.0, 1.0)],
        vec![tone(20.0, 30.0, -1.0)],
        vec![tone(20.0, 30.0, f64::NAN)],
        vec![tone(20.0, 30.0, 1.0); 2],
    ] {
        assert!(Plan::new(n, dt, &tones, budget().fit).is_err());
    }
    let truth = [tone(120.0, 150.0, 2.0), tone(180.0, 140.0, 3.0)];
    let samples = signal(
        n,
        dt,
        &truth,
        &[C64::new(0.8, -0.4), C64::new(-0.2, 0.65)],
        0.7,
    );
    let overlap = [Bounds {
        start_frequency_hz: [100.0, 200.0],
        end_frequency_hz: [100.0, 200.0],
        decay_rate_per_s: [0.0, 10.0],
    }; 2];
    assert!(matches!(
        refine_damped_chirped_tones(&samples, dt, &truth, &overlap, budget()),
        Err(Error::InvalidParameter { .. })
    ));
    let mut bad_domain = bounds(&truth);
    bad_domain[0].decay_rate_per_s[0] = -0.1;
    assert!(matches!(
        refine_damped_chirped_tones(&samples, dt, &truth, &bad_domain, budget()),
        Err(Error::InvalidParameter { .. })
    ));
    assert!(matches!(
        refine_damped_chirped_tones_channels(
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
