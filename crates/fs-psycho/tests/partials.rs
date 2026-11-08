use core::{f64::consts::TAU, ops::ControlFlow};
use fs_fft::C64;
use fs_fft::harmonic::{
    DampedToneBounds, DampedTonePole, DampedToneRefinementBudget, HarmonicFitBudget,
    HarmonicFitError, ToneFrequencyBudget, ToneFrequencyError,
};
use fs_psycho::partials::{
    PartialFitFailure, PartialFrame, PartialFrameSchedule, PartialTrackingError,
    PartialTrackingModel, track_partials,
};

fn tone_budget() -> ToneFrequencyBudget {
    ToneFrequencyBudget {
        absolute_tolerance_hz: 1e-8,
        max_evaluations: 80,
        min_frequency_reciprocal_condition: 1e-9,
        fit: HarmonicFitBudget {
            max_workspace_bytes: 8 * 1024 * 1024,
            min_reciprocal_condition: 1e-11,
        },
    }
}

fn pole_budget() -> DampedToneRefinementBudget {
    DampedToneRefinementBudget {
        absolute_frequency_tolerance_hz: 1e-8,
        absolute_decay_tolerance_per_s: 1e-7,
        max_evaluations: 80,
        min_pole_reciprocal_condition: 1e-10,
        fit: tone_budget().fit,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn moving_inharmonic_partials_keep_identity_phase_and_exact_source_indices() {
    let dt = 1.0 / 8192.0;
    let n = 1025;
    let frequencies = [[221.3, 603.2], [222.1, 605.0], [220.7, 604.3]];
    let amplitudes = [[1.0, 0.43], [0.7, 0.51], [0.9, 0.3]];
    // Piecewise constant physical frequencies, with continuous carrier phase
    // through the changes. The two partials are deliberately not harmonic.
    let mut phases = [0.27, -0.41];
    let mut initial_phasors = Vec::new();
    let mut samples = Vec::new();
    for (f, a) in frequencies.iter().zip(amplitudes) {
        initial_phasors.push(core::array::from_fn::<_, 2, _>(|k| {
            let (s, c) = f64::sin_cos(phases[k]);
            C64::new(a[k] * c, -a[k] * s)
        }));
        for i in 0..n {
            let t = i as f64 * dt;
            samples.push(
                0.12 + (0..2)
                    .map(|k| a[k] * (TAU * f[k] * t + phases[k]).cos())
                    .sum::<f64>(),
            );
        }
        for k in 0..2 {
            phases[k] += TAU * f[k] * n as f64 * dt;
        }
    }
    let schedule = PartialFrameSchedule {
        frame_samples: n,
        hop_samples: n,
        max_frames: 3,
        first_sample_index: (1_u64 << 53) + 37,
    };
    let model = PartialTrackingModel::Stationary {
        initial_frequencies_hz: &[221.0, 604.0],
        frequency_bounds_hz: &[[218.0, 225.0], [599.0, 609.0]],
        budget: tone_budget(),
    };
    let mut frames = Vec::new();
    let summary = track_partials(&samples, dt, schedule, model, |frame| {
        frames.push(frame.clone());
        ControlFlow::Continue(())
    })
    .unwrap();
    assert!(summary.completed);
    assert_eq!(summary.planned_frames, 3);
    assert_eq!(summary.visited_frames, 3);
    assert_eq!(summary.resolved_frames, 3);
    assert_eq!(summary.maximum_candidate_evaluations, 240);
    for (j, frame) in frames.iter().enumerate() {
        assert_eq!(
            frame.first_sample_index,
            schedule.first_sample_index + (j * n) as u64
        );
        assert_eq!(frame.sample_count, n);
        let estimate = frame.estimate.as_ref().unwrap();
        close(estimate.fit().dc, 0.12, 2e-10);
        assert!(estimate.fit().residual_rms < 2e-10);
        assert!(estimate.partial(2).is_none());
        for k in 0..2 {
            let partial = estimate.partial(k).unwrap();
            close(partial.frequency_hz, frequencies[j][k], 2e-8);
            assert_eq!(partial.decay_rate_per_s, None);
            close(partial.peak_phasor.re, initial_phasors[j][k].re, 2e-8);
            close(partial.peak_phasor.im, initial_phasors[j][k].im, 2e-8);
        }
    }
    let mut replay = Vec::new();
    track_partials(&samples, dt, schedule, model, |frame| {
        replay.push(frame.clone());
        ControlFlow::Continue(())
    })
    .unwrap();
    assert_eq!(
        frames, replay,
        "identical declarations replay deterministically"
    );
}

#[test]
fn overlapping_frames_resolve_double_decay_and_cover_the_final_sample() {
    let dt = 1.0 / 2048.0;
    let poles = [
        DampedTonePole {
            frequency_hz: 211.37,
            decay_rate_per_s: 1.1,
        },
        DampedTonePole {
            frequency_hz: 211.37,
            decay_rate_per_s: 4.6,
        },
    ];
    let phasors = [C64::new(0.8, 0.3), C64::new(0.4, -0.25)];
    let samples: Vec<_> = (0..1299)
        .map(|i| {
            let t = i as f64 * dt;
            0.23 + poles
                .iter()
                .zip(phasors)
                .map(|(pole, c)| {
                    let (s, co) = (TAU * pole.frequency_hz * t).sin_cos();
                    (-pole.decay_rate_per_s * t).exp() * (c.re * co + c.im * s)
                })
                .sum::<f64>()
        })
        .collect();
    let initial = poles.map(|p| DampedTonePole {
        frequency_hz: p.frequency_hz + 0.02,
        decay_rate_per_s: p.decay_rate_per_s + 0.1,
    });
    let domains = [
        DampedToneBounds {
            frequency_hz: [211.0, 211.7],
            decay_rate_per_s: [0.5, 2.0],
        },
        DampedToneBounds {
            frequency_hz: [211.0, 211.7],
            decay_rate_per_s: [3.5, 6.0],
        },
    ];
    let schedule = PartialFrameSchedule {
        frame_samples: 769,
        hop_samples: 400,
        max_frames: 3,
        first_sample_index: 127,
    };
    let model = PartialTrackingModel::Ringdown {
        initial_poles: &initial,
        bounds: &domains,
        budget: pole_budget(),
    };
    let mut frames = Vec::new();
    let summary = track_partials(&samples, dt, schedule, model, |frame| {
        frames.push(frame.clone());
        ControlFlow::Continue(())
    })
    .unwrap();
    assert_eq!(summary.resolved_frames, 3);
    assert!(summary.completed);
    let starts: Vec<_> = frames
        .iter()
        .map(|frame| frame.first_sample_index)
        .collect();
    assert_eq!(starts, [127, 527, 657]);
    assert_eq!(
        frames.last().unwrap().first_sample_index + 769,
        127 + samples.len() as u64
    );
    for frame in &frames {
        let t = (frame.first_sample_index - 127) as f64 * dt;
        let estimate = frame.estimate.as_ref().unwrap();
        close(estimate.fit().dc, 0.23, 1e-9);
        assert!(estimate.fit().residual_rms < 1e-9);
        for k in 0..2 {
            let partial = estimate.partial(k).unwrap();
            close(partial.frequency_hz, poles[k].frequency_hz, 1e-7);
            close(
                partial.decay_rate_per_s.unwrap(),
                poles[k].decay_rate_per_s,
                2e-6,
            );
            let (s, c) = (TAU * poles[k].frequency_hz * t).sin_cos();
            let envelope = (-poles[k].decay_rate_per_s * t).exp();
            // Local t=0 shifts BOTH amplitude and carrier phase. No averaging
            // or incoherent addition may erase the two-polarization beat/decay.
            close(
                partial.peak_phasor.re,
                envelope * (phasors[k].re * c + phasors[k].im * s),
                5e-7,
            );
            close(
                partial.peak_phasor.im,
                envelope * (phasors[k].im * c - phasors[k].re * s),
                5e-7,
            );
        }
    }
}

fn dropout_fixture() -> (Vec<f64>, PartialFrameSchedule) {
    let n = 513;
    let mut record: Vec<_> = (0..n)
        .map(|i| 0.6 * (TAU * 213.2 * i as f64 / 4096.0).cos())
        .collect();
    record.extend(core::iter::repeat_n(3.0, n));
    record.extend((0..n).map(|i| 0.8 * (TAU * 215.0 * i as f64 / 4096.0 + 0.37).cos()));
    (
        record,
        PartialFrameSchedule {
            frame_samples: n,
            hop_samples: n,
            max_frames: 3,
            first_sample_index: 11,
        },
    )
}

#[test]
fn unresolved_frames_are_explicit_and_do_not_destroy_tracking_seeds() {
    let (record, schedule) = dropout_fixture();
    let model = PartialTrackingModel::Stationary {
        initial_frequencies_hz: &[211.0],
        frequency_bounds_hz: &[[205.0, 220.0]],
        budget: tone_budget(),
    };
    let mut frames: Vec<PartialFrame> = Vec::new();
    let summary = track_partials(&record, 1.0 / 4096.0, schedule, model, |frame| {
        frames.push(frame.clone());
        ControlFlow::Continue(())
    })
    .unwrap();
    assert_eq!(summary.resolved_frames, 2);
    assert!(summary.completed);
    assert!(matches!(
        frames[1].estimate,
        Err(PartialFitFailure::Stationary(
            ToneFrequencyError::Unidentifiable { .. }
        ))
    ));
    close(
        frames[0]
            .estimate
            .as_ref()
            .unwrap()
            .partial(0)
            .unwrap()
            .frequency_hz,
        213.2,
        2e-8,
    );
    let last = frames[2].estimate.as_ref().unwrap().partial(0).unwrap();
    close(last.frequency_hz, 215.0, 2e-8);
    close(last.peak_phasor.re, 0.8 * 0.37_f64.cos(), 1e-8);
    close(last.peak_phasor.im, -0.8 * 0.37_f64.sin(), 1e-8);
}

#[test]
fn work_storage_and_source_clock_are_admitted_before_observation() {
    let (mut record, schedule) = dropout_fixture();
    let model = PartialTrackingModel::Stationary {
        initial_frequencies_hz: &[211.0],
        frequency_bounds_hz: &[[205.0, 220.0]],
        budget: tone_budget(),
    };
    let mut called = 0;
    let error = track_partials(
        &record,
        1.0 / 4096.0,
        PartialFrameSchedule {
            max_frames: 2,
            ..schedule
        },
        model,
        |_| {
            called += 1;
            ControlFlow::Continue(())
        },
    )
    .unwrap_err();
    assert_eq!(
        error,
        PartialTrackingError::FrameBudgetExceeded {
            required_frames: 3,
            max_frames: 2
        }
    );
    assert_eq!(called, 0);
    let small = PartialTrackingModel::Stationary {
        initial_frequencies_hz: &[211.0],
        frequency_bounds_hz: &[[205.0, 220.0]],
        budget: ToneFrequencyBudget {
            fit: HarmonicFitBudget {
                max_workspace_bytes: 7,
                ..tone_budget().fit
            },
            ..tone_budget()
        },
    };
    assert!(matches!(
        track_partials(&record, 1.0 / 4096.0, schedule, small, |_| panic!(
            "storage admission must precede observations"
        )),
        Err(PartialTrackingError::Fit {
            failure: PartialFitFailure::Stationary(ToneFrequencyError::Fit(
                HarmonicFitError::BudgetExceeded {
                    required_bytes: 8,
                    max_bytes: 7
                }
            )),
            ..
        })
    ));
    assert!(matches!(
        track_partials(
            &record,
            1.0 / 4096.0,
            PartialFrameSchedule {
                first_sample_index: u64::MAX - 2,
                ..schedule
            },
            model,
            |_| panic!("overflow must be refused")
        ),
        Err(PartialTrackingError::InvalidParameter { .. })
    ));
    let last = record.len() - 1;
    record[last] = f64::NAN;
    assert!(matches!(
        track_partials(&record, 1.0 / 4096.0, schedule, model, |_| panic!(
            "all source samples admitted first"
        )),
        Err(PartialTrackingError::InvalidParameter { .. })
    ));
}

#[test]
fn frame_callback_stops_work_and_unconverged_frames_never_publish_partials() {
    let (record, schedule) = dropout_fixture();
    let model = PartialTrackingModel::Stationary {
        initial_frequencies_hz: &[211.0],
        frequency_bounds_hz: &[[205.0, 220.0]],
        budget: ToneFrequencyBudget {
            max_evaluations: 1,
            ..tone_budget()
        },
    };
    let mut called = 0;
    let summary = track_partials(&record, 1.0 / 4096.0, schedule, model, |frame| {
        called += 1;
        assert!(matches!(
            frame.estimate,
            Err(PartialFitFailure::Stationary(
                ToneFrequencyError::NotConverged { evaluations: 1, .. }
            ))
        ));
        ControlFlow::Break(())
    })
    .unwrap();
    assert_eq!(called, 1);
    assert_eq!(summary.visited_frames, 1);
    assert_eq!(summary.resolved_frames, 0);
    assert!(!summary.completed);
    assert_eq!(summary.maximum_candidate_evaluations, 3);
}

#[test]
fn continuous_pitch_glide_has_a_measured_frame_length_ladder() {
    let dt = 1.0 / 4096.0;
    let start_hz = 173.37;
    let glide_hz_per_s = 20.0;
    // Instantaneous frequency is the derivative of the integrated phase.
    // This is deliberately not an exactly stationary fit model.
    let record: Vec<_> = (0..1025)
        .map(|i| {
            let t = i as f64 * dt;
            (TAU * (start_hz * t + 0.5 * glide_hz_per_s * t * t) + 0.31).cos()
        })
        .collect();
    let model = PartialTrackingModel::Stationary {
        initial_frequencies_hz: &[175.0],
        frequency_bounds_hz: &[[170.0, 185.0]],
        budget: ToneFrequencyBudget {
            absolute_tolerance_hz: 1e-7,
            ..tone_budget()
        },
    };
    let mut previous_residual = f64::INFINITY;
    for frame_samples in [513, 257, 129] {
        let schedule = PartialFrameSchedule {
            frame_samples,
            hop_samples: frame_samples - 1,
            max_frames: 8,
            first_sample_index: 0,
        };
        let mut max_residual: f64 = 0.0;
        let summary = track_partials(&record, dt, schedule, model, |frame| {
            let estimate = frame.estimate.as_ref().unwrap();
            let center_time =
                (frame.first_sample_index as f64 + 0.5 * (frame.sample_count - 1) as f64) * dt;
            close(
                estimate.partial(0).unwrap().frequency_hz,
                start_hz + glide_hz_per_s * center_time,
                0.03,
            );
            max_residual = max_residual.max(estimate.fit().residual_rms);
            ControlFlow::Continue(())
        })
        .unwrap();
        assert_eq!(summary.resolved_frames, summary.planned_frames);
        assert!(summary.completed);
        assert!(
            max_residual < 0.4 * previous_residual,
            "halving a smooth chirp's frame should reduce local-model residual: {max_residual:e} vs {previous_residual:e}"
        );
        previous_residual = max_residual;
    }
}
