use core::{f64::consts::TAU, ops::ControlFlow};
use fs_fft::C64;
use fs_fft::harmonic::{
    ChirpedTone, ChirpedToneBounds, ChirpedToneRefinementBudget, ChirpedToneRefinementError,
    DampedChirpedTone, DampedChirpedToneBounds, DampedChirpedToneRefinementBudget,
    DampedChirpedToneRefinementError, DampedToneBounds, DampedTonePole, DampedToneRefinementBudget,
    HarmonicFitBudget, HarmonicFitError, ToneFrequencyBudget, ToneFrequencyError,
};
use fs_psycho::partials::{
    PartialFitFailure, PartialFrame, PartialFrameEstimate, PartialFrameSchedule,
    PartialTrackingError, PartialTrackingModel, track_partials,
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

fn chirp_budget() -> ChirpedToneRefinementBudget {
    ChirpedToneRefinementBudget {
        absolute_frequency_tolerance_hz: 1e-8,
        max_evaluations: 32,
        min_chirp_reciprocal_condition: 1e-9,
        fit: tone_budget().fit,
    }
}

fn damped_chirp_budget() -> DampedChirpedToneRefinementBudget {
    DampedChirpedToneRefinementBudget {
        absolute_frequency_tolerance_hz: 1e-8,
        absolute_decay_tolerance_per_s: 1e-7,
        max_evaluations: 80,
        min_trajectory_reciprocal_condition: 1e-9,
        fit: tone_budget().fit,
    }
}

fn chirp_sample(t: f64, first_hz: f64, drift_hz_per_s: f64, phasor: C64) -> f64 {
    let (s, c) = (TAU * (first_hz * t + 0.5 * drift_hz_per_s * t * t)).sin_cos();
    phasor.re * c + phasor.im * s
}

fn translated_chirp_phasor(t: f64, first_hz: f64, drift_hz_per_s: f64, phasor: C64) -> C64 {
    let (s, c) = (TAU * (first_hz * t + 0.5 * drift_hz_per_s * t * t)).sin_cos();
    C64::new(phasor.re * c + phasor.im * s, phasor.im * c - phasor.re * s)
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
            assert_eq!(partial.end_frequency_hz, None);
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
            assert_eq!(partial.end_frequency_hz, None);
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

#[test]
fn chirped_partials_follow_actual_frame_endpoints_and_short_terminal_hops() {
    let dt = 1.0 / 2048.0;
    let frame_samples = 257;
    let first = [173.37, 605.2];
    let drift = [20.0, -9.0];
    let phasors = [C64::new(0.8, -0.3), C64::new(-0.2, 0.5)];
    let record: Vec<_> = (0..1100)
        .map(|j| {
            let t = j as f64 * dt;
            0.12 + (0..2)
                .map(|k| chirp_sample(t, first[k], drift[k], phasors[k]))
                .sum::<f64>()
        })
        .collect();
    let initial = core::array::from_fn::<_, 2, _>(|k| ChirpedTone {
        start_frequency_hz: first[k],
        end_frequency_hz: first[k] + drift[k] * (frame_samples - 1) as f64 * dt,
    });
    let bounds = [
        ChirpedToneBounds {
            start_frequency_hz: [170.0, 190.0],
            end_frequency_hz: [170.0, 190.0],
        },
        ChirpedToneBounds {
            start_frequency_hz: [595.0, 610.0],
            end_frequency_hz: [595.0, 610.0],
        },
    ];
    let schedule = PartialFrameSchedule {
        frame_samples,
        hop_samples: 160,
        max_frames: 7,
        first_sample_index: (1_u64 << 53) + 113,
    };
    let model = PartialTrackingModel::Chirped {
        initial_tones: &initial,
        bounds: &bounds,
        budget: ChirpedToneRefinementBudget {
            // Every warm start must already match the physical trajectory.
            // A wrong N vs N-1 duration or nominal final hop cannot be repaired.
            max_evaluations: 1,
            ..chirp_budget()
        },
    };
    let mut frames = Vec::new();
    let summary = track_partials(&record, dt, schedule, model, |frame| {
        frames.push(frame.clone());
        ControlFlow::Continue(())
    })
    .unwrap();
    assert!(summary.completed);
    assert_eq!(summary.planned_frames, 7);
    assert_eq!(summary.resolved_frames, 7);
    assert_eq!(summary.maximum_candidate_evaluations, 7);
    let offsets = [0, 160, 320, 480, 640, 800, 843];
    for (frame, offset) in frames.iter().zip(offsets) {
        assert_eq!(
            frame.first_sample_index,
            schedule.first_sample_index + offset
        );
        assert_eq!(frame.sample_count, frame_samples);
        let estimate = frame.estimate.as_ref().unwrap();
        let PartialFrameEstimate::Chirped(result) = estimate else {
            panic!("chirped model must retain its full numerical diagnostics");
        };
        assert_eq!(result.evaluations, 1);
        close(estimate.fit().dc, 0.12, 2e-10);
        assert!(estimate.fit().residual_rms < 3e-10);
        assert!(estimate.partial(2).is_none());
        let t = offset as f64 * dt;
        let last_t = (offset + (frame_samples - 1) as u64) as f64 * dt;
        for k in 0..2 {
            let partial = estimate.partial(k).unwrap();
            close(partial.frequency_hz, first[k] + drift[k] * t, 2e-8);
            close(
                partial.end_frequency_hz.unwrap(),
                first[k] + drift[k] * last_t,
                2e-8,
            );
            assert_eq!(partial.decay_rate_per_s, None);
            let expected = translated_chirp_phasor(t, first[k], drift[k], phasors[k]);
            close(partial.peak_phasor.re, expected.re, 2e-8);
            close(partial.peak_phasor.im, expected.im, 2e-8);
        }
    }
    let mut replay = Vec::new();
    track_partials(&record, dt, schedule, model, |frame| {
        replay.push(frame.clone());
        ControlFlow::Continue(())
    })
    .unwrap();
    assert_eq!(frames, replay);
}

#[test]
fn chirp_dropout_retains_successful_anchor_and_full_sample_gap() {
    let n = 513;
    let dt = 1.0 / 4096.0;
    let first = 213.2;
    let drift = 20.0;
    let phasor = C64::new(0.7, 0.2);
    let mut record: Vec<_> = (0..3 * n)
        .map(|j| 0.3 + chirp_sample(j as f64 * dt, first, drift, phasor))
        .collect();
    record[n..2 * n].fill(3.0);
    let initial = [ChirpedTone {
        start_frequency_hz: first,
        end_frequency_hz: first + drift * (n - 1) as f64 * dt,
    }];
    let bounds = [ChirpedToneBounds {
        start_frequency_hz: [205.0, 225.0],
        end_frequency_hz: [205.0, 225.0],
    }];
    let mut frames = Vec::new();
    let summary = track_partials(
        &record,
        dt,
        PartialFrameSchedule {
            frame_samples: n,
            hop_samples: n,
            max_frames: 3,
            first_sample_index: 11,
        },
        PartialTrackingModel::Chirped {
            initial_tones: &initial,
            bounds: &bounds,
            budget: ChirpedToneRefinementBudget {
                max_evaluations: 1,
                ..chirp_budget()
            },
        },
        |frame| {
            frames.push(frame.clone());
            ControlFlow::Continue(())
        },
    )
    .unwrap();
    assert!(summary.completed);
    assert_eq!(summary.resolved_frames, 2);
    assert_eq!(summary.maximum_candidate_evaluations, 3);
    assert!(matches!(
        frames[1].estimate,
        Err(PartialFitFailure::Chirped(
            ChirpedToneRefinementError::Unidentifiable { .. }
        ))
    ));
    let last = frames[2].estimate.as_ref().unwrap().partial(0).unwrap();
    let t = (2 * n) as f64 * dt;
    close(last.frequency_hz, first + drift * t, 2e-8);
    close(
        last.end_frequency_hz.unwrap(),
        first + drift * (3 * n - 1) as f64 * dt,
        2e-8,
    );
    let expected = translated_chirp_phasor(t, first, drift, phasor);
    close(last.peak_phasor.re, expected.re, 2e-8);
    close(last.peak_phasor.im, expected.im, 2e-8);
}

#[test]
fn chirp_warm_prediction_clamps_only_initial_guesses_to_fixed_domains() {
    let n = 257;
    let dt = 1.0 / 2048.0;
    let phasor = C64::new(0.8, -0.3);
    let second_phasor = translated_chirp_phasor(n as f64 * dt, 170.0, 72.0, phasor);
    let mut record: Vec<_> = (0..n)
        .map(|j| 0.1 + chirp_sample(j as f64 * dt, 170.0, 72.0, phasor))
        .collect();
    // This frame has a different, shallower local frequency law. Extrapolating
    // the successful first law predicts an end above the fixed 180 Hz bound.
    record.extend((0..n).map(|j| 0.1 + chirp_sample(j as f64 * dt, 179.0, 4.0, second_phasor)));
    let initial = [ChirpedTone {
        start_frequency_hz: 170.0,
        end_frequency_hz: 179.0,
    }];
    let bounds = [ChirpedToneBounds {
        start_frequency_hz: [169.0, 190.0],
        end_frequency_hz: [170.0, 180.0],
    }];
    let mut frames = Vec::new();
    let summary = track_partials(
        &record,
        dt,
        PartialFrameSchedule {
            frame_samples: n,
            hop_samples: n,
            max_frames: 2,
            first_sample_index: 0,
        },
        PartialTrackingModel::Chirped {
            initial_tones: &initial,
            bounds: &bounds,
            budget: chirp_budget(),
        },
        |frame| {
            frames.push(frame.clone());
            ControlFlow::Continue(())
        },
    )
    .unwrap();
    assert_eq!(summary.resolved_frames, 2);
    let last = frames[1].estimate.as_ref().unwrap();
    let partial = last.partial(0).unwrap();
    close(partial.frequency_hz, 179.0, 2e-8);
    close(partial.end_frequency_hz.unwrap(), 179.5, 2e-8);
    close(partial.peak_phasor.re, second_phasor.re, 2e-8);
    close(partial.peak_phasor.im, second_phasor.im, 2e-8);
    assert!(last.fit().residual_rms < 3e-10);
}

#[test]
fn chirp_tracker_admits_both_seed_arrays_and_keeps_typed_refusals_and_stop() {
    let n = 257;
    let dt = 1.0 / 2048.0;
    let record: Vec<_> = (0..2 * n)
        .map(|j| chirp_sample(j as f64 * dt, 170.0, 4.0, C64::new(0.8, -0.3)))
        .collect();
    let initial = [ChirpedTone {
        start_frequency_hz: 170.0,
        end_frequency_hz: 170.5,
    }];
    let bounds = [ChirpedToneBounds {
        start_frequency_hz: [165.0, 180.0],
        end_frequency_hz: [165.0, 180.0],
    }];
    let schedule = PartialFrameSchedule {
        frame_samples: n,
        hop_samples: n,
        max_frames: 2,
        first_sample_index: 7,
    };
    let model = PartialTrackingModel::Chirped {
        initial_tones: &initial,
        bounds: &bounds,
        budget: chirp_budget(),
    };
    let required_bytes = 2 * core::mem::size_of::<ChirpedTone>();
    let max_bytes = required_bytes - 1;
    let small = PartialTrackingModel::Chirped {
        initial_tones: &initial,
        bounds: &bounds,
        budget: ChirpedToneRefinementBudget {
            fit: HarmonicFitBudget {
                max_workspace_bytes: max_bytes,
                ..chirp_budget().fit
            },
            ..chirp_budget()
        },
    };
    assert_eq!(
        track_partials(&record, dt, schedule, small, |_| panic!(
            "both seed arrays admitted before observation"
        )),
        Err(PartialTrackingError::Fit {
            first_sample_index: 7,
            failure: PartialFitFailure::Chirped(ChirpedToneRefinementError::Fit(
                HarmonicFitError::BudgetExceeded {
                    required_bytes,
                    max_bytes
                }
            )),
        })
    );
    let invalid_bounds = [ChirpedToneBounds {
        start_frequency_hz: [180.0, 165.0],
        ..bounds[0]
    }];
    assert!(matches!(
        track_partials(
            &record,
            dt,
            schedule,
            PartialTrackingModel::Chirped {
                initial_tones: &initial,
                bounds: &invalid_bounds,
                budget: chirp_budget(),
            },
            |_| panic!("initial declarations must refuse before prediction")
        ),
        Err(PartialTrackingError::Fit {
            failure: PartialFitFailure::Chirped(
                ChirpedToneRefinementError::InvalidParameter { .. }
            ),
            ..
        })
    ));
    let mut called = 0;
    let summary = track_partials(&record, dt, schedule, model, |frame| {
        called += 1;
        assert!(frame.estimate.is_ok());
        ControlFlow::Break(())
    })
    .unwrap();
    assert_eq!(called, 1);
    assert_eq!(summary.visited_frames, 1);
    assert_eq!(summary.resolved_frames, 1);
    assert_eq!(summary.maximum_candidate_evaluations, 64);
    assert!(!summary.completed);
}

#[test]
fn damped_chirp_frames_retain_initial_amplitude_phase_and_physical_endpoints() {
    let n = 513;
    let dt = 1.0 / 4096.0;
    let first = [211.37, 603.2];
    let drift = [20.0, -14.0];
    let decay = [1.1, 4.6];
    let phasors = [C64::new(0.8, 0.3), C64::new(0.4, -0.25)];
    let record: Vec<_> = (0..1301)
        .map(|j| {
            let t = j as f64 * dt;
            0.23 + (0..2)
                .map(|k| (-decay[k] * t).exp() * chirp_sample(t, first[k], drift[k], phasors[k]))
                .sum::<f64>()
        })
        .collect();
    let initial = core::array::from_fn::<_, 2, _>(|k| DampedChirpedTone {
        start_frequency_hz: first[k],
        end_frequency_hz: first[k] + drift[k] * (n - 1) as f64 * dt,
        decay_rate_per_s: decay[k],
    });
    let bounds = [
        DampedChirpedToneBounds {
            start_frequency_hz: [205.0, 225.0],
            end_frequency_hz: [205.0, 225.0],
            decay_rate_per_s: [0.5, 2.0],
        },
        DampedChirpedToneBounds {
            start_frequency_hz: [595.0, 610.0],
            end_frequency_hz: [595.0, 610.0],
            decay_rate_per_s: [3.5, 6.0],
        },
    ];
    let schedule = PartialFrameSchedule {
        frame_samples: n,
        hop_samples: 256,
        max_frames: 5,
        first_sample_index: (1_u64 << 53) + 73,
    };
    let mut frames = Vec::new();
    let summary = track_partials(
        &record,
        dt,
        schedule,
        PartialTrackingModel::DampedChirped {
            initial_tones: &initial,
            bounds: &bounds,
            budget: DampedChirpedToneRefinementBudget {
                // Incorrect endpoint span, nominal terminal hop or decay
                // prediction cannot be repaired by a second candidate.
                max_evaluations: 1,
                ..damped_chirp_budget()
            },
        },
        |frame| {
            frames.push(frame.clone());
            ControlFlow::Continue(())
        },
    )
    .unwrap();
    assert!(summary.completed);
    assert_eq!(summary.planned_frames, 5);
    assert_eq!(summary.resolved_frames, 5);
    assert_eq!(summary.maximum_candidate_evaluations, 5);
    for (frame, offset) in frames.iter().zip([0, 256, 512, 768, 788]) {
        assert_eq!(
            frame.first_sample_index,
            schedule.first_sample_index + offset
        );
        let estimate = frame.estimate.as_ref().unwrap();
        let PartialFrameEstimate::DampedChirped(result) = estimate else {
            panic!("combined model retains its full numerical diagnostics");
        };
        assert_eq!(result.evaluations, 1);
        close(estimate.fit().dc, 0.23, 2e-10);
        assert!(estimate.fit().residual_rms < 3e-10);
        assert!(estimate.partial(2).is_none());
        let t = offset as f64 * dt;
        let last_t = (offset + (n - 1) as u64) as f64 * dt;
        for k in 0..2 {
            let partial = estimate.partial(k).unwrap();
            close(partial.frequency_hz, first[k] + drift[k] * t, 2e-8);
            close(
                partial.end_frequency_hz.unwrap(),
                first[k] + drift[k] * last_t,
                2e-8,
            );
            close(partial.decay_rate_per_s.unwrap(), decay[k], 1e-7);
            let rotated = translated_chirp_phasor(t, first[k], drift[k], phasors[k]);
            let envelope = (-decay[k] * t).exp();
            close(partial.peak_phasor.re, envelope * rotated.re, 2e-8);
            close(partial.peak_phasor.im, envelope * rotated.im, 2e-8);
        }
    }
}

#[test]
fn damped_chirp_frames_refine_changed_decay_after_bounded_endpoint_prediction() {
    let n = 257;
    let dt = 1.0 / 2048.0;
    let phasor = C64::new(0.8, -0.3);
    let transition = n as f64 * dt;
    let rotated = translated_chirp_phasor(transition, 170.0, 72.0, phasor);
    let envelope = (-1.1 * transition).exp();
    let second_phasor = C64::new(envelope * rotated.re, envelope * rotated.im);
    let mut record: Vec<_> = (0..n)
        .map(|j| {
            let t = j as f64 * dt;
            0.1 + (-1.1 * t).exp() * chirp_sample(t, 170.0, 72.0, phasor)
        })
        .collect();
    // Continuous amplitude and phase at the transition, but a new local
    // drift and decay. The first law predicts an endpoint outside the fixed
    // domain, so the guess is bounded before all three coordinates are refit.
    record.extend((0..2 * n).map(|j| {
        let t = j as f64 * dt;
        0.1 + (-4.6 * t).exp() * chirp_sample(t, 179.0, 4.0, second_phasor)
    }));
    let mut frames = Vec::new();
    let summary = track_partials(
        &record,
        dt,
        PartialFrameSchedule {
            frame_samples: n,
            hop_samples: n,
            max_frames: 3,
            first_sample_index: 19,
        },
        PartialTrackingModel::DampedChirped {
            initial_tones: &[DampedChirpedTone {
                start_frequency_hz: 170.0,
                end_frequency_hz: 179.0,
                decay_rate_per_s: 1.1,
            }],
            bounds: &[DampedChirpedToneBounds {
                start_frequency_hz: [169.0, 190.0],
                end_frequency_hz: [170.0, 180.7],
                decay_rate_per_s: [0.5, 6.0],
            }],
            budget: damped_chirp_budget(),
        },
        |frame| {
            frames.push(frame.clone());
            ControlFlow::Continue(())
        },
    )
    .unwrap();
    assert_eq!(summary.resolved_frames, 3);
    for (j, frame) in frames[1..].iter().enumerate() {
        let estimate = frame.estimate.as_ref().unwrap();
        let PartialFrameEstimate::DampedChirped(result) = estimate else {
            panic!("combined model");
        };
        if j == 0 {
            assert!(result.evaluations > 1, "new decay must be measured");
        } else {
            assert_eq!(
                result.evaluations, 1,
                "accepted decay must warm-start the next frame"
            );
        }
        let t = (j * n) as f64 * dt;
        let partial = estimate.partial(0).unwrap();
        close(partial.frequency_hz, 179.0 + 4.0 * t, 2e-8);
        close(partial.end_frequency_hz.unwrap(), 179.5 + 4.0 * t, 2e-8);
        close(partial.decay_rate_per_s.unwrap(), 4.6, 1e-7);
        let rotated = translated_chirp_phasor(t, 179.0, 4.0, second_phasor);
        let envelope = (-4.6 * t).exp();
        close(partial.peak_phasor.re, envelope * rotated.re, 2e-8);
        close(partial.peak_phasor.im, envelope * rotated.im, 2e-8);
        assert!(estimate.fit().residual_rms < 3e-10);
    }
}

#[test]
fn damped_chirp_dropout_retains_decay_and_the_last_successful_clock_anchor() {
    let n = 513;
    let dt = 1.0 / 4096.0;
    let phasor = C64::new(0.7, 0.2);
    let mut record: Vec<_> = (0..3 * n)
        .map(|j| {
            let t = j as f64 * dt;
            0.3 + (-2.3 * t).exp() * chirp_sample(t, 213.2, 20.0, phasor)
        })
        .collect();
    record[n..2 * n].fill(3.0);
    let mut frames = Vec::new();
    let summary = track_partials(
        &record,
        dt,
        PartialFrameSchedule {
            frame_samples: n,
            hop_samples: n,
            max_frames: 3,
            first_sample_index: 11,
        },
        PartialTrackingModel::DampedChirped {
            initial_tones: &[DampedChirpedTone {
                start_frequency_hz: 213.2,
                end_frequency_hz: 215.7,
                decay_rate_per_s: 2.3,
            }],
            bounds: &[DampedChirpedToneBounds {
                start_frequency_hz: [205.0, 225.0],
                end_frequency_hz: [205.0, 225.0],
                decay_rate_per_s: [1.0, 4.0],
            }],
            budget: DampedChirpedToneRefinementBudget {
                max_evaluations: 1,
                ..damped_chirp_budget()
            },
        },
        |frame| {
            frames.push(frame.clone());
            ControlFlow::Continue(())
        },
    )
    .unwrap();
    assert!(summary.completed);
    assert_eq!(summary.resolved_frames, 2);
    assert!(matches!(
        frames[1].estimate,
        Err(PartialFitFailure::DampedChirped(
            DampedChirpedToneRefinementError::Unidentifiable { .. }
        ))
    ));
    let last = frames[2].estimate.as_ref().unwrap().partial(0).unwrap();
    let t = (2 * n) as f64 * dt;
    close(last.frequency_hz, 213.2 + 20.0 * t, 2e-8);
    close(last.end_frequency_hz.unwrap(), 215.7 + 20.0 * t, 2e-8);
    close(last.decay_rate_per_s.unwrap(), 2.3, 1e-7);
    let rotated = translated_chirp_phasor(t, 213.2, 20.0, phasor);
    let envelope = (-2.3 * t).exp();
    close(last.peak_phasor.re, envelope * rotated.re, 2e-8);
    close(last.peak_phasor.im, envelope * rotated.im, 2e-8);
}

#[test]
fn damped_chirp_tracker_admits_storage_and_retains_typed_refusals_and_stop() {
    let n = 257;
    let dt = 1.0 / 2048.0;
    let record: Vec<_> = (0..2 * n)
        .map(|j| {
            let t = j as f64 * dt;
            (-1.1 * t).exp() * chirp_sample(t, 170.0, 4.0, C64::new(0.8, -0.3))
        })
        .collect();
    let initial = [DampedChirpedTone {
        start_frequency_hz: 170.0,
        end_frequency_hz: 170.5,
        decay_rate_per_s: 1.1,
    }];
    let bounds = [DampedChirpedToneBounds {
        start_frequency_hz: [165.0, 180.0],
        end_frequency_hz: [165.0, 180.0],
        decay_rate_per_s: [0.5, 3.0],
    }];
    let schedule = PartialFrameSchedule {
        frame_samples: n,
        hop_samples: n,
        max_frames: 2,
        first_sample_index: 7,
    };
    let required_bytes = 2 * core::mem::size_of::<DampedChirpedTone>();
    let max_bytes = required_bytes - 1;
    assert_eq!(
        track_partials(
            &record,
            dt,
            schedule,
            PartialTrackingModel::DampedChirped {
                initial_tones: &initial,
                bounds: &bounds,
                budget: DampedChirpedToneRefinementBudget {
                    fit: HarmonicFitBudget {
                        max_workspace_bytes: max_bytes,
                        ..damped_chirp_budget().fit
                    },
                    ..damped_chirp_budget()
                },
            },
            |_| panic!("both three-coordinate arrays must be admitted first")
        ),
        Err(PartialTrackingError::Fit {
            first_sample_index: 7,
            failure: PartialFitFailure::DampedChirped(DampedChirpedToneRefinementError::Fit(
                HarmonicFitError::BudgetExceeded {
                    required_bytes,
                    max_bytes
                }
            )),
        })
    );
    assert!(matches!(
        track_partials(
            &record,
            dt,
            schedule,
            PartialTrackingModel::DampedChirped {
                initial_tones: &[DampedChirpedTone {
                    decay_rate_per_s: 4.0,
                    ..initial[0]
                }],
                bounds: &bounds,
                budget: damped_chirp_budget(),
            },
            |_| panic!("initial decay must reach admission without clamping")
        ),
        Err(PartialTrackingError::Fit {
            failure: PartialFitFailure::DampedChirped(
                DampedChirpedToneRefinementError::InvalidParameter { .. }
            ),
            ..
        })
    ));
    let summary = track_partials(
        &record,
        dt,
        schedule,
        PartialTrackingModel::DampedChirped {
            initial_tones: &initial,
            bounds: &bounds,
            budget: damped_chirp_budget(),
        },
        |frame| {
            assert!(frame.estimate.is_ok());
            ControlFlow::Break(())
        },
    )
    .unwrap();
    assert_eq!(summary.visited_frames, 1);
    assert_eq!(summary.resolved_frames, 1);
    assert_eq!(summary.maximum_candidate_evaluations, 160);
    assert!(!summary.completed);
    let summary = track_partials(
        &record[..n],
        dt,
        schedule,
        PartialTrackingModel::DampedChirped {
            initial_tones: &[DampedChirpedTone {
                decay_rate_per_s: 1.6,
                ..initial[0]
            }],
            bounds: &bounds,
            budget: DampedChirpedToneRefinementBudget {
                max_evaluations: 1,
                ..damped_chirp_budget()
            },
        },
        |frame| {
            assert!(matches!(
                frame.estimate,
                Err(PartialFitFailure::DampedChirped(
                    DampedChirpedToneRefinementError::NotConverged { .. }
                ))
            ));
            ControlFlow::Continue(())
        },
    )
    .unwrap();
    assert!(summary.completed);
    assert_eq!(summary.resolved_frames, 0);
}
