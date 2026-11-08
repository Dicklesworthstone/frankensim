use core::{f64::consts::TAU, mem::size_of, ops::ControlFlow};
use fs_fft::C64;
use fs_fft::harmonic::{
    ChirpedTone, ChirpedToneBounds, ChirpedToneRefinementBudget, ChirpedToneRefinementError,
    DampedChirpedTone, DampedChirpedToneBounds, DampedChirpedToneRefinementBudget,
    DampedChirpedToneRefinementError, DampedToneBounds, DampedTonePole, DampedToneRefinementBudget,
    HarmonicFitBudget, HarmonicFitError, ToneFrequencyBudget, ToneFrequencyError,
};
use fs_psycho::partials::{
    PartialChannelsFrame, PartialChannelsFrameEstimate, PartialFitFailure, PartialFrameSchedule,
    PartialTrackingError, PartialTrackingModel, track_partial_channels, track_partials,
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

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn shifted(phasor: C64, cycles: f64, scale: f64) -> C64 {
    let (s, c) = (TAU * cycles).sin_cos();
    C64::new(
        scale * (phasor.re * c + phasor.im * s),
        scale * (phasor.im * c - phasor.re * s),
    )
}

fn observe(
    records: &[Vec<f64>],
    dt: f64,
    schedule: PartialFrameSchedule,
    model: PartialTrackingModel<'_>,
) -> Vec<PartialChannelsFrame> {
    let channels: Vec<_> = records.iter().map(Vec::as_slice).collect();
    let mut frames = Vec::new();
    let summary = track_partial_channels(&channels, dt, schedule, model, |frame| {
        frames.push(frame.clone());
        ControlFlow::Continue(())
    })
    .unwrap();
    assert!(summary.completed);
    assert_eq!(summary.visited_frames, summary.planned_frames);
    frames
}

const FREQUENCIES: [f64; 2] = [131.37, 371.83];
const PHASORS: [C64; 2] = [C64::new(0.8, 0.3), C64::new(0.4, -0.25)];
const SHAPES: [[f64; 3]; 2] = [[0.0, 1.0, -1.0], [1.0, 0.0, 0.3]];
const OFFSETS: [f64; 3] = [-2.0, 1.0, 0.25];

fn modal_records(count: usize, dt: f64) -> Vec<Vec<f64>> {
    (0..3)
        .map(|c| {
            (0..count)
                .map(|i| {
                    let t = i as f64 * dt;
                    OFFSETS[c]
                        + (0..2)
                            .map(|k| {
                                let (s, co) = (TAU * FREQUENCIES[k] * t).sin_cos();
                                SHAPES[k][c] * (PHASORS[k].re * co + PHASORS[k].im * s)
                            })
                            .sum::<f64>()
                })
                .collect()
        })
        .collect()
}

fn stationary_model() -> PartialTrackingModel<'static> {
    PartialTrackingModel::Stationary {
        initial_frequencies_hz: &[131.6, 371.6],
        frequency_bounds_hz: &[[130.0, 133.0], [370.0, 374.0]],
        budget: tone_budget(),
    }
}

#[test]
fn array_recovers_modes_at_sensor_nodes_without_antiphase_cancellation() {
    let dt = 1.0 / 4096.0;
    let records = modal_records(1537, dt);
    let schedule = PartialFrameSchedule {
        frame_samples: 769,
        hop_samples: 512,
        max_frames: 3,
        first_sample_index: (1_u64 << 53) + 71,
    };
    let frames = observe(&records, dt, schedule, stationary_model());
    let starts = [0, 512, 768];
    for (frame, offset) in frames.iter().zip(starts) {
        assert_eq!(
            frame.first_sample_index,
            schedule.first_sample_index + offset
        );
        let estimate = frame.estimate.as_ref().unwrap();
        assert_eq!(estimate.fits().len(), 3);
        assert!(estimate.partial(3, 0).is_none());
        assert!(estimate.partial(0, 2).is_none());
        let t = offset as f64 * dt;
        for c in 0..3 {
            close(estimate.fits()[c].dc, OFFSETS[c], 2e-10);
            assert!(estimate.fits()[c].residual_rms < 2e-10);
            for k in 0..2 {
                let partial = estimate.partial(c, k).unwrap();
                close(partial.frequency_hz, FREQUENCIES[k], 2e-8);
                assert_eq!(partial.end_frequency_hz, None);
                assert_eq!(partial.decay_rate_per_s, None);
                let expected = shifted(PHASORS[k], FREQUENCIES[k] * t, SHAPES[k][c]);
                close(partial.peak_phasor.re, expected.re, 2e-8);
                close(partial.peak_phasor.im, expected.im, 2e-8);
            }
        }
    }
    // Sensor zero is at the first mode's node; sensors one and two observe it
    // with opposite sign. Either selecting sensor zero or averaging loses it.
    let summary = track_partials(&records[0], dt, schedule, stationary_model(), |frame| {
        assert!(
            frame.estimate.is_err(),
            "the nodal scalar observation must remain unresolved: {:?}",
            frame.estimate
        );
        ControlFlow::Continue(())
    })
    .unwrap();
    assert_eq!(summary.resolved_frames, 0);
    assert_eq!(frames, observe(&records, dt, schedule, stationary_model()));
}

#[test]
fn rotating_common_unit_components_preserves_joint_fit_and_complex_mode_shapes() {
    let dt = 1.0 / 4096.0;
    let mut records = modal_records(769, dt);
    // Unmodeled physical motion makes the objective nonzero; rotation must
    // preserve its channel weights as well as the exact-mode limit.
    for (c, record) in records.iter_mut().enumerate() {
        for (i, sample) in record.iter_mut().enumerate() {
            *sample += [0.017, -0.023, 0.009][c] * (TAU * 517.6 * i as f64 * dt + 0.4).cos();
        }
    }
    let q = core::f64::consts::FRAC_1_SQRT_2;
    let rotation = [[q, q, 0.0], [-q, q, 0.0], [0.0, 0.0, 1.0]];
    let rotated: Vec<Vec<_>> = rotation
        .iter()
        .map(|row| {
            (0..769)
                .map(|i| (0..3).map(|c| row[c] * records[c][i]).sum())
                .collect()
        })
        .collect();
    let schedule = PartialFrameSchedule {
        frame_samples: 769,
        hop_samples: 769,
        max_frames: 1,
        first_sample_index: 19,
    };
    let original = observe(&records, dt, schedule, stationary_model());
    let transformed = observe(&rotated, dt, schedule, stationary_model());
    let a = original[0].estimate.as_ref().unwrap();
    let b = transformed[0].estimate.as_ref().unwrap();
    for k in 0..2 {
        close(
            a.partial(0, k).unwrap().frequency_hz,
            b.partial(0, k).unwrap().frequency_hz,
            2e-8,
        );
    }
    for c in 0..3 {
        let dc: f64 = (0..3).map(|j| rotation[c][j] * a.fits()[j].dc).sum();
        close(b.fits()[c].dc, dc, 2e-8);
        for k in 0..2 {
            let re: f64 = (0..3)
                .map(|j| rotation[c][j] * a.partial(j, k).unwrap().peak_phasor.re)
                .sum();
            let im: f64 = (0..3)
                .map(|j| rotation[c][j] * a.partial(j, k).unwrap().peak_phasor.im)
                .sum();
            close(b.partial(c, k).unwrap().peak_phasor.re, re, 2e-8);
            close(b.partial(c, k).unwrap().peak_phasor.im, im, 2e-8);
        }
    }
    let energy_a: f64 = a.fits().iter().map(|f| f.residual_rms.powi(2)).sum();
    let energy_b: f64 = b.fits().iter().map(|f| f.residual_rms.powi(2)).sum();
    assert!(energy_a > 1e-5);
    close(energy_a, energy_b, 2e-12);
}

#[test]
fn spatial_array_tracks_equal_frequency_double_decay_with_separate_sensor_phases() {
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
    let shapes = [[1.0, 0.0, -1.0], [0.0, 1.0, 0.25]];
    let records: Vec<Vec<_>> = (0..3)
        .map(|c| {
            (0..1299)
                .map(|i| {
                    let t = i as f64 * dt;
                    OFFSETS[c]
                        + (0..2)
                            .map(|k| {
                                let (s, co) = (TAU * poles[k].frequency_hz * t).sin_cos();
                                shapes[k][c]
                                    * (-poles[k].decay_rate_per_s * t).exp()
                                    * (PHASORS[k].re * co + PHASORS[k].im * s)
                            })
                            .sum::<f64>()
                })
                .collect()
        })
        .collect();
    let initial = poles.map(|p| DampedTonePole {
        frequency_hz: p.frequency_hz + 0.02,
        decay_rate_per_s: p.decay_rate_per_s + 0.1,
    });
    let bounds = [
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
    let frames = observe(
        &records,
        dt,
        schedule,
        PartialTrackingModel::Ringdown {
            initial_poles: &initial,
            bounds: &bounds,
            budget: DampedToneRefinementBudget {
                absolute_frequency_tolerance_hz: 1e-8,
                absolute_decay_tolerance_per_s: 1e-7,
                max_evaluations: 80,
                min_pole_reciprocal_condition: 1e-10,
                fit: tone_budget().fit,
            },
        },
    );
    assert_eq!(frames.len(), 3);
    for frame in frames {
        let t = (frame.first_sample_index - 127) as f64 * dt;
        let result = frame.estimate.unwrap();
        for c in 0..3 {
            close(result.fits()[c].dc, OFFSETS[c], 1e-9);
            assert!(result.fits()[c].residual_rms < 1e-9);
            for k in 0..2 {
                let partial = result.partial(c, k).unwrap();
                close(partial.frequency_hz, poles[k].frequency_hz, 1e-7);
                close(
                    partial.decay_rate_per_s.unwrap(),
                    poles[k].decay_rate_per_s,
                    2e-6,
                );
                let expected = shifted(
                    PHASORS[k],
                    poles[k].frequency_hz * t,
                    shapes[k][c] * (-poles[k].decay_rate_per_s * t).exp(),
                );
                close(partial.peak_phasor.re, expected.re, 5e-7);
                close(partial.peak_phasor.im, expected.im, 5e-7);
            }
        }
    }
}

#[test]
fn shared_chirp_keeps_successful_anchor_across_dropout_and_short_terminal_hop() {
    let n = 257;
    let dt = 1.0 / 2048.0;
    let first_hz = 173.37;
    let drift = 20.0;
    let phasors = [C64::new(0.7, 0.2), C64::new(-0.7, -0.2), C64::new(0.0, 0.0)];
    let mut records: Vec<Vec<_>> = (0..3)
        .map(|c| {
            (0..3 * n + 73)
                .map(|i| {
                    let t = i as f64 * dt;
                    OFFSETS[c] + shifted(phasors[c], first_hz * t + 0.5 * drift * t * t, 1.0).re
                })
                .collect()
        })
        .collect();
    for (c, record) in records.iter_mut().enumerate() {
        record[n..2 * n].fill(OFFSETS[c]);
    }
    let initial = [ChirpedTone {
        start_frequency_hz: first_hz,
        end_frequency_hz: first_hz + drift * (n - 1) as f64 * dt,
    }];
    let bounds = [ChirpedToneBounds {
        start_frequency_hz: [170.0, 185.0],
        end_frequency_hz: [170.0, 185.0],
    }];
    let schedule = PartialFrameSchedule {
        frame_samples: n,
        hop_samples: n,
        max_frames: 4,
        first_sample_index: (1_u64 << 53) + 11,
    };
    let model = PartialTrackingModel::Chirped {
        initial_tones: &initial,
        bounds: &bounds,
        budget: ChirpedToneRefinementBudget {
            absolute_frequency_tolerance_hz: 1e-8,
            // Exact seeds plus the shared predictor must suffice at EVERY gap.
            max_evaluations: 1,
            min_chirp_reciprocal_condition: 1e-9,
            fit: tone_budget().fit,
        },
    };
    let frames = observe(&records, dt, schedule, model);
    assert_eq!(frames.len(), 4);
    assert!(matches!(
        frames[1].estimate,
        Err(PartialFitFailure::Chirped(
            ChirpedToneRefinementError::Unidentifiable { .. }
        ))
    ));
    for (frame, offset) in frames.iter().zip([0, n, 2 * n, 2 * n + 73]) {
        assert_eq!(
            frame.first_sample_index,
            schedule.first_sample_index + offset as u64
        );
        if offset == n {
            continue;
        }
        let result = frame.estimate.as_ref().unwrap();
        let PartialChannelsFrameEstimate::Chirped(fit) = result else {
            panic!("chirp diagnostics must be retained")
        };
        assert_eq!(fit.evaluations, 1);
        let t = offset as f64 * dt;
        for c in 0..3 {
            let partial = result.partial(c, 0).unwrap();
            close(partial.frequency_hz, first_hz + drift * t, 2e-8);
            close(
                partial.end_frequency_hz.unwrap(),
                first_hz + drift * (offset + n - 1) as f64 * dt,
                2e-8,
            );
            assert_eq!(partial.decay_rate_per_s, None);
            let expected = shifted(phasors[c], first_hz * t + 0.5 * drift * t * t, 1.0);
            close(partial.peak_phasor.re, expected.re, 2e-8);
            close(partial.peak_phasor.im, expected.im, 2e-8);
            close(result.fits()[c].dc, OFFSETS[c], 2e-10);
            assert!(result.fits()[c].residual_rms < 3e-10);
        }
    }
    assert_eq!(frames, observe(&records, dt, schedule, model));
}

#[test]
fn channels_are_admitted_together_and_their_views_are_budgeted_before_observation() {
    let records = modal_records(1537, 1.0 / 4096.0);
    let schedule = PartialFrameSchedule {
        frame_samples: 769,
        hop_samples: 512,
        max_frames: 3,
        first_sample_index: 17,
    };
    for channels in [vec![], vec![&records[0][..], &records[1][..100]]] {
        assert!(matches!(
            track_partial_channels(
                &channels,
                1.0 / 4096.0,
                schedule,
                stationary_model(),
                |_| panic!("invalid channel shapes must precede callbacks")
            ),
            Err(PartialTrackingError::InvalidParameter { .. })
        ));
    }
    let mut bad = records[2].clone();
    *bad.last_mut().unwrap() = f64::NAN;
    assert!(matches!(
        track_partial_channels(
            &[&records[0], &bad],
            1.0 / 4096.0,
            schedule,
            stationary_model(),
            |_| panic!("all input samples admitted together")
        ),
        Err(PartialTrackingError::InvalidParameter { .. })
    ));
    let required_bytes = 2 * size_of::<f64>() + 3 * size_of::<&[f64]>();
    let limit = required_bytes - 1;
    let small = PartialTrackingModel::Stationary {
        initial_frequencies_hz: &[131.6, 371.6],
        frequency_bounds_hz: &[[130.0, 133.0], [370.0, 374.0]],
        budget: ToneFrequencyBudget {
            fit: HarmonicFitBudget {
                max_workspace_bytes: limit,
                ..tone_budget().fit
            },
            ..tone_budget()
        },
    };
    let channels: Vec<_> = records.iter().map(Vec::as_slice).collect();
    assert_eq!(
        track_partial_channels(&channels, 1.0 / 4096.0, schedule, small, |_| panic!(
            "seed and view storage must be admitted first"
        )),
        Err(PartialTrackingError::Fit {
            first_sample_index: 17,
            failure: PartialFitFailure::Stationary(ToneFrequencyError::Fit(
                HarmonicFitError::BudgetExceeded {
                    required_bytes,
                    max_bytes: limit
                }
            )),
        })
    );
    let mut called = 0;
    let stopped = track_partial_channels(
        &channels,
        1.0 / 4096.0,
        schedule,
        stationary_model(),
        |frame| {
            called += 1;
            assert!(frame.estimate.is_ok());
            ControlFlow::Break(())
        },
    )
    .unwrap();
    assert_eq!(called, 1);
    assert_eq!(stopped.visited_frames, 1);
    assert_eq!(stopped.resolved_frames, 1);
    assert!(!stopped.completed);
    assert_eq!(stopped.maximum_candidate_evaluations, 240);
}
#[test]
fn shared_damped_chirp_refits_decaying_channel_phases_across_dropout_and_short_hop() {
    let n = 257;
    let dt = 1.0 / 2048.0;
    let first_hz = 173.37;
    let drift = 20.0;
    let decay: f64 = 2.3;
    let phasors = [C64::new(0.7, 0.2), C64::new(-0.7, -0.2), C64::new(0.0, 0.0)];
    let mut records: Vec<Vec<_>> = (0..3)
        .map(|c| {
            (0..3 * n + 73)
                .map(|i| {
                    let t = i as f64 * dt;
                    OFFSETS[c]
                        + shifted(
                            phasors[c],
                            first_hz * t + 0.5 * drift * t * t,
                            (-decay * t).exp(),
                        )
                        .re
                })
                .collect()
        })
        .collect();
    for (c, record) in records.iter_mut().enumerate() {
        record[n..2 * n].fill(OFFSETS[c]);
    }
    let initial = [DampedChirpedTone {
        start_frequency_hz: first_hz,
        end_frequency_hz: first_hz + drift * (n - 1) as f64 * dt,
        decay_rate_per_s: decay,
    }];
    let bounds = [DampedChirpedToneBounds {
        start_frequency_hz: [170.0, 185.0],
        end_frequency_hz: [170.0, 185.0],
        decay_rate_per_s: [1.0, 4.0],
    }];
    let schedule = PartialFrameSchedule {
        frame_samples: n,
        hop_samples: n,
        max_frames: 4,
        first_sample_index: (1_u64 << 53) + 11,
    };
    let model = PartialTrackingModel::DampedChirped {
        initial_tones: &initial,
        bounds: &bounds,
        budget: DampedChirpedToneRefinementBudget {
            absolute_frequency_tolerance_hz: 1e-8,
            absolute_decay_tolerance_per_s: 1e-8,
            // Exact seeds plus the shared predictor must suffice at EVERY gap.
            max_evaluations: 1,
            min_trajectory_reciprocal_condition: 1e-9,
            fit: tone_budget().fit,
        },
    };
    let frames = observe(&records, dt, schedule, model);
    assert_eq!(frames.len(), 4);
    assert!(matches!(
        frames[1].estimate,
        Err(PartialFitFailure::DampedChirped(
            DampedChirpedToneRefinementError::Unidentifiable { .. }
        ))
    ));
    for (frame, offset) in frames.iter().zip([0, n, 2 * n, 2 * n + 73]) {
        assert_eq!(
            frame.first_sample_index,
            schedule.first_sample_index + offset as u64
        );
        if offset == n {
            continue;
        }
        let result = frame.estimate.as_ref().unwrap();
        let PartialChannelsFrameEstimate::DampedChirped(fit) = result else {
            panic!("chirp diagnostics must be retained")
        };
        assert_eq!(fit.evaluations, 1);
        let t = offset as f64 * dt;
        for c in 0..3 {
            let partial = result.partial(c, 0).unwrap();
            close(partial.frequency_hz, first_hz + drift * t, 2e-8);
            close(
                partial.end_frequency_hz.unwrap(),
                first_hz + drift * (offset + n - 1) as f64 * dt,
                2e-8,
            );
            close(partial.decay_rate_per_s.unwrap(), decay, 2e-8);
            let expected = shifted(
                phasors[c],
                first_hz * t + 0.5 * drift * t * t,
                (-decay * t).exp(),
            );
            close(partial.peak_phasor.re, expected.re, 2e-8);
            close(partial.peak_phasor.im, expected.im, 2e-8);
            close(result.fits()[c].dc, OFFSETS[c], 2e-10);
            assert!(result.fits()[c].residual_rms < 3e-10);
        }
    }
    assert_eq!(frames, observe(&records, dt, schedule, model));
}
