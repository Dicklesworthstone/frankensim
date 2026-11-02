//! G1/G3 accepted physical histories, retarded clocks and the original PCM path.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, TimedCompactSource,
    UniformSourceMotion,
};
use fs_couple::pcm_wav::decimate::Decimator;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16};
use fs_couple::render::RenderError;
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;
use std::io::Cursor;

fn observer(position: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        position,
        MovingCompactLimits {
            maximum_mach: 0.9,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap()
}

fn budget() -> MovingPressureBudget {
    MovingPressureBudget {
        maximum_sources: 8,
        maximum_source_samples: 2048,
        maximum_output_samples: 100_000,
        maximum_payload_bytes: 1_000_000,
        maximum_source_receiver_pairs: 1_000_000,
    }
}

fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.05,
        sample_rate_hz: 1000,
        samples: 81,
        max_block: 81,
        receiver_bandwidth_hz: 10.0,
        maximum_time_residual_s: 1e-12,
    }
}

// Cubics are represented exactly by endpoint values and derivatives on any
// admitted knot partition. The oracle evaluates this original polynomial.
fn strength(t: f64, scale: f64) -> CompactSourceStrength {
    CompactSourceStrength {
        volume_flow_m3_s: scale * (0.001 + 0.002 * t - 0.003 * t * t + 0.004 * t * t * t),
        volume_flow_rate_m3_s2: scale * (0.002 - 0.006 * t + 0.012 * t * t),
        force_on_fluid_n: [
            scale * (1.0 + 2.0 * t - 3.0 * t * t + 4.0 * t * t * t),
            scale * (-2.0 + 3.0 * t + t * t - 2.0 * t * t * t),
            scale * (0.5 - t + 2.0 * t * t + 3.0 * t * t * t),
        ],
        force_rate_on_fluid_n_s: [
            scale * (2.0 - 6.0 * t + 12.0 * t * t),
            scale * (3.0 + 2.0 * t - 6.0 * t * t),
            scale * (-1.0 + 4.0 * t + 9.0 * t * t),
        ],
    }
}

fn history(times: &[f64], shift: f64, scale: f64) -> Vec<TimedCompactSource> {
    times
        .iter()
        .map(|&source_time_s| TimedCompactSource {
            source_time_s,
            strength: strength(source_time_s + shift, scale),
        })
        .collect()
}

fn source(
    samples: &[TimedCompactSource],
    position: [f64; 3],
    velocity: [f64; 3],
) -> CompactSourceHistory<'_> {
    CompactSourceHistory {
        motion: UniformSourceMotion::new(0.0, position, velocity).unwrap(),
        samples,
        source_bandwidth_hz: 1.0,
    }
}

fn build<'a>(
    sources: &[CompactSourceHistory<'a>],
    config: MovingPressureConfig,
) -> MovingCompactRenderer<'a> {
    MovingCompactRenderer::new(
        observer([12.0, 0.0, 0.0]),
        sources,
        config,
        budget(),
        &mut || false,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-12 * expected.abs().max(1e-10),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn cubic_source_and_its_derivative_follow_stationary_and_moving_receiver_clocks() {
    let accepted = history(&[-0.25, 0.0, 0.25, 0.5], 0.0, 1.0);
    for velocity in [-60.0, 0.0, 60.0] {
        let sources = [source(&accepted, [0.0; 3], [velocity, 0.0, 0.0])];
        let mut renderer = build(&sources, config());
        let mut actual = [0.0; 81];
        renderer.block(&mut actual).unwrap();
        for (i, pressure) in actual.into_iter().enumerate() {
            let time = renderer.receiver_time_at(i as u64).unwrap();
            let kappa = 1.0 - velocity / 300.0;
            let tau = (time - 12.0 / 300.0) / kappa;
            let radius = 12.0 - velocity * tau;
            let s = strength(tau, 1.0);
            // Independent one-dimensional reduction of the moving Green law.
            // Constant Q also contributes when its source is translating.
            let monopole = 1.2
                * (s.volume_flow_rate_m3_s2 / radius
                    + s.volume_flow_m3_s * velocity / (radius * radius));
            let force = s.force_on_fluid_n[0] / (radius * radius)
                + s.force_rate_on_fluid_n_s[0] / (300.0 * radius);
            close(pressure, (monopole + force) / (4.0 * PI * kappa * kappa));
        }
        assert_eq!(renderer.samples_rendered(), 81);
    }
}

#[test]
fn different_source_grids_sum_at_the_same_arrival_time_without_normalization() {
    let a = history(&[-0.25, 0.0, 0.25, 0.5], 0.0, 1.0);
    let b = history(&[-0.3, -0.1, 0.1, 0.4], 0.03, -0.7);
    let sources = [
        source(&a, [0.0; 3], [0.0; 3]),
        source(&b, [12.0, -9.0, 0.0], [0.0; 3]),
    ];
    let mut renderer = build(&sources, config());
    let mut actual = [0.0; 81];
    renderer.block(&mut actual).unwrap();
    for (i, pressure) in actual.into_iter().enumerate() {
        let t = renderer.receiver_time_at(i as u64).unwrap();
        let first = strength(t - 12.0 / 300.0, 1.0);
        let second = strength(t - 9.0 / 300.0 + 0.03, -0.7);
        let p = |s: CompactSourceStrength, radius: f64, axis: usize| {
            (1.2 * s.volume_flow_rate_m3_s2 / radius
                + s.force_on_fluid_n[axis] / (radius * radius)
                + s.force_rate_on_fluid_n_s[axis] / (300.0 * radius))
                / (4.0 * PI)
        };
        close(pressure, p(first, 12.0, 0) + p(second, 9.0, 1));
    }
}

#[test]
fn callback_partitions_preserve_pressure_bits_and_global_timestamps() {
    let accepted = history(&[-0.25, 0.0, 0.25, 0.5], 0.0, 1.0);
    let sources = [source(&accepted, [0.0; 3], [60.0, 0.0, 0.0])];
    let mut whole = build(&sources, config());
    let mut expected = [0.0; 81];
    whole.block(&mut expected).unwrap();
    for block in [1, 7, 37] {
        let mut split = build(&sources, config());
        let mut actual = [0.0; 81];
        for chunk in actual.chunks_mut(block) {
            split.block(chunk).unwrap();
        }
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        for i in 0..81 {
            assert_eq!(
                split.receiver_time_at(i).unwrap(),
                whole.receiver_time_at(i).unwrap()
            );
        }
        assert_eq!(split.samples_rendered(), 81);
        assert_eq!(split.remaining_samples(), 0);
    }
}

#[test]
fn exact_accepted_endpoints_work_but_missing_prehistory_and_tails_refuse() {
    let accepted = history(&[0.0, 0.5], 0.0, 1.0);
    let mut track = source(&accepted, [0.0; 3], [0.0; 3]);
    track.source_bandwidth_hz = 0.5;
    let sources = [track];
    let config = MovingPressureConfig {
        first_receiver_time_s: 0.25,
        sample_rate_hz: 4,
        samples: 3,
        max_block: 3,
        receiver_bandwidth_hz: 1.0,
        ..config()
    };
    let prepare = |config| {
        MovingCompactRenderer::new(
            observer([75.0, 0.0, 0.0]),
            &sources,
            config,
            budget(),
            &mut || false,
        )
    };
    let mut renderer = prepare(config).unwrap();
    let mut actual = [0.0; 3];
    renderer.block(&mut actual).unwrap();
    for (i, pressure) in actual.into_iter().enumerate() {
        let s = strength(i as f64 / 4.0, 1.0);
        close(
            pressure,
            (1.2 * s.volume_flow_rate_m3_s2 / 75.0
                + s.force_on_fluid_n[0] / (75.0 * 75.0)
                + s.force_rate_on_fluid_n_s[0] / (300.0 * 75.0))
                / (4.0 * PI),
        );
    }
    let mut untouched = [123.0];
    assert!(renderer.block(&mut untouched).is_err());
    assert_eq!(untouched, [123.0]);
    assert_eq!(renderer.samples_rendered(), 3);
    assert!(renderer.receiver_time_at(3).is_err());
    assert!(
        prepare(MovingPressureConfig {
            first_receiver_time_s: 0.25_f64.next_down(),
            ..config
        })
        .is_err()
    );
    assert!(
        prepare(MovingPressureConfig {
            samples: 4,
            ..config
        })
        .is_err()
    );

    // An exact retained emission is still admissible when the analytic inverse
    // accumulates a tiny residual that exceeds a zero residual budget.
    let microphone = MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, 10.0).unwrap(),
        [10.0, 0.0, 0.0],
        MovingCompactLimits {
            maximum_mach: 0.9,
            minimum_distance_m: 0.01,
        },
    )
    .unwrap();
    let motion = UniformSourceMotion::new(0.0, [0.0; 3], [3.0, 0.0, 0.0]).unwrap();
    let constant = CompactSourceStrength {
        volume_flow_m3_s: 0.001,
        force_on_fluid_n: [1.0, 0.0, 0.0],
        ..CompactSourceStrength::default()
    };
    let accepted = [0.0, 0.1].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: constant,
    });
    let track = CompactSourceHistory {
        motion,
        samples: &accepted,
        source_bandwidth_hz: 0.0,
    };
    let config = MovingPressureConfig {
        first_receiver_time_s: 1.0,
        sample_rate_hz: 10,
        samples: 1,
        max_block: 1,
        receiver_bandwidth_hz: 1.0,
        maximum_time_residual_s: 0.0,
    };
    let expected = microphone
        .prepare_emission(0.0, motion.kinematics_at(0.0).unwrap())
        .unwrap()
        .evaluate(constant)
        .unwrap();
    assert_eq!(expected.receiver_time_s, 1.0);
    assert_eq!(expected.retarded_time_residual_s, 0.0);
    let mut renderer =
        MovingCompactRenderer::new(microphone, &[track], config, budget(), &mut || false).unwrap();
    let mut actual = [0.0];
    renderer.block(&mut actual).unwrap();
    assert_eq!(actual[0], expected.pressure_pa);
}

#[test]
fn a_close_pass_between_receiver_samples_is_admitted_continuously() {
    let accepted = [-0.1, 0.1].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: CompactSourceStrength {
            volume_flow_m3_s: 0.001,
            ..CompactSourceStrength::default()
        },
    });
    let attempt = |impact: f64| {
        let microphone = observer([0.0; 3]);
        let motion = UniformSourceMotion::new(0.0, [0.0, impact, 0.0], [100.0, 0.0, 0.0]).unwrap();
        let first = microphone
            .prepare_emission(-0.01, motion.kinematics_at(-0.01).unwrap())
            .unwrap();
        let config = MovingPressureConfig {
            first_receiver_time_s: first.receiver_time_s(),
            sample_rate_hz: 50,
            samples: 2,
            max_block: 2,
            ..config()
        };
        // Both requested samples clear the microphone. The closest point lies
        // strictly between their emissions, so endpoint-only admission is wrong.
        for t in [
            config.first_receiver_time_s,
            config.first_receiver_time_s + 0.02,
        ] {
            assert!(
                microphone
                    .prepare_uniform(t, &motion, 1e-12)
                    .unwrap()
                    .distance_m()
                    > 1.0
            );
        }
        let track = CompactSourceHistory {
            motion,
            samples: &accepted,
            source_bandwidth_hz: 0.0,
        };
        MovingCompactRenderer::new(microphone, &[track], config, budget(), &mut || false).is_ok()
    };
    assert!(!attempt(0.005));
    assert!(attempt(0.02));
}

#[test]
fn sub_ulp_closest_pass_respects_clearance_at_a_large_physical_time_origin() {
    let epoch = 1e12;
    let microphone = MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        [0.0; 3],
        MovingCompactLimits {
            maximum_mach: 0.9,
            minimum_distance_m: 0.10001,
        },
    )
    .unwrap();
    let motion = UniformSourceMotion::new(epoch, [-0.003, 0.1, 0.0], [100.0, 0.0, 0.0]).unwrap();
    let accepted = [epoch - 1.0, epoch + 1.0].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: CompactSourceStrength {
            volume_flow_m3_s: 0.001,
            ..CompactSourceStrength::default()
        },
    });
    let config = MovingPressureConfig {
        first_receiver_time_s: epoch - 0.05,
        sample_rate_hz: 10,
        samples: 2,
        max_block: 2,
        receiver_bandwidth_hz: 1.0,
        maximum_time_residual_s: 0.005,
    };
    // The actual closest event is epoch + 3e-5, with 0.1 m separation.
    // Rounding that event first incorrectly evaluates the trajectory at epoch
    // and reports sqrt(0.1^2 + 0.003^2) > 0.10001 m, a false admission.
    assert_eq!(epoch + 3e-5, epoch);
    for t in [
        config.first_receiver_time_s,
        config.first_receiver_time_s + 0.1,
    ] {
        assert!(
            microphone
                .prepare_uniform(t, &motion, 0.005)
                .unwrap()
                .distance_m()
                > 1.0
        );
    }
    let track = CompactSourceHistory {
        motion,
        samples: &accepted,
        source_bandwidth_hz: 0.0,
    };
    assert!(
        MovingCompactRenderer::new(microphone, &[track], config, budget(), &mut || false).is_err()
    );
}

#[test]
fn bandwidth_and_complete_window_work_are_admitted_before_output() {
    let accepted = history(&[-0.25, 0.0, 0.25, 0.5], 0.0, 1.0);
    let track = source(&accepted, [0.0; 3], [30.0, 0.0, 0.0]);
    let prepare = |track, config, budget| {
        MovingCompactRenderer::new(
            observer([12.0, 0.0, 0.0]),
            &[track],
            config,
            budget,
            &mut || false,
        )
        .is_ok()
    };
    assert!(prepare(track, config(), budget()));
    // Recession compresses the admitted source band by 2/3, rather than the
    // global speed-only bound of 2. Keep both windows before any close pass so
    // the approaching refusal isolates its actual Doppler expansion.
    for (velocity, admitted) in [(-150.0, true), (150.0, false)] {
        assert_eq!(
            prepare(
                source(&accepted, [0.0; 3], [velocity, 0.0, 0.0]),
                MovingPressureConfig {
                    samples: 11,
                    receiver_bandwidth_hz: 0.8,
                    ..config()
                },
                budget(),
            ),
            admitted,
        );
    }
    assert!(!prepare(
        CompactSourceHistory {
            source_bandwidth_hz: 2.0,
            ..track
        },
        config(),
        budget()
    ));
    assert!(!prepare(
        track,
        MovingPressureConfig {
            receiver_bandwidth_hz: 1.05,
            ..config()
        },
        budget()
    ));
    assert!(!prepare(
        track,
        MovingPressureConfig {
            receiver_bandwidth_hz: 500.0,
            ..config()
        },
        budget()
    ));
    assert!(!prepare(
        track,
        config(),
        MovingPressureBudget {
            maximum_output_samples: 80,
            ..budget()
        }
    ));
    assert!(!prepare(
        track,
        config(),
        MovingPressureBudget {
            maximum_source_receiver_pairs: 80,
            ..budget()
        }
    ));
    assert!(!prepare(
        track,
        config(),
        MovingPressureBudget {
            maximum_source_samples: 3,
            ..budget()
        }
    ));
    assert!(!prepare(
        track,
        config(),
        MovingPressureBudget {
            maximum_payload_bytes: 1,
            ..budget()
        }
    ));
    let renderer = build(&[track], config());
    assert!(renderer.validate_sample_rate(999).is_err());
    assert!(renderer.validate_sample_count(82).is_err());
    assert_eq!(renderer.samples_rendered(), 0);
}

#[test]
fn cancellation_inside_a_block_preserves_output_and_counter_then_retries_exactly() {
    let accepted = history(&[-0.25, 0.0, 0.25, 0.5], 0.0, 1.0);
    let sources = [source(&accepted, [0.0; 3], [60.0, 0.0, 0.0])];
    let mut baseline = build(&sources, config());
    let mut expected = [0.0; 22];
    baseline.block(&mut expected).unwrap();
    let mut renderer = build(&sources, config());
    let mut prefix = [0.0; 5];
    renderer.block(&mut prefix).unwrap();
    let mut untouched = [123.0; 17];
    let mut polls = 0;
    assert!(matches!(
        renderer.block_controlled(&mut untouched, &mut || {
            polls += 1;
            polls == 8
        }),
        Err(RenderError::Acoustic(AeroacError::Cancelled { .. }))
    ));
    assert_eq!(untouched, [123.0; 17]);
    assert_eq!(renderer.samples_rendered(), 5);
    assert_eq!(renderer.remaining_samples(), 76);
    renderer.block(&mut untouched).unwrap();
    assert_eq!(prefix.as_slice(), &expected[..5]);
    assert_eq!(untouched.as_slice(), &expected[5..]);
    assert_eq!(renderer.samples_rendered(), 22);
}

#[test]
fn decimated_pcm_cancellation_and_resume_keep_the_same_physical_pressure_history() {
    let accepted = history(&[-0.25, 0.0, 0.25, 0.5], 0.0, 1.0);
    let sources = [source(&accepted, [0.0; 3], [60.0, 0.0, 0.0])];
    let config = MovingPressureConfig {
        sample_rate_hz: 96_000,
        samples: 514,
        max_block: 3,
        ..config()
    };
    let mut direct = build(&sources, config);
    let mut input = [0.0; 514];
    for chunk in input.chunks_mut(3) {
        direct.block(chunk).unwrap();
    }
    let mut filter = Decimator::new(2, 1).unwrap();
    let expected: Vec<_> = input
        .chunks_exact(2)
        .map(|group| {
            let p = filter.preview(group).unwrap()[0];
            filter.commit();
            p
        })
        .collect();
    let (wav, clips) = encode_pcm16_wav(&expected, 48_000, 0.01).unwrap();
    let mut renderer = DecimatedRenderer::new(build(&sources, config), 96_000, 48_000, 37).unwrap();
    assert_eq!(
        renderer.info().delay_output_samples,
        filter.delay_output_frames()
    );
    assert!(renderer.info().delay_output_samples > 0.0);
    assert_eq!(renderer.output_samples_for(514).unwrap(), 257);
    assert!(renderer.output_samples_for(513).is_err());
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 48_000, 0.01, 37).unwrap();
    let mut scratch = [0.0; 37];
    render_pressure_pcm16(
        &mut renderer,
        &mut stream,
        &CancelGate::new(),
        &mut scratch,
        37,
    )
    .unwrap();
    assert_eq!(renderer.source().samples_rendered(), 74);
    let bytes = stream.sink().get_ref().clone();
    let gate = CancelGate::new();
    gate.request();
    scratch.fill(123.0);
    assert_eq!(
        render_pressure_pcm16(&mut renderer, &mut stream, &gate, &mut scratch, 220).unwrap(),
        ScheduledWavProgress::Cancelled { samples: 0 }
    );
    assert_eq!(scratch, [123.0; 37]);
    assert_eq!(stream.sink().get_ref(), &bytes);
    assert_eq!(renderer.source().samples_rendered(), 74);
    render_pressure_pcm16(
        &mut renderer,
        &mut stream,
        &CancelGate::new(),
        &mut scratch,
        220,
    )
    .unwrap();
    assert_eq!(renderer.source().samples_rendered(), 514);
    assert_eq!(renderer.source().remaining_samples(), 0);
    let (output, summary) = stream.finish().unwrap();
    assert_eq!(output.into_inner(), wav);
    assert_eq!(summary.samples, 257);
    assert_eq!(summary.clipped_samples, clips as u64);
}
