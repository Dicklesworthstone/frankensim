//! G1/G3 moving microphones on the original pressure clock and PCM pipeline.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, PreparedKinematicTrajectory,
    PreparedMovingReceiverSource, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{CompactSourceStrength, SourceKinematics, TimedCompactSource};
use fs_couple::pcm_wav::decimate::Decimator;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::DecimatedRenderer;
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, render_pressure_pcm16};
use fs_couple::render::RenderError;
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::{f64::consts::PI, io::Cursor};

type Track<'a> = CompactSourceHistory<'a, &'a PreparedMovingReceiverSource<'a>>;
type Renderer<'a> = MovingCompactRenderer<'a, &'a PreparedMovingReceiverSource<'a>>;

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.2, 300.0).unwrap()
}
fn path(
    times: &[f64],
    kinematics: impl Fn(f64) -> SourceKinematics,
    maximum_mach: f64,
) -> PreparedKinematicTrajectory {
    let samples: Vec<_> = times
        .iter()
        .map(|&source_time_s| TimedSourceKinematics {
            source_time_s,
            kinematics: kinematics(source_time_s),
        })
        .collect();
    PreparedKinematicTrajectory::new(
        medium(),
        &samples,
        maximum_mach,
        TrajectoryBudget {
            maximum_segments: 8,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 16,
            maximum_root_iterations: 64,
        },
        &mut || false,
    )
    .unwrap()
}
fn pair<'a>(
    source: &'a PreparedKinematicTrajectory,
    receiver: &'a PreparedKinematicTrajectory,
) -> PreparedMovingReceiverSource<'a> {
    PreparedMovingReceiverSource::new(
        source,
        receiver,
        MovingReceiverConfig {
            receiver_time_interval_s: [0.0, 0.16],
            minimum_distance_m: 0.1,
            maximum_time_residual_s: 1e-12,
        },
        MovingReceiverBudget {
            maximum_segment_pairs: 64,
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 16,
            maximum_root_iterations: 64,
        },
        &mut || false,
    )
    .unwrap()
}
fn budget() -> MovingPressureBudget {
    MovingPressureBudget {
        maximum_sources: 4,
        maximum_source_samples: 128,
        maximum_output_samples: 4096,
        maximum_payload_bytes: 131_072,
        maximum_source_receiver_pairs: 16_384,
    }
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.04,
        sample_rate_hz: 1000,
        samples: 101,
        max_block: 101,
        receiver_bandwidth_hz: 10.0,
        maximum_time_residual_s: 1e-12,
    }
}
fn build<'a>(sources: &[Track<'a>], config: MovingPressureConfig) -> Renderer<'a> {
    MovingCompactRenderer::with_moving_receiver(sources, config, budget(), &mut || false).unwrap()
}
fn strength(t: f64, scale: f64) -> CompactSourceStrength {
    CompactSourceStrength {
        volume_flow_m3_s: scale * (0.001 + 0.002 * t + 0.003 * t * t - 0.001 * t * t * t),
        volume_flow_rate_m3_s2: scale * (0.002 + 0.006 * t - 0.003 * t * t),
        force_on_fluid_n: [
            scale * (2.0 + 3.0 * t - 2.0 * t * t + t * t * t),
            scale * (-1.0 + 2.0 * t * t),
            scale * (0.5 - 0.4 * t),
        ],
        force_rate_on_fluid_n_s: [
            scale * (3.0 - 4.0 * t + 3.0 * t * t),
            scale * 4.0 * t,
            scale * -0.4,
        ],
    }
}
fn history(times: &[f64], scale: f64) -> Vec<TimedCompactSource> {
    times
        .iter()
        .map(|&source_time_s| TimedCompactSource {
            source_time_s,
            strength: strength(source_time_s, scale),
        })
        .collect()
}
fn track<'a>(
    motion: &'a PreparedMovingReceiverSource<'a>,
    samples: &'a [TimedCompactSource],
) -> Track<'a> {
    CompactSourceHistory {
        motion,
        samples,
        source_bandwidth_hz: 1.0,
    }
}
fn straight_receiver(t: f64, velocity: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [12.0 + velocity * t, 0.0, 0.0],
        velocity_m_s: [velocity, 0.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    }
}
fn source_a(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            0.3 + 18.0 * t + 8.0 * t * t,
            -0.2 + 4.0 * t - 20.0 * t * t + 5.0 * t * t * t,
            0.1 - 3.0 * t + 4.0 * t * t,
        ],
        velocity_m_s: [
            18.0 + 16.0 * t,
            4.0 - 40.0 * t + 15.0 * t * t,
            -3.0 + 8.0 * t,
        ],
        acceleration_m_s2: [16.0, -40.0 + 30.0 * t, 8.0],
    }
}
fn source_b(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [-2.0 + 4.0 * t, -3.0 - 2.0 * t, 2.0 + 0.5 * t * t],
        velocity_m_s: [4.0, -2.0, t],
        acceleration_m_s2: [0.0, 0.0, 1.0],
    }
}
fn receiver(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            12.0 - 25.0 * t + 8.0 * t * t,
            4.0 + 10.0 * t - 3.0 * t * t,
            1.5 + 2.0 * t * t,
        ],
        velocity_m_s: [-25.0 + 16.0 * t, 10.0 - 6.0 * t, 4.0 * t],
        acceleration_m_s2: [16.0, -6.0, 4.0],
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn reference_pressure(time: f64, source: fn(f64) -> SourceKinematics, scale: f64) -> f64 {
    let x = receiver(time).position_m;
    let mut lo = -0.1;
    let mut hi = 0.3;
    // Independent bisection of the original analytic path, not the prepared
    // quintic controls or the production safeguarded Newton implementation.
    for _ in 0..90 {
        let tau = 0.5 * (lo + hi);
        let y = source(tau).position_m;
        let radius = norm(core::array::from_fn(|i| x[i] - y[i]));
        if tau + radius / 300.0 < time {
            lo = tau;
        } else {
            hi = tau;
        }
    }
    let tau = 0.5 * (lo + hi);
    let motion = source(tau);
    let s = strength(tau, scale);
    let offset = core::array::from_fn(|i| x[i] - motion.position_m[i]);
    let r = norm(offset);
    let n = offset.map(|v| v / r);
    let m = motion.velocity_m_s.map(|v| v / 300.0);
    let mr = dot(n, m);
    let kappa = 1.0 - mr;
    let m2 = dot(m, m);
    let ar = dot(n, motion.acceleration_m_s2);
    let q = 1.2
        * (s.volume_flow_rate_m3_s2 / (r * kappa * kappa)
            + s.volume_flow_m3_s * (r * ar / 300.0 + 300.0 * (mr - m2)) / (r * r * kappa.powi(3)));
    let f = dot(n, s.force_rate_on_fluid_n_s) / (300.0 * r * kappa * kappa)
        - dot(s.force_on_fluid_n, m) / (r * r * kappa * kappa)
        + dot(n, s.force_on_fluid_n) * (1.0 - m2 + r * ar / 90_000.0) / (r * r * kappa.powi(3));
    (q + f) / (4.0 * PI)
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-9 * expected.abs().max(1e-10),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn moving_microphone_changes_retarded_source_clock_without_artificial_pressure_gain() {
    let source = path(&[-0.1, 0.3], |_| SourceKinematics::default(), 0.0);
    let accepted = history(&[-0.1, 0.0, 0.1, 0.3], 1.0);
    for velocity in [-60.0, 0.0, 60.0] {
        let receiver = path(&[0.0, 0.16], |t| straight_receiver(t, velocity), 0.3);
        let motion = pair(&source, &receiver);
        let mut renderer = build(&[track(&motion, &accepted)], config());
        let mut pressure = [0.0; 101];
        renderer.block(&mut pressure).unwrap();
        for (index, actual) in pressure.into_iter().enumerate() {
            let time = renderer.receiver_time_at(index as u64).unwrap();
            let radius = 12.0 + velocity * time;
            let tau = (1.0 - velocity / 300.0) * time - 12.0 / 300.0;
            let s = strength(tau, 1.0);
            let expected = (1.2 * s.volume_flow_rate_m3_s2 / radius
                + s.force_on_fluid_n[0] / (radius * radius)
                + s.force_rate_on_fluid_n_s[0] / (300.0 * radius))
                / (4.0 * PI);
            close(actual, expected);
        }
        assert_eq!(renderer.maximum_retarded_iterations(), 101 * 64);
    }
}

#[test]
fn two_curved_sources_share_the_actual_curved_microphone_and_block_clock() {
    let a = path(&[-0.1, 0.07, 0.3], source_a, 0.6);
    let b = path(&[-0.1, 0.19, 0.3], source_b, 0.6);
    let receiver = path(&[0.0, 0.06, 0.2], receiver, 0.6);
    let first = pair(&a, &receiver);
    let second = pair(&b, &receiver);
    let first_history = history(&[-0.1, 0.0, 0.11, 0.3], 1.0);
    let second_history = history(&[-0.1, -0.01, 0.2, 0.3], -0.7);
    let sources = [
        track(&first, &first_history),
        track(&second, &second_history),
    ];
    let mut whole = build(&sources, config());
    let mut expected = [0.0; 101];
    whole.block(&mut expected).unwrap();
    for (index, actual) in expected.into_iter().enumerate() {
        let time = whole.receiver_time_at(index as u64).unwrap();
        close(
            actual,
            reference_pressure(time, source_a, 1.0) + reference_pressure(time, source_b, -0.7),
        );
    }
    for block in [1, 7, 37] {
        let mut split = build(&sources, config());
        let mut actual = [0.0; 101];
        for chunk in actual.chunks_mut(block) {
            split.block(chunk).unwrap();
        }
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        assert_eq!(split.maximum_retarded_iterations(), 2 * 101 * 64);
        assert_eq!(split.samples_rendered(), 101);
        assert_eq!(split.remaining_samples(), 0);
    }
}

#[test]
fn receiver_doppler_and_original_microphone_identity_are_admission_gates() {
    let source = path(&[-0.1, 0.3], |_| SourceKinematics::default(), 0.0);
    let receiver = path(&[0.0, 0.16], |t| straight_receiver(t, -60.0), 0.3);
    let other_receiver = path(&[0.0, 0.16], |t| straight_receiver(t, -60.0), 0.3);
    let motion = pair(&source, &receiver);
    let other_motion = pair(&source, &other_receiver);
    let accepted = history(&[-0.1, 0.0, 0.1, 0.3], 1.0);
    let source = track(&motion, &accepted);
    let cfg = MovingPressureConfig {
        receiver_bandwidth_hz: 1.1,
        ..config()
    };
    // A stationary emitter has no source-side Doppler. The approaching
    // microphone alone raises the source-clock rate from one to 1.2.
    assert!(
        MovingCompactRenderer::with_moving_receiver(&[source], cfg, budget(), &mut || false)
            .is_err()
    );
    assert!(
        MovingCompactRenderer::with_moving_receiver(
            &[source],
            MovingPressureConfig {
                receiver_bandwidth_hz: 1.5,
                ..cfg
            },
            budget(),
            &mut || false,
        )
        .is_ok()
    );
    // Equal coordinates do not establish one original accepted receiver owner.
    assert!(
        MovingCompactRenderer::with_moving_receiver(
            &[source, track(&other_motion, &accepted)],
            config(),
            budget(),
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn finite_histories_and_one_shot_cancellation_preserve_accepted_output() {
    let source = path(&[-0.1, 0.07, 0.3], source_a, 0.6);
    let receiver = path(&[0.0, 0.06, 0.2], receiver, 0.6);
    let motion = pair(&source, &receiver);
    let accepted = history(&[-0.1, 0.0, 0.1, 0.3], 1.0);
    let sources = [track(&motion, &accepted)];
    let mut expected = [0.0; 101];
    build(&sources, config()).block(&mut expected).unwrap();
    for cancel_at in [1, 3, 9] {
        let mut renderer = build(&sources, config());
        let mut prefix = [0.0; 3];
        renderer.block(&mut prefix).unwrap();
        for (actual, reference) in prefix.iter().zip(&expected[..3]) {
            assert_eq!(actual.to_bits(), reference.to_bits());
        }
        let mut output = [123.0; 8];
        let mut polls = 0;
        assert!(matches!(
            renderer.block_controlled(&mut output, &mut || {
                polls += 1;
                polls == cancel_at
            }),
            Err(RenderError::Acoustic(AeroacError::Cancelled { .. }))
        ));
        assert_eq!(polls, cancel_at);
        assert_eq!(output, [123.0; 8]);
        assert_eq!(renderer.samples_rendered(), 3);
        renderer.block(&mut output).unwrap();
        for (actual, reference) in output.iter().zip(&expected[3..11]) {
            assert_eq!(actual.to_bits(), reference.to_bits());
        }
    }
    let incomplete = history(&[0.03, 0.1, 0.3], 1.0);
    assert!(
        MovingCompactRenderer::with_moving_receiver(
            &[track(&motion, &incomplete)],
            config(),
            budget(),
            &mut || false,
        )
        .is_err()
    );
    for cfg in [
        MovingPressureConfig {
            first_receiver_time_s: -0.001,
            ..config()
        },
        MovingPressureConfig {
            samples: 162,
            max_block: 162,
            ..config()
        },
        MovingPressureConfig {
            maximum_time_residual_s: 1e-13,
            ..config()
        },
    ] {
        assert!(
            MovingCompactRenderer::with_moving_receiver(&sources, cfg, budget(), &mut || false)
                .is_err()
        );
    }

    // The exact source endpoint tau=0 arrives at t=12/300=0.04 with zero
    // residual. Endpoint reconciliation must still respect the prepared
    // pair's configured clock policy, even for a single-sample observation.
    let endpoint_source = path(&[-0.1, 0.3], |_| SourceKinematics::default(), 0.0);
    let endpoint_receiver = path(&[0.0, 0.16], |t| straight_receiver(t, 0.0), 0.0);
    let endpoint_motion = pair(&endpoint_source, &endpoint_receiver);
    let endpoint_history = history(&[0.0, 0.1], 1.0);
    let endpoint_tracks = [track(&endpoint_motion, &endpoint_history)];
    let endpoint_config = MovingPressureConfig {
        samples: 1,
        max_block: 1,
        ..config()
    };
    assert!(
        MovingCompactRenderer::with_moving_receiver(
            &endpoint_tracks,
            endpoint_config,
            budget(),
            &mut || false,
        )
        .is_ok()
    );
    assert!(
        MovingCompactRenderer::with_moving_receiver(
            &endpoint_tracks,
            MovingPressureConfig {
                maximum_time_residual_s: 1e-13,
                ..endpoint_config
            },
            budget(),
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn moving_microphone_uses_original_decimator_and_pcm_without_a_second_sample_clock() {
    let source = path(&[-0.1, 0.07, 0.3], source_a, 0.6);
    let receiver = path(&[0.0, 0.06, 0.2], receiver, 0.6);
    let motion = pair(&source, &receiver);
    let accepted = history(&[-0.1, 0.0, 0.1, 0.3], 1.0);
    let sources = [track(&motion, &accepted)];
    let cfg = MovingPressureConfig {
        sample_rate_hz: 96_000,
        samples: 128,
        max_block: 26,
        ..config()
    };
    let mut direct = build(&sources, cfg);
    let mut pressure = [0.0; 128];
    for chunk in pressure.chunks_mut(26) {
        direct.block(chunk).unwrap();
    }
    let mut filter = Decimator::new(2, 1).unwrap();
    let expected: Vec<_> = pressure
        .chunks_exact(2)
        .map(|group| {
            let value = filter.preview(group).unwrap()[0];
            filter.commit();
            value
        })
        .collect();
    let (expected_wav, expected_clips) = encode_pcm16_wav(&expected, 48_000, 0.02).unwrap();
    let mut renderer = DecimatedRenderer::new(build(&sources, cfg), 96_000, 48_000, 13).unwrap();
    assert_eq!(
        renderer.info().delay_output_samples,
        filter.delay_output_frames()
    );
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 48_000, 0.02, 13).unwrap();
    render_pressure_pcm16(
        &mut renderer,
        &mut stream,
        &CancelGate::new(),
        &mut [0.0; 13],
        64,
    )
    .unwrap();
    assert_eq!(renderer.source().samples_rendered(), 128);
    let (sink, summary) = stream.finish().unwrap();
    assert_eq!(sink.into_inner(), expected_wav);
    assert_eq!(summary.clipped_samples, expected_clips as u64);
    assert_eq!(summary.samples, 64);
}
