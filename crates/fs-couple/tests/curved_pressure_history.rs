//! G1/G3 physical curved trajectories through the existing receiver/PCM path.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    PreparedSourceTrajectory, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    TimedCompactSource, UniformSourceMotion,
};
use fs_couple::pcm_wav::decimate::Decimator;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::DecimatedRenderer;
use fs_couple::pcm_wav::observation::moving::{
    CompactMotion, CompactSourceHistory, CompactSourceMotion, MovingCompactRenderer,
    MovingPressureBudget, MovingPressureConfig,
};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, render_pressure_pcm16};
use fs_couple::render::RenderError;
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::{f64::consts::PI, io::Cursor};

const RECEIVER: [f64; 3] = [8.0, 5.0, 2.0];

fn observer() -> MovingCompactObserver {
    MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        RECEIVER,
        MovingCompactLimits {
            maximum_mach: 0.8,
            minimum_distance_m: 0.1,
        },
    )
    .unwrap()
}
fn trajectory_budget() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 16,
        maximum_payload_bytes: 131_072,
        maximum_bound_nodes: 4096,
        maximum_subdivision_depth: 16,
        maximum_root_iterations: 64,
    }
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

// Original polynomial motion: quintic endpoint interpolation must reproduce
// these lower-degree components, including all three accelerations.
fn kinematics(t: f64) -> SourceKinematics {
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
fn trajectory() -> PreparedSourceTrajectory {
    let samples = [-0.1, 0.1, 0.3].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: kinematics(source_time_s),
    });
    PreparedSourceTrajectory::new(observer(), &samples, trajectory_budget(), &mut || false).unwrap()
}
fn strength(t: f64) -> CompactSourceStrength {
    CompactSourceStrength {
        volume_flow_m3_s: 0.001 + 0.002 * t + 0.003 * t * t - 0.001 * t * t * t,
        volume_flow_rate_m3_s2: 0.002 + 0.006 * t - 0.003 * t * t,
        force_on_fluid_n: [
            2.0 + 3.0 * t - 2.0 * t * t + t * t * t,
            -1.0 + 2.0 * t * t,
            0.5 - 0.4 * t,
        ],
        force_rate_on_fluid_n_s: [3.0 - 4.0 * t + 3.0 * t * t, 4.0 * t, -0.4],
    }
}
fn history(times: &[f64]) -> Vec<TimedCompactSource> {
    times
        .iter()
        .map(|&source_time_s| TimedCompactSource {
            source_time_s,
            strength: strength(source_time_s),
        })
        .collect()
}
fn build<'a, M: CompactMotion<Observer = MovingCompactObserver>>(
    sources: &[CompactSourceHistory<'a, M>],
    cfg: MovingPressureConfig,
) -> MovingCompactRenderer<'a, M> {
    MovingCompactRenderer::new(observer(), sources, cfg, budget(), &mut || false).unwrap()
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}
fn retarded(time: f64) -> f64 {
    let mut lo = -0.1;
    let mut hi = 0.3;
    // Independent scalar bisection of the original polynomial, without using
    // the prepared trajectory's coefficients, search, or Newton iteration.
    for _ in 0..90 {
        let tau = 0.5 * (lo + hi);
        let y = kinematics(tau).position_m;
        let r = norm(core::array::from_fn(|i| RECEIVER[i] - y[i]));
        if tau + r / 300.0 < time {
            lo = tau;
        } else {
            hi = tau;
        }
    }
    0.5 * (lo + hi)
}
fn reference_pressure(time: f64) -> f64 {
    let tau = retarded(time);
    let motion = kinematics(tau);
    let s = strength(tau);
    let offset = core::array::from_fn(|i| RECEIVER[i] - motion.position_m[i]);
    let r = norm(offset);
    let n = offset.map(|x| x / r);
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
        (actual - expected).abs() <= 5e-9 * expected.abs().max(1e-9),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn curved_motion_and_source_derivatives_reach_the_physical_receiver_clock() {
    let path = trajectory();
    let accepted = history(&[-0.1, 0.0, 0.1, 0.3]);
    let sources = [CompactSourceHistory {
        motion: &path,
        samples: &accepted,
        source_bandwidth_hz: 1.0,
    }];
    let mut renderer = build(&sources, config());
    assert_eq!(renderer.maximum_retarded_iterations(), 101 * 64);
    let mut pressure = [0.0; 101];
    renderer.block(&mut pressure).unwrap();
    for (index, actual) in pressure.into_iter().enumerate() {
        let time = renderer.receiver_time_at(index as u64).unwrap();
        close(actual, reference_pressure(time));
        let prepared = path.prepare_receiver(time, 1e-12, &mut || false).unwrap();
        assert!((prepared.source_time_s() - retarded(time)).abs() <= 2e-12);
        assert_eq!(prepared.receiver_time_s(), time);
        assert!(prepared.retarded_time_residual_s() <= 1e-12);
        let source = prepared.source_time_s();
        close(
            path.kinematics_at(source).unwrap().acceleration_m_s2[1],
            kinematics(source).acceleration_m_s2[1],
        );
    }
    assert_eq!(renderer.remaining_samples(), 0);
}

#[test]
fn mixed_motion_owners_sum_coherently_and_preserve_partition_bits() {
    let path = trajectory();
    let accepted = history(&[-0.1, 0.0, 0.1, 0.3]);
    let uniform = UniformSourceMotion::new(0.0, [1.0, 0.0, -1.0], [-20.0, 0.0, 0.0]).unwrap();
    let curved = [CompactSourceHistory {
        motion: &path,
        samples: &accepted,
        source_bandwidth_hz: 1.0,
    }];
    let straight = [CompactSourceHistory {
        motion: uniform,
        samples: &accepted,
        source_bandwidth_hz: 1.0,
    }];
    let mut original: MovingCompactRenderer<'_> = build(&straight, config());
    assert_eq!(original.maximum_retarded_iterations(), 0);
    let mut straight_pressure = [0.0; 101];
    original.block(&mut straight_pressure).unwrap();
    let mut curved_pressure = [0.0; 101];
    build(&curved, config())
        .block(&mut curved_pressure)
        .unwrap();
    let mixed = [
        CompactSourceHistory {
            motion: CompactSourceMotion::Trajectory(&path),
            samples: &accepted,
            source_bandwidth_hz: 1.0,
        },
        CompactSourceHistory {
            motion: CompactSourceMotion::Uniform(uniform),
            samples: &accepted,
            source_bandwidth_hz: 1.0,
        },
    ];
    let expected: [f64; 101] = core::array::from_fn(|i| curved_pressure[i] + straight_pressure[i]);
    for block in [1, 7, 37, 101] {
        let mut renderer = build(&mixed, config());
        let mut actual = [0.0; 101];
        for chunk in actual.chunks_mut(block) {
            renderer.block(chunk).unwrap();
        }
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        assert_eq!(renderer.maximum_retarded_iterations(), 101 * 64);
        assert_eq!(renderer.samples_rendered(), 101);
    }
}

#[test]
fn trajectory_and_strength_coverage_and_original_observer_are_independent_gates() {
    let path = trajectory();
    let accepted = history(&[-0.3, -0.1, 0.1, 0.3, 0.4]);
    let track = CompactSourceHistory {
        motion: &path,
        samples: &accepted,
        source_bandwidth_hz: 1.0,
    };
    let attempt = |microphone, track, cfg| {
        MovingCompactRenderer::new(microphone, &[track], cfg, budget(), &mut || false).is_ok()
    };
    assert!(attempt(observer(), track, config()));
    let [first, last] = path.receiver_time_interval_s();
    for time in [first - 0.001, last + 0.001] {
        assert!(!attempt(
            observer(),
            track,
            MovingPressureConfig {
                first_receiver_time_s: time,
                samples: 1,
                ..config()
            }
        ));
    }
    let short = history(&[0.05, 0.07]);
    assert!(!attempt(
        observer(),
        CompactSourceHistory {
            samples: &short,
            ..track
        },
        config()
    ));
    assert!(!attempt(
        observer(),
        track,
        MovingPressureConfig {
            samples: 401,
            ..config()
        }
    ));
    for microphone in [
        MovingCompactObserver::new(
            FluidAcousticState::try_new(1.3, 300.0).unwrap(),
            RECEIVER,
            observer().limits(),
        )
        .unwrap(),
        MovingCompactObserver::new(observer().medium(), [9.0, 5.0, 2.0], observer().limits())
            .unwrap(),
        MovingCompactObserver::new(
            observer().medium(),
            RECEIVER,
            MovingCompactLimits {
                maximum_mach: 0.7,
                ..observer().limits()
            },
        )
        .unwrap(),
    ] {
        assert!(!attempt(microphone, track, config()));
    }
}

#[test]
fn a_turn_toward_the_microphone_requires_the_whole_trajectory_doppler_bound() {
    let microphone =
        MovingCompactObserver::new(observer().medium(), [30.0, 0.0, 0.0], observer().limits())
            .unwrap();
    let samples = [0.0, 1.0].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: SourceKinematics {
            position_m: [120.0 * source_time_s * (source_time_s - 1.0), 0.0, 0.0],
            velocity_m_s: [-120.0 + 240.0 * source_time_s, 0.0, 0.0],
            acceleration_m_s2: [240.0, 0.0, 0.0],
        },
    });
    let path =
        PreparedSourceTrajectory::new(microphone, &samples, trajectory_budget(), &mut || false)
            .unwrap();
    let accepted = history(&[0.0, 0.25, 0.5, 0.75, 1.0]);
    let track = CompactSourceHistory {
        motion: &path,
        samples: &accepted,
        source_bandwidth_hz: 1.0,
    };
    let cfg = MovingPressureConfig {
        first_receiver_time_s: 0.1,
        sample_rate_hz: 10,
        samples: 11,
        max_block: 11,
        receiver_bandwidth_hz: 0.9,
        ..config()
    };
    let first = path.prepare_receiver(0.1, 1e-12, &mut || false).unwrap();
    assert!(first.source_time_rate() < 0.9);
    assert!(path.maximum_source_time_rate() > 1.6);
    // This exact quadratic has velocity controls [-120,-60,0,60,120]
    // m/s. Its outward control-hull bound has ample margin below factor 2.
    assert!(path.maximum_source_time_rate() < 2.0);
    assert!(
        MovingCompactRenderer::new(microphone, &[track], cfg, budget(), &mut || false).is_err()
    );
    assert!(
        MovingCompactRenderer::new(
            microphone,
            &[track],
            MovingPressureConfig {
                receiver_bandwidth_hz: 2.0,
                ..cfg
            },
            budget(),
            &mut || false
        )
        .is_ok()
    );
}

#[test]
fn cancellation_in_the_motion_owner_is_atomic_even_at_a_retained_endpoint() {
    let path = trajectory();
    let accepted = history(&[-0.1, 0.0, 0.1, 0.3]);
    let sources = [CompactSourceHistory {
        motion: &path,
        samples: &accepted,
        source_bandwidth_hz: 1.0,
    }];
    // Two renderer polls precede the trajectory call. A one-shot cancellation
    // inside that owner must not be consumed by exact-endpoint reconciliation.
    for (time, cancel_at) in [(path.receiver_time_interval_s()[0], 3), (0.04, 4)] {
        let cfg = MovingPressureConfig {
            first_receiver_time_s: time,
            samples: 1,
            ..config()
        };
        let mut baseline = build(&sources, cfg);
        let mut expected = [0.0];
        baseline.block(&mut expected).unwrap();
        let mut renderer = build(&sources, cfg);
        let mut out = [123.0];
        let mut polls = 0;
        assert!(matches!(
            renderer.block_controlled(&mut out, &mut || {
                polls += 1;
                polls == cancel_at
            }),
            Err(RenderError::Acoustic(AeroacError::Cancelled { .. }))
        ));
        assert_eq!(polls, cancel_at);
        assert_eq!(out, [123.0]);
        assert_eq!(renderer.samples_rendered(), 0);
        renderer.block(&mut out).unwrap();
        assert_eq!(out.map(f64::to_bits), expected.map(f64::to_bits));
    }
}

#[test]
fn curved_histories_use_the_original_decimator_and_pcm_writer() {
    let path = trajectory();
    let accepted = history(&[-0.1, 0.0, 0.1, 0.3]);
    let sources = [CompactSourceHistory {
        motion: &path,
        samples: &accepted,
        source_bandwidth_hz: 1.0,
    }];
    let cfg = MovingPressureConfig {
        sample_rate_hz: 96_000,
        samples: 128,
        max_block: 16,
        ..config()
    };
    let mut direct = build(&sources, cfg);
    let mut pressure = [0.0; 128];
    for chunk in pressure.chunks_mut(16) {
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
    assert_eq!(renderer.source().remaining_samples(), 0);
    let (sink, summary) = stream.finish().unwrap();
    assert_eq!(sink.into_inner(), expected_wav);
    assert_eq!(summary.clipped_samples, expected_clips as u64);
    assert_eq!(summary.samples, 64);
}
