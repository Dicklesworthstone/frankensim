//! G1/G3 reflected source histories sampled on the original moving microphone.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, MovingReflectionObserver,
    PreparedKinematicTrajectory, PreparedReflectedMovingReceiverSource, ReflectedTrajectoryBudget,
    TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{CompactSourceStrength, SourceKinematics, TimedCompactSource};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_couple::render::RenderError;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;

type Track<'a> = CompactSourceHistory<'a, &'a PreparedReflectedMovingReceiverSource<'a>>;
type Renderer<'a> = MovingCompactRenderer<'a, &'a PreparedReflectedMovingReceiverSource<'a>>;

fn path(times: &[f64], motion: impl Fn(f64) -> SourceKinematics) -> PreparedKinematicTrajectory {
    let samples: Vec<_> = times
        .iter()
        .map(|&source_time_s| TimedSourceKinematics {
            source_time_s,
            kinematics: motion(source_time_s),
        })
        .collect();
    PreparedKinematicTrajectory::new(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        &samples,
        0.8,
        TrajectoryBudget {
            maximum_segments: 8,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}
fn source(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            3.0 + 6.0 * t + 2.0 * t * t,
            2.0 + 3.0 * t + t * t * t,
            1.0 + t * t,
        ],
        velocity_m_s: [6.0 + 4.0 * t, 3.0 + 3.0 * t * t, 2.0 * t],
        acceleration_m_s2: [4.0, 6.0 * t, 2.0],
    }
}
fn microphone(t: f64) -> SourceKinematics {
    SourceKinematics {
        position_m: [
            12.0 - 20.0 * t + 4.0 * t * t,
            8.0 + 5.0 * t - 0.5 * t * t,
            5.0 + t * t,
        ],
        velocity_m_s: [-20.0 + 8.0 * t, 5.0 - t, 2.0 * t],
        acceleration_m_s2: [8.0, -1.0, 2.0],
    }
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
fn history(scale: f64) -> Vec<TimedCompactSource> {
    [-0.3, -0.1, 0.0, 0.1, 0.4]
        .into_iter()
        .map(|source_time_s| TimedCompactSource {
            source_time_s,
            strength: strength(source_time_s, scale),
        })
        .collect()
}
fn bounds() -> ReflectedTrajectoryBudget {
    ReflectedTrajectoryBudget {
        maximum_bound_nodes: 4096,
        maximum_subdivision_depth: 20,
    }
}
fn observer(receiver: &PreparedKinematicTrajectory) -> MovingReflectionObserver<'_> {
    MovingReflectionObserver::new(
        receiver,
        &[
            IdealAcousticPlane::new([0.0; 3], [1.0, 0.0, 0.0], IdealPlaneBoundary::Rigid).unwrap(),
            IdealAcousticPlane::new(
                [0.0; 3],
                [0.0, 1.0, 0.0],
                IdealPlaneBoundary::PressureRelease,
            )
            .unwrap(),
        ],
        4,
        bounds(),
        &mut || false,
    )
    .unwrap()
}
fn reflected<'a>(
    source: &'a PreparedKinematicTrajectory,
    observer: &'a MovingReflectionObserver<'a>,
) -> PreparedReflectedMovingReceiverSource<'a> {
    PreparedReflectedMovingReceiverSource::new(
        source,
        observer,
        MovingReceiverConfig {
            receiver_time_interval_s: [0.0, 0.2],
            minimum_distance_m: 0.1,
            maximum_time_residual_s: 1e-13,
        },
        bounds(),
        MovingReceiverBudget {
            maximum_segment_pairs: 1024,
            maximum_bound_nodes: 16_384,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap()
}
fn track<'a>(
    motion: &'a PreparedReflectedMovingReceiverSource<'a>,
    samples: &'a [TimedCompactSource],
) -> Track<'a> {
    CompactSourceHistory {
        motion,
        samples,
        source_bandwidth_hz: 1.0,
    }
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.04,
        sample_rate_hz: 1000,
        samples: 101,
        max_block: 101,
        receiver_bandwidth_hz: 10.0,
        maximum_time_residual_s: 1e-13,
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
fn build<'a>(sources: &[Track<'a>], config: MovingPressureConfig) -> Renderer<'a> {
    MovingCompactRenderer::with_moving_receiver_reflections(sources, config, budget(), &mut || {
        false
    })
    .unwrap()
}
fn image(v: [f64; 3], mask: usize) -> [f64; 3] {
    core::array::from_fn(|i| if mask & (1 << i) == 0 { v[i] } else { -v[i] })
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn reference_pressure(time: f64) -> f64 {
    let x = microphone(time).position_m;
    (0..4)
        .map(|mask| {
            let mut lo: f64 = -0.3;
            let mut hi: f64 = 0.4;
            // Independent source-image polynomial and bisection, never the
            // production virtual-receiver path or its source interpolant.
            for _ in 0..90 {
                let tau = lo.midpoint(hi);
                let y = image(source(tau).position_m, mask);
                let distance = norm(core::array::from_fn(|i| x[i] - y[i]));
                if tau + distance / 300.0 < time {
                    lo = tau;
                } else {
                    hi = tau;
                }
            }
            let tau = lo.midpoint(hi);
            let state = source(tau);
            let y = image(state.position_m, mask);
            let v = image(state.velocity_m_s, mask);
            let acceleration = image(state.acceleration_m_s2, mask);
            let s = strength(tau, 1.0);
            let force = image(s.force_on_fluid_n, mask);
            let force_rate = image(s.force_rate_on_fluid_n_s, mask);
            let displacement = core::array::from_fn(|i| x[i] - y[i]);
            let r = norm(displacement);
            let n = displacement.map(|v| v / r);
            let m = v.map(|v| v / 300.0);
            let mr = dot(n, m);
            let m2 = dot(m, m);
            let kappa = 1.0 - mr;
            let ar = dot(n, acceleration);
            let q = 1.2
                * (s.volume_flow_rate_m3_s2 / (r * kappa.powi(2))
                    + s.volume_flow_m3_s * (r * ar / 300.0 + 300.0 * (mr - m2))
                        / (r * r * kappa.powi(3)));
            let f = dot(n, force_rate) / (300.0 * r * kappa.powi(2))
                - dot(force, m) / (r * r * kappa.powi(2))
                + dot(n, force) * (1.0 - m2 + r * ar / 90_000.0) / (r * r * kappa.powi(3));
            (if mask & 2 == 0 { 1.0 } else { -1.0 }) * (q + f) / (4.0 * PI)
        })
        .sum()
}

#[test]
fn original_strength_interpolation_reaches_the_moving_reflected_pressure_sink() {
    let a = path(&[-0.3, -0.04, 0.11, 0.4], source);
    let b = path(&[0.0, 0.07, 0.2], microphone);
    let microphone = observer(&b);
    let motion = reflected(&a, &microphone);
    let samples = history(1.0);
    let tracks = [track(&motion, &samples)];
    let mut renderer = MovingCompactRenderer::with_moving_receiver_reflections(
        &tracks,
        config(),
        MovingPressureBudget {
            maximum_sources: 1,
            maximum_source_samples: samples.len(),
            maximum_source_receiver_pairs: config().samples * 4,
            ..budget()
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(renderer.source_count(), 1);
    assert_eq!(renderer.source_path_count(), 4);
    assert_eq!(
        renderer.maximum_retarded_iterations(),
        config().samples * 4 * 96
    );
    let sink: &mut dyn PressureRenderer = &mut renderer;
    sink.validate_sample_rate(1000).unwrap();
    sink.validate_sample_count(101).unwrap();
    let mut output = [0.0; 101];
    sink.block(&mut output).unwrap();
    for (index, value) in output.into_iter().enumerate() {
        let time = renderer.receiver_time_at(index as u64).unwrap();
        let expected = reference_pressure(time);
        assert!(
            (value - expected).abs() <= 2e-9 * expected.abs().max(1e-10),
            "{index} {value:e} != {expected:e}"
        );
        assert_eq!(
            renderer.observer_at(time).unwrap().observer_position_m(),
            b.kinematics_at(time).unwrap().position_m
        );
    }
}

#[test]
fn moving_reflected_callbacks_keep_partition_bits_coherence_and_cancellation_rollback() {
    let a = path(&[-0.3, -0.04, 0.11, 0.4], source);
    let b = path(&[0.0, 0.07, 0.2], microphone);
    let observer = observer(&b);
    let motion = reflected(&a, &observer);
    let samples = history(1.0);
    let tracks = [track(&motion, &samples)];
    let mut whole = build(&tracks, config());
    let mut expected = [0.0; 101];
    whole.block(&mut expected).unwrap();
    for block in [1, 7, 37] {
        let mut split = build(
            &tracks,
            MovingPressureConfig {
                max_block: block,
                ..config()
            },
        );
        let mut actual = [0.0; 101];
        for chunk in actual.chunks_mut(block) {
            split.block(chunk).unwrap();
        }
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
    }
    let mut cancelled = build(&tracks, config());
    let mut output = [17.5; 17];
    let mut polls = 0;
    assert!(matches!(
        cancelled.block_controlled(&mut output, &mut || {
            polls += 1;
            polls == 20
        }),
        Err(RenderError::Acoustic(AeroacError::Cancelled { .. }))
    ));
    assert_eq!(output, [17.5; 17]);
    assert_eq!(cancelled.samples_rendered(), 0);
    cancelled.block(&mut output).unwrap();
    assert_eq!(
        output.map(f64::to_bits),
        expected[..17]
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>()
            .as_slice()
    );
    let opposite = history(-1.0);
    let mut coherent = build(
        &[track(&motion, &samples), track(&motion, &opposite)],
        config(),
    );
    assert_eq!(coherent.source_count(), 2);
    assert_eq!(coherent.source_path_count(), 8);
    let mut pressure = [0.0; 101];
    coherent.block(&mut pressure).unwrap();
    assert!(pressure.iter().all(|v| v.abs() < 1e-20));
}

#[test]
fn complete_image_history_work_and_shared_moving_boundary_are_required() {
    let a = path(&[-0.3, -0.04, 0.11, 0.4], source);
    let b = path(&[0.0, 0.07, 0.2], microphone);
    let walls = observer(&b);
    let other_walls = observer(&b);
    let motion = reflected(&a, &walls);
    let other = reflected(&a, &other_walls);
    let samples = history(1.0);
    let tracks = [track(&motion, &samples)];
    assert!(
        MovingCompactRenderer::with_moving_receiver_reflections(
            &tracks,
            config(),
            MovingPressureBudget {
                maximum_source_receiver_pairs: config().samples * 4 - 1,
                ..budget()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(
        MovingCompactRenderer::with_moving_receiver_reflections(
            &[track(&motion, &samples), track(&other, &samples)],
            config(),
            budget(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        MovingCompactRenderer::with_moving_receiver_reflections(
            &tracks,
            MovingPressureConfig {
                maximum_time_residual_s: 1e-14,
                ..config()
            },
            budget(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        MovingCompactRenderer::with_moving_receiver_reflections(
            &tracks,
            MovingPressureConfig {
                first_receiver_time_s: 0.18,
                ..config()
            },
            budget(),
            &mut || false
        )
        .is_err()
    );
    let truncated: Vec<_> = [0.0, 0.1, 0.4]
        .into_iter()
        .map(|source_time_s| TimedCompactSource {
            source_time_s,
            strength: strength(source_time_s, 1.0),
        })
        .collect();
    // At 0.04 the direct path has a positive source time, but a reflected
    // path still needs pre-zero history. The whole render must refuse.
    assert!(
        motion
            .prepare_path(0, 0.04, &mut || false)
            .unwrap()
            .pressure
            .source_time_s()
            > 0.0
    );
    assert!(
        motion
            .prepare_path(1, 0.04, &mut || false)
            .unwrap()
            .pressure
            .source_time_s()
            < 0.0
    );
    assert!(
        MovingCompactRenderer::with_moving_receiver_reflections(
            &[track(&motion, &truncated)],
            config(),
            budget(),
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn reflected_receiver_doppler_is_included_in_audio_bandwidth_admission() {
    let straight = |origin, t| SourceKinematics {
        position_m: [origin - 90.0 * t, 1.0, 1.0],
        velocity_m_s: [-90.0, 0.0, 0.0],
        acceleration_m_s2: [0.0; 3],
    };
    let a = path(&[-0.3, 0.03], |t| straight(3.0, t));
    let b = path(&[0.0, 0.02], |t| straight(20.0, t));
    let wall =
        IdealAcousticPlane::new([0.0; 3], [1.0, 0.0, 0.0], IdealPlaneBoundary::Rigid).unwrap();
    let observer = MovingReflectionObserver::new(&b, &[wall], 2, bounds(), &mut || false).unwrap();
    let motion = PreparedReflectedMovingReceiverSource::new(
        &a,
        &observer,
        MovingReceiverConfig {
            receiver_time_interval_s: [0.0, 0.02],
            minimum_distance_m: 0.1,
            maximum_time_residual_s: 1e-13,
        },
        bounds(),
        MovingReceiverBudget {
            maximum_segment_pairs: 2,
            maximum_bound_nodes: 16_384,
            maximum_subdivision_depth: 20,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap();
    let samples: Vec<_> = [-0.3, 0.0, 0.03]
        .into_iter()
        .map(|source_time_s| TimedCompactSource {
            source_time_s,
            strength: strength(source_time_s, 1.0),
        })
        .collect();
    let tracks = [track(&motion, &samples)];
    let config = MovingPressureConfig {
        first_receiver_time_s: 0.0,
        samples: 21,
        receiver_bandwidth_hz: 1.1,
        ..config()
    };
    assert!(motion.path_maximum_source_time_rate(0).unwrap() < 1.000_001);
    assert!(
        MovingCompactRenderer::with_moving_receiver_reflections(
            &tracks,
            config,
            budget(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        MovingCompactRenderer::with_moving_receiver_reflections(
            &tracks,
            MovingPressureConfig {
                receiver_bandwidth_hz: 2.0,
                ..config
            },
            budget(),
            &mut || false
        )
        .is_ok()
    );
}
