//! G1/G3 actual reflected audio from one original accepted source history.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    IdealReflectionObserver, PreparedReflectedSourceTrajectory, PreparedSourceTrajectory,
    ReflectedTrajectoryBudget, TimedSourceKinematics, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, SourceKinematics,
    TimedCompactSource,
};
use fs_aeroac::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_couple::render::RenderError;
use fs_material::fluid_acoustics::FluidAcousticState;
use std::f64::consts::PI;

const RECEIVER: [f64; 3] = [12.0, 8.0, 5.0];

fn observer() -> MovingCompactObserver {
    MovingCompactObserver::new(
        FluidAcousticState::try_new(1.2, 300.0).unwrap(),
        RECEIVER,
        MovingCompactLimits {
            maximum_mach: 0.5,
            minimum_distance_m: 0.1,
        },
    )
    .unwrap()
}
fn motion(t: f64) -> SourceKinematics {
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
fn trajectory() -> PreparedSourceTrajectory {
    let knots = [-0.2, -0.03, 0.11, 0.4].map(|source_time_s| TimedSourceKinematics {
        source_time_s,
        kinematics: motion(source_time_s),
    });
    PreparedSourceTrajectory::new(
        observer(),
        &knots,
        TrajectoryBudget {
            maximum_segments: 16,
            maximum_payload_bytes: 131_072,
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 16,
            maximum_root_iterations: 64,
        },
        &mut || false,
    )
    .unwrap()
}
fn walls() -> IdealReflectionObserver {
    IdealReflectionObserver::new(
        observer(),
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
    )
    .unwrap()
}
fn reflected<'a>(
    trajectory: &'a PreparedSourceTrajectory,
    walls: &'a IdealReflectionObserver,
) -> PreparedReflectedSourceTrajectory<'a> {
    PreparedReflectedSourceTrajectory::new(
        trajectory,
        walls,
        ReflectedTrajectoryBudget {
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 16,
        },
        &mut || false,
    )
    .unwrap()
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
fn history() -> Vec<TimedCompactSource> {
    [-0.2, -0.05, 0.07, 0.2, 0.4]
        .into_iter()
        .map(|source_time_s| TimedCompactSource {
            source_time_s,
            strength: strength(source_time_s),
        })
        .collect()
}
fn config() -> MovingPressureConfig {
    MovingPressureConfig {
        first_receiver_time_s: 0.05,
        sample_rate_hz: 2000,
        samples: 241,
        max_block: 241,
        receiver_bandwidth_hz: 10.0,
        maximum_time_residual_s: 1e-13,
    }
}
fn budget() -> MovingPressureBudget {
    MovingPressureBudget {
        maximum_sources: 2,
        maximum_source_samples: 10,
        maximum_output_samples: 4096,
        maximum_payload_bytes: 131_072,
        maximum_source_receiver_pairs: 16_384,
    }
}
fn source<'a, 'b>(
    path: &'a PreparedReflectedSourceTrajectory<'b>,
    samples: &'a [TimedCompactSource],
) -> CompactSourceHistory<'a, &'a PreparedReflectedSourceTrajectory<'b>> {
    CompactSourceHistory {
        motion: path,
        samples,
        source_bandwidth_hz: 2.0,
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}
fn image(v: [f64; 3], mask: usize) -> [f64; 3] {
    core::array::from_fn(|i| if mask & (1 << i) == 0 { v[i] } else { -v[i] })
}

// Independent original-polynomial bisection and explicit SOURCE reflection.
// Neither the production quintics, virtual microphones nor Hermite strength
// coefficients enter this oracle.
fn reference_pressure(time: f64) -> f64 {
    let mut pressure = 0.0;
    for mask in 0..4 {
        let mut lo = -0.2_f64;
        let mut hi = 0.4_f64;
        for _ in 0..100 {
            let tau = lo.midpoint(hi);
            let y = image(motion(tau).position_m, mask);
            let r = norm(core::array::from_fn(|i| RECEIVER[i] - y[i]));
            if tau + r / 300.0 < time {
                lo = tau;
            } else {
                hi = tau;
            }
        }
        let tau = lo.midpoint(hi);
        let kinematics = motion(tau);
        let s = strength(tau);
        let y = image(kinematics.position_m, mask);
        let delta = core::array::from_fn(|i| RECEIVER[i] - y[i]);
        let r = norm(delta);
        let n = delta.map(|v| v / r);
        let m = image(kinematics.velocity_m_s, mask).map(|v| v / 300.0);
        let a = image(kinematics.acceleration_m_s2, mask);
        let force = image(s.force_on_fluid_n, mask);
        let force_rate = image(s.force_rate_on_fluid_n_s, mask);
        let mr = dot(n, m);
        let m2 = dot(m, m);
        let kappa = 1.0 - mr;
        let ar = dot(n, a);
        let q = 1.2
            * (s.volume_flow_rate_m3_s2 / (r * kappa * kappa)
                + s.volume_flow_m3_s * (r * ar / 300.0 + 300.0 * (mr - m2))
                    / (r * r * kappa.powi(3)));
        let f = dot(n, force_rate) / (300.0 * r * kappa * kappa)
            - dot(force, m) / (r * r * kappa * kappa)
            + dot(n, force) * (1.0 - m2 + r * ar / 90_000.0) / (r * r * kappa.powi(3));
        pressure += if mask & 2 == 0 { q + f } else { -q - f } / (4.0 * PI);
    }
    pressure
}

fn through_pressure_sink(renderer: &mut impl PressureRenderer, output: &mut [f64]) {
    renderer.validate_sample_rate(2000).unwrap();
    renderer.validate_sample_count(output.len() as u64).unwrap();
    renderer.block(output).unwrap();
}

#[test]
fn original_strength_history_reaches_actual_reflected_pressure_and_the_shared_sink() {
    let trajectory = trajectory();
    let walls = walls();
    let path = reflected(&trajectory, &walls);
    let samples = history();
    let mut renderer = MovingCompactRenderer::with_ideal_reflections(
        &[source(&path, &samples)],
        config(),
        MovingPressureBudget {
            maximum_sources: 1,
            maximum_source_samples: samples.len(),
            maximum_source_receiver_pairs: 4 * config().samples,
            ..budget()
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(renderer.source_count(), 1);
    assert_eq!(renderer.source_path_count(), 4);
    assert_eq!(
        renderer.maximum_retarded_iterations(),
        4 * 64 * config().samples
    );
    assert_eq!(
        renderer.observer_at(0.1).unwrap().observer_position_m(),
        RECEIVER
    );
    let mut pressure = vec![0.0; config().samples as usize];
    through_pressure_sink(&mut renderer, &mut pressure);
    for (index, &actual) in pressure.iter().enumerate() {
        let expected = reference_pressure(renderer.receiver_time_at(index as u64).unwrap());
        assert!(
            (actual - expected).abs() <= 3e-10 * expected.abs().max(1e-7),
            "{index}: actual={actual:e}, expected={expected:e}"
        );
    }
    assert_eq!(renderer.samples_rendered(), config().samples);
}

#[test]
fn reflected_callbacks_preserve_partition_bits_and_coherent_source_cancellation() {
    let trajectory = trajectory();
    let walls = walls();
    let path = reflected(&trajectory, &walls);
    let samples = history();
    let tracks = [source(&path, &samples)];
    let mut whole =
        MovingCompactRenderer::with_ideal_reflections(&tracks, config(), budget(), &mut || false)
            .unwrap();
    let mut split =
        MovingCompactRenderer::with_ideal_reflections(&tracks, config(), budget(), &mut || false)
            .unwrap();
    let mut expected = vec![0.0; 241];
    whole.block(&mut expected).unwrap();
    let mut observed = vec![0.0; 241];
    let mut offset = 0;
    for count in [1, 7, 32, 89, 112] {
        split.block(&mut observed[offset..offset + count]).unwrap();
        offset += count;
    }
    assert!(
        expected
            .iter()
            .zip(&observed)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    let opposite = samples
        .iter()
        .map(|sample| {
            let s = sample.strength;
            TimedCompactSource {
                source_time_s: sample.source_time_s,
                strength: CompactSourceStrength {
                    volume_flow_m3_s: -s.volume_flow_m3_s,
                    volume_flow_rate_m3_s2: -s.volume_flow_rate_m3_s2,
                    force_on_fluid_n: s.force_on_fluid_n.map(|v| -v),
                    force_rate_on_fluid_n_s: s.force_rate_on_fluid_n_s.map(|v| -v),
                },
            }
        })
        .collect::<Vec<_>>();
    let mut cancelled = MovingCompactRenderer::with_ideal_reflections(
        &[source(&path, &samples), source(&path, &opposite)],
        config(),
        budget(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(cancelled.source_count(), 2);
    assert_eq!(cancelled.source_path_count(), 8);
    cancelled.block(&mut observed).unwrap();
    assert!(observed.iter().all(|value| value.abs() < 1e-20));
}

#[test]
fn every_reflected_history_and_all_path_work_must_be_admitted() {
    let trajectory = trajectory();
    let walls = walls();
    let path = reflected(&trajectory, &walls);
    let samples = history();
    let cfg = config();
    assert!(
        MovingCompactRenderer::with_ideal_reflections(
            &[source(&path, &samples)],
            cfg,
            MovingPressureBudget {
                maximum_source_receiver_pairs: cfg.samples * 4 - 1,
                ..budget()
            },
            &mut || false
        )
        .is_err()
    );
    let partial = [0.0, 0.2, 0.4].map(|source_time_s| TimedCompactSource {
        source_time_s,
        strength: strength(source_time_s),
    });
    let cfg = MovingPressureConfig {
        first_receiver_time_s: observer()
            .prepare_emission(0.0, motion(0.0))
            .unwrap()
            .receiver_time_s()
            + 1e-6,
        ..cfg
    };
    assert!(
        MovingCompactRenderer::new(
            observer(),
            &[CompactSourceHistory {
                motion: &trajectory,
                samples: &partial,
                source_bandwidth_hz: 2.0,
            }],
            cfg,
            budget(),
            &mut || false
        )
        .is_ok()
    );
    assert!(
        MovingCompactRenderer::with_ideal_reflections(
            &[source(&path, &partial)],
            cfg,
            budget(),
            &mut || false
        )
        .is_err()
    );
    let another_wall_owner = walls.clone();
    let another_path = reflected(&trajectory, &another_wall_owner);
    assert!(
        MovingCompactRenderer::with_ideal_reflections(
            &[source(&path, &samples), source(&another_path, &samples)],
            config(),
            budget(),
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn cancelling_inside_an_image_path_leaves_output_and_physical_clock_unchanged() {
    let trajectory = trajectory();
    let walls = walls();
    let path = reflected(&trajectory, &walls);
    let samples = history();
    let tracks = [source(&path, &samples)];
    let mut renderer =
        MovingCompactRenderer::with_ideal_reflections(&tracks, config(), budget(), &mut || false)
            .unwrap();
    let mut output = [17.5; 17];
    let mut polls = 0;
    assert!(matches!(
        renderer.block_controlled(&mut output, &mut || {
            polls += 1;
            polls == 14
        }),
        Err(RenderError::Acoustic(AeroacError::Cancelled { .. }))
    ));
    assert_eq!(output, [17.5; 17]);
    assert_eq!(renderer.samples_rendered(), 0);
    let mut fresh =
        MovingCompactRenderer::with_ideal_reflections(&tracks, config(), budget(), &mut || false)
            .unwrap();
    let mut expected = [0.0; 17];
    fresh.block(&mut expected).unwrap();
    renderer.block(&mut output).unwrap();
    assert_eq!(output.map(f64::to_bits), expected.map(f64::to_bits));
}
