//! Native accepted rigid motion and rotating forces through physical receivers.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::trajectory::{
    MovingReceiverBudget, MovingReceiverConfig, PreparedKinematicTrajectory,
    PreparedMovingReceiverSource, PreparedSourceTrajectory, TrajectoryBudget,
};
use fs_aeroac::curle3d::moving::{MovingCompactLimits, MovingCompactObserver};
use fs_couple::pcm_wav::observation::moving::rigid::{
    BodyFrameCompactSourceStrength, RigidAcousticAttachment, RigidAcousticError,
    RigidAcousticHistory, RigidHistoryBudget,
};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingCompactRenderer, MovingPressureBudget, MovingPressureConfig,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};
use std::f64::consts::PI;

const SOUND_SPEED: f64 = 300.0;
const DENSITY: f64 = 1.2;
const OMEGA: f64 = 18.0;
const RADIUS: f64 = 0.7;
const DURATION: f64 = 0.24;
const RECEIVER: [f64; 3] = [8.0, 5.0, 2.0];

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(DENSITY, SOUND_SPEED).unwrap()
}
fn properties() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap()
}
fn state(pose: Pose, velocity_world: Vec3, omega_body: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        pose,
        velocity_world.scale(2.0),
        properties()
            .angular_momentum_body_checked(omega_body)
            .unwrap(),
    )
    .unwrap()
}
fn observer(position: [f64; 3]) -> MovingCompactObserver {
    MovingCompactObserver::new(
        medium(),
        position,
        MovingCompactLimits {
            maximum_mach: 0.2,
            minimum_distance_m: 0.1,
        },
    )
    .unwrap()
}
fn trajectory_budget() -> TrajectoryBudget {
    TrajectoryBudget {
        maximum_segments: 64,
        maximum_payload_bytes: 262_144,
        maximum_bound_nodes: 16_384,
        maximum_subdivision_depth: 20,
        maximum_root_iterations: 96,
    }
}
fn pressure_budget() -> MovingPressureBudget {
    MovingPressureBudget {
        maximum_sources: 1,
        maximum_source_samples: 128,
        maximum_output_samples: 256,
        maximum_payload_bytes: 262_144,
        maximum_source_receiver_pairs: 256,
    }
}
fn history(arm: Vec3, maximum_samples: usize) -> RigidAcousticHistory {
    RigidAcousticHistory::new(
        RigidAcousticAttachment::new(arm).unwrap(),
        RigidHistoryBudget {
            maximum_samples,
            maximum_payload_bytes: 262_144,
        },
    )
    .unwrap()
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn close_vector(actual: [f64; 3], expected: [f64; 3]) {
    for i in 0..3 {
        close(actual[i], expected[i]);
    }
}

#[test]
fn rotating_constant_body_force_at_com_has_a_radiating_transport_derivative() {
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let initial = state(Pose::identity(), Vec3::ZERO, Vec3::new(0.0, 0.0, OMEGA));
    let time = 0.7 / OMEGA;
    let accepted = integrator
        .step(initial, properties(), Wrench::ZERO, time)
        .unwrap();
    let prepared = integrator
        .prepare_point_motion(accepted.state_after, properties(), Wrench::ZERO)
        .unwrap();
    let sample = RigidAcousticAttachment::new(Vec3::ZERO)
        .unwrap()
        .sample(
            time,
            &prepared,
            BodyFrameCompactSourceStrength {
                force_on_fluid_body_n: [2.0, 0.0, 0.0],
                ..Default::default()
            },
        )
        .unwrap();
    let (sin, cos) = 0.7_f64.sin_cos();
    close_vector(
        sample.source.strength.force_on_fluid_n,
        [2.0 * cos, 2.0 * sin, 0.0],
    );
    close_vector(
        sample.source.strength.force_rate_on_fluid_n_s,
        [-2.0 * OMEGA * sin, 2.0 * OMEGA * cos, 0.0],
    );
    assert_eq!(sample.kinematics.kinematics.velocity_m_s, [0.0; 3]);
    assert_eq!(sample.kinematics.kinematics.acceleration_m_s2, [0.0; 3]);
    assert_eq!(sample.kinematics.source_time_s, sample.source.source_time_s);
    let pressure = observer([0.0, 6.0, 0.0])
        .prepare_emission(time, sample.kinematics.kinematics)
        .unwrap()
        .evaluate(sample.source.strength)
        .unwrap();
    let near_field = 2.0 * sin / (4.0 * PI * 36.0);
    let radiating = 2.0 * OMEGA * cos / (4.0 * PI * SOUND_SPEED * 6.0);
    close(pressure.pressure_pa, near_field + radiating);
    assert!(radiating.abs() > 0.1 * near_field.abs());
    assert_eq!(pressure.monopole_pressure_pa, 0.0);
}

#[test]
fn changing_world_axes_rotates_geometry_and_force_but_preserves_source_scalars() {
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 0.5), 0.83).unwrap();
    let rotate = |v| rotation.rotate_body_to_world(v);
    let position = Vec3::new(0.4, -0.2, 0.9);
    let velocity = Vec3::new(0.8, -0.6, 0.2);
    let omega = Vec3::new(1.0, -2.0, 3.0);
    let force = Vec3::new(1.5, -2.0, 3.0);
    let gravity = Vec3::new(0.0, -9.0, 0.0);
    let torque = Vec3::new(0.2, 0.4, -0.1);
    let base = RigidBodyIntegrator::new(Gravity::new(gravity).unwrap())
        .prepare_point_motion(
            state(
                Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
                velocity,
                omega,
            ),
            properties(),
            Wrench {
                force_world: force,
                torque_body: torque,
            },
        )
        .unwrap();
    let rotated = RigidBodyIntegrator::new(Gravity::new(rotate(gravity)).unwrap())
        .prepare_point_motion(
            state(
                Pose::new(rotate(position), rotation).unwrap(),
                rotate(velocity),
                omega,
            ),
            properties(),
            Wrench {
                force_world: rotate(force),
                torque_body: torque,
            },
        )
        .unwrap();
    let strength = BodyFrameCompactSourceStrength {
        volume_flow_m3_s: 0.03,
        volume_flow_rate_m3_s2: -0.02,
        force_on_fluid_body_n: [0.6, -0.8, 0.4],
        force_rate_on_fluid_body_n_s: [-0.1, 0.3, 0.2],
    };
    let attachment = RigidAcousticAttachment::new(Vec3::new(0.3, 0.2, -0.1)).unwrap();
    let a = attachment.sample(0.12, &base, strength).unwrap();
    let b = attachment.sample(0.12, &rotated, strength).unwrap();
    for (original, transformed) in [
        (
            a.kinematics.kinematics.position_m,
            b.kinematics.kinematics.position_m,
        ),
        (
            a.kinematics.kinematics.velocity_m_s,
            b.kinematics.kinematics.velocity_m_s,
        ),
        (
            a.kinematics.kinematics.acceleration_m_s2,
            b.kinematics.kinematics.acceleration_m_s2,
        ),
        (
            a.source.strength.force_on_fluid_n,
            b.source.strength.force_on_fluid_n,
        ),
        (
            a.source.strength.force_rate_on_fluid_n_s,
            b.source.strength.force_rate_on_fluid_n_s,
        ),
    ] {
        close_vector(
            transformed,
            components(rotate(Vec3::new(original[0], original[1], original[2]))),
        );
    }
    close_vector(a.source.strength.force_rate_on_fluid_n_s, [1.5, 1.7, 0.6]);
    for sample in [a, b] {
        assert_eq!(sample.source.source_time_s, 0.12);
        assert_eq!(sample.kinematics.source_time_s, 0.12);
        assert_eq!(
            sample.source.strength.volume_flow_m3_s,
            strength.volume_flow_m3_s
        );
        assert_eq!(
            sample.source.strength.volume_flow_rate_m3_s2,
            strength.volume_flow_rate_m3_s2
        );
    }
    // A mechanical load and gravity alone do not invent a force on the fluid.
    let silent = attachment.sample(0.12, &base, Default::default()).unwrap();
    assert_eq!(silent.source.strength, Default::default());
}

fn rotor_history(steps: usize) -> RigidAcousticHistory {
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut accepted = state(Pose::identity(), Vec3::ZERO, Vec3::new(0.0, 0.0, OMEGA));
    let mut result = history(Vec3::new(RADIUS, 0.0, 0.0), steps + 1);
    let dt = DURATION / steps as f64;
    for index in 0..=steps {
        let time = index as f64 * dt;
        let prepared = integrator
            .prepare_point_motion(accepted, properties(), Wrench::ZERO)
            .unwrap();
        result
            .append(
                time,
                &prepared,
                BodyFrameCompactSourceStrength {
                    volume_flow_m3_s: 0.001 + 0.003 * time,
                    volume_flow_rate_m3_s2: 0.003,
                    force_on_fluid_body_n: [1.3, -0.4, 0.2],
                    ..Default::default()
                },
                &mut || false,
            )
            .unwrap();
        if index < steps {
            accepted = integrator
                .step(accepted, properties(), Wrench::ZERO, dt)
                .unwrap()
                .state_after;
        }
    }
    result
}

fn circle_position(time: f64) -> [f64; 3] {
    let (sin, cos) = (OMEGA * time).sin_cos();
    [RADIUS * cos, RADIUS * sin, 0.0]
}

// Independent original circular geometry and scalar bisection: no native
// quaternion, trajectory polynomial, retarded solver or force adapter is used.
fn circular_pressure(receiver_time: f64) -> f64 {
    let mut lo = 0.0;
    let mut hi = DURATION;
    for _ in 0..90 {
        let time = (lo + hi) * 0.5;
        let position = circle_position(time);
        let offset = core::array::from_fn(|i| RECEIVER[i] - position[i]);
        if time + dot(offset, offset).sqrt() / SOUND_SPEED < receiver_time {
            lo = time;
        } else {
            hi = time;
        }
    }
    let time = (lo + hi) * 0.5;
    let (sin, cos) = (OMEGA * time).sin_cos();
    let position = circle_position(time);
    let velocity = [-OMEGA * RADIUS * sin, OMEGA * RADIUS * cos, 0.0];
    let acceleration = [
        -OMEGA * OMEGA * RADIUS * cos,
        -OMEGA * OMEGA * RADIUS * sin,
        0.0,
    ];
    let force = [1.3 * cos + 0.4 * sin, 1.3 * sin - 0.4 * cos, 0.2];
    let derivative = [-OMEGA * force[1], OMEGA * force[0], 0.0];
    let offset = core::array::from_fn(|i| RECEIVER[i] - position[i]);
    let distance = dot(offset, offset).sqrt();
    let n = offset.map(|x| x / distance);
    let mach = velocity.map(|v| v / SOUND_SPEED);
    let radial_mach = dot(n, mach);
    let mach_squared = dot(mach, mach);
    let kappa = 1.0 - radial_mach;
    let radial_acceleration = dot(n, acceleration);
    let monopole = DENSITY
        * (0.003 / (distance * kappa.powi(2))
            + (0.001 + 0.003 * time)
                * (distance * radial_acceleration / SOUND_SPEED
                    + SOUND_SPEED * (radial_mach - mach_squared))
                / (distance * distance * kappa.powi(3)));
    let dipole = dot(n, derivative) / (SOUND_SPEED * distance * kappa.powi(2))
        - dot(force, mach) / (distance * distance * kappa.powi(2))
        + dot(n, force)
            * (1.0 - mach_squared + distance * radial_acceleration / SOUND_SPEED.powi(2))
            / (distance * distance * kappa.powi(3));
    (monopole + dipole) / (4.0 * PI)
}

#[test]
fn accepted_native_rotor_history_converges_on_circular_retarded_pressure() {
    let config = MovingPressureConfig {
        first_receiver_time_s: 0.07,
        sample_rate_hz: 1000,
        samples: 121,
        max_block: 121,
        receiver_bandwidth_hz: 30.0,
        maximum_time_residual_s: 1e-12,
    };
    let mut errors = [0.0_f64; 3];
    for (level, steps) in [12, 24, 48].into_iter().enumerate() {
        let accepted = rotor_history(steps);
        let path = PreparedSourceTrajectory::new(
            observer(RECEIVER),
            accepted.kinematics(),
            trajectory_budget(),
            &mut || false,
        )
        .unwrap();
        let tracks = [CompactSourceHistory {
            motion: &path,
            samples: accepted.sources(),
            source_bandwidth_hz: 4.0,
        }];
        let mut renderer = MovingCompactRenderer::new(
            observer(RECEIVER),
            &tracks,
            config,
            pressure_budget(),
            &mut || false,
        )
        .unwrap();
        let mut pressure = [0.0; 121];
        renderer.block(&mut pressure).unwrap();
        for (index, actual) in pressure.into_iter().enumerate() {
            let time = renderer.receiver_time_at(index as u64).unwrap();
            errors[level] = errors[level].max((actual - circular_pressure(time)).abs());
        }
        for sample in accepted.kinematics() {
            let expected = circle_position(sample.source_time_s);
            close_vector(sample.kinematics.position_m, expected);
            close_vector(
                sample.kinematics.acceleration_m_s2,
                expected.map(|x| -OMEGA * OMEGA * x),
            );
        }
    }
    assert!(
        errors[0] > 1e-10,
        "fixture must resolve interpolation error: {errors:?}"
    );
    assert!(errors[1] < errors[0] / 3.5, "first refinement: {errors:?}");
    assert!(errors[2] < errors[1] / 3.5, "second refinement: {errors:?}");
    assert!(
        errors[2] < 1e-7,
        "finest absolute pressure error [Pa]: {errors:?}"
    );
}

#[test]
fn native_attached_microphone_uses_its_receiver_clock_in_the_shared_renderer() {
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let source = state(Pose::identity(), Vec3::ZERO, Vec3::ZERO);
    let source_motion = integrator
        .prepare_point_motion(source, properties(), Wrench::ZERO)
        .unwrap();
    let mut sources = history(Vec3::ZERO, 3);
    let mut receiver = state(
        Pose::new(Vec3::new(12.0, -2.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(-30.0, 0.0, 0.0),
        Vec3::ZERO,
    );
    let microphone = RigidAcousticAttachment::new(Vec3::new(0.0, 2.0, 0.0)).unwrap();
    let mut receiver_samples = Vec::new();
    for index in 0..3 {
        let time = index as f64 * 0.1;
        sources
            .append(
                time,
                &source_motion,
                BodyFrameCompactSourceStrength {
                    volume_flow_m3_s: 0.001 + 0.003 * time,
                    volume_flow_rate_m3_s2: 0.003,
                    force_on_fluid_body_n: [1.0 + 2.0 * time, 0.2, 0.0],
                    force_rate_on_fluid_body_n_s: [2.0, 0.0, 0.0],
                },
                &mut || false,
            )
            .unwrap();
        let motion = integrator
            .prepare_point_motion(receiver, properties(), Wrench::ZERO)
            .unwrap();
        receiver_samples.push(microphone.kinematics(time, &motion).unwrap());
        if index < 2 {
            receiver = integrator
                .step(receiver, properties(), Wrench::ZERO, 0.1)
                .unwrap()
                .state_after;
        }
    }
    let source_path = PreparedKinematicTrajectory::new(
        medium(),
        sources.kinematics(),
        0.0,
        trajectory_budget(),
        &mut || false,
    )
    .unwrap();
    let receiver_path = PreparedKinematicTrajectory::new(
        medium(),
        &receiver_samples,
        0.2,
        trajectory_budget(),
        &mut || false,
    )
    .unwrap();
    let pair = PreparedMovingReceiverSource::new(
        &source_path,
        &receiver_path,
        MovingReceiverConfig {
            receiver_time_interval_s: [0.05, 0.14],
            minimum_distance_m: 0.1,
            maximum_time_residual_s: 1e-12,
        },
        MovingReceiverBudget {
            maximum_segment_pairs: 16,
            maximum_bound_nodes: 4096,
            maximum_subdivision_depth: 16,
            maximum_root_iterations: 96,
        },
        &mut || false,
    )
    .unwrap();
    let tracks = [CompactSourceHistory {
        motion: &pair,
        samples: sources.sources(),
        source_bandwidth_hz: 1.0,
    }];
    let mut renderer = MovingCompactRenderer::with_moving_receiver(
        &tracks,
        MovingPressureConfig {
            first_receiver_time_s: 0.05,
            sample_rate_hz: 1000,
            samples: 91,
            max_block: 91,
            receiver_bandwidth_hz: 10.0,
            maximum_time_residual_s: 1e-12,
        },
        pressure_budget(),
        &mut || false,
    )
    .unwrap();
    let mut pressure = [0.0; 91];
    renderer.block(&mut pressure).unwrap();
    for (index, actual) in pressure.into_iter().enumerate() {
        let time = renderer.receiver_time_at(index as u64).unwrap();
        let distance = 12.0 - 30.0 * time;
        let tau = time - distance / SOUND_SPEED;
        let expected = (DENSITY * 0.003 / distance
            + (1.0 + 2.0 * tau) / (distance * distance)
            + 2.0 / (SOUND_SPEED * distance))
            / (4.0 * PI);
        close(actual, expected);
        let prepared = pair.prepare_receiver(time, &mut || false).unwrap();
        close(prepared.pressure.source_time_s(), tau);
        close(prepared.observed_source_time_rate, 1.1);
        close_vector(
            prepared.receiver_kinematics.position_m,
            [distance, 0.0, 0.0],
        );
    }
}

#[test]
fn nonfinite_source_inputs_and_transport_overflow_leave_native_state_unchanged() {
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let native = state(Pose::identity(), Vec3::ZERO, Vec3::new(0.0, 0.0, OMEGA));
    let before = native;
    let prepared = integrator
        .prepare_point_motion(native, properties(), Wrench::ZERO)
        .unwrap();
    let attachment = RigidAcousticAttachment::new(Vec3::new(RADIUS, 0.0, 0.0)).unwrap();
    let valid = attachment
        .sample(0.0, &prepared, Default::default())
        .unwrap();
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(RigidAcousticAttachment::new(Vec3::new(bad, 0.0, 0.0)).is_err());
        assert!(attachment.kinematics(bad, &prepared).is_err());
        assert!(
            attachment
                .sample(bad, &prepared, Default::default())
                .is_err()
        );
        for field in 0..8 {
            let mut strength = BodyFrameCompactSourceStrength::default();
            match field {
                0 => strength.volume_flow_m3_s = bad,
                1 => strength.volume_flow_rate_m3_s2 = bad,
                2..=4 => strength.force_on_fluid_body_n[field - 2] = bad,
                _ => strength.force_rate_on_fluid_body_n_s[field - 5] = bad,
            }
            assert!(attachment.sample(0.0, &prepared, strength).is_err());
        }
    }
    assert!(
        attachment
            .sample(
                0.0,
                &prepared,
                BodyFrameCompactSourceStrength {
                    force_on_fluid_body_n: [f64::MAX, 0.0, 0.0],
                    ..Default::default()
                }
            )
            .is_err()
    );
    assert_eq!(
        attachment
            .sample(0.0, &prepared, Default::default())
            .unwrap(),
        valid
    );
    assert_eq!(native, before);
}

#[test]
fn bounded_history_publishes_paired_events_only_after_a_successful_append() {
    let attachment = RigidAcousticAttachment::new(Vec3::ZERO).unwrap();
    for budget in [
        RigidHistoryBudget {
            maximum_samples: 1,
            maximum_payload_bytes: 262_144,
        },
        RigidHistoryBudget {
            maximum_samples: 2,
            maximum_payload_bytes: 0,
        },
        RigidHistoryBudget {
            maximum_samples: usize::MAX,
            maximum_payload_bytes: usize::MAX,
        },
    ] {
        assert!(RigidAcousticHistory::new(attachment, budget).is_err());
    }
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let native = state(Pose::identity(), Vec3::ZERO, Vec3::ZERO);
    let prepared = integrator
        .prepare_point_motion(native, properties(), Wrench::ZERO)
        .unwrap();
    let mut accepted = history(Vec3::ZERO, 3);
    assert!(accepted.is_empty());
    accepted
        .append(0.0, &prepared, Default::default(), &mut || false)
        .unwrap();
    let motion_pointer = accepted.kinematics().as_ptr();
    let source_pointer = accepted.sources().as_ptr();
    let before_motion = accepted.kinematics().to_vec();
    let before_source = accepted.sources().to_vec();
    for time in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            accepted
                .append(time, &prepared, Default::default(), &mut || false)
                .is_err()
        );
    }
    assert!(
        accepted
            .append(
                1.0,
                &prepared,
                BodyFrameCompactSourceStrength {
                    volume_flow_m3_s: f64::NAN,
                    ..Default::default()
                },
                &mut || false
            )
            .is_err()
    );
    for cancel_at in [1, 2] {
        let mut calls = 0;
        let result = accepted.append(1.0, &prepared, Default::default(), &mut || {
            calls += 1;
            calls == cancel_at
        });
        assert!(matches!(
            result,
            Err(RigidAcousticError::Acoustic(AeroacError::Cancelled { .. }))
        ));
        assert_eq!(calls, cancel_at);
        assert_eq!(accepted.kinematics(), before_motion);
        assert_eq!(accepted.sources(), before_source);
        assert_eq!(accepted.len(), 1);
        assert_eq!(accepted.remaining_capacity(), 2);
    }
    for time in [1.0, 2.0] {
        accepted
            .append(time, &prepared, Default::default(), &mut || false)
            .unwrap();
    }
    assert_eq!(accepted.remaining_capacity(), 0);
    assert!(
        accepted
            .append(3.0, &prepared, Default::default(), &mut || false)
            .is_err()
    );
    assert_eq!(accepted.len(), 3);
    assert_eq!(accepted.kinematics().as_ptr(), motion_pointer);
    assert_eq!(accepted.sources().as_ptr(), source_pointer);
    for (motion, source) in accepted.kinematics().iter().zip(accepted.sources()) {
        assert_eq!(motion.source_time_s, source.source_time_s);
    }
}
