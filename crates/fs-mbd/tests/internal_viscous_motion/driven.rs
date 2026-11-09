//! G1 moving driven-damped trajectories; G3 frame and source limits.

use super::*;
use fs_mbd::spatial_loads::internal_motion::{
    PrescribedInternalWrench, step_internal_wrench,
    viscous::driven::{
        DrivenViscousMotionError, HeldInternalEffort, step_driven_viscous_connector,
        step_driven_viscous_connector_with_loads, step_driven_viscous_connector_with_spatial_loads,
    },
};

fn near_states(actual: [RigidBodyState; 2], expected: [RigidBodyState; 2], tolerance: f64) {
    for i in 0..2 {
        near_vec(
            actual[i].pose().position_world(),
            expected[i].pose().position_world(),
            tolerance,
        );
        near_vec(
            actual[i].linear_momentum_world(),
            expected[i].linear_momentum_world(),
            tolerance,
        );
        near_vec(
            actual[i].angular_momentum_body(),
            expected[i].angular_momentum_body(),
            tolerance,
        );
        for (a, b) in actual[i]
            .pose()
            .orientation()
            .components()
            .into_iter()
            .zip(expected[i].pose().orientation().components())
        {
            near(a, b, tolerance);
        }
    }
}

#[test]
fn held_force_has_exact_terminal_slip_work_heat_and_refining_position() {
    let (initial, properties, law) = axial();
    let duration: f64 = 0.7;
    let force = 1.7;
    let damping = 2.3;
    let rate = damping * (1.0 / 2.0 + 1.0 / 3.0);
    let terminal = force / damping;
    let transient = -1.5 - terminal;
    let integrated_decay = -(-rate * duration).exp_m1() / rate;
    let expected_slip = terminal + transient * (-rate * duration).exp();
    let travel = terminal * duration + transient * integrated_decay;
    let work = force * travel;
    let heat = damping
        * (terminal * terminal * duration
            + 2.0 * terminal * transient * integrated_decay
            + transient * transient * -(-2.0 * rate * duration).exp_m1() / (2.0 * rate));
    let effort = HeldInternalEffort::World([force, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let mut errors = [0.0; 3];
    for (rung, steps) in [1, 2, 4].into_iter().enumerate() {
        let mut states = initial;
        let mut observed_heat = 0.0;
        let mut observed_work = 0.0;
        for _ in 0..steps {
            let result = step_driven_viscous_connector(
                states,
                properties,
                &law,
                effort,
                duration / steps as f64,
                controls(),
                || false,
            )
            .unwrap();
            near(result.work_heat_residual_j, 0.0, 5e-13);
            near_vec(result.linear_momentum_residual_ns, Vec3::ZERO, 5e-14);
            near_vec(result.angular_momentum_residual_nms, Vec3::ZERO, 5e-14);
            states = result.states_after;
            observed_heat += result.dissipated_energy_j;
            observed_work += result.actuator_work_j;
        }
        let velocity = core::array::from_fn::<_, 2, _>(|i| {
            states[i]
                .center_of_mass_velocity_world(properties[i])
                .unwrap()
                .x
        });
        near(velocity[1] - velocity[0], expected_slip, 2e-13);
        near((2.0 * velocity[0] + 3.0 * velocity[1]) / 5.0, 0.3, 2e-14);
        near(observed_heat, heat, 8e-13);
        near(observed_work, work, 8e-13);
        assert!(observed_heat > 0.0 && observed_work < 0.0);
        errors[rung] = (states[1].pose().position_world().x
            - states[0].pose().position_world().x
            - 2.0
            - travel)
            .abs();
    }
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.26 * errors[1], "{errors:?}");

    let long =
        step_driven_viscous_connector(initial, properties, &law, effort, 12.0, controls(), || {
            false
        })
        .unwrap();
    let slip = long.states_after[1]
        .center_of_mass_velocity_world(properties[1])
        .unwrap()
        .x
        - long.states_after[0]
            .center_of_mass_velocity_world(properties[0])
            .unwrap()
            .x;
    near(slip, terminal + transient * (-rate * 12.0).exp(), 2e-12);
    near(long.work_heat_residual_j, 0.0, 2e-11);
}

// Extend the existing independent continuous Newton-Euler derivative with
// physical actuator forces/couples. Neither impulse nor source kernels are used.
#[derive(Clone, Copy)]
struct DrivenReference {
    motion: Reference,
    work: f64,
}
impl DrivenReference {
    fn shifted(self, slope: Self, dt: f64) -> Self {
        Self {
            motion: self.motion.shifted(slope.motion, dt),
            work: self.work + dt * slope.work,
        }
    }
}
fn driven_derivative(
    s: DrivenReference,
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    c: [[f64; 6]; 6],
    effort: HeldInternalEffort,
) -> DrivenReference {
    let mut motion = derivative(s.motion, properties, law, c);
    let q = s
        .motion
        .q
        .map(|q| UnitQuaternion::new(q[0], q[1], q[2], q[3]).unwrap());
    let authored = match effort {
        HeldInternalEffort::MaterialFirst(value) | HeldInternalEffort::World(value) => value,
    };
    let mut force = Vec3::new(authored[0], authored[1], authored[2]);
    let mut couple = Vec3::new(authored[3], authored[4], authored[5]);
    if matches!(effort, HeldInternalEffort::MaterialFirst(_)) {
        force = q[0].rotate_body_to_world(force);
        couple = q[0].rotate_body_to_world(couple);
    }
    let arms = core::array::from_fn::<_, 2, _>(|i| q[i].rotate_body_to_world(law.arms_body_m[i]));
    let separation = s.motion.x[1].add(arms[1]).sub(s.motion.x[0].add(arms[0]));
    let torque = [
        arms[0]
            .cross(force)
            .scale(-1.0)
            .sub(couple)
            .sub(separation.cross(force)),
        arms[1].cross(force).add(couple),
    ];
    let omega = core::array::from_fn::<_, 2, _>(|i| {
        q[i].rotate_body_to_world(
            properties[i]
                .angular_velocity_body_checked(s.motion.h[i])
                .unwrap(),
        )
    });
    let relative = motion.x[1].add(omega[1].cross(arms[1])).sub(
        motion.x[0]
            .add(omega[0].cross(arms[0]))
            .add(omega[0].cross(separation)),
    );
    motion.p[0] = motion.p[0].sub(force);
    motion.p[1] = motion.p[1].add(force);
    for i in 0..2 {
        motion.h[i] = motion.h[i].add(q[i].rotate_world_to_body(torque[i]));
    }
    DrivenReference {
        motion,
        work: force.dot(relative) + couple.dot(omega[1].sub(omega[0])),
    }
}
fn driven_reference(
    initial: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    c: [[f64; 6]; 6],
    effort: HeldInternalEffort,
    duration: f64,
) -> DrivenReference {
    let mut s = DrivenReference {
        motion: Reference::from_states(initial),
        work: 0.0,
    };
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let k1 = driven_derivative(s, properties, law, c, effort);
        let k2 = driven_derivative(s.shifted(k1, dt / 2.0), properties, law, c, effort);
        let k3 = driven_derivative(s.shifted(k2, dt / 2.0), properties, law, c, effort);
        let k4 = driven_derivative(s.shifted(k3, dt), properties, law, c, effort);
        s = s
            .shifted(k1, dt / 6.0)
            .shifted(k2, dt / 3.0)
            .shifted(k3, dt / 3.0)
            .shifted(k4, dt / 6.0);
    }
    s
}

#[test]
fn both_actuator_frames_refine_to_independent_full_spatial_dynamics() {
    let (initial, properties, law, c) = coupled();
    let command = [0.35, -0.2, 0.15, 0.08, -0.12, 0.1];
    let mut endpoints = [initial; 2];
    for (frame, effort) in [
        HeldInternalEffort::MaterialFirst(command),
        HeldInternalEffort::World(command),
    ]
    .into_iter()
    .enumerate()
    {
        let expected = driven_reference(initial, properties, &law, c, effort, 0.4);
        let mut errors = [0.0; 3];
        for (rung, steps) in [2, 4, 8].into_iter().enumerate() {
            let mut states = initial;
            let mut heat = 0.0;
            let mut work = 0.0;
            for _ in 0..steps {
                let result = step_driven_viscous_connector(
                    states,
                    properties,
                    &law,
                    effort,
                    0.4 / steps as f64,
                    controls(),
                    || false,
                )
                .unwrap();
                near(result.work_heat_residual_j, 0.0, 2e-11);
                near_vec(result.linear_momentum_residual_ns, Vec3::ZERO, 2e-12);
                near_vec(result.angular_momentum_residual_nms, Vec3::ZERO, 2e-11);
                assert!(result.dissipated_energy_j >= 0.0);
                let endpoint_poses = result.half_sources[1]
                    .impulse
                    .bodies
                    .map(|b| b.state_before.pose());
                assert_eq!(endpoint_poses, result.states_after.map(|s| s.pose()));
                assert_ne!(endpoint_poses, states.map(|s| s.pose()));
                for source in result.half_sources {
                    let mut components = command;
                    if matches!(effort, HeldInternalEffort::World(_)) {
                        let q = source.impulse.bodies[0].state_before.pose().orientation();
                        let f =
                            q.rotate_world_to_body(Vec3::new(command[0], command[1], command[2]));
                        let t =
                            q.rotate_world_to_body(Vec3::new(command[3], command[4], command[5]));
                        components = [f.x, f.y, f.z, t.x, t.y, t.z];
                    }
                    let power: f64 = components
                        .into_iter()
                        .zip(source.average_relative_twist)
                        .map(|(force, speed)| force * speed)
                        .sum();
                    near(source.actuator_work_j, power * 0.2 / steps as f64, 3e-13);
                }
                states = result.states_after;
                heat += result.dissipated_energy_j;
                work += result.actuator_work_j;
            }
            errors[rung] =
                reference_error(states, heat, expected.motion) + (work - expected.work).abs();
            endpoints[frame] = states;
        }
        assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.3 * errors[1], "{errors:?}");
    }
    assert!(
        norm(
            endpoints[0][1]
                .linear_momentum_world()
                .sub(endpoints[1][1].linear_momentum_world())
        ) > 1e-3
    );
}

#[test]
fn zero_damping_and_zero_actuation_recover_original_motion_owners() {
    let (initial, properties, mut law, _) = coupled();
    let passive =
        step_viscous_connector(initial, properties, &law, 0.1, controls(), || false).unwrap();
    let unforced = step_driven_viscous_connector(
        initial,
        properties,
        &law,
        HeldInternalEffort::MaterialFirst([0.0; 6]),
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    near_states(unforced.states_after, passive.states_after, 2e-12);
    near(
        unforced.dissipated_energy_j,
        passive.dissipated_energy_j,
        2e-12,
    );
    assert_eq!(unforced.actuator_work_j, 0.0);

    law.damping = matrix([[0.0; 6]; 6]);
    let f = Vec3::new(0.2, -0.1, 0.3);
    let t = Vec3::new(-0.05, 0.08, 0.02);
    let native = step_internal_wrench(
        initial,
        properties,
        PrescribedInternalWrench {
            arms_body_m: law.arms_body_m,
            force_on_second_world_n: f,
            couple_on_second_world_nm: t,
        },
        0.1,
        controls().impulse,
        || false,
    )
    .unwrap();
    let driven = step_driven_viscous_connector(
        initial,
        properties,
        &law,
        HeldInternalEffort::World([f.x, f.y, f.z, t.x, t.y, t.z]),
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    near_states(driven.states_after, native.states_after, 2e-12);
    near(driven.actuator_work_j, native.actuator_work_j, 2e-12);
    assert_eq!(driven.dissipated_energy_j, 0.0);
}

#[test]
fn gravity_and_external_work_stay_separate_from_internal_drive_and_heat() {
    let (initial, properties, law) = axial();
    let effort = HeldInternalEffort::World([1.7, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let duration = 0.1;
    let acceleration = Vec3::new(0.2, 0.1, -0.3);
    let gravity = Gravity::new(Vec3::new(0.0, -9.81, 0.0)).unwrap();
    let wrenches = properties.map(|p| Wrench {
        force_world: acceleration.scale(p.mass()),
        torque_body: Vec3::ZERO,
    });
    let free = step_driven_viscous_connector(
        initial,
        properties,
        &law,
        effort,
        duration,
        controls(),
        || false,
    )
    .unwrap();
    let loaded = step_driven_viscous_connector_with_loads(
        initial,
        properties,
        &law,
        effort,
        HeldExternalLoads {
            wrenches: &wrenches,
            gravity,
        },
        duration,
        controls(),
        || false,
    )
    .unwrap();
    near(loaded.dissipated_energy_j, free.dissipated_energy_j, 3e-13);
    near(loaded.actuator_work_j, free.actuator_work_j, 3e-13);
    let total_acceleration = acceleration.add(gravity.acceleration_world());
    let mut force_work = 0.0;
    for i in 0..2 {
        near_vec(
            loaded.states_after[i].linear_momentum_world(),
            free.states_after[i]
                .linear_momentum_world()
                .add(total_acceleration.scale(duration * properties[i].mass())),
            3e-13,
        );
        near_vec(
            loaded.states_after[i].pose().position_world(),
            free.states_after[i]
                .pose()
                .position_world()
                .add(total_acceleration.scale(0.5 * duration * duration)),
            3e-13,
        );
        force_work += wrenches[i].force_world.dot(
            loaded.states_after[i]
                .pose()
                .position_world()
                .sub(initial[i].pose().position_world()),
        );
    }
    near(
        loaded.external_loads.applied_force_work_j,
        force_work,
        3e-13,
    );
    near(
        loaded.gravitational_potential_energy_change_j,
        -loaded.external_loads.gravity_work_j,
        0.0,
    );
    near(
        loaded.mechanical_energy_change_j + loaded.dissipated_energy_j,
        loaded.actuator_work_j + loaded.external_loads.applied_work_j,
        3e-13,
    );
}

#[test]
fn rotating_material_point_loads_keep_their_actual_source_work() {
    let (initial, properties, law, _) = coupled();
    let points = [
        WorldPointForce {
            arm_body: Vec3::new(0.2, -0.1, 0.15),
            force_world: Vec3::new(0.1, 0.3, -0.2),
        },
        WorldPointForce {
            arm_body: Vec3::new(-0.3, 0.2, 0.1),
            force_world: Vec3::new(-0.2, 0.1, 0.05),
        },
    ];
    let couples = [Vec3::new(0.03, -0.02, 0.04), Vec3::new(-0.01, 0.03, 0.02)];
    let loads = core::array::from_fn::<_, 2, _>(|i| {
        PreparedSpatialLoads::new(Wrench::ZERO, couples[i], &[points[i]], || false).unwrap()
    });
    let result = step_driven_viscous_connector_with_spatial_loads(
        initial,
        properties,
        &law,
        HeldInternalEffort::MaterialFirst([0.35, -0.2, 0.15, 0.08, -0.12, 0.1]),
        HeldSpatialLoads {
            bodies: &loads,
            gravity: Gravity::ZERO,
            midpoint: SpatialMidpointConfig {
                maximum_iterations: 24,
                maximum_backtracks: 12,
                angular_momentum_tolerance_nms: 1e-13,
                relative_tolerance: 1e-13,
                maximum_cayley_increment: 0.5,
            },
        },
        0.07,
        controls(),
        || false,
    )
    .unwrap();
    let mut work = 0.0;
    for i in 0..2 {
        let before = result.half_sources[0].impulse.bodies[i].state_after;
        let after = result.half_sources[1].impulse.bodies[i].state_before;
        let point = |s: RigidBodyState| {
            s.pose().position_world().add(
                s.pose()
                    .orientation()
                    .rotate_body_to_world(points[i].arm_body),
            )
        };
        work += points[i].force_world.dot(point(after).sub(point(before)));
    }
    near(result.external_loads.applied_force_work_j, work, 3e-13);
    near(
        result.mechanical_energy_change_j + result.dissipated_energy_j,
        result.actuator_work_j + result.external_loads.applied_work_j,
        2e-11,
    );
    assert!(result.dissipated_energy_j > 0.0);
}

#[test]
fn source_motion_and_final_cancellation_publish_nothing_and_retry_exactly() {
    let (initial, properties, law, _) = coupled();
    let saved_law = law.clone();
    let effort = HeldInternalEffort::World([0.35, -0.2, 0.15, 0.08, -0.12, 0.1]);
    let mut count = 0;
    let expected =
        step_driven_viscous_connector(initial, properties, &law, effort, 0.1, controls(), || {
            count += 1;
            false
        })
        .unwrap();
    for target in [1, count / 2, count] {
        let mut polls = 0;
        let refused = step_driven_viscous_connector(
            initial,
            properties,
            &law,
            effort,
            0.1,
            controls(),
            || {
                polls += 1;
                polls == target
            },
        );
        assert!(refused.is_err());
        assert_eq!(polls, target);
        if target == count {
            assert!(matches!(
                refused,
                Err(DrivenViscousMotionError::Motion(
                    InternalWrenchStepError::Cancelled
                ))
            ));
        }
    }
    for duration in [0.0, -0.1, f64::NAN, f64::from_bits(1)] {
        assert!(
            step_driven_viscous_connector(
                initial,
                properties,
                &law,
                effort,
                duration,
                controls(),
                || false,
            )
            .is_err()
        );
    }
    assert!(
        step_driven_viscous_connector(
            initial,
            properties,
            &law,
            HeldInternalEffort::World([f64::NAN; 6]),
            0.1,
            controls(),
            || false,
        )
        .is_err()
    );
    assert!(matches!(
        step_driven_viscous_connector_with_loads(
            initial,
            properties,
            &law,
            effort,
            HeldExternalLoads {
                wrenches: &[Wrench::ZERO],
                gravity: Gravity::ZERO
            },
            0.1,
            controls(),
            || false,
        ),
        Err(DrivenViscousMotionError::Motion(
            InternalWrenchStepError::ExternalLoadCount { .. }
        ))
    ));
    assert_eq!(law, saved_law);
    assert_eq!(step_driven_viscous_connector(
        initial, properties, &law, effort, 0.1, controls(), || false,
    ).unwrap(), expected);
}
