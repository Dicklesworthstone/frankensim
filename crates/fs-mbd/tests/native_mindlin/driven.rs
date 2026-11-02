//! Sustained loading, motor work and recovery of retained Mindlin energy.

use super::*;
use fs_mbd::Wrench;

fn driven(
    law: &SourceHertzMindlinContact,
    state: &SourceHertzMindlinState,
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    external: [Wrench; 2],
    dt: f64,
) -> RigidMindlinCandidate {
    advance_rigid_mindlin_driven(
        bodies,
        properties,
        contact(law, arms),
        state,
        external,
        dt,
        controls(),
        &mut || false,
    )
    .unwrap()
}

fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.sub(expected).norm_squared().sqrt(), 0.0, tolerance);
}

// Independent full-body work/momentum check. No effective mass, reduced free
// displacement or reported energy residual supplies the expected values.
fn native_balance(
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    external: [Wrench; 2],
    dt: f64,
    result: &RigidMindlinCandidate,
) {
    let mut external_work = 0.0;
    let mut contact_work = 0.0;
    let mut kinetic_change = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    for (i, sign) in [(0, -1.0), (1, 1.0)] {
        let old = before[i].point_kinematics(properties[i], arms[i]).unwrap();
        let new = result.bodies_after[i]
            .point_kinematics(properties[i], arms[i])
            .unwrap();
        let mean_v = old
            .center_of_mass_velocity_world
            .add(new.center_of_mass_velocity_world)
            .scale(0.5);
        let mean_w = old
            .angular_velocity_body
            .add(new.angular_velocity_body)
            .scale(0.5);
        external_work +=
            dt * (external[i].force_world.dot(mean_v) + external[i].torque_body.dot(mean_w));
        let impulse = result.report.impulse_on_b_world_ns.scale(sign);
        contact_work += impulse.dot(
            old.point_velocity_world
                .add(new.point_velocity_world)
                .scale(0.5),
        );
        kinetic_change += result.bodies_after[i]
            .kinetic_energy(properties[i])
            .unwrap()
            - before[i].kinetic_energy(properties[i]).unwrap();
        let commanded_linear = external[i].force_world.scale(dt);
        let commanded_angular = external[i].torque_body.scale(dt);
        vector_close(
            result.bodies_after[i]
                .linear_momentum_world()
                .sub(before[i].linear_momentum_world()),
            impulse.add(commanded_linear),
            2e-13,
        );
        vector_close(
            result.bodies_after[i]
                .angular_momentum_body()
                .sub(before[i].angular_momentum_body()),
            arms[i]
                .cross(before[i].pose().orientation().rotate_world_to_body(impulse))
                .add(commanded_angular),
            2e-13,
        );
        vector_close(
            result.report.body_displacement_world_m[i],
            mean_v.scale(dt),
            2e-15,
        );
        linear = linear.add(commanded_linear);
        angular = angular
            .add(before[i].pose().position_world().cross(commanded_linear))
            .add(
                before[i]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(commanded_angular),
            );
    }
    close(result.report.external_wrench_work_j, external_work, 2e-13);
    close(result.report.midpoint_impulse_work_j, contact_work, 2e-13);
    close(kinetic_change - contact_work - external_work, 0.0, 2e-13);
    close(
        kinetic_change
            + result.contact.normal_stored_change_j
            + result.contact.tangential.stored_change_j
            + result.contact.tangential.dissipated_j
            - external_work,
        0.0,
        2e-13,
    );
    vector_close(result.report.external_linear_impulse_ns, linear, 2e-13);
    vector_close(result.report.external_angular_impulse_nms, angular, 2e-13);
    vector_close(result.report.linear_momentum_change_ns, linear, 2e-13);
    vector_close(result.report.angular_momentum_change_nms, angular, 2e-13);
}

#[test]
fn compression_and_motor_torque_drive_partial_slip_from_rest() {
    let law = source(1e7, 0.6);
    let state = law.initial_state();
    let (bodies, properties, arms) = pair(state.indentation_m(), 0.0, 0.0);
    let load = state.normal_force_n();
    let external = [
        Wrench {
            force_world: Vec3::new(2.0, load, 0.0),
            torque_body: Vec3::new(0.0, 0.0, -0.05),
        },
        Wrench {
            force_world: Vec3::new(-0.4, -load, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 0.02),
        },
    ];
    let dt = 1e-3;
    let result = driven(&law, &state, bodies, properties, arms, external, dt);
    close(
        result.contact.tangential.state.indentation_m(),
        state.indentation_m(),
        1e-16,
    );
    close(result.report.normal_average_force_n, load, 1e-13);
    close(result.report.relative_velocity_after_m_s[0], 0.0, 1e-13);
    assert!(result.report.relative_velocity_after_m_s[1] > 0.0);
    assert!(result.report.external_wrench_work_j > 0.0);
    assert!(result.report.kinetic_after_j > result.report.kinetic_before_j);
    assert!(result.contact.tangential.stored_change_j > 0.0);
    assert!(result.contact.tangential.dissipated_j > 0.0);
    // Classical virgin fixed-load Cattaneo--Mindlin curve and its independent
    // displacement integral, before gross sliding begins.
    let ratio = result.contact.tangential.slip_increment_m[0]
        / law.fixed_source().law().gross_slip_displacement_m();
    assert!(ratio > 0.01 && ratio < 1.0, "{ratio:e}");
    let capacity = 0.6 * load;
    let endpoint = capacity * (1.0 - (1.0 - ratio).powf(1.5));
    let average = capacity * (1.0 - 0.4 / ratio * (1.0 - (1.0 - ratio).powf(2.5)));
    close(
        result.contact.tangential.state.endpoint_force_n(),
        endpoint,
        1e-12,
    );
    close(result.contact.tangential.traction_n[0], average, 1e-12);
    native_balance(bodies, properties, arms, external, dt, &result);
}

#[test]
fn held_preload_then_reduced_drive_returns_elastic_work_to_bodies_and_actuator() {
    let law = source(1e7, 0.6);
    let initial = law.initial_state();
    let state = law
        .advance(
            &initial,
            &TangentialStepInput {
                normal_n: initial.normal_force_n(),
                slip_increment_m: [
                    0.9 * law.fixed_source().law().gross_slip_displacement_m(),
                    0.0,
                ],
                dt_s: 1e-3,
                contact_temperature_k: Some(law.temperature_kelvin()),
                nominal_pressure_pa: None,
                compliance_m_per_n: 0.0,
                solve: controls().solve,
            },
        )
        .unwrap()
        .state;
    let (bodies, properties, arms) = pair(state.indentation_m(), 0.0, 0.0);
    let opposing = Vec3::new(state.endpoint_force_n(), state.normal_force_n(), 0.0);
    let held = [0, 1].map(|i| {
        let force = opposing.scale(if i == 0 { 1.0 } else { -1.0 });
        Wrench {
            force_world: force,
            torque_body: arms[i].cross(force),
        }
    });
    let dt = 1e-3;
    let balanced = driven(&law, &state, bodies, properties, arms, held, dt);
    close(balanced.report.relative_displacement_m[0], 0.0, 1e-16);
    close(balanced.report.relative_displacement_m[1], 0.0, 1e-16);
    close(balanced.contact.tangential.dissipated_j, 0.0, 1e-16);
    close(balanced.report.kinetic_after_j, 0.0, 1e-22);
    native_balance(bodies, properties, arms, held, dt, &balanced);
    let reduced = held.map(|w| Wrench {
        force_world: w.force_world.scale(0.25),
        torque_body: w.torque_body.scale(0.25),
    });
    let result = driven(
        &law,
        &balanced.contact.tangential.state,
        balanced.bodies_after,
        properties,
        arms,
        reduced,
        dt,
    );
    assert!(result.contact.normal_stored_change_j < 0.0);
    assert!(result.contact.tangential.stored_change_j < 0.0);
    assert!(result.contact.tangential.dissipated_j >= 0.0);
    assert!(result.report.kinetic_after_j > balanced.report.kinetic_after_j);
    assert!(result.report.external_wrench_work_j < 0.0);
    native_balance(
        balanced.bodies_after,
        properties,
        arms,
        reduced,
        dt,
        &result,
    );
}

#[test]
fn full_inertia_world_rotation_common_acceleration_and_zero_drive_are_consistent() {
    let law = source(1e7, 0.6);
    let state = law.initial_state();
    let (original, _, _) = pair(state.indentation_m(), 0.005, 0.001);
    let arms = [Vec3::new(0.0, 0.2, 0.0), Vec3::new(0.0, -0.2, 0.0)];
    let bodies = [0, 1].map(|i| {
        RigidBodyState::new(
            Pose::new(arms[i].scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
            original[i].linear_momentum_world(),
            Vec3::ZERO,
        )
        .unwrap()
    });
    // Mirrored nonprincipal inertia has real x/z coupling. Equal z impulses
    // give opposing x rotation, retaining zero relative transverse contact slip.
    let properties = [0, 1].map(|i| {
        MassProperties::from_principal_frame(
            (i + 1) as f64,
            Vec3::ZERO,
            Vec3::new(0.012, 0.017, 0.025),
            UnitQuaternion::from_axis_angle(
                Vec3::new(0.0, 1.0, 0.0),
                if i == 0 { 0.37 } else { -0.37 },
            )
            .unwrap(),
        )
        .unwrap()
    });
    let load = state.normal_force_n();
    let external = [
        Wrench {
            force_world: Vec3::new(0.6, load, 0.0),
            torque_body: Vec3::new(0.0, 0.0, -0.01),
        },
        Wrench {
            force_world: Vec3::new(-0.2, -load, 0.0),
            torque_body: Vec3::new(0.0, 0.0, -0.01),
        },
    ];
    let dt = 1e-3;
    let result = driven(&law, &state, bodies, properties, arms, external, dt);
    native_balance(bodies, properties, arms, external, dt, &result);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.73).unwrap();
    let rotated = bodies.map(|body| {
        RigidBodyState::new(
            Pose::new(
                rotation.rotate_body_to_world(body.pose().position_world()),
                rotation,
            )
            .unwrap(),
            rotation.rotate_body_to_world(body.linear_momentum_world()),
            body.angular_momentum_body(),
        )
        .unwrap()
    });
    let rotated_drive = external.map(|w| Wrench {
        force_world: rotation.rotate_body_to_world(w.force_world),
        torque_body: w.torque_body,
    });
    let array = |v: Vec3| [v.x, v.y, v.z];
    let mut rotated_contact = contact(&law, arms);
    rotated_contact.frame = TangentFrame::new(
        array(rotation.rotate_body_to_world(Vec3::new(0.0, 1.0, 0.0))),
        array(rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
    )
    .unwrap();
    let covariant = advance_rigid_mindlin_driven(
        rotated,
        properties,
        rotated_contact,
        &state,
        rotated_drive,
        dt,
        controls(),
        &mut || false,
    )
    .unwrap();
    for i in 0..2 {
        vector_close(
            covariant.bodies_after[i].linear_momentum_world(),
            rotation.rotate_body_to_world(result.bodies_after[i].linear_momentum_world()),
            2e-13,
        );
        vector_close(
            covariant.bodies_after[i].angular_momentum_body(),
            result.bodies_after[i].angular_momentum_body(),
            2e-13,
        );
    }
    close(
        covariant.contact.tangential.dissipated_j,
        result.contact.tangential.dissipated_j,
        1e-15,
    );
    close(
        covariant.report.external_wrench_work_j,
        result.report.external_wrench_work_j,
        1e-13,
    );
    native_balance(rotated, properties, arms, rotated_drive, dt, &covariant);
    let common = Vec3::new(0.2, -0.4, 0.3);
    let accelerated = [0, 1].map(|i| Wrench {
        force_world: external[i]
            .force_world
            .add(common.scale(properties[i].mass())),
        torque_body: external[i].torque_body,
    });
    let changed = driven(&law, &state, bodies, properties, arms, accelerated, dt);
    for i in 0..2 {
        vector_close(
            changed.bodies_after[i].linear_momentum_world(),
            result.bodies_after[i]
                .linear_momentum_world()
                .add(common.scale(dt * properties[i].mass())),
            2e-13,
        );
        vector_close(
            changed.bodies_after[i].angular_momentum_body(),
            result.bodies_after[i].angular_momentum_body(),
            2e-13,
        );
        close(
            changed.report.relative_displacement_m[i],
            result.report.relative_displacement_m[i],
            1e-16,
        );
    }
    close(
        changed.contact.tangential.dissipated_j,
        result.contact.tangential.dissipated_j,
        1e-15,
    );
    native_balance(bodies, properties, arms, accelerated, dt, &changed);
    assert_eq!(
        advance(&law, &state, bodies, properties, arms, dt),
        driven(
            &law,
            &state,
            bodies,
            properties,
            arms,
            [Wrench::ZERO; 2],
            dt
        )
    );
}

#[test]
fn driven_unowned_motion_source_limits_and_late_cancellation_publish_nothing() {
    let law = source(1e7, 0.6);
    let state = law.initial_state();
    let (bodies, properties, arms) = pair(state.indentation_m(), 0.0, 0.0);
    let saved = bodies;
    let saved_state = state.clone();
    for wrench in [
        Wrench {
            force_world: Vec3::new(0.0, 0.0, 1.0),
            torque_body: Vec3::ZERO,
        },
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(0.0, 0.01, 0.0),
        },
        Wrench {
            force_world: Vec3::new(f64::NAN, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
    ] {
        assert!(
            advance_rigid_mindlin_driven(
                bodies,
                properties,
                contact(&law, arms),
                &state,
                [wrench, Wrench::ZERO],
                1e-3,
                controls(),
                &mut || false
            )
            .is_err()
        );
    }
    // Mean normal speed 0.15 is admissible, but the driven native endpoint 0.30
    // exceeds the source's 0.20 m/s envelope and must still refuse.
    let dt = 1e-4;
    let increment = 0.15 * dt;
    let average = law
        .fixed_source()
        .normal_source()
        .elastic_path_average(state.indentation_m(), state.indentation_m() + increment, dt)
        .unwrap();
    let drive_n = 2.0 * increment / (dt * dt) + 1.5 * average.mean_force_n;
    assert!(matches!(
        advance_rigid_mindlin_driven(
            bodies,
            properties,
            contact(&law, arms),
            &state,
            [
                Wrench {
                    force_world: Vec3::new(0.0, drive_n, 0.0),
                    torque_body: Vec3::ZERO
                },
                Wrench::ZERO
            ],
            dt,
            controls(),
            &mut || false
        ),
        Err(RigidMindlinError::Normal(_))
    ));
    let load = state.normal_force_n();
    let external = [
        Wrench {
            force_world: Vec3::new(1.0, load, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench {
            force_world: Vec3::new(0.0, -load, 0.0),
            torque_body: Vec3::ZERO,
        },
    ];
    let mut polls = 0;
    let accepted = advance_rigid_mindlin_driven(
        bodies,
        properties,
        contact(&law, arms),
        &state,
        external,
        1e-3,
        controls(),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 5);
    let mut cancelled_polls = 0;
    assert!(matches!(
        advance_rigid_mindlin_driven(
            bodies,
            properties,
            contact(&law, arms),
            &state,
            external,
            1e-3,
            controls(),
            &mut || {
                cancelled_polls += 1;
                cancelled_polls == polls
            }
        ),
        Err(RigidMindlinError::Cancelled)
    ));
    assert_eq!(bodies, saved);
    assert_eq!(state, saved_state);
    assert_eq!(
        accepted,
        driven(&law, &state, bodies, properties, arms, external, 1e-3)
    );
}
