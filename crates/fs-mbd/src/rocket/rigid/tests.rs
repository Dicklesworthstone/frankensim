use super::*;
use staging::RocketStage;

fn inertia() -> RigidRocketInertia {
    RigidRocketInertia {
        dry_principal_kg_m2: Vec3::new(2.0, 3.0, 4.0),
        propellant_specific_principal_m2: Vec3::new(0.1, 0.2, 0.2),
    }
}
fn point() -> RocketState {
    RocketState {
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::ZERO,
        dry_mass_kg: 10.0,
        propellant_mass_kg: 4.0,
    }
}
fn state() -> RigidRocketState {
    RigidRocketState::new(point(), UnitQuaternion::IDENTITY, Vec3::ZERO, inertia()).unwrap()
}
fn forcing() -> RigidRocketForcing {
    RigidRocketForcing {
        mass_flow_kg_s: 1.0,
        effective_exhaust_velocity_body_m_s: Vec3::new(100.0, 0.0, 0.0),
        thrust_application_arm_body_m: Vec3::ZERO,
        co_rotating_exhaust_specific_inertia_m2: Vec3::ZERO,
        exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::ZERO,
        gravity_world_m_s2: Vec3::ZERO,
        external_force_world_n: Vec3::ZERO,
        external_torque_body_n_m: Vec3::ZERO,
    }
}
fn controls(step: f64) -> RigidRocketSolve {
    RigidRocketSolve {
        maximum_substep_s: step,
        maximum_substeps: 10000,
    }
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute,
        "{a:.16e} != {b:.16e}, error {}",
        (a - b).abs()
    );
}
fn vector_close(a: Vec3, b: Vec3, absolute: f64) {
    close(a.x, b.x, absolute);
    close(a.y, b.y, absolute);
    close(a.z, b.z, absolute);
}

#[test]
fn g1_native_analytic_translation_survives_body_rotation_and_burnout() {
    let orientation =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), core::f64::consts::FRAC_PI_2)
            .unwrap();
    let body = RigidRocketState::new(point(), orientation, Vec3::ZERO, inertia()).unwrap();
    let input = RigidRocketForcing {
        gravity_world_m_s2: Vec3::new(0.0, 0.0, -2.0),
        external_force_world_n: Vec3::new(3.0, 0.0, 0.0),
        ..forcing()
    };
    let result = body.advance(input, 6.0, controls(6.0), || false).unwrap();
    let expected = point()
        .advance(
            RocketForcing {
                mass_flow_kg_s: input.mass_flow_kg_s,
                effective_exhaust_velocity_m_s: orientation
                    .rotate_body_to_world(input.effective_exhaust_velocity_body_m_s),
                gravity_m_s2: input.gravity_world_m_s2,
                external_force_n: input.external_force_world_n,
            },
            6.0,
            || false,
        )
        .unwrap();
    vector_close(
        result.state_after.translation.velocity_m_s,
        expected.state_after.velocity_m_s,
        2e-13,
    );
    vector_close(
        result.state_after.translation.position_m,
        expected.state_after.position_m,
        2e-13,
    );
    assert_eq!(result.state_after.translation.propellant_mass_kg, 0.0);
    assert_eq!(result.powered_duration_s, 4.0);
    assert_eq!(result.coast_duration_s, 2.0);
    assert_eq!(result.axial_thrust_impulse_n_s, 400.0);
    assert_eq!(result.substeps, 2);
    vector_close(
        result.state_after.angular_momentum_body_kg_m2_s,
        Vec3::ZERO,
        0.0,
    );
}

#[test]
fn g0_original_free_rotation_preserves_world_spin_for_asymmetric_body() {
    let body = RigidRocketState::new(
        point(),
        UnitQuaternion::IDENTITY,
        Vec3::new(1.2, -0.7, 2.3),
        inertia(),
    )
    .unwrap();
    let result = body
        .advance(
            RigidRocketForcing {
                mass_flow_kg_s: 0.0,
                ..forcing()
            },
            20.0,
            controls(0.02),
            || false,
        )
        .unwrap();
    let spin = result
        .state_after
        .orientation
        .rotate_body_to_world(result.state_after.angular_momentum_body_kg_m2_s);
    vector_close(spin, body.angular_momentum_body_kg_m2_s, 3e-11);
    vector_close(
        result.angular_momentum_residual_world_kg_m2_s,
        Vec3::ZERO,
        3e-11,
    );
    close(
        result
            .state_after
            .orientation
            .components()
            .iter()
            .map(|x| x * x)
            .sum::<f64>(),
        1.0,
        1e-15,
    );
}

#[test]
fn g1_corotating_mass_loss_preserves_spin_rate_instead_of_inventing_spin_up() {
    let model = RigidRocketInertia {
        dry_principal_kg_m2: Vec3::new(2.0, 2.0, 2.0),
        propellant_specific_principal_m2: Vec3::new(0.5, 0.5, 0.5),
    };
    let body = RigidRocketState::new(
        point(),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 0.0, 8.0),
        model,
    )
    .unwrap();
    let input = RigidRocketForcing {
        effective_exhaust_velocity_body_m_s: Vec3::ZERO,
        co_rotating_exhaust_specific_inertia_m2: model.propellant_specific_principal_m2,
        ..forcing()
    };
    let result = body.advance(input, 2.0, controls(2.0), || false).unwrap();
    let properties = result.state_after.mass_properties().unwrap();
    close(
        properties
            .angular_velocity_body(result.state_after.angular_momentum_body_kg_m2_s)
            .z,
        2.0,
        2e-14,
    );
    close(
        result.state_after.angular_momentum_body_kg_m2_s.z,
        6.0,
        2e-14,
    );
    vector_close(
        result.exported_angular_momentum_world_kg_m2_s,
        Vec3::new(0.0, 0.0, 2.0),
        2e-14,
    );
    let expected = UnitQuaternion::IDENTITY
        .right_exp(Vec3::new(0.0, 0.0, 4.0))
        .unwrap();
    vector_close(
        result
            .state_after
            .orientation
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        expected.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        3e-14,
    );
    vector_close(
        result.angular_momentum_residual_world_kg_m2_s,
        Vec3::ZERO,
        3e-14,
    );
    let no_flux = body
        .advance(
            RigidRocketForcing {
                co_rotating_exhaust_specific_inertia_m2: Vec3::ZERO,
                ..input
            },
            2.0,
            controls(0.02),
            || false,
        )
        .unwrap();
    close(
        no_flux.state_after.angular_momentum_body_kg_m2_s.z,
        8.0,
        1e-13,
    );
    assert!(
        no_flux
            .state_after
            .mass_properties()
            .unwrap()
            .angular_velocity_body(no_flux.state_after.angular_momentum_body_kg_m2_s)
            .z
            > 2.6
    );
}

#[test]
fn g1_off_axis_thrust_and_external_torque_use_the_original_rotation_map() {
    let model = RigidRocketInertia {
        dry_principal_kg_m2: Vec3::new(2.0, 2.0, 2.0),
        propellant_specific_principal_m2: Vec3::ZERO,
    };
    let body = RigidRocketState::new(point(), UnitQuaternion::IDENTITY, Vec3::ZERO, model).unwrap();
    let input = RigidRocketForcing {
        thrust_application_arm_body_m: Vec3::new(0.0, 0.01, 0.0),
        external_torque_body_n_m: Vec3::new(0.0, 0.0, -1.0),
        ..forcing()
    };
    let result = body.advance(input, 1.0, controls(0.01), || false).unwrap();
    close(
        result.state_after.angular_momentum_body_kg_m2_s.z,
        -2.0,
        5e-14,
    );
    let expected = UnitQuaternion::IDENTITY
        .right_exp(Vec3::new(0.0, 0.0, -0.5))
        .unwrap();
    vector_close(
        result
            .state_after
            .orientation
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        expected.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        2e-13,
    );
    assert!(result.state_after.translation.velocity_m_s.y < 0.0);
    assert!(
        result.thrust_impulse_world_n_s.norm_squared().sqrt() < result.axial_thrust_impulse_n_s
    );
    vector_close(
        result.angular_momentum_residual_world_kg_m2_s,
        Vec3::ZERO,
        2e-13,
    );
}

#[test]
fn g1_rotating_thrust_refines_position_velocity_and_attitude_at_second_order() {
    let body = RigidRocketState::new(
        point(),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.8, -0.4, 1.7),
        inertia(),
    )
    .unwrap();
    let input = RigidRocketForcing {
        thrust_application_arm_body_m: Vec3::new(-0.1, 0.001, 0.002),
        co_rotating_exhaust_specific_inertia_m2: Vec3::new(0.2, 0.4, 0.4),
        external_torque_body_n_m: Vec3::new(0.3, -0.2, 0.1),
        ..forcing()
    };
    let reference = body
        .advance(input, 1.0, controls(1.0 / 2048.0), || false)
        .unwrap()
        .state_after;
    let error = |steps| {
        let actual = body
            .advance(input, 1.0, controls(1.0 / f64::from(steps)), || false)
            .unwrap()
            .state_after;
        [
            actual
                .translation
                .position_m
                .sub(reference.translation.position_m)
                .norm_squared()
                .sqrt(),
            actual
                .translation
                .velocity_m_s
                .sub(reference.translation.velocity_m_s)
                .norm_squared()
                .sqrt(),
            actual
                .orientation
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))
                .sub(
                    reference
                        .orientation
                        .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
                )
                .norm_squared()
                .sqrt(),
        ]
    };
    let a = error(8);
    let b = error(16);
    let c = error(32);
    for i in 0..3 {
        assert!(
            a[i] / b[i] > 3.7 && b[i] / c[i] > 3.7,
            "component {i}: {a:?} {b:?} {c:?}"
        );
    }
}

#[test]
fn g3_world_rotation_and_galilean_shift_preserve_body_motion() {
    let first = RigidRocketState::new(
        point(),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.3, 0.2, 1.0),
        inertia(),
    )
    .unwrap();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap();
    let drift = Vec3::new(11.0, -4.0, 2.0);
    let shift = Vec3::new(100.0, 50.0, -12.0);
    let second = RigidRocketState::new(
        RocketState {
            position_m: shift,
            velocity_m_s: drift,
            ..point()
        },
        rotation,
        first.angular_momentum_body_kg_m2_s,
        inertia(),
    )
    .unwrap();
    let input = RigidRocketForcing {
        gravity_world_m_s2: Vec3::new(0.0, 0.0, -1.0),
        external_force_world_n: Vec3::new(2.0, 1.0, 0.0),
        thrust_application_arm_body_m: Vec3::new(0.0, 0.01, 0.0),
        ..forcing()
    };
    let a = first
        .advance(input, 1.0, controls(0.02), || false)
        .unwrap()
        .state_after;
    let b = second
        .advance(
            RigidRocketForcing {
                gravity_world_m_s2: rotation.rotate_body_to_world(input.gravity_world_m_s2),
                external_force_world_n: rotation.rotate_body_to_world(input.external_force_world_n),
                ..input
            },
            1.0,
            controls(0.02),
            || false,
        )
        .unwrap()
        .state_after;
    vector_close(
        b.translation.position_m,
        shift
            .add(drift)
            .add(rotation.rotate_body_to_world(a.translation.position_m)),
        2e-12,
    );
    vector_close(
        b.translation.velocity_m_s,
        drift.add(rotation.rotate_body_to_world(a.translation.velocity_m_s)),
        2e-12,
    );
    vector_close(
        a.angular_momentum_body_kg_m2_s,
        b.angular_momentum_body_kg_m2_s,
        2e-13,
    );
}

#[test]
fn g4_invalid_inertia_step_budget_and_late_cancellation_publish_no_state() {
    let initial = state();
    assert!(
        RigidRocketState::new(
            point(),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            RigidRocketInertia {
                dry_principal_kg_m2: Vec3::new(1.0, 1.0, 4.0),
                ..inertia()
            }
        )
        .is_err()
    );
    assert!(
        initial
            .advance(
                RigidRocketForcing {
                    co_rotating_exhaust_specific_inertia_m2: Vec3::new(0.0, 0.0, 1.0),
                    ..forcing()
                },
                1.0,
                controls(0.1),
                || false
            )
            .is_err()
    );
    assert!(matches!(
        initial.advance(
            forcing(),
            1.0,
            RigidRocketSolve {
                maximum_substep_s: 0.1,
                maximum_substeps: 2
            },
            || false
        ),
        Err(RigidRocketError::SubstepBudget {
            completed: 2,
            maximum: 2
        })
    ));
    let mut polls = 0;
    assert!(
        initial
            .advance(forcing(), 1.0, controls(0.1), || {
                polls += 1;
                polls > 20
            })
            .is_err()
    );
    assert!(polls > 20);
    assert_eq!(initial, state());
}

fn separation_fixture() -> (RigidRocketState, RocketStage, RocketStage) {
    let piece_inertia = RigidRocketInertia {
        dry_principal_kg_m2: Vec3::new(1.0, 1.0, 1.0),
        propellant_specific_principal_m2: Vec3::ZERO,
    };
    let retained = RocketStage {
        dry_mass_kg: 2.0,
        propellant_mass_kg: 1.0,
        center_of_mass_body_m: Vec3::new(0.0, 0.0, -1.0),
        inertia: piece_inertia,
    };
    let released = RocketStage {
        center_of_mass_body_m: Vec3::new(0.0, 0.0, 1.0),
        ..retained
    };
    let parent = RigidRocketState::new(
        RocketState {
            position_m: Vec3::new(5.0, 7.0, 9.0),
            velocity_m_s: Vec3::new(10.0, -2.0, 4.0),
            dry_mass_kg: 4.0,
            propellant_mass_kg: 2.0,
        },
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.4).unwrap(),
        Vec3::new(1.6, -2.4, 1.4),
        RigidRocketInertia {
            dry_principal_kg_m2: Vec3::new(6.0, 6.0, 2.0),
            propellant_specific_principal_m2: Vec3::new(1.0, 1.0, 0.0),
        },
    )
    .unwrap();
    (parent, retained, released)
}

#[test]
fn g0_staging_preserves_all_momenta_and_accounts_for_separation_work() {
    let (parent, retained, released) = separation_fixture();
    for impulse in [Vec3::ZERO, Vec3::new(1.0, -0.3, 0.5)] {
        let event = parent
            .separate(
                retained,
                released,
                impulse,
                Vec3::new(0.2, 0.1, 0.0),
                1e-12,
                || false,
            )
            .unwrap();
        close(event.mass_residual_kg, 0.0, 0.0);
        vector_close(
            event.linear_momentum_residual_world_kg_m_s,
            Vec3::ZERO,
            3e-14,
        );
        vector_close(
            event.angular_momentum_residual_world_kg_m2_s,
            Vec3::ZERO,
            3e-14,
        );
        close(event.energy_balance_residual_j, 0.0, 3e-14);
        if impulse == Vec3::ZERO {
            assert_eq!(event.separation_work_j, 0.0);
        } else {
            assert!(event.separation_work_j > 0.0);
        }
        assert_eq!(event.retained.translation.propellant_mass_kg, 1.0);
        assert_eq!(event.released.translation.propellant_mass_kg, 1.0);
        assert!(event.retained.translation.velocity_m_s != event.released.translation.velocity_m_s);
    }
}

#[test]
fn g4_staging_refuses_inconsistent_geometry_or_mass_and_is_atomic() {
    let (parent, retained, released) = separation_fixture();
    assert!(
        parent
            .separate(
                RocketStage {
                    dry_mass_kg: 3.0,
                    ..retained
                },
                released,
                Vec3::ZERO,
                Vec3::ZERO,
                1e-12,
                || false
            )
            .is_err()
    );
    assert!(
        parent
            .separate(
                RocketStage {
                    center_of_mass_body_m: Vec3::new(0.1, 0.0, -1.0),
                    ..retained
                },
                released,
                Vec3::ZERO,
                Vec3::ZERO,
                1e-12,
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        parent.separate(retained, released, Vec3::ZERO, Vec3::ZERO, 1e-12, || {
            polls += 1;
            polls == 3
        }),
        Err(RigidRocketError::Cancelled)
    ));
    assert_eq!(parent, separation_fixture().0);
}
