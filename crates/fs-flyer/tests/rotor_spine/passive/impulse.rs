//! Physical event, work and continuation tests on the retained canonical state.
use super::*;

fn no_impulses<const N: usize>() -> RotorImpulses<N> {
    RotorImpulses {
        carrier: CarrierImpulse::default(),
        motor_nms: [0.0; N],
        aerodynamic_nms: [0.0; N],
    }
}

#[test]
fn full_tensor_event_has_independent_recoil_work_and_fixed_spring_storage() {
    let assembly = RotorAssembly {
        reduced_body: RigidBody::from_inertia_tensor(
            3.0,
            Mat3 {
                m: [2.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 4.0],
            },
        )
        .unwrap(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let initial = assembly
        .state_from_relative_spin(rigid(), [0.37], [0.0])
        .unwrap();
    let mut motion = assembly
        .shaft_motion(
            &initial,
            ShaftProperties {
                stiffness_nm_per_rad: [2.0],
                rest_phase_rad: [0.1],
                damping_nms_per_rad: [0.1],
            },
            0.025,
        )
        .unwrap();
    let input = RotorImpulses {
        carrier: CarrierImpulse {
            linear_world_ns: [3.0, -1.5, 0.75],
            angular_body_nms: [0.1, -0.2, 0.3],
        },
        motor_nms: [0.2],
        aerodynamic_nms: [-0.05],
    };
    let origin = [0.3, 0.4, -0.2];
    let event = motion.apply_impulse(input, origin).unwrap();
    // Carrier receives [0.1,-0.2,0.3] - 0.2*x. Invert the actual
    // [[2,1],[1,3]] xy block and the independent z inertia 4.
    assert!((v(event.state.rigid.omega_body) - v([-0.02, -0.06, 0.075])).norm() < 2e-14);
    assert_eq!(event.state.rigid.vel_mps, [1.0, -0.5, 0.25]);
    assert!((event.state.axial_momentum[0] - 0.15).abs() < 1e-15);
    assert!((event.relative_speed_rad_s[0] - 0.77).abs() < 2e-14);
    assert!((event.motor_work_j - 0.077).abs() < 2e-14);
    assert!((event.external_work_j - 1.96625).abs() < 2e-14);
    assert!((event.kinetic_energy_change_j - 2.04325).abs() < 2e-14);
    assert!((event.elastic_energy_j - 0.0729).abs() < 1e-14);
    assert!(event.work_energy_residual_j.abs() < 2e-14);
    let expected_angular = v([0.05, -0.2, 0.3])
        + (v(initial.rigid.pos_m) - v(origin)).cross(v(input.carrier.linear_world_ns));
    assert!((v(event.external_angular_impulse_world) - expected_angular).norm() < 2e-14);
    let (p0, h0) = total(&assembly, &initial, origin);
    let (p1, h1) = total(&assembly, &event.state, origin);
    assert!((p1 - p0 - v(event.external_linear_impulse_world)).norm() < 1e-14);
    assert!((h1 - h0 - expected_angular).norm() < 2e-14);
    assert_eq!(event.state.rigid.pos_m, initial.rigid.pos_m);
    assert_eq!(event.state.rigid.quat, initial.rigid.quat);
    assert_eq!(event.state.relative_phase_rad, initial.relative_phase_rad);
    assert_eq!(motion.duration_s(), 0.025);
    assert_eq!(motion.state(), &event.state);
    assert_eq!(motion.native_state(), &event.native.state_after);
}

#[test]
fn off_center_event_uses_material_point_work_and_is_world_frame_covariant() {
    let assembly = RotorAssembly {
        reduced_body: RigidBody::from_inertia_tensor(
            3.0,
            Mat3 {
                m: [3.0, 0.2, 0.1, 0.2, 4.0, 0.3, 0.1, 0.3, 5.0],
            },
        )
        .unwrap(),
        rotors: [rotor([1.0, 0.0, 0.0]), rotor([0.0, 1.0, 0.0])],
    };
    let mut carrier = rigid();
    carrier.vel_mps = [2.0, -1.0, 3.0];
    carrier.omega_body = [0.2, -0.3, 0.4];
    set_rotation(
        &mut carrier,
        So3::exp(So3Tangent::new(v([0.3, -0.4, 0.2]))).unwrap(),
    );
    let state = assembly
        .state_from_relative_spin(carrier, [0.1, -0.2], [7.0, -4.0])
        .unwrap();
    let properties = ShaftProperties {
        stiffness_nm_per_rad: [0.8, 1.2],
        rest_phase_rad: [0.2, -0.1],
        damping_nms_per_rad: [0.1, 0.2],
    };
    let mut motion = assembly.shaft_motion(&state, properties, 0.03125).unwrap();
    let before = *motion.state();
    let arm = [0.7, 0.2, -0.3];
    let impulse = [1.2, -0.8, 0.5];
    let origin = [4.0, 5.0, 6.0];
    let input = RotorImpulses {
        carrier: CarrierImpulse::at_body_point(&before.rigid, arm, impulse).unwrap(),
        motor_nms: [0.12, -0.07],
        aerodynamic_nms: [-0.03, 0.02],
    };
    let event = motion.apply_impulse(input, origin).unwrap();
    let r = rotation(&before.rigid);
    let point = v(before.rigid.pos_m) + r.rotate(v(arm)).unwrap();
    let expected_angular =
        (point - v(origin)).cross(v(impulse)) + r.rotate(v([-0.03, 0.02, 0.0])).unwrap();
    assert!((v(event.external_angular_impulse_world) - expected_angular).norm() < 2e-12);
    let mean_omega = (v(before.rigid.omega_body) + v(event.state.rigid.omega_body)).scale(0.5);
    let point_mean_velocity = (v(before.rigid.vel_mps) + v(event.state.rigid.vel_mps)).scale(0.5)
        + r.rotate(mean_omega.cross(v(arm))).unwrap();
    let mut external_work = v(impulse).dot(point_mean_velocity);
    for i in 0..2 {
        external_work += input.aerodynamic_nms[i]
            * (0.5 * before.axial_momentum[i] + 0.5 * event.state.axial_momentum[i])
            / assembly.rotors[i].axial_inertia_kgm2;
    }
    assert!((event.external_work_j - external_work).abs() < 2e-12);
    assert!(event.work_energy_residual_j.abs() < 2e-12);
    assert!(v(event.angular_momentum_residual).norm() < 2e-12);

    let q = So3::exp(So3Tangent::new(v([-0.2, 0.1, 0.3]))).unwrap();
    let mut transformed = before;
    transformed.rigid.pos_m = array(q.rotate(v(before.rigid.pos_m)).unwrap());
    transformed.rigid.vel_mps = array(q.rotate(v(before.rigid.vel_mps)).unwrap());
    set_rotation(&mut transformed.rigid, q.compose(r).unwrap());
    let mut rotated_motion = assembly
        .shaft_motion(&transformed, properties, 0.03125)
        .unwrap();
    let rotated = rotated_motion
        .apply_impulse(
            RotorImpulses {
                carrier: CarrierImpulse::at_body_point(
                    &transformed.rigid,
                    arm,
                    array(q.rotate(v(impulse)).unwrap()),
                )
                .unwrap(),
                ..input
            },
            array(q.rotate(v(origin)).unwrap()),
        )
        .unwrap();
    assert!((v(rotated.state.rigid.omega_body) - v(event.state.rigid.omega_body)).norm() < 2e-12);
    assert!(
        (v(rotated.state.rigid.vel_mps) - q.rotate(v(event.state.rigid.vel_mps)).unwrap()).norm()
            < 2e-12
    );
    assert!(
        (v(rotated.external_angular_impulse_world) - q.rotate(expected_angular).unwrap()).norm()
            < 2e-12
    );
    assert!((rotated.motor_work_j - event.motor_work_j).abs() < 2e-12);
    assert!((rotated.external_work_j - event.external_work_j).abs() < 2e-12);
    assert_eq!(rotated.state.relative_phase_rad, before.relative_phase_rad);
}

#[test]
fn repeated_events_continue_the_same_native_springs_drag_and_canonical_state() {
    let assembly = RotorAssembly {
        reduced_body: RigidBody::from_inertia_tensor(
            3.0,
            Mat3 {
                m: [2.0, 0.3, 0.2, 0.3, 3.0, 0.1, 0.2, 0.1, 4.0],
            },
        )
        .unwrap(),
        rotors: [rotor([1.0, 0.0, 0.0]), rotor([0.0, 1.0, 0.0])],
    };
    let mut carrier = rigid();
    carrier.vel_mps = [0.5, -0.25, 0.125];
    carrier.omega_body = [0.2, -0.1, 0.3];
    let initial = assembly
        .state_from_relative_spin(carrier, [0.4, -0.1], [3.0, -2.0])
        .unwrap();
    let properties = ShaftProperties {
        stiffness_nm_per_rad: [1.7, 0.8],
        rest_phase_rad: [0.1, -0.2],
        damping_nms_per_rad: [0.05, 0.08],
    };
    let dt = 0.03125;
    let mut motion = assembly.shaft_motion(&initial, properties, dt).unwrap();
    let bank = RotorBank::from_reduced(
        assembly.reduced_body.inertia_tensor(),
        assembly
            .rotors
            .map(|r| NativeRotor::from_unit_axis(v(r.axis_body), r.axial_inertia_kgm2).unwrap()),
    )
    .unwrap();
    let native_stepper = bank
        .torsional_stepper(
            properties.stiffness_nm_per_rad,
            properties.rest_phase_rad,
            properties.damping_nms_per_rad,
            dt,
        )
        .unwrap();
    let mut expected_state = *motion.native_state();
    let input = RotorImpulses {
        carrier: CarrierImpulse {
            linear_world_ns: [0.09375, -0.046875, 0.0234375],
            angular_body_nms: [0.0625, -0.03125, 0.015625],
        },
        motor_nms: [0.0625, -0.03125],
        aerodynamic_nms: [-0.015625, 0.0078125],
    };
    let mut heat = 0.0;
    for _ in 0..12 {
        let phase_before = motion.state().relative_phase_rad;
        let expected_event = bank
            .apply_impulse(
                &expected_state,
                [0.046875, -0.0234375],
                v([0.046875, -0.0234375, 0.015625]),
            )
            .unwrap();
        let event = motion.apply_impulse(input, ZERO).unwrap();
        assert_eq!(event.native, expected_event);
        assert_eq!(event.state.relative_phase_rad, phase_before);
        assert_eq!(
            event.elastic_energy_j,
            native_stepper.potential_energy(&expected_state).unwrap()
        );
        assert!(event.work_energy_residual_j.abs() < 1e-12);
        expected_state = expected_event.state_after;

        let expected_step = native_stepper
            .step(
                &expected_state,
                [0.375, -0.1875],
                HeldTorque::Body(v([0.125, -0.0625, 0.0625])),
            )
            .unwrap();
        let step = motion
            .advance(
                Loads {
                    force_n: [0.75, -0.375, 0.1875],
                    moment_nm: [0.25, -0.125, 0.0625],
                },
                RotorTorques {
                    motor_nm: [0.5, -0.25],
                    aerodynamic_nm: [-0.125, 0.0625],
                },
                ZERO,
            )
            .unwrap();
        assert_eq!(step.native, expected_step);
        assert!(step.work_energy_residual_j.abs() < 1e-12);
        heat += step.native.dissipated_energy_j.iter().sum::<f64>();
        expected_state = expected_step.motion.state_after;
    }
    assert!(
        heat > 0.0,
        "the same prepared drag must act after every physical event"
    );
    assert_eq!(motion.duration_s(), dt);
    assert_eq!(motion.native_state(), &expected_state);
}

#[test]
fn event_work_resolves_large_background_energy_without_advancing_phase() {
    let assembly = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let mut carrier = rigid();
    carrier.vel_mps[0] = 1e10;
    let initial = assembly
        .state_from_relative_spin(carrier, [1e30], [0.0])
        .unwrap();
    let mut motion = assembly
        .shaft_motion(&initial, ShaftProperties::default(), 0.05)
        .unwrap();
    let event = motion
        .apply_impulse(
            RotorImpulses {
                carrier: CarrierImpulse {
                    linear_world_ns: [3.0, 0.0, 0.0],
                    ..CarrierImpulse::default()
                },
                ..no_impulses()
            },
            ZERO,
        )
        .unwrap();
    let expected_change = 3.0 * (1e10 + 0.5);
    assert_eq!(event.kinetic_energy_change_j, expected_change);
    assert_eq!(event.external_work_j, expected_change);
    assert_eq!(event.work_energy_residual_j, 0.0);
    let naive_difference = 0.5 * 3.0 * event.state.rigid.vel_mps[0].powi(2)
        - 0.5 * 3.0 * initial.rigid.vel_mps[0].powi(2);
    assert_ne!(naive_difference, expected_change);
    let drive = motion
        .apply_impulse(
            RotorImpulses {
                motor_nms: [0.25],
                ..no_impulses()
            },
            ZERO,
        )
        .unwrap();
    assert!((drive.kinetic_energy_change_j - 0.1875).abs() < 1e-14);
    assert!((drive.motor_work_j - 0.1875).abs() < 1e-14);
    assert!(drive.work_energy_residual_j.abs() < 1e-14);
    assert_eq!(drive.state.rigid.pos_m, initial.rigid.pos_m);
    assert_eq!(drive.state.rigid.quat, initial.rigid.quat);
    assert_eq!(drive.state.relative_phase_rad, [1e30]);
    let accepted = *motion.state();
    let canonical = *motion.native_state();
    assert!(motion.advance(no_loads(), zero_torques(), ZERO).is_err());
    assert_eq!(*motion.state(), accepted);
    assert_eq!(*motion.native_state(), canonical);
}

#[test]
fn unresolved_events_and_late_accounting_refusals_keep_the_accepted_native_state() {
    let assembly = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let mut initial = assembly
        .state_from_relative_spin(rigid(), [0.3], [5.0])
        .unwrap();
    initial.rigid.vel_mps[0] = 1e30;
    let properties = ShaftProperties {
        stiffness_nm_per_rad: [2.0],
        rest_phase_rad: [0.0],
        damping_nms_per_rad: [0.1],
    };
    let mut motion = assembly.shaft_motion(&initial, properties, 0.05).unwrap();
    let mut fresh = assembly.shaft_motion(&initial, properties, 0.05).unwrap();
    let first = motion.advance(no_loads(), zero_torques(), ZERO).unwrap();
    assert_eq!(
        first,
        fresh.advance(no_loads(), zero_torques(), ZERO).unwrap()
    );
    assert!(first.native.dissipated_energy_j[0] > 0.0);
    let accepted = *motion.state();
    let canonical = *motion.native_state();
    let unresolved = RotorImpulses {
        carrier: CarrierImpulse {
            linear_world_ns: [1.0, 0.0, 0.0],
            ..CarrierImpulse::default()
        },
        motor_nms: [0.2],
        ..no_impulses()
    };
    assert_eq!(
        motion.apply_impulse(unresolved, ZERO).unwrap_err().code,
        "impulse-unrepresentable"
    );
    assert_eq!(*motion.state(), accepted);
    assert_eq!(*motion.native_state(), canonical);
    let overflow = RotorImpulses {
        carrier: CarrierImpulse {
            linear_world_ns: [0.0, 3.0, 0.0],
            ..CarrierImpulse::default()
        },
        motor_nms: [0.2],
        ..no_impulses()
    };
    assert!(
        motion
            .apply_impulse(overflow, [f64::MAX, 0.0, 0.0])
            .is_err()
    );
    assert_eq!(*motion.state(), accepted);
    assert_eq!(*motion.native_state(), canonical);
    assert_eq!(
        motion.advance(no_loads(), zero_torques(), ZERO).unwrap(),
        fresh.advance(no_loads(), zero_torques(), ZERO).unwrap(),
    );

    let enormous = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0]); 2],
    };
    let mut state = enormous
        .state_from_relative_spin(rigid(), [0.0; 2], [0.0; 2])
        .unwrap();
    state.axial_momentum[0] = 1e30;
    let mut prepared = enormous
        .shaft_motion(&state, ShaftProperties::default(), 0.05)
        .unwrap();
    let input = RotorImpulses {
        motor_nms: [0.0, 1.0],
        ..no_impulses()
    };
    assert!(prepared.apply_impulse(input, ZERO).is_err());
    assert_eq!(*prepared.state(), state);
    // The distinct velocity-based impulse API retains its wider high-spin domain.
    assert_eq!(
        apply_impulse(&enormous, &state, input, ZERO)
            .unwrap()
            .state
            .rigid
            .omega_body[0],
        -1.0
    );
    assert_eq!(
        prepared.apply_impulse(no_impulses(), ZERO).unwrap().state,
        state
    );
}
