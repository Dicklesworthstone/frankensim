//! Physical integration oracles for the original native shaft owner and Flyer state.
use super::*;
use fs_flyer::spine::rotors::passive::ShaftProperties;
use fs_time::lie::gyrostat::{
    HeldTorque,
    rotors::{AxialRotor as NativeRotor, RotorBank},
};

fn zero_torques<const N: usize>() -> RotorTorques<N> {
    RotorTorques {
        motor_nm: [0.0; N],
        aerodynamic_nm: [0.0; N],
    }
}

#[test]
fn g1_native_bearing_drag_recoils_the_carrier_and_closes_actual_heat() {
    let assembly = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let initial = assembly
        .state_from_relative_spin(rigid(), [0.0], [5.0])
        .unwrap();
    let dt = 0.05;
    let damping = 0.2;
    let mut motion = assembly
        .shaft_motion(
            &initial,
            ShaftProperties {
                damping_nms_per_rad: [damping],
                ..ShaftProperties::default()
            },
            dt,
        )
        .unwrap();
    let step = motion.advance(no_loads(), zero_torques(), ZERO).unwrap();
    // Two native midpoint half-stages solve ds/dt = -d*(1/J+1/I)*s.
    let mobility = 1.0 / 0.2 + 1.0;
    let ratio = (1.0 - 0.25 * dt * damping * mobility) / (1.0 + 0.25 * dt * damping * mobility);
    let speed = 5.0 * ratio * ratio;
    let momentum = (speed + 1.0) / mobility;
    let omega = 1.0 - momentum;
    let phase = 0.25 * dt * 5.0 * (1.0 + ratio) * (1.0 + ratio);
    let energy_after = 0.5 * omega * omega + 0.5 * momentum * momentum / 0.2;
    assert!((step.relative_speed_rad_s[0] - speed).abs() < 1e-13);
    assert!((step.state.axial_momentum[0] - momentum).abs() < 1e-14);
    assert!((step.state.rigid.omega_body[0] - omega).abs() < 1e-14);
    assert!((step.state.relative_phase_rad[0] - phase).abs() < 1e-14);
    assert!(
        step.state.rigid.omega_body[0] > 0.0,
        "bearing reaction must move the carrier"
    );
    assert!((step.native.dissipated_energy_j[0] - (2.5 - energy_after)).abs() < 1e-13);
    assert!(step.native.dissipated_energy_j[0] > 0.0);
    assert_eq!(step.native.elastic_energy_change_j, [0.0]);
    assert!(step.work_energy_residual_j.abs() < 1e-13);
    assert!(
        (total(&assembly, &step.state, ZERO).1 - total(&assembly, &initial, ZERO).1).norm() < 1e-13
    );
    assert_eq!(motion.state(), &step.state);
    assert_eq!(motion.native_state(), &step.native.motion.state_after);
}

#[test]
fn g1_torsional_state_and_storage_match_the_coupled_two_inertia_oscillator() {
    let assembly = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let initial = assembly
        .state_from_relative_spin(rigid(), [0.4], [-0.2])
        .unwrap();
    let dt = 0.05;
    let k = 2.0;
    let d = 0.1;
    let mut motion = assembly
        .shaft_motion(
            &initial,
            ShaftProperties {
                stiffness_nm_per_rad: [k],
                rest_phase_rad: [0.1],
                damping_nms_per_rad: [d],
            },
            dt,
        )
        .unwrap();
    let mut strain = 0.3;
    let mut speed = -0.2;
    let mobility = 6.0;
    for _ in 0..2 {
        let h = 0.5 * dt;
        let next_speed = ((1.0 - 0.5 * h * mobility * d - 0.25 * h * h * mobility * k) * speed
            - h * mobility * k * strain)
            / (1.0 + 0.5 * h * mobility * d + 0.25 * h * h * mobility * k);
        strain += 0.5 * h * (speed + next_speed);
        speed = next_speed;
    }
    let step = motion.advance(no_loads(), zero_torques(), ZERO).unwrap();
    assert!((step.relative_speed_rad_s[0] - speed).abs() < 1e-13);
    assert!((step.state.relative_phase_rad[0] - (0.1 + strain)).abs() < 1e-14);
    assert!((step.elastic_energy_before_j - 0.5 * k * 0.3 * 0.3).abs() < 1e-14);
    assert!((step.elastic_energy_after_j - 0.5 * k * strain * strain).abs() < 1e-14);
    assert!(step.native.elastic_energy_change_j[0] < 0.0);
    assert!(step.native.dissipated_energy_j[0] > 0.0);
    assert!(step.work_energy_residual_j.abs() < 1e-13);
    assert!(
        (total(&assembly, &step.state, ZERO).1 - total(&assembly, &initial, ZERO).1).norm() < 1e-13
    );
}

#[test]
fn g3_full_tensor_flyer_motion_retains_native_momentum_and_stage_receipts() {
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
    carrier.vel_mps = [1.0, -0.5, 0.25];
    carrier.omega_body = [0.2, -0.1, 0.3];
    set_rotation(
        &mut carrier,
        So3::exp(So3Tangent::new(v([0.1, 0.2, -0.1]))).unwrap(),
    );
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
    let native = bank
        .torsional_stepper(
            properties.stiffness_nm_per_rad,
            properties.rest_phase_rad,
            properties.damping_nms_per_rad,
            dt,
        )
        .unwrap();
    let mut expected_state = *motion.native_state();
    let torques = RotorTorques {
        motor_nm: [0.5, -0.25],
        aerodynamic_nm: [-0.125, 0.0625],
    };
    let loads = Loads {
        force_n: [0.75, -0.375, 0.1875],
        moment_nm: [0.25, -0.125, 0.0625],
    };
    let origin = [0.5, -0.25, 0.75];
    for _ in 0..12 {
        let before = *motion.state();
        let expected = native
            .step(
                &expected_state,
                [0.375, -0.1875],
                HeldTorque::Body(v([0.125, -0.0625, 0.0625])),
            )
            .unwrap();
        let accepted = motion.advance(loads, torques, origin).unwrap();
        assert_eq!(accepted.native, expected);
        assert_eq!(
            accepted.state.axial_momentum,
            expected.motion.state_after.axial_momentum_nms
        );
        assert_eq!(
            accepted.state.rigid.omega_body,
            array(bank.angular_velocity(&expected.motion.state_after).unwrap())
        );
        let (p0, h0) = total(&assembly, &before, origin);
        let (p1, h1) = total(&assembly, &accepted.state, origin);
        assert!((p1 - p0 - v(accepted.external_linear_impulse_world)).norm() < 1e-12);
        assert!((h1 - h0 - v(accepted.external_angular_impulse_world)).norm() < 1e-11);
        assert!(v(accepted.angular_momentum_residual).norm() < 1e-11);
        assert!(accepted.work_energy_residual_j.abs() < 1e-12);
        expected_state = expected.motion.state_after;
    }
    assert_ne!(motion.state().rigid.pos_m, initial.rigid.pos_m);
}

#[test]
fn g1_native_motor_and_aerodynamic_ports_keep_distinct_recoil_and_work() {
    let assembly = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let initial = assembly
        .state_from_relative_spin(rigid(), [0.0], [5.0])
        .unwrap();
    let dt = 0.05;
    let motor = 0.25;
    let aero = -0.125;
    let mut motion = assembly
        .shaft_motion(&initial, ShaftProperties::default(), dt)
        .unwrap();
    let step = motion
        .advance(
            no_loads(),
            RotorTorques {
                motor_nm: [motor],
                aerodynamic_nm: [aero],
            },
            ZERO,
        )
        .unwrap();
    assert!((step.state.rigid.omega_body[0] + motor * dt).abs() < 1e-14);
    assert!((step.state.axial_momentum[0] - (1.0 + (motor + aero) * dt)).abs() < 1e-14);
    assert!((step.external_angular_impulse_world[0] - aero * dt).abs() < 1e-14);
    assert!((step.motor_work_j - motor * step.state.relative_phase_rad[0]).abs() < 1e-14);
    let absolute_mean_speed =
        (initial.axial_momentum[0] + step.state.axial_momentum[0]) / (2.0 * 0.2);
    assert!((step.external_work_j - aero * dt * absolute_mean_speed).abs() < 1e-13);
    assert!(step.work_energy_residual_j.abs() < 1e-13);
}

#[test]
fn g0_native_admission_and_late_refusal_preserve_the_complete_accepted_state() {
    let assembly = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let mut initial = assembly
        .state_from_relative_spin(rigid(), [0.0], [5.0])
        .unwrap();
    let properties = ShaftProperties {
        damping_nms_per_rad: [0.1],
        ..ShaftProperties::default()
    };
    initial.rigid.vel_mps[0] = 1e30;
    let mut motion = assembly.shaft_motion(&initial, properties, 0.05).unwrap();
    let mut fresh = assembly.shaft_motion(&initial, properties, 0.05).unwrap();
    let before = *motion.state();
    let canonical = *motion.native_state();
    let error = motion
        .advance(
            Loads {
                force_n: [1.0, 0.0, 0.0],
                moment_nm: ZERO,
            },
            zero_torques(),
            ZERO,
        )
        .unwrap_err();
    assert_eq!(error.code, "impulse-unrepresentable");
    assert_eq!(*motion.state(), before);
    assert_eq!(*motion.native_state(), canonical);
    assert_eq!(
        motion.advance(no_loads(), zero_torques(), ZERO).unwrap(),
        fresh.advance(no_loads(), zero_torques(), ZERO).unwrap(),
    );

    initial.rigid.vel_mps = ZERO;
    initial.rigid.omega_body[0] = 1.0;
    initial.axial_momentum[0] = 1e30;
    assert!(assembly.shaft_motion(&initial, properties, 0.05).is_err());

    let regular = assembly
        .state_from_relative_spin(rigid(), [0.0], [5.0])
        .unwrap();
    assert!(
        assembly
            .shaft_motion(
                &regular,
                ShaftProperties {
                    stiffness_nm_per_rad: [1e12],
                    ..ShaftProperties::default()
                },
                0.05,
            )
            .is_err()
    );
    let at_cap = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0]); 64],
    };
    let state = at_cap
        .state_from_relative_spin(rigid(), [0.0; 64], [0.0; 64])
        .unwrap();
    let mut capped = at_cap
        .shaft_motion(&state, ShaftProperties::default(), 0.01)
        .unwrap();
    assert_eq!(
        capped
            .advance(no_loads(), zero_torques(), ZERO)
            .unwrap()
            .state,
        state
    );
}
