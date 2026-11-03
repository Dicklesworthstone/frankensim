//! G1/G3 stage coupling, independent ODE refinement and transactional publication.
use super::*;

fn physical_energy<const N: usize>(
    assembly: &RotorAssembly<N>,
    state: &RotorState<N>,
    properties: ShaftProperties<N>,
) -> f64 {
    let omega = v(state.rigid.omega_body);
    let velocity = v(state.rigid.vel_mps);
    let inertia = assembly.reduced_body.inertia_tensor().matrix();
    let mut energy = 0.5 * assembly.reduced_body.mass_kg * velocity.dot(velocity)
        + 0.5 * omega.dot(inertia.apply(omega));
    for i in 0..N {
        let strain = state.relative_phase_rad[i] - properties.rest_phase_rad[i];
        energy += 0.5 * state.axial_momentum[i].powi(2) / assembly.rotors[i].axial_inertia_kgm2
            + 0.5 * properties.stiffness_nm_per_rad[i] * strain * strain;
    }
    energy
}

fn feedback(time: f64, state: &RotorState<2>) -> (Loads, RotorTorques<2>) {
    let relative = [
        state.axial_momentum[0] / 0.2 - state.rigid.omega_body[0],
        state.axial_momentum[1] / 0.2 - state.rigid.omega_body[1],
    ];
    (
        Loads {
            force_n: [
                -0.125 * state.rigid.pos_m[0] + 0.0625 * time,
                -0.0625 * state.rigid.vel_mps[1],
                -0.25 * state.rigid.pos_m[2],
            ],
            moment_nm: [
                0.03125 * state.rigid.omega_body[1],
                -0.0625 * state.rigid.omega_body[2],
                0.125 * state.rigid.omega_body[0],
            ],
        },
        RotorTorques {
            motor_nm: [
                0.5 - 0.125 * relative[0] + 0.0625 * time,
                -0.25 - 0.0625 * relative[1],
            ],
            aerodynamic_nm: [
                -0.015625 * state.axial_momentum[0] / 0.2,
                -0.03125 * state.axial_momentum[1] / 0.2,
            ],
        },
    )
}

#[test]
fn g3_native_stage_oracle_retains_full_tensor_momentum_across_impulses() {
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
    let half = 0.5 * dt;
    let origin = [0.5, -0.25, 0.75];
    let mut motion = assembly.shaft_motion(&initial, properties, dt).unwrap();
    let bank = RotorBank::from_reduced(
        assembly.reduced_body.inertia_tensor(),
        assembly
            .rotors
            .map(|r| NativeRotor::from_unit_axis(v(r.axis_body), r.axial_inertia_kgm2).unwrap()),
    )
    .unwrap();
    let prepare = |duration| {
        bank.torsional_stepper(
            properties.stiffness_nm_per_rad,
            properties.rest_phase_rad,
            properties.damping_nms_per_rad,
            duration,
        )
        .unwrap()
    };
    let predictor_owner = prepare(half);
    let accepted_owner = prepare(dt);
    let native_loads = |loads: Loads, torques: RotorTorques<2>| {
        (
            core::array::from_fn(|i| torques.motor_nm[i] + torques.aerodynamic_nm[i]),
            HeldTorque::Body(v([
                loads.moment_nm[0] + torques.aerodynamic_nm[0],
                loads.moment_nm[1] + torques.aerodynamic_nm[1],
                loads.moment_nm[2],
            ])),
        )
    };
    for tick in 0..8 {
        if tick == 0 || tick == 4 {
            motion
                .apply_impulse(
                    RotorImpulses {
                        carrier: CarrierImpulse {
                            linear_world_ns: [0.3, -0.15, 0.075],
                            angular_body_nms: [0.025, -0.0125, 0.00625],
                        },
                        motor_nms: [0.02, -0.01],
                        aerodynamic_nms: [-0.005, 0.0025],
                    },
                    origin,
                )
                .unwrap();
        }
        let before = *motion.state();
        let canonical = *motion.native_state();
        let time = 2.0 + f64::from(tick) * dt;
        let (force0, torque0) = feedback(time, &before);
        let (drive0, external0) = native_loads(force0, torque0);
        let prediction = predictor_owner.step(&canonical, drive0, external0).unwrap();
        let predicted = prediction.motion.state_after;
        let mut expected_midpoint = before;
        for j in 0..3 {
            expected_midpoint.rigid.vel_mps[j] =
                before.rigid.vel_mps[j] + (force0.force_n[j] * half) / 3.0;
            expected_midpoint.rigid.pos_m[j] = before.rigid.pos_m[j]
                + before.rigid.vel_mps[j] * (0.5 * half)
                + expected_midpoint.rigid.vel_mps[j] * (0.5 * half);
        }
        set_rotation(&mut expected_midpoint.rigid, predicted.rotation);
        expected_midpoint.rigid.omega_body = array(bank.angular_velocity(&predicted).unwrap());
        expected_midpoint.axial_momentum = predicted.axial_momentum_nms;
        expected_midpoint.relative_phase_rad = predicted.relative_angle_rad;
        let (force_mid, torque_mid) = feedback(time + half, &expected_midpoint);
        let (drive_mid, external_mid) = native_loads(force_mid, torque_mid);
        let expected = accepted_owner
            .step(&canonical, drive_mid, external_mid)
            .unwrap();
        let mut calls = 0;
        let accepted = motion
            .advance_with(time, origin, |at, state| {
                assert_eq!(at, time + f64::from(calls) * half);
                assert_eq!(
                    *state,
                    if calls == 0 {
                        before
                    } else {
                        expected_midpoint
                    }
                );
                calls += 1;
                Ok(feedback(at, state))
            })
            .unwrap();
        assert_eq!(calls, 2);
        assert_eq!(accepted.native, expected);
        assert_eq!(*motion.native_state(), expected.motion.state_after);
        let (p0, h0) = total(&assembly, &before, origin);
        let (p1, h1) = total(&assembly, &accepted.state, origin);
        assert!((p1 - p0 - v(accepted.external_linear_impulse_world)).norm() < 2e-12);
        assert!((h1 - h0 - v(accepted.external_angular_impulse_world)).norm() < 2e-11);
        let storage = physical_energy(&assembly, &accepted.state, properties)
            - physical_energy(&assembly, &before, properties);
        let heat: f64 = accepted.native.dissipated_energy_j.into_iter().sum();
        assert!(heat > 0.0);
        let independent_balance = storage + heat
            - accepted.motor_work_j
            - accepted.external_work_j
            - accepted.native.motion.drift_energy_defect_j;
        assert!(independent_balance.abs() < 2e-12);
        // The rejected predictor's heat must not be included in accepted storage.
        assert!(
            (independent_balance + prediction.dissipated_energy_j.iter().sum::<f64>()).abs() > 1e-6
        );
    }
}

#[test]
fn g3_constant_law_is_exactly_the_original_held_step() {
    let assembly = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let initial = assembly
        .state_from_relative_spin(rigid(), [0.4], [3.0])
        .unwrap();
    let properties = ShaftProperties {
        stiffness_nm_per_rad: [1.7],
        rest_phase_rad: [0.1],
        damping_nms_per_rad: [0.08],
    };
    let mut held = assembly.shaft_motion(&initial, properties, 0.05).unwrap();
    let mut coupled = assembly.shaft_motion(&initial, properties, 0.05).unwrap();
    let force = Loads {
        force_n: [0.3, -0.2, 0.1],
        moment_nm: [0.01, 0.02, -0.03],
    };
    let torque = RotorTorques {
        motor_nm: [0.2],
        aerodynamic_nm: [-0.1],
    };
    for step in 0..10 {
        let expected = held.advance(force, torque, ZERO).unwrap();
        let actual = coupled
            .advance_with(f64::from(step) * 0.05, ZERO, |_, _| Ok((force, torque)))
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(coupled.native_state(), held.native_state());
    }
}

// Independent coaxial equations in [spring strain, shaft momentum, total
// angular momentum, COM position, COM velocity]. This oracle has no gyrostat
// time-step or Flyer state-conversion calls.
fn rhs(time: f64, y: [f64; 5]) -> [f64; 5] {
    let omega = y[2] - y[1]; // reduced I = 1
    let relative = y[1] / 0.2 - omega;
    let motor = 0.3 - 0.11 * relative + 0.04 * time;
    let aerodynamic = -0.07 * y[1] / 0.2;
    [
        relative,
        -1.7 * y[0] - 0.08 * relative + motor + aerodynamic,
        -0.09 * omega + aerodynamic,
        y[4],
        (-1.3 * y[3] - 0.2 * y[4] + 0.1 * time) / 3.0,
    ]
}

#[test]
fn g1_state_dependent_shaft_and_carrier_loads_converge_to_independent_ode() {
    let assembly = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let mut carrier = rigid();
    carrier.omega_body[0] = 0.2;
    carrier.vel_mps[0] = 0.4;
    let initial = assembly
        .state_from_relative_spin(carrier, [0.35], [1.2])
        .unwrap();
    let properties = ShaftProperties {
        stiffness_nm_per_rad: [1.7],
        rest_phase_rad: [0.1],
        damping_nms_per_rad: [0.08],
    };
    let initial_y = [
        0.25,
        initial.axial_momentum[0],
        0.2 + initial.axial_momentum[0],
        2.0,
        0.4,
    ];
    let mut reference = initial_y;
    let h = 0.5 / 8192.0;
    for i in 0..8192 {
        let time = f64::from(i) * h;
        let stage = |k: [f64; 5], scale: f64| core::array::from_fn(|j| reference[j] + scale * k[j]);
        let k1 = rhs(time, reference);
        let k2 = rhs(time + 0.5 * h, stage(k1, 0.5 * h));
        let k3 = rhs(time + 0.5 * h, stage(k2, 0.5 * h));
        let k4 = rhs(time + h, stage(k3, h));
        reference = core::array::from_fn(|j| {
            reference[j] + h * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0
        });
    }
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let dt = 0.5 / f64::from(steps);
        let mut motion = assembly.shaft_motion(&initial, properties, dt).unwrap();
        for i in 0..steps {
            motion
                .advance_with(f64::from(i) * dt, ZERO, |time, state| {
                    let omega = state.rigid.omega_body[0];
                    let absolute = state.axial_momentum[0] / 0.2;
                    let mut loads = no_loads();
                    loads.force_n[0] =
                        -1.3 * state.rigid.pos_m[0] - 0.2 * state.rigid.vel_mps[0] + 0.1 * time;
                    loads.moment_nm[0] = -0.09 * omega;
                    Ok((
                        loads,
                        RotorTorques {
                            motor_nm: [0.3 - 0.11 * (absolute - omega) + 0.04 * time],
                            aerodynamic_nm: [-0.07 * absolute],
                        },
                    ))
                })
                .unwrap();
        }
        let state = motion.state();
        let actual = [
            state.relative_phase_rad[0] - 0.1,
            state.axial_momentum[0],
            state.rigid.omega_body[0] + state.axial_momentum[0],
            state.rigid.pos_m[0],
            state.rigid.vel_mps[0],
        ];
        errors.push(
            actual
                .into_iter()
                .zip(reference)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt(),
        );
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.6..4.4).contains(&ratio),
            "second-order refinement: {errors:?}"
        );
    }
    assert!(errors[2] < 5e-5, "independent ODE discrepancy: {errors:?}");
}

#[test]
fn g0_stage_refusals_and_unresolved_clock_preserve_impulsed_state_and_retry() {
    let assembly = RotorAssembly {
        reduced_body: body(),
        rotors: [rotor([1.0, 0.0, 0.0])],
    };
    let initial = assembly
        .state_from_relative_spin(rigid(), [0.4], [3.0])
        .unwrap();
    let properties = ShaftProperties {
        stiffness_nm_per_rad: [1.7],
        rest_phase_rad: [0.1],
        damping_nms_per_rad: [0.08],
    };
    let mut motion = assembly.shaft_motion(&initial, properties, 0.05).unwrap();
    let mut fresh = assembly.shaft_motion(&initial, properties, 0.05).unwrap();
    let impulse = RotorImpulses {
        carrier: CarrierImpulse::default(),
        motor_nms: [0.125],
        aerodynamic_nms: [-0.0625],
    };
    motion.apply_impulse(impulse, ZERO).unwrap();
    fresh.apply_impulse(impulse, ZERO).unwrap();
    let state = *motion.state();
    let native = *motion.native_state();
    for time in [f64::NAN, f64::INFINITY, 1e30] {
        let mut calls = 0;
        assert!(
            motion
                .advance_with(time, ZERO, |_, _| {
                    calls += 1;
                    Ok((no_loads(), zero_torques()))
                })
                .is_err()
        );
        assert_eq!(calls, 0);
        assert_eq!(*motion.state(), state);
        assert_eq!(*motion.native_state(), native);
    }
    for failing_stage in [0, 1] {
        let mut calls = 0;
        let error = motion
            .advance_with(0.0, ZERO, |_, _| {
                let stage = calls;
                calls += 1;
                if stage == failing_stage {
                    return Err(fs_flyer::Refusal {
                        code: "upstream-load-cancelled",
                        message: "cancelled load evaluation".into(),
                        ranked_repairs: vec![],
                    });
                }
                Ok((no_loads(), zero_torques()))
            })
            .unwrap_err();
        assert_eq!(error.code, "upstream-load-cancelled");
        assert_eq!(calls, failing_stage + 1);
        assert_eq!(*motion.state(), state);
        assert_eq!(*motion.native_state(), native);
    }
    for failed_force in [f64::NAN, 1e308] {
        let mut calls = 0;
        assert!(
            motion
                .advance_with(0.0, ZERO, |_, _| {
                    calls += 1;
                    let mut force = no_loads();
                    if calls == 2 {
                        force.force_n[0] = failed_force;
                    }
                    Ok((force, zero_torques()))
                })
                .is_err()
        );
        assert_eq!(calls, 2);
        assert_eq!(*motion.state(), state);
        assert_eq!(*motion.native_state(), native);
    }
    assert_eq!(
        motion
            .advance_with(0.0, ZERO, |_, _| Ok((no_loads(), zero_torques())))
            .unwrap(),
        fresh.advance(no_loads(), zero_torques(), ZERO).unwrap(),
    );
}
