//! G0/G1/G3 ideal-gear mechanics and real coupled aerodynamic consumers.
use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::propcoupling::CANDIDATE_B;
use fs_flyer::spine::rotors::geared::{
    self, GearedAssembly, GearedLoads, GearedRotor, GearedState,
};
use fs_flyer::spine::rotors::propulsion::{
    GearedPropulsion, PropellerHandedness, PropellerMount, PropulsionInputs,
};
use fs_flyer::spine::{Loads, SixDofState};
use fs_ga::{Mat3, Quat, So3, Vec3};

const ZERO: [f64; 3] = [0.0; 3];
fn v(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn diagonal(d: [f64; 3]) -> Mat3 {
    Mat3 {
        m: [d[0], 0.0, 0.0, 0.0, d[1], 0.0, 0.0, 0.0, d[2]],
    }
}
fn rigid() -> SixDofState {
    SixDofState {
        pos_m: [2.0, -3.0, 1.0],
        vel_mps: ZERO,
        quat: [1.0, 0.0, 0.0, 0.0],
        omega_body: ZERO,
    }
}
fn rotor(j: f64, g: f64, axis: [f64; 3]) -> GearedRotor {
    GearedRotor {
        axis_body: axis,
        axial_inertia_kgm2: j,
        ratio: g,
    }
}
fn no_loads<const N: usize>() -> GearedLoads<N> {
    GearedLoads {
        carrier: Loads {
            force_n: ZERO,
            moment_nm: ZERO,
        },
        motor_nm: 0.0,
        aerodynamic_nm: [0.0; N],
    }
}
fn rotation(state: &GearedState) -> So3 {
    let q = state.rigid.quat;
    So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })
    .unwrap()
}
fn direct_energy<const N: usize>(assembly: &GearedAssembly<N>, state: &GearedState) -> f64 {
    let omega = v(state.rigid.omega_body);
    let speed = assembly.speed(state).unwrap();
    let mut energy = 0.5 * assembly.mass_kg() * v(state.rigid.vel_mps).dot(v(state.rigid.vel_mps))
        + 0.5 * omega.dot(assembly.locked_inertia_kgm2().apply(omega));
    for rotor in assembly.rotors() {
        let relative = rotor.ratio * speed;
        energy += rotor.axial_inertia_kgm2 * relative * v(rotor.axis_body).dot(omega)
            + 0.5 * rotor.axial_inertia_kgm2 * relative * relative;
    }
    energy
}
fn direct_momentum<const N: usize>(
    assembly: &GearedAssembly<N>,
    state: &GearedState,
    origin: [f64; 3],
) -> Vec3 {
    let mut spin = assembly
        .locked_inertia_kgm2()
        .apply(v(state.rigid.omega_body));
    let speed = assembly.speed(state).unwrap();
    for rotor in assembly.rotors() {
        spin = spin + v(rotor.axis_body).scale(rotor.axial_inertia_kgm2 * rotor.ratio * speed);
    }
    rotation(state).rotate(spin).unwrap()
        + (v(state.rigid.pos_m) - v(origin)).cross(v(state.rigid.vel_mps).scale(assembly.mass_kg()))
}

#[test]
fn g1_nonunit_gear_ratio_motor_recoil_phase_and_work_are_analytic() {
    let a = GearedAssembly::new(
        3.0,
        diagonal([3.0, 4.0, 5.0]),
        [rotor(0.5, 2.0, [1.0, 0.0, 0.0])],
    )
    .unwrap();
    assert_eq!(a.equivalent_inertia_kgm2(), 2.0);
    assert_eq!(a.cross_inertia_kgm2(), [1.0, 0.0, 0.0]);
    let state = a.state_from_speed(rigid(), 0.0, 0.0).unwrap();
    let load = GearedLoads {
        motor_nm: 1.6,
        ..no_loads()
    };
    let result = geared::step_held(&a, &state, load, ZERO, 0.05).unwrap();
    assert!((result.state.momentum - 0.08).abs() < 1e-15);
    assert!((result.state.rigid.omega_body[0] + 0.016).abs() < 1e-15);
    assert!((result.speed_rad_s - 0.048).abs() < 1e-15);
    assert!((result.relative_speed_rad_s[0] - 0.096).abs() < 1e-15);
    assert!((result.absolute_speed_rad_s[0] - 0.08).abs() < 1e-15);
    assert!((result.state.phase_rad - 0.0012).abs() < 1e-15);
    assert!((result.motor_work_j - 0.00192).abs() < 1e-15);
    assert!(result.kinetic_work_residual_j.abs() < 1e-15);
    assert!(v(result.angular_momentum_residual).norm() < 1e-15);
    assert_eq!(result.external_angular_impulse_world, ZERO);
}

#[test]
fn g0_counterrotation_cancels_motor_recoil_but_unequal_aero_is_external_once() {
    let a = GearedAssembly::new(
        3.0,
        diagonal([3.0, 4.0, 5.0]),
        [
            rotor(0.5, 1.0, [1.0, 0.0, 0.0]),
            rotor(0.5, -1.0, [1.0, 0.0, 0.0]),
        ],
    )
    .unwrap();
    assert_eq!(a.cross_inertia_kgm2(), ZERO);
    let state = a.state_from_speed(rigid(), 0.4, 20.0).unwrap();
    let motor = geared::step_held(
        &a,
        &state,
        GearedLoads {
            motor_nm: 2.0,
            ..no_loads()
        },
        ZERO,
        0.05,
    )
    .unwrap();
    assert_eq!(motor.state.rigid.omega_body, ZERO);
    assert_eq!(
        motor.relative_speed_rad_s[0],
        -motor.relative_speed_rad_s[1]
    );
    let load = GearedLoads {
        aerodynamic_nm: [-0.4, 0.7],
        ..no_loads()
    };
    let result = geared::step_held(&a, &state, load, ZERO, 0.05).unwrap();
    assert!((result.state.momentum - (20.0 - 1.1 * 0.05)).abs() < 1e-14);
    assert!((result.state.rigid.omega_body[0] - 0.3 * 0.05 / 3.0).abs() < 1e-14);
    assert!((result.external_angular_impulse_world[0] - 0.3 * 0.05).abs() < 1e-14);
    let absolute_mid_left = 0.5 * (20.0 + result.absolute_speed_rad_s[0]);
    let absolute_mid_right = 0.5 * (-20.0 + result.absolute_speed_rad_s[1]);
    let external = 0.05 * (-0.4 * absolute_mid_left + 0.7 * absolute_mid_right);
    assert!((result.external_work_j - external).abs() < 1e-14);
    assert!(result.kinetic_work_residual_j.abs() < 1e-13);
}

#[test]
fn g3_arbitrary_axes_retain_products_of_inertia_and_physical_energy_momentum() {
    let a = GearedAssembly::new(
        3.5,
        Mat3 {
            m: [4.0, 0.2, -0.1, 0.2, 5.0, 0.3, -0.1, 0.3, 6.0],
        },
        [
            rotor(0.3, 2.0, [0.6, 0.8, 0.0]),
            rotor(0.7, -1.5, [0.0, 0.0, 1.0]),
        ],
    )
    .unwrap();
    assert!(a.reduced_inertia_kgm2().m[1].abs() > 0.1);
    let mut r = rigid();
    r.vel_mps = [0.7, 0.3, -0.2];
    r.omega_body = [0.4, -0.2, 0.3];
    let mut state = a.state_from_speed(r, 0.1, 12.0).unwrap();
    let origin = [-3.0, 1.0, -2.0];
    let initial = direct_momentum(&a, &state, origin);
    let mut accumulated = Vec3::default();
    for _ in 0..200 {
        let load = GearedLoads {
            carrier: Loads {
                force_n: [0.2, -0.1, 0.3],
                moment_nm: [0.01, -0.02, 0.03],
            },
            motor_nm: 0.04,
            aerodynamic_nm: [-0.05, 0.04],
        };
        let result = geared::step_held(&a, &state, load, origin, 0.003).unwrap();
        assert!((result.kinetic_energy_before_j - direct_energy(&a, &state)).abs() < 2e-12);
        assert!((result.kinetic_energy_after_j - direct_energy(&a, &result.state)).abs() < 2e-12);
        assert!(v(result.linear_momentum_residual).norm() < 2e-15);
        assert!(v(result.angular_momentum_residual).norm() < 2e-12);
        accumulated = accumulated + v(result.external_angular_impulse_world);
        state = result.state;
    }
    assert!((direct_momentum(&a, &state, origin) - initial - accumulated).norm() < 2e-10);
}

#[test]
fn g1_state_dependent_motor_and_force_have_second_order_midpoint_refinement() {
    let a = GearedAssembly::new(
        2.0,
        diagonal([4.0, 4.0, 4.0]),
        [
            rotor(0.5, 1.0, [1.0, 0.0, 0.0]),
            rotor(0.5, -1.0, [1.0, 0.0, 0.0]),
        ],
    )
    .unwrap();
    let mut r = rigid();
    r.vel_mps = [3.0, 0.0, 0.0];
    let initial = a.state_from_speed(r, 0.0, 10.0).unwrap();
    let exact_speed = 10.0 * f64::exp(-0.7 * 0.4);
    let exact_velocity = 3.0 * f64::exp(-1.2 * 0.4);
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let dt = 0.4 / f64::from(count);
        let mut state = initial;
        for k in 0..count {
            state = geared::step(&a, &state, f64::from(k) * dt, dt, ZERO, |_, s| {
                Ok(GearedLoads {
                    motor_nm: -0.7 * s.momentum,
                    carrier: Loads {
                        force_n: [-2.0 * 1.2 * s.rigid.vel_mps[0], 0.0, 0.0],
                        moment_nm: ZERO,
                    },
                    aerodynamic_nm: [0.0; 2],
                })
            })
            .unwrap()
            .state;
        }
        errors.push(
            (a.speed(&state).unwrap() - exact_speed).abs()
                + (state.rigid.vel_mps[0] - exact_velocity).abs(),
        );
    }
    eprintln!("midpoint analytic errors {errors:?}");
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])));
    }
    let mut calls = Vec::new();
    geared::step(&a, &initial, 2.0, 1e-6, ZERO, |t, s| {
        calls.push((t, s.momentum));
        Ok(no_loads())
    })
    .unwrap();
    assert_eq!(calls.len(), 2);
    assert_eq!(calls[0].0, 2.0);
    assert_eq!(calls[1].0, 2.0 + 0.5e-6);
}

#[test]
fn g0_invalid_gear_data_and_states_refuse_before_load_callbacks() {
    let definition = rotor(0.5, 1.0, [1.0, 0.0, 0.0]);
    assert!(GearedAssembly::<0>::new(1.0, diagonal([2.0; 3]), []).is_err());
    assert!(GearedAssembly::new(1.0, diagonal([100.0; 3]), [definition; 64]).is_ok());
    assert!(GearedAssembly::new(1.0, diagonal([100.0; 3]), [definition; 65]).is_err());
    for bad in [
        rotor(0.5, 0.0, [1.0, 0.0, 0.0]),
        rotor(-0.5, 1.0, [1.0, 0.0, 0.0]),
        rotor(0.5, 1.0, [2.0, 0.0, 0.0]),
    ] {
        assert!(GearedAssembly::new(1.0, diagonal([2.0; 3]), [bad]).is_err());
    }
    assert!(GearedAssembly::new(1.0, diagonal([0.4, 2.0, 3.0]), [definition]).is_err());
    let a = GearedAssembly::new(1.0, diagonal([2.0; 3]), [definition]).unwrap();
    let mut state = a.state_from_speed(rigid(), 0.0, -5.0).unwrap();
    assert_eq!(a.speed(&state).unwrap(), -5.0);
    state.rigid.quat = [2.0, 0.0, 0.0, 0.0];
    let mut calls = 0;
    assert!(
        geared::step(&a, &state, 0.0, 0.01, ZERO, |_, _| {
            calls += 1;
            Ok(no_loads())
        })
        .is_err()
    );
    assert_eq!(calls, 0);
}

fn mounts() -> [PropellerMount; 2] {
    [
        PropellerMount {
            rotor_index: 0,
            handedness: PropellerHandedness::Positive,
            thrust_direction_body: [1.0, 0.0, 0.0],
        },
        PropellerMount {
            rotor_index: 1,
            handedness: PropellerHandedness::Negative,
            thrust_direction_body: [1.0, 0.0, 0.0],
        },
    ]
}
fn aero_assembly(mass: f64) -> GearedAssembly<2> {
    // Synthetic mechanical fixture values, never historical rotor calibration.
    GearedAssembly::new(
        mass,
        diagonal([1000.0, 400.0, 1200.0]),
        [
            rotor(1.5, 1.0, [1.0, 0.0, 0.0]),
            rotor(2.5, -1.0, [1.0, 0.0, 0.0]),
        ],
    )
    .unwrap()
}
fn flight_state(a: &GearedAssembly<2>) -> GearedState {
    let mut r = rigid();
    r.pos_m = [0.0, 0.0, -10.0];
    r.vel_mps = [12.0 * f64::cos(0.08), 0.0, 12.0 * f64::sin(0.08)];
    a.state_from_speed(r, 0.0, 35.0).unwrap()
}
fn inputs() -> PropulsionInputs {
    PropulsionInputs {
        density_kg_m3: 1.225,
        canard_rad: 0.03,
        motor_nm: 80.0,
        gravity_world_mps2: ZERO,
    }
}

#[test]
fn g0_real_coupled_force_matches_existing_symmetric_aerodynamics_and_signed_ports() {
    let mut design = wright_openloop_v1();
    design.coupling = CANDIDATE_B;
    let a = aero_assembly(design.gross_mass_kg);
    let model = GearedPropulsion::new(&design, &a, mounts()).unwrap();
    let state = flight_state(&a);
    let evaluation = model
        .evaluate(&state, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    let legacy = design
        .force_buildup_with_incident(12.0, 0.08, 0.03, 35.0, 0.0, 1.225, None, &|_| {
            Ok(state.rigid.vel_mps)
        })
        .unwrap();
    let weight = design.gross_mass_kg * 9.80665;
    let expected = [
        legacy.force_n[0] + weight * f64::sin(0.08),
        legacy.force_n[1],
        legacy.force_n[2] - weight * f64::cos(0.08),
    ];
    assert!((v(expected) - v(evaluation.aerodynamic_force_body_n)).norm() < 1e-9);
    assert!((legacy.moment_y_nm - evaluation.loads.carrier.moment_nm[1]).abs() < 1e-9);
    assert!(evaluation.loads.carrier.moment_nm[0].abs() < 1e-9);
    assert!(evaluation.loads.carrier.moment_nm[2].abs() < 1e-9);
    assert_eq!(
        evaluation.loads.aerodynamic_nm[0],
        -evaluation.coupled.torque_nm[0]
    );
    assert_eq!(
        evaluation.loads.aerodynamic_nm[1],
        evaluation.coupled.torque_nm[1]
    );
    assert_eq!(evaluation.propeller_absolute_speed_rad_s, [35.0, -35.0]);
}

#[test]
fn g3_real_asymmetric_aerodynamics_drive_all_motion_and_close_external_momentum() {
    let mut design = wright_openloop_v1();
    design.coupling = CANDIDATE_B;
    let a = aero_assembly(design.gross_mass_kg);
    let model = GearedPropulsion::new(&design, &a, mounts()).unwrap();
    let mut rigid = flight_state(&a).rigid;
    rigid.omega_body = [0.2, 0.03, -0.04];
    let state = a.state_from_speed(rigid, 0.1, 35.0).unwrap();
    let original = state;
    let mut times = Vec::new();
    let result = model
        .step(
            &state,
            1.0,
            0.002,
            ZERO,
            None,
            |time, s| {
                times.push(time);
                Ok(PropulsionInputs {
                    motor_nm: 2.0 * fs_airscrew::engine_torque_at_prop_nm(a.speed(s)?),
                    gravity_world_mps2: [0.0, 0.0, 9.80665],
                    ..inputs()
                })
            },
            |_, point| Ok([0.2 * point[1], 0.0, 0.0]),
        )
        .unwrap();
    assert_eq!(state, original);
    assert_eq!(times, [1.0, 1.001]);
    assert!(
        (result.midpoint.coupled.thrust_n[0] - result.midpoint.coupled.thrust_n[1]).abs() > 0.1
    );
    assert_ne!(
        result.initial.coupled.thrust_n,
        result.midpoint.coupled.thrust_n
    );
    assert_ne!(result.motion.state.rigid.vel_mps, state.rigid.vel_mps);
    for i in 0..3 {
        assert_ne!(
            result.motion.state.rigid.omega_body[i],
            state.rigid.omega_body[i]
        );
    }
    let load = result.midpoint.loads;
    let expected_p =
        state.momentum + 0.002 * (load.motor_nm + load.aerodynamic_nm[0] - load.aerodynamic_nm[1]);
    assert!((expected_p - result.motion.state.momentum).abs() < 1e-12);
    assert!(v(result.motion.linear_momentum_residual).norm() < 1e-10);
    assert!(v(result.motion.angular_momentum_residual).norm() < 1e-9);
    let angular_change =
        direct_momentum(&a, &result.motion.state, ZERO) - direct_momentum(&a, &state, ZERO);
    assert!((angular_change - v(result.motion.external_angular_impulse_world)).norm() < 1e-9);
    assert!(
        (result.motion.motor_work_j
            - load.motor_nm * (result.motion.state.phase_rad - state.phase_rad))
            .abs()
            < 1e-12
    );
}

#[test]
fn g0_actual_absolute_spin_handedness_and_mount_domain_are_explicit() {
    let mut design = wright_openloop_v1();
    design.coupling = CANDIDATE_B;
    let a = aero_assembly(design.gross_mass_kg);
    let model = GearedPropulsion::new(&design, &a, mounts()).unwrap();
    let mut rigid = flight_state(&a).rigid;
    rigid.omega_body = [0.4, 0.0, 0.0];
    let state = a.state_from_speed(rigid, 0.0, 35.0).unwrap();
    let evaluation = model
        .evaluate(&state, inputs(), &|_| Ok(ZERO), None)
        .unwrap();
    assert_eq!(evaluation.propeller_absolute_speed_rad_s, [35.4, -34.6]);
    for speed in [0.0, -1.0] {
        let state = a
            .state_from_speed(flight_state(&a).rigid, 0.0, speed)
            .unwrap();
        assert_eq!(
            model
                .evaluate(&state, inputs(), &|_| Ok(ZERO), None)
                .unwrap_err()
                .code,
            "propeller-spin-outside-domain"
        );
    }
    let mut bad = mounts();
    bad[1].handedness = PropellerHandedness::Positive;
    let wrong = GearedPropulsion::new(&design, &a, bad).unwrap();
    assert_eq!(
        wrong
            .evaluate(&flight_state(&a), inputs(), &|_| Ok(ZERO), None)
            .unwrap_err()
            .code,
        "propeller-spin-outside-domain"
    );
    let mut bad = mounts();
    bad[1].thrust_direction_body = [-1.0, 0.0, 0.0];
    assert!(GearedPropulsion::new(&design, &a, bad).is_err());
    assert_eq!(
        model
            .evaluate(&flight_state(&a), inputs(), &|_| Ok([20.0, 0.0, 0.0]), None)
            .unwrap_err()
            .code,
        "propulsion-state-invalid"
    );
}

#[test]
fn g1_real_coupled_aerodynamic_midpoint_refines_with_step_size() {
    let mut design = wright_openloop_v1();
    design.coupling = CANDIDATE_B;
    design.coupling.tol = 1e-10;
    design.coupling.cap = 64;
    let a = aero_assembly(design.gross_mass_kg);
    let model = GearedPropulsion::new(&design, &a, mounts()).unwrap();
    let initial = flight_state(&a);
    let mut endpoints = Vec::new();
    for count in [2, 4, 8, 16] {
        let dt = 0.04 / f64::from(count);
        let mut state = initial;
        let mut warm = None;
        for k in 0..count {
            let result = model
                .step(
                    &state,
                    f64::from(k) * dt,
                    dt,
                    ZERO,
                    warm,
                    |_, _| Ok(inputs()),
                    |_, point| Ok([0.1 * point[1], 0.0, 0.0]),
                )
                .unwrap();
            warm = Some(result.midpoint.coupled.w_slip);
            state = result.motion.state;
        }
        endpoints.push(state);
    }
    let difference = |a: &GearedState, b: &GearedState| {
        (v(a.rigid.vel_mps) - v(b.rigid.vel_mps)).norm()
            + (v(a.rigid.omega_body) - v(b.rigid.omega_body)).norm()
            + (a.momentum - b.momentum).abs()
    };
    let errors: Vec<_> = endpoints
        .windows(2)
        .map(|pair| difference(&pair[0], &pair[1]))
        .collect();
    eprintln!("real coupled successive differences {errors:?}");
    for pair in errors.windows(2) {
        assert!((3.5..4.5).contains(&(pair[0] / pair[1])));
    }
}
