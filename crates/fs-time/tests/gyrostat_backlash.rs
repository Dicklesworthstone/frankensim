use fs_ga::{Mat3, So3, Vec3};
use fs_la::factor::lu;
use fs_time::lie::gyrostat::rotors::{
    AxialRotor, BacklashContact, BacklashCoupling, BacklashFlank, BacklashStepError,
    MAX_BACKLASH_ROOT_ITERATIONS, RotorBank, RotorState, TorsionalCoupling,
};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn zero() -> HeldTorque {
    HeldTorque::Body(Vec3::default())
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn flank(k: f64, c: f64) -> BacklashFlank {
    BacklashFlank {
        stiffness_nm_per_rad: k,
        compression_damping_nms_per_rad: c,
    }
}
fn law(kp: f64, kn: f64, cp: f64, cn: f64) -> BacklashCoupling {
    BacklashCoupling {
        first_rotor: 0,
        second_rotor: 1,
        ratio: 1.0,
        center_angle_rad: 0.0,
        half_gap_rad: 0.1,
        positive_flank: flank(kp, cp),
        negative_flank: flank(kn, cn),
    }
}
fn bank() -> RotorBank<2> {
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(3.0, 3.0, 2.0), So3::identity()).unwrap(),
        [0.4, 0.6].map(|j| AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), j).unwrap()),
    )
    .unwrap()
}
fn state(bank: &RotorBank<2>, x: f64, v: f64) -> RotorState<2> {
    bank.state_from_velocities(
        So3::identity(),
        Vec3::new(0.0, 0.0, 0.3),
        [0.6 * v, -0.4 * v],
        [0.6 * x, -0.4 * x],
    )
    .unwrap()
}
fn coordinate(law: BacklashCoupling, state: &RotorState<2>) -> f64 {
    state.relative_angle_rad[0] - law.ratio * state.relative_angle_rad[1] - law.center_angle_rad
}
fn potential(law: BacklashCoupling, state: &RotorState<2>) -> f64 {
    let x = coordinate(law, state);
    0.5 * law.positive_flank.stiffness_nm_per_rad * (x - law.half_gap_rad).max(0.0).powi(2)
        + 0.5 * law.negative_flank.stiffness_nm_per_rad * (-x - law.half_gap_rad).max(0.0).powi(2)
}

#[test]
fn g1_free_gap_and_zero_gap_linear_limit_keep_original_shaft_mechanics() {
    let bank = bank();
    let coupling = law(3.0, 5.0, 0.2, 0.4);
    let initial = state(&bank, 0.0, 0.5);
    let prepared = bank.backlash_stepper(coupling, 0.04).unwrap();
    let step = prepared.step(&initial, [0.0; 2], zero()).unwrap();
    assert_eq!(step.contact_after, BacklashContact::Gap);
    assert_eq!(step.root_iterations, [0; 2]);
    assert_eq!(step.compressive_angular_impulse_nms, [0.0; 2]);
    assert_eq!(step.dissipated_energy_j, [0.0; 2]);
    assert_eq!(step.elastic_energy_change_j, [0.0; 2]);
    assert_eq!(step.motion.axial_momentum_change_nms, [0.0; 2]);
    near(coordinate(coupling, &step.motion.state_after), 0.02, 1e-16);
    near(
        bank.angular_velocity(&step.motion.state_after).unwrap().z,
        0.3,
        1e-15,
    );

    // Removing clearance with equal elastic flanks gives the existing native
    // reciprocal spring, with the same negative ratio and source torques.
    let mut coupling = law(2.0, 2.0, 0.0, 0.0);
    coupling.half_gap_rad = 0.0;
    coupling.ratio = -1.5;
    coupling.center_angle_rad = 0.05;
    let initial = state(&bank, 0.2, 0.7);
    let backlash = bank.backlash_stepper(coupling, 0.01).unwrap();
    let spring = bank
        .torsional_network_stepper(
            [TorsionalCoupling {
                first_rotor: 0,
                second_rotor: 1,
                ratio: coupling.ratio,
                stiffness_nm_per_rad: 2.0,
                damping_nms_per_rad: 0.0,
                rest_angle_rad: 0.05,
            }],
            0.01,
        )
        .unwrap();
    let drive = [0.3, -0.2];
    let external = HeldTorque::Body(Vec3::new(0.0, 0.0, 0.1));
    let a = backlash.step(&initial, drive, external).unwrap();
    let b = spring.step(&initial, drive, external).unwrap();
    for i in 0..2 {
        near(
            a.motion.state_after.axial_momentum_nms[i],
            b.motion.state_after.axial_momentum_nms[i],
            2e-14,
        );
        near(
            a.motion.state_after.relative_angle_rad[i],
            b.motion.state_after.relative_angle_rad[i],
            2e-14,
        );
    }
    near(a.work_energy_residual_j, 0.0, 2e-14);
    near(
        a.elastic_energy_change_j.iter().sum(),
        b.elastic_energy_change_j[0],
        2e-14,
    );
}

#[test]
fn g1_gap_crossings_and_both_flanks_refine_to_freeflight_and_exact_contact_arcs() {
    let bank = bank();
    let mobility = 1.0 / 0.4 + 1.0 / 0.6;
    // Unit incident slip: entry at t=.1, positive-flank half-period 1 s,
    // exit at1.1, gap traversal .2 s, negative-flank entry at1.3.
    let coupling = law(
        core::f64::consts::PI.powi(2) / mobility,
        4.0 * core::f64::consts::PI.powi(2) / mobility,
        0.0,
        0.0,
    );
    let initial = state(&bank, 0.0, 1.0);
    let time = 1.5;
    let omega = 2.0 * core::f64::consts::PI;
    let expected_x = -0.1 - (omega * (time - 1.3)).sin() / omega;
    let expected_v = -(omega * (time - 1.3)).cos();
    let initial_energy = bank.kinetic_energy(&initial).unwrap();
    let mut previous = None;
    for count in [64, 128, 256] {
        let prepared = bank
            .backlash_stepper(coupling, time / f64::from(count))
            .unwrap();
        let mut current = initial;
        let mut compressive = [0.0; 2];
        let mut crossed = false;
        for index in 0..count {
            let step = prepared
                .step(&current, [0.0; 2], zero())
                .unwrap_or_else(|error| {
                    panic!(
                        "elastic crossing count={count}, step={index}, state={current:?}: {error}"
                    )
                });
            assert_eq!(step.dissipated_energy_j, [0.0; 2]);
            assert!(
                step.compressive_angular_impulse_nms
                    .iter()
                    .all(|v| *v >= 0.0)
            );
            crossed |= step.root_iterations.iter().any(|n| *n > 0);
            for side in 0..2 {
                compressive[side] += step.compressive_angular_impulse_nms[side];
            }
            near(step.work_energy_residual_j, 0.0, 2e-13);
            assert!(step.motion.angular_momentum_residual_world.norm() < 3e-14);
            current = step.motion.state_after;
            near(
                bank.kinetic_energy(&current).unwrap() + potential(coupling, &current),
                initial_energy,
                3e-12,
            );
            near(bank.angular_velocity(&current).unwrap().z, 0.3, 2e-14);
        }
        assert!(crossed && compressive[0] > 0.0 && compressive[1] > 0.0);
        assert_eq!(
            prepared.contact(&current).unwrap(),
            BacklashContact::NegativeFlank
        );
        let rates = bank.relative_rates(&current).unwrap();
        let error = (coordinate(coupling, &current) - expected_x).abs()
            + (rates[0] - rates[1] - expected_v).abs();
        if let Some(coarse) = previous {
            assert!(
                coarse / error > 3.5 && coarse / error < 4.5,
                "crossing refinement {coarse:e}/{error:e}"
            );
        }
        previous = Some(error);
    }
}

#[test]
fn g1_compression_damping_dissipates_without_attractive_unloading_or_gap_heat() {
    let bank = bank();
    let coupling = law(2.0, 3.0, 10.0, 12.0);
    let prepared = bank.backlash_stepper(coupling, 0.01).unwrap();
    // A conventional k*z+c*v would be negative here. This authored law returns
    // elastic energy without tensile tooth force or unloading dashpot heat.
    let unloading = state(&bank, 0.25, -1.0);
    let out = prepared.step(&unloading, [0.0; 2], zero()).unwrap();
    assert_eq!(out.root_iterations, [0; 2]);
    assert_eq!(out.dissipated_energy_j, [0.0; 2]);
    assert!(out.compressive_angular_impulse_nms[0] > 0.0);
    assert_eq!(out.compressive_angular_impulse_nms[1], 0.0);
    assert!(out.motion.axial_momentum_change_nms[0] < 0.0);
    assert!(out.elastic_energy_change_j[0] < 0.0);
    assert!(out.motion.kinetic_energy_change_j > 0.0);
    near(out.work_energy_residual_j, 0.0, 2e-14);
    let loading = state(&bank, 0.25, 1.0);
    let incoming = prepared.step(&loading, [0.0; 2], zero()).unwrap();
    assert!(incoming.dissipated_energy_j[0] > 0.0);
    assert_eq!(incoming.dissipated_energy_j[1], 0.0);
    assert!(incoming.elastic_energy_change_j[0] > 0.0);
    near(incoming.work_energy_residual_j, 0.0, 2e-13);

    let coupling = law(3.0, 5.0, 0.2, 0.4);
    let prepared = bank.backlash_stepper(coupling, 3.0 / 257.0).unwrap();
    let initial = state(&bank, 0.0, 1.0);
    let initial_energy = bank.kinetic_energy(&initial).unwrap();
    let mut current = initial;
    let mut heat = [0.0; 2];
    for index in 0..257 {
        let step = prepared
            .step(&current, [0.0; 2], zero())
            .unwrap_or_else(|error| {
                panic!("damped crossing step={index}, state={current:?}: {error}")
            });
        for side in 0..2 {
            assert!(step.dissipated_energy_j[side] >= 0.0);
            heat[side] += step.dissipated_energy_j[side];
        }
        current = step.motion.state_after;
        near(
            bank.kinetic_energy(&current).unwrap()
                + potential(coupling, &current)
                + heat.iter().sum::<f64>(),
            initial_energy,
            3e-12,
        );
    }
    assert!(heat[0] > 0.0 && heat[1] > 0.0);
}

#[test]
fn g3_oblique_loaded_mesh_refines_against_original_euler_and_shaft_equations() {
    let matrix = Mat3 {
        m: [2.4, 0.3, -0.2, 0.3, 3.1, 0.25, -0.2, 0.25, 4.0],
    };
    let rotors = [
        AxialRotor::new(Vec3::new(1.0, 0.2, 0.0), 0.4).unwrap(),
        AxialRotor::new(Vec3::new(0.3, 1.0, 0.4), 0.7).unwrap(),
    ];
    let bank = RotorBank::from_reduced(InertiaTensor3::new(matrix).unwrap(), rotors).unwrap();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(0.4, -0.3, 0.7),
            [1.2, 0.2],
            [0.3, 0.0],
        )
        .unwrap();
    let mut coupling = law(1.5, 2.5, 0.2, 0.35);
    coupling.ratio = -0.7;
    coupling.center_angle_rad = 0.02;
    let drive = [0.4, -0.1];
    let external = Vec3::new(0.07, -0.04, 0.03);
    let inverse = lu(&matrix.m, 3).unwrap();
    let pack = |state: &RotorState<2>| {
        [
            state.total_momentum_body.x,
            state.total_momentum_body.y,
            state.total_momentum_body.z,
            state.axial_momentum_nms[0],
            state.axial_momentum_nms[1],
            state.relative_angle_rad[0],
            state.relative_angle_rad[1],
        ]
    };
    // Continuous original-coordinate equations, independent of the stage's
    // mobility, discrete gradient, root solve and native gyrostat split.
    let rhs = |y: [f64; 7]| {
        let total = Vec3::new(y[0], y[1], y[2]);
        let carrier = total - rotors[0].axis_body().scale(y[3]) - rotors[1].axis_body().scale(y[4]);
        let mut w = [carrier.x, carrier.y, carrier.z];
        inverse.solve(&mut w);
        let omega = Vec3::new(w[0], w[1], w[2]);
        let rates = [
            y[3] / 0.4 - rotors[0].axis_body().dot(omega),
            y[4] / 0.7 - rotors[1].axis_body().dot(omega),
        ];
        let x = y[5] - coupling.ratio * y[6] - coupling.center_angle_rad;
        let v = rates[0] - coupling.ratio * rates[1];
        assert!(
            x > coupling.half_gap_rad && v > 0.0,
            "fixture retains a smooth compressing positive flank"
        );
        let effort = 1.5 * (x - 0.1) + 0.2 * v;
        let dl = total.cross(omega) + external;
        [
            dl.x,
            dl.y,
            dl.z,
            drive[0] - effort,
            drive[1] + coupling.ratio * effort,
            rates[0],
            rates[1],
        ]
    };
    let time = 0.2;
    let mut reference = pack(&initial);
    let dt = time / 4096.0;
    for _ in 0..4096 {
        let a = rhs(reference);
        let b = rhs(core::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]));
        let c = rhs(core::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]));
        let d = rhs(core::array::from_fn(|i| reference[i] + dt * c[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    let mut previous = None;
    for count in [8, 16, 32] {
        let prepared = bank
            .backlash_stepper(coupling, time / f64::from(count))
            .unwrap();
        let mut current = initial;
        let mut source = 0.0;
        let mut heat = 0.0;
        for _ in 0..count {
            let step = prepared
                .step(&current, drive, HeldTorque::Body(external))
                .unwrap();
            assert!(step.dissipated_energy_j[0] > 0.0);
            assert_eq!(step.dissipated_energy_j[1], 0.0);
            assert_eq!(step.root_iterations, [0; 2]);
            near(step.work_energy_residual_j, 0.0, 3e-13);
            assert!(step.motion.angular_momentum_residual_world.norm() < 3e-13);
            source += step.drive_work_j.iter().sum::<f64>()
                + step.motion.external_work_j
                + step.motion.drift_energy_defect_j;
            heat += step.dissipated_energy_j.iter().sum::<f64>();
            current = step.motion.state_after;
        }
        near(
            bank.kinetic_energy(&current).unwrap() + potential(coupling, &current) + heat
                - bank.kinetic_energy(&initial).unwrap()
                - potential(coupling, &initial),
            source,
            2e-12,
        );
        let error = pack(&current)
            .iter()
            .zip(reference)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>();
        if let Some(coarse) = previous {
            assert!(
                coarse / error > 3.8 && coarse / error < 4.2,
                "oblique refinement {coarse:e}/{error:e}"
            );
        }
        previous = Some(error);
    }
}

#[test]
fn g0_bounded_crossing_solve_cancels_atomically_and_zero_law_delegates() {
    let bank = bank();
    let initial = state(&bank, 0.095, 1.0);
    let prepared = bank
        .backlash_stepper(law(3.0, 5.0, 0.2, 0.4), 0.04)
        .unwrap();
    let mut calls = 0;
    let result = prepared.step_controlled(&initial, [0.0; 2], zero(), || {
        calls += 1;
        calls == 3
    });
    assert_eq!(result, Err(BacklashStepError::Cancelled));
    assert_eq!(calls, 3);
    let retry = prepared.step(&initial, [0.0; 2], zero()).unwrap();
    assert!(retry.root_iterations[0] > 0);
    assert!(
        retry
            .root_iterations
            .iter()
            .all(|n| *n <= MAX_BACKLASH_ROOT_ITERATIONS)
    );
    assert_eq!(retry, prepared.step(&initial, [0.0; 2], zero()).unwrap());
    let zero_law = bank
        .backlash_stepper(law(0.0, 0.0, 0.0, 0.0), 0.04)
        .unwrap();
    let drive = [0.2, -0.3];
    let external = HeldTorque::Space(Vec3::new(0.1, -0.03, 0.2));
    assert_eq!(
        zero_law.step(&initial, drive, external).unwrap().motion,
        bank.step(&initial, drive, external, 0.04).unwrap()
    );
    for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, 100.0] {
        assert!(bank.backlash_stepper(law(3.0, 5.0, 0.2, 0.4), dt).is_err());
    }
    let original = law(3.0, 5.0, 0.2, 0.4);
    for invalid in [
        BacklashCoupling {
            half_gap_rad: -0.1,
            ..original
        },
        BacklashCoupling {
            ratio: f64::NAN,
            ..original
        },
        BacklashCoupling {
            center_angle_rad: f64::INFINITY,
            ..original
        },
        BacklashCoupling {
            second_rotor: 0,
            ..original
        },
        BacklashCoupling {
            second_rotor: 2,
            ..original
        },
        BacklashCoupling {
            positive_flank: flank(-1.0, 0.0),
            ..original
        },
        BacklashCoupling {
            negative_flank: flank(1.0, -0.2),
            ..original
        },
    ] {
        assert!(bank.backlash_stepper(invalid, 0.01).is_err());
    }
    let before = initial;
    assert!(prepared.step(&initial, [f64::NAN, 0.0], zero()).is_err());
    assert!(
        prepared
            .step(
                &initial,
                [0.0; 2],
                HeldTorque::Body(Vec3::new(f64::INFINITY, 0.0, 0.0))
            )
            .is_err()
    );
    assert_eq!(initial, before);
    // A nominal grid endpoint can sit only a few ulps beyond first contact.
    // Its tiny physical reaction then cannot change the original stored p.
    // Retain that native refusal; the damped 3 s trajectory above uses a grid
    // with representable crossings, rather than silently discarding the force.
    let boundary_stepper = bank.backlash_stepper(original, 0.01).unwrap();
    let mut boundary = state(&bank, 0.0, 1.0);
    for _ in 0..9 {
        boundary = boundary_stepper
            .step(&boundary, [0.0; 2], zero())
            .unwrap()
            .motion
            .state_after;
    }
    let boundary_before = boundary;
    assert!(matches!(
        boundary_stepper.step(&boundary, [0.0; 2], zero()),
        Err(BacklashStepError::Native(
            fs_time::lie::LieStepError::InvalidParameter {
                context: "rotor axial impulse"
            }
        ))
    ));
    assert_eq!(boundary, boundary_before);
    let empty = RotorBank::from_reduced(bank.reduced_inertia(), []).unwrap();
    assert!(empty.backlash_stepper(original, 0.01).is_err());
    let full = RotorBank::from_reduced(bank.reduced_inertia(), [bank.rotors()[0]; 64]).unwrap();
    let coupling = BacklashCoupling {
        second_rotor: 63,
        ..original
    };
    let full_stepper = full.backlash_stepper(coupling, 0.01).unwrap();
    let rest = full
        .state_from_velocities(So3::identity(), Vec3::default(), [0.0; 64], [0.0; 64])
        .unwrap();
    assert_eq!(
        full_stepper
            .step(&rest, [0.0; 64], zero())
            .unwrap()
            .motion
            .state_after,
        rest
    );
}
