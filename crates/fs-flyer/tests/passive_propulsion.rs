//! G1/G3 real wing/BEMT consumers of independent native spring/drag shafts.
use fs_airscrew::swirl::Controls;
use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::propcoupling::CANDIDATE_B;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor};
use fs_flyer::spine::rotors::passive::ShaftProperties;
use fs_flyer::spine::rotors::propulsion::passive::{
    ShaftInputs, ShaftPropulsion, SwirlShaftPropulsion,
};
use fs_flyer::spine::rotors::propulsion::swirl::SwirlPropulsion;
use fs_flyer::spine::rotors::propulsion::{
    GearedPropulsion, PropellerHandedness, PropellerMount, PropulsionInputs,
};
use fs_flyer::spine::rotors::{
    AxialRotor, CarrierImpulse, RotorAssembly, RotorImpulses, RotorState, RotorTorques,
};
use fs_flyer::spine::{Loads, RigidBody, SixDofState};
use fs_ga::{Mat3, Quat, So3, Vec3};

const ZERO: [f64; 3] = [0.0; 3];
const AXIAL: [[f64; 3]; 2] = [[1.0, 0.0, 0.0]; 2];

#[path = "passive_propulsion/dynamic.rs"]
mod dynamic;

fn v(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn rotation(state: &RotorState<2>) -> So3 {
    let q = state.rigid.quat;
    So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })
    .unwrap()
}
fn assembly(mass: f64, axes: [[f64; 3]; 2]) -> RotorAssembly<2> {
    // Explicit synthetic inertia; no historical propeller calibration claim.
    RotorAssembly {
        reduced_body: RigidBody::from_inertia_tensor(
            mass,
            Mat3 {
                m: [1000.0, 10.0, 20.0, 10.0, 400.0, 5.0, 20.0, 5.0, 1200.0],
            },
        )
        .unwrap(),
        rotors: [
            AxialRotor {
                axis_body: axes[0],
                axial_inertia_kgm2: 2.0,
            },
            AxialRotor {
                axis_body: axes[1],
                axial_inertia_kgm2: 4.0,
            },
        ],
    }
}
fn mounts(axes: [[f64; 3]; 2]) -> [PropellerMount; 2] {
    [
        PropellerMount {
            rotor_index: 0,
            handedness: PropellerHandedness::Positive,
            thrust_direction_body: axes[0],
        },
        PropellerMount {
            rotor_index: 1,
            handedness: PropellerHandedness::Negative,
            thrust_direction_body: axes[1],
        },
    ]
}
fn initial(a: &RotorAssembly<2>, relative: [f64; 2]) -> RotorState<2> {
    a.state_from_relative_spin(
        SixDofState {
            pos_m: [1.0, -2.0, -10.0],
            vel_mps: [12.0 * 0.08_f64.cos(), 0.0, 12.0 * 0.08_f64.sin()],
            quat: [1.0, 0.0, 0.0, 0.0],
            omega_body: [0.2, 0.03, -0.04],
        },
        [0.2, -0.1],
        relative,
    )
    .unwrap()
}
fn inputs(a: &RotorAssembly<2>, state: &RotorState<2>) -> ShaftInputs<2> {
    let rate = core::array::from_fn::<_, 2, _>(|i| {
        state.axial_momentum[i] / a.rotors[i].axial_inertia_kgm2
            - v(a.rotors[i].axis_body).dot(v(state.rigid.omega_body))
    });
    ShaftInputs {
        density_kg_m3: 1.225,
        canard_rad: 0.03,
        motor_nm: [
            fs_airscrew::engine_torque_at_prop_nm(rate[0]),
            -fs_airscrew::engine_torque_at_prop_nm(-rate[1]),
        ],
        gravity_world_mps2: [0.0, 0.0, 9.80665],
    }
}
fn properties() -> ShaftProperties<2> {
    ShaftProperties {
        stiffness_nm_per_rad: [0.8, 1.2],
        rest_phase_rad: [0.05, -0.05],
        damping_nms_per_rad: [0.12, 0.18],
    }
}
fn energy(a: &RotorAssembly<2>, s: &RotorState<2>) -> f64 {
    let omega = v(s.rigid.omega_body);
    let velocity = v(s.rigid.vel_mps);
    let laws = properties();
    let mut result = 0.5 * a.reduced_body.mass_kg * velocity.dot(velocity)
        + 0.5 * omega.dot(a.reduced_body.inertia_tensor().matrix().apply(omega));
    for i in 0..2 {
        result += 0.5 * s.axial_momentum[i].powi(2) / a.rotors[i].axial_inertia_kgm2
            + 0.5
                * laws.stiffness_nm_per_rad[i]
                * (s.relative_phase_rad[i] - laws.rest_phase_rad[i]).powi(2);
    }
    result
}
fn momentum(a: &RotorAssembly<2>, s: &RotorState<2>, origin: [f64; 3]) -> (Vec3, Vec3) {
    let linear = v(s.rigid.vel_mps).scale(a.reduced_body.mass_kg);
    let mut spin = a
        .reduced_body
        .inertia_tensor()
        .matrix()
        .apply(v(s.rigid.omega_body));
    for i in 0..2 {
        spin = spin + v(a.rotors[i].axis_body).scale(s.axial_momentum[i]);
    }
    (
        linear,
        rotation(s).rotate(spin).unwrap() + (v(s.rigid.pos_m) - v(origin)).cross(linear),
    )
}

#[test]
fn g3_matching_snapshots_share_axial_spatial_and_swirling_aerodynamics() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::PropCouplingSolverSpec {
        cap: 64,
        tol: 1e-10,
        ..CANDIDATE_B
    };
    for axes in [AXIAL, [[0.96, 0.0, 0.28], [0.96, 0.0, -0.28]]] {
        let a = assembly(design.gross_mass_kg, axes);
        let mut locked = a.reduced_body.inertia_tensor().matrix();
        for rotor in a.rotors {
            for i in 0..3 {
                for j in 0..3 {
                    locked.m[3 * i + j] +=
                        rotor.axial_inertia_kgm2 * rotor.axis_body[i] * rotor.axis_body[j];
                }
            }
        }
        let gears = GearedAssembly::new(
            design.gross_mass_kg,
            locked,
            core::array::from_fn(|i| GearedRotor {
                axis_body: axes[i],
                axial_inertia_kgm2: a.rotors[i].axial_inertia_kgm2,
                ratio: if i == 0 { 1.0 } else { -1.0 },
            }),
        )
        .unwrap();
        let m = mounts(axes);
        let (passive, geared) = if axes == AXIAL {
            (
                ShaftPropulsion::new(&design, &a, m).unwrap(),
                GearedPropulsion::new(&design, &gears, m).unwrap(),
            )
        } else {
            (
                ShaftPropulsion::new_spatial(&design, &a, m).unwrap(),
                GearedPropulsion::new_spatial(&design, &gears, m).unwrap(),
            )
        };
        let state = initial(&a, [35.0, -35.0]);
        let geared_state = gears.state_from_speed(state.rigid, 0.2, 35.0).unwrap();
        let control = inputs(&a, &state);
        let wind = |point: [f64; 3]| Ok([0.02 * point[1], 0.0, 0.0]);
        let actual = passive.evaluate(&state, control, &wind, None).unwrap();
        let expected = geared
            .evaluate(
                &geared_state,
                PropulsionInputs {
                    density_kg_m3: control.density_kg_m3,
                    canard_rad: control.canard_rad,
                    motor_nm: 0.0,
                    gravity_world_mps2: control.gravity_world_mps2,
                },
                &wind,
                None,
            )
            .unwrap();
        assert_eq!(actual.coupled, expected.coupled);
        assert_eq!(actual.carrier, expected.loads.carrier);
        assert_eq!(
            actual.aerodynamic_force_body_n,
            expected.aerodynamic_force_body_n
        );
        assert_eq!(actual.torques.aerodynamic_nm, expected.loads.aerodynamic_nm);
        assert_eq!(
            actual.propeller_absolute_speed_rad_s,
            expected.propeller_absolute_speed_rad_s
        );
        assert_eq!(
            actual.torques.aerodynamic_nm,
            [-actual.coupled.torque_nm[0], actual.coupled.torque_nm[1]]
        );
        let actual_swirl = SwirlShaftPropulsion::new(&passive, Controls::default())
            .evaluate(&state, control, &wind, None)
            .unwrap();
        let expected_swirl = SwirlPropulsion::new(&geared, Controls::default())
            .evaluate(
                &geared_state,
                PropulsionInputs {
                    density_kg_m3: control.density_kg_m3,
                    canard_rad: control.canard_rad,
                    motor_nm: 0.0,
                    gravity_world_mps2: control.gravity_world_mps2,
                },
                &wind,
                None,
            )
            .unwrap();
        assert_eq!(actual_swirl.coupled, expected_swirl.coupled);
        assert_eq!(actual_swirl.carrier, expected_swirl.loads.carrier);
        assert_eq!(
            actual_swirl.torques.aerodynamic_nm,
            expected_swirl.loads.aerodynamic_nm
        );
        assert_ne!(
            actual_swirl.torques.aerodynamic_nm,
            actual.torques.aerodynamic_nm
        );
        for k in 0..2 {
            assert_eq!(
                actual_swirl.coupled.rotors[k].power.shaft_w,
                -actual_swirl.torques.aerodynamic_nm[k]
                    * actual_swirl.propeller_absolute_speed_rad_s[k]
            );
            assert!(actual_swirl.coupled.rotors[k].power.swirl_w > 0.0);
        }
    }
}

#[test]
fn g1_real_passive_propulsion_uses_both_stages_and_closes_physical_work_after_impulse() {
    let mut design = wright_openloop_v1();
    design.coupling = CANDIDATE_B;
    let a = assembly(design.gross_mass_kg, AXIAL);
    let model = ShaftPropulsion::new(&design, &a, mounts(AXIAL)).unwrap();
    let dt = 0.002;
    let mut motion = a
        .shaft_motion(&initial(&a, [32.0, -39.0]), properties(), dt)
        .unwrap();
    let origin = [0.3, -0.4, 0.5];
    motion
        .apply_impulse(
            RotorImpulses {
                carrier: CarrierImpulse::default(),
                motor_nms: [0.25, -0.125],
                aerodynamic_nms: [-0.0625, 0.03125],
            },
            origin,
        )
        .unwrap();
    let before = *motion.state();
    let mut stages = Vec::new();
    let wind = |point: [f64; 3]| Ok([0.02 * point[1], 0.0, 0.0]);
    let result = model
        .step(
            &mut motion,
            1.0,
            origin,
            None,
            |time, state| {
                stages.push((time, *state));
                Ok(inputs(&a, state))
            },
            |_, point| wind(point),
        )
        .unwrap();
    assert_eq!(stages.len(), 2);
    assert_eq!(stages[0], (1.0, before));
    assert_eq!(stages[1].0, 1.001);
    assert_ne!(stages[1].1.axial_momentum, before.axial_momentum);
    assert_ne!(stages[1].1.relative_phase_rad, before.relative_phase_rad);
    assert_eq!(
        result.initial,
        model
            .evaluate(&before, inputs(&a, &before), &wind, None)
            .unwrap()
    );
    assert_eq!(
        result.midpoint,
        model
            .evaluate(
                &stages[1].1,
                inputs(&a, &stages[1].1),
                &wind,
                Some(result.initial.coupled.w_slip)
            )
            .unwrap()
    );
    assert_ne!(
        result.initial.propeller_absolute_speed_rad_s,
        result.midpoint.propeller_absolute_speed_rad_s
    );
    assert!(
        (result.midpoint.coupled.thrust_n[0] - result.midpoint.coupled.thrust_n[1]).abs() > 1.0
    );
    let next = &result.motion.state;
    let (p0, l0) = momentum(&a, &before, origin);
    let (p1, l1) = momentum(&a, next, origin);
    let force_impulse = v(result.midpoint.carrier.force_n).scale(dt);
    let mut torque = v(result.midpoint.carrier.moment_nm);
    for i in 0..2 {
        torque = torque + v(a.rotors[i].axis_body).scale(result.midpoint.torques.aerodynamic_nm[i]);
    }
    let angular_impulse = (rotation(&before).rotate(torque).unwrap()
        + rotation(next).rotate(torque).unwrap())
    .scale(0.5 * dt)
        + (v(before.rigid.pos_m) + v(before.rigid.vel_mps).scale(0.5 * dt) - v(origin))
            .cross(force_impulse);
    assert!((p1 - p0 - force_impulse).norm() < 1e-10);
    assert!((l1 - l0 - angular_impulse).norm() < 1e-9);
    assert!((v(result.motion.external_angular_impulse_world) - angular_impulse).norm() < 1e-10);
    let motor_work: f64 = (0..2)
        .map(|i| {
            result.midpoint.torques.motor_nm[i]
                * (next.relative_phase_rad[i] - before.relative_phase_rad[i])
        })
        .sum();
    assert!((result.motion.motor_work_j - motor_work).abs() < 1e-11);
    let heat: f64 = result.motion.native.dissipated_energy_j.into_iter().sum();
    assert!(heat > 0.0);
    let balance = energy(&a, next) - energy(&a, &before) + heat
        - motor_work
        - result.motion.external_work_j
        - result.motion.native.motion.drift_energy_defect_j;
    assert!(
        balance.abs() < 1e-8,
        "native physical storage/work balance {balance}"
    );
    assert_eq!(motion.state(), next);
    assert_eq!(
        motion.native_state(),
        &result.motion.native.motion.state_after
    );
}

#[test]
fn g3_independent_spin_mapping_is_invariant_to_rotor_inventory_order() {
    let mut design = wright_openloop_v1();
    design.coupling = CANDIDATE_B;
    let a = assembly(design.gross_mass_kg, AXIAL);
    let state = initial(&a, [29.0, -42.0]);
    let model = ShaftPropulsion::new(&design, &a, mounts(AXIAL)).unwrap();
    let original = model
        .evaluate(&state, inputs(&a, &state), &|_| Ok(ZERO), None)
        .unwrap();
    let swapped = RotorAssembly {
        rotors: [a.rotors[1], a.rotors[0]],
        ..a
    };
    let reordered = RotorState {
        axial_momentum: [state.axial_momentum[1], state.axial_momentum[0]],
        relative_phase_rad: [state.relative_phase_rad[1], state.relative_phase_rad[0]],
        ..state
    };
    let mut mapping = mounts(AXIAL);
    mapping[0].rotor_index = 1;
    mapping[1].rotor_index = 0;
    let reordered_model = ShaftPropulsion::new(&design, &swapped, mapping).unwrap();
    let mut control = inputs(&a, &state);
    control.motor_nm.swap(0, 1);
    let actual = reordered_model
        .evaluate(&reordered, control, &|_| Ok(ZERO), None)
        .unwrap();
    assert_eq!(actual.coupled, original.coupled);
    assert_eq!(actual.carrier, original.carrier);
    assert_eq!(
        actual.propeller_absolute_speed_rad_s,
        original.propeller_absolute_speed_rad_s
    );
    assert_eq!(
        actual.torques.aerodynamic_nm,
        [
            original.torques.aerodynamic_nm[1],
            original.torques.aerodynamic_nm[0]
        ]
    );
    let mut changed_carrier = state;
    changed_carrier.rigid.omega_body[0] += 0.1;
    let same_spin = model
        .evaluate(
            &changed_carrier,
            inputs(&a, &changed_carrier),
            &|_| Ok(ZERO),
            None,
        )
        .unwrap();
    assert_eq!(
        same_spin.propeller_absolute_speed_rad_s, original.propeller_absolute_speed_rad_s,
        "absolute axial momentum already includes carrier rotation"
    );
}

#[test]
fn g0_inventory_spin_and_midpoint_wind_refusals_leave_native_motion_unchanged() {
    let mut design = wright_openloop_v1();
    design.coupling = CANDIDATE_B;
    let a = assembly(design.gross_mass_kg, AXIAL);
    let model = ShaftPropulsion::new(&design, &a, mounts(AXIAL)).unwrap();
    let state = initial(&a, [35.0, -35.0]);
    let mut stopped = state;
    stopped.axial_momentum[0] = 0.0;
    assert_eq!(
        model
            .evaluate(&stopped, inputs(&a, &state), &|_| Ok(ZERO), None)
            .unwrap_err()
            .code,
        "propeller-spin-outside-domain"
    );
    let mut other = a;
    other.rotors[0].axial_inertia_kgm2 *= 2.0;
    let mut incompatible = other.shaft_motion(&state, properties(), 0.002).unwrap();
    let mut calls = 0;
    assert_eq!(
        model
            .step(
                &mut incompatible,
                0.0,
                ZERO,
                None,
                |_, _| {
                    calls += 1;
                    Ok(inputs(&a, &state))
                },
                |_, _| Ok(ZERO)
            )
            .unwrap_err()
            .code,
        "propulsion-assembly-invalid"
    );
    assert_eq!(calls, 0);
    let mut motion = a.shaft_motion(&state, properties(), 0.002).unwrap();
    let before = *motion.state();
    let canonical = *motion.native_state();
    let mut controls = 0;
    let error = model
        .step(
            &mut motion,
            0.0,
            ZERO,
            None,
            |_, stage| {
                controls += 1;
                Ok(inputs(&a, stage))
            },
            |time, _| {
                if time > 0.0 {
                    Err(fs_flyer::Refusal {
                        code: "wind-domain",
                        message: "midpoint outside material wind domain".into(),
                        ranked_repairs: vec![],
                    })
                } else {
                    Ok(ZERO)
                }
            },
        )
        .unwrap_err();
    assert_eq!(error.code, "wind-domain");
    assert_eq!(controls, 2);
    assert_eq!(*motion.state(), before);
    assert_eq!(*motion.native_state(), canonical);
    let mut fresh = a.shaft_motion(&state, properties(), 0.002).unwrap();
    let retry = model
        .step(
            &mut motion,
            0.0,
            ZERO,
            None,
            |_, stage| Ok(inputs(&a, stage)),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    let expected = model
        .step(
            &mut fresh,
            0.0,
            ZERO,
            None,
            |_, stage| Ok(inputs(&a, stage)),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(motion.native_state(), fresh.native_state());
}

#[test]
fn g0_continuous_force_retains_rounding_without_relaxing_explicit_impulse_admission() {
    let assembly = RotorAssembly {
        reduced_body: RigidBody::from_inertia_tensor(
            3.0,
            Mat3 {
                m: [1.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 2.0],
            },
        )
        .unwrap(),
        rotors: [AxialRotor {
            axis_body: [1.0, 0.0, 0.0],
            axial_inertia_kgm2: 0.2,
        }],
    };
    let state = assembly
        .state_from_relative_spin(
            SixDofState {
                pos_m: ZERO,
                vel_mps: [12.0, 0.0, 0.0],
                quat: [1.0, 0.0, 0.0, 0.0],
                omega_body: ZERO,
            },
            [0.0],
            [0.0],
        )
        .unwrap();
    let dt = 0.002;
    let force = 0.001;
    let requested = force * dt;
    let mut motion = assembly
        .shaft_motion(&state, ShaftProperties::default(), dt)
        .unwrap();
    assert_eq!(
        motion
            .apply_impulse(
                RotorImpulses {
                    carrier: CarrierImpulse {
                        linear_world_ns: [requested, 0.0, 0.0],
                        angular_body_nms: ZERO,
                    },
                    motor_nms: [0.0],
                    aerodynamic_nms: [0.0],
                },
                ZERO,
            )
            .unwrap_err()
            .code,
        "impulse-unrepresentable"
    );
    assert_eq!(motion.state(), &state);
    let result = motion
        .advance_with(0.0, ZERO, |_, _| {
            Ok((
                Loads {
                    force_n: [force, 0.0, 0.0],
                    moment_nm: ZERO,
                },
                RotorTorques {
                    motor_nm: [0.0],
                    aerodynamic_nm: [0.0],
                },
            ))
        })
        .unwrap();
    let actual =
        assembly.reduced_body.mass_kg * (result.state.rigid.vel_mps[0] - state.rigid.vel_mps[0]);
    let defect = actual - requested;
    assert_ne!(actual, 0.0);
    assert!(defect.abs() > 4096.0 * f64::EPSILON * requested);
    assert!(defect.abs() < 64.0 * f64::EPSILON * 36.0);
    assert_eq!(result.linear_momentum_residual, [defect, 0.0, 0.0]);
    let mean_velocity = 0.5 * (result.state.rigid.vel_mps[0] + state.rigid.vel_mps[0]);
    assert!((result.work_energy_residual_j - defect * mean_velocity).abs() < 1e-18);
    let before = *motion.state();
    let canonical = *motion.native_state();
    assert_eq!(
        motion
            .advance_with(dt, ZERO, |_, _| {
                Ok((
                    Loads {
                        force_n: [1e-300, 0.0, 0.0],
                        moment_nm: ZERO,
                    },
                    RotorTorques {
                        motor_nm: [0.0],
                        aerodynamic_nm: [0.0],
                    },
                ))
            })
            .unwrap_err()
            .code,
        "impulse-unrepresentable"
    );
    assert_eq!(motion.state(), &before);
    assert_eq!(motion.native_state(), &canonical);
}

#[test]
fn g1_swirl_drives_independent_springs_with_physical_momentum_and_work_after_impulse() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::PropCouplingSolverSpec {
        cap: 64,
        tol: 1e-10,
        ..CANDIDATE_B
    };
    let a = assembly(design.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new(&design, &a, mounts(AXIAL)).unwrap();
    let model = SwirlShaftPropulsion::new(&base, Controls::default());
    let dt = 0.002;
    let mut motion = a
        .shaft_motion(&initial(&a, [32.0, -39.0]), properties(), dt)
        .unwrap();
    let origin = [0.3, -0.4, 0.5];
    motion
        .apply_impulse(
            RotorImpulses {
                carrier: CarrierImpulse::default(),
                motor_nms: [0.25, -0.125],
                aerodynamic_nms: [-0.0625, 0.03125],
            },
            origin,
        )
        .unwrap();
    let before = *motion.state();
    let mut stages = Vec::new();
    let wind = |point: [f64; 3]| Ok([0.02 * point[1], 0.0, 0.0]);
    let result = model
        .step(
            &mut motion,
            1.0,
            origin,
            None,
            |time, state| {
                stages.push((time, *state));
                Ok(inputs(&a, state))
            },
            |_, point| wind(point),
        )
        .unwrap();
    assert_eq!(stages.len(), 2);
    assert_eq!(stages[0], (1.0, before));
    assert_eq!(stages[1].0, 1.001);
    assert_ne!(stages[1].1.axial_momentum, before.axial_momentum);
    assert_eq!(
        result.initial,
        model
            .evaluate(&before, inputs(&a, &before), &wind, None)
            .unwrap()
    );
    let seed = result.initial.coupled.warm_start();
    assert_eq!(
        result.midpoint,
        model
            .evaluate(&stages[1].1, inputs(&a, &stages[1].1), &wind, Some(&seed),)
            .unwrap()
    );
    assert_ne!(
        result.initial.propeller_absolute_speed_rad_s,
        result.midpoint.propeller_absolute_speed_rad_s
    );
    let next = &result.motion.state;
    let (p0, l0) = momentum(&a, &before, origin);
    let (p1, l1) = momentum(&a, next, origin);
    let impulse = v(result.midpoint.carrier.force_n).scale(dt);
    let mut torque = v(result.midpoint.carrier.moment_nm);
    for i in 0..2 {
        torque = torque + v(a.rotors[i].axis_body).scale(result.midpoint.torques.aerodynamic_nm[i]);
        assert_eq!(
            result.midpoint.coupled.rotors[i].power.shaft_w,
            -result.midpoint.torques.aerodynamic_nm[i]
                * result.midpoint.propeller_absolute_speed_rad_s[i]
        );
    }
    let angular_impulse = (rotation(&before).rotate(torque).unwrap()
        + rotation(next).rotate(torque).unwrap())
    .scale(0.5 * dt)
        + (v(before.rigid.pos_m) + v(before.rigid.vel_mps).scale(0.5 * dt) - v(origin))
            .cross(impulse);
    assert!((p1 - p0 - impulse).norm() < 1e-10);
    assert!((l1 - l0 - angular_impulse).norm() < 1e-9);
    let motor_work: f64 = (0..2)
        .map(|i| {
            result.midpoint.torques.motor_nm[i]
                * (next.relative_phase_rad[i] - before.relative_phase_rad[i])
        })
        .sum();
    assert!((result.motion.motor_work_j - motor_work).abs() < 1e-11);
    let heat: f64 = result.motion.native.dissipated_energy_j.into_iter().sum();
    assert!(heat > 0.0);
    assert_ne!(
        result.motion.elastic_energy_before_j,
        result.motion.elastic_energy_after_j
    );
    let balance = energy(&a, next) - energy(&a, &before) + heat
        - motor_work
        - result.motion.external_work_j
        - result.motion.native.motion.drift_energy_defect_j;
    assert!(
        balance.abs() < 1e-8,
        "physical swirling-shaft work balance {balance}"
    );
    assert_eq!(motion.state(), next);
    assert_eq!(
        motion.native_state(),
        &result.motion.native.motion.state_after
    );
}

/// Same native attached windmilling construction as the geared regression:
/// w=-2, u=-0.15 at V=10 and Omega=60, satisfying both annular momentum laws.
fn windmilling_rotor() -> fs_airscrew::Rotor {
    let stations = [0.3, 0.5, 0.7, 0.9]
        .into_iter()
        .map(|r| {
            let va = 8.0_f64;
            let vt = 60.0 * r + 0.15;
            let phi = va.atan2(vt);
            let ratio = (phi.cos() + 0.075 * phi.sin()) / (0.075 * phi.cos() - phi.sin());
            let b = std::f64::consts::TAU;
            let alpha = 0.024 * ratio / (b + (b * b - 4.0 * 0.8 * 0.012 * ratio * ratio).sqrt());
            let loss = |distance: f64| {
                2.0 / std::f64::consts::PI * (-1.5 * distance / (r * phi.sin())).exp().acos()
            };
            let f = (loss(1.0 - r) * loss(r - 0.24)).clamp(1e-3, 1.0);
            let target = 4.0 * std::f64::consts::PI * 1.225 * r * f * va * (-2.0);
            let unit_force = 1.5
                * 1.225
                * (va * va + vt * vt)
                * (b * alpha * phi.cos() - (0.012 + 0.8 * alpha * alpha) * phi.sin());
            fs_airscrew::BladeStation {
                r_over_r: r,
                chord_m: target / unit_force,
                beta_rad: phi + alpha,
            }
        })
        .collect();
    fs_airscrew::Rotor {
        radius_m: 1.0,
        n_blades: 3,
        camber_ratio: 0.0,
        stations,
    }
}

#[test]
fn g1_swirl_windmilling_accelerates_both_free_shafts_without_motor_work() {
    let mut design = wright_openloop_v1();
    design.camber_ratio = 0.0;
    design.rotor = windmilling_rotor();
    design.coupling = fs_flyer::propcoupling::PropCouplingSolverSpec {
        cap: 64,
        tol: 1e-10,
        ..CANDIDATE_B
    };
    let a = assembly(design.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new(&design, &a, mounts(AXIAL)).unwrap();
    let model = SwirlShaftPropulsion::new(&base, Controls::default());
    let mut rigid = initial(&a, [60.0, -60.0]).rigid;
    rigid.vel_mps = [10.0, 0.0, 0.0];
    rigid.omega_body = ZERO;
    let start = a
        .state_from_relative_spin(rigid, [0.0; 2], [60.0, -60.0])
        .unwrap();
    let dt = 0.001;
    let mut motion = a
        .shaft_motion(&start, ShaftProperties::default(), dt)
        .unwrap();
    let result = model
        .step(
            &mut motion,
            0.0,
            start.rigid.pos_m,
            None,
            |_, _| {
                Ok(ShaftInputs {
                    density_kg_m3: 1.225,
                    canard_rad: 0.0,
                    motor_nm: [0.0; 2],
                    gravity_world_mps2: ZERO,
                })
            },
            |_, _| Ok(ZERO),
        )
        .unwrap();
    let mut spin_energy_gain = 0.0;
    for k in 0..2 {
        for stage in [&result.initial, &result.midpoint] {
            let rotor = &stage.coupled.rotors[k];
            assert!(rotor.thrust_n < 0.0 && rotor.torque_nm < 0.0);
            assert!(rotor.power.shaft_w < 0.0 && rotor.power.swirl_w > 0.0);
            assert!(
                stage.torques.aerodynamic_nm[k] * stage.propeller_absolute_speed_rad_s[k] > 0.0
            );
        }
        let p0 = start.axial_momentum[k];
        let p1 = result.motion.state.axial_momentum[k];
        let applied = result.midpoint.torques.aerodynamic_nm[k] * dt;
        assert!((p1 - p0 - applied).abs() < 1e-12);
        assert!(p1.abs() > p0.abs());
        spin_energy_gain += (p1 - p0) * (p1 + p0) / (2.0 * a.rotors[k].axial_inertia_kgm2);
    }
    assert!(spin_energy_gain > 0.0);
    assert_eq!(result.motion.motor_work_j, 0.0);
    assert_eq!(result.motion.native.dissipated_energy_j, [0.0; 2]);
    assert!(result.motion.state.rigid.vel_mps[0] < start.rigid.vel_mps[0]);
    assert!(v(result.motion.angular_momentum_residual).norm() < 1e-10);
}

#[test]
fn g0_swirl_budget_and_midpoint_refusals_preserve_both_states_and_retry() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::PropCouplingSolverSpec {
        cap: 64,
        tol: 1e-10,
        ..CANDIDATE_B
    };
    let a = assembly(design.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new(&design, &a, mounts(AXIAL)).unwrap();
    let model = SwirlShaftPropulsion::new(&base, Controls::default());
    let start = initial(&a, [35.0, -35.0]);
    let mut motion = a.shaft_motion(&start, properties(), 0.002).unwrap();
    let prepared = *motion.state();
    let canonical = *motion.native_state();
    let bounded = SwirlShaftPropulsion::new(
        &base,
        Controls {
            maximum_station_evaluations: 3,
            ..Controls::default()
        },
    );
    assert_eq!(
        bounded
            .step(
                &mut motion,
                0.0,
                ZERO,
                None,
                |_, state| Ok(inputs(&a, state)),
                |_, _| Ok(ZERO),
            )
            .unwrap_err()
            .code,
        "swirl-nonconvergence"
    );
    assert_eq!(motion.state(), &prepared);
    assert_eq!(motion.native_state(), &canonical);
    let seed = model
        .evaluate(&start, inputs(&a, &start), &|_| Ok(ZERO), None)
        .unwrap()
        .coupled
        .warm_start();
    let saved_seed = seed.clone();
    let mut calls = 0;
    let error = model
        .step(
            &mut motion,
            0.0,
            ZERO,
            Some(&seed),
            |_, state| {
                calls += 1;
                Ok(inputs(&a, state))
            },
            |time, _| {
                if time > 0.0 {
                    Err(fs_flyer::Refusal {
                        code: "wind-domain",
                        message: "midpoint outside material wind domain".into(),
                        ranked_repairs: vec![],
                    })
                } else {
                    Ok(ZERO)
                }
            },
        )
        .unwrap_err();
    assert_eq!(error.code, "wind-domain");
    assert_eq!(calls, 2);
    assert_eq!(motion.state(), &prepared);
    assert_eq!(motion.native_state(), &canonical);
    assert_eq!(seed, saved_seed);
    let mut fresh = a.shaft_motion(&start, properties(), 0.002).unwrap();
    let retry = model
        .step(
            &mut motion,
            0.0,
            ZERO,
            Some(&seed),
            |_, state| Ok(inputs(&a, state)),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    let expected = model
        .step(
            &mut fresh,
            0.0,
            ZERO,
            Some(&seed),
            |_, state| Ok(inputs(&a, state)),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(motion.native_state(), fresh.native_state());
}
