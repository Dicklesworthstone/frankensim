//! Native modal/skew wake memory through independent compliant shafts.
use super::*;
use fs_airscrew::dynamic_inflow::first_harmonic::HarmonicInduction;
use fs_airscrew::dynamic_inflow::first_harmonic::skew::SkewInflowControls;
use fs_flyer::spine::rotors::propulsion::dynamic::first_harmonic::{
    HarmonicPropulsion, HarmonicState,
};
use fs_flyer::spine::rotors::propulsion::passive::dynamic::first_harmonic::{
    HarmonicShaftPropulsion, HarmonicShaftState,
};

fn modes() -> [HarmonicInduction; 2] {
    [
        HarmonicInduction {
            mean_mps: 0.4,
            first_mps: [0.10, -0.05],
        },
        HarmonicInduction {
            mean_mps: 0.8,
            first_mps: [-0.03, 0.07],
        },
    ]
}

fn control(a: &RotorAssembly<2>, t: f64, state: &HarmonicShaftState<2>) -> ShaftInputs<2> {
    let mut value = inputs(a, &state.motion);
    value.motor_nm[0] += 20.0 * t;
    value.motor_nm[1] -= 10.0 * t;
    value
}

fn model<'p, 'a>(
    base: &'p ShaftPropulsion<'a, 2>,
    skew: bool,
) -> HarmonicShaftPropulsion<'p, 'a, 2> {
    let model = if skew {
        HarmonicShaftPropulsion::new_skewed(base, RHO, 16, SkewInflowControls::default())
    } else {
        HarmonicShaftPropulsion::new(base, RHO, 16)
    };
    model
        .unwrap()
        .with_blade_incident(BladeIncident::WindAndWingField)
        .unwrap()
}

#[test]
fn g3_independent_harmonic_and_skew_snapshots_match_the_original_geared_owner() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, CANTED);
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
        core::array::from_fn(|k| GearedRotor {
            axis_body: CANTED[k],
            axial_inertia_kgm2: a.rotors[k].axial_inertia_kgm2,
            ratio: [1.0, -1.0][k],
        }),
    )
    .unwrap();
    let passive = ShaftPropulsion::new_spatial(&design, &a, mounts(CANTED)).unwrap();
    let geared = GearedPropulsion::new_spatial(&design, &gears, mounts(CANTED)).unwrap();
    let state = HarmonicShaftState {
        motion: initial(&a, [35.0, -35.0]),
        induction: modes(),
    };
    let expected_state = HarmonicState {
        motion: gears
            .state_from_speed(state.motion.rigid, 0.2, 35.0)
            .unwrap(),
        induction: state.induction,
    };
    for skew in [false, true] {
        let actual_model = model(&passive, skew);
        let expected_model = if skew {
            HarmonicPropulsion::new_skewed(&geared, RHO, 16, SkewInflowControls::default())
        } else {
            HarmonicPropulsion::new(&geared, RHO, 16)
        }
        .unwrap()
        .with_blade_incident(BladeIncident::WindAndWingField)
        .unwrap();
        assert_eq!(actual_model.disk_frames(), expected_model.disk_frames());
        assert_eq!(actual_model.inflow_model(), expected_model.inflow_model());
        let values = control(&a, 0.0, &state);
        let actual = actual_model
            .evaluate(&state, values, &|p| wind(0.0, p))
            .unwrap();
        let expected = expected_model
            .evaluate(
                &expected_state,
                PropulsionInputs {
                    density_kg_m3: values.density_kg_m3,
                    canard_rad: values.canard_rad,
                    motor_nm: 0.0,
                    gravity_world_mps2: values.gravity_world_mps2,
                },
                &|p| wind(0.0, p),
            )
            .unwrap();
        assert_eq!(actual.coupled, expected.aerodynamics);
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
            [
                -actual.coupled.propellers[0].torque_nm,
                actual.coupled.propellers[1].torque_nm
            ]
        );
        assert_eq!(actual.coupled.skew_transport.is_some(), skew);
    }
}

#[test]
fn g1_asymmetric_modal_memory_and_independent_shafts_share_real_stages_and_source_accounts() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, CANTED);
    let base = ShaftPropulsion::new_spatial(&design, &a, mounts(CANTED)).unwrap();
    let start = initial(&a, [32.0, -39.0]);
    let dt = 0.004;
    let origin = [0.3, -0.4, 0.5];
    for skew in [false, true] {
        let model = model(&base, skew);
        let mut motion = a.shaft_motion(&start, properties(), dt).unwrap();
        let mut half = a.shaft_motion(&start, properties(), 0.5 * dt).unwrap();
        let mut full = a.shaft_motion(&start, properties(), dt).unwrap();
        let event = RotorImpulses {
            carrier: CarrierImpulse::default(),
            motor_nms: [0.25, -0.125],
            aerodynamic_nms: [-0.0625, 0.03125],
        };
        for owner in [&mut motion, &mut half, &mut full] {
            owner.apply_impulse(event, origin).unwrap();
        }
        let before = *motion.state();
        let before_modes = modes();
        let mut induction = before_modes;
        let mut stages = Vec::new();
        let result = model
            .step(
                &mut motion,
                &mut induction,
                0.0,
                origin,
                |time, state| {
                    stages.push((time, *state));
                    Ok(control(&a, time, state))
                },
                wind,
            )
            .unwrap();
        assert_eq!(stages.len(), 2);
        assert_eq!(
            stages[0].1,
            HarmonicShaftState {
                motion: before,
                induction: before_modes
            }
        );
        assert_eq!(stages[1].0, 0.5 * dt);
        let predicted = half
            .advance(result.initial.carrier, result.initial.torques, origin)
            .unwrap();
        assert_eq!(stages[1].1.motion, predicted.state);
        assert_ne!(stages[1].1.motion.axial_momentum, before.axial_momentum);
        assert_ne!(stages[1].1.induction, before_modes);
        assert_eq!(
            result.midpoint,
            model
                .evaluate(
                    &stages[1].1,
                    control(&a, 0.5 * dt, &stages[1].1),
                    &|p| wind(0.5 * dt, p),
                )
                .unwrap()
        );
        let expected = full
            .advance(result.midpoint.carrier, result.midpoint.torques, origin)
            .unwrap();
        assert_eq!(result.motion, expected);
        assert_eq!(motion.native_state(), full.native_state());
        assert_eq!(
            result.state,
            HarmonicShaftState {
                motion: *motion.state(),
                induction
            }
        );

        let native = model.inflow_model();
        let radius = design.rotor.radius_m;
        let masses = [
            (8.0 / 3.0) * RHO * radius.powi(3),
            (16.0 / 45.0) * RHO * radius.powi(3),
            (16.0 / 45.0) * RHO * radius.powi(3),
        ];
        for k in 0..2 {
            let held = |stage: &fs_flyer::spine::rotors::propulsion::dynamic::first_harmonic::HarmonicAerodynamics, duration| {
                if skew {
                    native.step_skewed_held(
                        before_modes[k], stage.disk_incident_mps[k], stage.wake_force_n[k],
                        stage.frame_rate_rad_s[k], duration, SkewInflowControls::default(),
                    ).unwrap().modal
                } else {
                    native.step_held(
                        before_modes[k], stage.disk_axial_mps[k], stage.wake_force_n[k],
                        stage.frame_rate_rad_s[k], duration,
                    ).unwrap()
                }
            };
            assert_eq!(
                stages[1].1.induction[k],
                held(&result.initial.coupled, 0.5 * dt).induction
            );
            assert_eq!(result.inflow[k], held(&result.midpoint.coupled, dt));
            let receipt = result.inflow[k];
            let u0 = before_modes[k].modes();
            let u1 = induction[k].modes();
            let mut kinetic_change = 0.0;
            for j in 0..3 {
                let balance = masses[j] * (u1[j] - u0[j])
                    - dt * result.midpoint.coupled.wake_force_n[k][j]
                    + receipt.transport_impulse_ns[j]
                    + receipt.frame_impulse_ns[j];
                assert!(balance.abs() < 1e-11, "mode {k}/{j}: {balance}");
                kinetic_change += 0.5 * masses[j] * (u1[j].powi(2) - u0[j].powi(2));
            }
            let balance = kinetic_change - receipt.induced_work_j
                + receipt.transport_work_j
                + receipt.frame_work_j;
            assert!(balance.abs() < 1e-11, "modal storage work: {balance}");
            assert!(receipt.frame_work_j.abs() < 1e-13);
            assert_eq!(
                result.coupling_work_residual_j[k],
                receipt.induced_work_j - dt * result.midpoint.coupled.propellers[k].induced_work_w
            );
            if skew {
                let detail = result.skew_inflow.unwrap()[k];
                assert_eq!(detail.modal, receipt);
                assert!(detail.cross_transport_work_j.abs() > 1e-12);
                assert!(
                    (detail.diagonal_transport_work_j + detail.cross_transport_work_j
                        - receipt.transport_work_j)
                        .abs()
                        < 1e-12
                );
            }
        }
        assert_eq!(result.skew_inflow.is_some(), skew);
        let (p0, l0) = momentum(&a, &before, origin);
        let (p1, l1) = momentum(&a, motion.state(), origin);
        let impulse = v(result.midpoint.carrier.force_n).scale(dt);
        let mut torque = v(result.midpoint.carrier.moment_nm);
        for k in 0..2 {
            torque =
                torque + v(a.rotors[k].axis_body).scale(result.midpoint.torques.aerodynamic_nm[k]);
        }
        let angular = (rotation(&before).rotate(torque).unwrap()
            + rotation(motion.state()).rotate(torque).unwrap())
        .scale(0.5 * dt)
            + (v(before.rigid.pos_m) + v(before.rigid.vel_mps).scale(0.5 * dt) - v(origin))
                .cross(impulse);
        assert!((p1 - p0 - impulse).norm() < 1e-10);
        assert!((l1 - l0 - angular).norm() < 1e-9);
        let motor_work: f64 = (0..2)
            .map(|k| {
                result.midpoint.torques.motor_nm[k]
                    * (motion.state().relative_phase_rad[k] - before.relative_phase_rad[k])
            })
            .sum();
        let heat: f64 = result.motion.native.dissipated_energy_j.into_iter().sum();
        assert!(heat > 0.0);
        assert!((result.motion.motor_work_j - motor_work).abs() < 1e-11);
        let balance = energy(&a, motion.state()) - energy(&a, &before) + heat
            - motor_work
            - result.motion.external_work_j
            - result.motion.native.motion.drift_energy_defect_j;
        assert!(balance.abs() < 1e-8, "mechanical source balance: {balance}");
    }
}

#[test]
fn g0_harmonic_and_skew_late_refusals_preserve_every_mode_and_canonical_state() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, CANTED);
    let base = ShaftPropulsion::new_spatial(&design, &a, mounts(CANTED)).unwrap();
    let start = initial(&a, [32.0, -39.0]);
    for skew in [false, true] {
        let model = model(&base, skew);
        for late_native in [false, true] {
            let mut motion = a.shaft_motion(&start, properties(), 0.004).unwrap();
            let before = *motion.state();
            let canonical = *motion.native_state();
            let mut induction = modes();
            let mut calls = 0;
            let error = model
                .step(
                    &mut motion,
                    &mut induction,
                    0.0,
                    ZERO,
                    |time, state| {
                        calls += 1;
                        let mut values = control(&a, time, state);
                        if late_native && time > 0.0 {
                            values.motor_nm[0] = f64::MAX / 32.0;
                        }
                        Ok(values)
                    },
                    |time, point| {
                        if !late_native && time > 0.0 {
                            Err(fs_flyer::Refusal {
                                code: "wind-domain",
                                message: "midpoint blade field unavailable".into(),
                                ranked_repairs: vec![],
                            })
                        } else {
                            wind(time, point)
                        }
                    },
                )
                .unwrap_err();
            assert_eq!(calls, 2);
            assert_eq!(
                error.code,
                if late_native {
                    "rotor-shaft-refused"
                } else {
                    "wind-domain"
                }
            );
            assert_eq!(motion.state(), &before);
            assert_eq!(motion.native_state(), &canonical);
            assert_eq!(induction, modes());
            let mut fresh = a.shaft_motion(&start, properties(), 0.004).unwrap();
            let mut fresh_induction = modes();
            let retry = model
                .step(
                    &mut motion,
                    &mut induction,
                    0.0,
                    ZERO,
                    |t, state| Ok(control(&a, t, state)),
                    wind,
                )
                .unwrap();
            let expected = model
                .step(
                    &mut fresh,
                    &mut fresh_induction,
                    0.0,
                    ZERO,
                    |t, state| Ok(control(&a, t, state)),
                    wind,
                )
                .unwrap();
            assert_eq!(retry, expected);
            assert_eq!(motion.native_state(), fresh.native_state());
            assert_eq!(induction, fresh_induction);
        }
    }
}
