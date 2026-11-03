//! G0/G1/G3 joint physical wake memory and independent compliant shafts.
use super::*;
use fs_flyer::spine::rotors::propulsion::dynamic::{
    BladeIncident, DynamicPropulsion, DynamicState,
};
use fs_flyer::spine::rotors::propulsion::passive::dynamic::{
    DynamicShaftPropulsion, DynamicShaftState,
};

const RHO: f64 = 1.225;
const CANTED: [[f64; 3]; 2] = [[0.96, 0.0, 0.28], [0.96, 0.0, -0.28]];

fn wind(time: f64, point: [f64; 3]) -> Result<[f64; 3], fs_flyer::Refusal> {
    Ok([0.02 * point[1] + 0.03 * time, 0.01 * point[2], 0.0])
}

fn control(a: &RotorAssembly<2>, time: f64, state: &DynamicShaftState<2>) -> ShaftInputs<2> {
    let mut values = inputs(a, &state.motion);
    values.motor_nm[0] += 20.0 * time;
    values.motor_nm[1] -= 10.0 * time;
    values
}

#[test]
fn g3_dynamic_snapshots_share_uniform_oblique_and_material_air_owners() {
    let design = wright_openloop_v1();
    for axes in [AXIAL, CANTED] {
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
            core::array::from_fn(|k| GearedRotor {
                axis_body: axes[k],
                axial_inertia_kgm2: a.rotors[k].axial_inertia_kgm2,
                ratio: [1.0, -1.0][k],
            }),
        )
        .unwrap();
        let passive = ShaftPropulsion::new_spatial(&design, &a, mounts(axes)).unwrap();
        let geared = GearedPropulsion::new_spatial(&design, &gears, mounts(axes)).unwrap();
        let state = DynamicShaftState {
            motion: initial(&a, [35.0, -35.0]),
            induced_mps: [0.4, 0.8],
        };
        let geared_state = DynamicState {
            motion: gears
                .state_from_speed(state.motion.rigid, 0.2, 35.0)
                .unwrap(),
            induced_mps: state.induced_mps,
        };
        let inputs = control(&a, 0.0, &state);
        for oblique in [false, true] {
            let (actual, expected) = if oblique {
                (
                    DynamicShaftPropulsion::new_oblique(&passive, RHO, 16)
                        .unwrap()
                        .with_blade_incident(BladeIncident::WindAndWingField)
                        .unwrap(),
                    DynamicPropulsion::new_oblique(&geared, RHO, 16)
                        .unwrap()
                        .with_blade_incident(BladeIncident::WindAndWingField)
                        .unwrap(),
                )
            } else {
                (
                    DynamicShaftPropulsion::new(&passive, RHO).unwrap(),
                    DynamicPropulsion::new(&geared, RHO).unwrap(),
                )
            };
            let actual = actual.evaluate(&state, inputs, &|p| wind(0.0, p)).unwrap();
            let expected = expected
                .evaluate(
                    &geared_state,
                    PropulsionInputs {
                        density_kg_m3: inputs.density_kg_m3,
                        canard_rad: inputs.canard_rad,
                        motor_nm: 0.0,
                        gravity_world_mps2: inputs.gravity_world_mps2,
                    },
                    &|p| wind(0.0, p),
                )
                .unwrap();
            assert_eq!(actual.coupled, expected.aerodynamics);
            assert_eq!(actual.carrier, expected.loads.carrier);
            assert_eq!(actual.torques.aerodynamic_nm, expected.loads.aerodynamic_nm);
            assert_eq!(
                actual.propeller_absolute_speed_rad_s,
                expected.propeller_absolute_speed_rad_s
            );
            assert_eq!(actual.torques.motor_nm, inputs.motor_nm);
            assert_eq!(
                actual.torques.aerodynamic_nm,
                [
                    -actual.coupled.propellers[0].torque_nm,
                    actual.coupled.propellers[1].torque_nm
                ]
            );
        }
    }
}

#[test]
fn g1_native_joint_midpoint_after_impulse_closes_physical_work_and_wake_balances() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, CANTED);
    let base = ShaftPropulsion::new_spatial(&design, &a, mounts(CANTED)).unwrap();
    let model = DynamicShaftPropulsion::new_oblique(&base, RHO, 16)
        .unwrap()
        .with_blade_incident(BladeIncident::WindAndWingField)
        .unwrap();
    let start = initial(&a, [32.0, -39.0]);
    let dt = 0.004;
    let origin = [0.3, -0.4, 0.5];
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
    assert_eq!(motion.native_state(), half.native_state());
    let before = *motion.state();
    let before_wake = [0.2, 0.8];
    let mut induced = before_wake;
    let mut stages = Vec::new();
    let result = model
        .step(
            &mut motion,
            &mut induced,
            0.0,
            origin,
            |t, state| {
                stages.push((t, *state));
                Ok(control(&a, t, state))
            },
            wind,
        )
        .unwrap();
    assert_eq!(stages.len(), 2);
    assert_eq!(stages[0].0, 0.0);
    assert_eq!(stages[0].1.motion, before);
    assert_eq!(stages[0].1.induced_mps, before_wake);
    assert_eq!(stages[1].0, 0.5 * dt);
    let prediction = half
        .advance(result.initial.carrier, result.initial.torques, origin)
        .unwrap();
    assert_eq!(stages[1].1.motion, prediction.state);
    for k in 0..2 {
        let prediction = model
            .inflow_model()
            .step_held(
                before_wake[k],
                result.initial.coupled.disk_axial_mps[k],
                result.initial.coupled.propellers[k].thrust_n,
                0.5 * dt,
            )
            .unwrap();
        assert_eq!(stages[1].1.induced_mps[k], prediction.induced_mps);
        assert_ne!(stages[1].1.induced_mps[k], before_wake[k]);
    }
    assert_eq!(
        result.midpoint,
        model
            .evaluate(
                &stages[1].1,
                control(&a, 0.5 * dt, &stages[1].1),
                &|p| wind(0.5 * dt, p)
            )
            .unwrap()
    );
    let expected = full
        .advance(result.midpoint.carrier, result.midpoint.torques, origin)
        .unwrap();
    assert_eq!(
        result.motion, expected,
        "predictor heat/work cannot enter the accepted account"
    );
    assert_eq!(motion.native_state(), full.native_state());
    assert_eq!(result.state.motion, *motion.state());
    assert_eq!(result.state.induced_mps, induced);

    // Independent world momentum from physical reduced inertia, each actual
    // rotor momentum and the COM orbital contribution; internal motors cancel.
    let (p0, l0) = momentum(&a, &before, origin);
    let (p1, l1) = momentum(&a, motion.state(), origin);
    let force_impulse = v(result.midpoint.carrier.force_n).scale(dt);
    let mut torque = v(result.midpoint.carrier.moment_nm);
    for k in 0..2 {
        torque = torque + v(a.rotors[k].axis_body).scale(result.midpoint.torques.aerodynamic_nm[k]);
    }
    let angular_impulse = (rotation(&before).rotate(torque).unwrap()
        + rotation(motion.state()).rotate(torque).unwrap())
    .scale(0.5 * dt)
        + (v(before.rigid.pos_m) + v(before.rigid.vel_mps).scale(0.5 * dt) - v(origin))
            .cross(force_impulse);
    assert!((p1 - p0 - force_impulse).norm() < 1e-10);
    assert!((l1 - l0 - angular_impulse).norm() < 1e-9);
    let motor_work: f64 = (0..2)
        .map(|k| {
            result.midpoint.torques.motor_nm[k]
                * (motion.state().relative_phase_rad[k] - before.relative_phase_rad[k])
        })
        .sum();
    assert!((result.motion.motor_work_j - motor_work).abs() < 1e-11);
    let heat: f64 = result.motion.native.dissipated_energy_j.into_iter().sum();
    assert!(heat > 0.0);
    let mechanical_balance = energy(&a, motion.state()) - energy(&a, &before) + heat
        - motor_work
        - result.motion.external_work_j
        - result.motion.native.motion.drift_energy_defect_j;
    assert!(
        mechanical_balance.abs() < 1e-8,
        "physical mechanical balance {mechanical_balance:e}"
    );

    // Reconstruct the apparent-fluid mass, transported momentum and work
    // directly from disk geometry and accepted velocities, not receipt totals.
    let radius = design.rotor.radius_m;
    let mass = (8.0 / 3.0) * RHO * radius.powi(3);
    let area = core::f64::consts::PI * radius.powi(2);
    for k in 0..2 {
        let midpoint = 0.5 * (before_wake[k] + induced[k]);
        let thrust = result.midpoint.coupled.propellers[k].thrust_n;
        let transport =
            2.0 * RHO * area * midpoint * (result.midpoint.coupled.disk_axial_mps[k] + midpoint);
        let momentum = mass * (induced[k] - before_wake[k]) - dt * thrust + dt * transport;
        let energy = 0.5 * mass * (induced[k].powi(2) - before_wake[k].powi(2))
            - dt * thrust * midpoint
            + dt * transport * midpoint;
        assert!(momentum.abs() < 1e-12);
        assert!(energy.abs() < 1e-12);
        assert!((result.inflow[k].thrust_impulse_ns - dt * thrust).abs() < 1e-13);
        assert!((result.inflow[k].transport_work_j - dt * transport * midpoint).abs() < 1e-13);
    }
}

#[test]
fn g1_body_independent_shafts_and_wake_trajectory_refine_together() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, CANTED);
    let base = ShaftPropulsion::new_spatial(&design, &a, mounts(CANTED)).unwrap();
    let model = DynamicShaftPropulsion::new(&base, RHO).unwrap();
    let start = initial(&a, [32.0, -39.0]);
    let mut endpoints = Vec::new();
    for count in [12, 24, 48, 96] {
        let dt = 0.12 / f64::from(count);
        let mut motion = a.shaft_motion(&start, properties(), dt).unwrap();
        let mut induced = [0.2, 0.8];
        for tick in 0..count {
            model
                .step(
                    &mut motion,
                    &mut induced,
                    f64::from(tick) * dt,
                    ZERO,
                    |t, stage| Ok(control(&a, t, stage)),
                    wind,
                )
                .unwrap();
        }
        let s = motion.state();
        endpoints.push([
            s.rigid
                .pos_m
                .into_iter()
                .chain(s.rigid.vel_mps)
                .chain(s.rigid.omega_body)
                .chain(s.rigid.quat)
                .collect::<Vec<_>>(),
            s.relative_phase_rad
                .into_iter()
                .chain(core::array::from_fn::<_, 2, _>(|k| {
                    s.axial_momentum[k] / a.rotors[k].axial_inertia_kgm2
                }))
                .collect(),
            induced.to_vec(),
        ]);
    }
    for group in 0..3 {
        let difference: Vec<f64> = endpoints
            .windows(2)
            .map(|pair| {
                pair[0][group]
                    .iter()
                    .zip(&pair[1][group])
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .collect();
        println!("joint passive body/shaft/wake group {group} differences: {difference:?}");
        for pair in difference.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!(
                (3.5..4.5).contains(&ratio),
                "group {group} refinement ratio {ratio}"
            );
        }
    }
}

#[test]
fn g3_reversing_one_shaft_coordinate_preserves_oblique_wake_motion() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, CANTED);
    let mut b = a;
    b.rotors[0].axis_body = b.rotors[0].axis_body.map(|x| -x);
    let mut mapping = mounts(CANTED);
    mapping[0].handedness = PropellerHandedness::Negative;
    let base_a = ShaftPropulsion::new_spatial(&design, &a, mounts(CANTED)).unwrap();
    let base_b = ShaftPropulsion::new_spatial(&design, &b, mapping).unwrap();
    let model_a = DynamicShaftPropulsion::new_oblique(&base_a, RHO, 16).unwrap();
    let model_b = DynamicShaftPropulsion::new_oblique(&base_b, RHO, 16).unwrap();
    let first = initial(&a, [32.0, -39.0]);
    let mut reversed = first;
    reversed.axial_momentum[0] = -reversed.axial_momentum[0];
    reversed.relative_phase_rad[0] = -reversed.relative_phase_rad[0];
    let mut law = properties();
    law.rest_phase_rad[0] = -law.rest_phase_rad[0];
    let mut motion_a = a.shaft_motion(&first, properties(), 0.004).unwrap();
    let mut motion_b = b.shaft_motion(&reversed, law, 0.004).unwrap();
    let mut wake_a = [0.2, 0.8];
    let mut wake_b = wake_a;
    for tick in 0..8 {
        let time = f64::from(tick) * 0.004;
        let result_a = model_a
            .step(
                &mut motion_a,
                &mut wake_a,
                time,
                ZERO,
                |t, state| Ok(control(&a, t, state)),
                wind,
            )
            .unwrap();
        let result_b = model_b
            .step(
                &mut motion_b,
                &mut wake_b,
                time,
                ZERO,
                |t, state| {
                    let mut physical = *state;
                    physical.motion.axial_momentum[0] = -physical.motion.axial_momentum[0];
                    physical.motion.relative_phase_rad[0] = -physical.motion.relative_phase_rad[0];
                    let mut value = control(&a, t, &physical);
                    value.motor_nm[0] = -value.motor_nm[0];
                    Ok(value)
                },
                wind,
            )
            .unwrap();
        assert_eq!(motion_a.state().rigid, motion_b.state().rigid);
        assert_eq!(wake_a, wake_b);
        assert_eq!(result_a.inflow, result_b.inflow);
        assert_eq!(result_a.midpoint.coupled, result_b.midpoint.coupled);
        assert_eq!(result_a.motion.motor_work_j, result_b.motion.motor_work_j);
        assert_eq!(
            result_a.motion.native.dissipated_energy_j,
            result_b.motion.native.dissipated_energy_j
        );
        assert_eq!(
            motion_a.state().axial_momentum[0],
            -motion_b.state().axial_momentum[0]
        );
        assert_eq!(
            motion_a.state().relative_phase_rad[0],
            -motion_b.state().relative_phase_rad[0]
        );
    }
}

#[test]
fn g0_midpoint_density_and_final_native_refusals_preserve_all_states_and_retry() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, CANTED);
    let base = ShaftPropulsion::new_spatial(&design, &a, mounts(CANTED)).unwrap();
    let model = DynamicShaftPropulsion::new(&base, RHO).unwrap();
    let start = initial(&a, [32.0, -39.0]);
    let mut other = a;
    other.rotors[0].axial_inertia_kgm2 *= 2.0;
    let mut incompatible = other.shaft_motion(&start, properties(), 0.004).unwrap();
    let mut incompatible_wake = [0.2, 0.8];
    let mut calls = 0;
    assert_eq!(
        model
            .step(
                &mut incompatible,
                &mut incompatible_wake,
                0.0,
                ZERO,
                |t, state| {
                    calls += 1;
                    Ok(control(&a, t, state))
                },
                wind,
            )
            .unwrap_err()
            .code,
        "propulsion-assembly-invalid"
    );
    assert_eq!(calls, 0);
    assert_eq!(incompatible_wake, [0.2, 0.8]);
    let mut motion = a.shaft_motion(&start, properties(), 0.004).unwrap();
    let before = *motion.state();
    let canonical = *motion.native_state();
    let mut induced = [0.2, 0.8];
    for mechanical_failure in [false, true] {
        let mut calls = 0;
        let error = model
            .step(
                &mut motion,
                &mut induced,
                0.0,
                ZERO,
                |time, state| {
                    calls += 1;
                    let mut values = control(&a, time, state);
                    if time > 0.0 {
                        if mechanical_failure {
                            values.motor_nm[0] = f64::MAX / 32.0;
                        } else {
                            values.density_kg_m3 *= 1.01;
                        }
                    }
                    Ok(values)
                },
                wind,
            )
            .unwrap_err();
        assert_eq!(calls, 2);
        assert_eq!(
            error.code,
            if mechanical_failure {
                "rotor-shaft-refused"
            } else {
                "dynamic-propulsion-invalid"
            }
        );
        assert_eq!(motion.state(), &before);
        assert_eq!(motion.native_state(), &canonical);
        assert_eq!(induced, [0.2, 0.8]);
    }
    let mut fresh = a.shaft_motion(&start, properties(), 0.004).unwrap();
    let mut fresh_wake = induced;
    let retry = model
        .step(
            &mut motion,
            &mut induced,
            0.0,
            ZERO,
            |t, state| Ok(control(&a, t, state)),
            wind,
        )
        .unwrap();
    let expected = model
        .step(
            &mut fresh,
            &mut fresh_wake,
            0.0,
            ZERO,
            |t, state| Ok(control(&a, t, state)),
            wind,
        )
        .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(motion.native_state(), fresh.native_state());
    assert_eq!(induced, fresh_wake);
}

#[test]
fn g0_late_full_wake_refusal_cannot_publish_native_mechanics() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new(&design, &a, mounts(AXIAL)).unwrap();
    let model = DynamicShaftPropulsion::new(&base, RHO).unwrap();
    let start = initial(&a, [35.0, -35.0]);
    let dt = 0.1;
    let mut motion = a.shaft_motion(&start, properties(), dt).unwrap();
    let before = *motion.state();
    let canonical = *motion.native_state();
    let before_wake = [30.0, 0.8];
    let mut induced = before_wake;
    let mut stages = Vec::new();
    let error = model
        .step(
            &mut motion,
            &mut induced,
            0.0,
            ZERO,
            |time, state| {
                stages.push((time, *state));
                Ok(control(&a, time, state))
            },
            |_, _| Ok(ZERO),
        )
        .unwrap_err();
    println!(
        "late wake refusal: {}, completed aerodynamic stages={}",
        error.code,
        stages.len()
    );
    assert_eq!(stages.len(), 2);
    assert_eq!(error.code, "dynamic-inflow-wake-domain");
    let midpoint = model
        .evaluate(&stages[1].1, control(&a, 0.5 * dt, &stages[1].1), &|_| {
            Ok(ZERO)
        })
        .unwrap();
    assert_eq!(
        model
            .inflow_model()
            .step_held(
                before_wake[0],
                midpoint.coupled.disk_axial_mps[0],
                midpoint.coupled.propellers[0].thrust_n,
                dt
            )
            .unwrap_err()
            .code,
        "dynamic-inflow-wake-domain"
    );
    assert_eq!(motion.state(), &before);
    assert_eq!(motion.native_state(), &canonical);
    assert_eq!(induced, before_wake);
}
