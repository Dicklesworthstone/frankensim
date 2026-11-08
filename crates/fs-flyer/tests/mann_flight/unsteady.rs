//! Actual lift/moment, shaft/wake memory and trajectories under coherent gusts.
use super::*;
use fs_airfoil::unsteady::UNSTEADY_SECTION_V1;
use fs_flyer::Refusal;
use fs_flyer::spine::rotors::geared::GearedState;
use fs_flyer::spine::rotors::propulsion::dynamic::unsteady::DynamicUnsteadyPropulsion;
use fs_flyer::spine::rotors::propulsion::unsteady::UnsteadyPropulsion;
use fs_ga::{Quat, So3, Vec3};
use fs_math::det;

fn v(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn a(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
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

// Independent kinematic construction: the difference is measured at the
// original bound-row midpoint, with the same rigid-point velocity in each air.
fn expected_gust(
    design: &OpenLoopDesign,
    motion: &GearedState,
    control: f64,
    strip: usize,
    total: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
    reference: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
) -> f64 {
    let (panels, strips) = design.layout(control).unwrap();
    let panel = &panels[strips[strip].panel_indices[0]];
    let point = (v(panel.a) + v(panel.b)).scale(0.5);
    let arm = point - v(design.cg_m);
    let r = rotation(motion);
    let world = a(v(motion.rigid.pos_m) + r.rotate(arm).unwrap());
    let body_velocity = |wind| {
        r.inverse()
            .rotate(v(motion.rigid.vel_mps) - v(wind))
            .unwrap()
            + v(motion.rigid.omega_body).cross(arm)
    };
    let full = body_velocity(total(world).unwrap());
    let mean = body_velocity(reference(world).unwrap());
    det::atan2(full.z, full.x) - det::atan2(mean.z, mean.x)
}

#[test]
fn native_reference_partition_matches_actual_material_point_kinematics() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    let assembly = assembly(&design);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts()).unwrap();
    let model = UnsteadyPropulsion::new(&base, UNSTEADY_SECTION_V1).unwrap();
    let mut motion = state(&assembly).motion;
    // Nonzero pose and all three angular rates make CG-only and unrotated
    // gust decompositions observably incorrect.
    let pose = Quat::from_axis_angle(v([0.0, 1.0, 0.0]), 0.04);
    motion.rigid.quat = [pose.w, pose.x, pose.y, pose.z];
    let reference = |point: [f64; 3]| Ok([-1.0 + 0.01 * point[2], 0.02 * point[0], 0.0]);
    let total = |point: [f64; 3]| {
        let mean = reference(point)?;
        Ok([mean[0] + 0.03 * point[1], mean[1], -0.18 - 0.01 * point[1]])
    };
    let zero = vec![0.0; model.section_count()];
    let accepted = model
        .trimmed_state(motion, inputs(), &reference, &zero, None)
        .unwrap();
    let before = accepted.clone();
    let calls = Cell::new(0);
    let resolved = model
        .evaluate_with_reference_wind(
            &accepted,
            inputs(),
            &total,
            &|point| {
                calls.set(calls.get() + 1);
                reference(point)
            },
            None,
        )
        .unwrap();
    assert_eq!(calls.get(), model.section_count());
    let expected = (0..model.section_count())
        .map(|s| expected_gust(&design, &motion, inputs().canard_rad, s, &total, &reference))
        .collect::<Vec<_>>();
    for (actual, expected) in resolved.aerodynamics.gust_angle_rad.iter().zip(&expected) {
        assert!((actual - expected).abs() < 2e-15);
    }
    assert!(
        expected
            .windows(2)
            .any(|pair| (pair[0] - pair[1]).abs() > 1e-6)
    );
    let explicit = model
        .evaluate(&accepted, inputs(), &total, &expected, None)
        .unwrap();
    assert_eq!(resolved, explicit);
    let as_motion = model
        .evaluate(&accepted, inputs(), &total, &zero, None)
        .unwrap();
    assert!(
        (resolved.aerodynamic_force_body_n[2] - as_motion.aerodynamic_force_body_n[2]).abs() > 0.1
    );
    assert_ne!(
        resolved.loads.carrier.moment_nm,
        as_motion.loads.carrier.moment_nm
    );
    assert_eq!(accepted, before);
}

#[test]
fn identical_total_and_reference_preserve_both_original_zero_gust_owners() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    let assembly = assembly(&design);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts()).unwrap();
    let model = UnsteadyPropulsion::new(&base, UNSTEADY_SECTION_V1).unwrap();
    let field = |point: [f64; 3]| Ok([-1.0, 0.0, 0.01 * point[1]]);
    let controls = |time: f64| {
        let mut result = inputs();
        if time > 0.0 {
            result.canard_rad += 0.01;
        }
        result
    };
    let zero = vec![0.0; model.section_count()];
    let accepted = model
        .trimmed_state(state(&assembly).motion, inputs(), &field, &zero, None)
        .unwrap();
    let reference = model
        .trimmed_state_with_reference_wind(state(&assembly).motion, inputs(), &field, &field, None)
        .unwrap();
    assert_eq!(reference, accepted);
    let explicit = model
        .step(
            &accepted,
            0.0,
            0.002,
            ORIGIN,
            None,
            |time, _| Ok(controls(time)),
            |_, point| field(point),
            |_, _, _| Ok(0.0),
        )
        .unwrap();
    let actual = model
        .step_with_reference_wind(
            &accepted,
            0.0,
            0.002,
            ORIGIN,
            None,
            |time, _| Ok(controls(time)),
            |_, point| field(point),
            |_, point| field(point),
        )
        .unwrap();
    assert_eq!(actual, explicit);
    for stage in [&actual.initial, &actual.midpoint] {
        assert!(
            stage
                .aerodynamics
                .gust_angle_rad
                .iter()
                .all(|gust| *gust == 0.0)
        );
    }

    let dynamic = DynamicPropulsion::new(&base, RHO).unwrap();
    let model = DynamicUnsteadyPropulsion::new(&dynamic, UNSTEADY_SECTION_V1).unwrap();
    let accepted = model
        .trimmed_state(state(&assembly), inputs(), &field, &zero)
        .unwrap();
    let reference = model
        .trimmed_state_with_reference_wind(state(&assembly), inputs(), &field, &field)
        .unwrap();
    assert_eq!(reference, accepted);
    let explicit = model
        .step(
            &accepted,
            0.0,
            0.002,
            ORIGIN,
            |time, _| Ok(controls(time)),
            |_, point| field(point),
            |_, _, _| Ok(0.0),
        )
        .unwrap();
    let actual = model
        .step_with_reference_wind(
            &accepted,
            0.0,
            0.002,
            ORIGIN,
            |time, _| Ok(controls(time)),
            |_, point| field(point),
            |_, point| field(point),
        )
        .unwrap();
    assert_eq!(actual, explicit);
    for stage in [&actual.initial, &actual.midpoint] {
        assert!(
            stage
                .aerodynamics
                .gust_angle_rad
                .iter()
                .all(|gust| *gust == 0.0)
        );
    }
}

#[test]
fn zero_and_reversed_reference_air_refuse_without_fabricating_gust() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    let assembly = assembly(&design);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts()).unwrap();
    let model = UnsteadyPropulsion::new(&base, UNSTEADY_SECTION_V1).unwrap();
    let mut motion = state(&assembly).motion;
    motion.rigid.omega_body = [0.0; 3];
    let total = |_| Ok([-1.0, 0.0, 0.0]);
    let zero = vec![0.0; model.section_count()];
    let accepted = model
        .trimmed_state(motion, inputs(), &total, &zero, None)
        .unwrap();
    let before = accepted.clone();
    for forward_wind_increment in [0.0, 1.0] {
        let error = model
            .evaluate_with_reference_wind(
                &accepted,
                inputs(),
                &total,
                &|_| {
                    let mut wind = motion.rigid.vel_mps;
                    wind[0] += forward_wind_increment;
                    Ok(wind)
                },
                None,
            )
            .unwrap_err();
        assert_eq!(error.code, "incident-flow-invalid");
        assert_eq!(accepted, before);
    }
}

#[test]
fn mann_dynamic_unsteady_loads_use_real_half_step_fields_and_replay_all_histories() {
    let design = wright_openloop_v1();
    let assembly = assembly(&design);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts()).unwrap();
    let dynamic = DynamicPropulsion::new_oblique(&base, RHO, 16)
        .unwrap()
        .with_blade_incident(BladeIncident::WindField)
        .unwrap();
    let model = DynamicUnsteadyPropulsion::new(&dynamic, UNSTEADY_SECTION_V1).unwrap();
    let mut atmosphere = wind();
    let accepted = atmosphere
        .trimmed_dynamic_unsteady_state(&model, state(&assembly), inputs())
        .unwrap();
    let before_state = accepted.clone();
    let checkpoint = atmosphere.checkpoint();
    let start = atmosphere.time_s();
    let dt = atmosphere.dt_s();
    let held = atmosphere.clone();
    let mut source = FixedStepMannField::resume(checkpoint.field.clone(), 88 * MODES).unwrap();
    source
        .advance_to(source.step() + 1, WORK / 2, || false)
        .unwrap();
    let mut mid_checkpoint = checkpoint.clone();
    mid_checkpoint.field = source.checkpoint();
    let midpoint = MannFlightWind::resume(mid_checkpoint, BUDGET).unwrap();
    let controls = |time: f64| {
        let mut result = inputs();
        if time != start {
            result.canard_rad += 0.002;
        }
        result
    };
    let direct = model
        .step_with_reference_wind(
            &accepted,
            start,
            dt,
            ORIGIN,
            |time, _| Ok(controls(time)),
            |time, point| {
                if time == start {
                    held.sample(point)
                } else {
                    assert_eq!(time, start + dt * 0.5);
                    midpoint.sample(point)
                }
                .map(|flow| flow.u)
            },
            |_, point| held.sample_mean(point).map(|flow| flow.u),
        )
        .unwrap();
    let actual = atmosphere
        .step_dynamic_unsteady(
            &model,
            &accepted,
            ORIGIN,
            WORK,
            || false,
            |time, _| Ok(controls(time)),
        )
        .unwrap();
    assert_eq!(actual, direct);
    assert_eq!(accepted, before_state);
    source
        .advance_to(source.step() + 1, WORK / 2, || false)
        .unwrap();
    assert_eq!(atmosphere.checkpoint().field, source.checkpoint());
    assert_ne!(actual.state.sections, accepted.sections);
    assert_ne!(
        actual.state.dynamic.induced_mps,
        accepted.dynamic.induced_mps
    );
    assert_ne!(actual.state.dynamic.motion, accepted.dynamic.motion);
    assert!(
        actual
            .initial
            .aerodynamics
            .gust_angle_rad
            .iter()
            .any(|a| a.abs() > 1e-6)
    );
    assert_ne!(
        actual.initial.aerodynamics.gust_angle_rad,
        actual.midpoint.aerodynamics.gust_angle_rad
    );
    let frozen = model
        .step_with_reference_wind(
            &accepted,
            start,
            dt,
            ORIGIN,
            |time, _| Ok(controls(time)),
            |_, point| held.sample(point).map(|flow| flow.u),
            |_, point| held.sample_mean(point).map(|flow| flow.u),
        )
        .unwrap();
    assert_eq!(actual.initial, frozen.initial);
    assert_ne!(
        actual.midpoint.aerodynamic_force_body_n,
        frozen.midpoint.aerodynamic_force_body_n
    );
    assert_ne!(
        actual.midpoint.loads.aerodynamic_nm,
        frozen.midpoint.loads.aerodynamic_nm
    );
    assert_ne!(actual.state, frozen.state);
    let mut replay = MannFlightWind::resume(checkpoint, BUDGET).unwrap();
    let repeated = replay
        .step_dynamic_unsteady(
            &model,
            &accepted,
            ORIGIN,
            WORK,
            || false,
            |time, _| Ok(controls(time)),
        )
        .unwrap();
    assert_eq!(actual, repeated);
    assert_eq!(atmosphere.checkpoint(), replay.checkpoint());
}

#[test]
fn mann_equilibrium_wake_unsteady_loads_share_the_same_reference_and_atomic_clock() {
    let mut design = wright_openloop_v1();
    design.coupling = fs_flyer::propcoupling::CANDIDATE_B;
    let assembly = assembly(&design);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts()).unwrap();
    let model = UnsteadyPropulsion::new(&base, UNSTEADY_SECTION_V1).unwrap();
    let mut atmosphere = wind();
    let accepted = atmosphere
        .trimmed_unsteady_state(&model, state(&assembly).motion, inputs(), None)
        .unwrap();
    let direct = model
        .trimmed_state_with_reference_wind(
            state(&assembly).motion,
            inputs(),
            &|point| atmosphere.sample(point).map(|flow| flow.u),
            &|point| atmosphere.sample_mean(point).map(|flow| flow.u),
            None,
        )
        .unwrap();
    assert_eq!(accepted, direct);
    let before = atmosphere.checkpoint();
    let mut replay = MannFlightWind::resume(before.clone(), BUDGET).unwrap();
    let actual = atmosphere
        .step_unsteady(
            &model,
            &accepted,
            ORIGIN,
            None,
            WORK,
            || false,
            |_, _| Ok(inputs()),
        )
        .unwrap();
    let repeated = replay
        .step_unsteady(
            &model,
            &accepted,
            ORIGIN,
            None,
            WORK,
            || false,
            |_, _| Ok(inputs()),
        )
        .unwrap();
    assert_eq!(actual, repeated);
    assert_eq!(atmosphere.checkpoint(), replay.checkpoint());
    assert_ne!(actual.state.sections, accepted.sections);
    assert_ne!(actual.state.motion, accepted.motion);
    assert_ne!(atmosphere.checkpoint().field, before.field);
}

#[test]
fn late_reference_or_mann_refusal_withholds_wind_wake_and_section_memories() {
    let design = wright_openloop_v1();
    let assembly = assembly(&design);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts()).unwrap();
    let dynamic = DynamicPropulsion::new(&base, RHO).unwrap();
    let model = DynamicUnsteadyPropulsion::new(&dynamic, UNSTEADY_SECTION_V1).unwrap();
    let mut atmosphere = wind();
    let accepted = atmosphere
        .trimmed_dynamic_unsteady_state(&model, state(&assembly), inputs())
        .unwrap();
    let before = accepted.clone();
    let wind_before = atmosphere.checkpoint();
    let start = atmosphere.time_s();
    let stage_calls = Cell::new(0);
    let failed = model.step_with_reference_wind(
        &accepted,
        start,
        atmosphere.dt_s(),
        ORIGIN,
        |_, _| {
            stage_calls.set(stage_calls.get() + 1);
            Ok(inputs())
        },
        |_, point| atmosphere.sample(point).map(|flow| flow.u),
        |time, point| {
            if time == start {
                atmosphere.sample_mean(point).map(|flow| flow.u)
            } else {
                Ok([f64::NAN, 0.0, 0.0])
            }
        },
    );
    assert!(failed.is_err());
    assert_eq!(stage_calls.get(), 2);
    assert_eq!(accepted, before);

    let calls = Cell::new(0);
    assert!(
        atmosphere
            .step_dynamic_unsteady(
                &model,
                &accepted,
                ORIGIN,
                WORK,
                || false,
                |time, _| {
                    calls.set(calls.get() + 1);
                    let mut input = inputs();
                    if time != start {
                        input.density_kg_m3 *= 1.1;
                    }
                    Ok(input)
                }
            )
            .is_err()
    );
    assert_eq!(calls.get(), 2);
    assert_eq!(atmosphere.checkpoint(), wind_before);
    assert_eq!(accepted, before);
    let finished = Cell::new(false);
    let error = atmosphere
        .step_dynamic_unsteady(
            &model,
            &accepted,
            ORIGIN,
            WORK,
            || finished.get(),
            |time, _| {
                if time != start {
                    finished.set(true);
                }
                Ok(inputs())
            },
        )
        .unwrap_err();
    assert_eq!(error.code, "mann-flight-cancelled");
    assert_eq!(atmosphere.checkpoint(), wind_before);
    assert_eq!(accepted, before);
}
