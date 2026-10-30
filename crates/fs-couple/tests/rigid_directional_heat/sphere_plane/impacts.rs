//! Original sourced impact, supported arrest and static capture in one owner.
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::SourceFrictionConditions;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_couple::friction_heat::directional::sphere_plane::*;
use fs_tribo::tangential::elliptic::impulse::EllipticImpulseControls;

fn conditions(query: &QueryPoint) -> SourceFrictionConditions {
    SourceFrictionConditions {
        query: query.clone(),
        contact_temperature_k: Some(query.axes()["Ti"]),
        nominal_pressure_pa: None,
    }
}

fn support_conditions(
    query: &QueryPoint,
    phase: &SpherePlaneDirectionalPhase,
) -> SourceSpherePlaneSupportConditions {
    SourceSpherePlaneSupportConditions {
        friction: conditions(query),
        motion: MotionQuerySolve {
            distance_bracket_m: [0.0, 12.0 * phase.input.dt_s],
            absolute_distance_tolerance_m: 1e-12,
            relative_distance_tolerance: 1e-12,
            maximum_evaluations: 128,
        },
    }
}

fn controls() -> SourceSpherePlaneControls {
    SourceSpherePlaneControls {
        motion: settings(0.015625),
        impulse: EllipticImpulseControls {
            maximum_iterations: 128,
            maximum_impulse_ns: 1e5,
            impulse_absolute_tolerance_ns: 1e-11,
            impulse_relative_tolerance: 1e-13,
            velocity_tolerance_m_s: 1e-9,
            energy_tolerance_j: 1e-9,
        },
        impact_coefficient_role: SourceImpulseCoefficientRole::Kinetic,
    }
}

fn incoming(height: f64, speed: f64, normal_speed: f64) -> SpherePlaneState {
    SpherePlaneState::new(
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, height), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(2.0 * speed, 0.0, 2.0 * normal_speed),
            Vec3::ZERO,
        )
        .unwrap(),
        0.0,
    )
    .unwrap()
}

#[test]
fn sourced_preimpact_speed_drives_actual_impulse_spin_and_bounce_without_travel() {
    let (card, law, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.7)]));
    let mut scene = sphere(Vec3::ZERO);
    scene.restitution = 0.5;
    scene.friction_coefficient = 0.8;
    let before = incoming(0.5, 1.0, -1.0);
    let original = law.initial_state();
    let mut history = original.clone();
    let report = advance_source_sphere_plane(
        &scene,
        before,
        &law,
        &mut history,
        frame(),
        0.125,
        controls(),
        |_, phase| Ok::<_, &'static str>(support_conditions(&query, phase)),
        |_, phase| {
            close(phase.time_s, 0.0);
            close(phase.input.normal_impulse_ns, 3.0);
            Ok::<_, &'static str>(conditions(&query))
        },
        || false,
    )
    .unwrap();
    let motion = &report.trajectory;
    assert_eq!(motion.support_phases, 0);
    assert_eq!(report.support_source_evaluations, 0);
    assert_eq!(report.impacts.len(), 1);
    assert_eq!(motion.motion.state_after.mode(), SpherePlaneMode::Flight);
    close(motion.support_distance_m, 0.0);
    close(motion.sustained_friction_heat_j, 0.0);
    let receipt = &report.impacts[0];
    close(receipt.response.queried_speed_m_s, 1.0);
    close(receipt.response.coefficients_mu[0], 0.15);
    close(receipt.response.impulse.impulse_ns[0], -0.45);
    close(receipt.response.impulse.velocity_after_m_s[0], 0.2125);
    close(receipt.response.impulse.dissipated_j, 0.2728125);
    close(motion.impact_loss_j, 0.75 + 0.2728125);
    close(motion.motion.contact_impulse_world.x, -0.45);
    close(motion.motion.contact_impulse_world.z, 3.0);
    let after = motion.motion.state_after.body();
    close(after.pose().position_world().x, 0.775 * 0.125);
    close(after.pose().position_world().z, 0.5 + 0.5 * 0.125);
    close(
        after
            .center_of_mass_velocity_world(scene.properties)
            .unwrap()
            .x,
        0.775,
    );
    close(
        after
            .point_kinematics(scene.properties, Vec3::ZERO)
            .unwrap()
            .angular_velocity_world
            .y,
        1.125,
    );
    assert_eq!(history.constitutive_state().slip_m(), [0.0; 2]);
    close(history.constitutive_law().static_coefficients()[0], 0.9);
    close(history.constitutive_law().kinetic_coefficients()[0], 0.15);
    for source in history.friction_sources() {
        close(source.state_coordinate("path").unwrap(), 0.25);
        close(source.state_coordinate("speed").unwrap(), 1.0);
    }
    original_claims(&card, &original, &receipt.response.state);
    original_claims(&card, &original, &history);
    balance(&scene, before, motion);
}

#[test]
fn impact_then_sourced_arrest_retains_receipts_and_full_static_hold_capacity() {
    let (card, law, query) = source_fixture(scalar(0.2));
    let scene = sphere(Vec3::new(0.0, 0.0, -20.0));
    let before = incoming(0.5, 1.0, -1.0);
    let original = law.initial_state();
    let mut history = original.clone();
    let report = advance_source_sphere_plane(
        &scene,
        before,
        &law,
        &mut history,
        frame(),
        0.125,
        controls(),
        |_, phase| Ok::<_, &'static str>(support_conditions(&query, phase)),
        |_, _| Ok::<_, &'static str>(conditions(&query)),
        || false,
    )
    .unwrap();
    let motion = &report.trajectory;
    assert_eq!(report.impacts.len(), 1);
    assert!(report.support_source_evaluations > 0);
    assert_eq!(motion.motion.state_after.mode(), SpherePlaneMode::Rolling);
    // Jn=2, pt=-.4, v_after=1-1.75*.4=.3; support decelerates at7 m/s².
    let distance = 0.3 * 0.3 / 14.0;
    close(motion.impact_loss_j, 1.0 + 0.26);
    close(motion.support_distance_m, distance);
    close(motion.sustained_friction_heat_j, 4.0 * distance);
    close(history.constitutive_state().slip_m()[0], distance);
    let receipt = &report.impacts[0];
    close(receipt.response.impulse.impulse_ns[0], -0.4);
    for source in receipt.response.state.friction_sources() {
        close(source.state_coordinate("path").unwrap(), 0.25);
        close(source.state_coordinate("speed").unwrap(), 1.0);
    }
    for source in history.friction_sources() {
        close(source.state_coordinate("path").unwrap(), 0.25 + distance);
    }
    close(history.constitutive_law().static_coefficients()[0], 0.9);
    close(history.constitutive_law().kinetic_coefficients()[0], 0.2);
    original_claims(&card, &original, &receipt.response.state);
    original_claims(&card, &original, &history);
    balance(&scene, before, motion);

    // Holding needs40/7 N, above the4 N kinetic capacity and below18 N static.
    let driven_scene = sphere(Vec3::new(20.0, 0.0, -20.0));
    let supported = motion.motion.state_after;
    let held = advance_source_sphere_plane(
        &driven_scene,
        supported,
        &law,
        &mut history,
        frame(),
        0.125,
        controls(),
        |_, phase| Ok::<_, &'static str>(support_conditions(&query, phase)),
        |_, _| Ok::<_, &'static str>(conditions(&query)),
        || false,
    )
    .unwrap();
    assert!(held.impacts.is_empty());
    assert_eq!(
        held.trajectory.motion.state_after.mode(),
        SpherePlaneMode::Rolling
    );
    close(held.trajectory.support_distance_m, 0.0);
    close(held.trajectory.sustained_friction_heat_j, 0.0);
    close(
        physical_slip(held.trajectory.motion.state_after.body(), scene.properties),
        0.0,
    );
    for source in history.friction_sources() {
        close(source.state_coordinate("path").unwrap(), 0.25 + distance);
    }
    balance(&driven_scene, supported, &held.trajectory);
}

#[test]
fn source_refusal_or_cancellation_after_impact_rolls_back_the_whole_trajectory() {
    let (_, law, query) = source_fixture(scalar(0.2));
    let scene = sphere(Vec3::new(0.0, 0.0, -20.0));
    let before = incoming(0.5, 1.0, -1.0);
    let original = law.initial_state();
    for cancel_late in [false, true] {
        let mut history = original.clone();
        let impacts = Cell::new(0);
        let supports = Cell::new(0);
        let stop = Cell::new(false);
        let result = advance_source_sphere_plane(
            &scene,
            before,
            &law,
            &mut history,
            frame(),
            0.125,
            controls(),
            |accepted, phase| {
                assert_eq!(impacts.get(), 1);
                supports.set(supports.get() + 1);
                if supports.get() == 2 {
                    assert!(accepted.constitutive_state().slip_m()[0] > 0.0);
                    if cancel_late {
                        stop.set(true);
                    } else {
                        let invalid = query.clone().with("Ti", 450.0).unwrap();
                        return Ok::<_, &'static str>(support_conditions(&invalid, phase));
                    }
                }
                Ok::<_, &'static str>(support_conditions(&query, phase))
            },
            |_, _| {
                impacts.set(impacts.get() + 1);
                Ok::<_, &'static str>(conditions(&query))
            },
            || stop.get(),
        );
        assert!(matches!(
            result,
            Err(SpherePlaneDirectionalError::Source(
                SourceSpherePlaneError::Support(_)
            ))
        ));
        assert_eq!(impacts.get(), 1);
        assert_eq!(supports.get(), 2);
        assert_eq!(history, original);
        assert_eq!(before.body(), incoming(0.5, 1.0, -1.0).body());
        assert_eq!(before.time_s(), 0.0);
    }

    let mut accepted = None;
    for _ in 0..2 {
        let mut history = original.clone();
        let report = advance_source_sphere_plane(
            &scene,
            before,
            &law,
            &mut history,
            frame(),
            0.125,
            controls(),
            |_, phase| Ok::<_, &'static str>(support_conditions(&query, phase)),
            |_, _| Ok::<_, &'static str>(conditions(&query)),
            || false,
        )
        .unwrap();
        let values = (
            report.trajectory.motion.state_after,
            history,
            report.impacts[0].response.impulse,
            report.trajectory.sustained_friction_heat_j,
        );
        if let Some(previous) = &accepted {
            assert_eq!(&values, previous);
        } else {
            accepted = Some(values);
        }
    }
}

#[test]
fn failed_collision_query_after_ballistic_travel_preserves_original_body_and_sources() {
    let (_, law, query) = source_fixture(scalar(0.2));
    let mut scene = sphere(Vec3::ZERO);
    scene.restitution = 0.5;
    let before = incoming(1.0, 1.0, -1.0);
    let original = law.initial_state();
    let mut history = original.clone();
    let calls = Cell::new(0);
    let result = advance_source_sphere_plane(
        &scene,
        before,
        &law,
        &mut history,
        frame(),
        1.0,
        controls(),
        |_, phase| Ok::<_, &'static str>(support_conditions(&query, phase)),
        |accepted, phase| {
            calls.set(calls.get() + 1);
            assert_eq!(accepted, &original);
            close(phase.time_s, 0.5);
            close(phase.body.pose().position_world().x, 0.5);
            Err::<SourceFrictionConditions, _>("external interface source unavailable")
        },
        || false,
    );
    assert!(matches!(
        result,
        Err(SpherePlaneDirectionalError::Source(
            SourceSpherePlaneError::ImpactConditions("external interface source unavailable")
        ))
    ));
    assert_eq!(calls.get(), 1);
    assert_eq!(history, original);
    assert_eq!(before.time_s(), 0.0);
    close(before.body().pose().position_world().z, 1.0);
}
