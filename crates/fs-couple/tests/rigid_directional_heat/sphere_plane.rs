//! G1/G3/G4: original material claims through moving poses and normal events.
//! This point-contact trajectory reports friction heat and unallocated impact
//! loss. It supplies neither a contact area nor a moving thermal receiver.
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::SourceThermalEllipticState;
use fs_contact::interface_binding::thermal_friction::elliptic::driven::{
    SourceDrivenError, SourceDrivenReport, SourceDrivenSegmentError,
};
use fs_contact::interface_binding::thermal_friction::elliptic::motion::MotionQueryError;
use fs_mbd::sphere_plane::directional::*;
use fs_mbd::sphere_plane::{
    SpherePlane, SpherePlaneConfig, SpherePlaneError, SpherePlaneEventKind, SpherePlaneMode,
    SpherePlaneState,
};
use fs_tribo::tangential::elliptic::driven::{
    EllipticDrivenControls, EllipticDrivenError, EllipticDrivenMode, EllipticDrivenSegment,
    EllipticDrivenSegmentError,
};
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;
use std::cell::Cell;

#[path = "sphere_plane/impacts.rs"]
mod impacts;

/// Pressure-free source: an ideal rigid point contact does not define an area.
pub(super) fn source_fixture(
    kinetic_x: PropertyValue,
) -> (
    InterfaceSystemCard,
    SourceThermalEllipticCoulomb,
    QueryPoint,
) {
    // The bulk sphere fixture agrees with m=2, R=.5 and I=.4*m*R².
    // No caloric curve or thermal receiver is needed for this mechanical path.
    let cards = [0, 1].map(|side| {
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new("density", Density::DIMS),
                value: PropertyValue::Scalar {
                    value: if side == 0 {
                        12.0 / core::f64::consts::PI
                    } else {
                        3.0
                    },
                    dims: Density::DIMS,
                },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("synthetic-moving-contact-{side}"),
                phase: "solid".into(),
                process: "analytical fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap()
    });
    let q = QueryPoint::new()
        .with_quantity("Ti", absolute_temperature(), 300.0)
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.0)
        .unwrap()
        .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.25)
        .unwrap();
    let mut claims = ClaimSet::new();
    for (name, value) in [
        ("static-x", scalar(0.9)),
        ("static-y", scalar(0.9)),
        ("kinetic-x", kinetic_x),
        ("kinetic-y", scalar(0.4)),
    ] {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity("Ti", absolute_temperature(), 200.0, 400.0)
                    .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.0, 12.0)
                    .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.0, 30.0),
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "sphere-surface".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "fixed-plane-texture".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic".into(),
            history: "original".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let axes = vec![
        FrictionQueryAxis {
            name: "speed".into(),
            kind: FrictionQueryAxisKind::SlidingSpeed,
        },
        FrictionQueryAxis {
            name: "path".into(),
            kind: FrictionQueryAxisKind::AccumulatedDistance,
        },
    ];
    let sources = [0, 1].map(|i| {
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
            q.clone(),
            ScalarPropertyRequirement::try_new(
                ["static-x", "static-y"][i],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                ["kinetic-x", "kinetic-y"][i],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "Ti".into(),
                nominal_pressure: None,
            },
            axes.clone(),
        )
        .unwrap()
    });
    let law = SourceThermalEllipticCoulomb::from_principal_sources(sources, [1.0, 0.0]).unwrap();
    (card, law, q)
}

fn sphere(force: Vec3) -> SpherePlane {
    SpherePlane {
        properties: MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.2, 0.2, 0.2)).unwrap(),
        radius: 0.5,
        point_world: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        velocity_world: Vec3::ZERO,
        force_world: force,
        restitution: 0.0,
        friction_coefficient: 0.0,
    }
}
fn resting(speed: f64) -> SpherePlaneState {
    SpherePlaneState::new(
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, 0.0, 0.5), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(2.0 * speed, 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap(),
        0.0,
    )
    .unwrap()
}
fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}
fn settings(phase_s: f64) -> SpherePlaneDirectionalControls {
    SpherePlaneDirectionalControls {
        motion: SpherePlaneConfig::default(),
        maximum_support_phase_s: phase_s,
        maximum_support_phases: 256,
        driven: EllipticDrivenControls {
            maximum_step_s: 0.002,
            minimum_step_s: 1e-10,
            maximum_substeps: 100_000,
            rest_velocity_tolerance_m_s: 1e-10,
            solve: TangentialSolveConfig {
                max_iterations: 128,
                maximum_force_n: 1e5,
                force_absolute_tolerance_n: 1e-11,
                force_relative_tolerance: 1e-13,
            },
            relaxation: EllipticRelaxationControls {
                time_tolerance_s: 1e-11,
                travel_tolerance_m: 1e-11,
                maximum_evaluations: 40_000,
                maximum_root_iterations: 96,
            },
        },
        velocity_tolerance_m_s: 1e-9,
        work_tolerance_j: 1e-9,
    }
}
type SourceError = SourceDrivenSegmentError<SpherePlaneError>;

fn sourced(
    law: &SourceThermalEllipticCoulomb,
    history: &SourceThermalEllipticState,
    phase: &SpherePlaneDirectionalPhase,
    query: &QueryPoint,
    controls: SpherePlaneDirectionalControls,
    visit: &mut dyn FnMut(&EllipticDrivenSegment) -> Result<(), SpherePlaneError>,
    cancel: &mut dyn FnMut() -> bool,
) -> Result<SourceDrivenReport, SourceError> {
    assert_eq!(phase.input.contact_temperature_k, None);
    assert_eq!(phase.input.nominal_pressure_pa, None);
    let mut input = phase.input;
    input.contact_temperature_k = Some(query.axes()["Ti"]);
    law.advance_driven_motion_consistent_with_segments(
        history,
        &input,
        query,
        MotionQuerySolve {
            distance_bracket_m: [0.0, 12.0 * input.dt_s],
            absolute_distance_tolerance_m: 1e-12,
            relative_distance_tolerance: 1e-12,
            maximum_evaluations: 128,
        },
        controls.driven,
        &mut || cancel(),
        &mut |segment| visit(segment),
    )
}

pub(super) fn original_claims(
    card: &InterfaceSystemCard,
    original: &SourceThermalEllipticState,
    current: &SourceThermalEllipticState,
) {
    for (axis, source) in current.friction_sources().iter().enumerate() {
        assert!(!source.query().axes().contains_key("pressure"));
        close(source.state_coordinate("Ti").unwrap(), 300.0);
        for property in source.properties() {
            let old = original.friction_sources()[axis]
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                old.answer().receipt.selected
            );
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}
fn physical_slip(body: RigidBodyState, properties: MassProperties) -> f64 {
    let arm = body
        .pose()
        .orientation()
        .rotate_world_to_body(Vec3::new(0.0, 0.0, -0.5));
    body.point_kinematics(properties, arm)
        .unwrap()
        .point_velocity_world
        .x
}
fn balance(scene: &SpherePlane, before: SpherePlaneState, report: &SpherePlaneDirectionalAdvance) {
    let after = report.motion.state_after.body();
    close(
        after.kinetic_energy(scene.properties).unwrap()
            - before.body().kinetic_energy(scene.properties).unwrap()
            + report.sustained_friction_heat_j
            + report.impact_loss_j,
        report.motion.external_force_work + report.motion.plane_work,
    );
    close(
        report.motion.dissipated_energy,
        report.sustained_friction_heat_j + report.impact_loss_j,
    );
    close(report.motion.energy_balance_residual, 0.0);
    for r in [
        report.motion.linear_momentum_residual,
        report.motion.angular_momentum_residual,
    ] {
        for x in [r.x, r.y, r.z] {
            close(x, 0.0);
        }
    }
}

#[test]
fn sourced_speed_feedback_moves_actual_poses_with_a_three_rung_continuum_ladder() {
    let (card, law, q) = source_fixture(curve("speed", SPEED, [(0.0, 0.2), (12.0, 0.8)]));
    let scene = sphere(Vec3::new(80.0, 0.0, -20.0));
    let before = resting(0.0);
    let original = law.initial_state();
    let duration = 0.5;
    // Actual slip s satisfies s'=40-W*20*(.2+.05s)=33-1.75s.
    // COM motion follows x=(100/7)t²+(2/7) integral(s), including spin recoil.
    let rate: f64 = 1.75;
    let terminal = 33.0 / rate;
    let decay = (-rate * duration).exp();
    let slip = terminal * (1.0 - decay);
    let distance = terminal * (duration - (1.0 - decay) / rate);
    let position = (100.0 / 7.0) * duration * duration + (2.0 / 7.0) * distance;
    let velocity = (200.0 / 7.0) * duration + (2.0 / 7.0) * slip;
    let squared_slip = terminal
        * terminal
        * (duration - 2.0 * (1.0 - decay) / rate + (1.0 - decay * decay) / (2.0 * rate));
    let heat = 4.0 * distance + squared_slip;
    let work = 80.0 * position;
    let mut errors = Vec::new();
    for phases in [4, 8, 16] {
        let controls = settings(duration / phases as f64);
        let mut history = original.clone();
        let mut calls = 0;
        let mut source_distance = 0.0;
        let report = advance_sphere_plane_directional_with(
            &scene,
            before,
            &mut history,
            frame(),
            duration,
            controls,
            |history, phase, visit, cancel| {
                calls += 1;
                close(phase.input.mobility_per_kg[0][0], 1.75);
                close(phase.input.mobility_per_kg[1][1], 1.75);
                assert_eq!(phase.input.mobility_per_kg[0][1], 0.0);
                close(phase.input.normal_n, 20.0);
                close(phase.input.driving_acceleration_m_s2[0], 40.0);
                let response = sourced(&law, history, phase, &q, controls, visit, cancel)?;
                original_claims(&card, &original, &response.state);
                assert!(response.distance_residual_m <= 0.0);
                assert!(response.distance_residual_m.abs() < 2e-11);
                source_distance += response.actual_distance_m;
                for source in response.state.friction_sources() {
                    close(
                        source.state_coordinate("speed").unwrap(),
                        response.actual_distance_m / phase.input.dt_s,
                    );
                    close(
                        source.state_coordinate("path").unwrap(),
                        0.25 + source_distance,
                    );
                }
                close(
                    response.state.constitutive_law().kinetic_coefficients()[0],
                    0.2 + 0.05 * response.actual_distance_m / phase.input.dt_s,
                );
                Ok::<_, SourceError>((response.driven, response.state))
            },
            || false,
        )
        .unwrap();
        assert_eq!(calls, phases);
        assert_eq!(report.support_phases, phases);
        assert_eq!(report.motion.state_after.mode(), SpherePlaneMode::Sliding);
        let after = report.motion.state_after.body();
        let x = after.pose().position_world().x;
        let v = after
            .center_of_mass_velocity_world(scene.properties)
            .unwrap()
            .x;
        let actual_slip = physical_slip(after, scene.properties);
        close(report.support_distance_m, source_distance);
        close(report.support_slip_m[0], source_distance);
        close(report.support_slip_m[1], 0.0);
        close(history.constitutive_state().slip_m()[0], source_distance);
        close(report.motion.external_force_work, 80.0 * x);
        close(report.motion.contact_impulse_world.z, 20.0 * duration);
        close(report.motion.gap_after, 0.0);
        assert_eq!(report.impact_loss_j, 0.0);
        balance(&scene, before, &report);
        original_claims(&card, &original, &history);
        errors.push([
            (actual_slip - slip).abs(),
            (x - position).abs(),
            (v - velocity).abs(),
            (source_distance - distance).abs(),
            (report.sustained_friction_heat_j - heat).abs(),
            (report.motion.external_force_work - work).abs(),
        ]);
    }
    for coordinate in 0..6 {
        assert!(
            errors[2][coordinate] < 0.05,
            "source/pose refinement {errors:?}"
        );
        for pair in errors.windows(2) {
            let ratio = pair[0][coordinate] / pair[1][coordinate];
            assert!(
                ratio > 3.5 && ratio < 4.5,
                "mean-speed source continuum ladder {errors:?}"
            );
        }
    }
}

#[test]
fn original_sources_survive_liftoff_and_recontact_while_normal_reaction_changes_friction() {
    let (card, law, q) = source_fixture(scalar(0.2));
    let controls = settings(0.03125);
    let scene = sphere(Vec3::new(0.0, 0.0, -20.0));
    let before = resting(1.0);
    let original = law.initial_state();
    let mut history = original.clone();
    let supported = advance_sphere_plane_directional_with(
        &scene,
        before,
        &mut history,
        frame(),
        0.0625,
        controls,
        |history, phase, visit, cancel| {
            close(phase.input.normal_n, 20.0);
            let response = sourced(
                &law,
                history,
                phase,
                &q,
                controls,
                &mut |segment| {
                    close(segment.traction_n[0], 4.0);
                    visit(segment)
                },
                cancel,
            )?;
            Ok::<_, SourceError>((response.driven, response.state))
        },
        || false,
    )
    .unwrap();
    close(supported.support_distance_m, 0.0625 - 3.5 * 0.0625 * 0.0625);
    close(
        supported.sustained_friction_heat_j,
        4.0 * supported.support_distance_m,
    );
    let slip_before_flight = physical_slip(supported.motion.state_after.body(), scene.properties);
    close(slip_before_flight, 0.5625);
    balance(&scene, before, &supported);
    let accepted = history.clone();
    let airborne_scene = sphere(Vec3::new(0.0, 0.0, 20.0));
    let mut queries_in_flight = 0;
    let flight = advance_sphere_plane_directional_with(
        &airborne_scene,
        supported.motion.state_after,
        &mut history,
        frame(),
        0.125,
        controls,
        |history, phase, visit, cancel| {
            queries_in_flight += 1;
            let response = sourced(&law, history, phase, &q, controls, visit, cancel)?;
            Ok::<_, SourceError>((response.driven, response.state))
        },
        || false,
    )
    .unwrap();
    assert_eq!(queries_in_flight, 0);
    assert_eq!(history, accepted);
    assert_eq!(flight.motion.state_after.mode(), SpherePlaneMode::Flight);
    assert_eq!(flight.support_phases, 0);
    assert_eq!(flight.sustained_friction_heat_j, 0.0);
    assert_eq!(flight.impact_loss_j, 0.0);
    close(flight.motion.gap_after, 0.078125);
    close(
        flight
            .motion
            .state_after
            .body()
            .center_of_mass_velocity_world(scene.properties)
            .unwrap()
            .z,
        1.25,
    );
    let lift = flight
        .motion
        .events
        .iter()
        .find(|e| e.kind == SpherePlaneEventKind::LiftOff)
        .unwrap();
    close(lift.time_s, 0.0625);
    balance(&airborne_scene, supported.motion.state_after, &flight);

    let return_scene = sphere(Vec3::new(0.0, 0.0, -40.0));
    let impact_speed = 4.6875_f64.sqrt();
    let flight_time = (1.25 + impact_speed) / 20.0;
    let mut loaded_segments = 0;
    let returned = advance_sphere_plane_directional_with(
        &return_scene,
        flight.motion.state_after,
        &mut history,
        frame(),
        0.25,
        controls,
        |history, phase, visit, cancel| {
            close(phase.input.normal_n, 40.0);
            let response = sourced(
                &law,
                history,
                phase,
                &q,
                controls,
                &mut |segment| {
                    if segment.mode_before == EllipticDrivenMode::Sliding {
                        close(segment.traction_n[0], 8.0);
                        loaded_segments += 1;
                    }
                    visit(segment)
                },
                cancel,
            )?;
            original_claims(&card, &original, &response.state);
            Ok::<_, SourceError>((response.driven, response.state))
        },
        || false,
    )
    .unwrap();
    assert!(loaded_segments > 0);
    let stop_time = slip_before_flight / 14.0;
    let stopped_distance = 0.5 * slip_before_flight * stop_time;
    close(returned.support_distance_m, stopped_distance);
    close(returned.sustained_friction_heat_j, 8.0 * stopped_distance);
    close(returned.impact_loss_j, 4.6875);
    assert!(returned.impact_loss_j > 40.0 * returned.sustained_friction_heat_j);
    let impact = returned
        .motion
        .events
        .iter()
        .find(|e| e.kind == SpherePlaneEventKind::Impact)
        .unwrap();
    close(impact.time_s, 0.1875 + flight_time);
    close(impact.normal_impulse, 2.0 * impact_speed);
    let arrest = returned
        .motion
        .events
        .iter()
        .find(|e| e.kind == SpherePlaneEventKind::Rolling)
        .unwrap();
    close(arrest.time_s, impact.time_s + stop_time);
    close(
        returned.motion.contact_impulse_world.z,
        2.0 * impact_speed + 40.0 * (0.25 - flight_time),
    );
    assert_eq!(returned.motion.state_after.mode(), SpherePlaneMode::Rolling);
    close(
        physical_slip(returned.motion.state_after.body(), scene.properties),
        0.0,
    );
    close(returned.motion.gap_after, 0.0);
    for source in history.friction_sources() {
        close(
            source.state_coordinate("path").unwrap(),
            0.25 + supported.support_distance_m + stopped_distance,
        );
    }
    balance(&return_scene, flight.motion.state_after, &returned);
    original_claims(&card, &original, &history);
}

#[test]
fn late_source_failure_and_cancellation_discard_private_poses_and_history_then_retry() {
    let (card, law, q) = source_fixture(curve("speed", SPEED, [(0.0, 0.2), (12.0, 0.8)]));
    let scene = sphere(Vec3::new(80.0, 0.0, -20.0));
    let input = resting(0.0);
    let original_body = input.body();
    let original_history = law.initial_state();
    let controls = settings(0.03125);
    for cancel_late in [false, true] {
        let mut history = original_history.clone();
        let mut calls = 0;
        let mut visits = 0;
        let cancelled = Cell::new(false);
        let mut private_position = 0.0_f64;
        let result = advance_sphere_plane_directional_with(
            &scene,
            input,
            &mut history,
            frame(),
            0.125,
            controls,
            |history, phase, visit, cancel| {
                calls += 1;
                private_position = private_position.max(phase.body.pose().position_world().x);
                let query = if !cancel_late && calls == 3 {
                    q.clone().with("Ti", 450.0).unwrap()
                } else {
                    q.clone()
                };
                let response = sourced(
                    &law,
                    history,
                    phase,
                    &query,
                    controls,
                    &mut |segment| {
                        visit(segment)?;
                        visits += 1;
                        if cancel_late && visits == 3 {
                            cancelled.set(true);
                        }
                        Ok(())
                    },
                    cancel,
                )?;
                Ok::<_, SourceError>((response.driven, response.state))
            },
            || cancelled.get(),
        );
        if cancel_late {
            assert!(matches!(
                result,
                Err(SpherePlaneDirectionalError::Source(
                    SourceDrivenSegmentError::Segments(EllipticDrivenSegmentError::Driven(
                        EllipticDrivenError::Cancelled
                    ))
                ))
            ));
            assert_eq!(visits, 3);
        } else {
            assert!(matches!(
                result,
                Err(SpherePlaneDirectionalError::Source(
                    SourceDrivenSegmentError::Source(SourceDrivenError::Query(
                        MotionQueryError::Source(_)
                    ))
                ))
            ));
            assert_eq!(visits, 2);
        }
        assert_eq!(calls, 3);
        assert!(private_position > 0.0);
        assert_eq!(history, original_history);
        assert_eq!(input.body(), original_body);
        assert_eq!(input.time_s(), 0.0);
    }
    let mut accepted = None;
    for _ in 0..2 {
        let mut history = original_history.clone();
        let report = advance_sphere_plane_directional_with(
            &scene,
            input,
            &mut history,
            frame(),
            0.125,
            controls,
            |history, phase, visit, cancel| {
                let response = sourced(&law, history, phase, &q, controls, visit, cancel)?;
                Ok::<_, SourceError>((response.driven, response.state))
            },
            || false,
        )
        .unwrap();
        assert_eq!(report.support_phases, 4);
        assert!(report.motion.state_after.body().pose().position_world().x > 0.0);
        assert_ne!(history, original_history);
        balance(&scene, input, &report);
        original_claims(&card, &original_history, &history);
        let values = (
            report.motion.state_after,
            history,
            report.support_distance_m,
            report.sustained_friction_heat_j,
        );
        if let Some(previous) = &accepted {
            assert_eq!(&values, previous);
        } else {
            accepted = Some(values);
        }
    }
}
