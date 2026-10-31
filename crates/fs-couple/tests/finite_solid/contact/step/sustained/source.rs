//! Native source travel, retry rollback, retained feature identity and texture.
use super::super::super::elliptic::{conditions, source as interface_source};
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::history::SourcePostcontactMotionError;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticCoulomb,
};
use fs_couple::finite_solid::contact::mixed::SourceFiniteSolidMixedContactBinding;
use fs_couple::finite_solid::contact::step::sustained::source::*;
use fs_couple::finite_solid::contact::step::sustained::*;
use fs_couple::finite_solid::contact::texture::FiniteSolidContactTexture;
use std::cell::Cell;

const SOURCE_STEP: f64 = STEP;

fn source_controls(step: f64) -> FiniteSolidSphereTrajectoryControls {
    let mut controls = sustained_controls(step);
    controls.maximum_steps = 512;
    controls
}

fn binding(law: &SourceThermalEllipticCoulomb) -> FiniteSolidSourceContactBinding<'_> {
    FiniteSolidSourceContactBinding {
        source: SourceFiniteSolidMixedContactBinding {
            law,
            coefficient_role: SourceImpulseCoefficientRole::Kinetic,
        },
        texture: FiniteSolidContactTexture::FirstRigid {
            direction_body: Vec3::new(1.0, 0.0, 0.0),
        },
    }
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance<'a>(
    solid: &FiniteSolidBody<'a>,
    surface: &FiniteSolidContactSurface,
    scene: &FiniteSolidSphereScene,
    bodies: &[RigidBodyState],
    duration: f64,
    gravity: Gravity,
    controls: FiniteSolidSphereTrajectoryControls,
    history: &FiniteSolidSourceTrajectoryHistory<'a>,
    law: &SourceThermalEllipticCoulomb,
    boundary: impl Fn(
        &FiniteSolidSourceConditionContext<'_>,
    ) -> Result<SourceFrictionConditions, &'static str>,
    cx: &Cx<'_>,
) -> Result<
    (
        FiniteSolidBody<'a>,
        Vec<RigidBodyState>,
        FiniteSolidSourceTrajectoryHistory<'a>,
        FiniteSolidSustainedSourceTrajectoryReport,
    ),
    FiniteSolidSustainedSourceTrajectoryError<&'static str>,
> {
    solid.sustained_source_elliptic_sphere_trajectory_candidate(
        surface,
        scene,
        bodies,
        &vec![homogeneous_sphere(); bodies.len()],
        &[],
        [0.0; 3],
        &vec![Wrench::ZERO; bodies.len()],
        gravity,
        duration,
        ORIGIN,
        controls,
        history,
        Default::default(),
        |_, _| Ok(binding(law)),
        boundary,
        cx,
    )
}

#[test]
fn loaded_sliding_records_postcontact_point_slip_and_original_source_path_once() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 0.0, 0.0)], &[], 0).unwrap();
        let initial = body_at(
            [0.25, 0.25, -0.125 - GAP],
            Vec3::ZERO,
            Vec3::new(1.0, 0.0, 0.0),
        );
        let (_, law) = interface_source([1.0, 0.0]);
        let history = FiniteSolidSourceTrajectoryHistory::new(&solid, &scene);
        let (next, bodies, accepted, report) = advance(
            &solid,
            &surface,
            &scene,
            &[initial],
            4.0 * SOURCE_STEP,
            upward_gravity(),
            source_controls(SOURCE_STEP),
            &history,
            &law,
            |_| Ok(conditions(300.0)),
            cx,
        )
        .unwrap();
        let mut speed = 1.0;
        let mut distance = 0.0;
        for step in &report.steps {
            let event = step.event_data.as_ref().unwrap();
            assert_eq!(event.sources.len(), 1);
            assert_eq!(event.motion.len(), 1);
            assert_eq!(event.work.len(), 1);
            let impact = step.impact.as_ref().unwrap();
            let jn = 2.0 * ACCELERATION * step.duration_s;
            let coefficient = 0.2 + 0.05 * speed;
            let jt = -coefficient * jn;
            // m=2, I=1/80, R=1/8: Wtt=1/m+R²/I=7/4. The
            // moving contact point includes actual rigid spin, not just COM travel.
            near(event.sources[0].sample.queried_speed_m_s, speed, 2e-9);
            near(
                event.sources[0].sample.coefficients_mu[0],
                coefficient,
                2e-9,
            );
            near(impact.contact.impulses[0][0], jn, 2e-9);
            near(impact.contact.impulses[0][1], jt, 2e-9);
            let after_tangent = speed + 1.75 * jt;
            let work = &event.work[0];
            assert_eq!(work.contact_id, event.sources[0].kinematics.contact_id);
            assert_eq!(work.authority, FiniteSolidFrictionHeatAuthority::UniqueSustainedImpulse);
            near(
                work.tangential_impulse_work_j,
                jt * 0.5 * (speed + after_tangent),
                2e-10,
            );
            near(
                work.friction_heat_candidate_j.unwrap(),
                -work.tangential_impulse_work_j,
                2e-10,
            );
            assert!(work.friction_heat_candidate_j.unwrap() > 0.0);
            speed = after_tangent;
            let travel = speed * step.duration_s;
            distance += travel;
            let motion = event.motion[0].motion.as_ref().unwrap();
            near(motion.slip_increment_m[0], travel, 2e-12);
            near(motion.slip_increment_m[1], 0.0, 2e-12);
            near(motion.path_increment_m, travel, 2e-12);
            near(motion.queried_speed_m_s, speed, 2e-9);
            assert!(
                event.motion[0].maximum_clearance_m.unwrap()
                    <= source_controls(SOURCE_STEP).surface.activation_tolerance_m
            );
            for source in motion.state.friction_sources() {
                near(
                    source.state_coordinate("path").unwrap(),
                    0.25 + distance,
                    2e-12,
                );
                near(source.state_coordinate("speed").unwrap(), speed, 2e-9);
            }
        }
        let entry = accepted.entries().next().unwrap();
        assert_eq!(accepted.entries().len(), 1);
        near(
            entry.state.constitutive_state().slip_m()[0],
            distance,
            2e-12,
        );
        near(entry.last_contact_time_s, next.native_state().t, 0.0);
        near(accepted.time_s(), 4.0 * SOURCE_STEP, 0.0);
        let arm = bodies[0]
            .pose()
            .orientation()
            .rotate_world_to_body(Vec3::new(0.0, 0.0, 0.125));
        near(
            bodies[0]
                .point_kinematics(homogeneous_sphere(), arm)
                .unwrap()
                .point_velocity_world
                .x,
            speed,
            2e-9,
        );
        assert!(distance < bodies[0].pose().position_world().x - initial.pose().position_world().x);
        assert_eq!(history.entries().len(), 0);
        assert_eq!(history.time_s(), 0.0);
        assert_eq!(solid.native_state().steps, 0);
        near(report.work_energy_residual_j, 0.0, 2e-10);
    });
}

#[test]
fn force_split_refinement_advances_only_accepted_source_travel() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(
            &[sphere(1, 0, 0.0, 0.0), sphere(2, 1, 0.0, 0.0)],
            &[SpherePairLaw {
                first: 1,
                second: 2,
                friction_coefficient: 0.0,
                restitution: 0.0,
            }],
            1,
        )
        .unwrap();
        let bodies = [
            body_at(
                [0.15, 0.15, -0.125 - GAP],
                Vec3::ZERO,
                Vec3::new(0.1, 0.0, 0.0),
            ),
            body_at(
                [0.5, 0.1, -0.125 - GAP],
                Vec3::ZERO,
                Vec3::new(1.0, 0.0, 0.0),
            ),
        ];
        let (_, law) = interface_source([1.0, 0.0]);
        let history = FiniteSolidSourceTrajectoryHistory::new(&solid, &scene);
        let first_samples = Cell::new(0usize);
        let first_motion_calls = Cell::new(0usize);
        let mut controls = source_controls(4.0 * SOURCE_STEP);
        // Both resting normal supports contribute g²*T*h. This allowance
        // admits h=STEP but refuses 2*STEP and 4*STEP independently of friction.
        controls.integration_tolerance_j = 3.0 * ACCELERATION * ACCELERATION * DURATION * STEP;
        let (_, _, accepted, report) = advance(
            &solid,
            &surface,
            &scene,
            &bodies,
            4.0 * SOURCE_STEP,
            upward_gravity(),
            controls,
            &history,
            &law,
            |context| {
                if context.time_s == 0.0
                    && context.phase == FiniteSolidSourceConditionPhase::PreContact
                {
                    near(
                        context.state.friction_sources()[0]
                            .state_coordinate("path")
                            .unwrap(),
                        0.25,
                        0.0,
                    );
                }
                if matches!(
                    context.key,
                    FiniteSolidSphereSceneKey::Solid { sphere: 1, .. }
                ) {
                    match context.phase {
                        FiniteSolidSourceConditionPhase::PreContact => {
                            first_samples.set(first_samples.get() + 1)
                        }
                        FiniteSolidSourceConditionPhase::AcceptedMotion => {
                            first_motion_calls.set(first_motion_calls.get() + 1)
                        }
                    }
                }
                Ok(conditions(300.0))
            },
            cx,
        )
        .unwrap();
        assert!(report.rejected_steps >= 2);
        assert!(first_samples.get() > report.steps.len());
        assert_eq!(first_motion_calls.get(), report.steps.len());
        for entry in accepted.entries() {
            let mut travel = 0.0;
            let mut slip = [0.0; 2];
            let mut speed = match entry.key {
                FiniteSolidSphereSceneKey::Solid { sphere: 1, .. } => 0.1,
                FiniteSolidSphereSceneKey::Solid { sphere: 2, .. } => 1.0,
                other => panic!("unexpected active source feature: {other:?}"),
            };
            for step in &report.steps {
                assert!(step.duration_s <= STEP);
                let jt = -(0.2 + 0.05 * speed) * 2.0 * ACCELERATION * step.duration_s;
                speed += 1.75 * jt;
                for event in &step.event_data.as_ref().unwrap().motion {
                    if event.contact_id == entry.contact_id {
                        let motion = event.motion.as_ref().unwrap();
                        near(motion.path_increment_m, speed * step.duration_s, 2e-12);
                        travel += motion.path_increment_m;
                        for axis in 0..2 {
                            slip[axis] += motion.slip_increment_m[axis];
                        }
                    }
                }
            }
            near(
                entry.state.friction_sources()[0]
                    .state_coordinate("path")
                    .unwrap(),
                0.25 + travel,
                2e-12,
            );
            near(entry.state.constitutive_state().slip_m()[0], slip[0], 2e-12);
            near(entry.state.constitutive_state().slip_m()[1], slip[1], 2e-12);
        }
        assert_eq!(accepted.entries().len(), 2);
        assert_eq!(history.entries().len(), 0);
        assert_eq!(history.time_s(), 0.0);
    });
}

#[test]
fn source_endpoint_refusal_and_scene_clock_mismatch_preserve_accepted_history() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 0.0, 0.0)], &[], 0).unwrap();
        let body = body_at(
            [0.25, 0.25, -0.125 - GAP],
            Vec3::ZERO,
            Vec3::new(1.0, 0.0, 0.0),
        );
        let (_, law) = interface_source([1.0, 0.0]);
        let history = FiniteSolidSourceTrajectoryHistory::new(&solid, &scene);
        let (next, bodies, accepted, _) = advance(
            &solid,
            &surface,
            &scene,
            &[body],
            SOURCE_STEP,
            upward_gravity(),
            source_controls(SOURCE_STEP),
            &history,
            &law,
            |_| Ok(conditions(300.0)),
            cx,
        )
        .unwrap();
        let saved = accepted.entries().cloned().collect::<Vec<_>>();
        let saved_native = next.native_state().clone();
        let result = advance(
            &next,
            &surface,
            &scene,
            &bodies,
            2.0 * SOURCE_STEP,
            upward_gravity(),
            source_controls(SOURCE_STEP),
            &accepted,
            &law,
            |context| {
                Ok(conditions(
                    if context.phase == FiniteSolidSourceConditionPhase::AcceptedMotion
                        && context.time_s > 2.0 * SOURCE_STEP
                    {
                        340.0
                    } else {
                        300.0
                    },
                ))
            },
            cx,
        );
        assert!(matches!(
            result,
            Err(FiniteSolidSustainedTrajectoryError::Impact(
                FiniteSolidSustainedSourceFailure::Motion {
                    error: SourcePostcontactMotionError::Source(_),
                    ..
                }
            ))
        ));
        assert_eq!(accepted.entries().cloned().collect::<Vec<_>>(), saved);
        assert_eq!(accepted.time_s(), SOURCE_STEP);
        assert_eq!(next.native_state(), &saved_native);

        assert!(matches!(
            advance(
                &solid,
                &surface,
                &scene,
                &[body],
                SOURCE_STEP,
                upward_gravity(),
                source_controls(SOURCE_STEP),
                &accepted,
                &law,
                |_| panic!("stale clock must fail before source query"),
                cx
            ),
            Err(FiniteSolidSustainedTrajectoryError::Invalid(_))
        ));
        let changed = FiniteSolidSphereScene::new(&[sphere(10, 0, 0.0, 0.0)], &[], 0).unwrap();
        assert!(matches!(
            advance(
                &next,
                &surface,
                &changed,
                &bodies,
                SOURCE_STEP,
                upward_gravity(),
                source_controls(SOURCE_STEP),
                &accepted,
                &law,
                |_| panic!("changed scene must fail before source query"),
                cx
            ),
            Err(FiniteSolidSustainedTrajectoryError::Invalid(_))
        ));
        assert_eq!(accepted.entries().cloned().collect::<Vec<_>>(), saved);
    });
}

#[test]
fn separation_and_return_keep_feature_history_without_free_flight_travel() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 0.0, 0.0)], &[], 0).unwrap();
        let body = body_at(
            [0.25, 0.25, -0.125 - GAP],
            Vec3::ZERO,
            Vec3::new(1.0, 0.0, 0.0),
        );
        let (_, law) = interface_source([1.0, 0.0]);
        let history = FiniteSolidSourceTrajectoryHistory::new(&solid, &scene);
        let (s1, b1, h1, _) = advance(
            &solid,
            &surface,
            &scene,
            &[body],
            SOURCE_STEP,
            upward_gravity(),
            source_controls(SOURCE_STEP),
            &history,
            &law,
            |_| Ok(conditions(300.0)),
            cx,
        )
        .unwrap();
        let original_entry = h1.entries().next().unwrap().clone();
        let path_before = original_entry.state.friction_sources()[0]
            .state_coordinate("path")
            .unwrap();
        // An explicit native external kick separates the body. The trajectory
        // then admits its full free path without counting tangential free flight.
        let departing = b1[0]
            .apply_momentum_impulse(homogeneous_sphere(), Vec3::new(0.0, 0.0, -2.0), Vec3::ZERO)
            .unwrap()
            .state_after;
        let (s2, b2, h2, outward) = advance(
            &s1,
            &surface,
            &scene,
            &[departing],
            SOURCE_STEP,
            Gravity::ZERO,
            source_controls(SOURCE_STEP),
            &h1,
            &law,
            |_| Ok(conditions(300.0)),
            cx,
        )
        .unwrap();
        assert!(
            outward
                .steps
                .iter()
                .filter_map(|s| s.event_data.as_ref())
                .flat_map(|e| &e.motion)
                .all(|m| m.motion.is_none())
        );
        let separated = h2.get(original_entry.contact_id).unwrap();
        near(
            separated.state.friction_sources()[0]
                .state_coordinate("path")
                .unwrap(),
            path_before,
            0.0,
        );
        assert_eq!(
            separated.state.constitutive_state().slip_m(),
            original_entry.state.constitutive_state().slip_m()
        );
        let returning = b2[0]
            .apply_momentum_impulse(homogeneous_sphere(), Vec3::new(0.0, 0.0, 4.0), Vec3::ZERO)
            .unwrap()
            .state_after;
        let (s3, b3, h3, inward) = advance(
            &s2,
            &surface,
            &scene,
            &[returning],
            SOURCE_STEP,
            Gravity::ZERO,
            source_controls(SOURCE_STEP),
            &h2,
            &law,
            |_| Ok(conditions(300.0)),
            cx,
        )
        .unwrap();
        assert!(inward.steps.iter().all(|s| s.event_data.is_none()));
        near(
            h3.get(original_entry.contact_id)
                .unwrap()
                .state
                .friction_sources()[0]
                .state_coordinate("path")
                .unwrap(),
            path_before,
            0.0,
        );
        let (s4, _, h4, reentry) = advance(
            &s3,
            &surface,
            &scene,
            &b3,
            SOURCE_STEP,
            Gravity::ZERO,
            source_controls(SOURCE_STEP),
            &h3,
            &law,
            |_| Ok(conditions(300.0)),
            cx,
        )
        .unwrap();
        let added: f64 = reentry
            .steps
            .iter()
            .filter_map(|s| s.event_data.as_ref())
            .flat_map(|e| &e.motion)
            .filter_map(|m| m.motion.as_ref())
            .map(|m| m.path_increment_m)
            .sum();
        assert!(added > 0.0);
        let retained = h4.get(original_entry.contact_id).unwrap();
        assert_eq!(retained.key, original_entry.key);
        assert_eq!(h4.entries().len(), 1);
        near(
            retained.state.friction_sources()[0]
                .state_coordinate("path")
                .unwrap(),
            path_before + added,
            2e-12,
        );
        near(retained.last_contact_time_s, s4.native_state().t, 0.0);
        near(h4.time_s(), 4.0 * SOURCE_STEP, 0.0);
        near(h1.time_s(), SOURCE_STEP, 0.0);
    });
}

#[test]
fn loaded_rebound_keeps_instantaneous_source_sample_without_free_flight_travel() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        // This larger friction impulse still has a real finite-gap couple;
        // keep that couple within the original tight angular admission gate.
        let body = body_at(
            [0.25, 0.25, -0.125 - 0.25 * GAP],
            Vec3::ZERO,
            Vec3::new(4.0, 0.0, 1.0),
        );
        let (_, law) = interface_source([1.0, 0.0]);
        // Any positive physical restitution target is rebound, including a
        // target smaller than the solver's numerical velocity tolerance.
        for restitution in [1.0, 1e-12] {
            let scene =
                FiniteSolidSphereScene::new(&[sphere(9, 0, restitution, 0.0)], &[], 0).unwrap();
            let history = FiniteSolidSourceTrajectoryHistory::new(&solid, &scene);
            let (_, _, accepted, report) = advance(
                &solid,
                &surface,
                &scene,
                &[body],
                SOURCE_STEP,
                Gravity::ZERO,
                source_controls(SOURCE_STEP),
                &history,
                &law,
                |_| Ok(conditions(300.0)),
                cx,
            )
            .unwrap();
            let step = &report.steps[0];
            assert!(step.normal_targets_m_s[0] > 0.0);
            assert!(step.impact.as_ref().unwrap().contact.impulses[0][0] > 0.0);
            assert!(
                step.impact
                    .as_ref()
                    .unwrap()
                    .contact
                    .relative_velocities_after[0][1]
                    > 1.0
            );
            let event = step.event_data.as_ref().unwrap();
            near(event.sources[0].sample.queried_speed_m_s, 4.0, 2e-10);
            assert_eq!(event.work.len(), 1);
            assert_eq!(
                event.work[0].authority,
                FiniteSolidFrictionHeatAuthority::InstantaneousOrOpening
            );
            assert!(event.work[0].friction_heat_candidate_j.is_none());
            assert!(event.motion[0].motion.is_none());
            assert!(event.motion[0].maximum_clearance_m.is_none());
            let entry = accepted.entries().next().unwrap();
            assert_eq!(entry.state.constitutive_state().slip_m(), [0.0; 2]);
            for source in entry.state.friction_sources() {
                near(source.state_coordinate("path").unwrap(), 0.25, 0.0);
                near(source.state_coordinate("speed").unwrap(), 4.0, 2e-10);
            }
            near(accepted.time_s(), SOURCE_STEP, 0.0);
            assert_eq!(history.entries().len(), 0);
        }
    });
}

#[test]
fn material_spin_convects_source_axis_components_without_inventing_slip() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 0.0, 0.0)], &[], 0).unwrap();
        let omega = 256.0;
        let body = RigidBodyState::new(
            Pose::new(
                Vec3::new(0.25, 0.25, -0.125 - GAP),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, 0.0125 * omega),
        )
        .unwrap();
        let (_, law) = interface_source([0.6, 0.8]);
        let history = FiniteSolidSourceTrajectoryHistory::new(&solid, &scene);
        let (_, bodies, accepted, report) = advance(
            &solid,
            &surface,
            &scene,
            &[body],
            4.0 * SOURCE_STEP,
            upward_gravity(),
            source_controls(SOURCE_STEP),
            &history,
            &law,
            |_| Ok(conditions(300.0)),
            cx,
        )
        .unwrap();
        let entry = accepted.entries().next().unwrap();
        let axis = entry.state.constitutive_state().principal_axis();
        near(axis[0], 0.6, 2e-12);
        near(axis[1], 0.8, 2e-12);
        same(
            vec(entry.frame.longitudinal()),
            bodies[0]
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            2e-12,
        );
        assert!(norm(vec(entry.frame.longitudinal()).sub(Vec3::new(1.0, 0.0, 0.0))) > 1e-3);
        assert_eq!(entry.state.constitutive_state().slip_m(), [0.0; 2]);
        near(
            entry.state.friction_sources()[0]
                .state_coordinate("path")
                .unwrap(),
            0.25,
            0.0,
        );
        let mut spin = 0.0;
        for event in report.steps.iter().filter_map(|s| s.event_data.as_ref()) {
            for motion in &event.motion {
                near(motion.motion.as_ref().unwrap().path_increment_m, 0.0, 0.0);
                spin += motion.material_spin_angle_rad.unwrap();
            }
        }
        // Contact normal is -z, opposite the physical +z material rotation.
        near(spin, -omega * 4.0 * SOURCE_STEP, 2e-10);
    });
}
