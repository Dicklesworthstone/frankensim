//! Actual moving native participants, splitting refinement and atomic refusal.
use super::*;
use fs_couple::finite_solid::contact::step::free::FreeFlightControls;
use fs_couple::finite_solid::gravity::reciprocal::trajectory::{
    ReciprocalGravityStepControls, ReciprocalGravityStepError,
};
use fs_mbd::{Gravity, Wrench};

fn step_controls(h: f64) -> ReciprocalGravityStepControls {
    ReciprocalGravityStepControls {
        gravity: ReciprocalGravityKickControls {
            duration_s: h,
            ..kick_controls()
        },
        middle: FreeFlightControls {
            solid: controls(1e-8),
            maximum_rigid_bodies: 1,
            ..Default::default()
        },
        maximum_duration_s: 0.1,
        maximum_nodes: 4,
        minimum_separation_m: 0.2,
        splitting_work_tolerance_j: 0.1,
        energy_tolerance_j: 0.1,
    }
}

fn total_momentum(
    solid: &FiniteSolidBody<'_>,
    source: RigidBodyState,
    origin: [f64; 3],
    cx: &Cx<'_>,
) -> (Vec3, Vec3) {
    let solid = solid.momentum_about(origin, cx).unwrap();
    (
        vector(solid.linear_kg_m_s).add(source.linear_momentum_world()),
        vector(solid.angular_kg_m2_s)
            .add(
                source
                    .pose()
                    .position_world()
                    .sub(vector(origin))
                    .cross(source.linear_momentum_world()),
            )
            .add(
                source
                    .pose()
                    .orientation()
                    .rotate_body_to_world(source.angular_momentum_body()),
            ),
    )
}

#[test]
fn rotating_extended_source_matches_native_subflows_and_exchanges_actual_spin() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (source, properties) = extended_source();
        let c = step_controls(0.02);
        let before = solid.native_state().clone();
        let (next, moved, report) = solid
            .reciprocal_quadrupole_gravity_step_candidate(
                source,
                properties,
                SOURCE_EXTENT,
                &[],
                [0.0; 3],
                None,
                c,
                cx,
            )
            .unwrap();
        assert_ne!(moved.pose().orientation(), source.pose().orientation());
        assert!(
            report
                .first_kick
                .source
                .angular_increment_body
                .x
                .hypot(report.first_kick.source.angular_increment_body.y)
                .hypot(report.first_kick.source.angular_increment_body.z)
                > 1e-7
        );
        assert!(
            report.path_separation_lower_bound_m > SOURCE_EXTENT.minimum_sample_radius_m().unwrap()
        );
        let half = ReciprocalGravityKickControls {
            duration_s: 0.5 * c.gravity.duration_s,
            ..c.gravity
        };
        let (a, s0, _) = solid
            .reciprocal_quadrupole_gravity_kick_candidate(
                source,
                properties,
                SOURCE_EXTENT,
                half,
                cx,
            )
            .unwrap();
        let (b, s1, _) = a
            .rigid_free_flight_candidate(
                &[s0],
                &[properties],
                &[],
                [0.0; 3],
                &[Wrench::ZERO],
                Gravity::ZERO,
                c.gravity.duration_s,
                c.gravity.balance_origin_m,
                c.middle,
                cx,
            )
            .unwrap();
        let (expected, expected_source, _) = b
            .reciprocal_quadrupole_gravity_kick_candidate(
                s1[0],
                properties,
                SOURCE_EXTENT,
                half,
                cx,
            )
            .unwrap();
        assert_eq!(next.native_state(), expected.native_state());
        assert_eq!(moved, expected_source);
        assert_eq!(next.native_state().a, b.native_state().a);
        assert_eq!(next.native_state().steps, 1);
        let (p0, l0) = total_momentum(&solid, source, c.gravity.balance_origin_m, cx);
        let (p1, l1) = total_momentum(&next, moved, c.gravity.balance_origin_m, cx);
        near_vector(components(p1.sub(p0)), [0.0; 3], 2e-10);
        near_vector(components(l1.sub(l0)), [0.0; 3], 2e-10);
        near(
            next.kinetic_energy_j().unwrap() - solid.kinetic_energy_j().unwrap()
                + moved.kinetic_energy(properties).unwrap()
                - source.kinetic_energy(properties).unwrap()
                + report.middle.strain_change_j
                + report.interaction_potential_change_j,
            report.integration_defect_j + report.middle.source_work_j,
            3e-12,
        );
        near(report.work_energy_residual_j, 0.0, 3e-12);
        let (shifted, shifted_source, _) = solid
            .reciprocal_quadrupole_gravity_step_candidate(
                source,
                properties,
                SOURCE_EXTENT,
                &[],
                [0.0; 3],
                None,
                ReciprocalGravityStepControls {
                    gravity: ReciprocalGravityKickControls {
                        balance_origin_m: [-1.3, 2.2, -0.6],
                        ..c.gravity
                    },
                    ..c
                },
                cx,
            )
            .unwrap();
        assert_eq!(shifted.native_state(), next.native_state());
        assert_eq!(shifted_source, moved);
        assert!(matches!(
            solid.reciprocal_quadrupole_gravity_step_candidate(
                source,
                properties,
                QuadrupoleSourceExtent {
                    radius_m: 1.2,
                    ..SOURCE_EXTENT
                },
                &[],
                [0.0; 3],
                None,
                c,
                cx
            ),
            Err(ReciprocalGravityStepError::Separation { .. })
        ));
        assert_eq!(solid.native_state(), &before);
    });
}

#[test]
fn actual_moving_source_step_matches_native_strang_subflows_and_pair_balances() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.4,
            Default::default(),
            cx,
        )
        .unwrap();
        let before = solid.native_state().clone();
        let (source, inertia) = source_at([-2.0, 0.3, -0.7]);
        let source =
            RigidBodyState::new(source.pose(), Vec3::ZERO, source.angular_momentum_body()).unwrap();
        let c = step_controls(0.02);
        let (next, moved, report) = solid
            .reciprocal_gravity_step_candidate(source, inertia, &[], [0.0; 3], None, c, cx)
            .unwrap();
        assert!(moved.pose().position_world().x > source.pose().position_world().x);
        assert!(next.displacement_m().iter().map(|u| u[0]).sum::<f64>() < 0.0);
        assert_ne!(
            report.first_kick.potential_energy_j,
            report.last_kick.potential_energy_j
        );
        near(next.native_state().t, 0.42, 1e-15);
        assert_eq!(next.native_state().steps, 1);
        assert_eq!(next.native_state().history.len(), 1);
        assert_eq!(solid.native_state(), &before);
        let half = ReciprocalGravityKickControls {
            duration_s: 0.01,
            ..c.gravity
        };
        let (a, s0, _) = solid
            .reciprocal_gravity_kick_candidate(source, inertia, half, cx)
            .unwrap();
        let (b, s1, middle) = a
            .rigid_free_flight_candidate(
                &[s0],
                &[inertia],
                &[],
                [0.0; 3],
                &[Wrench::ZERO],
                Gravity::ZERO,
                0.02,
                c.gravity.balance_origin_m,
                c.middle,
                cx,
            )
            .unwrap();
        let (expected, expected_source, _) = b
            .reciprocal_gravity_kick_candidate(s1[0], inertia, half, cx)
            .unwrap();
        assert_eq!(next.native_state(), expected.native_state());
        assert_eq!(moved, expected_source);
        assert_eq!(next.native_state().a, b.native_state().a);
        assert_eq!(report.middle.solid.time, middle.solid.time);
        let (p0, l0) = total_momentum(&solid, source, c.gravity.balance_origin_m, cx);
        let (p1, l1) = total_momentum(&next, moved, c.gravity.balance_origin_m, cx);
        near_vector(components(p1.sub(p0)), [0.0; 3], 2e-10);
        near_vector(components(l1.sub(l0)), [0.0; 3], 2e-10);
        near_vector(
            components(report.linear_momentum_change_kg_m_s),
            components(p1.sub(p0)),
            2e-10,
        );
        near_vector(
            components(report.angular_momentum_change_kg_m2_s),
            components(l1.sub(l0)),
            2e-10,
        );
        let physical_change = next.kinetic_energy_j().unwrap() - solid.kinetic_energy_j().unwrap()
            + moved.kinetic_energy(inertia).unwrap()
            - source.kinetic_energy(inertia).unwrap()
            + report.middle.strain_change_j
            + report.last_kick.potential_energy_j
            - report.first_kick.potential_energy_j;
        near(
            physical_change,
            report.integration_defect_j + report.middle.source_work_j,
            2e-12,
        );
        near(report.work_energy_residual_j, 0.0, 2e-12);
        assert!(
            report.absolute_integration_defect_j >= report.gravitational_splitting_defect_j.abs()
        );
        let shifted = ReciprocalGravityStepControls {
            gravity: ReciprocalGravityKickControls {
                balance_origin_m: [-1.1, 2.3, -0.6],
                ..c.gravity
            },
            ..c
        };
        let (translated, translated_source, _) = solid
            .reciprocal_gravity_step_candidate(source, inertia, &[], [0.0; 3], None, shifted, cx)
            .unwrap();
        assert_eq!(translated.native_state(), next.native_state());
        assert_eq!(translated_source, moved);
    });
}

#[test]
fn coupled_motion_refines_at_second_order_while_preserving_actual_pair_momentum() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let initial = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (source, inertia) = source_at([-3.0, 0.3, -1.0]);
        let source = RigidBodyState::new(source.pose(), Vec3::ZERO, Vec3::ZERO).unwrap();
        let run = |steps: usize, extended: bool| {
            let mut solid = initial.clone();
            let mut source = source;
            let mut c = step_controls(0.08 / steps as f64);
            c.gravity.gravitational_constant_m3_kg_s2 = 20.0;
            let initial_momentum = total_momentum(&solid, source, c.gravity.balance_origin_m, cx);
            for _ in 0..steps {
                let (next, next_source, report) = if extended {
                    solid.reciprocal_quadrupole_gravity_step_candidate(
                        source,
                        inertia,
                        SOURCE_EXTENT,
                        &[],
                        [0.0; 3],
                        None,
                        c,
                        cx,
                    )
                } else {
                    solid.reciprocal_gravity_step_candidate(
                        source,
                        inertia,
                        &[],
                        [0.0; 3],
                        None,
                        c,
                        cx,
                    )
                }
                .unwrap();
                near(report.work_energy_residual_j, 0.0, 2e-10);
                solid = next;
                source = next_source;
            }
            let final_momentum = total_momentum(&solid, source, c.gravity.balance_origin_m, cx);
            near_vector(
                components(final_momentum.0.sub(initial_momentum.0)),
                [0.0; 3],
                2e-9,
            );
            near_vector(
                components(final_momentum.1.sub(initial_momentum.1)),
                [0.0; 3],
                2e-9,
            );
            assert_eq!(solid.native_state().steps, steps);
            near(solid.native_state().t, 0.08, 2e-15);
            let mut state = solid.native_state().q.clone();
            state.extend(&solid.native_state().v);
            state.extend(components(source.pose().position_world()));
            state.extend(components(source.linear_momentum_world()));
            state.extend(components(source.angular_momentum_body()));
            state.extend(source.pose().orientation().components());
            state
        };
        let difference = |a: &[f64], b: &[f64]| {
            a.iter()
                .zip(b)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt()
        };
        for extended in [false, true] {
            let coarse = run(4, extended);
            let medium = run(8, extended);
            let fine = run(16, extended);
            let coarse_error = difference(&coarse, &medium);
            let fine_error = difference(&medium, &fine);
            assert!(
                fine_error > 1e-12,
                "refinement difference must resolve the physical trajectory"
            );
            let ratio = coarse_error / fine_error;
            assert!(
                (2.8..5.2).contains(&ratio),
                "symmetric splitting refinement ratio {ratio:e}, extended={extended}"
            );
        }
    });
}

#[test]
fn supported_motion_accounts_for_both_kick_reactions_and_native_middle_loads() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0)];
        let solid = FiniteSolidBody::new(
            problem(&material, &fixed, &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (source, inertia) = source_at([-2.0, 0.3, -0.7]);
        let c = step_controls(0.01);
        let (p0, l0) = total_momentum(&solid, source, c.gravity.balance_origin_m, cx);
        let (next, moved, report) = solid
            .reciprocal_gravity_step_candidate(source, inertia, &[(4, 0.2)], [0.0; 3], None, c, cx)
            .unwrap();
        let (p1, l1) = total_momentum(&next, moved, c.gravity.balance_origin_m, cx);
        assert_eq!(next.velocity_m_s()[0], [0.0; 3]);
        near_vector(
            components(p1.sub(p0)),
            components(
                report
                    .support_impulse_kg_m_s
                    .add(report.middle.external_impulse_kg_m_s),
            ),
            2e-10,
        );
        near_vector(
            components(l1.sub(l0)),
            components(
                report
                    .support_angular_impulse_kg_m2_s
                    .add(report.middle.external_angular_impulse_kg_m2_s),
            ),
            2e-10,
        );
        assert!(
            report
                .first_kick
                .solid
                .support
                .linear_kg_m_s
                .iter()
                .any(|x| x.abs() > 1e-6)
        );
        assert!(
            report
                .last_kick
                .solid
                .support
                .linear_kg_m_s
                .iter()
                .any(|x| x.abs() > 1e-6)
        );
        near_vector(
            components(report.linear_momentum_change_kg_m_s),
            components(p1.sub(p0)),
            2e-10,
        );
        near_vector(
            components(report.angular_momentum_change_kg_m2_s),
            components(l1.sub(l0)),
            2e-10,
        );
    });
}

#[test]
fn source_tunneling_late_energy_refusal_and_cancellation_are_atomic() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let before = solid.native_state().clone();
        let (source, inertia) = source_at([-2.0, 0.25, 0.25]);
        let flying =
            RigidBodyState::new(source.pose(), Vec3::new(150.0, 0.0, 0.0), Vec3::ZERO).unwrap();
        let mut c = step_controls(0.1);
        c.gravity.gravitational_constant_m3_kg_s2 = 1e-6;
        // Both endpoints lie well outside the solid, but the source crosses its
        // declared linear numerical path. Endpoint-only field admission misses it.
        assert!(matches!(
            solid.reciprocal_gravity_step_candidate(flying, inertia, &[], [0.0; 3], None, c, cx),
            Err(ReciprocalGravityStepError::Separation { .. })
        ));
        let (source, inertia) = source_at([-2.0, 0.3, -0.7]);
        let mut c = step_controls(0.02);
        c.splitting_work_tolerance_j = 1e-30;
        assert!(matches!(
            solid.reciprocal_gravity_step_candidate(source, inertia, &[], [0.0; 3], None, c, cx),
            Err(ReciprocalGravityStepError::Balance(_))
        ));
        assert_eq!(solid.native_state(), &before);
        let mut bounded = step_controls(0.02);
        bounded.maximum_nodes = 3;
        assert!(matches!(
            solid.reciprocal_gravity_step_candidate(
                source,
                inertia,
                &[],
                [0.0; 3],
                None,
                bounded,
                cx
            ),
            Err(ReciprocalGravityStepError::Invalid(_))
        ));
        gate.request();
        assert!(
            solid
                .reciprocal_gravity_step_candidate(
                    source,
                    inertia,
                    &[],
                    [0.0; 3],
                    None,
                    step_controls(0.02),
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.native_state(), &before);
    });
}
