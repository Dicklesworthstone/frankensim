use super::*;

#[test]
fn static_capture_and_kinetic_breakaway_follow_original_force_motion() {
    with_cx(|cx, _| {
        let h = 0.125;
        for load in [0.0, 9.0, 12.0, 14.0] {
            let mut solid = point(h, [0.0; 3], cx);
            let step = solid
                .advance_sustained_static_kinetic_contacts(
                    &[contact(h, 0.6, [0.0; 3])],
                    &[0.3],
                    &[[0.0; 3]],
                    &[load, -20.0, 0.0],
                    controls(),
                    cx,
                    admit_point,
                )
                .unwrap();
            let sliding = load > 12.0;
            let tangent_force = if sliding { -6.0 } else { -load };
            let acceleration = (load + tangent_force) / 2.0;
            near(step.solve.forces_n[0][0], 20.0, 1e-8);
            near(step.solve.forces_n[0][1], tangent_force, 1e-8);
            near(solid.velocity_mps()[0], acceleration * h, 1e-9);
            near(solid.displacement_m()[0], 0.5 * acceleration * h * h, 1e-10);
            near(
                step.frictional_dissipation_j,
                -tangent_force * 0.5 * acceleration * h * h,
                1e-9,
            );
            assert_eq!(
                step.friction_coefficients,
                [if sliding { 0.3 } else { 0.6 }]
            );
            assert_eq!(step.friction_regime_resolves, usize::from(sliding));
            assert_eq!(solid.accepted_steps(), 1);
            near(solid.time_s(), h, 0.0);
            energy(&step);
        }
    });
}

#[test]
fn incoming_slip_and_moving_surface_select_kinetic_friction_immediately() {
    with_cx(|cx, _| {
        let h = 0.125;
        // The same relative trajectory viewed in two different inertial frames.
        for (velocity, surface) in [([1.0, 0.0, 0.0], [0.0; 3]), ([0.0; 3], [0.0, -1.0, 0.0])] {
            let mut solid = point(h, velocity, cx);
            let step = solid
                .advance_sustained_static_kinetic_contacts(
                    &[contact(h, 0.6, surface)],
                    &[0.3],
                    &[[20.0, -12.0, 0.0]],
                    &[0.0, -20.0, 0.0],
                    controls(),
                    cx,
                    admit_point,
                )
                .unwrap();
            near(step.solve.forces_n[0][1], -6.0, 1e-8);
            near(solid.velocity_mps()[0], velocity[0] - 3.0 * h, 1e-9);
            near(step.frictional_dissipation_j, 6.0 * (h - 1.5 * h * h), 1e-9);
            assert_eq!(step.friction_coefficients, [0.3]);
            assert_eq!(step.friction_regime_resolves, 0);
            energy(&step);
        }
    });
}

#[test]
fn even_tiny_incoming_slip_cannot_borrow_static_strength() {
    with_cx(|cx, _| {
        let h = 0.125;
        let mut solid = point(h, [1e-12, 0.0, 0.0], cx);
        let step = solid
            .advance_sustained_static_kinetic_contacts(
                &[contact(h, 0.6, [0.0; 3])],
                &[0.0],
                &[[20.0, 0.0, 0.0]],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                admit_point,
            )
            .unwrap();
        assert_eq!(step.friction_coefficients, [0.0]);
        near(step.solve.forces_n[0][1], 0.0, 0.0);
        near(solid.velocity_mps()[0], 1e-12, 0.0);
        near(step.frictional_dissipation_j, 0.0, 0.0);
    });
}

#[test]
fn coupled_breakaway_cascade_rechecks_every_contact_and_keeps_input_order_irrelevant() {
    with_cx(|cx, _| {
        let h = 0.125;
        let mut mass = [0.0; 36];
        for a in 0..2 {
            for b in 0..2 {
                for axis in 0..3 {
                    mass[(3 * a + axis) * 6 + 3 * b + axis] = if a == b { 2.0 } else { 1.0 };
                }
            }
        }
        let pencil = assembly(&[0.0; 36], &mass, 6);
        let rows = [
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ];
        let contacts = [
            SolidSustainedContact {
                id: 3,
                velocity_weights: [&rows[0], &rows[1], &rows[2]],
                ..contact(h, 1.2, [0.0; 3])
            },
            SolidSustainedContact {
                id: 9,
                velocity_weights: [&rows[3], &rows[4], &rows[5]],
                ..contact(h, 0.25, [0.0; 3])
            },
        ];
        let mut baseline = None;
        for reverse in [false, true] {
            let mut solid =
                ExactLinearSolid::prepare(&pencil, h, &[0.0; 6], &[0.0; 6], native_controls(), cx)
                    .unwrap();
            let input = if reverse {
                [contacts[1], contacts[0]]
            } else {
                contacts
            };
            let kinetic = if reverse { [0.1, 0.6] } else { [0.6, 0.1] };
            let step = solid
                .advance_sustained_static_kinetic_contacts(
                    &input,
                    &kinetic,
                    &[[0.0; 3]; 2],
                    &[14.0, -10.0, 0.0, 0.0, -10.0, 0.0],
                    controls(),
                    cx,
                    |p| {
                        // With M=[[2,1],[1,2]] on each axis, the final forces give
                        // constant accelerations +5 and -2. Both normal paths are
                        // identically zero and neither tangent can reverse.
                        for (i, acceleration) in [(0, 5.0), (3, -2.0)] {
                            near(p.velocity_after_mps[i], h * acceleration, 1e-8);
                            near(p.displacement_change_m[i], 0.5 * h * h * acceleration, 1e-9);
                        }
                        for i in [1, 2, 4, 5] {
                            near(p.velocity_after_mps[i], 0.0, 1e-8);
                        }
                        Ok::<_, &'static str>(())
                    },
                )
                .unwrap();
            assert_eq!(step.contact_ids, [3, 9]);
            assert_eq!(step.friction_coefficients, [0.6, 0.1]);
            assert_eq!(step.friction_regime_resolves, 2);
            near(step.solve.forces_n[0][1], -6.0, 1e-7);
            near(step.solve.forces_n[1][1], 1.0, 1e-7);
            near(step.frictional_dissipation_j, 16.0 * h * h, 1e-8);
            assert_eq!(solid.accepted_steps(), 1);
            energy(&step);
            if let Some(ref expected) = baseline {
                assert_eq!(bits(&solid), *expected);
            } else {
                baseline = Some(bits(&solid));
            }
        }
    });
}

#[test]
fn equal_static_and_kinetic_bounds_retain_single_cone_results() {
    with_cx(|cx, _| {
        let h = 0.125;
        for (velocity, load) in [([0.0; 3], 5.0), ([0.0; 3], 14.0), ([1.0, 0.0, 0.0], 0.0)] {
            let contacts = [contact(h, 0.3, [0.0; 3])];
            let mut original = point(h, velocity, cx);
            let mut split = point(h, velocity, cx);
            let old = original
                .advance_sustained_contacts(
                    &contacts,
                    &[[0.0; 3]],
                    &[load, -20.0, 0.0],
                    controls(),
                    cx,
                    admit_point,
                )
                .unwrap();
            let new = split
                .advance_sustained_static_kinetic_contacts(
                    &contacts,
                    &[0.3],
                    &[[0.0; 3]],
                    &[load, -20.0, 0.0],
                    controls(),
                    cx,
                    admit_point,
                )
                .unwrap();
            assert_eq!(bits(&original), bits(&split));
            assert_eq!(old.solve.forces_n, new.solve.forces_n);
            assert_eq!(
                old.solve.response_evaluations,
                new.solve.response_evaluations
            );
            assert_eq!(new.friction_regime_resolves, 0);
        }
    });
}

#[test]
fn global_response_budget_invalid_pairs_and_path_refusal_leave_state_unchanged() {
    with_cx(|cx, _| {
        let h = 0.125;
        let contacts = [contact(h, 0.6, [0.0; 3])];
        let load = [14.0, -20.0, 0.0];
        let mut static_only = point(h, [0.0; 3], cx);
        let static_step = static_only
            .advance_sustained_contacts(&contacts, &[[0.0; 3]], &load, controls(), cx, admit_point)
            .unwrap();
        let mut solid = point(h, [0.0; 3], cx);
        let before = bits(&solid);
        let mut limited = controls();
        limited.solve.maximum_response_evaluations = static_step.solve.response_evaluations;
        let result = solid.advance_sustained_static_kinetic_contacts(
            &contacts,
            &[0.3],
            &[[0.0; 3]],
            &load,
            limited,
            cx,
            |_| -> Result<(), &str> { panic!("exhausted candidate cannot reach path admission") },
        );
        assert!(matches!(
            result,
            Err(SolidSustainedError::Budget {
                resource: "sustained contact response evaluations",
                ..
            })
        ));
        assert_eq!(bits(&solid), before);
        assert_eq!(solid.accepted_steps(), 0);
        for kinetic in [
            vec![],
            vec![-0.1],
            vec![0.7],
            vec![f64::NAN],
            vec![f64::INFINITY],
        ] {
            let result = solid.advance_sustained_static_kinetic_contacts(
                &contacts,
                &kinetic,
                &[[0.0; 3]],
                &load,
                controls(),
                cx,
                |_| -> Result<(), &str> { panic!("invalid friction cannot reach path admission") },
            );
            assert!(matches!(result, Err(SolidSustainedError::Invalid(_))));
            assert_eq!(bits(&solid), before);
        }
        let refusal = solid.advance_sustained_static_kinetic_contacts(
            &contacts,
            &[0.3],
            &[[0.0; 3]],
            &load,
            controls(),
            cx,
            |_| Err::<(), _>("whole path refused after breakaway"),
        );
        assert!(matches!(
            refusal,
            Err(SolidSustainedError::Path(
                "whole path refused after breakaway"
            ))
        ));
        assert_eq!(bits(&solid), before);
    });
}
