//! G1/G3 impulse and momentum oracles, independent of the production mass solve.
use super::*;
use fs_couple::finite_solid::impulse::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn vector(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for c in 0..3 {
        near(a[c], b[c], tolerance);
    }
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|c| a[c] + b[c])
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|c| a[c] - b[c])
}
fn scale(a: [f64; 3], b: f64) -> [f64; 3] {
    a.map(|v| v * b)
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|c| a[c] * b[c]).sum()
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn mbd(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

#[test]
fn off_center_impulse_uses_consistent_mass_current_geometry_and_conjugate_work() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false)); // total mass 1 kg
        let mut u = [[0.0; 3]; 4];
        for (a, x) in NODES.iter().enumerate() {
            // Finite shear/stretch and a translation, not a reference-origin fixture.
            u[a] = [
                0.2 + 0.3 * x[0] + 0.2 * x[1],
                -0.3 - 0.1 * x[1],
                0.4 + 0.1 * x[2],
            ];
        }
        let x: [[f64; 3]; 4] = std::array::from_fn(|a| add(NODES[a], u[a]));
        let v: [[f64; 3]; 4] =
            std::array::from_fn(|a| add([0.2, -0.1, 0.3], cross([0.3, -0.4, 0.5], x[a])));
        let body = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let before = body.native_state().clone();
        let weights = [0.1, 0.2, 0.3, 0.4];
        let point = FiniteSolidPoint::new(0, weights).unwrap();
        let j = [0.2, -0.3, 0.15];
        let origin = [0.7, -0.8, 0.4];
        let xp: [f64; 3] = std::array::from_fn(|c| (0..4).map(|a| weights[a] * x[a][c]).sum());
        let vp: [f64; 3] = std::array::from_fn(|c| (0..4).map(|a| weights[a] * v[a][c]).sum());
        vector(body.point_position_m(point, cx).unwrap(), xp, 2e-15);
        vector(body.point_velocity_m_s(point, cx).unwrap(), vp, 2e-15);
        // For one tet M=(m/20)(I+11^T), M^-1=(20/m)(I-11^T/5).
        let inverse_mass = 20.0 * (weights.iter().map(|w| w * w).sum::<f64>() - 0.2);
        near(
            body.point_directional_inverse_mass_per_kg(
                point,
                [0.0, 1.0, 0.0],
                Default::default(),
                cx,
            )
            .unwrap(),
            inverse_mass,
            3e-13,
        );
        near(
            body.point_directional_inverse_mass_per_kg(
                point,
                [0.0, 3.0, 0.0],
                Default::default(),
                cx,
            )
            .unwrap(),
            9.0 * inverse_mass,
            3e-12,
        );
        let (next, report) = body
            .point_impulse_candidate(point, j, origin, Default::default(), cx)
            .unwrap();
        for a in 0..4 {
            vector(
                next.velocity_m_s()[a],
                add(v[a], scale(j, 20.0 * (weights[a] - 0.2))),
                2e-13,
            );
        }
        vector(
            sub(report.after.linear_kg_m_s, report.before.linear_kg_m_s),
            j,
            2e-14,
        );
        vector(
            sub(report.after.angular_kg_m2_s, report.before.angular_kg_m2_s),
            cross(sub(xp, origin), j),
            3e-14,
        );
        vector(
            report.applied.angular_kg_m2_s,
            cross(sub(xp, origin), j),
            2e-15,
        );
        let expected_work = dot(j, vp) + 0.5 * inverse_mass * dot(j, j);
        near(report.kinetic_change_j, expected_work, 4e-14);
        near(report.impulse_work_j, expected_work, 4e-14);
        near(
            report.kinetic_after_j - report.kinetic_before_j,
            expected_work,
            4e-14,
        );
        assert_eq!(body.native_state(), &before);
        assert_eq!(next.native_state().q, before.q);
        assert_eq!(next.native_state().a, before.a);
        assert_eq!(next.native_state().t, before.t);
        assert_eq!(next.native_state().steps, before.steps);
        assert_eq!(next.native_state().history, before.history);
        let shifted_origin = add(origin, [1.0, -0.2, 0.8]);
        let shifted = next.momentum_about(shifted_origin, cx).unwrap();
        vector(
            shifted.angular_kg_m2_s,
            sub(
                report.after.angular_kg_m2_s,
                cross(sub(shifted_origin, origin), report.after.linear_kg_m_s),
            ),
            3e-14,
        );
        // Same physical trajectory continues through the original time owner.
        let (continued, _) = next
            .step_candidate(1e-4, &[], [0.0; 3], controls(1e-6), cx)
            .unwrap();
        let observed = continued.momentum_about(origin, cx).unwrap();
        vector(observed.linear_kg_m_s, report.after.linear_kg_m_s, 3e-11);
        vector(
            observed.angular_kg_m2_s,
            report.after.angular_kg_m2_s,
            3e-11,
        );
    });
}

#[test]
fn stationary_support_impulses_close_both_momenta_without_doing_work() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0)];
        let body = FiniteSolidBody::new(
            problem(&material, &fixed, &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let point = FiniteSolidPoint::new(0, [0.25; 4]).unwrap();
        let j = [0.4, -0.2, 0.3];
        let (next, report) = body
            .point_impulse_candidate(point, j, [0.2, 0.3, -0.4], Default::default(), cx)
            .unwrap();
        vector(next.velocity_m_s()[0], [0.0; 3], 0.0);
        // Free 3x3 block has inverse 20(I-11^T/4), giving each free v=1.25 J.
        for v in &next.velocity_m_s()[1..] {
            vector(*v, scale(j, 1.25), 1e-13);
        }
        vector(report.support.linear_kg_m_s, scale(j, -1.0 / 16.0), 1e-14);
        vector(
            report.after.linear_kg_m_s,
            add(report.applied.linear_kg_m_s, report.support.linear_kg_m_s),
            2e-14,
        );
        vector(
            report.after.angular_kg_m2_s,
            add(
                report.applied.angular_kg_m2_s,
                report.support.angular_kg_m2_s,
            ),
            2e-14,
        );
        near(report.impulse_work_j, 0.5 * 0.9375 * dot(j, j), 2e-14);
        near(report.work_residual_j, 0.0, 1e-14);
        near(
            body.point_directional_inverse_mass_per_kg(
                FiniteSolidPoint::new(0, [1.0, 0.0, 0.0, 0.0]).unwrap(),
                [1.0, 0.0, 0.0],
                Default::default(),
                cx,
            )
            .unwrap(),
            0.0,
            0.0,
        );
    });
}

#[test]
fn equal_opposite_rigid_solid_contact_impulse_preserves_total_angular_momentum() {
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
        let point = FiniteSolidPoint::new(0, [0.0, 0.2, 0.3, 0.5]).unwrap();
        let x = solid.point_position_m(point, cx).unwrap();
        let properties = MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(0.3, 0.4, 0.5),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap(),
        )
        .unwrap();
        let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.3, -0.1, 0.7), 0.8).unwrap();
        let center = Vec3::new(0.7, -0.2, 0.9);
        let rigid = RigidBodyState::new(
            Pose::new(center, rotation).unwrap(),
            Vec3::new(-0.2, 0.1, 0.0),
            Vec3::new(0.1, 0.2, 0.3),
        )
        .unwrap();
        let arm = rotation.rotate_world_to_body(mbd(x).sub(center));
        let direction = [0.0, 0.0, 1.0];
        let rigid_mass = rigid
            .directional_effective_mass(properties, arm, mbd(direction))
            .unwrap();
        let solid_mass = solid
            .point_directional_inverse_mass_per_kg(point, direction, Default::default(), cx)
            .unwrap();
        let relative_before = -rigid
            .point_kinematics(properties, arm)
            .unwrap()
            .point_velocity_world
            .z;
        // A bilateral sticking impulse computed from the two original mass owners.
        let j = scale(
            direction,
            -relative_before / (solid_mass + rigid_mass.inverse_mass_kg_inverse),
        );
        let origin = [-0.3, 0.2, 0.1];
        let (after_solid, s) = solid
            .point_impulse_candidate(point, j, origin, Default::default(), cx)
            .unwrap();
        let r = rigid
            .apply_impulse_at_body_point(properties, arm, mbd(scale(j, -1.0)))
            .unwrap();
        let rigid_l = |b: RigidBodyState| {
            add(
                array(
                    b.pose()
                        .orientation()
                        .rotate_body_to_world(b.angular_momentum_body()),
                ),
                cross(
                    sub(array(b.pose().position_world()), origin),
                    array(b.linear_momentum_world()),
                ),
            )
        };
        vector(
            add(s.after.angular_kg_m2_s, rigid_l(r.state_after)),
            add(s.before.angular_kg_m2_s, rigid_l(rigid)),
            4e-14,
        );
        vector(
            add(
                s.after.linear_kg_m_s,
                array(r.state_after.linear_momentum_world()),
            ),
            add(s.before.linear_kg_m_s, array(rigid.linear_momentum_world())),
            4e-14,
        );
        near(
            after_solid.point_velocity_m_s(point, cx).unwrap()[2]
                - r.point_kinematics_after.point_velocity_world.z,
            0.0,
            2e-14,
        );
        let loss = -0.5 * relative_before * relative_before
            / (solid_mass + rigid_mass.inverse_mass_kg_inverse);
        near(
            s.kinetic_change_j + r.work.kinetic_energy_change_j,
            loss,
            5e-14,
        );
    });
}

#[test]
fn full_point_mobility_matches_analytic_consistent_mass_with_component_supports() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false)); // total mass 1 kg
        let rows = [
            (
                FiniteSolidPoint::new(0, [1.0, 0.0, 0.0, 0.0]).unwrap(),
                [1.0, 0.0, 0.0],
            ),
            (
                FiniteSolidPoint::new(0, [0.0, 1.0, 0.0, 0.0]).unwrap(),
                [1.0, 1.0, 0.0],
            ),
            (
                FiniteSolidPoint::new(0, [0.1, 0.2, 0.3, 0.4]).unwrap(),
                [0.3, -0.4, 0.5],
            ),
            (
                FiniteSolidPoint::new(0, [0.25; 4]).unwrap(),
                [-0.5, 0.7, 0.2],
            ),
        ];
        // Different supports per component produce anisotropic point mobility.
        for fixed in [vec![], vec![(0, 0.0), (2, 0.0), (4, 0.0), (8, 0.0)]] {
            let body = FiniteSolidBody::new(
                problem(&material, &fixed, &[]),
                &[[0.0; 3]; 4],
                &[[0.0; 3]; 4],
                0.0,
                Default::default(),
                cx,
            )
            .unwrap();
            let original = body.native_state().clone();
            let mobility = body
                .point_mobility_matrix_per_kg(&rows, rows.len(), Default::default(), cx)
                .unwrap();
            assert_eq!(mobility.dimension, rows.len());
            assert_eq!(mobility.matrix_per_kg.len(), rows.len() * rows.len());
            assert_eq!(mobility.solves.len(), rows.len());
            for (i, &(point_i, direction_i)) in rows.iter().enumerate() {
                for (j, &(point_j, direction_j)) in rows.iter().enumerate() {
                    let wi = point_i.weights();
                    let wj = point_j.weights();
                    let mut expected = 0.0;
                    for c in 0..3 {
                        let free: Vec<usize> = (0..4)
                            .filter(|a| !fixed.iter().any(|&(dof, _)| dof == 3 * a + c))
                            .collect();
                        let sum_i: f64 = free.iter().map(|&a| wi[a]).sum();
                        let sum_j: f64 = free.iter().map(|&a| wj[a]).sum();
                        let product: f64 = free.iter().map(|&a| wi[a] * wj[a]).sum();
                        // For k free vertices, M_ff^-1 = 20(I - 11^T/(k+1)).
                        expected += 20.0
                            * direction_i[c]
                            * direction_j[c]
                            * (product - sum_i * sum_j / (free.len() + 1) as f64);
                    }
                    near(mobility.matrix_per_kg[i * rows.len() + j], expected, 3e-12);
                }
            }
            if fixed.is_empty() {
                // The consistent mass transmits a NEGATIVE cross-vertex response;
                // diagonal point approximations or nodal lumping cannot pass.
                near(mobility.matrix_per_kg[1], -4.0, 1e-12);
            } else {
                for j in 0..rows.len() {
                    near(mobility.matrix_per_kg[j], 0.0, 0.0);
                    near(mobility.matrix_per_kg[j * rows.len()], 0.0, 0.0);
                }
            }
            assert_eq!(body.native_state(), &original);
        }
    });
}

#[test]
fn simultaneous_point_impulses_retain_cross_point_work_and_current_moments() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let u: [[f64; 3]; 4] = std::array::from_fn(|a| {
            [
                0.1 + 0.2 * NODES[a][0],
                -0.2 + 0.1 * NODES[a][1],
                0.3 - 0.1 * NODES[a][2],
            ]
        });
        let x: [[f64; 3]; 4] = std::array::from_fn(|a| add(NODES[a], u[a]));
        let v: [[f64; 3]; 4] =
            std::array::from_fn(|a| add([0.2, -0.1, 0.3], cross([0.3, -0.4, 0.5], x[a])));
        let body = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let original = body.native_state().clone();
        let points = [
            FiniteSolidPoint::new(0, [1.0, 0.0, 0.0, 0.0]).unwrap(),
            FiniteSolidPoint::new(0, [0.0, 1.0, 0.0, 0.0]).unwrap(),
            FiniteSolidPoint::new(0, [0.1, 0.2, 0.3, 0.4]).unwrap(),
        ];
        let directions = [[1.0, 0.2, -0.3], [0.5, -0.4, 0.6], [-0.3, 0.7, 0.1]];
        let amplitudes = [0.3, -0.2, 0.4];
        let rows: [(FiniteSolidPoint, [f64; 3]); 3] =
            std::array::from_fn(|i| (points[i], directions[i]));
        let j: [[f64; 3]; 3] = std::array::from_fn(|i| scale(directions[i], amplitudes[i]));
        let mobility = body
            .point_mobility_matrix_per_kg(&rows, 3, Default::default(), cx)
            .unwrap();
        let origin = [-0.3, 0.2, 0.1];
        let (next, report) = body
            .point_impulses_candidate(
                &[
                    (points[0], scale(j[0], 0.5)),
                    (points[1], j[1]),
                    (points[0], scale(j[0], 0.5)), // repeated point accumulates
                    (points[2], j[2]),
                ],
                origin,
                Default::default(),
                cx,
            )
            .unwrap();
        let mut expected_work = 0.0;
        let mut diagonal_work = 0.0;
        let mut expected_linear = [0.0; 3];
        let mut expected_angular = [0.0; 3];
        for i in 0..3 {
            let wi = points[i].weights();
            let xp = std::array::from_fn(|c| (0..4).map(|a| wi[a] * x[a][c]).sum());
            let vp = std::array::from_fn(|c| (0..4).map(|a| wi[a] * v[a][c]).sum());
            expected_linear = add(expected_linear, j[i]);
            expected_angular = add(expected_angular, cross(sub(xp, origin), j[i]));
            expected_work += dot(vp, j[i]);
            diagonal_work += dot(vp, j[i]);
            let mut expected_jump = [0.0; 3];
            for k in 0..3 {
                let wk = points[k].weights();
                let coupling = 20.0 * ((0..4).map(|a| wi[a] * wk[a]).sum::<f64>() - 0.2);
                expected_jump = add(expected_jump, scale(j[k], coupling));
                let work = 0.5 * coupling * dot(j[i], j[k]);
                expected_work += work;
                if i == k {
                    diagonal_work += work;
                }
            }
            let actual_jump = sub(next.point_velocity_m_s(points[i], cx).unwrap(), vp);
            vector(actual_jump, expected_jump, 3e-12);
            let matrix_jump: f64 = (0..3)
                .map(|k| mobility.matrix_per_kg[3 * i + k] * amplitudes[k])
                .sum();
            near(dot(directions[i], actual_jump), matrix_jump, 3e-12);
        }
        assert!((expected_work - diagonal_work).abs() > 0.05);
        near(report.kinetic_change_j, expected_work, 3e-12);
        near(report.impulse_work_j, expected_work, 3e-12);
        vector(
            sub(report.after.linear_kg_m_s, report.before.linear_kg_m_s),
            expected_linear,
            3e-12,
        );
        vector(
            sub(report.after.angular_kg_m2_s, report.before.angular_kg_m2_s),
            expected_angular,
            3e-12,
        );
        vector(report.support.linear_kg_m_s, [0.0; 3], 0.0);
        assert_eq!(body.native_state(), &original);
        assert_eq!(next.native_state().q, original.q);
        assert_eq!(next.native_state().a, original.a);
        assert_eq!(next.native_state().t, original.t);
        assert_eq!(next.native_state().steps, original.steps);
        assert_eq!(next.native_state().history, original.history);
    });
}

#[test]
fn invalid_unresolved_and_cancelled_impulses_leave_the_original_state_unchanged() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let body = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[1e20, 0.0, 0.0]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let original = body.native_state().clone();
        let point = FiniteSolidPoint::new(0, [0.25; 4]).unwrap();
        let empty = body
            .point_mobility_matrix_per_kg(&[], 0, Default::default(), cx)
            .unwrap();
        assert_eq!(empty.dimension, 0);
        assert!(empty.matrix_per_kg.is_empty());
        assert!(empty.solves.is_empty());
        assert!(
            body.point_mobility_matrix_per_kg(
                &[(point, [1.0, 0.0, 0.0])],
                0,
                Default::default(),
                cx
            )
            .is_err()
        );
        for direction in [[0.0; 3], [f64::NAN, 0.0, 0.0]] {
            assert!(
                body.point_mobility_matrix_per_kg(&[(point, direction)], 1, Default::default(), cx)
                    .is_err()
            );
        }
        assert!(
            body.point_impulses_candidate(
                &[
                    (point, [1.0, 0.0, 0.0]),
                    (FiniteSolidPoint::new(2, [0.25; 4]).unwrap(), [0.0; 3])
                ],
                [0.0; 3],
                Default::default(),
                cx,
            )
            .is_err()
        );
        assert!(FiniteSolidPoint::new(0, [0.2; 4]).is_err());
        assert!(FiniteSolidPoint::new(0, [-0.1, 0.1, 0.5, 0.5]).is_err());
        assert!(
            body.point_position_m(FiniteSolidPoint::new(2, [0.25; 4]).unwrap(), cx)
                .is_err()
        );
        assert!(
            body.nodal_impulse_candidate(&[(4, [0.0; 3])], [0.0; 3], Default::default(), cx)
                .is_err()
        );
        let mut too_small = FiniteSolidImpulseControls::default();
        too_small.mass_solve.maximum_workspace_components = 0;
        assert!(matches!(
            body.point_mobility_matrix_per_kg(
                &[(point, [1.0, 0.0, 0.0])],
                1,
                too_small.mass_solve,
                cx
            ),
            Err(FiniteSolidError::Initialization(_))
        ));
        assert!(matches!(
            body.point_impulse_candidate(point, [1.0, 0.0, 0.0], [0.0; 3], too_small, cx),
            Err(FiniteSolidError::Initialization(_))
        ));
        // Adding a small impulse to huge existing velocity must not silently lose it.
        assert!(matches!(
            body.point_impulse_candidate(point, [1.0, 0.0, 0.0], [0.0; 3], Default::default(), cx),
            Err(FiniteSolidError::ImpulseBalance { .. })
        ));
        assert_eq!(body.native_state(), &original);
        gate.request();
        assert!(
            body.point_mobility_matrix_per_kg(
                &[(point, [1.0, 0.0, 0.0])],
                1,
                Default::default(),
                cx
            )
            .is_err()
        );
        assert!(
            body.point_impulses_candidate(
                &[(point, [1.0, 0.0, 0.0])],
                [0.0; 3],
                Default::default(),
                cx
            )
            .is_err()
        );
        assert!(
            body.point_impulse_candidate(point, [0.0, 1.0, 0.0], [0.0; 3], Default::default(), cx)
                .is_err()
        );
        assert_eq!(body.native_state(), &original);
    });
}
