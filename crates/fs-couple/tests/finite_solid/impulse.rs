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
            body.point_impulse_candidate(point, [0.0, 1.0, 0.0], [0.0; 3], Default::default(), cx)
                .is_err()
        );
        assert_eq!(body.native_state(), &original);
    });
}
