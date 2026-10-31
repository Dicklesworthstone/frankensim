//! G1/G3/G4: native free-flight work, spin/orbital balance and atomic refusal.
use super::super::*;
use fs_couple::finite_solid::contact::step::free::*;
use fs_mbd::{Gravity, Wrench};

fn flight_controls(energy_tolerance_j: f64) -> FreeFlightControls {
    FreeFlightControls {
        solid: controls(energy_tolerance_j),
        energy_tolerance_j,
        ..FreeFlightControls::default()
    }
}
fn total_momentum(
    solid: &FiniteSolidBody<'_>,
    bodies: &[RigidBodyState],
    origin: [f64; 3],
    cx: &Cx<'_>,
) -> (Vec3, Vec3) {
    let s = solid.momentum_about(origin, cx).unwrap();
    let mut p = vec(s.linear_kg_m_s);
    let mut l = vec(s.angular_kg_m2_s);
    for body in bodies {
        p = p.add(body.linear_momentum_world());
        l = l
            .add(
                body.pose()
                    .position_world()
                    .sub(vec(origin))
                    .cross(body.linear_momentum_world()),
            )
            .add(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(body.angular_momentum_body()),
            );
    }
    (p, l)
}
fn flight_audit(report: &FreeFlightReport) {
    near(norm(report.linear_momentum_residual_kg_m_s), 0.0, 2e-10);
    near(norm(report.angular_momentum_residual_kg_m2_s), 0.0, 2e-10);
    near(report.work_energy_residual_j, 0.0, 2e-10);
    near(report.absolute_linear_residual_kg_m_s, 0.0, 2e-10);
    near(report.absolute_angular_residual_kg_m2_s, 0.0, 2e-10);
    near(report.absolute_work_residual_j, 0.0, 2e-10);
    near(
        report.kinetic_change_j + report.strain_change_j,
        report.source_work_j + report.integration_defect_j,
        2e-10,
    );
}

#[test]
fn held_native_forces_and_axial_torque_match_constant_acceleration_and_rotation() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let mut prepared = problem(&material, &[], &[]);
        prepared.body_force_n_m3 = [0.0, 6.0, 0.0];
        let velocity = [0.5, -0.25, 0.125];
        let solid = FiniteSolidBody::new(
            prepared,
            &[[0.0; 3]; 4],
            &[velocity; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let property = properties();
        let x0 = Vec3::new(1.0, -0.5, 0.75);
        let v0 = Vec3::new(-0.2, 0.1, 0.3);
        let body = RigidBodyState::new(
            Pose::new(
                x0,
                UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.2).unwrap(),
            )
            .unwrap(),
            v0.scale(2.0),
            Vec3::new(0.0, 0.0, 0.4),
        )
        .unwrap();
        let duration = 0.125;
        let origin = [0.2, -0.3, 0.4];
        let wrench = Wrench {
            force_world: Vec3::new(1.0, 2.0, 3.0),
            torque_body: Vec3::new(0.0, 0.0, 0.5),
        };
        let gravity = Gravity::new(Vec3::new(0.0, 0.0, -3.0)).unwrap();
        let (next, rigid, report) = solid
            .rigid_free_flight_candidate(
                &[body],
                &[property],
                &[],
                [0.0, 6.0, 0.0],
                &[wrench],
                gravity,
                duration,
                origin,
                flight_controls(1e-8),
                cx,
            )
            .unwrap();
        // rho*V=1 kg, so the solid accelerates at +1 m/s² in y. The rigid
        // total force includes its own m*g once, hence F=(1,2,-3) N.
        let solid_shift =
            vec(velocity)
                .scale(duration)
                .add(Vec3::new(0.0, 0.5 * duration * duration, 0.0));
        for (q, v) in next.displacement_m().into_iter().zip(next.velocity_m_s()) {
            same(vec(q), solid_shift, 2e-11);
            same(
                vec(v),
                vec(velocity).add(Vec3::new(0.0, duration, 0.0)),
                2e-11,
            );
        }
        let force = Vec3::new(1.0, 2.0, -3.0);
        let rigid_shift = v0
            .scale(duration)
            .add(force.scale(0.25 * duration * duration));
        same(rigid[0].pose().position_world(), x0.add(rigid_shift), 2e-12);
        same(
            rigid[0].linear_momentum_world(),
            body.linear_momentum_world().add(force.scale(duration)),
            2e-12,
        );
        near(
            rigid[0].angular_momentum_body().z,
            0.4 + 0.5 * duration,
            2e-12,
        );
        let angle = 0.2 + duration * (0.4 + 0.25 * duration) / 0.5;
        same(
            rigid[0]
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            Vec3::new(angle.cos(), angle.sin(), 0.0),
            2e-12,
        );
        let torque_work = duration * 0.5 * (0.4 + 0.25 * duration) / 0.5;
        near(
            report.source_work_j,
            solid_shift.y + force.dot(rigid_shift) + torque_work,
            2e-11,
        );
        let solid_impulse = Vec3::new(0.0, duration, 0.0);
        let rigid_impulse = force.scale(duration);
        same(
            report.external_impulse_kg_m_s,
            solid_impulse.add(rigid_impulse),
            2e-12,
        );
        let solid_midpoint = Vec3::new(0.25, 0.25, 0.25).add(solid_shift.scale(0.5));
        let expected_angular = solid_midpoint
            .sub(vec(origin))
            .cross(solid_impulse)
            .add(
                x0.add(rigid_shift.scale(0.5))
                    .sub(vec(origin))
                    .cross(rigid_impulse),
            )
            .add(Vec3::new(0.0, 0.0, 0.5 * duration));
        same(
            report.external_angular_impulse_kg_m2_s,
            expected_angular,
            2e-11,
        );
        near(report.rigid_drift_energy_defect_j, 0.0, 2e-12);
        near(norm(report.support_impulse_kg_m_s), 0.0, 1e-14);
        flight_audit(&report);
        assert_eq!(solid.native_state().t, 0.0);
        assert_eq!(next.native_state().t, duration);
        assert_eq!(next.native_state().steps, solid.native_state().steps + 1);
    });
}

#[test]
fn deforming_consistent_solid_and_asymmetric_rigid_preserve_actual_total_momentum() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let u = [
            [0.0; 3],
            [0.08, 0.025, -0.01],
            [0.015, -0.04, 0.02],
            [-0.005, 0.01, 0.05],
        ];
        let v = [
            [0.2, -0.1, 0.3],
            [0.0, 0.25, 0.1],
            [-0.15, 0.1, -0.1],
            [0.05, -0.2, 0.2],
        ];
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let property = MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(0.3, 0.4, 0.5),
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.35).unwrap(),
        )
        .unwrap();
        let body = RigidBodyState::new(
            Pose::new(
                Vec3::new(1.25, -0.5, 0.75),
                UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.4).unwrap(),
            )
            .unwrap(),
            Vec3::new(0.3, -0.4, 0.2),
            Vec3::new(0.2, 0.35, -0.1),
        )
        .unwrap();
        let origin = [-0.2, 0.4, 0.7];
        let before = total_momentum(&solid, &[body], origin, cx);
        let (next, rigid, report) = solid
            .rigid_free_flight_candidate(
                &[body],
                &[property],
                &[],
                [0.0; 3],
                &[Wrench::ZERO],
                Gravity::ZERO,
                0.001,
                origin,
                flight_controls(1e-7),
                cx,
            )
            .unwrap();
        let after = total_momentum(&next, &rigid, origin, cx);
        same(after.0, before.0, 2e-10);
        same(after.1, before.1, 2e-10);
        same(
            report.angular_momentum_change_kg_m2_s,
            after.1.sub(before.1),
            2e-11,
        );
        assert_ne!(rigid[0].pose().orientation(), body.pose().orientation());
        assert_ne!(next.native_state().q, solid.native_state().q);
        assert!(report.absolute_integration_defect_j >= report.integration_defect_j.abs());
        assert!(report.absolute_integration_defect_j <= 1e-7);
        near(report.source_work_j, 0.0, 1e-14);
        flight_audit(&report);
    });
}

#[test]
fn stationary_support_reactions_close_linear_and_angular_balance_without_work() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let force = [(3, 0.4)];
        let mut u = [[0.0; 3]; 4];
        u[1][0] = 0.1;
        let mut v = [[0.0; 3]; 4];
        v[1][0] = 0.2;
        let solid = FiniteSolidBody::new(
            problem(&material, &fixed, &force),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let origin = [0.3, -0.2, 0.5];
        let before = total_momentum(&solid, &[], origin, cx);
        let (next, rigid, report) = solid
            .rigid_free_flight_candidate(
                &[],
                &[],
                &force,
                [0.0; 3],
                &[],
                Gravity::ZERO,
                0.001,
                origin,
                flight_controls(1e-7),
                cx,
            )
            .unwrap();
        assert!(rigid.is_empty());
        let external = Vec3::new(0.0004, 0.0, 0.0);
        same(report.external_impulse_kg_m_s, external, 1e-14);
        same(
            report.external_angular_impulse_kg_m2_s,
            Vec3::new(0.0, -0.0002, -0.00008),
            1e-14,
        );
        let after = total_momentum(&next, &[], origin, cx);
        same(
            after.0.sub(before.0),
            external.add(report.support_impulse_kg_m_s),
            2e-10,
        );
        same(
            after.1.sub(before.1),
            report
                .external_angular_impulse_kg_m2_s
                .add(report.support_angular_impulse_kg_m2_s),
            2e-10,
        );
        assert!(norm(report.support_impulse_kg_m_s) > 1e-6);
        near(
            report.source_work_j,
            0.4 * (next.displacement_m()[1][0] - u[1][0]),
            2e-12,
        );
        for &(dof, _) in &fixed {
            near(next.velocity_m_s()[dof / 3][dof % 3], 0.0, 1e-14);
        }
        flight_audit(&report);
    });
}

#[test]
fn a_nonrepresentable_native_kick_is_rejected_instead_of_credited_as_applied_work() {
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
        let original = solid.native_state().clone();
        let body =
            RigidBodyState::new(Pose::identity(), Vec3::new(1e20, 0.0, 0.0), Vec3::ZERO).unwrap();
        let result = solid.rigid_free_flight_candidate(
            &[body],
            &[properties()],
            &[],
            [0.0; 3],
            &[Wrench {
                force_world: Vec3::new(1.0, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            }],
            Gravity::ZERO,
            1e-4,
            [0.0; 3],
            flight_controls(1e-8),
            cx,
        );
        let Err(FreeFlightError::Audit(report)) = result else {
            panic!("lost native kick must fail physical balance")
        };
        near(report.linear_momentum_change_kg_m_s.x, 0.0, 1e-14);
        near(report.linear_momentum_residual_kg_m_s.x, -1e-4, 1e-14);
        assert!(report.source_work_j > 1e15);
        assert!(report.work_energy_residual_j < -1e15);
        assert_eq!(solid.native_state(), &original);
        assert_eq!(body.linear_momentum_world(), Vec3::new(1e20, 0.0, 0.0));
        let opposite = solid.rigid_free_flight_candidate(
            &[body; 2],
            &[properties(); 2],
            &[],
            [0.0; 3],
            &[1.0, -1.0].map(|x| Wrench {
                force_world: Vec3::new(x, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            }),
            Gravity::ZERO,
            1e-4,
            [0.0; 3],
            flight_controls(1e-8),
            cx,
        );
        let Err(FreeFlightError::Audit(report)) = opposite else {
            panic!("opposite lost kicks must not cancel the admission residual")
        };
        near(norm(report.linear_momentum_residual_kg_m_s), 0.0, 1e-14);
        near(report.work_energy_residual_j, 0.0, 1e-14);
        near(report.absolute_linear_residual_kg_m_s, 2e-4, 1e-14);
        assert!(report.absolute_work_residual_j > 9e15);
        assert_eq!(solid.native_state(), &original);
        gate.request();
        assert!(matches!(
            solid.rigid_free_flight_candidate(
                &[],
                &[],
                &[],
                [0.0; 3],
                &[],
                Gravity::ZERO,
                0.001,
                [0.0; 3],
                flight_controls(1e-8),
                cx,
            ),
            Err(FreeFlightError::Cancelled)
        ));
        assert_eq!(solid.native_state(), &original);
    });
}
