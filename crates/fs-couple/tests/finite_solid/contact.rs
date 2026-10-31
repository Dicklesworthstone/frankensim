//! Actual simultaneous native rigid/solid Coulomb impacts and independent oracles.
#[path = "contact/elliptic.rs"]
mod elliptic;

use super::*;
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_couple::finite_solid::contact::*;
use fs_couple::finite_solid::impulse::FiniteSolidPoint;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn vec(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn unit(v: Vec3) -> Vec3 {
    v.scale(norm(v).recip())
}
fn same(a: Vec3, b: Vec3, tolerance: f64) {
    near(norm(a.sub(b)), 0.0, tolerance);
}
fn properties() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.5, 0.5, 0.5)).unwrap()
}
fn body_at(point: [f64; 3], arm: Vec3, velocity: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(vec(point).sub(arm), UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(2.0),
        Vec3::ZERO,
    )
    .unwrap()
}
fn witness(id: u64, rigid_body: usize, solid_point: FiniteSolidPoint) -> FiniteSolidRigidContact {
    FiniteSolidRigidContact {
        id,
        rigid_body,
        arm_body_m: Vec3::ZERO,
        solid_point,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance_m: 0.0,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}
fn audit<R>(report: &FiniteSolidRigidContactReport<R>) {
    near(norm(report.linear_momentum_residual_kg_m_s), 0.0, 2e-11);
    near(norm(report.angular_momentum_residual_kg_m2_s), 0.0, 2e-11);
    near(report.work_energy_residual_j, 0.0, 2e-11);
    near(report.dissipation_residual_j, 0.0, 2e-11);
    near(report.endpoint_velocity_residual_m_s, 0.0, 2e-10);
    assert!(report.kinetic_change_j <= 1e-12);
}

#[test]
fn source_mass_contact_resolves_stick_slide_and_restitution_on_native_bodies() {
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
        let point = FiniteSolidPoint::new(0, [0.0, 0.25, 0.25, 0.5]).unwrap();
        let x = solid.point_position_m(point, cx).unwrap();
        let n = unit(Vec3::new(1.0, 1.0, 1.0));
        let t = unit(Vec3::new(1.0, -1.0, 0.0));
        // M^-1 = 20(I-11^T/5) for this 1 kg reference tet.
        // The point's inverse mass is 20*(3/8-1/5)=3.5; rigid adds 0.5.
        let mobility = 4.0;
        for (vt, mu, restitution) in [(0.2, 0.5, 0.0), (4.0, 0.1, 0.0), (0.0, 0.0, 0.5)] {
            let rigid = body_at(x, Vec3::ZERO, n.scale(-2.0).add(t.scale(vt)));
            let mut c = witness(17, 0, point);
            c.normal_world = n;
            c.tangent_world = t;
            c.friction_coefficient = mu;
            c.normal_law = NormalVelocityLaw::NewtonRestitution {
                coefficient: restitution,
            };
            let (next, rigid_next, report) = solid
                .rigid_contact_candidate(
                    &[rigid],
                    &[properties()],
                    &[c],
                    [0.6, -0.3, 0.2],
                    Default::default(),
                    cx,
                )
                .unwrap();
            let jn = (1.0 + restitution) * 2.0 / mobility;
            let jt = -(vt / mobility).min(mu * jn);
            near(report.contact.impulses[0][0], jn, 2e-10);
            near(report.contact.impulses[0][1], jt, 2e-10);
            near(report.contact.impulses[0][2], 0.0, 2e-11);
            let expected_work =
                0.5 * (-2.0 + 2.0 * restitution) * jn + (vt + 0.5 * mobility * jt) * jt;
            near(report.kinetic_change_j, expected_work, 2e-10);
            for axis in 0..3 {
                near(
                    report.contact_mobility_per_kg[3 * axis + axis],
                    mobility,
                    3e-12,
                );
            }
            audit(&report);
            assert_eq!(solid.native_state().steps, next.native_state().steps);
            assert_eq!(solid.native_state().q, next.native_state().q);
            assert_eq!(solid.native_state().t, next.native_state().t);
            assert_eq!(rigid.pose(), rigid_next[0].pose());
            assert!(next.native_state().v.iter().any(|v| v.abs() > 0.0));
        }
    });
}

#[test]
fn two_rigids_striking_one_solid_retain_negative_cross_contact_mobility() {
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
        let a = FiniteSolidPoint::new(0, [0.0, 1.0, 0.0, 0.0]).unwrap();
        let b = FiniteSolidPoint::new(0, [0.0, 0.0, 1.0, 0.0]).unwrap();
        let rigid = [
            body_at(NODES[1], Vec3::ZERO, Vec3::new(0.0, 0.0, -1.0)),
            body_at(NODES[2], Vec3::ZERO, Vec3::new(0.0, 0.0, -1.0)),
        ];
        let contacts = [witness(20, 0, a), witness(10, 1, b)];
        let origin = [2.0, -3.0, 4.0];
        let (next, rigid_next, report) = solid
            .rigid_contact_candidate(
                &rigid,
                &[properties(); 2],
                &contacts,
                origin,
                Default::default(),
                cx,
            )
            .unwrap();
        // Normal block [[16.5,-4],[-4,16.5]] has (1,1) eigenvalue 12.5.
        // A diagonal or pairwise approximation gives a different impulse.
        assert_eq!(report.contact_ids, vec![10, 20]);
        near(report.contact_mobility_per_kg[0], 16.5, 2e-12);
        near(report.contact_mobility_per_kg[3], -4.0, 2e-12);
        for impulse in &report.contact.impulses {
            near(impulse[0], 0.08, 2e-10);
        }
        for body in &rigid_next {
            near(body.linear_momentum_world().z / 2.0, -0.96, 2e-10);
        }
        near(next.point_velocity_m_s(a, cx).unwrap()[2], -0.96, 2e-10);
        near(next.point_velocity_m_s(b, cx).unwrap()[2], -0.96, 2e-10);
        near(report.kinetic_change_j, -0.08, 3e-10);
        audit(&report);
        let (permuted, permuted_rigid, permuted_report) = solid
            .rigid_contact_candidate(
                &rigid,
                &[properties(); 2],
                &[contacts[1], contacts[0]],
                origin,
                Default::default(),
                cx,
            )
            .unwrap();
        assert_eq!(permuted.native_state(), next.native_state());
        assert_eq!(permuted_rigid, rigid_next);
        assert_eq!(permuted_report.contact.impulses, report.contact.impulses);
    });
}

#[test]
fn off_center_full_inertia_impact_conserves_world_spin_and_orbital_momentum() {
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
        let point = FiniteSolidPoint::new(0, [0.0, 0.25, 0.25, 0.5]).unwrap();
        let arm = Vec3::new(0.125, -0.25, 0.375);
        let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap();
        let property = MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(0.3, 0.4, 0.5),
            principal,
        )
        .unwrap();
        let rigid = body_at(
            solid.point_position_m(point, cx).unwrap(),
            arm,
            Vec3::new(0.0, 0.0, -2.0),
        );
        let mut c = witness(1, 0, point);
        c.arm_body_m = arm;
        let origin = [0.75, -1.25, 0.5];
        let (next, rigid_next, report) = solid
            .rigid_contact_candidate(&[rigid], &[property], &[c], origin, Default::default(), cx)
            .unwrap();
        let torque_principal = principal.rotate_world_to_body(arm.cross(Vec3::new(0.0, 0.0, 1.0)));
        let inv_mass = 4.0
            + torque_principal.x * torque_principal.x / 0.3
            + torque_principal.y * torque_principal.y / 0.4
            + torque_principal.z * torque_principal.z / 0.5;
        let j = 2.0 / inv_mass;
        near(report.contact.impulses[0][0], j, 2e-10);
        same(
            rigid_next[0].angular_momentum_body(),
            arm.cross(Vec3::new(0.0, 0.0, j)),
            2e-10,
        );
        let before = solid.momentum_about(origin, cx).unwrap();
        let after = next.momentum_about(origin, cx).unwrap();
        let total = |body: RigidBodyState, solid_l: [f64; 3]| {
            body.pose()
                .position_world()
                .sub(vec(origin))
                .cross(body.linear_momentum_world())
                .add(
                    body.pose()
                        .orientation()
                        .rotate_body_to_world(body.angular_momentum_body()),
                )
                .add(vec(solid_l))
        };
        same(
            total(rigid, before.angular_kg_m2_s),
            total(rigid_next[0], after.angular_kg_m2_s),
            3e-11,
        );
        audit(&report);
    });
}

#[test]
fn fixed_solid_supports_receive_impact_reaction_and_do_zero_work() {
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
        let point = FiniteSolidPoint::new(0, [1.0, 0.0, 0.0, 0.0]).unwrap();
        let rigid = body_at(NODES[0], Vec3::ZERO, Vec3::new(0.0, 0.0, -1.0));
        let (next, rigid_next, report) = solid
            .rigid_contact_candidate(
                &[rigid],
                &[properties()],
                &[witness(1, 0, point)],
                [0.3, -0.4, 0.2],
                Default::default(),
                cx,
            )
            .unwrap();
        near(report.contact.impulses[0][0], 2.0, 1e-10);
        same(rigid_next[0].linear_momentum_world(), Vec3::ZERO, 2e-10);
        near(report.solid_impulse.support.linear_kg_m_s[2], 2.0, 1e-10);
        near(report.solid_impulse.impulse_work_j, 0.0, 0.0);
        near(report.kinetic_change_j, -1.0, 1e-10);
        assert_eq!(next.native_state(), solid.native_state());
        audit(&report);
    });
}

#[test]
fn invalid_frames_gaps_and_late_angular_refusal_leave_both_inputs_unchanged() {
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
        let point = FiniteSolidPoint::new(0, [0.0, 0.25, 0.25, 0.5]).unwrap();
        let x = solid.point_position_m(point, cx).unwrap();
        let rigid = body_at(x, Vec3::ZERO, Vec3::new(1.0, 0.0, -1.0));
        let c = witness(1, 0, point);
        let mut collinear = c;
        collinear.normal_world = Vec3::new(1.0, 1.0, 1.0);
        collinear.tangent_world = collinear.normal_world;
        assert!(matches!(
            solid.rigid_contact_candidate(
                &[rigid],
                &[properties()],
                &[collinear],
                [0.0; 3],
                Default::default(),
                cx
            ),
            Err(FiniteSolidRigidContactError::Invalid(_))
        ));
        let penetrated = body_at(
            [x[0], x[1], x[2] - 0.125],
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, -1.0),
        );
        assert!(matches!(
            solid.rigid_contact_candidate(
                &[penetrated],
                &[properties()],
                &[c],
                [0.0; 3],
                Default::default(),
                cx
            ),
            Err(FiniteSolidRigidContactError::Contact(_))
        ));
        let mut gap_contact = c;
        gap_contact.activation_tolerance_m = 0.25;
        gap_contact.friction_coefficient = 0.5;
        let separated = body_at(
            [x[0], x[1], x[2] + 0.125],
            Vec3::ZERO,
            Vec3::new(1.0, 0.0, -1.0),
        );
        // Activation is no permission to move witnesses to a fictitious point.
        // Tangential impulse across the finite gap leaves a real angular couple.
        assert!(matches!(
            solid.rigid_contact_candidate(
                &[separated],
                &[properties()],
                &[gap_contact],
                [0.0; 3],
                Default::default(),
                cx
            ),
            Err(FiniteSolidRigidContactError::Audit(_))
        ));
        let mut capped = FiniteSolidRigidContactControls::default();
        capped.contact.maximum_contacts = 1;
        let mut second = c;
        second.id = 2;
        assert!(
            solid
                .rigid_contact_candidate(
                    &[rigid],
                    &[properties()],
                    &[c, second],
                    [0.0; 3],
                    capped,
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.native_state(), &original);
        gate.request();
        assert!(matches!(
            solid.rigid_contact_candidate(
                &[rigid],
                &[properties()],
                &[c],
                [0.0; 3],
                Default::default(),
                cx
            ),
            Err(FiniteSolidRigidContactError::Cancelled)
        ));
        assert_eq!(solid.native_state(), &original);
        same(
            rigid.linear_momentum_world(),
            Vec3::new(2.0, 0.0, -2.0),
            0.0,
        );
    });
}
