//! Independent simultaneous chain, full-tensor, source and refusal witnesses.
use super::*;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_couple::finite_solid::contact::mixed::*;

fn rigid_pair(
    id: u64,
    first: RigidImpactPoint,
    second: RigidImpactPoint,
) -> FiniteSolidMixedContact {
    FiniteSolidMixedContact {
        id,
        point: first,
        counterpart: FiniteSolidContactCounterpart::Rigid(second),
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance_m: 0.0,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}

fn stack(
    point: FiniteSolidPoint,
    position: [f64; 3],
    velocity: Vec3,
) -> ([RigidBodyState; 2], [FiniteSolidMixedContact; 2]) {
    let bottom = Vec3::new(0.0, 0.0, -0.5);
    let bodies = [
        body_at(position, bottom, velocity),
        body_at(
            [position[0], position[1], position[2] + 1.0],
            bottom,
            velocity,
        ),
    ];
    let mut wall = witness(10, 0, point);
    wall.arm_body_m = bottom;
    let pair = rigid_pair(
        20,
        RigidImpactPoint {
            body: 1,
            arm_body: bottom,
        },
        RigidImpactPoint {
            body: 0,
            arm_body: bottom.scale(-1.0),
        },
    );
    (bodies, [wall.into(), pair])
}

#[test]
fn mixed_stack_resolves_shared_body_and_consistent_mass_in_one_solve() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        for supported in [false, true] {
            let fixed = if supported {
                vec![(0, 0.0), (1, 0.0), (2, 0.0)]
            } else {
                vec![]
            };
            let solid = FiniteSolidBody::new(
                problem(&material, &fixed, &[]),
                &[[0.0; 3]; 4],
                &[[0.0; 3]; 4],
                0.0,
                Default::default(),
                cx,
            )
            .unwrap();
            let point = FiniteSolidPoint::new(
                0,
                if supported {
                    [1.0, 0.0, 0.0, 0.0]
                } else {
                    [0.0, 0.25, 0.25, 0.5]
                },
            )
            .unwrap();
            let (bodies, contacts) = stack(
                point,
                solid.point_position_m(point, cx).unwrap(),
                Vec3::new(0.0, 0.0, -1.0),
            );
            let origin = [0.75, -0.25, 1.5];
            let (next, rigid, report) = solid
                .mixed_contact_candidate(
                    &bodies,
                    &[properties(); 2],
                    &contacts,
                    origin,
                    Default::default(),
                    cx,
                )
                .unwrap();
            // The normal block is [[s+1/2,-1/2],[-1/2,1]], with s=7/2
            // from the 1 kg consistent tet or zero at a stationary support.
            // W J=(1,0) => J_wall=1/(s+1/4), J_pair=J_wall/2.
            let s = if supported { 0.0 } else { 3.5 };
            let impulse = 1.0 / (s + 0.25);
            near(report.contact_mobility_per_kg[0], s + 0.5, 3e-12);
            near(report.contact_mobility_per_kg[3], -0.5, 3e-12);
            near(report.contact_mobility_per_kg[21], 1.0, 3e-12);
            near(report.contact.impulses[0][0], impulse, 3e-10);
            near(report.contact.impulses[1][0], impulse * 0.5, 3e-10);
            for body in &rigid {
                near(
                    body.center_of_mass_velocity_world(properties()).unwrap().z,
                    -1.0 + impulse * 0.25,
                    3e-10,
                );
                same(body.angular_momentum_body(), Vec3::ZERO, 3e-11);
            }
            near(
                next.point_velocity_m_s(point, cx).unwrap()[2],
                -s * impulse,
                3e-10,
            );
            near(report.kinetic_change_j, -0.5 * impulse, 3e-10);
            if supported {
                near(
                    report.solid_impulse.support.linear_kg_m_s[2],
                    impulse,
                    3e-10,
                );
                near(report.solid_impulse.impulse_work_j, 0.0, 0.0);
                assert_eq!(next.native_state(), solid.native_state());
            }
            assert_eq!(report.solid_mass_solves.len(), 3);
            audit(&report);
            let (permuted, rigid_permuted, receipt) = solid
                .mixed_contact_candidate(
                    &bodies,
                    &[properties(); 2],
                    &[contacts[1], contacts[0]],
                    origin,
                    Default::default(),
                    cx,
                )
                .unwrap();
            assert_eq!(permuted.native_state(), next.native_state());
            assert_eq!(rigid_permuted, rigid);
            assert_eq!(receipt.contact_ids, vec![10, 20]);
            assert_eq!(receipt.contact.impulses, report.contact.impulses);
        }
    });
}

#[test]
fn rigid_only_chain_has_no_solid_mass_columns_and_preserves_solid_state() {
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
        let property = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(0.5, 0.5, 0.5)).unwrap();
        let bodies: Vec<_> = [2.0, 1.0, 0.0]
            .into_iter()
            .enumerate()
            .map(|(i, z)| {
                RigidBodyState::new(
                    Pose::new(Vec3::new(0.0, 0.0, z), UnitQuaternion::IDENTITY).unwrap(),
                    Vec3::new(0.0, 0.0, if i == 0 { -1.0 } else { 0.0 }),
                    Vec3::ZERO,
                )
                .unwrap()
            })
            .collect();
        let pair = |id, a, b| {
            rigid_pair(
                id,
                RigidImpactPoint {
                    body: a,
                    arm_body: Vec3::new(0.0, 0.0, -0.5),
                },
                RigidImpactPoint {
                    body: b,
                    arm_body: Vec3::new(0.0, 0.0, 0.5),
                },
            )
        };
        let (next, rigid, report) = solid
            .mixed_contact_candidate(
                &bodies,
                &[property; 3],
                &[pair(1, 0, 1), pair(2, 1, 2)],
                [1.5, -0.5, 0.25],
                Default::default(),
                cx,
            )
            .unwrap();
        // Independent chain normal matrix [[2,-1],[-1,2]]. A pairwise
        // event cannot produce this simultaneous common outgoing velocity.
        near(report.contact_mobility_per_kg[0], 2.0, 1e-13);
        near(report.contact_mobility_per_kg[3], -1.0, 1e-13);
        near(report.contact_mobility_per_kg[21], 2.0, 1e-13);
        near(report.contact.impulses[0][0], 2.0 / 3.0, 3e-10);
        near(report.contact.impulses[1][0], 1.0 / 3.0, 3e-10);
        for body in rigid {
            near(body.linear_momentum_world().z, -1.0 / 3.0, 3e-10);
        }
        near(report.kinetic_change_j, -1.0 / 3.0, 3e-10);
        assert!(report.solid_mass_solves.is_empty());
        assert_eq!(next.native_state(), solid.native_state());
        audit(&report);
    });
}

#[test]
fn directional_mixed_stack_uses_physical_sliding_impulse_and_spin() {
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
        let (bodies, contacts) = stack(
            point,
            solid.point_position_m(point, cx).unwrap(),
            Vec3::new(2.0, 0.0, -1.0),
        );
        let directional = contacts.map(|c| {
            c.with_friction(EllipticFriction {
                coefficients_mu: [0.1, 0.4],
                principal_axis: [1.0, 0.0],
            })
        });
        let (next, rigid, report) = solid
            .mixed_elliptic_contact_candidate(
                &bodies,
                &[properties(); 2],
                &directional,
                [0.0; 3],
                Default::default(),
                cx,
            )
            .unwrap();
        let jn = 4.0 / 15.0;
        let jt = -2.0 / 75.0;
        near(report.contact.impulses[0][0], jn, 3e-10);
        near(report.contact.impulses[0][1], jt, 3e-10);
        near(report.contact.impulses[1][0], jn * 0.5, 3e-10);
        near(report.contact.impulses[1][1], 0.0, 3e-10);
        // The shared body's translational and rotational tangential cross
        // terms cancel: (-1/2) + (+1/2)=0, although its normal cross is -1/2.
        near(report.contact_mobility_per_kg[10], 0.0, 2e-13);
        near(rigid[0].angular_momentum_body().y, -jt * 0.5, 3e-10);
        near(rigid[1].angular_momentum_body().y, 0.0, 3e-10);
        near(
            report.kinetic_change_j,
            -0.5 * jn + (2.0 + 2.25 * jt) * jt,
            3e-10,
        );
        audit(&report);
        let (circular_solid, circular_rigid, circular) = solid
            .mixed_contact_candidate(
                &bodies,
                &[properties(); 2],
                &contacts.map(|c| c.with_friction(0.1)),
                [0.0; 3],
                Default::default(),
                cx,
            )
            .unwrap();
        for (a, b) in circular_solid
            .native_state()
            .v
            .iter()
            .zip(&next.native_state().v)
        {
            near(*a, *b, 3e-10);
        }
        for (a, b) in circular_rigid.iter().zip(&rigid) {
            same(a.linear_momentum_world(), b.linear_momentum_world(), 3e-10);
            same(a.angular_momentum_body(), b.angular_momentum_body(), 3e-10);
        }
        near(circular.kinetic_change_j, report.kinetic_change_j, 3e-10);
    });
}

#[test]
fn mixed_self_contacts_and_finite_gap_couples_refuse_without_publishing() {
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
        let bodies = [
            body_at([0.0, 0.0, 0.125], Vec3::ZERO, Vec3::new(2.0, 0.0, -1.0)),
            body_at([0.0; 3], Vec3::ZERO, Vec3::ZERO),
        ];
        let a = RigidImpactPoint {
            body: 0,
            arm_body: Vec3::ZERO,
        };
        let b = RigidImpactPoint {
            body: 1,
            arm_body: Vec3::ZERO,
        };
        assert!(matches!(
            solid.mixed_contact_candidate(
                &bodies,
                &[properties(); 2],
                &[rigid_pair(1, a, a)],
                [0.0; 3],
                Default::default(),
                cx
            ),
            Err(FiniteSolidMixedContactError::Invalid(_))
        ));
        let mut separated = rigid_pair(1, a, b);
        separated.activation_tolerance_m = 0.25;
        separated.friction_coefficient = 0.5;
        // Tangential action and reaction at distinct real witnesses produce
        // a nonzero couple. Activation tolerance does not authorize hiding it.
        let Err(FiniteSolidMixedContactError::Audit(report)) = solid.mixed_contact_candidate(
            &bodies,
            &[properties(); 2],
            &[separated],
            [0.0; 3],
            Default::default(),
            cx,
        ) else {
            panic!("a finite-gap friction couple must fail angular admission")
        };
        assert!(norm(report.angular_momentum_residual_kg_m2_s) > 0.01);
        assert_eq!(solid.native_state(), &original);
        gate.request();
        assert!(matches!(
            solid.mixed_contact_candidate(
                &bodies,
                &[properties(); 2],
                &[separated],
                [0.0; 3],
                Default::default(),
                cx
            ),
            Err(FiniteSolidMixedContactError::Cancelled)
        ));
        assert_eq!(solid.native_state(), &original);
    });
}
