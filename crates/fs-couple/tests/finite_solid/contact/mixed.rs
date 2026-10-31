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

#[test]
fn full_tensor_mixed_event_matches_independent_world_velocity_operator() {
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
        let x = vec(solid.point_position_m(point, cx).unwrap());
        let orientations = [
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 1.0), 0.4).unwrap(),
            UnitQuaternion::from_axis_angle(Vec3::new(2.0, 1.0, -1.0), -0.7).unwrap(),
        ];
        let principals = [
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.6).unwrap(),
            UnitQuaternion::from_axis_angle(Vec3::new(-1.0, 1.0, 2.0), 0.3).unwrap(),
        ];
        let diagonals = [Vec3::new(0.3, 0.4, 0.5), Vec3::new(0.6, 0.8, 0.9)];
        let properties = [0, 1].map(|i| {
            MassProperties::from_principal_frame(2.0, Vec3::ZERO, diagonals[i], principals[i])
                .unwrap()
        });
        let arms = [
            Vec3::new(0.125, -0.25, 0.375),
            Vec3::new(-0.25, 0.125, 0.25),
            Vec3::new(0.2, 0.1, -0.125),
        ];
        let world_arms = [
            orientations[0].rotate_body_to_world(arms[0]),
            orientations[0].rotate_body_to_world(arms[1]),
            orientations[1].rotate_body_to_world(arms[2]),
        ];
        // A declared dyadic activation gap keeps rotation/reconstruction
        // roundoff from accidentally making this tensor fixture penetrating.
        // It does not change either point's velocity/mobility oracle.
        let gap = 1.0 / 1_099_511_627_776.0;
        let separation = Vec3::new(0.0, 0.0, gap);
        let positions = [
            x.sub(world_arms[0]).add(separation),
            x.sub(world_arms[0])
                .add(world_arms[1])
                .sub(world_arms[2])
                .add(separation.scale(2.0)),
        ];
        let axes = [
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        // Independent world-coordinate virtual impulse construction. No native
        // contact mobility or angular-velocity routine supplies this oracle.
        let inverse_world = |body: usize, torque: Vec3| {
            let principal = principals[body]
                .rotate_world_to_body(orientations[body].rotate_world_to_body(torque));
            orientations[body].rotate_body_to_world(principals[body].rotate_body_to_world(
                Vec3::new(
                    principal.x / diagonals[body].x,
                    principal.y / diagonals[body].y,
                    principal.z / diagonals[body].z,
                ),
            ))
        };
        let mut matrix = [0.0; 36];
        for column in 0..6 {
            let j = axes[column % 3];
            let (forces, moments, solid_velocity) = if column < 3 {
                (
                    [j, Vec3::ZERO],
                    [world_arms[0].cross(j), Vec3::ZERO],
                    j.scale(-3.5),
                )
            } else {
                (
                    [j.scale(-1.0), j],
                    [world_arms[1].cross(j.scale(-1.0)), world_arms[2].cross(j)],
                    Vec3::ZERO,
                )
            };
            let angular = [inverse_world(0, moments[0]), inverse_world(1, moments[1])];
            let wall_velocity = forces[0]
                .scale(0.5)
                .add(angular[0].cross(world_arms[0]))
                .sub(solid_velocity);
            let pair_velocity = forces[1]
                .scale(0.5)
                .add(angular[1].cross(world_arms[2]))
                .sub(forces[0].scale(0.5).add(angular[0].cross(world_arms[1])));
            for axis in 0..3 {
                matrix[axis * 6 + column] = axes[axis].dot(wall_velocity);
                matrix[(axis + 3) * 6 + column] = axes[axis].dot(pair_velocity);
            }
        }
        // Prescribe an interior sticking impulse and construct incoming point
        // velocities -WJ. The known answer requires no second contact solver.
        let expected = [1.0, 0.05, -0.03, 0.8, -0.04, 0.02];
        let incoming: [f64; 6] = std::array::from_fn(|row| {
            -(0..6)
                .map(|column| matrix[6 * row + column] * expected[column])
                .sum::<f64>()
        });
        assert!(incoming[0] < 0.0 && incoming[3] < 0.0);
        let expand = |offset: usize| {
            axes[0]
                .scale(incoming[offset])
                .add(axes[1].scale(incoming[offset + 1]))
                .add(axes[2].scale(incoming[offset + 2]))
        };
        let velocities = [expand(0), expand(0).add(expand(3))];
        let bodies = [0, 1].map(|i| {
            RigidBodyState::new(
                Pose::new(positions[i], orientations[i]).unwrap(),
                velocities[i].scale(2.0),
                Vec3::ZERO,
            )
            .unwrap()
        });
        let mut wall = witness(10, 0, point);
        wall.arm_body_m = arms[0];
        let pair = rigid_pair(
            20,
            RigidImpactPoint {
                body: 1,
                arm_body: arms[2],
            },
            RigidImpactPoint {
                body: 0,
                arm_body: arms[1],
            },
        );
        let contacts = [FiniteSolidMixedContact::from(wall), pair].map(|mut c| {
            c.activation_tolerance_m = 2.0 * gap;
            c.with_friction(EllipticFriction {
                coefficients_mu: [0.8, 0.4],
                principal_axis: [0.6, 0.8],
            })
        });
        let (_, rigid, report) = solid
            .mixed_elliptic_contact_candidate(
                &bodies,
                &properties,
                &contacts,
                [0.75, -1.25, 0.5],
                Default::default(),
                cx,
            )
            .unwrap();
        for row in 0..6 {
            near(
                report.contact.impulses[row / 3][row % 3],
                expected[row],
                2e-8,
            );
            near(
                report.contact.relative_velocities_after[row / 3][row % 3],
                0.0,
                2e-8,
            );
            for column in 0..6 {
                near(
                    report.contact_mobility_per_kg[6 * row + column],
                    matrix[6 * row + column],
                    2e-11,
                );
            }
        }
        let impulses = [
            Vec3::new(expected[1], expected[2], expected[0]),
            Vec3::new(expected[4], expected[5], expected[3]),
        ];
        same(
            rigid[0].linear_momentum_world(),
            bodies[0]
                .linear_momentum_world()
                .add(impulses[0])
                .sub(impulses[1]),
            2e-8,
        );
        same(
            rigid[1].linear_momentum_world(),
            bodies[1].linear_momentum_world().add(impulses[1]),
            2e-8,
        );
        same(
            rigid[0].angular_momentum_body(),
            arms[0]
                .cross(orientations[0].rotate_world_to_body(impulses[0]))
                .sub(arms[1].cross(orientations[0].rotate_world_to_body(impulses[1]))),
            2e-8,
        );
        same(
            rigid[1].angular_momentum_body(),
            arms[2].cross(orientations[1].rotate_world_to_body(impulses[1])),
            2e-8,
        );
        near(
            report.kinetic_change_j,
            0.5 * (0..6).map(|i| incoming[i] * expected[i]).sum::<f64>(),
            2e-8,
        );
        audit(&report);
    });
}
