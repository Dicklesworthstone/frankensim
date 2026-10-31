//! Material direction convection, physical anisotropic recoil and objectivity.
use super::*;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_couple::finite_solid::contact::mixed::*;
use fs_couple::finite_solid::contact::texture::*;

fn rigid(position: Vec3, orientation: UnitQuaternion, velocity: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, orientation).unwrap(),
        velocity.scale(2.0),
        Vec3::ZERO,
    )
    .unwrap()
}

#[test]
fn rigid_material_spin_rotates_anisotropic_friction_and_actual_angular_impulse() {
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
        let texture = FiniteSolidContactTexture::FirstRigid {
            direction_body: Vec3::new(7.0, 0.0, 0.0),
        };
        for (angle, coefficient) in [(0.0, 0.1), (core::f64::consts::FRAC_PI_2, 0.4)] {
            let orientation =
                UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), angle).unwrap();
            let body = rigid(
                Vec3::new(0.0, 0.0, 1.0),
                orientation,
                Vec3::new(4.0, 0.0, -1.0),
            );
            let mut geometry = witness(1, 0, point);
            geometry.arm_body_m = Vec3::new(0.0, 0.0, -1.0);
            let contact = FiniteSolidMixedContact::from(geometry).with_friction(EllipticFriction {
                coefficients_mu: [0.1, 0.4],
                principal_axis: [1.0, 0.0],
            });
            let frame = solid
                .contact_texture_frame(&[body], &contact, texture, cx)
                .unwrap();
            same(
                vec(frame.frame.longitudinal()),
                orientation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
                3e-14,
            );
            near(frame.tangential_fraction, 1.0, 3e-14);
            let contact = frame.bind(contact);
            let (_, next, report) = solid
                .mixed_elliptic_contact_candidate(
                    &[body],
                    &[properties()],
                    &[contact],
                    [0.0; 3],
                    Default::default(),
                    cx,
                )
                .unwrap();
            // Wnn=1/2; Wtt=Wbb=1/2+1/(1/2)=5/2. Normal impulse is 2.
            // Physical +x slip samples the longitudinal coefficient before a
            // quarter-turn, and the lateral coefficient after material spin.
            let jt = -2.0 * coefficient;
            same(
                next[0]
                    .linear_momentum_world()
                    .sub(body.linear_momentum_world()),
                Vec3::new(jt, 0.0, 2.0),
                3e-8,
            );
            same(
                orientation.rotate_body_to_world(next[0].angular_momentum_body()),
                Vec3::new(0.0, -jt, 0.0),
                3e-8,
            );
            near(report.kinetic_change_j, -1.0 + (4.0 + 1.25 * jt) * jt, 3e-8);
            audit(&report);
            assert_eq!(contact.friction_coefficient.coefficients_mu, [0.1, 0.4]);
        }
    });
}

#[test]
fn material_frames_follow_full_world_rotations_and_select_the_declared_rigid_side() {
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
        let point = FiniteSolidPoint::new(0, [1.0, 0.0, 0.0, 0.0]).unwrap();
        let angle = 0.37;
        let spin = Vec3::new(0.0, 0.0, angle);
        let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -3.0), 1.1).unwrap();
        let original_orientation = UnitQuaternion::IDENTITY.right_exp(spin).unwrap();
        let rotated_orientation = rotation.right_exp(spin).unwrap();
        let first = rigid(Vec3::ZERO, original_orientation, Vec3::ZERO);
        let second = rigid(Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::ZERO);
        let mut contact = FiniteSolidMixedContact::from(witness(3, 0, point));
        contact.counterpart = FiniteSolidContactCounterpart::Rigid(RigidImpactPoint {
            body: 1,
            arm_body: Vec3::ZERO,
        });
        let reference = Vec3::new(1.0, 0.25, 0.5);
        let first_texture = FiniteSolidContactTexture::FirstRigid {
            direction_body: reference,
        };
        let old = solid
            .contact_texture_frame(&[first, second], &contact, first_texture, cx)
            .unwrap();
        let counterpart = solid
            .contact_texture_frame(
                &[first, second],
                &contact,
                FiniteSolidContactTexture::SecondRigid {
                    direction_body: reference,
                },
                cx,
            )
            .unwrap();
        same(counterpart.material_direction_world, unit(reference), 3e-14);
        assert!(
            norm(
                counterpart
                    .material_direction_world
                    .sub(old.material_direction_world)
            ) > 0.1
        );

        let rotated = [
            rigid(Vec3::ZERO, rotated_orientation, Vec3::ZERO),
            rigid(Vec3::ZERO, rotation, Vec3::ZERO),
        ];
        contact.normal_world = rotation.rotate_body_to_world(contact.normal_world);
        let new = solid
            .contact_texture_frame(&rotated, &contact, first_texture, cx)
            .unwrap();
        same(
            new.material_direction_world,
            rotation.rotate_body_to_world(old.material_direction_world),
            3e-14,
        );
        same(
            vec(new.frame.normal()),
            rotation.rotate_body_to_world(vec(old.frame.normal())),
            3e-14,
        );
        same(
            vec(new.frame.longitudinal()),
            rotation.rotate_body_to_world(vec(old.frame.longitudinal())),
            3e-14,
        );
        same(
            vec(new.frame.lateral()),
            rotation.rotate_body_to_world(vec(old.frame.lateral())),
            3e-14,
        );
        near(new.tangential_fraction, old.tangential_fraction, 3e-14);
    });
}

#[test]
fn solid_texture_uses_total_affine_shear_fixed_dofs_and_superposed_world_motion() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let point = FiniteSolidPoint::new(0, [0.0, 0.25, 0.25, 0.5]).unwrap();
        let reference = Vec3::new(0.0, 1.0, 1.0);
        let expected_line = unit(Vec3::new(0.5, 1.0, 1.0));
        let expected_tangent = unit(Vec3::new(0.5, 1.0, 0.0));
        let rotation = UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 3.0), 0.8).unwrap();
        for (world_rotation, translation) in [
            (UnitQuaternion::IDENTITY, Vec3::ZERO),
            (
                UnitQuaternion::IDENTITY,
                Vec3::new(1_073_741_824.0, -536_870_912.0, 268_435_456.0),
            ),
            (rotation, Vec3::new(0.5, -0.25, 0.75)),
        ] {
            let displacement = NODES.map(|x| {
                let current = world_rotation
                    .rotate_body_to_world(Vec3::new(x[0] + 0.5 * x[1], x[1], x[2]))
                    .add(translation);
                [current.x - x[0], current.y - x[1], current.z - x[2]]
            });
            // A nonzero prescribed DOF supplies the actual shear contribution.
            let fixed = [(6, displacement[2][0])];
            let solid = FiniteSolidBody::new(
                problem(&material, &fixed, &[]),
                &displacement,
                &[[0.0; 3]; 4],
                0.0,
                Default::default(),
                cx,
            )
            .unwrap();
            let original = solid.native_state().clone();
            let body = rigid(
                vec(solid.point_position_m(point, cx).unwrap()),
                UnitQuaternion::IDENTITY,
                Vec3::ZERO,
            );
            let mut contact = FiniteSolidMixedContact::from(witness(1, 0, point));
            contact.normal_world = world_rotation.rotate_body_to_world(contact.normal_world);
            let texture = FiniteSolidContactTexture::Solid {
                direction_reference: reference,
            };
            let frame = solid
                .contact_texture_frame(&[body], &contact, texture, cx)
                .unwrap();
            // F=[[1,.5,0],[0,1,0],[0,0,1]], so F(0,1,1)=(.5,1,1).
            // Projection is (.5,1,0), rather than a polar-rotation-only line.
            same(
                frame.material_direction_world,
                world_rotation.rotate_body_to_world(expected_line),
                3e-13,
            );
            same(
                vec(frame.frame.longitudinal()),
                world_rotation.rotate_body_to_world(expected_tangent),
                3e-13,
            );
            near(frame.tangential_fraction, 5.0_f64.sqrt() / 3.0, 3e-13);
            // P1 carries the same reference material line throughout its tet.
            contact.counterpart = FiniteSolidContactCounterpart::Solid(
                FiniteSolidPoint::new(0, [1.0, 0.0, 0.0, 0.0]).unwrap(),
            );
            assert_eq!(
                solid
                    .contact_texture_frame(&[body], &contact, texture, cx)
                    .unwrap(),
                frame
            );
            assert_eq!(solid.native_state(), &original);
        }
    });
}

#[test]
fn ambiguous_or_unresolved_material_attachments_refuse_without_state_changes() {
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
        let point = FiniteSolidPoint::new(0, [1.0, 0.0, 0.0, 0.0]).unwrap();
        let body = rigid(Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::ZERO);
        let mut contact = FiniteSolidMixedContact::from(witness(1, 0, point));
        for direction in [
            Vec3::ZERO,
            Vec3::new(f64::NAN, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1e-16, 0.0, 1.0),
        ] {
            assert!(
                solid
                    .contact_texture_frame(
                        &[body],
                        &contact,
                        FiniteSolidContactTexture::FirstRigid {
                            direction_body: direction
                        },
                        cx
                    )
                    .is_err()
            );
        }
        assert!(
            solid
                .contact_texture_frame(
                    &[body],
                    &contact,
                    FiniteSolidContactTexture::SecondRigid {
                        direction_body: Vec3::new(1.0, 0.0, 0.0)
                    },
                    cx
                )
                .is_err()
        );
        contact.counterpart = FiniteSolidContactCounterpart::Rigid(contact.point);
        assert!(
            solid
                .contact_texture_frame(
                    &[body],
                    &contact,
                    FiniteSolidContactTexture::Solid {
                        direction_reference: Vec3::new(1.0, 0.0, 0.0)
                    },
                    cx
                )
                .is_err()
        );
        contact.counterpart = FiniteSolidContactCounterpart::Solid(
            FiniteSolidPoint::new(7, [1.0, 0.0, 0.0, 0.0]).unwrap(),
        );
        assert!(matches!(
            solid.contact_texture_frame(
                &[body],
                &contact,
                FiniteSolidContactTexture::Solid {
                    direction_reference: Vec3::new(1.0, 0.0, 0.0)
                },
                cx
            ),
            Err(FiniteSolidContactTextureError::Solid(_))
        ));
        gate.request();
        assert!(matches!(
            solid.contact_texture_frame(
                &[body],
                &contact,
                FiniteSolidContactTexture::FirstRigid {
                    direction_body: Vec3::new(1.0, 0.0, 0.0)
                },
                cx
            ),
            Err(FiniteSolidContactTextureError::Cancelled)
        ));
        assert_eq!(solid.native_state(), &original);
    });
}
