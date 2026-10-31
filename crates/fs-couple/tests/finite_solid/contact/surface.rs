//! Current-geometry witnesses feeding the original angular-momentum response.
use super::*;
use fs_couple::finite_solid::contact::surface::*;

fn sphere(id: u64, rigid_body: usize, center_body_m: Vec3, radius_m: f64) -> FiniteSolidSphere {
    FiniteSolidSphere {
        id,
        rigid_body,
        center_body_m,
        radius_m,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}

#[test]
fn deformed_exterior_retains_original_material_weights_and_separate_witnesses() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        // Current x = (1.5 X, .75 Y, 1.25 Z) + (.25,-.125,.5).
        let displacement = [
            [0.25, -0.125, 0.5],
            [0.75, -0.125, 0.5],
            [0.25, -0.375, 0.5],
            [0.25, -0.125, 0.75],
        ];
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &displacement,
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let surface = FiniteSolidContactSurface::new(&solid, 4, cx).unwrap();
        let center = [0.625, 0.0625, 0.125];
        let offset = Vec3::new(0.125, 0.0, 0.0);
        let rigid = body_at(center, offset, Vec3::ZERO);
        let shape = sphere(11, 0, offset, 0.25);
        let controls = FiniteSolidSurfaceControls {
            activation_tolerance_m: 0.125,
            ..Default::default()
        };
        let query = surface
            .sphere_separations(&solid, &[rigid], &[shape], controls, cx)
            .unwrap();
        assert_eq!(query.contacts.len(), 1);
        let w = query.contacts[0];
        assert_eq!(w.contact.solid_point.tetrahedron(), 0);
        for (actual, expected) in w
            .contact
            .solid_point
            .weights()
            .into_iter()
            .zip([0.5, 0.25, 0.25, 0.0])
        {
            near(actual, expected, 2e-14);
        }
        same(w.solid_point_world_m, Vec3::new(0.625, 0.0625, 0.5), 2e-14);
        same(
            w.rigid_point_world_m,
            Vec3::new(0.625, 0.0625, 0.375),
            2e-14,
        );
        same(w.contact.normal_world, Vec3::new(0.0, 0.0, -1.0), 2e-14);
        same(w.contact.arm_body_m, Vec3::new(0.125, 0.0, 0.25), 2e-14);
        near(w.gap_m, 0.125, 2e-14);
        same(
            vec(solid.point_position_m(w.contact.solid_point, cx).unwrap()),
            w.solid_point_world_m,
            2e-14,
        );
        assert_ne!(w.rigid_point_world_m, w.solid_point_world_m);
        assert_eq!(query.admitted_contacts().unwrap(), vec![w.contact]);
    });
}

#[test]
fn shared_edges_and_vertices_deduplicate_without_losing_material_identity() {
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
        let surface = FiniteSolidContactSurface::new(&solid, 4, cx).unwrap();
        let rigid = [
            body_at([-0.125, 0.0, 0.0], Vec3::ZERO, Vec3::ZERO),
            body_at([0.25, -0.125, 0.0], Vec3::ZERO, Vec3::ZERO),
        ];
        let spheres = [
            sphere(9, 0, Vec3::ZERO, 0.125),
            sphere(3, 1, Vec3::ZERO, 0.125),
        ];
        let query = surface
            .sphere_separations(&solid, &rigid, &spheres, Default::default(), cx)
            .unwrap();
        assert_eq!(query.contacts.len(), 2);
        assert!(query.candidate_triangles >= 5);
        let edge_contact = query.contacts.iter().find(|w| w.sphere_id == 3).unwrap();
        let vertex_contact = query.contacts.iter().find(|w| w.sphere_id == 9).unwrap();
        assert_eq!(
            edge_contact.contact.solid_point.weights(),
            [0.75, 0.25, 0.0, 0.0]
        );
        assert_eq!(
            vertex_contact.contact.solid_point.weights(),
            [1.0, 0.0, 0.0, 0.0]
        );
        let permuted = surface
            .sphere_separations(
                &solid,
                &rigid,
                &[spheres[1], spheres[0]],
                Default::default(),
                cx,
            )
            .unwrap();
        assert_eq!(query.contacts, permuted.contacts);
        assert_eq!(query.closest, permuted.closest);
    });
}

#[test]
fn generated_offset_sphere_strike_drives_full_inertia_and_closes_world_angular_momentum() {
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
        let surface = FiniteSolidContactSurface::new(&solid, 4, cx).unwrap();
        let offset = Vec3::new(0.125, 0.0, 0.0);
        let rigid = body_at([0.25, 0.25, -0.125], offset, Vec3::new(0.0, 0.0, 2.0));
        let principal =
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), core::f64::consts::FRAC_PI_4)
                .unwrap();
        let property = MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(0.3, 0.4, 0.5),
            principal,
        )
        .unwrap();
        let contacts = surface
            .sphere_contacts(
                &solid,
                &[rigid],
                &[sphere(4, 0, offset, 0.125)],
                Default::default(),
                cx,
            )
            .unwrap();
        assert_eq!(contacts.len(), 1);
        same(contacts[0].arm_body_m, Vec3::new(0.125, 0.0, 0.125), 2e-14);
        let origin = [0.75, -1.25, 0.5];
        let (next, rigid_next, report) = solid
            .rigid_contact_candidate(
                &[rigid],
                &[property],
                &contacts,
                origin,
                Default::default(),
                cx,
            )
            .unwrap();
        // A one-kilogram P1 tet gives 20*(3/8-1/5)=3.5 at this point.
        // Rigid translation adds .5. The generated offset arm adds rotational
        // inverse mass .125^2*(1/.3+1/.4)/2 for the 45-degree principal frame.
        let mobility = 4.0 + 0.125_f64.powi(2) * (1.0 / 0.3 + 1.0 / 0.4) / 2.0;
        let j = 2.0 / mobility;
        near(report.contact.impulses[0][0], j, 2e-10);
        same(
            rigid_next[0].angular_momentum_body(),
            Vec3::new(0.0, 0.125 * j, 0.0),
            2e-10,
        );
        let before_solid = solid.momentum_about(origin, cx).unwrap();
        let after_solid = next.momentum_about(origin, cx).unwrap();
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
            total(rigid, before_solid.angular_kg_m2_s),
            total(rigid_next[0], after_solid.angular_kg_m2_s),
            3e-11,
        );
        audit(&report);
    });
}

#[test]
fn generated_rough_sphere_contact_converts_tangential_motion_into_physical_spin() {
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
        let surface = FiniteSolidContactSurface::new(&solid, 4, cx).unwrap();
        let rigid = body_at([0.25, 0.25, -0.125], Vec3::ZERO, Vec3::new(2.0, 0.0, 2.0));
        // Homogeneous 2 kg sphere: I = (2/5) m R^2 = .0125 kg m^2.
        let property =
            MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.0125, 0.0125, 0.0125)).unwrap();
        let mut shape = sphere(8, 0, Vec3::ZERO, 0.125);
        shape.friction_coefficient = 0.1;
        let contacts = surface
            .sphere_contacts(&solid, &[rigid], &[shape], Default::default(), cx)
            .unwrap();
        let (next, rigid_next, report) = solid
            .rigid_contact_candidate(
                &[rigid],
                &[property],
                &contacts,
                [0.0; 3],
                Default::default(),
                cx,
            )
            .unwrap();
        // Wn=4 gives Jn=.5. Coulomb sliding gives Jx=-mu*Jn=-.05;
        // Wt=4+R^2/I=5.25 leaves positive slip, so sticking is inadmissible.
        same(
            rigid_next[0]
                .linear_momentum_world()
                .sub(rigid.linear_momentum_world()),
            Vec3::new(-0.05, 0.0, -0.5),
            2e-10,
        );
        same(
            rigid_next[0].angular_momentum_body(),
            Vec3::new(0.0, -0.00625, 0.0),
            2e-10,
        );
        let angular_change = vec(next.momentum_about([0.0; 3], cx).unwrap().angular_kg_m2_s);
        same(angular_change, Vec3::new(0.125, -0.125, -0.0125), 2e-10);
        near(report.kinetic_change_j, -0.5934375, 2e-10);
        audit(&report);
    });
}

#[test]
fn penetration_inside_centres_and_exhausted_work_refuse_without_repair() {
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
        let surface = FiniteSolidContactSurface::new(&solid, 4, cx).unwrap();
        let shape = sphere(1, 0, Vec3::ZERO, 0.25);
        let overlapping = body_at([0.25, 0.25, -0.125], Vec3::ZERO, Vec3::ZERO);
        let query = surface
            .sphere_separations(&solid, &[overlapping], &[shape], Default::default(), cx)
            .unwrap();
        near(query.closest[0].gap_m, -0.125, 2e-14);
        assert!(matches!(
            query.admitted_contacts(),
            Err(FiniteSolidSurfaceError::Penetration { .. })
        ));
        assert!(matches!(
            surface.sphere_contacts(&solid, &[overlapping], &[shape], Default::default(), cx),
            Err(FiniteSolidSurfaceError::Penetration { .. })
        ));
        let inside = body_at([0.125, 0.125, 0.125], Vec3::ZERO, Vec3::ZERO);
        assert!(matches!(
            surface.sphere_separations(&solid, &[inside], &[shape], Default::default(), cx),
            Err(FiniteSolidSurfaceError::SphereCenterInside { .. })
        ));
        let touching = body_at([0.25, 0.25, -0.25], Vec3::ZERO, Vec3::ZERO);
        let no_work = FiniteSolidSurfaceControls {
            maximum_query_visits: 0,
            ..Default::default()
        };
        assert!(matches!(
            surface.sphere_contacts(&solid, &[touching], &[shape], no_work, cx),
            Err(FiniteSolidSurfaceError::WorkBudgetExhausted { .. })
                | Err(FiniteSolidSurfaceError::Query(
                    fs_query::TriangleQueryError::WorkBudgetExhausted { .. }
                ))
        ));
        let no_contacts = FiniteSolidSurfaceControls {
            maximum_contacts: 0,
            ..Default::default()
        };
        assert!(matches!(
            surface.sphere_contacts(&solid, &[touching], &[shape], no_contacts, cx),
            Err(FiniteSolidSurfaceError::ContactBudgetExhausted { .. })
        ));
        gate.request();
        assert!(matches!(
            surface.sphere_contacts(&solid, &[touching], &[shape], Default::default(), cx),
            Err(FiniteSolidSurfaceError::Cancelled)
        ));
    });
}
