//! A native graded-memory body and unequal spheres share one detected impact.
use super::*;
use fs_couple::inelastic_impact::surface::{
    FiniteSolidContactSurface, FiniteSolidSphere, FiniteSolidSurfaceControls, scene::*,
};

pub(super) fn check_visco_scene(motion: &FiniteViscoTetMotion<'_>, cx: &Cx<'_>) {
    let point = TetImpulsePoint {
        tetrahedron: 0,
        weights: [0.5, 0.25, 0.25, 0.0],
    };
    let trace = motion.point_kinematics(point, cx).unwrap();
    let surface = FiniteSolidContactSurface::for_visco(motion, 4, cx).unwrap();
    let masses = [2.0, 3.0];
    // Uniform spheres of radius 1/8 have I=(2/5)*m*r²=m/160. Mechanics still
    // receives each complete native inertia; no scalar contact proxy is used.
    let properties = masses.map(|m| {
        MassProperties::new(m, Vec3::ZERO, Vec3::new(m / 160.0, m / 160.0, m / 160.0)).unwrap()
    });
    let bodies = [0, 1].map(|i| {
        RigidBodyState::new(
            Pose::new(
                vector(trace.position_m).sub(Vec3::new(0.0, 0.0, 0.125 + 0.25 * i as f64)),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::new(0.2 * masses[i], 0.0, (i + 1) as f64 * masses[i]),
            Vec3::ZERO,
        )
        .unwrap()
    });
    let spheres = [0, 1].map(|i| FiniteSolidSphere {
        id: [11, 22][i],
        rigid_body: i,
        center_body_m: Vec3::ZERO,
        radius_m: 0.125,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.25 },
    });
    let pair_law = SpherePairLaw {
        first: 11,
        second: 22,
        friction_coefficient: 0.0,
        restitution: 0.5,
    };
    let scene = FiniteSolidSphereScene::new(&spheres, &[pair_law], 1).unwrap();
    let controls = FiniteSolidSurfaceControls::default();
    let query = scene
        .visco_contacts(&surface, motion, &bodies, controls, cx)
        .unwrap();
    assert_eq!(query.contacts.len(), 2);
    assert_eq!(query.sphere_pair_candidates, 1);
    assert_eq!(
        query.keys[0],
        FiniteSolidSphereSceneKey::Spheres {
            first: 11,
            second: 22
        }
    );
    let original = query.surface.contacts[0];
    assert_eq!(
        query.keys[1],
        FiniteSolidSphereSceneKey::Solid {
            sphere: 11,
            contact_id: original.contact.id,
            face: original.face,
        }
    );
    assert_eq!(query.contacts[0].id, 1);
    assert_eq!(query.contacts[1].id, 2 * original.contact.id);
    assert_eq!(query.contacts[1].point.body, original.contact.rigid_body);
    assert_eq!(
        query.contacts[1].point.arm_body,
        original.contact.arm_body_m
    );
    assert_eq!(
        query.contacts[1].normal_world,
        original.contact.normal_world
    );
    assert_eq!(
        query.contacts[1].tangent_world,
        original.contact.tangent_world
    );
    assert_eq!(query.contacts[1].normal_law, spheres[0].normal_law);
    assert_eq!(
        query.contacts[0].normal_law,
        NormalVelocityLaw::NewtonRestitution {
            coefficient: pair_law.restitution
        }
    );
    let FiniteSolidContactCounterpart::Solid(mapped) = query.contacts[1].counterpart else {
        panic!("scene lost the original solid point")
    };
    assert_eq!(
        mapped.tetrahedron,
        original.contact.solid_point.tetrahedron()
    );
    assert_eq!(mapped.weights, original.contact.solid_point.weights());
    for (actual, expected) in mapped.weights.into_iter().zip(point.weights) {
        near(actual, expected, 1e-13);
    }
    let reordered = FiniteSolidSphereScene::new(&[spheres[1], spheres[0]], &[pair_law], 1)
        .unwrap()
        .visco_contacts(&surface, motion, &bodies, controls, cx)
        .unwrap();
    assert_eq!(query.keys, reordered.keys);
    assert_eq!(query.contacts, reordered.contacts);

    let event = visco_mixed_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &query.contacts,
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    // For rho=3+12x, this face trace has bᵀ M⁻¹ b=293/71. In (solid,pair)
    // order the normal matrix is [[293/71+1/2,-1/2],[-1/2,1/2+1/3]].
    // Incoming speeds are both -1; independent targets 1/4 and 1/2 give:
    let js = 3053.0 / 6144.0;
    let jp = 4297.0 / 2048.0;
    let loss = 11647.0 / 16384.0;
    near(event.report.contact_mobility_per_kg[0], 5.0 / 6.0, 2e-11);
    near(event.report.contact_mobility_per_kg[3], -0.5, 2e-11);
    near(
        event.report.contact_mobility_per_kg[21],
        293.0 / 71.0 + 0.5,
        2e-11,
    );
    assert_eq!(event.report.solid_mass_solves.len(), 3);
    near(event.report.contact.impulses[0][0], jp, 2e-10);
    near(event.report.contact.impulses[1][0], js, 2e-10);
    near(
        event.report.contact.relative_velocities_after[0][0],
        0.5,
        2e-10,
    );
    near(
        event.report.contact.relative_velocities_after[1][0],
        0.25,
        2e-10,
    );
    near(
        event.rigid_bodies[0].linear_momentum_world().z,
        2.0 - js + jp,
        2e-10,
    );
    near(
        event.rigid_bodies[1].linear_momentum_world().z,
        6.0 - jp,
        2e-10,
    );
    near(
        event
            .solid
            .point_kinematics(point, cx)
            .unwrap()
            .velocity_mps[2],
        293.0 / 71.0 * js,
        2e-10,
    );
    near(event.report.contact.dissipated_energy, loss, 2e-10);
    near(event.report.kinetic_change_j, -loss, 2e-10);
    mixed_audit(&event.report);
    assert_eq!(
        event.solid.solid().element_states(),
        motion.solid().element_states()
    );
    assert_eq!(
        event.solid.solid().displacement_m(),
        motion.solid().displacement_m()
    );
    assert_eq!(event.solid.solid().time_s(), motion.solid().time_s());
    assert_eq!(
        event.solid.reference_nodal_density_kg_m3(),
        Some(RHO.as_slice())
    );
    assert_eq!(event.rigid_bodies[0].pose(), bodies[0].pose());
    assert_eq!(event.rigid_bodies[1].pose(), bodies[1].pose());

    // Native geometry keeps the original refusals. A neighboring sphere is
    // never silently omitted or assigned the surface's different impact law.
    assert!(matches!(
        FiniteSolidSphereScene::new(&spheres, &[], 1),
        Err(FiniteSolidSphereSceneError::Invalid(_))
    ));
    assert!(matches!(
        FiniteSolidSphereScene::new(
            &spheres,
            &[SpherePairLaw {
                restitution: 1.1,
                ..pair_law
            }],
            1
        ),
        Err(FiniteSolidSphereSceneError::Invalid(_))
    ));
    assert!(FiniteSolidSphereScene::new(&spheres, &[pair_law], 0).is_err());
    let limited = FiniteSolidSurfaceControls {
        maximum_contacts: 1,
        ..controls
    };
    assert!(matches!(
        scene.visco_contacts(&surface, motion, &bodies, limited, cx),
        Err(FiniteSolidSphereSceneError::ContactBudget {
            maximum_contacts: 1
        })
    ));
    let mut penetrating = bodies;
    penetrating[1] = RigidBodyState::new(
        Pose::new(
            bodies[1]
                .pose()
                .position_world()
                .add(Vec3::new(0.0, 0.0, 1e-6)),
            bodies[1].pose().orientation(),
        )
        .unwrap(),
        bodies[1].linear_momentum_world(),
        bodies[1].angular_momentum_body(),
    )
    .unwrap();
    assert!(
        matches!(scene.visco_contacts(&surface, motion, &penetrating, controls, cx),
        Err(FiniteSolidSphereSceneError::Penetration { first: 11, second: 22, gap_m }) if gap_m < 0.0)
    );
    let retry = scene
        .visco_contacts(&surface, motion, &bodies, controls, cx)
        .unwrap();
    assert_eq!(retry.contacts, query.contacts);
    assert_eq!(retry.keys, query.keys);
}
