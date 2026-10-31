//! Automatic mixed geometry, simultaneous chain response and native path motion.
use super::*;
use fs_couple::finite_solid::contact::scene::{
    FiniteSolidSphereScene, FiniteSolidSphereSceneError, FiniteSolidSphereSceneKey, SpherePairLaw,
};

fn pair(first: u64, second: u64, restitution: f64) -> SpherePairLaw {
    SpherePairLaw {
        first,
        second,
        friction_coefficient: 0.0,
        restitution,
    }
}
fn sphere_properties() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.0125, 0.0125, 0.0125)).unwrap()
}

#[test]
fn previously_unsupported_rigid_pair_advances_through_detected_elastic_collision() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = [
            body_at([-0.625, -2.0, -2.0], Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)),
            body_at([0.625, -2.0, -2.0], Vec3::ZERO, Vec3::new(-1.0, 0.0, 0.0)),
        ];
        // Geometry input order and sphere identities differ from native addressing.
        let scene = FiniteSolidSphereScene::new(
            &[sphere(9, 0, 0.0, 0.0), sphere(3, 1, 0.0, 0.0)],
            &[pair(3, 9, 1.0)],
            1,
        )
        .unwrap();
        let (next, bodies, report) = solid
            .mixed_sphere_trajectory_candidate(
                &surface,
                &scene,
                &initial,
                &[sphere_properties(); 2],
                &[],
                [0.0; 3],
                &[Wrench::ZERO; 2],
                Gravity::ZERO,
                1.0,
                [0.75, -0.5, 0.25],
                settings(),
                cx,
            )
            .unwrap();
        near(next.native_state().t, 1.0, 2e-15);
        near(bodies[0].pose().position_world().x, -0.625, 2.1e-8);
        near(bodies[1].pose().position_world().x, 0.625, 2.1e-8);
        same(
            bodies[0].linear_momentum_world(),
            Vec3::new(-2.0, 0.0, 0.0),
            2e-10,
        );
        same(
            bodies[1].linear_momentum_world(),
            Vec3::new(2.0, 0.0, 0.0),
            2e-10,
        );
        let events: Vec<_> = report
            .steps
            .iter()
            .filter(|s| {
                s.impact
                    .as_ref()
                    .is_some_and(|i| i.contact.impulses.iter().any(|j| j[0] > 1e-8))
            })
            .collect();
        assert_eq!(events.len(), 1);
        assert_eq!(
            events[0].contact_keys,
            vec![FiniteSolidSphereSceneKey::Spheres {
                first: 3,
                second: 9
            }]
        );
        near(
            events[0].impact.as_ref().unwrap().contact.impulses[0][0],
            4.0,
            2e-10,
        );
        same(report.support_impulse_kg_m_s, Vec3::ZERO, 2e-10);
        same(report.linear_momentum_residual_kg_m_s, Vec3::ZERO, 2e-10);
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
        near(report.unallocated_contact_loss_j, 0.0, 2e-10);
        assert!(report.rejected_steps > 0);
        assert!(
            report
                .steps
                .iter()
                .all(|s| s.path.sphere_pairs_checked == 1)
        );
    });
}

#[test]
fn initially_resting_solid_contact_responds_to_a_simultaneous_rigid_neighbor() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = [
            // Positive activation gaps retain strict path admission despite
            // representable solver residuals; no negative gap is clamped.
            body_at([0.25, 0.25, -0.1250000001], Vec3::ZERO, Vec3::ZERO),
            body_at(
                [0.25, 0.25, -0.3750000002],
                Vec3::ZERO,
                Vec3::new(0.0, 0.0, 1.0),
            ),
        ];
        let spheres = [sphere(11, 0, 0.0, 0.0), sphere(22, 1, 0.0, 0.0)];
        let scene = FiniteSolidSphereScene::new(&spheres, &[pair(11, 22, 0.0)], 1).unwrap();
        let query = scene
            .contacts(&surface, &solid, &initial, settings().surface, cx)
            .unwrap();
        assert_eq!(query.contacts.len(), 2);
        assert_eq!(query.contacts[0].id % 2, 1);
        assert_eq!(query.contacts[1].id % 2, 0);
        // Reordering input leaves identity and the physical incidence unchanged.
        let reordered =
            FiniteSolidSphereScene::new(&[spheres[1], spheres[0]], &[pair(11, 22, 0.0)], 1)
                .unwrap()
                .contacts(&surface, &solid, &initial, settings().surface, cx)
                .unwrap();
        assert_eq!(query.contacts, reordered.contacts);
        assert_eq!(query.keys, reordered.keys);
        let (next, bodies, report) = solid
            .mixed_sphere_trajectory_candidate(
                &surface,
                &scene,
                &initial,
                &[sphere_properties(); 2],
                &[],
                [0.0; 3],
                &[Wrench::ZERO; 2],
                Gravity::ZERO,
                0.01,
                [0.5, -0.125, 0.3],
                settings(),
                cx,
            )
            .unwrap();
        let event = report.steps[0].impact.as_ref().unwrap();
        assert_eq!(event.contact.impulses.len(), 2);
        for impulse in &event.contact.impulses {
            near(impulse[0], 2.0, 2e-10);
        }
        for body in &bodies {
            same(body.linear_momentum_world(), Vec3::ZERO, 2e-10);
        }
        near(next.native_state().t, 0.01, 2e-16);
        near(report.unallocated_contact_loss_j, 1.0, 2e-10);
        same(
            report.support_impulse_kg_m_s,
            Vec3::new(0.0, 0.0, -2.0),
            2e-10,
        );
        same(
            report.support_angular_impulse_kg_m2_s,
            Vec3::new(-0.75, -0.5, 0.0),
            2e-10,
        );
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
        near(report.work_energy_residual_j, 0.0, 2e-10);
        assert_eq!(solid.native_state().steps, 0);
        assert_eq!(solid.native_state().q, next.native_state().q);
    });
}

#[test]
fn explicit_pair_law_coverage_and_nonpenetrating_geometry_are_required() {
    with_cx(|_, cx| {
        let spheres = [sphere(1, 0, 0.0, 0.1), sphere(2, 1, 0.0, 0.7)];
        assert!(FiniteSolidSphereScene::new(&spheres, &[], 1).is_err());
        assert!(FiniteSolidSphereScene::new(&spheres, &[pair(2, 1, 0.0)], 1).is_err());
        assert!(FiniteSolidSphereScene::new(&spheres, &[pair(1, 2, 0.0)], 0).is_err());
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let bodies = [
            body_at([0.0, -2.0, -2.0], Vec3::ZERO, Vec3::ZERO),
            body_at([0.25_f64.next_down(), -2.0, -2.0], Vec3::ZERO, Vec3::ZERO),
        ];
        let scene = FiniteSolidSphereScene::new(&spheres, &[pair(1, 2, 0.0)], 1).unwrap();
        assert!(
            matches!(scene.contacts(&surface, &solid, &bodies, settings().surface, cx),
            Err(FiniteSolidSphereSceneError::Penetration { first: 1, second: 2, gap_m }) if gap_m < 0.0)
        );
        assert_eq!(solid.native_state().steps, 0);
    });
}
