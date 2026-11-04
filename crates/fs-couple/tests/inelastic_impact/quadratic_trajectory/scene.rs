//! G1/G3 declared pair encounters with the original native material trajectory.

use super::*;
use fs_couple::inelastic_impact::surface::quadratic::scene::{
    FiniteSolidSphereScene, FiniteSolidSphereSceneKey, SpherePairLaw,
};

#[test]
fn unequal_pair_collision_localizes_and_retains_native_maxwell_history() {
    with_cx(|_, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = material();
        let motion = relaxed(&nodes, &elements, &law, cx);
        let original = motion.solid().quadrature_states().to_vec();
        let pressure = motion.solid().pressure_pa().to_vec();
        let held: Vec<_> = (0..3 * nodes.len()).collect();
        let surface = QuadraticSolidContactSurface::for_mixed_visco(&motion, 4, cx).unwrap();
        let properties = [2.0, 3.0].map(|mass| {
            MassProperties::new(mass, Vec3::ZERO, Vec3::new(0.01, 0.02, 0.03)).unwrap()
        });
        let bodies = [(-3.0, 2.0), (-1.0, -3.0)].map(|(x, p)| {
            RigidBodyState::new(
                Pose::new(Vec3::new(x, -3.0, -3.0), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::new(p, 0.0, 0.0),
                Vec3::ZERO,
            )
            .unwrap()
        });
        let spheres = [0, 1].map(|i| FiniteSolidSphere {
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 1.0 },
            ..sphere([7, 9][i], i)
        });
        let pair = SpherePairLaw {
            first: 7,
            second: 9,
            friction_coefficient: 0.0,
            restitution: 0.5,
        };
        let scene = FiniteSolidSphereScene::new(&spheres, &[pair], 1).unwrap();
        let mut c = controls();
        c.maximum_step_s = 1.5;
        let run = |scene: &FiniteSolidSphereScene, c| {
            mixed_visco_sphere_scene_trajectory_candidate(
                &motion,
                &surface,
                scene,
                &bodies,
                &properties,
                &held,
                &[],
                [0.0; 3],
                None,
                &[Wrench::ZERO; 2],
                Gravity::ZERO,
                1.5,
                ORIGIN,
                c,
                cx,
            )
        };
        let result = run(&scene, c).unwrap();
        // u=(1,-1), m=(2,3), pair e=1/2 give J=18/5, outgoing
        // momenta (-8/5,3/5), loss9/5 J. Both solid-interface laws have e=1.
        for (index, expected) in [-1.6, 0.6].into_iter().enumerate() {
            near(
                norm(
                    result.rigid_bodies[index]
                        .linear_momentum_world()
                        .sub(Vec3::new(expected, 0.0, 0.0)),
                ),
                0.0,
                2e-8,
            );
        }
        // Zero-gap time .95 s. Early activation retains its actual positive gap.
        for (index, x) in [-2.49, -1.84].into_iter().enumerate() {
            near(
                norm(
                    result.rigid_bodies[index]
                        .pose()
                        .position_world()
                        .sub(Vec3::new(x, -3.0, -3.0)),
                ),
                0.0,
                2e-5,
            );
        }
        near(result.report.unallocated_contact_loss_j, 1.8, 2e-8);
        near(result.report.kinetic_change_j, -1.8, 2e-8);
        near(norm(result.report.linear_momentum_change_kg_m_s), 0.0, 1e-8);
        near(norm(result.report.support_impulse_kg_m_s), 0.0, 1e-8);
        assert!(result.report.rejected_steps > 0);
        assert_balanced(&result.report);
        let mut replay = original.clone();
        let mut collisions = 0;
        for step in &result.report.steps {
            if let Some(impact) = &step.impact {
                assert_eq!(
                    step.path.contact_keys,
                    [FiniteSolidSphereSceneKey::Spheres {
                        first: 7,
                        second: 9
                    }]
                );
                assert_eq!(impact.contact_ids, [1]);
                if impact.contact.impulses[0][0] > 1e-8 {
                    collisions += 1;
                    near(impact.contact.impulses[0][0], 3.6, 2e-8);
                }
            }
            for state in &mut replay {
                *state = law
                    .step_remainder_response(state, &state.deformation_gradient, step.duration_s)
                    .unwrap()
                    .state;
            }
        }
        assert_eq!(collisions, 1);
        assert_eq!(result.solid.solid().quadrature_states(), replay);
        assert_eq!(
            result.solid.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        near(result.solid.solid().time_s(), 1.507, 4e-15);
        let reordered = FiniteSolidSphereScene::new(&[spheres[1], spheres[0]], &[pair], 1).unwrap();
        let replayed = run(&reordered, c).unwrap();
        assert_eq!(
            replayed.solid.solid().quadrature_states(),
            result.solid.solid().quadrature_states()
        );
        assert_eq!(
            replayed.solid.solid().pressure_pa(),
            result.solid.solid().pressure_pa()
        );
        for i in 0..2 {
            assert_eq!(
                replayed.rigid_bodies[i].pose(),
                result.rigid_bodies[i].pose()
            );
            assert_eq!(
                replayed.rigid_bodies[i].linear_momentum_world(),
                result.rigid_bodies[i].linear_momentum_world()
            );
        }
        let mut limited = c;
        limited.maximum_step_s = 0.6;
        limited.maximum_rejections = 0;
        assert!(matches!(run(&scene, limited),
            Err(MixedViscoSphereTrajectoryError::RejectionBudget { elapsed_s, last })
            if elapsed_s > 0.0 && matches!(*last, MixedViscoSpherePathRefusal::ScenePath(_))));
        assert_eq!(motion.solid().quadrature_states(), original);
        assert_eq!(motion.solid().pressure_pa(), pressure);
        assert_eq!(motion.solid().time_s(), 0.007);
    });
}

#[test]
fn single_curved_sphere_keeps_identical_native_motion_through_both_entries() {
    with_cx(|_, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = material();
        let motion = relaxed(&nodes, &elements, &law, cx);
        let held: Vec<_> = (0..3 * nodes.len()).collect();
        let c = controls();
        let ordinary = bounce(&motion, &held, 0.1, 0.3, c, cx).unwrap();
        let surface = QuadraticSolidContactSurface::for_mixed_visco(&motion, 4, cx).unwrap();
        let (point, normal) = curved_point_and_normal();
        let body = RigidBodyState::new(
            Pose::new(
                point.add(normal.scale(0.05 + 0.1)),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            normal.scale(-2.0),
            Vec3::ZERO,
        )
        .unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(7, 0)], &[], 0).unwrap();
        let complete = mixed_visco_sphere_scene_trajectory_candidate(
            &motion,
            &surface,
            &scene,
            &[body],
            &[properties()],
            &held,
            &[],
            [0.0; 3],
            None,
            &[Wrench::ZERO],
            Gravity::ZERO,
            0.3,
            ORIGIN,
            c,
            cx,
        )
        .unwrap();
        assert_eq!(
            complete.solid.solid().displacement_m(),
            ordinary.solid.solid().displacement_m()
        );
        assert_eq!(
            complete.solid.solid().pressure_pa(),
            ordinary.solid.solid().pressure_pa()
        );
        assert_eq!(
            complete.solid.solid().quadrature_states(),
            ordinary.solid.solid().quadrature_states()
        );
        assert_eq!(
            complete.rigid_bodies[0].pose(),
            ordinary.rigid_bodies[0].pose()
        );
        assert_eq!(
            complete.rigid_bodies[0].linear_momentum_world(),
            ordinary.rigid_bodies[0].linear_momentum_world()
        );
        assert_eq!(complete.report.steps.len(), ordinary.report.steps.len());
        assert_eq!(
            complete.report.rejected_steps,
            ordinary.report.rejected_steps
        );
        for (scene_step, raw_step) in complete.report.steps.iter().zip(&ordinary.report.steps) {
            assert_eq!(scene_step.duration_s, raw_step.duration_s);
            assert_eq!(scene_step.material, raw_step.material);
            assert_eq!(scene_step.path.contact_keys, raw_step.path.contact_keys);
            assert_eq!(scene_step.path.minimum_pair_clearance_m, None);
        }
        assert_balanced(&complete.report);
    });
}
