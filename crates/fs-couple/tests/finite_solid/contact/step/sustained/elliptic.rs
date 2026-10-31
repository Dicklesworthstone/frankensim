//! Original directional contact through the identical native force/drift path.
use super::*;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContactError,
};
use fs_couple::finite_solid::contact::step::sustained::FiniteSolidSustainedTrajectoryError;
use fs_couple::finite_solid::contact::step::sustained::elliptic::FiniteSolidSustainedEllipticFailure;

#[test]
fn equal_principal_coefficients_recover_circular_trajectory_with_full_inertia() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 1.0, 0.3)], &[], 0).unwrap();
        let property = MassProperties::from_inertia_tensor(
            2.0,
            Vec3::ZERO,
            fs_ga::Mat3 {
                m: [
                    0.016, 0.002, 0.001, 0.002, 0.014, -0.001, 0.001, -0.001, 0.012,
                ],
            },
        )
        .unwrap();
        let initial = RigidBodyState::new(
            Pose::new(
                Vec3::new(0.25, 0.25, -0.125 - GAP),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::new(0.8, -0.6, 0.0),
            Vec3::new(0.003, -0.004, 0.002),
        )
        .unwrap();
        let wrench = Wrench {
            force_world: Vec3::new(0.2, -0.1, 0.0),
            torque_body: Vec3::new(0.01, 0.02, -0.015),
        };
        let (c_solid, c_bodies, circular) = solid
            .sustained_sphere_trajectory_candidate(
                &surface,
                &scene,
                &[initial],
                &[property],
                &[],
                [0.0; 3],
                &[wrench],
                upward_gravity(),
                STEP,
                ORIGIN,
                sustained_controls(STEP),
                cx,
            )
            .unwrap();
        let mut calls = 0;
        let (e_solid, e_bodies, elliptic) = solid
            .sustained_elliptic_sphere_trajectory_candidate(
                &surface,
                &scene,
                &[initial],
                &[property],
                &[],
                [0.0; 3],
                &[wrench],
                upward_gravity(),
                STEP,
                ORIGIN,
                sustained_controls(STEP),
                |context, key, contact| {
                    calls += 1;
                    assert!(matches!(
                        key,
                        FiniteSolidSphereSceneKey::Solid { sphere: 9, .. }
                    ));
                    near(context.elapsed_s, 0.0, 0.0);
                    near(context.duration_s, STEP, 0.0);
                    // The real force kick has happened, but normal restitution
                    // still refers to the unforced resting incoming speed.
                    near(
                        context.bodies[0].linear_momentum_world().z,
                        2.0 * ACCELERATION * STEP,
                        2e-14,
                    );
                    assert_eq!(contact.normal_law, NormalVelocityLaw::Target { speed: 0.0 });
                    Ok::<_, core::convert::Infallible>(EllipticFriction {
                        coefficients_mu: [0.3; 2],
                        principal_axis: [0.3, 0.4],
                    })
                },
                cx,
            )
            .unwrap();
        assert_eq!(calls, 1);
        assert_eq!(c_solid.native_state(), e_solid.native_state());
        same(
            c_bodies[0].linear_momentum_world(),
            e_bodies[0].linear_momentum_world(),
            2e-11,
        );
        same(
            c_bodies[0].angular_momentum_body(),
            e_bodies[0].angular_momentum_body(),
            2e-11,
        );
        same(
            c_bodies[0].pose().position_world(),
            e_bodies[0].pose().position_world(),
            2e-13,
        );
        same(
            circular.support_impulse_kg_m_s,
            elliptic.support_impulse_kg_m_s,
            2e-11,
        );
        same(
            circular.support_angular_impulse_kg_m2_s,
            elliptic.support_angular_impulse_kg_m2_s,
            2e-11,
        );
        near(circular.source_work_j, elliptic.source_work_j, 2e-11);
        near(
            circular.unallocated_contact_loss_j,
            elliptic.unallocated_contact_loss_j,
            2e-11,
        );
        near(elliptic.work_energy_residual_j, 0.0, 2e-10);
        let disk = &circular.steps[0].impact.as_ref().unwrap().contact;
        let ellipse = &elliptic.steps[0].impact.as_ref().unwrap().contact;
        assert_eq!(disk.impulses, ellipse.impulses);
        assert_eq!(disk.impulses, ellipse.transformed.impulses);
        assert_eq!(
            ellipse.coordinates[0].tangent_impulse_map,
            [[1.0, 0.0], [0.0, 1.0]]
        );
        assert!(elliptic.steps[0].rigid[0].absolute_angular_residual_kg_m2_s < 1e-9);
    });
}

#[test]
fn turning_the_ellipse_changes_actual_sliding_impulse_and_spin_by_the_principal_law() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        // The scene's scalar friction is intentionally irrelevant: this entry
        // requires an explicit ellipse for each original current scene row.
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 1.0, 0.0)], &[], 0).unwrap();
        let initial = body_at(
            [0.25, 0.25, -0.125 - GAP],
            Vec3::ZERO,
            Vec3::new(1.0, 0.0, 0.0),
        );
        let mut impulses = Vec::new();
        // n=-z, t=-y and b=-x. Rotating the principal ellipse by 90 degrees
        // switches world-x sliding from mu=.8 to mu=.2.
        for (axis, mu) in [([1.0, 0.0], 0.8), ([0.0, 1.0], 0.2)] {
            let (_, bodies, report) = solid
                .sustained_elliptic_sphere_trajectory_candidate(
                    &surface,
                    &scene,
                    &[initial],
                    &[homogeneous_sphere()],
                    &[],
                    [0.0; 3],
                    &[Wrench::ZERO],
                    upward_gravity(),
                    STEP,
                    ORIGIN,
                    sustained_controls(STEP),
                    |_, _, _| {
                        Ok::<_, core::convert::Infallible>(EllipticFriction {
                            coefficients_mu: [0.2, 0.8],
                            principal_axis: axis,
                        })
                    },
                    cx,
                )
                .unwrap();
            let normal_impulse = 2.0 * ACCELERATION * STEP;
            let tangent_impulse = mu * normal_impulse;
            let impact = &report.steps[0].impact.as_ref().unwrap().contact;
            near(impact.impulses[0][0], normal_impulse, 2e-10);
            near(impact.impulses[0][1], 0.0, 2e-10);
            near(impact.impulses[0][2], tangent_impulse, 2e-10);
            same(
                bodies[0].linear_momentum_world(),
                Vec3::new(2.0 - tangent_impulse, 0.0, 0.0),
                2e-10,
            );
            let world_spin = bodies[0]
                .pose()
                .orientation()
                .rotate_body_to_world(bodies[0].angular_momentum_body());
            same(
                world_spin,
                Vec3::new(0.0, -0.125 * tangent_impulse, 0.0),
                2e-10,
            );
            near(
                bodies[0].pose().position_world().x,
                0.25 + STEP * (1.0 - 0.5 * tangent_impulse),
                2e-13,
            );
            // m=2, I=1/80 and r=1/8 give tangential inverse mass 7/4.
            // Both cases remain sliding, so the independent boundary formula
            // Jt=mu*Jn applies rather than a zero-speed/sticking approximation.
            assert!(1.0 - 1.75 * tangent_impulse > 0.9);
            near(
                impact.relative_velocities_after[0][2],
                -1.0 + 1.75 * tangent_impulse,
                2e-10,
            );
            same(report.linear_momentum_residual_kg_m_s, Vec3::ZERO, 2e-10);
            same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
            near(report.work_energy_residual_j, 0.0, 2e-10);
            assert!(report.unallocated_contact_loss_j > 0.0);
            impulses.push(tangent_impulse);
        }
        near(impulses[0], 4.0 * impulses[1], 2e-14);
    });
}

#[test]
fn directional_coefficient_and_degenerate_law_refusals_keep_original_native_state() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let before = solid.native_state().clone();
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 0.0, 0.0)], &[], 0).unwrap();
        let initial = body_at(
            [0.25, 0.25, -0.125 - GAP],
            Vec3::ZERO,
            Vec3::new(1.0, 0.0, 0.0),
        );
        let refused = solid.sustained_elliptic_sphere_trajectory_candidate(
            &surface,
            &scene,
            &[initial],
            &[homogeneous_sphere()],
            &[],
            [0.0; 3],
            &[Wrench::ZERO],
            upward_gravity(),
            STEP,
            ORIGIN,
            sustained_controls(STEP),
            |_, _, _| Err::<EllipticFriction, _>("missing directional source"),
            cx,
        );
        assert!(matches!(
            refused,
            Err(FiniteSolidSustainedTrajectoryError::Impact(
                FiniteSolidSustainedEllipticFailure::Coefficients {
                    original_index: 0,
                    error: "missing directional source",
                    ..
                }
            ))
        ));
        let singular = solid.sustained_elliptic_sphere_trajectory_candidate(
            &surface,
            &scene,
            &[initial],
            &[homogeneous_sphere()],
            &[],
            [0.0; 3],
            &[Wrench::ZERO],
            upward_gravity(),
            STEP,
            ORIGIN,
            sustained_controls(STEP),
            |_, _, _| {
                Ok::<_, core::convert::Infallible>(EllipticFriction {
                    coefficients_mu: [0.0, 0.8],
                    principal_axis: [1.0, 0.0],
                })
            },
            cx,
        );
        assert!(matches!(
            singular,
            Err(FiniteSolidSustainedTrajectoryError::Impact(
                FiniteSolidSustainedEllipticFailure::Contact(
                    FiniteSolidRigidContactError::Contact(
                        EllipticFrictionalContactError::UnsupportedDegenerateEllipse { .. }
                    )
                )
            ))
        ));
        assert_eq!(solid.native_state(), &before);
    });
}
