//! Pressure-loaded sustained contact through every original mechanical dispatch.
use super::*;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_couple::finite_solid::contact::step::sustained::FiniteSolidSustainedTrajectoryError;
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;
use std::sync::Arc;

pub(super) fn pressure_body<'a>(
    material: &'a FiniteSolidMaterial,
    cx: &Cx<'_>,
) -> FiniteSolidBody<'a> {
    let mut u = [[0.0; 3]; 4];
    u[3][2] = 0.125;
    let pressure = Arc::new(
        ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap(),
    );
    FiniteSolidBody::new_with_pressure(
        problem(material, &FLOOR, &[]),
        pressure,
        3.0,
        &u,
        &[[0.0; 3]; 4],
        0.0,
        Default::default(),
        cx,
    )
    .unwrap()
}
fn volume(body: &FiniteSolidBody<'_>) -> f64 {
    let u = body.displacement_m();
    let f: [f64; 9] =
        core::array::from_fn(|ij| NODES[ij % 3 + 1][ij / 3] + u[ij % 3 + 1][ij / 3] - u[0][ij / 3]);
    (f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6]))
        / 6.0
}

#[test]
fn circular_and_elliptic_sustained_pressure_share_load_work_and_support_balance() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = pressure_body(&material, cx);
        let before = solid.native_state().clone();
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 1.0, 0.0)], &[], 0).unwrap();
        let initial = body_at([0.25, 0.25, -0.125 - GAP], Vec3::ZERO, Vec3::ZERO);
        let (next, bodies, report) = solid
            .sustained_sphere_trajectory_candidate_with_pressure(
                &surface,
                &scene,
                &[initial],
                &[homogeneous_sphere()],
                &[],
                [0.0; 3],
                3.0,
                &[Wrench::ZERO],
                upward_gravity(),
                DURATION,
                ORIGIN,
                sustained_controls(STEP),
                cx,
            )
            .unwrap();
        let mut pressure_work = 0.0;
        let mut expected_contact_loss = 0.0;
        for step in &report.steps {
            pressure_work += step.solid_drift.pressure_work_j;
            let impulse = step.impact.as_ref().unwrap().contact.impulses[0][0];
            near(impulse, 2.0 * ACCELERATION * step.duration_s, 2e-10);
            expected_contact_loss += ACCELERATION.powi(2) * step.duration_s.powi(2);
            assert!(
                step.solid_drift.absolute_integration_defect_j
                    >= step.solid_drift.pressure_splitting_work_j.abs()
            );
        }
        assert!(pressure_work > 1e-6);
        near(
            pressure_work,
            -3.0 * (volume(&next) - volume(&solid)),
            2e-12,
        );
        near(report.source_work_j, pressure_work, 2e-12);
        near(
            report.unallocated_contact_loss_j,
            expected_contact_loss,
            2e-10,
        );
        near(report.work_energy_residual_j, 0.0, 2e-10);
        same(report.linear_momentum_residual_kg_m_s, Vec3::ZERO, 2e-10);
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
        same(bodies[0].linear_momentum_world(), Vec3::ZERO, 2e-10);
        let (e_solid, e_bodies, elliptic) = solid
            .sustained_elliptic_sphere_trajectory_candidate_with_pressure(
                &surface,
                &scene,
                &[initial],
                &[homogeneous_sphere()],
                &[],
                [0.0; 3],
                3.0,
                &[Wrench::ZERO],
                upward_gravity(),
                DURATION,
                ORIGIN,
                sustained_controls(STEP),
                |_, _, _| {
                    Ok::<_, core::convert::Infallible>(EllipticFriction {
                        coefficients_mu: [0.0; 2],
                        principal_axis: [1.0, 0.0],
                    })
                },
                cx,
            )
            .unwrap();
        assert_eq!(next.native_state(), e_solid.native_state());
        assert_eq!(bodies, e_bodies);
        near(elliptic.source_work_j, report.source_work_j, 2e-12);
        near(
            elliptic.unallocated_contact_loss_j,
            report.unallocated_contact_loss_j,
            2e-12,
        );
        assert_eq!(solid.native_state(), &before);
    });
}

#[test]
fn pressure_defect_refines_whole_intervals_and_late_budget_failure_is_atomic() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = pressure_body(&material, cx);
        let before = solid.native_state().clone();
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 1.0, 0.0)], &[], 0).unwrap();
        let initial = body_at([0.25, 0.25, -0.125 - GAP], Vec3::ZERO, Vec3::ZERO);
        let run = |controls| {
            solid.sustained_sphere_trajectory_candidate_with_pressure(
                &surface,
                &scene,
                &[initial],
                &[homogeneous_sphere()],
                &[],
                [0.0; 3],
                3.0,
                &[Wrench::ZERO],
                upward_gravity(),
                DURATION,
                ORIGIN,
                controls,
                cx,
            )
        };
        let coarse = run(sustained_controls(STEP)).unwrap();
        let mut controls = sustained_controls(STEP);
        controls.integration_tolerance_j = 0.25 * coarse.2.absolute_integration_defect_j;
        let refined = run(controls).unwrap();
        assert!(refined.2.rejected_steps > 0);
        assert!(refined.2.steps.len() > coarse.2.steps.len());
        assert!(refined.2.absolute_integration_defect_j <= controls.integration_tolerance_j);
        let pressure_work: f64 = refined
            .2
            .steps
            .iter()
            .map(|step| step.solid_drift.pressure_work_j)
            .sum();
        near(
            pressure_work,
            -3.0 * (volume(&refined.0) - volume(&solid)),
            2e-12,
        );
        assert!(pressure_work > 1e-6);
        controls.maximum_steps = 1;
        assert!(matches!(
            run(controls),
            Err(FiniteSolidSustainedTrajectoryError::Trajectory(
                FiniteSolidSphereTrajectoryError::StepBudget { .. }
            ))
        ));
        assert_eq!(solid.native_state(), &before);
        near(solid.native_state().t, 0.0, 0.0);
    });
}
