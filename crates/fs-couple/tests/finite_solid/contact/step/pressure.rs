//! Actual pressure-loaded deformable impacts and pressure-preserving retries.
use super::*;
use fs_couple::finite_solid::contact::scene::FiniteSolidSphereScene;
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;
use std::sync::Arc;

fn pressure_surface(cx: &Cx<'_>) -> Arc<ClosedTetPressureSurface> {
    Arc::new(
        ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap(),
    )
}
fn volume(solid: &FiniteSolidBody<'_>) -> f64 {
    let u = solid.displacement_m();
    let f: [f64; 9] =
        core::array::from_fn(|ij| NODES[ij % 3 + 1][ij / 3] + u[ij % 3 + 1][ij / 3] - u[0][ij / 3]);
    (f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6]))
        / 6.0
}
fn audit(report: &FiniteSolidSphereTrajectoryReport) {
    same(report.linear_momentum_residual_kg_m_s, Vec3::ZERO, 1e-9);
    same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 1e-9);
    near(report.work_energy_residual_j, 0.0, 1e-9);
    near(report.integration_defect_j, 0.0, 1e-9);
    assert!(report.absolute_integration_defect_j < 1e-9);
    near(
        report.kinetic_change_j + report.strain_change_j - report.source_work_j
            + report.unallocated_contact_loss_j,
        0.0,
        1e-9,
    );
}

#[test]
fn freely_recoiling_pressurized_solid_and_rigid_impact_share_one_work_account() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let u = NODES.map(|x| x.map(|v| 0.125 * v));
        let solid = FiniteSolidBody::new_with_pressure(
            problem(&material, &[], &[]),
            pressure_surface(cx),
            3.0,
            &u,
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = body_at(
            [0.28125, 0.28125, -0.125],
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, 1.0),
        );
        let before = solid.native_state().clone();
        let mut c = settings();
        c.maximum_step_s = 0.005;
        let (next, bodies, report) = solid
            .sphere_trajectory_candidate_with_pressure(
                &surface,
                &[initial],
                &[properties()],
                &[sphere(17, 0, 0.5, 0.0)],
                &[],
                [0.0; 3],
                3.0,
                &[Wrench::ZERO],
                Gravity::ZERO,
                0.005,
                [0.0; 3],
                c,
                cx,
            )
            .unwrap();
        // Reference mass is 1 kg even at stretch 1.125. At weights (.5,.25,.25,0),
        // the exact consistent point mobility is 20*(.375-.2)=3.5/kg. Adding
        // the 2 kg rigid's radial mobility gives 4/kg: Jn=1.5/4 and loss=.75/8.
        let impact = report.steps[0].impact.as_ref().unwrap();
        near(impact.contact.impulses[0][0], 0.375, 2e-10);
        near(report.unallocated_contact_loss_j, 0.09375, 3e-10);
        near(bodies[0].linear_momentum_world().z, 1.625, 3e-10);
        let pressure_work: f64 = report
            .steps
            .iter()
            .map(|s| s.free.solid.pressure_work_j)
            .sum();
        near(
            pressure_work,
            -3.0 * (volume(&next) - volume(&solid)),
            1e-10,
        );
        assert!(pressure_work > 1e-6);
        near(report.source_work_j, pressure_work, 1e-12);
        same(report.support_impulse_kg_m_s, Vec3::ZERO, 1e-12);
        same(report.external_impulse_kg_m_s, Vec3::ZERO, 1e-12);
        audit(&report);
        assert_eq!(solid.native_state(), &before);
        assert_eq!(next.native_state().steps, report.steps.len());
    });
}

#[test]
fn mixed_pressure_trajectory_retains_pressure_through_geometric_retries_and_zero_parity() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let mut u = [[0.0; 3]; 4];
        u[3][2] = 0.125;
        let problem = problem(&material, &FLOOR, &[]);
        let plain =
            FiniteSolidBody::new(problem, &u, &[[0.0; 3]; 4], 0.0, Default::default(), cx).unwrap();
        let bound = FiniteSolidBody::new_with_pressure(
            problem,
            pressure_surface(cx),
            0.0,
            &u,
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let surface = FiniteSolidContactSurface::new(&bound, 64, cx).unwrap();
        let plain_surface = FiniteSolidContactSurface::new(&plain, 64, cx).unwrap();
        let initial = body_at([0.25, 0.25, -0.127], Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0));
        let scene = FiniteSolidSphereScene::new(&[sphere(5, 0, 1.0, 0.0)], &[], 1).unwrap();
        let mut c = settings();
        c.maximum_step_s = 0.005;
        let (next, bodies, report) = bound
            .mixed_sphere_trajectory_candidate_with_pressure(
                &surface,
                &scene,
                &[initial],
                &[properties()],
                &[],
                [0.0; 3],
                3.0,
                &[Wrench::ZERO],
                Gravity::ZERO,
                0.005,
                [0.0; 3],
                c,
                cx,
            )
            .unwrap();
        assert!(report.rejected_steps > 0);
        near(bodies[0].pose().position_world().z, -0.128, 2.1e-8);
        near(bodies[0].linear_momentum_world().z, -2.0, 3e-10);
        let work: f64 = report
            .steps
            .iter()
            .map(|s| s.free.solid.pressure_work_j)
            .sum();
        near(work, -3.0 * (volume(&next) - volume(&bound)), 1e-10);
        assert!(work > 1e-6);
        audit(&report);
        let zero = bound
            .mixed_sphere_trajectory_candidate_with_pressure(
                &surface,
                &scene,
                &[initial],
                &[properties()],
                &[],
                [0.0; 3],
                0.0,
                &[Wrench::ZERO],
                Gravity::ZERO,
                0.005,
                [0.0; 3],
                c,
                cx,
            )
            .unwrap();
        let old = plain
            .mixed_sphere_trajectory_candidate(
                &plain_surface,
                &scene,
                &[initial],
                &[properties()],
                &[],
                [0.0; 3],
                &[Wrench::ZERO],
                Gravity::ZERO,
                0.005,
                [0.0; 3],
                c,
                cx,
            )
            .unwrap();
        assert_eq!(zero.0.native_state(), old.0.native_state());
        assert_eq!(zero.1, old.1);
        assert_eq!(zero.2.rejected_steps, old.2.rejected_steps);
        assert_eq!(
            zero.2.source_work_j.to_bits(),
            old.2.source_work_j.to_bits()
        );
    });
}

#[test]
fn pressure_binding_and_cancellation_refuse_before_any_native_state_is_published() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let plain = stored_solid(&material, cx);
        let bound = FiniteSolidBody::new_with_pressure(
            problem(&material, &FLOOR, &[]),
            pressure_surface(cx),
            0.0,
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let surface = FiniteSolidContactSurface::new(&bound, 64, cx).unwrap();
        let before = bound.native_state().clone();
        assert!(
            plain
                .sphere_trajectory_candidate_with_pressure(
                    &surface,
                    &[],
                    &[],
                    &[],
                    &[],
                    [0.0; 3],
                    3.0,
                    &[],
                    Gravity::ZERO,
                    0.001,
                    [0.0; 3],
                    settings(),
                    cx
                )
                .is_err()
        );
        assert!(
            bound
                .sphere_trajectory_candidate_with_pressure(
                    &surface,
                    &[],
                    &[],
                    &[],
                    &[],
                    [0.0; 3],
                    f64::NAN,
                    &[],
                    Gravity::ZERO,
                    0.001,
                    [0.0; 3],
                    settings(),
                    cx
                )
                .is_err()
        );
        gate.request();
        assert!(
            bound
                .sphere_trajectory_candidate_with_pressure(
                    &surface,
                    &[],
                    &[],
                    &[],
                    &[],
                    [0.0; 3],
                    3.0,
                    &[],
                    Gravity::ZERO,
                    0.001,
                    [0.0; 3],
                    settings(),
                    cx
                )
                .is_err()
        );
        assert_eq!(bound.native_state(), &before);
    });
}
