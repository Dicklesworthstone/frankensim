//! G1/G3 automatic contact trajectories on actual native mechanical states.
#[path = "step/free.rs"]
mod free;
#[path = "step/scene.rs"]
mod scene;

use super::*;
use fs_couple::finite_solid::contact::step::*;
use fs_couple::finite_solid::contact::surface::*;
use fs_mbd::{Gravity, Wrench};

const FLOOR: [(usize, f64); 9] = [
    (0, 0.0),
    (1, 0.0),
    (2, 0.0),
    (3, 0.0),
    (4, 0.0),
    (5, 0.0),
    (6, 0.0),
    (7, 0.0),
    (8, 0.0),
];
fn sphere(id: u64, body: usize, restitution: f64, friction: f64) -> FiniteSolidSphere {
    FiniteSolidSphere {
        id,
        rigid_body: body,
        center_body_m: Vec3::ZERO,
        radius_m: 0.125,
        friction_coefficient: friction,
        normal_law: NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        },
    }
}
fn settings() -> FiniteSolidSphereTrajectoryControls {
    let mut c = FiniteSolidSphereTrajectoryControls {
        maximum_step_s: 0.3,
        maximum_steps: 128,
        maximum_rejections: 512,
        ..Default::default()
    };
    c.surface.activation_tolerance_m = 1e-8;
    c.impact.geometry_tolerance_m = 2e-8;
    c.free.solid = controls(1e-10);
    c
}
fn stored_solid<'a>(material: &'a FiniteSolidMaterial, cx: &Cx<'_>) -> FiniteSolidBody<'a> {
    FiniteSolidBody::new(
        problem(material, &FLOOR, &[]),
        &[[0.0; 3]; 4],
        &[[0.0; 3]; 4],
        0.0,
        Default::default(),
        cx,
    )
    .unwrap()
}

#[test]
fn detected_bounce_advances_native_time_without_tunneling_or_position_repair() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = body_at([0.25, 0.25, -0.625], Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0));
        let (next, bodies, report) = solid
            .sphere_trajectory_candidate(
                &surface,
                &[initial],
                &[properties()],
                &[sphere(7, 0, 1.0, 0.0)],
                &[],
                [0.0; 3],
                &[Wrench::ZERO],
                Gravity::ZERO,
                1.0,
                [0.0; 3],
                settings(),
                cx,
            )
            .unwrap();
        // Exact zero-gap bounce t=.5; activation moves it by at most 1e-8 m.
        near(bodies[0].pose().position_world().z, -0.625, 2.1e-8);
        near(bodies[0].linear_momentum_world().z, -2.0, 2e-10);
        near(next.native_state().t, 1.0, 2e-15);
        assert_eq!(next.native_state().steps, report.steps.len());
        assert_eq!(solid.native_state().steps, 0);
        assert_eq!(solid.native_state().q, next.native_state().q);
        assert!(report.rejected_steps > 0);
        let events: Vec<_> = report
            .steps
            .iter()
            .filter_map(|s| s.impact.as_ref())
            .filter(|i| i.contact.impulses.iter().any(|j| j[0] > 1e-8))
            .collect();
        assert_eq!(events.len(), 1);
        near(events[0].contact.impulses[0][0], 4.0, 2e-10);
        near(report.unallocated_contact_loss_j, 0.0, 2e-10);
        same(
            report.support_impulse_kg_m_s,
            Vec3::new(0.0, 0.0, -4.0),
            2e-10,
        );
        same(
            report.support_angular_impulse_kg_m2_s,
            Vec3::new(-1.0, 1.0, 0.0),
            2e-10,
        );
        same(report.linear_momentum_residual_kg_m_s, Vec3::ZERO, 2e-10);
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
        near(report.work_energy_residual_j, 0.0, 2e-10);
        for step in &report.steps {
            assert!(step.maximum_contact_gap_m <= settings().surface.activation_tolerance_m);
            assert_eq!(step.path.tetrahedra_checked, 1);
            assert!(step.path.minimum_clearance_m >= 0.0);
            for receipt in &step.free.rigid {
                assert!(receipt.state_before.pose().position_world().z <= -0.125);
                assert!(receipt.state_after.pose().position_world().z <= -0.125);
            }
        }
    });
}

#[test]
fn rough_contact_spin_and_support_moment_continue_through_native_rotation() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        // Homogeneous physical sphere: m=2, R=1/8, I=2mR^2/5=1/80.
        let property =
            MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.0125, 0.0125, 0.0125)).unwrap();
        let initial = body_at([0.25, 0.25, -0.125], Vec3::ZERO, Vec3::new(1.0, 0.0, 1.0));
        let origin = [0.75, -0.5, 0.25];
        let (next, bodies, report) = solid
            .sphere_trajectory_candidate(
                &surface,
                &[initial],
                &[property],
                &[sphere(9, 0, 0.5, 0.1)],
                &[],
                [0.0; 3],
                &[Wrench::ZERO],
                Gravity::ZERO,
                0.02,
                origin,
                settings(),
                cx,
            )
            .unwrap();
        // Wn=.5 => Jn=3. Wt=.5+R^2/I=1.75; mu*Jn=.3 < 1/Wt.
        // Sliding world impulse is (-.3,0,-3), spin increment (0,-.0375,0).
        let event = report.steps[0].impact.as_ref().unwrap();
        near(event.contact.impulses[0][0], 3.0, 2e-10);
        same(
            bodies[0].linear_momentum_world(),
            Vec3::new(1.7, 0.0, -1.0),
            2e-10,
        );
        same(
            bodies[0]
                .pose()
                .orientation()
                .rotate_body_to_world(bodies[0].angular_momentum_body()),
            Vec3::new(0.0, -0.0375, 0.0),
            2e-10,
        );
        assert_ne!(bodies[0].pose().orientation(), UnitQuaternion::IDENTITY);
        // .75 J normal loss plus .3 - .5*1.75*.3^2 tangential loss.
        near(report.unallocated_contact_loss_j, 0.97125, 3e-10);
        near(report.kinetic_change_j, -0.97125, 3e-10);
        same(
            report.support_impulse_kg_m_s,
            Vec3::new(-0.3, 0.0, -3.0),
            2e-10,
        );
        same(
            report.support_angular_impulse_kg_m2_s,
            Vec3::new(0.25 - origin[0], 0.25 - origin[1], -origin[2])
                .cross(Vec3::new(-0.3, 0.0, -3.0)),
            2e-10,
        );
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
        near(report.work_energy_residual_j, 0.0, 3e-10);
        assert_eq!(next.native_state().q, solid.native_state().q);
    });
}

#[test]
fn complete_path_refuses_unassembled_rigid_pair_with_clear_endpoints() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let bodies = [
            body_at([-1.0, 0.0, -2.0], Vec3::ZERO, Vec3::new(4.0, 0.0, 0.0)),
            body_at([1.0, 0.0, -2.0], Vec3::ZERO, Vec3::new(-4.0, 0.0, 0.0)),
        ];
        let c = FiniteSolidSphereTrajectoryControls {
            maximum_step_s: 0.5,
            ..settings()
        };
        let result = solid.sphere_trajectory_candidate(
            &surface,
            &bodies,
            &[properties(); 2],
            &[sphere(20, 0, 1.0, 0.0), sphere(10, 1, 1.0, 0.0)],
            &[],
            [0.0; 3],
            &[Wrench::ZERO; 2],
            Gravity::ZERO,
            0.5,
            [0.0; 3],
            c,
            cx,
        );
        // At t=.5 they have exchanged positions; first touch is .21875 s.
        // The returned fractional window is conservative, not an exact root.
        assert!(matches!(result,
            Err(FiniteSolidSphereTrajectoryError::UnsupportedSpherePair { first: 10, second: 20, fraction })
            if fraction <= 0.4375 && fraction > 0.4375 - 2e-8));
        assert_eq!(solid.native_state().steps, 0);
    });
}

#[test]
fn late_budget_offset_geometry_and_cancellation_publish_no_partial_motion() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = body_at([0.25, 0.25, -0.625], Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0));
        let c = FiniteSolidSphereTrajectoryControls {
            maximum_steps: 1,
            ..settings()
        };
        assert!(
            matches!(solid.sphere_trajectory_candidate(&surface, &[initial], &[properties()],
            &[sphere(7, 0, 1.0, 0.0)], &[], [0.0; 3], &[Wrench::ZERO], Gravity::ZERO,
            1.0, [0.0; 3], c, cx), Err(FiniteSolidSphereTrajectoryError::StepBudget { elapsed_s }) if elapsed_s == 0.3)
        );
        assert_eq!(solid.native_state().steps, 0);
        assert_eq!(initial.pose().position_world().z, -0.625);
        let mut offset = sphere(7, 0, 1.0, 0.0);
        offset.center_body_m = Vec3::new(0.01, 0.0, 0.0);
        assert!(matches!(
            solid.sphere_trajectory_candidate(
                &surface,
                &[initial],
                &[properties()],
                &[offset],
                &[],
                [0.0; 3],
                &[Wrench::ZERO],
                Gravity::ZERO,
                1.0,
                [0.0; 3],
                settings(),
                cx
            ),
            Err(FiniteSolidSphereTrajectoryError::Invalid(_))
        ));
        gate.request();
        assert!(matches!(
            solid.sphere_trajectory_candidate(
                &surface,
                &[initial],
                &[properties()],
                &[sphere(7, 0, 1.0, 0.0)],
                &[],
                [0.0; 3],
                &[Wrench::ZERO],
                Gravity::ZERO,
                1.0,
                [0.0; 3],
                settings(),
                cx
            ),
            Err(FiniteSolidSphereTrajectoryError::Cancelled)
        ));
    });
}

#[test]
fn large_native_clock_cannot_silently_change_the_requested_duration() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &FLOOR, &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            2f64.powi(40),
            Default::default(),
            cx,
        )
        .unwrap();
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let body = body_at([0.0, 0.0, -2.0], Vec3::ZERO, Vec3::ZERO);
        let c = FiniteSolidSphereTrajectoryControls {
            maximum_step_s: 0.001,
            ..settings()
        };
        assert!(matches!(
            solid.sphere_trajectory_candidate(
                &surface,
                &[body],
                &[properties()],
                &[sphere(0, 0, 1.0, 0.0)],
                &[],
                [0.0; 3],
                &[Wrench::ZERO],
                Gravity::ZERO,
                0.01,
                [0.0; 3],
                c,
                cx
            ),
            Err(FiniteSolidSphereTrajectoryError::TimeResolution { .. })
        ));
        assert_eq!(solid.native_state().t, 2f64.powi(40));
        assert_eq!(solid.native_state().steps, 0);
    });
}
