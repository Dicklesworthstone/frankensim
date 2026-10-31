//! Independent loaded-support and angular oracles for the actual native split.
use super::*;
use fs_couple::finite_solid::contact::scene::{
    FiniteSolidSphereScene, FiniteSolidSphereSceneKey, SpherePairLaw,
};

const GAP: f64 = 1.0 / 1_073_741_824.0;
const STEP: f64 = 1.0 / 1024.0;
const DURATION: f64 = 4.0 * STEP;
const ACCELERATION: f64 = 8.0;
const ORIGIN: [f64; 3] = [0.75, -0.5, 0.25];

fn sustained_controls(step: f64) -> FiniteSolidSphereTrajectoryControls {
    let mut controls = settings();
    controls.maximum_step_s = step;
    // The first-order force/contact split creates a KNOWN O(h) support loss.
    // Admit its explicitly reported integration defect, retaining the normal
    // tight momentum and storage-minus-work residual gates.
    controls.free.energy_tolerance_j = 1e-3;
    controls.free.solid.energy_tolerance_j = 1e-3;
    controls.integration_tolerance_j = 1e-3;
    controls
}

fn homogeneous_sphere() -> MassProperties {
    // m=2, R=1/8 gives I=2mR^2/5=1/80 about every principal axis.
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.0125, 0.0125, 0.0125)).unwrap()
}

fn upward_gravity() -> Gravity {
    Gravity::new(Vec3::new(0.0, 0.0, ACCELERATION)).unwrap()
}

fn expected_support(count: f64) -> Vec3 {
    Vec3::new(0.0, 0.0, -count * 2.0 * ACCELERATION * DURATION)
}

fn support_moment(impulse: Vec3) -> Vec3 {
    Vec3::new(0.25 - ORIGIN[0], 0.25 - ORIGIN[1], -ORIGIN[2]).cross(impulse)
}

#[test]
fn loaded_resting_sphere_advances_time_and_exposes_refinable_support_splitting_loss() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let input = solid.native_state().clone();
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = body_at([0.25, 0.25, -0.125 - GAP], Vec3::ZERO, Vec3::ZERO);
        // An elastic impact coefficient must not turn this held load into a
        // fresh bounce: restitution sees the pre-force resting velocity.
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 1.0, 0.0)], &[], 0).unwrap();
        let mut losses = Vec::new();
        for (step, expected_steps) in [(STEP, 4), (0.5 * STEP, 8)] {
            let (next, bodies, report) = solid
                .sustained_sphere_trajectory_candidate(
                    &surface,
                    &scene,
                    &[initial],
                    &[homogeneous_sphere()],
                    &[],
                    [0.0; 3],
                    &[Wrench::ZERO],
                    upward_gravity(),
                    DURATION,
                    ORIGIN,
                    sustained_controls(step),
                    cx,
                )
                .unwrap();
            near(next.native_state().t, DURATION, 0.0);
            assert_eq!(report.steps.len(), expected_steps);
            assert_eq!(next.native_state().steps, expected_steps);
            assert_eq!(next.native_state().history.len(), expected_steps);
            assert_eq!(next.native_state().q, input.q);
            same(
                bodies[0].pose().position_world(),
                initial.pose().position_world(),
                3e-11,
            );
            same(bodies[0].linear_momentum_world(), Vec3::ZERO, 2e-10);
            same(bodies[0].angular_momentum_body(), Vec3::ZERO, 2e-10);
            same(report.support_impulse_kg_m_s, expected_support(1.0), 2e-10);
            same(
                report.support_angular_impulse_kg_m2_s,
                support_moment(expected_support(1.0)),
                2e-10,
            );
            same(
                report.external_impulse_kg_m_s,
                expected_support(-1.0),
                2e-10,
            );
            same(report.linear_momentum_residual_kg_m_s, Vec3::ZERO, 2e-10);
            same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
            near(report.source_work_j, 0.0, 2e-10);
            near(report.kinetic_change_j, 0.0, 2e-10);
            near(report.strain_change_j, 0.0, 2e-10);
            near(report.work_energy_residual_j, 0.0, 2e-10);
            // Each force kick creates 1/2 m (g h)^2 kinetic energy, removed
            // by zero-target support. It is a split defect, not continuum
            // friction heat; its whole-trajectory value is 1/2 m g^2 T h.
            let expected_loss = ACCELERATION * ACCELERATION * DURATION * step;
            near(report.unallocated_contact_loss_j, expected_loss, 2e-10);
            near(report.integration_defect_j, expected_loss, 2e-10);
            losses.push(report.unallocated_contact_loss_j);
            for (history, accepted) in next.native_state().history.iter().zip(&report.steps) {
                assert_eq!(history, &accepted.solid_drift.time);
                assert!(history.newton.converged);
                assert!(accepted.solid_drift.endpoint_acceleration.converged);
                assert!(accepted.path.minimum_clearance_m >= 0.0);
                assert_eq!(accepted.incoming_normal_speeds_m_s.len(), 1);
                assert_eq!(accepted.normal_targets_m_s.len(), 1);
                near(accepted.incoming_normal_speeds_m_s[0], 0.0, 2e-10);
                near(accepted.normal_targets_m_s[0], 0.0, 2e-10);
                let impact = accepted.impact.as_ref().unwrap();
                assert_eq!(impact.contact.impulses.len(), 1);
                near(
                    impact.contact.impulses[0][0],
                    2.0 * ACCELERATION * step,
                    2e-10,
                );
            }
        }
        near(losses[1], 0.5 * losses[0], 2e-10);
        assert_eq!(solid.native_state(), &input);
    });
}

#[test]
fn loaded_resting_chain_transmits_neighbor_force_through_one_mixed_solve() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = [
            body_at([0.25, 0.25, -0.125 - GAP], Vec3::ZERO, Vec3::ZERO),
            body_at([0.25, 0.25, -0.375 - 2.0 * GAP], Vec3::ZERO, Vec3::ZERO),
        ];
        // Canonical sphere identity runs opposite to native body addressing.
        let scene = FiniteSolidSphereScene::new(
            &[sphere(9, 0, 0.0, 0.0), sphere(3, 1, 0.0, 0.0)],
            &[SpherePairLaw {
                first: 3,
                second: 9,
                friction_coefficient: 0.0,
                restitution: 0.0,
            }],
            1,
        )
        .unwrap();
        let (next, bodies, report) = solid
            .sustained_sphere_trajectory_candidate(
                &surface,
                &scene,
                &initial,
                &[homogeneous_sphere(); 2],
                &[],
                [0.0; 3],
                &[Wrench::ZERO; 2],
                upward_gravity(),
                DURATION,
                ORIGIN,
                sustained_controls(STEP),
                cx,
            )
            .unwrap();
        near(next.native_state().t, DURATION, 0.0);
        assert_eq!(report.steps.len(), 4);
        assert_eq!(next.native_state().q, solid.native_state().q);
        for (before, after) in initial.iter().zip(&bodies) {
            same(
                after.pose().position_world(),
                before.pose().position_world(),
                3e-11,
            );
            same(after.linear_momentum_world(), Vec3::ZERO, 2e-10);
            same(after.angular_momentum_body(), Vec3::ZERO, 2e-10);
        }
        for accepted in &report.steps {
            let impact = accepted.impact.as_ref().unwrap();
            assert_eq!(impact.contact.impulses.len(), 2);
            assert_eq!(accepted.contact_keys.len(), 2);
            for (key, impulse) in accepted.contact_keys.iter().zip(&impact.contact.impulses) {
                let expected = match key {
                    FiniteSolidSphereSceneKey::Solid { sphere: 9, .. } => 4.0 * ACCELERATION * STEP,
                    FiniteSolidSphereSceneKey::Spheres {
                        first: 3,
                        second: 9,
                    } => 2.0 * ACCELERATION * STEP,
                    other => panic!("unexpected mixed support key: {other:?}"),
                };
                near(impulse[0], expected, 2e-10);
                near(impulse[1], 0.0, 2e-10);
                near(impulse[2], 0.0, 2e-10);
            }
            assert_eq!(accepted.path.sphere_pairs_checked, 1);
            assert!(accepted.path.minimum_clearance_m >= 0.0);
        }
        same(report.support_impulse_kg_m_s, expected_support(2.0), 2e-10);
        same(
            report.support_angular_impulse_kg_m2_s,
            support_moment(expected_support(2.0)),
            2e-10,
        );
        same(
            report.external_impulse_kg_m_s,
            expected_support(-2.0),
            2e-10,
        );
        same(report.linear_momentum_residual_kg_m_s, Vec3::ZERO, 2e-10);
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
        near(report.work_energy_residual_j, 0.0, 2e-10);
        near(
            report.unallocated_contact_loss_j,
            2.0 * ACCELERATION * ACCELERATION * DURATION * STEP,
            3e-10,
        );
        assert_eq!(solid.native_state().steps, 0);
    });
}

#[test]
fn applied_normal_axis_torque_preserves_spin_while_support_carries_gravity() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = body_at([0.25, 0.25, -0.125 - GAP], Vec3::ZERO, Vec3::ZERO);
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 0.0, 0.0)], &[], 0).unwrap();
        let torque = 0.25;
        let wrench = Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(0.0, 0.0, torque),
        };
        let (_, bodies, report) = solid
            .sustained_sphere_trajectory_candidate(
                &surface,
                &scene,
                &[initial],
                &[homogeneous_sphere()],
                &[],
                [0.0; 3],
                &[wrench],
                upward_gravity(),
                DURATION,
                ORIGIN,
                sustained_controls(STEP),
                cx,
            )
            .unwrap();
        let spin = Vec3::new(0.0, 0.0, torque * DURATION);
        let body = bodies[0];
        same(body.linear_momentum_world(), Vec3::ZERO, 2e-10);
        same(body.angular_momentum_body(), spin, 2e-10);
        same(
            body.pose()
                .orientation()
                .rotate_body_to_world(body.angular_momentum_body()),
            spin,
            2e-10,
        );
        same(report.angular_momentum_change_kg_m2_s, spin, 2e-10);
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
        same(report.support_impulse_kg_m_s, expected_support(1.0), 2e-10);
        same(
            report.support_angular_impulse_kg_m2_s,
            support_moment(expected_support(1.0)),
            2e-10,
        );
        // Four force/contact/drift intervals use successive post-kick speeds:
        // theta = (tau/I) h^2 (1+2+3+4). The native constant-axis drift is exact.
        let angle = (torque / 0.0125) * STEP * STEP * 10.0;
        same(
            body.pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            Vec3::new(angle.cos(), angle.sin(), 0.0),
            2e-10,
        );
        assert_ne!(body.pose().orientation(), UnitQuaternion::IDENTITY);
        near(report.work_energy_residual_j, 0.0, 3e-10);
    });
}

#[test]
fn loaded_contact_recoils_against_deforming_consistent_mass_and_keeps_support_reaction() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        // Only the three base-node z coordinates move. The apex and all x/y
        // coordinates stay prescribed, permitting a real uniform first-step
        // face displacement without replacing the tet by a rigid support.
        let fixed: Vec<_> = (0..12)
            .filter(|dof| ![2, 5, 8].contains(dof))
            .map(|dof| (dof, 0.0))
            .collect();
        let solid = FiniteSolidBody::new(
            problem(&material, &fixed, &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let initial = body_at([1.0 / 3.0, 1.0 / 3.0, -0.125 - GAP], Vec3::ZERO, Vec3::ZERO);
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 0.0, 0.0)], &[], 0).unwrap();
        let (next, bodies, report) = solid
            .sustained_sphere_trajectory_candidate(
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
                cx,
            )
            .unwrap();
        // The one-kilogram tet's free consistent mass is (I+11^T)/20.
        // Its uniform face mode has mass 3*(4/20)=0.6, hence the shared
        // outgoing speed is (2 g h)/(2+0.6), not a lumped-node approximation.
        let velocity = 2.0 * ACCELERATION * STEP / 2.6;
        let displacement = STEP * velocity;
        let impulse = 2.0 * ACCELERATION * STEP - 2.0 * velocity;
        assert_eq!(report.steps.len(), 1);
        near(next.native_state().t, STEP, 0.0);
        near(bodies[0].linear_momentum_world().z, 2.0 * velocity, 2e-10);
        near(
            bodies[0].pose().position_world().z,
            initial.pose().position_world().z + displacement,
            3e-12,
        );
        for node in 0..3 {
            near(next.velocity_m_s()[node][2], velocity, 2e-10);
            near(next.displacement_m()[node][2], displacement, 3e-12);
            near(next.displacement_m()[node][0], 0.0, 0.0);
            near(next.displacement_m()[node][1], 0.0, 0.0);
        }
        assert_eq!(next.displacement_m()[3], [0.0; 3]);
        let accepted = &report.steps[0];
        let impact = accepted.impact.as_ref().unwrap();
        near(impact.contact_mobility_per_kg[0], 0.5 + 5.0 / 3.0, 2e-11);
        near(impact.contact.impulses[0][0], impulse, 2e-10);
        // Consistent mass gives the fixed apex momentum 3v/20. Its constraint
        // supplies this external reaction even though the support does no work.
        let support = Vec3::new(0.0, 0.0, 0.15 * velocity);
        same(report.support_impulse_kg_m_s, support, 2e-10);
        same(
            report.support_angular_impulse_kg_m2_s,
            Vec3::new(-ORIGIN[0], -ORIGIN[1], 1.0 - ORIGIN[2]).cross(support),
            2e-10,
        );
        // This affine endpoint has F=diag(1,1,1-hv); compare with the parent
        // test's independently derived compressible neo-Hookean energy.
        near(
            accepted.solid_drift.strain_after_j,
            strain_energy(1.0 - displacement, 1.0),
            2e-13,
        );
        near(
            report.source_work_j,
            2.0 * ACCELERATION * displacement,
            2e-10,
        );
        same(report.linear_momentum_residual_kg_m_s, Vec3::ZERO, 2e-10);
        same(report.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
        near(report.work_energy_residual_j, 0.0, 2e-10);
        assert_eq!(solid.native_state().steps, 0);
        assert_eq!(solid.displacement_m(), [[0.0; 3]; 4]);
    });
}
