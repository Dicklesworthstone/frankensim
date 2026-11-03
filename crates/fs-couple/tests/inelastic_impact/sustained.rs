//! Loaded resting chain, refinement and accepted-only native Maxwell replay.
use super::*;
use fs_couple::finite_solid::contact::step::FiniteSolidSphereTrajectoryError;
use fs_couple::inelastic_impact::{
    surface::{scene::*, *},
    trajectory::sustained::*,
};
use fs_mbd::{Gravity, Wrench};
use fs_solid::hyper3::visco::dynamics::FiniteViscoContactDriftControls;

const GAP: f64 = 1.0 / 1_073_741_824.0;
const STEP: f64 = 1.0 / 1024.0;
const DURATION: f64 = 4.0 * STEP;
const G: f64 = 8.0;

fn settings_for(step: f64) -> ViscoSustainedSphereTrajectoryControls {
    let mut impact = mixed_controls();
    impact.geometry_tolerance_m = 2e-8;
    ViscoSustainedSphereTrajectoryControls {
        maximum_step_s: step,
        maximum_steps: 128,
        maximum_rejections: 512,
        path: Default::default(),
        surface: FiniteSolidSurfaceControls {
            activation_tolerance_m: 1e-8,
            ..Default::default()
        },
        impact,
        solid: FiniteViscoContactDriftControls {
            material_balance_tolerance_j: 1e-9,
            maximum_integration_defect_j: 1e-3,
        },
        linear_tolerance_kg_m_s: 1e-7,
        angular_tolerance_kg_m2_s: 1e-7,
        energy_tolerance_j: 1e-7,
        integration_tolerance_j: 1e-3,
    }
}

fn material() -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 1200.0,
                lambda: 2400.0,
            },
            4.0,
        )
        .unwrap(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 3600.0 },
            tau_s: 0.07,
        }],
    )
    .unwrap()
}

fn relaxed<'a>(law: &'a ViscoHyperelastic, cx: &Cx<'_>) -> FiniteViscoTetMotion<'a> {
    let mut solid = FiniteViscoTetSolid::try_new(
        &NODES,
        &TETS,
        ViscoelasticLawField::Uniform(law),
        TetAssemblyBudget::default(),
        cx,
    )
    .unwrap();
    // Relaxed in-plane shear keeps the supporting z=0 face planar and fixed.
    let fixed: Vec<_> = (0..12)
        .map(|d| (d, if d == 4 { 1e-4 } else { 0.0 }))
        .collect();
    for dt in [0.0, 0.007] {
        solid = solid
            .advance(
                dt,
                &fixed,
                &[],
                [0.0; 3],
                settings(),
                FiniteViscoStepControls {
                    equilibrium_work_tolerance_j: 1e-9,
                    maximum_endpoint_work_defect_j: 1e-2,
                },
                cx,
            )
            .unwrap()
            .solid;
    }
    assert!(solid.element_states()[0].dissipated > 0.0);
    FiniteViscoTetMotion::try_new_with_nodal_density(solid, &RHO, &[[0.0; 3]; 4], cx).unwrap()
}

fn run<'a>(
    motion: &FiniteViscoTetMotion<'a>,
    controls: ViscoSustainedSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<ViscoSustainedSphereTrajectoryCandidate<'a>, ViscoSustainedSphereTrajectoryError> {
    let surface = FiniteSolidContactSurface::for_visco(motion, 4, cx)
        .map_err(ViscoSustainedSphereTrajectoryError::Surface)?;
    let bodies = [-0.125 - GAP, -0.375 - 2.0 * GAP].map(|z| {
        RigidBodyState::new(
            Pose::new(Vec3::new(0.25, 0.25, z), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap()
    });
    // Canonical identity is deliberately opposite native body indexing.
    let spheres = [(9, 0), (3, 1)].map(|(id, rigid_body)| FiniteSolidSphere {
        id,
        rigid_body,
        center_body_m: Vec3::ZERO,
        radius_m: 0.125,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 1.0 },
    });
    let scene = FiniteSolidSphereScene::new(
        &spheres,
        &[SpherePairLaw {
            first: 3,
            second: 9,
            friction_coefficient: 0.0,
            restitution: 1.0,
        }],
        1,
    )
    .unwrap();
    visco_sustained_sphere_trajectory_candidate(
        motion,
        &surface,
        &scene,
        &bodies,
        &[properties(); 2],
        &(0..12).collect::<Vec<_>>(),
        &[],
        [0.0; 3],
        &[Wrench::ZERO; 2],
        Gravity::new(Vec3::new(0.0, 0.0, G)).unwrap(),
        DURATION,
        ORIGIN,
        controls,
        cx,
    )
}

fn replay(
    law: &ViscoHyperelastic,
    before: &FiniteViscoTetMotion<'_>,
    after: &ViscoSustainedSphereTrajectoryCandidate<'_>,
) {
    let mut state = before.solid().element_states()[0].clone();
    let mut elapsed = 0.0;
    for step in &after.report.steps {
        let response = law
            .step_response(&state, &state.deformation_gradient, step.duration_s)
            .unwrap();
        near(
            step.solid_drift.dissipated_j,
            response.dissipated_pa / 6.0,
            1e-16,
        );
        state = response.state;
        elapsed += step.duration_s;
        near(step.elapsed_s, elapsed, 0.0);
    }
    assert_eq!(after.solid.solid().element_states(), &[state]);
    assert_eq!(
        after.solid.solid().displacement_m(),
        before.solid().displacement_m()
    );
    assert_eq!(
        after.solid.reference_nodal_density_kg_m3(),
        Some(RHO.as_slice())
    );
    near(elapsed, DURATION, 0.0);
    near(
        after.solid.solid().time_s(),
        before.solid().time_s() + DURATION,
        2e-15,
    );
}

#[test]
fn loaded_viscoelastic_support_transmits_chain_force_without_repeated_restitution() {
    with_cx(|_, cx| {
        let law = material();
        let motion = relaxed(&law, cx);
        let original = motion.solid().element_states().to_vec();
        let mut losses = Vec::new();
        for (step, count) in [(STEP, 4), (0.5 * STEP, 8)] {
            let event = run(&motion, settings_for(step), cx).unwrap();
            assert_eq!(event.report.steps.len(), count);
            assert_eq!(event.report.rejected_steps, 0);
            for (body, z) in event
                .rigid_bodies
                .iter()
                .zip([-0.125 - GAP, -0.375 - 2.0 * GAP])
            {
                near(body.pose().position_world().z, z, 3e-11);
                near(norm(body.linear_momentum_world()), 0.0, 2e-10);
                near(norm(body.angular_momentum_body()), 0.0, 2e-10);
            }
            for accepted in &event.report.steps {
                near(accepted.duration_s, step, 0.0);
                assert_eq!(accepted.path.sphere_pairs_checked, 1);
                assert!(accepted.path.minimum_clearance_m >= 0.0);
                assert!(accepted.maximum_contact_gap_m <= 1e-8);
                assert_eq!(accepted.normal_targets_m_s.len(), 2);
                for (&incoming, &target) in accepted
                    .incoming_normal_speeds_m_s
                    .iter()
                    .zip(&accepted.normal_targets_m_s)
                {
                    near(incoming, 0.0, 2e-10);
                    near(target, 0.0, 2e-10); // e=1 still must not bounce a resting load.
                }
                let impact = accepted.impact.as_ref().unwrap();
                assert_eq!(impact.contact.impulses.len(), 2);
                for (key, impulse) in accepted.contact_keys.iter().zip(&impact.contact.impulses) {
                    let expected = match key {
                        FiniteSolidSphereSceneKey::Solid { sphere: 9, .. } => 4.0 * G * step,
                        FiniteSolidSphereSceneKey::Spheres {
                            first: 3,
                            second: 9,
                        } => 2.0 * G * step,
                        key => panic!("unexpected support row {key:?}"),
                    };
                    near(impulse[0], expected, 2e-10);
                    near(impulse[1], 0.0, 2e-10);
                    near(impulse[2], 0.0, 2e-10);
                }
            }
            let r = &event.report;
            let support = Vec3::new(0.0, 0.0, -4.0 * G * DURATION);
            near(norm(r.support_impulse_kg_m_s.sub(support)), 0.0, 2e-10);
            let moment = Vec3::new(0.25, 0.25, 0.0)
                .sub(vector(ORIGIN))
                .cross(support);
            near(
                norm(r.support_angular_impulse_kg_m2_s.sub(moment)),
                0.0,
                2e-10,
            );
            near(norm(r.external_impulse_kg_m_s.add(support)), 0.0, 2e-10);
            near(norm(r.linear_momentum_residual_kg_m_s), 0.0, 2e-10);
            near(norm(r.angular_momentum_residual_kg_m2_s), 0.0, 2e-10);
            near(r.source_work_j, 0.0, 2e-10);
            near(r.support_work_j, 0.0, 2e-10);
            near(r.kinetic_change_j, 0.0, 2e-10);
            near(r.stored_change_j + r.relaxation_loss_j, 0.0, 2e-10);
            near(r.work_energy_residual_j, 0.0, 2e-10);
            assert!(r.relaxation_loss_j > 0.0);
            // Two m=2 spheres: sum 1/2*m*g^2*T*h. This known numerical
            // support loss halves under refinement and is not Maxwell heat.
            let expected_loss = 2.0 * G * G * DURATION * step;
            near(r.unallocated_contact_loss_j, expected_loss, 2e-10);
            near(r.integration_defect_j, expected_loss, 2e-10);
            near(r.rigid_force_splitting_work_j, expected_loss, 2e-10);
            losses.push(r.unallocated_contact_loss_j);
            replay(&law, &motion, &event);
        }
        near(losses[1], 0.5 * losses[0], 2e-10);
        assert_eq!(motion.solid().element_states(), original);
        near(motion.solid().time_s(), 0.007, 0.0);
    });
}

#[test]
fn proportional_budget_replays_whole_trials_and_late_refusal_keeps_input_history() {
    with_cx(|gate, cx| {
        let law = material();
        let motion = relaxed(&law, cx);
        let original = motion.solid().element_states().to_vec();
        let stress = motion.solid().element_first_piola_pa().to_vec();
        let coarse_loss = 2.0 * G * G * DURATION * STEP;
        let mut tight = settings_for(STEP);
        tight.integration_tolerance_j = coarse_loss / 7.0;
        let refined = run(&motion, tight, cx).unwrap();
        assert!(refined.report.rejected_steps > 0);
        assert!(refined.report.steps.iter().all(|s| s.duration_s < STEP));
        assert!(refined.report.absolute_integration_defect_j <= tight.integration_tolerance_j);
        assert!(refined.report.unallocated_contact_loss_j < 0.2 * coarse_loss);
        near(
            refined.report.support_impulse_kg_m_s.z,
            -4.0 * G * DURATION,
            2e-10,
        );
        replay(&law, &motion, &refined);
        let mut no_retries = tight;
        no_retries.maximum_rejections = 0;
        assert!(matches!(
            run(&motion, no_retries, cx),
            Err(ViscoSustainedSphereTrajectoryError::RejectionBudget {
                elapsed_s: 0.0,
                last: ViscoSustainedSphereRefusal::Integration { .. }
            })
        ));
        let mut late = settings_for(STEP);
        late.maximum_steps = 1;
        assert!(matches!(run(&motion, late, cx),
            Err(ViscoSustainedSphereTrajectoryError::Path(
                FiniteSolidSphereTrajectoryError::StepBudget { elapsed_s })) if elapsed_s > 0.0));
        assert_eq!(motion.solid().element_states(), original);
        assert_eq!(motion.solid().element_first_piola_pa(), stress);
        near(motion.solid().time_s(), 0.007, 0.0);
        gate.request();
        assert!(run(&motion, tight, cx).is_err());
        assert_eq!(motion.solid().element_states(), original);
    });
}
