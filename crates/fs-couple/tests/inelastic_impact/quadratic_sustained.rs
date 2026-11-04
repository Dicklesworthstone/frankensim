//! Loaded curved dwell, refinable split loss and full-inertia friction transfer.

use super::*;
use fs_couple::finite_solid::contact::step::FiniteSolidSphereTrajectoryError;
use fs_couple::inelastic_impact::{
    surface::{
        FiniteSolidSphere, FiniteSolidSurfaceControls,
        quadratic::{
            QuadraticSolidContactSurface, QuadraticSolidSurfaceControls,
            path::QuadraticSolidPathControls,
            scene::{
                FiniteSolidSphereScene, FiniteSolidSphereSceneKey, QuadraticSphereScenePathVerdict,
                SpherePairLaw,
            },
            volume::QuadraticSolidVolumePathVerdict,
        },
    },
    trajectory::quadratic::sustained::*,
};
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_mbd::{Gravity, Wrench};
use fs_query::triangle::quadratic::{QuadraticQueryControls, sweep::QuadraticSweepControls};
use fs_solid::hyper3::mixed::{
    MixedUpSettings,
    visco::{
        MixedViscoLawField, MixedViscoTetSolid,
        dynamics::{FiniteViscoContactDriftControls, MixedViscoTetMotion},
    },
};
use fs_solid::linear3::tet10::promote_p1_to_p2;
use fs_solver::nonlinear::NewtonKrylovConfig;

const RADIUS: f64 = 0.125;
const GAP: f64 = 1.0 / 262_144.0;
const STEP: f64 = 1.0 / 1024.0;
const DURATION: f64 = 4.0 * STEP;
const G: f64 = 8.0;

fn material() -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6.0,
                c01: 0.0,
                kappa: 120.0,
            },
            4.0,
        )
        .unwrap(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 36.0 },
            tau_s: 0.07,
        }],
    )
    .unwrap()
}

fn settings_for(step: f64) -> MixedViscoSustainedSphereTrajectoryControls {
    let mut impact = mixed_controls();
    impact.geometry_tolerance_m = 2e-5;
    let distance = QuadraticQueryControls {
        distance_tolerance: 1e-10,
        witness_tolerance: 1e-5,
        ..Default::default()
    };
    MixedViscoSustainedSphereTrajectoryControls {
        maximum_step_s: step,
        maximum_steps: 128,
        maximum_rejections: 512,
        path: QuadraticSolidPathControls {
            initial: QuadraticSolidSurfaceControls {
                surface: FiniteSolidSurfaceControls {
                    activation_tolerance_m: 1e-5,
                    ..Default::default()
                },
                query: distance,
                maximum_witness_span_m: 3e-5,
            },
            sweep: QuadraticSweepControls {
                distance,
                ..Default::default()
            },
        },
        volume: Default::default(),
        impact,
        settings: MixedUpSettings {
            load_steps: 1,
            maximum_newton_iterations: 30,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-11,
                relative_tolerance: 0.0,
                ..Default::default()
            },
            length_scale_m: 1.0,
            stress_scale_pa: 12.0,
            maximum_assembly_bytes: 2 * 1024 * 1024,
            maximum_factor_bytes: 2 * 1024 * 1024,
        },
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

fn relaxed<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 10]],
    law: &'a ViscoHyperelastic,
    cx: &Cx<'_>,
) -> MixedViscoTetMotion<'a> {
    let mut solid = MixedViscoTetSolid::try_new(
        nodes,
        tets,
        MixedViscoLawField::Uniform(law),
        VolumetricTreatment::OriginalMaterial,
        TetAssemblyBudget::standard(),
        4,
        64,
        cx,
    )
    .unwrap();
    // A genuine saddle x=X+(Y-1/4)(Z-1/4)/4, with det F=1 everywhere.
    // At the material point (0,1/4,1/4), position and outward normal are
    // exactly (0,1/4,1/4) and -e_x. No corner plane supplies that tangent.
    let held: Vec<_> = nodes
        .iter()
        .flat_map(|p| [0.25 * (p[1] - 0.25) * (p[2] - 0.25), 0.0, 0.0])
        .enumerate()
        .collect();
    let c = settings_for(STEP);
    for dt in [0.0, 0.007] {
        solid = solid
            .advance(
                dt,
                &held,
                &[],
                [0.0; 3],
                None,
                c.settings,
                FiniteViscoStepControls {
                    equilibrium_work_tolerance_j: 1e-9,
                    maximum_endpoint_work_defect_j: 1.0,
                },
                cx,
            )
            .unwrap()
            .solid;
    }
    assert!(solid.quadrature_states().iter().any(|s| s.dissipated > 0.0));
    MixedViscoTetMotion::try_new(solid, &[6.0], &vec![[0.0; 3]; nodes.len()], 128 * 1024, cx)
        .unwrap()
}

fn descriptor(id: u64, rigid_body: usize, friction: f64) -> FiniteSolidSphere {
    FiniteSolidSphere {
        id,
        rigid_body,
        center_body_m: Vec3::ZERO,
        radius_m: RADIUS,
        friction_coefficient: friction,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 1.0 },
    }
}

fn run_chain<'a>(
    motion: &MixedViscoTetMotion<'a>,
    c: MixedViscoSustainedSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<
    MixedViscoSustainedSphereTrajectoryCandidate<'a>,
    MixedViscoSustainedSphereTrajectoryError,
> {
    let surface = QuadraticSolidContactSurface::for_mixed_visco(motion, 4, cx)
        .map_err(MixedViscoSustainedSphereTrajectoryError::Surface)?;
    let bodies = [-RADIUS - GAP, -3.0 * RADIUS - 2.0 * GAP].map(|x| {
        RigidBodyState::new(
            Pose::new(Vec3::new(x, 0.25, 0.25), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap()
    });
    // Identity order deliberately differs from the original body indexing.
    let scene = FiniteSolidSphereScene::new(
        &[descriptor(9, 0, 0.0), descriptor(3, 1, 0.0)],
        &[SpherePairLaw {
            first: 3,
            second: 9,
            friction_coefficient: 0.0,
            restitution: 1.0,
        }],
        1,
    )
    .unwrap();
    mixed_visco_sustained_sphere_trajectory_candidate(
        motion,
        &surface,
        &scene,
        &bodies,
        &[properties(); 2],
        &(0..3 * motion.velocity_mps().len()).collect::<Vec<_>>(),
        &[],
        [0.0; 3],
        &[Wrench::ZERO; 2],
        Gravity::new(Vec3::new(G, 0.0, 0.0)).unwrap(),
        DURATION,
        ORIGIN,
        c,
        cx,
    )
}

fn replay(
    law: &ViscoHyperelastic,
    before: &MixedViscoTetMotion<'_>,
    after: &MixedViscoSustainedSphereTrajectoryCandidate<'_>,
) {
    let mut states = before.solid().quadrature_states().to_vec();
    let mut elapsed = 0.0;
    for step in &after.report.steps {
        for state in &mut states {
            *state = law
                .step_remainder_response(state, &state.deformation_gradient, step.duration_s)
                .unwrap()
                .state;
        }
        elapsed += step.duration_s;
        near(step.elapsed_s, elapsed, 0.0);
    }
    assert_eq!(after.solid.solid().quadrature_states(), states);
    assert_eq!(
        after.solid.solid().displacement_m(),
        before.solid().displacement_m()
    );
    assert_eq!(
        after.solid.solid().pressure_pa(),
        before.solid().pressure_pa()
    );
    assert_eq!(after.solid.reference_density_kg_m3(), &[6.0]);
    near(elapsed, DURATION, 0.0);
    near(
        after.solid.solid().time_s(),
        before.solid().time_s() + DURATION,
        2e-15,
    );
}

fn balanced(r: &MixedViscoSustainedSphereTrajectoryReport) {
    near(norm(r.linear_momentum_residual_kg_m_s), 0.0, 1e-7);
    near(norm(r.angular_momentum_residual_kg_m2_s), 0.0, 1e-7);
    near(r.work_energy_residual_j, 0.0, 1e-7);
    assert!(r.accumulated_linear_residual_kg_m_s < 1e-7);
    assert!(r.accumulated_angular_residual_kg_m2_s < 1e-7);
    assert!(r.accumulated_energy_residual_j < 1e-7);
}

#[test]
fn curved_loaded_chain_uses_pre_force_restitution_and_replays_original_maxwell_history() {
    with_cx(|_, cx| {
        let (nodes, tets) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = material();
        let motion = relaxed(&nodes, &tets, &law, cx);
        let original = motion.solid().quadrature_states().to_vec();
        let mut losses = Vec::new();
        for (dt, count) in [(STEP, 4), (0.5 * STEP, 8)] {
            let event = run_chain(&motion, settings_for(dt), cx).unwrap();
            assert_eq!(event.report.steps.len(), count);
            assert_eq!(event.report.rejected_steps, 0);
            assert!(matches!(
                event.reference_volume.verdict,
                QuadraticSolidVolumePathVerdict::Positive { .. }
            ));
            assert!(matches!(
                event.initial_volume.verdict,
                QuadraticSolidVolumePathVerdict::Positive { .. }
            ));
            for (body, x) in event
                .rigid_bodies
                .iter()
                .zip([-RADIUS - GAP, -3.0 * RADIUS - 2.0 * GAP])
            {
                near(
                    norm(body.pose().position_world().sub(Vec3::new(x, 0.25, 0.25))),
                    0.0,
                    2e-9,
                );
                near(norm(body.linear_momentum_world()), 0.0, 2e-9);
                near(norm(body.angular_momentum_body()), 0.0, 2e-9);
            }
            for step in &event.report.steps {
                near(step.duration_s, dt, 0.0);
                assert_eq!(step.path.volume.tested_tetrahedra, 1);
                assert!(matches!(
                    step.path.volume.verdict,
                    QuadraticSolidVolumePathVerdict::Positive { .. }
                ));
                assert_eq!(step.path.contact.sphere_pairs_checked, 1);
                assert!(matches!(step.path.contact.verdict,
                    QuadraticSphereScenePathVerdict::Clear { minimum_clearance_m } if minimum_clearance_m >= 0.0));
                assert!(step.maximum_contact_gap_m <= 1e-5);
                assert_eq!(step.normal_targets_m_s.len(), 2);
                for (&speed, &target) in step
                    .incoming_normal_speeds_m_s
                    .iter()
                    .zip(&step.normal_targets_m_s)
                {
                    near(speed, 0.0, 2e-9);
                    near(target, 0.0, 2e-9); // Restitution one must not bounce a resting load.
                }
                let impact = step.impact.as_ref().unwrap();
                assert_eq!(impact.contact.impulses.len(), 2);
                for (key, impulse) in step.contact_keys.iter().zip(&impact.contact.impulses) {
                    let expected = match key {
                        FiniteSolidSphereSceneKey::Solid { sphere: 9, .. } => 4.0 * G * dt,
                        FiniteSolidSphereSceneKey::Spheres {
                            first: 3,
                            second: 9,
                        } => 2.0 * G * dt,
                        other => panic!("unexpected contact row {other:?}"),
                    };
                    near(impulse[0], expected, 2e-9);
                    near(impulse[1], 0.0, 2e-9);
                    near(impulse[2], 0.0, 2e-9);
                }
            }
            let r = &event.report;
            let support = Vec3::new(-4.0 * G * DURATION, 0.0, 0.0);
            let moment = Vec3::new(0.0, 0.25, 0.25)
                .sub(vector(ORIGIN))
                .cross(support);
            near(norm(r.support_impulse_kg_m_s.sub(support)), 0.0, 2e-9);
            near(
                norm(r.support_angular_impulse_kg_m2_s.sub(moment)),
                0.0,
                2e-9,
            );
            near(norm(r.external_impulse_kg_m_s.add(support)), 0.0, 2e-9);
            near(r.source_work_j, 0.0, 2e-9);
            near(r.support_work_j, 0.0, 2e-9);
            near(r.kinetic_change_j, 0.0, 2e-9);
            near(r.stored_change_j + r.relaxation_loss_j, 0.0, 2e-9);
            assert!(r.relaxation_loss_j > 0.0);
            balanced(r);
            // Two masses of 2 kg: sum .5*m*g²*T*h is numerical support split
            // loss. It halves under refinement and is separate from relaxation.
            let expected_loss = 2.0 * G * G * DURATION * dt;
            near(r.unallocated_contact_loss_j, expected_loss, 2e-9);
            near(r.rigid_force_splitting_work_j, expected_loss, 2e-9);
            near(r.integration_defect_j, expected_loss, 2e-9);
            losses.push(r.unallocated_contact_loss_j);
            replay(&law, &motion, &event);
        }
        near(losses[1], 0.5 * losses[0], 2e-9);
        assert_eq!(motion.solid().quadrature_states(), original);
        near(motion.solid().time_s(), 0.007, 0.0);
    });
}

#[test]
fn quadratic_split_refines_complete_trials_and_late_refusal_preserves_pressure_and_memory() {
    with_cx(|gate, cx| {
        let (nodes, tets) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = material();
        let motion = relaxed(&nodes, &tets, &law, cx);
        let states = motion.solid().quadrature_states().to_vec();
        let pressure = motion.solid().pressure_pa().to_vec();
        let stress = motion.solid().first_piola_pa().to_vec();
        let coarse_loss = 2.0 * G * G * DURATION * STEP;
        let mut tight = settings_for(STEP);
        tight.integration_tolerance_j = coarse_loss / 7.0;
        let refined = run_chain(&motion, tight, cx).unwrap();
        assert!(refined.report.rejected_steps > 0);
        assert!(refined.report.steps.iter().all(|s| s.duration_s < STEP));
        assert!(refined.report.absolute_integration_defect_j <= tight.integration_tolerance_j);
        assert!(refined.report.unallocated_contact_loss_j < 0.2 * coarse_loss);
        near(
            refined.report.support_impulse_kg_m_s.x,
            -4.0 * G * DURATION,
            2e-9,
        );
        replay(&law, &motion, &refined);
        let mut limited = tight;
        limited.maximum_rejections = 0;
        assert!(matches!(
            run_chain(&motion, limited, cx),
            Err(MixedViscoSustainedSphereTrajectoryError::RejectionBudget { elapsed_s: 0.0, .. })
        ));
        let mut late = settings_for(STEP);
        late.maximum_steps = 1;
        assert!(matches!(run_chain(&motion, late, cx),
            Err(MixedViscoSustainedSphereTrajectoryError::Path(FiniteSolidSphereTrajectoryError::StepBudget { elapsed_s })) if elapsed_s > 0.0));
        assert_eq!(motion.solid().quadrature_states(), states);
        assert_eq!(motion.solid().pressure_pa(), pressure);
        assert_eq!(motion.solid().first_piola_pa(), stress);
        near(motion.solid().time_s(), 0.007, 0.0);
        gate.request();
        assert!(run_chain(&motion, tight, cx).is_err());
        assert_eq!(motion.solid().quadrature_states(), states);
    });
}

#[test]
fn curved_sliding_support_transfers_linear_and_full_tensor_spin_momentum() {
    with_cx(|_, cx| {
        let (nodes, tets) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = material();
        let motion = relaxed(&nodes, &tets, &law, cx);
        let surface = QuadraticSolidContactSurface::for_mixed_visco(&motion, 4, cx).unwrap();
        let gap = 1.0 / 4_194_304.0;
        let dt = 1.0 / 131_072.0;
        let mu = 0.25;
        let body = RigidBodyState::new(
            Pose::new(
                Vec3::new(-RADIUS - gap, 0.25, 0.25),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::new(0.0, 0.0, 2.0),
            Vec3::ZERO,
        )
        .unwrap();
        let scene = FiniteSolidSphereScene::new(&[descriptor(7, 0, mu)], &[], 0).unwrap();
        let result = mixed_visco_sustained_sphere_trajectory_candidate(
            &motion,
            &surface,
            &scene,
            &[body],
            &[properties()],
            &(0..3 * nodes.len()).collect::<Vec<_>>(),
            &[],
            [0.0; 3],
            &[Wrench::ZERO],
            Gravity::new(Vec3::new(G, 0.0, 0.0)).unwrap(),
            dt,
            ORIGIN,
            settings_for(dt),
            cx,
        )
        .unwrap();
        assert_eq!(result.report.steps.len(), 1);
        assert_eq!(result.report.rejected_steps, 0);
        let after = result.rigid_bodies[0];
        let normal = 2.0 * G * dt;
        let tangent = mu * normal;
        let spin_y = RADIUS * tangent;
        near(
            norm(
                after
                    .linear_momentum_world()
                    .sub(Vec3::new(0.0, 0.0, 2.0 - tangent)),
            ),
            0.0,
            2e-9,
        );
        near(
            norm(
                after
                    .pose()
                    .orientation()
                    .rotate_body_to_world(after.angular_momentum_body())
                    .sub(Vec3::new(0.0, spin_y, 0.0)),
            ),
            0.0,
            2e-10,
        );
        // The original off-diagonal inertia produces omega_x even though the
        // applied world angular impulse is purely along y.
        let omega = properties()
            .angular_velocity_body_checked(after.angular_momentum_body())
            .unwrap();
        near(omega.x, -5.0 * spin_y / 12.0, 2e-10);
        near(omega.y, 25.0 * spin_y / 12.0, 2e-10);
        let r = &result.report;
        near(
            norm(r.external_impulse_kg_m_s.sub(Vec3::new(normal, 0.0, 0.0))),
            0.0,
            2e-10,
        );
        near(
            norm(
                r.support_impulse_kg_m_s
                    .sub(Vec3::new(-normal, 0.0, -tangent)),
            ),
            0.0,
            2e-9,
        );
        near(
            norm(
                r.linear_momentum_change_kg_m_s
                    .sub(Vec3::new(0.0, 0.0, -tangent)),
            ),
            0.0,
            2e-9,
        );
        // Analytic contact mobility in the sliding z direction is
        // 1/m + r²*(I^-1)_yy = 409/768 kg^-1. Normal/tangent decouple.
        let split_loss = G * G * dt * dt;
        let contact_loss = split_loss + tangent - 0.5 * (409.0 / 768.0) * tangent * tangent;
        near(r.unallocated_contact_loss_j, contact_loss, 2e-9);
        near(r.rigid_force_splitting_work_j, split_loss, 2e-10);
        near(r.source_work_j, 0.0, 2e-10);
        near(r.support_work_j, 0.0, 2e-10);
        // Actual paired witnesses retain the positive activation gap. Its tiny
        // tangential couple remains a measured angular residual, not erased by
        // moving either application point to fabricate exact conservation.
        near(r.angular_momentum_residual_kg_m2_s.y, -gap * tangent, 5e-13);
        balanced(r);
        assert!(r.relaxation_loss_j > 0.0);
        assert!(result.report.steps[0].maximum_contact_gap_m > 0.0);
    });
}
