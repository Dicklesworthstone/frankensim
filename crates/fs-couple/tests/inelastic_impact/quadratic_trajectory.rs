//! G1/G3 native curved bounce, retained Maxwell evolution and loaded momentum.

#[path = "quadratic_trajectory/scene.rs"]
mod scene;

use super::*;
use fs_couple::finite_solid::contact::step::FiniteSolidSphereTrajectoryError;
use fs_couple::inelastic_impact::{
    surface::{
        FiniteSolidSphere, FiniteSolidSurfaceControls,
        quadratic::{
            QuadraticSolidContactSurface, QuadraticSolidSurfaceControls,
            path::{QuadraticSolidPathControls, QuadraticSolidPathVerdict},
            volume::QuadraticSolidVolumePathVerdict,
        },
    },
    trajectory::quadratic::*,
};
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_mbd::{Gravity, Wrench};
use fs_query::triangle::quadratic::{QuadraticQueryControls, sweep::QuadraticSweepControls};
use fs_solid::hyper3::{
    mixed::{
        MixedUpSettings,
        visco::{MixedViscoLawField, MixedViscoTetSolid, dynamics::MixedViscoTetMotion},
    },
    visco::dynamics::FiniteViscoMotionControls,
};
use fs_solid::linear3::tet10::promote_p1_to_p2;
use fs_solver::nonlinear::NewtonKrylovConfig;

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

fn controls() -> MixedViscoSphereTrajectoryControls {
    let mut impact = mixed_controls();
    impact.geometry_tolerance_m = 2e-6;
    let distance = QuadraticQueryControls {
        distance_tolerance: 1e-9,
        witness_tolerance: 1e-5,
        ..Default::default()
    };
    MixedViscoSphereTrajectoryControls {
        maximum_step_s: 0.2,
        maximum_steps: 128,
        maximum_rejections: 512,
        path: QuadraticSolidPathControls {
            initial: QuadraticSolidSurfaceControls {
                surface: FiniteSolidSurfaceControls {
                    activation_tolerance_m: 1e-6,
                    maximum_contacts: 8,
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
        solid: FiniteViscoMotionControls {
            material: FiniteViscoStepControls {
                equilibrium_work_tolerance_j: 2e-9,
                maximum_endpoint_work_defect_j: 1.0,
            },
            maximum_numerical_kinetic_loss_j: 1e-3,
        },
        linear_tolerance_kg_m_s: 1e-7,
        angular_tolerance_kg_m2_s: 1e-7,
        energy_tolerance_j: 1e-7,
        integration_tolerance_j: 1e-3,
        angular_integration_tolerance_kg_m2_s: 1e-6,
    }
}

fn solid<'a>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    law: &'a ViscoHyperelastic,
    treatment: VolumetricTreatment,
    cx: &Cx<'_>,
) -> MixedViscoTetSolid<'a> {
    MixedViscoTetSolid::try_new(
        nodes,
        elements,
        MixedViscoLawField::Uniform(law),
        treatment,
        TetAssemblyBudget::standard(),
        4,
        64,
        cx,
    )
    .unwrap()
}

fn relaxed<'a>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    law: &'a ViscoHyperelastic,
    cx: &Cx<'_>,
) -> MixedViscoTetMotion<'a> {
    let mut state = solid(
        nodes,
        elements,
        law,
        VolumetricTreatment::OriginalMaterial,
        cx,
    );
    // Isochoric quadratic shear: x'=x+.2yz. All corners stay unchanged,
    // but the x=0 boundary is a genuine curved saddle with det F=1 everywhere.
    let held: Vec<_> = nodes
        .iter()
        .flat_map(|p| [0.2 * p[1] * p[2], 0.0, 0.0])
        .enumerate()
        .collect();
    for dt in [0.0, 0.007] {
        state = state
            .advance(
                dt,
                &held,
                &[],
                [0.0; 3],
                None,
                controls().settings,
                controls().solid.material,
                cx,
            )
            .unwrap()
            .solid;
    }
    assert!(state.quadrature_states().iter().any(|s| s.dissipated > 0.0));
    MixedViscoTetMotion::try_new(state, &[6.0], &vec![[0.0; 3]; nodes.len()], 128 * 1024, cx)
        .unwrap()
}

fn sphere(id: u64, rigid_body: usize) -> FiniteSolidSphere {
    FiniteSolidSphere {
        id,
        rigid_body,
        center_body_m: Vec3::ZERO,
        radius_m: 0.05,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
    }
}

fn curved_point_and_normal() -> (Vec3, Vec3) {
    let normal = Vec3::new(-1.0, 0.05, 0.05);
    (
        Vec3::new(0.0125, 0.25, 0.25),
        normal.scale(1.0 / norm(normal)),
    )
}

fn bounce<'a>(
    motion: &MixedViscoTetMotion<'a>,
    held: &[usize],
    gap: f64,
    duration: f64,
    c: MixedViscoSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<MixedViscoSphereTrajectoryCandidate<'a>, MixedViscoSphereTrajectoryError> {
    let surface = QuadraticSolidContactSurface::for_mixed_visco(motion, 4, cx)
        .map_err(MixedViscoSphereTrajectoryError::Surface)?;
    let (point, normal) = curved_point_and_normal();
    let body = RigidBodyState::new(
        Pose::new(
            point.add(normal.scale(0.05 + gap)),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        normal.scale(-2.0),
        Vec3::ZERO,
    )
    .unwrap();
    mixed_visco_sphere_trajectory_candidate(
        motion,
        &surface,
        &[body],
        &[properties()],
        &[sphere(7, 0)],
        held,
        &[],
        [0.0; 3],
        None,
        &[Wrench::ZERO],
        Gravity::ZERO,
        duration,
        ORIGIN,
        c,
        cx,
    )
}

fn assert_balanced(report: &MixedViscoSphereTrajectoryReport) {
    near(norm(report.linear_momentum_residual_kg_m_s), 0.0, 1e-7);
    near(norm(report.angular_momentum_residual_kg_m2_s), 0.0, 1e-7);
    near(report.work_energy_residual_j, 0.0, 1e-7);
    assert!(report.accumulated_linear_residual_kg_m_s < 1e-7);
    assert!(report.accumulated_angular_residual_kg_m2_s < 1e-7);
    assert!(report.accumulated_energy_residual_j < 1e-7);
}

#[test]
fn curved_bounce_shortens_native_paths_and_replays_only_accepted_maxwell_intervals() {
    with_cx(|gate, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = material();
        let motion = relaxed(&nodes, &elements, &law, cx);
        let original = motion.solid().quadrature_states().to_vec();
        let pressure = motion.solid().pressure_pa().to_vec();
        let held: Vec<_> = (0..3 * nodes.len()).collect();
        let c = controls();
        let result = bounce(&motion, &held, 0.1, 0.3, c, cx).unwrap();
        let (point, normal) = curved_point_and_normal();
        // Speed 1 m/s, exact zero-gap time .1 s, m=2 kg and e=.5 give
        // J=3 N.s and loss .75 J. Returned witnesses retain their authored
        // localization bound; activation remains a positive physical gap.
        near(
            norm(result.rigid_bodies[0].linear_momentum_world().sub(normal)),
            0.0,
            5e-4,
        );
        near(
            norm(
                result.rigid_bodies[0]
                    .pose()
                    .position_world()
                    .sub(point.add(normal.scale(0.15))),
            ),
            0.0,
            8e-5,
        );
        near(result.report.unallocated_contact_loss_j, 0.75, 1e-6);
        near(result.report.kinetic_change_j, -0.75, 1e-6);
        near(
            norm(result.report.support_impulse_kg_m_s.sub(normal.scale(3.0))),
            0.0,
            5e-4,
        );
        near(result.solid.solid().time_s(), 0.307, 2e-15);
        near(result.report.source_work_j, 0.0, 1e-13);
        near(result.report.support_work_j, 0.0, 1e-13);
        near(result.report.numerical_kinetic_loss_j, 0.0, 0.0);
        assert!(result.report.relaxation_loss_j > 0.0);
        assert!(result.report.rejected_steps > 0);
        assert_balanced(&result.report);
        assert!(matches!(
            result.reference_volume.verdict,
            QuadraticSolidVolumePathVerdict::Positive { .. }
        ));
        assert!(matches!(
            result.initial_volume.verdict,
            QuadraticSolidVolumePathVerdict::Positive { .. }
        ));
        let mut replay = original.clone();
        let mut impulses = 0;
        for interval in &result.report.steps {
            assert_eq!(interval.path.volume.tested_tetrahedra, 1);
            assert!(matches!(
                interval.path.volume.verdict,
                QuadraticSolidVolumePathVerdict::Positive { .. }
            ));
            assert!(matches!(interval.path.boundary.verdict,
                QuadraticSolidPathVerdict::Clear { minimum_clearance_m } if minimum_clearance_m >= 0.0));
            assert_eq!(interval.path.minimum_pair_clearance_m, None);
            assert!(
                interval.maximum_contact_gap_m <= c.path.initial.surface.activation_tolerance_m
            );
            if let Some(impact) = &interval.impact {
                if impact.contact.impulses.iter().any(|j| j[0] > 1e-8) {
                    impulses += 1;
                    near(impact.contact.impulses[0][0], 3.0, 1e-6);
                }
            }
            for state in &mut replay {
                *state = law
                    .step_remainder_response(
                        state,
                        &state.deformation_gradient,
                        interval.duration_s,
                    )
                    .unwrap()
                    .state;
            }
        }
        assert_eq!(impulses, 1);
        assert_eq!(result.solid.solid().quadrature_states(), replay);
        assert_eq!(
            result.solid.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        assert_eq!(result.solid.reference_density_kg_m3(), [6.0]);
        // A late refusal may have private accepted trials; none is published.
        let mut limited = c;
        limited.maximum_step_s = 0.06;
        limited.maximum_rejections = 0;
        assert!(matches!(bounce(&motion, &held, 0.1, 0.3, limited, cx),
            Err(MixedViscoSphereTrajectoryError::RejectionBudget { elapsed_s, .. }) if elapsed_s > 0.0));
        assert_eq!(motion.solid().quadrature_states(), original);
        assert_eq!(motion.solid().pressure_pa(), pressure);
        assert_eq!(motion.solid().time_s(), 0.007);
        gate.request();
        assert!(bounce(&motion, &held, 0.1, 0.3, c, cx).is_err());
    });
}

#[test]
fn curved_free_recoil_continues_original_pressure_and_relaxation_without_reconstruction() {
    with_cx(|_, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = material();
        let motion = relaxed(&nodes, &elements, &law, cx);
        let original = motion.solid().quadrature_states().to_vec();
        let mut c = controls();
        c.maximum_step_s = 5e-6;
        let gap = 0.5 * c.path.initial.surface.activation_tolerance_m;
        let result = bounce(&motion, &[], gap, 1e-5, c, cx).unwrap();
        assert!(
            result.report.steps[0]
                .impact
                .as_ref()
                .unwrap()
                .contact
                .impulses[0][0]
                > 0.0
        );
        assert_ne!(
            result.solid.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        assert_ne!(result.solid.solid().quadrature_states(), original);
        assert!(
            result
                .solid
                .solid()
                .pressure_pa()
                .iter()
                .all(|p| p.is_finite())
        );
        assert!(result.report.numerical_kinetic_loss_j > 0.0);
        assert_balanced(&result.report);
        let mut strict = c;
        strict.integration_tolerance_j = 0.5 * result.report.absolute_integration_defect_j;
        assert!(matches!(
            bounce(&motion, &[], gap, 1e-5, strict, cx),
            Err(MixedViscoSphereTrajectoryError::Audit(_))
        ));
        let next = result
            .solid
            .advance(1e-6, &[], &[], [0.0; 3], None, c.settings, c.solid, cx)
            .unwrap();
        for (before, after) in result
            .solid
            .solid()
            .quadrature_states()
            .iter()
            .zip(next.motion.solid().quadrature_states())
        {
            let expected = law
                .step_remainder_response(before, &after.deformation_gradient, 1e-6)
                .unwrap();
            assert_eq!(*after, expected.state);
        }
        assert_eq!(motion.solid().quadrature_states(), original);
        assert!(motion.velocity_mps().iter().all(|v| *v == [0.0; 3]));
        near(next.motion.solid().time_s(), 0.007011, 2e-17);
    });
}

#[test]
fn original_signed_body_load_and_full_rigid_wrench_close_the_complete_trajectory() {
    with_cx(|_, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = material();
        let state = solid(
            &nodes,
            &elements,
            &law,
            VolumetricTreatment::OriginalMaterial,
            cx,
        );
        let motion = MixedViscoTetMotion::try_new(
            state,
            &[6.0],
            &vec![[0.0; 3]; nodes.len()],
            128 * 1024,
            cx,
        )
        .unwrap();
        let surface = QuadraticSolidContactSurface::for_mixed_visco(&motion, 4, cx).unwrap();
        let body = RigidBodyState::new(
            Pose::new(Vec3::new(-3.0, -3.0, -3.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(0.3, 0.2, 0.0),
            Vec3::new(0.01, 0.02, 0.03),
        )
        .unwrap();
        let wrench = Wrench {
            force_world: Vec3::new(0.4, -0.1, 0.2),
            torque_body: Vec3::new(0.05, 0.02, -0.03),
        };
        let gravity = Gravity::new(Vec3::new(0.0, 0.0, -0.7)).unwrap();
        let mut c = controls();
        c.maximum_step_s = 5e-5;
        let result = mixed_visco_sphere_trajectory_candidate(
            &motion,
            &surface,
            &[body],
            &[properties()],
            &[sphere(7, 0)],
            &[],
            &[],
            [6.0, -12.0, 3.0],
            None,
            &[wrench],
            gravity,
            1e-4,
            ORIGIN,
            c,
            cx,
        )
        .unwrap();
        // V=1/6, rho=6: the original consistent solid has mass 1 kg.
        // Its force is [1,-2,.5] N; the rigid adds its own wrench and 2g.
        let expected_impulse = Vec3::new(1.4e-4, -2.1e-4, -0.7e-4);
        near(
            norm(result.report.external_impulse_kg_m_s.sub(expected_impulse)),
            0.0,
            2e-12,
        );
        near(
            norm(
                result
                    .report
                    .linear_momentum_change_kg_m_s
                    .sub(expected_impulse),
            ),
            0.0,
            2e-10,
        );
        let moment_rate = Vec3::new(0.25, 0.25, 0.25)
            .sub(vector(ORIGIN))
            .cross(Vec3::new(1.0, -2.0, 0.5));
        for interval in &result.report.steps {
            near(
                norm(
                    vector(interval.momentum.applied.linear_kg_m_s)
                        .sub(Vec3::new(1.0, -2.0, 0.5).scale(interval.duration_s)),
                ),
                0.0,
                2e-12,
            );
            near(
                norm(
                    vector(interval.momentum.applied.angular_kg_m2_s)
                        .sub(moment_rate.scale(interval.duration_s)),
                ),
                0.0,
                2e-12,
            );
            assert!(interval.impact.is_none());
        }
        // Two BE intervals: sum .5*m*|a*dt|². This remains numerical loss,
        // separate from the actual source work and any physical relaxation.
        near(result.report.numerical_kinetic_loss_j, 1.3125e-8, 2e-12);
        near(result.report.unallocated_contact_loss_j, 0.0, 0.0);
        near(norm(result.report.support_impulse_kg_m_s), 0.0, 2e-12);
        assert!(result.report.source_work_j != 0.0);
        assert_balanced(&result.report);
    });
}

#[test]
fn unsupported_exact_volume_velocity_and_hidden_sphere_pair_crossing_refuse_transactionally() {
    with_cx(|_, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        let law = material();
        let state = solid(
            &nodes,
            &elements,
            &law,
            VolumetricTreatment::ExactIncompressible,
            cx,
        );
        let exact = MixedViscoTetMotion::try_new(
            state,
            &[6.0],
            &vec![[0.0; 3]; nodes.len()],
            128 * 1024,
            cx,
        )
        .unwrap();
        let surface = QuadraticSolidContactSurface::for_mixed_visco(&exact, 4, cx).unwrap();
        assert!(matches!(
            mixed_visco_sphere_trajectory_candidate(
                &exact,
                &surface,
                &[],
                &[],
                &[],
                &[],
                &[],
                [0.0; 3],
                None,
                &[],
                Gravity::ZERO,
                0.1,
                ORIGIN,
                controls(),
                cx
            ),
            Err(MixedViscoSphereTrajectoryError::ExactIncompressibleChordVelocity)
        ));
        assert_eq!(exact.solid().time_s(), 0.0);
        let motion = relaxed(&nodes, &elements, &law, cx);
        let original = motion.solid().quadrature_states().to_vec();
        let surface = QuadraticSolidContactSurface::for_mixed_visco(&motion, 4, cx).unwrap();
        let bodies = [(-3.0, 4.0), (-1.0, -4.0)].map(|(x, p)| {
            RigidBodyState::new(
                Pose::new(Vec3::new(x, -3.0, -3.0), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::new(p, 0.0, 0.0),
                Vec3::ZERO,
            )
            .unwrap()
        });
        let held: Vec<_> = (0..3 * nodes.len()).collect();
        let mut c = controls();
        c.maximum_step_s = 1.0;
        // Both endpoints separate; the two centres cross halfway through.
        assert!(
            matches!(mixed_visco_sphere_trajectory_candidate(&motion, &surface, &bodies,
            &[properties(); 2], &[sphere(7, 0), sphere(9, 1)], &held, &[], [0.0; 3], None,
            &[Wrench::ZERO; 2], Gravity::ZERO, 1.0, ORIGIN, c, cx),
            Err(MixedViscoSphereTrajectoryError::Path(FiniteSolidSphereTrajectoryError::UnsupportedSpherePair {
                first: 7, second: 9, fraction,
            })) if fraction > 0.0 && fraction < 0.5)
        );
        assert_eq!(motion.solid().quadrature_states(), original);
        assert_eq!(motion.solid().time_s(), 0.007);
    });
}
