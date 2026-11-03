//! G1/G3 native impact localization, Maxwell replay and independent loss accounts.
use super::*;
use fs_couple::finite_solid::contact::step::FiniteSolidSphereTrajectoryError;
use fs_couple::inelastic_impact::{surface::*, trajectory::*};
use fs_mbd::{Gravity, Wrench};
use fs_solid::hyper3::visco::dynamics::FiniteViscoMotionControls;

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
fn trajectory_controls() -> ViscoSphereTrajectoryControls {
    let mut impact = mixed_controls();
    impact.geometry_tolerance_m = 2e-8;
    ViscoSphereTrajectoryControls {
        maximum_step_s: 0.3,
        maximum_steps: 128,
        maximum_rejections: 512,
        path: Default::default(),
        surface: FiniteSolidSurfaceControls {
            activation_tolerance_m: 1e-8,
            ..Default::default()
        },
        impact,
        settings: settings(),
        solid: FiniteViscoMotionControls {
            material: FiniteViscoStepControls {
                equilibrium_work_tolerance_j: 1e-9,
                maximum_endpoint_work_defect_j: 1e-2,
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
fn relaxed<'a>(law: &'a ViscoHyperelastic, cx: &Cx<'_>) -> FiniteViscoTetMotion<'a> {
    let mut solid = FiniteViscoTetSolid::try_new(
        &NODES,
        &TETS,
        ViscoelasticLawField::Uniform(law),
        TetAssemblyBudget::default(),
        cx,
    )
    .unwrap();
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
                trajectory_controls().solid.material,
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
    held: &[usize],
    z: f64,
    duration: f64,
    controls: ViscoSphereTrajectoryControls,
    cx: &Cx<'_>,
) -> Result<ViscoSphereTrajectoryCandidate<'a>, ViscoSphereTrajectoryError> {
    let surface = FiniteSolidContactSurface::for_visco(motion, 4, cx)
        .map_err(ViscoSphereTrajectoryError::Surface)?;
    let body = RigidBodyState::new(
        Pose::new(Vec3::new(0.25, 0.25, z), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.0, 0.0, 2.0),
        Vec3::ZERO,
    )
    .unwrap();
    let sphere = FiniteSolidSphere {
        id: 7,
        rigid_body: 0,
        center_body_m: Vec3::ZERO,
        radius_m: 0.125,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
    };
    visco_sphere_trajectory_candidate(
        motion,
        &surface,
        &[body],
        &[properties()],
        &[sphere],
        held,
        &[],
        [0.0; 3],
        &[Wrench::ZERO],
        Gravity::ZERO,
        duration,
        ORIGIN,
        controls,
        cx,
    )
}

#[test]
fn localized_bounce_replays_only_accepted_relaxation_and_keeps_loss_channels_separate() {
    with_cx(|gate, cx| {
        let law = material();
        let motion = relaxed(&law, cx);
        let original = motion.solid().element_states().to_vec();
        let fixed: Vec<_> = (0..12).collect();
        let c = trajectory_controls();
        let event = run(&motion, &fixed, -0.625, 1.0, c, cx).unwrap();
        // A stationary z=0 face: exact zero-gap t=.5; e=.5 gives J=3,
        // final z=-.375 and loss .5*m*(1-e^2)=.75 J. Activation stays explicit.
        near(
            event.rigid_bodies[0].pose().position_world().z,
            -0.375,
            2.1e-8,
        );
        near(event.rigid_bodies[0].linear_momentum_world().z, -1.0, 2e-10);
        near(event.solid.solid().time_s(), 1.007, 3e-15);
        assert!(event.report.rejected_steps > 0);
        let impacts: Vec<_> = event
            .report
            .steps
            .iter()
            .filter_map(|s| s.impact.as_ref())
            .filter(|i| i.contact.impulses.iter().any(|j| j[0] > 1e-8))
            .collect();
        assert_eq!(impacts.len(), 1);
        near(impacts[0].contact.impulses[0][0], 3.0, 2e-10);
        near(event.report.unallocated_contact_loss_j, 0.75, 2e-10);
        near(event.report.kinetic_change_j, -0.75, 2e-10);
        near(event.report.support_impulse_kg_m_s.z, -3.0, 2e-10);
        near(event.report.numerical_kinetic_loss_j, 0.0, 0.0);
        near(event.report.work_energy_residual_j, 0.0, 2e-10);
        assert!(event.report.relaxation_loss_j > 0.0);
        let mut replay = original[0].clone();
        for step in &event.report.steps {
            assert!(step.path.minimum_clearance_m >= 0.0);
            assert_eq!(step.path.tetrahedra_checked, 1);
            assert!(step.maximum_contact_gap_m <= c.surface.activation_tolerance_m);
            assert!(step.rigid[0].state_before.pose().position_world().z <= -0.125);
            assert!(step.rigid[0].state_after.pose().position_world().z <= -0.125);
            let exact_owner = law
                .step_response(&replay, &replay.deformation_gradient, step.duration_s)
                .unwrap();
            near(
                step.material.dissipated_j,
                exact_owner.dissipated_pa / 6.0,
                1e-16,
            );
            replay = exact_owner.state;
        }
        assert_eq!(event.solid.solid().element_states(), &[replay]);
        assert_eq!(
            event.solid.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        assert_eq!(
            event.solid.reference_nodal_density_kg_m3(),
            Some(RHO.as_slice())
        );
        let repeat = run(&motion, &fixed, -0.625, 1.0, c, cx).unwrap();
        assert_eq!(event.rigid_bodies, repeat.rigid_bodies);
        assert_eq!(
            event.solid.solid().element_states(),
            repeat.solid.solid().element_states()
        );
        let mut limited = c;
        limited.maximum_rejections = 0;
        assert!(matches!(run(&motion, &fixed, -0.625, 1.0, limited, cx),
            Err(ViscoSphereTrajectoryError::Path(FiniteSolidSphereTrajectoryError::RejectionBudget { elapsed_s, .. }))
            if elapsed_s > 0.0));
        assert_eq!(motion.solid().element_states(), original);
        assert_eq!(motion.solid().time_s(), 0.007);
        gate.request();
        assert!(run(&motion, &fixed, -0.625, 1.0, c, cx).is_err());
    });
}

#[test]
fn impact_recoil_advances_native_deformation_and_bounds_numerical_angular_drift() {
    with_cx(|_, cx| {
        let law = material();
        let motion = relaxed(&law, cx);
        let original = motion.solid().element_states().to_vec();
        let mut c = trajectory_controls();
        c.maximum_step_s = 5e-6;
        let z = -0.125 - 0.5 * c.surface.activation_tolerance_m;
        let event = run(&motion, &[], z, 1e-5, c, cx).unwrap();
        assert!(
            event.report.steps[0]
                .impact
                .as_ref()
                .unwrap()
                .contact
                .impulses[0][0]
                > 0.0
        );
        assert_ne!(
            event.solid.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        assert_ne!(event.solid.solid().element_states(), original);
        assert!(event.report.numerical_kinetic_loss_j > 0.0);
        near(
            norm(event.report.linear_momentum_residual_kg_m_s),
            0.0,
            1e-8,
        );
        near(
            norm(event.report.angular_momentum_residual_kg_m2_s),
            0.0,
            1e-8,
        );
        near(event.report.work_energy_residual_j, 0.0, 1e-8);
        let drift = event.report.absolute_angular_integration_defect_kg_m2_s;
        assert!(drift > 0.0);
        let mut strict = c;
        strict.angular_integration_tolerance_kg_m2_s = drift * 0.5;
        assert!(matches!(
            run(&motion, &[], z, 1e-5, strict, cx),
            Err(ViscoSphereTrajectoryError::Audit(_))
        ));
        // The ordinary original advance continues the final accepted branches.
        let next = event
            .solid
            .advance(1e-6, &[], &[], [0.0; 3], settings(), c.solid, cx)
            .unwrap();
        let expected = law
            .step_response(
                &event.solid.solid().element_states()[0],
                &next.motion.solid().element_states()[0].deformation_gradient,
                1e-6,
            )
            .unwrap();
        assert_eq!(next.motion.solid().element_states(), &[expected.state]);
        assert_eq!(motion.solid().element_states(), original);
        assert_eq!(motion.velocity_mps(), [[0.0; 3]; 4]);
    });
}
