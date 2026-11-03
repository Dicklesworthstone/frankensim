//! Current curved-boundary witnesses through original inertia and material state.

use super::*;
use fs_couple::inelastic_impact::surface::quadratic::{
    QuadraticSolidContactSurface, QuadraticSolidSurfaceControls,
};
use fs_couple::inelastic_impact::surface::{FiniteSolidSphere, FiniteSolidSurfaceControls};
use fs_material::hyper::discrete_gradient::DiscreteGradientSettings;
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_query::triangle::quadratic::QuadraticQueryControls;
use fs_solid::hyper3::mixed::conservative::{ConservativeFormulation, ConservativeTetDefinition};
use fs_solid::hyper3::mixed::visco::{
    MixedViscoLawField, MixedViscoTetSolid, dynamics::MixedViscoTetMotion,
};
use fs_solid::hyper3::mixed::{MixedLawField, MixedUpSettings};
use fs_solid::linear3::tet10::promote_p1_to_p2;
use fs_solver::nonlinear::NewtonKrylovConfig;

fn geometry_controls() -> QuadraticSolidSurfaceControls {
    QuadraticSolidSurfaceControls {
        surface: FiniteSolidSurfaceControls {
            activation_tolerance_m: 1e-6,
            maximum_contacts: 8,
            ..Default::default()
        },
        query: QuadraticQueryControls {
            distance_tolerance: 1e-9,
            witness_tolerance: 1e-5,
            ..Default::default()
        },
        maximum_witness_span_m: 3e-5,
    }
}

#[test]
fn detected_curved_contact_transfers_real_torque_through_both_native_solid_owners() {
    with_cx(|_, cx| {
        let (nodes, elements) = promote_p1_to_p2(&NODES, &TETS, None, cx).unwrap();
        // An exactly isochoric quadratic shear, F=[[1,.2z,.2y],[0,1,0],[0,0,1]].
        // Every corner is unchanged; the x=0 face is the actual saddle x=.2yz.
        let displacement: Vec<_> = nodes
            .iter()
            .map(|p| [0.2 * p[1] * p[2], 0.0, 0.0])
            .collect();
        let velocity = vec![[0.0; 3]; nodes.len()];
        let elastic = Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6.0,
                c01: 0.0,
                kappa: 120.0,
            },
            4.0,
        )
        .unwrap();
        let density = [6.0];
        let definition = ConservativeTetDefinition {
            nodes_m: &nodes,
            elements: &elements,
            material: MixedLawField::Uniform(&elastic),
            reference_density_kg_m3: &density,
            formulation: ConservativeFormulation::ExactIncompressible {
                maximum_relative_volume_residual: 1e-9,
                maximum_initial_relative_volume_rate_per_s: 1e-10,
                maximum_point_volume_work_defect: 1e-12,
            },
            budget: TetAssemblyBudget::standard(),
            quadrature_points_per_dir: 4,
            maximum_quadrature_points: 64,
            maximum_mass_assembly_bytes: 128 * 1024,
            material_settings: DiscreteGradientSettings {
                minimum_relative_strain_increment: 1e-6,
                absolute_work_tolerance_pa: 1e-10,
                relative_work_tolerance: 1e-12,
            },
        };
        let conservative = definition
            .prepare_with_pressure(&displacement, &velocity, &[0.3; 4], cx)
            .unwrap();
        let expected = Vec3::new(0.0125, 0.25, 0.25);
        let unscaled = Vec3::new(-1.0, 0.05, 0.05);
        let normal = unscaled.scale(1.0 / norm(unscaled));
        let radius = 0.05;
        let gap = 1e-7;
        let center_arm = Vec3::new(0.0, 0.0, 0.2);
        let center = expected.add(normal.scale(radius + gap));
        let rigid = RigidBodyState::new(
            Pose::new(center.sub(center_arm), UnitQuaternion::IDENTITY).unwrap(),
            normal.scale(-2.0),
            Vec3::ZERO,
        )
        .unwrap();
        let sphere = FiniteSolidSphere {
            id: 7,
            rigid_body: 0,
            center_body_m: center_arm,
            radius_m: radius,
            friction_coefficient: 0.0,
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
        };
        let surface = QuadraticSolidContactSurface::for_conservative(&conservative, 4, cx).unwrap();
        let event = surface
            .conservative_sphere_contact_candidate(
                &conservative,
                &[],
                &[rigid],
                &[properties()],
                &[sphere],
                None,
                ORIGIN,
                geometry_controls(),
                super::volume::volume_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(event.geometry.contacts.len(), 1);
        let witness = event.geometry.contacts[0];
        near(norm(witness.solid_point_world_m.sub(expected)), 0.0, 3e-5);
        assert!(witness.solid_point_world_m.x > 0.01); // corner-only face would be x=0
        assert!(witness.gap_lower_m >= 0.0);
        assert!(witness.gap_upper_m <= geometry_controls().surface.activation_tolerance_m);
        let report = &event.mechanical.report;
        let impulse = witness
            .contact
            .normal_world
            .scale(report.contact.impulses[0][0]);
        let torque = witness
            .rigid_point_world_m
            .sub(rigid.pose().position_world())
            .cross(impulse);
        assert!(norm(torque) > 1e-4);
        near(
            norm(
                event.mechanical.rigid_bodies[0]
                    .angular_momentum_body()
                    .sub(torque),
            ),
            0.0,
            2e-10,
        );
        near(
            norm(
                vector(report.solid_impulse.external.angular_kg_m2_s).sub(
                    witness
                        .solid_point_world_m
                        .sub(vector(ORIGIN))
                        .cross(impulse.scale(-1.0)),
                ),
            ),
            0.0,
            2e-10,
        );
        near(norm(report.linear_momentum_residual_kg_m_s), 0.0, 2e-10);
        near(norm(report.angular_momentum_residual_kg_m2_s), 0.0, 2e-10);
        near(report.work_energy_residual_j, 0.0, 2e-10);
        near(report.dissipation_residual_j, 0.0, 2e-10);
        near(report.solid_impulse.pressure_work_j, 0.0, 1e-11);
        assert!(report.solid_impulse.kkt_inertia.is_some());
        assert_eq!(
            event.mechanical.solid.displacement_m(),
            conservative.displacement_m()
        );
        assert_eq!(
            event.mechanical.solid.algorithmic_pressure_pa(),
            conservative.algorithmic_pressure_pa()
        );
        assert_eq!(
            event.mechanical.solid.stored_energy_j(),
            conservative.stored_energy_j()
        );
        assert_eq!(event.mechanical.solid.time_s(), conservative.time_s());

        let material = ViscoHyperelastic::new(
            elastic.clone(),
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 36.0 },
                tau_s: 0.07,
            }],
        )
        .unwrap();
        let mut accepted = MixedViscoTetSolid::try_new(
            &nodes,
            &elements,
            MixedViscoLawField::Uniform(&material),
            VolumetricTreatment::OriginalMaterial,
            TetAssemblyBudget::standard(),
            4,
            64,
            cx,
        )
        .unwrap();
        let prescribed: Vec<_> = displacement.iter().flatten().copied().enumerate().collect();
        for dt in [0.0, 0.007] {
            accepted = accepted
                .advance(
                    dt,
                    &prescribed,
                    &[],
                    [0.0; 3],
                    None,
                    MixedUpSettings {
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
                    FiniteViscoStepControls {
                        equilibrium_work_tolerance_j: 2e-9,
                        maximum_endpoint_work_defect_j: 1.0,
                    },
                    cx,
                )
                .unwrap()
                .solid;
        }
        assert!(
            accepted
                .quadrature_states()
                .iter()
                .any(|s| s.dissipated > 0.0)
        );
        let motion =
            MixedViscoTetMotion::try_new(accepted, &density, &velocity, 128 * 1024, cx).unwrap();
        let visco_surface = QuadraticSolidContactSurface::for_mixed_visco(&motion, 4, cx).unwrap();
        let visco = visco_surface
            .mixed_visco_sphere_contact_candidate(
                &motion,
                &[],
                &[rigid],
                &[properties()],
                &[sphere],
                ORIGIN,
                geometry_controls(),
                mixed_controls(),
                cx,
            )
            .unwrap();
        mixed_audit(&visco.mechanical.report);
        near(
            norm(visco.geometry.contacts[0].solid_point_world_m.sub(expected)),
            0.0,
            3e-5,
        );
        assert_eq!(
            visco.mechanical.solid.solid().quadrature_states(),
            motion.solid().quadrature_states()
        );
        assert_eq!(
            visco.mechanical.solid.solid().pressure_pa(),
            motion.solid().pressure_pa()
        );
        assert_eq!(
            visco.mechanical.solid.solid().time_s(),
            motion.solid().time_s()
        );

        // A raw penetration witness remains negative; the composed event refuses
        // it without modifying either accepted native material or rigid state.
        let penetrating = RigidBodyState::new(
            Pose::new(
                rigid.pose().position_world().sub(normal.scale(1e-4)),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            rigid.linear_momentum_world(),
            rigid.angular_momentum_body(),
        )
        .unwrap();
        let raw = visco_surface
            .mixed_visco_sphere_separations(
                &motion,
                &[penetrating],
                &[sphere],
                geometry_controls(),
                cx,
            )
            .unwrap();
        assert!(raw.closest[0].gap_m < 0.0);
        assert!(
            visco_surface
                .mixed_visco_sphere_contact_candidate(
                    &motion,
                    &[],
                    &[penetrating],
                    &[properties()],
                    &[sphere],
                    ORIGIN,
                    geometry_controls(),
                    mixed_controls(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(motion.velocity_mps(), velocity);
        assert_eq!(conservative.velocity_mps(), velocity);
    });
}
