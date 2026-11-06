//! Source natural strain reaches the actual conservative quadratic solid.
use super::*;
use fs_couple::finite_solid::thermal::conservative::{
    ConservativeThermalError, ConservativeThermalField, ConservativeThermalPreparation,
};
use fs_material::hyper::discrete_gradient::DiscreteGradientSettings;
use fs_material::hyper::discrete_gradient::bulk::DiscreteBulkSettings;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_solid::hyper3::mixed::MixedUpSettings;
use fs_solid::hyper3::mixed::conservative::ConservativeFormulation;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::promote_p1_to_p2;
use fs_solver::nonlinear::NewtonKrylovConfig;

fn preparation(formulation: ConservativeFormulation) -> ConservativeThermalPreparation {
    ConservativeThermalPreparation {
        formulation,
        budget: TetAssemblyBudget::standard(),
        quadrature_points_per_dir: 4,
        maximum_quadrature_points: 64,
        maximum_mass_assembly_bytes: 128 * 1024,
        material_settings: DiscreteGradientSettings {
            minimum_relative_strain_increment: 1e-6,
            absolute_work_tolerance_pa: 1e-10,
            relative_work_tolerance: 1e-12,
        },
    }
}

fn settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 30,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        length_scale_m: 1.0,
        stress_scale_pa: 15.0,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
    }
}

fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}

#[test]
fn source_reference_field_retains_order_original_density_and_exact_element_coverage() {
    with_cx(|_, cx| {
        let source = card(STRETCHES, "T", false);
        let hot = bind_expansion(&source, Q, 310.0);
        let warm = bind_expansion(&source, Q, 305.0);
        let manufacture =
            ReferenceDeformation::try_new([1.03, 0.16, 0.02, 0.0, 0.97, -0.01, 0.0, 0.0, 1.02])
                .unwrap();
        let thermal = match hot.thermal_reference().unwrap() {
            ThermalReference::Expansion(reference) => *reference.reference_deformation(),
            _ => panic!("expected authored directional source"),
        };
        let expected = manufacture.followed_by(thermal).unwrap();
        let reversed = thermal.followed_by(manufacture).unwrap();
        assert_ne!(
            expected.deformation_gradient(),
            reversed.deformation_gradient()
        );
        let uniform = ConservativeThermalField::uniform(hot.clone(), manufacture, 2, cx).unwrap();
        assert_eq!(uniform.element_count(), 2);
        assert_eq!(uniform.reference_at(0), Some(&expected));
        assert_eq!(uniform.reference_at(1), Some(&expected));
        assert_eq!(uniform.reference_before_thermal_at(1), Some(&manufacture));
        assert!(uniform.material_at(2).is_none());
        assert!(uniform.reference_at(2).is_none());
        assert!(uniform.reference_before_thermal_at(2).is_none());
        assert_eq!(uniform.material_at(0).unwrap(), &hot);
        assert_eq!(uniform.reference_density_kg_m3(), &[6.0, 6.0]);

        let references = [manufacture, ReferenceDeformation::default()];
        let field =
            ConservativeThermalField::per_element(vec![hot.clone(), warm.clone()], &references, cx)
                .unwrap();
        assert_eq!(field.material_at(0).unwrap(), &hot);
        assert_eq!(field.material_at(1).unwrap(), &warm);
        assert_eq!(field.reference_at(0), Some(&expected));
        assert_ne!(field.reference_at(0), field.reference_at(1));
        assert_eq!(field.reference_density_kg_m3(), &[6.0, 6.0]);
        assert!(ConservativeThermalField::uniform(hot.clone(), manufacture, 0, cx).is_err());
        assert!(ConservativeThermalField::per_element(vec![hot], &[], cx).is_err());
        let (nodes, elements) = mesh(cx);
        assert!(matches!(
            field.prepare(
                &nodes,
                &elements,
                &vec![[0.0; 3]; nodes.len()],
                &vec![[0.0; 3]; nodes.len()],
                &[],
                preparation(ConservativeFormulation::FullMaterial),
                cx
            ),
            Err(ConservativeThermalError::Invalid(_))
        ));
    });
}

#[test]
fn large_sourced_natural_expansion_moves_conservatively_with_original_mass_in_all_formulations() {
    with_cx(|_, cx| {
        // Existing authored source: 300->310 K has s=4 and rho0=6,
        // while endpoint density is 6/64. The original inventory must win.
        let source = super::super::thermal::thermal_card(0.3, 1.0, false);
        let material = super::super::thermal::prepared(&source);
        let reference = ReferenceDeformation::default();
        let field = ConservativeThermalField::uniform(material, reference, 1, cx).unwrap();
        let g = field.reference_at(0).unwrap();
        near(g.volume_ratio(), 64.0, 1e-12);
        assert!(
            field
                .material_at(0)
                .unwrap()
                .law()
                .law()
                .tangent(g.deformation_gradient())
                .is_err()
        );
        let (nodes, elements) = mesh(cx);
        let displacement: Vec<_> = nodes.iter().map(|x| x.map(|x| 3.0 * x)).collect();
        let speed = [0.25, -0.1, 0.05];
        let velocities = vec![speed; nodes.len()];
        let mass = 6.0 / 6.0;
        let kinetic = 0.5 * mass * speed.iter().map(|v| v * v).sum::<f64>();
        let exact = ConservativeFormulation::ExactIncompressible {
            maximum_relative_volume_residual: 1e-9,
            maximum_initial_relative_volume_rate_per_s: 1e-9,
            maximum_point_volume_work_defect: 1e-10,
        };
        let bulk = ConservativeFormulation::OriginalMaterial {
            bulk: DiscreteBulkSettings {
                minimum_volume_ratio: 0.5,
                maximum_volume_ratio: 1.5,
                maximum_iterations: 96,
                absolute_pressure_tolerance_pa: 1e-10,
                relative_pressure_tolerance: 1e-12,
                absolute_work_tolerance_pa: 1e-10,
                relative_work_tolerance: 1e-12,
            },
            maximum_relative_volume_residual: 1e-9,
            maximum_point_volume_work_defect: 1e-10,
        };
        for formulation in [ConservativeFormulation::FullMaterial, exact, bulk] {
            let mut body = field
                .prepare(
                    &nodes,
                    &elements,
                    &displacement,
                    &velocities,
                    &[],
                    preparation(formulation),
                    cx,
                )
                .unwrap();
            assert!(std::ptr::eq(body.source_field(), &field));
            near(body.native().reference_mass_kg(), mass, 2e-12);
            near(body.native().stored_energy_j(), 0.0, 1e-9);
            near(body.native().kinetic_energy_j(), kinetic, 1e-12);
            let acceleration = vec![[0.3, -0.2, 0.1]; nodes.len()];
            let load = body
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap();
            for axis in 0..3 {
                let total: f64 = load
                    .iter()
                    .filter(|(dof, _)| dof % 3 == axis)
                    .map(|(_, f)| f)
                    .sum();
                near(total, mass * acceleration[0][axis], 2e-12);
            }
            let mut time = 0.0;
            for dt in [0.01, 0.015] {
                let before_u = body.native().displacement_m().to_vec();
                let before_v = body.native().velocity_mps().to_vec();
                let before_t = body.native().time_s();
                let step = body
                    .advance(dt, &[], &[], [0.0; 3], settings(), 1e-8, cx)
                    .unwrap();
                assert_eq!(body.native().displacement_m(), before_u);
                assert_eq!(body.native().velocity_mps(), before_v);
                assert_eq!(body.native().time_s(), before_t);
                time += dt;
                for (node, (u, v)) in step
                    .solid
                    .native()
                    .displacement_m()
                    .iter()
                    .zip(step.solid.native().velocity_mps())
                    .enumerate()
                {
                    for axis in 0..3 {
                        near(u[axis], displacement[node][axis] + time * speed[axis], 2e-9);
                        near(v[axis], speed[axis], 2e-8);
                    }
                }
                near(step.inertia.numerical_loss_j, 0.0, 0.0);
                near(step.inertia.kinetic_after_j, kinetic, 1e-9);
                near(step.energy.balance_error_j, 0.0, 1e-8);
                near(step.energy.cumulative_balance_error_j, 0.0, 1e-8);
                assert!(std::ptr::eq(step.solid.source_field(), &field));
                assert_eq!(
                    step.solid
                        .reference_acceleration_load_n(&acceleration, cx)
                        .unwrap(),
                    load
                );
                near(step.solid.native().reference_mass_kg(), mass, 2e-12);
                body = step.solid;
            }
        }
    });
}
