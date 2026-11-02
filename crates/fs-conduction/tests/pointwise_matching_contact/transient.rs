//! Actual storage consumers retain the local contact rule and its full tangent.

use super::*;
use fs_blake3::ContentHash;
use fs_conduction::interface::temperature::ContactIterationConfig;
use fs_conduction::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::{LinearConfig, SolveConfig};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_solver::NewtonKrylovConfig;

fn transport(t: &[f64]) -> Vec<f64> {
    let mut result = independent_contact_load(t);
    // Both unit tetrahedra have volume 1/6 and this reference stiffness;
    // reflection in z changes gradients together, not their dot products.
    let stiffness = [
        [3., -1., -1., -1.],
        [-1., 1., 0., 0.],
        [-1., 0., 1., 0.],
        [-1., 0., 0., 1.],
    ];
    for base in [0, 4] {
        for a in 0..4 {
            for b in 0..4 {
                result[base + a] += 0.7 / 6. * stiffness[a][b] * t[base + b];
            }
        }
    }
    result
}

#[test]
fn capacity_newton_accepts_the_actual_pointwise_residual_and_preserves_closed_energy() {
    context(|_, cx| {
        let mesh = mesh(false);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = source_card(false);
        let set = contacts(cx, &mesh, &boundary, &card, false);
        let capacities = [
            VolumetricHeatCapacity::declared(10.).unwrap(),
            VolumetricHeatCapacity::declared(20.).unwrap(),
        ];
        let engine = BackwardEuler::per_element(cx, &mesh, &capacities).unwrap();
        let material = ConductivityModel::isotropic_declared(0.7).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let mut linear = SolveConfig::default().linear;
        linear.tolerance = 1e-12;
        linear.max_iterations = 2048;
        linear.restart = 8;
        let config = StepConfig {
            linear,
            energy_tolerance_j: 1e-8,
        };
        let nonlinear = NonlinearStepConfig {
            residual_rtol: 0.,
            residual_atol_j: 1e-11,
            ..Default::default()
        };
        let coupling = ContactIterationConfig {
            max_iterations: 64,
            relaxation: 0.5,
            residual_tolerance_j: 1e-9,
        };
        let old = [360., 300., 300., 340., 310., 325., 320., 315.];
        let before = old;
        let dt = 0.04;
        let accepted = set
            .advance_capacity_newton(
                cx, &engine, problem, &old, dt, config, nonlinear, coupling, 1,
            )
            .unwrap();
        assert_eq!(old, before);
        assert_eq!(accepted.physical_contacts.quadrature_states().len(), 3);
        assert!(accepted.physical_contacts.states().is_empty());
        let t = &accepted.conduction.temperature;
        let spatial = transport(t);
        let mut energy = 0.;
        for i in 0..8 {
            let storage = if i < 4 { 10. / 24. } else { 20. / 24. };
            near(storage * (t[i] - old[i]) + dt * spatial[i], 0., 1e-9);
            energy += storage * (t[i] - old[i]);
        }
        near(energy, 0., 1e-9);
        assert!(accepted.physical.residual_norm_j <= coupling.residual_tolerance_j);
        let legacy = TemperatureContacts::new(
            cx,
            &mesh,
            &boundary,
            vec![],
            vec![law(&mesh, &card, false, false)],
            1,
        )
        .unwrap();
        let mean = legacy
            .advance_capacity_newton(
                cx, &engine, problem, &old, dt, config, nonlinear, coupling, 1,
            )
            .unwrap();
        assert!(
            t.iter()
                .zip(&mean.conduction.temperature)
                .any(|(a, b)| (a - b).abs() > 1e-5),
            "transient must retain the changed local constitutive model"
        );
        let bound = engine
            .linearize_temperature_contact_endpoint(
                cx,
                problem,
                &set,
                &old,
                t,
                dt,
                config,
                coupling.residual_tolerance_j,
                &[],
                1,
            )
            .unwrap();
        let gradient = bound
            .pullback_with_resistance_controls(cx, &[0., 0., 0., 0., 1., 0., 0., 0.], &[], &[], 1)
            .unwrap();
        assert_eq!(gradient.resistance_multipliers.len(), 1);
        assert_eq!(gradient.resistance_multipliers[0].face_pairs, 1);
        assert!(gradient.resistance_multipliers[0].derivative.is_finite());
    });
}

fn chart() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([37; 32]),
        [
            (-640., 256., 0.),
            (0., 320., 0.),
            (100., 320., 1.),
            (740., 384., 1.),
        ]
        .into_iter()
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.,
        })
        .collect(),
    )
    .unwrap()
}

#[test]
fn coupled_enthalpy_newton_keeps_latent_storage_and_pointwise_contact_adjoint() {
    context(|_, cx| {
        let mesh = mesh(false);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = source_card(false);
        let set = contacts(cx, &mesh, &boundary, &card, false);
        let curve = chart();
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &curve,
                    reference_density_kg_m3: 1.,
                },
                ReferenceEnthalpyMaterial {
                    curve: &curve,
                    reference_density_kg_m3: 2.,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic_declared(0.7).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let step = EnthalpyStepConfig {
            newton: NewtonKrylovConfig {
                relative_tolerance: 0.,
                absolute_tolerance: 1e-11,
                linear_restart: 8,
                max_linear_cycles: 16,
                forcing_minimum: 1e-12,
                forcing_maximum: 1e-3,
                ..Default::default()
            },
            max_newton_iterations: 64,
            energy_tolerance_j: 1e-8,
        };
        let config = ContactEnthalpyNewtonConfig {
            step,
            max_residual_evaluations: 1024,
            max_face_pairs: 1,
            residual_tolerance_j: 1e-9,
        };
        let old = [500., 300., 400., 450., -150., -50., 50., -200.];
        let before = old;
        let dt = 0.04;
        let accepted = set
            .advance_enthalpy_newton(cx, &engine, problem, &old, dt, config)
            .unwrap();
        assert_eq!(old, before);
        assert_eq!(accepted.physical_contacts.quadrature_states().len(), 3);
        let spatial = transport(&accepted.conduction.temperature);
        let h = &accepted.conduction.specific_enthalpy_j_kg;
        let mut energy = 0.;
        for i in 0..8 {
            let mass = if i < 4 { 1. / 24. } else { 2. / 24. };
            near(mass * (h[i] - old[i]) + dt * spatial[i], 0., 2e-9);
            energy += mass * (h[i] - old[i]);
        }
        near(energy, 0., 2e-9);
        assert!(
            h[6] > 50. && h[6] < 100.,
            "contact heat must enter the latent plateau"
        );
        assert_eq!(accepted.conduction.temperature[6], 320.);
        let linear = set
            .linearize_accepted_enthalpy(
                cx,
                &engine,
                problem,
                &old,
                dt,
                step,
                config.residual_tolerance_j,
                accepted,
                1,
            )
            .unwrap();
        assert_eq!(
            linear
                .temperature_pullback(cx, &[0., 0., 0., 0., 0., 0., 1., 0.])
                .unwrap(),
            vec![0.; 8]
        );
        let gradient = linear
            .pullback_with_resistance_controls(
                cx,
                &[0., 0., 0., 0., 1., 0., 0., 0.],
                LinearConfig {
                    tolerance: 1e-10,
                    max_iterations: 512,
                    restart: 8,
                },
                1,
            )
            .unwrap();
        assert_eq!(gradient.resistance_multipliers.len(), 1);
        assert_eq!(gradient.resistance_multipliers[0].face_pairs, 1);
        assert!(gradient.resistance_multipliers[0].derivative.is_finite());
    });
}
