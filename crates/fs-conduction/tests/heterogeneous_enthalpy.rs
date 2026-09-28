//! Heterogeneous reference storage and phase ownership through actual solves.

mod support;

use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, HeterogeneousEnthalpyError, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig,
};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable, ScalarField,
    ThermalBoundaryBuilder,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

fn chart(offset_k: f64, upper_h: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x45; 32]),
        [(0.0, 300.0, 0.0), (100.0, 350.0, 0.0), (300.0, 350.0, 1.0),
         (upper_h, 450.0, 1.0)]
            .into_iter()
            .map(|(h, temperature, liquid)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: temperature + offset_k,
                liquid_mass_fraction: liquid,
                // Deliberately different from the second material's reference
                // density: phase-state density must never redefine storage.
                bulk_density_kg_m3: 2.0,
            })
            .collect(),
    )
    .unwrap()
}

fn mesh() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        vec![
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0], [3.0, 0.0, 0.0],
            [2.0, 1.0, 0.0], [2.0, 0.0, 1.0],
        ],
    )
    .unwrap()
}

fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 8,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-10,
    }
}

fn budget() -> EnthalpyBudget {
    EnthalpyBudget { max_vertices: 8, max_elements: 2 }
}

#[test]
fn each_material_uses_its_own_chart_and_frozen_density() {
    let mesh = mesh();
    let solid = chart(0.0, 700.0);
    let pcm = chart(100.0, 700.0);
    let materials = [
        ReferenceEnthalpyMaterial { curve: &solid, reference_density_kg_m3: 2.0 },
        ReferenceEnthalpyMaterial { curve: &pcm, reference_density_kg_m3: 6.0 },
    ];
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let conductivity = ConductivityModel::isotropic(ConductivityTable::declared(2.0).unwrap());
    let source = ScalarField::nodal("uniform deposition", 8, vec![12.0; 8]).unwrap();
    let old = [50.0, 50.0, 50.0, 50.0, 180.0, 180.0, 180.0, 180.0];
    with_cx(|cx| {
        let stepper = HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh, &materials, &[0, 1], budget(),
        ).unwrap();
        assert_eq!(stepper.vertex_material_ids(), &[0, 0, 0, 0, 1, 1, 1, 1]);
        assert_eq!(stepper.element_material_ids(), &[0, 1]);
        assert!(std::ptr::eq(stepper.phase_curve_at(0).unwrap(), &solid));
        assert!(std::ptr::eq(stepper.phase_curve_at(7).unwrap(), &pcm));
        assert!(stepper.phase_curve_at(8).is_none());
        for (i, &mass) in stepper.reference_nodal_masses_kg().iter().enumerate() {
            let expected = if i < 4 { 1.0 / 12.0 } else { 0.25 };
            assert!((mass - expected).abs() < 1e-15);
        }
        let result = stepper.advance(cx, ConductionProblem {
            mesh: &mesh, boundary: &boundary, material: &conductivity,
            element_materials: None, source: &source,
        }, None, &old, 0.5, config()).unwrap();
        for i in 0..8 {
            let (h, temperature, fraction) = if i < 4 {
                (53.0, 326.5, 0.0)
            } else {
                (181.0, 450.0, 0.405)
            };
            assert!((result.specific_enthalpy_j_kg[i] - h).abs() < 1e-9);
            assert!((result.temperature[i] - temperature).abs() < 1e-9);
            assert!((result.liquid_mass_fraction[i] - fraction).abs() < 1e-9);
        }
        assert!((result.stored_energy_change_j - 2.0).abs() < 1e-10);
        assert!(result.energy_residual_j.abs() < 1e-10);
        assert_eq!(old[4], 180.0);
    });
}

#[test]
fn one_material_reproduces_the_uniform_stepper() {
    let mesh = mesh();
    let curve = chart(0.0, 700.0);
    let materials = [ReferenceEnthalpyMaterial {
        curve: &curve, reference_density_kg_m3: 2.0,
    }];
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let conductivity = ConductivityModel::isotropic(ConductivityTable::declared(2.0).unwrap());
    let source = ScalarField::nodal("deposition", 8, vec![3.0; 8]).unwrap();
    let old = [40.0, 50.0, 60.0, 70.0, 65.0, 55.0, 45.0, 35.0];
    with_cx(|cx| {
        let uniform = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 2.0, budget()).unwrap();
        let heterogeneous = HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh, &materials, &[0, 0], budget(),
        ).unwrap();
        let problem = || ConductionProblem {
            mesh: &mesh, boundary: &boundary, material: &conductivity,
            element_materials: None, source: &source,
        };
        let a = uniform.advance(cx, problem(), None, &old, 0.03, config()).unwrap();
        let b = heterogeneous.advance(cx, problem(), None, &old, 0.03, config()).unwrap();
        assert_eq!(a.specific_enthalpy_j_kg, b.specific_enthalpy_j_kg);
        assert_eq!(a.temperature, b.temperature);
        assert_eq!(a.liquid_mass_fraction, b.liquid_mass_fraction);
        assert_eq!(uniform.reference_nodal_masses_kg(), heterogeneous.reference_nodal_masses_kg());
    });
}

#[test]
fn refuses_implicit_material_mixtures_at_a_shared_vertex() {
    let mesh = ConductionMesh::new(
        TetComplex::from_tets(5, vec![[0, 1, 2, 3], [0, 2, 1, 4]]),
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
             [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
    ).unwrap();
    let a = chart(0.0, 700.0);
    let b = chart(100.0, 700.0);
    let materials = [
        ReferenceEnthalpyMaterial { curve: &a, reference_density_kg_m3: 2.0 },
        ReferenceEnthalpyMaterial { curve: &b, reference_density_kg_m3: 6.0 },
    ];
    with_cx(|cx| {
        assert!(matches!(
            HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &materials, &[0, 1], budget()),
            Err(HeterogeneousEnthalpyError::SharedVertex {
                vertex: 0, first_material: 0, next_material: 1,
            })
        ));
    });
}

#[test]
fn assignment_density_and_budget_refusals_are_explicit() {
    let mesh = mesh();
    let curve = chart(0.0, 700.0);
    with_cx(|cx| {
        for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let materials = [ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: density }];
            assert!(matches!(HeterogeneousEnthalpyBackwardEuler::new(
                cx, &mesh, &materials, &[0, 0], budget(),
            ), Err(HeterogeneousEnthalpyError::Enthalpy(EnthalpyError::InvalidInput(_)))));
        }
        let materials = [ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: 2.0 }];
        assert!(matches!(HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh, &materials, &[0, 1], budget(),
        ), Err(HeterogeneousEnthalpyError::UnknownMaterial { element: 1, material: 1, material_count: 1 })));
        assert!(HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &materials, &[0], budget()).is_err());
        assert!(HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &[], &[0, 0], budget()).is_err());
        assert!(matches!(HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh, &materials, &[0, 0], EnthalpyBudget { max_vertices: 7, max_elements: 2 },
        ), Err(HeterogeneousEnthalpyError::Enthalpy(EnthalpyError::Budget { resource: "vertices", .. }))));
    });
    with_cancelled_cx(|cx| {
        let materials = [ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: 2.0 }];
        assert!(HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &materials, &[0, 0], budget()).is_err());
    });
}

#[test]
fn domain_refusal_names_the_vertex_of_the_assigned_chart() {
    let mesh = mesh();
    let a = chart(0.0, 700.0);
    let b = chart(100.0, 500.0);
    let materials = [
        ReferenceEnthalpyMaterial { curve: &a, reference_density_kg_m3: 2.0 },
        ReferenceEnthalpyMaterial { curve: &b, reference_density_kg_m3: 6.0 },
    ];
    let boundary = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let conductivity = ConductivityModel::isotropic(ConductivityTable::declared(2.0).unwrap());
    let source = ScalarField::nodal("zero", 8, vec![0.0; 8]).unwrap();
    with_cx(|cx| {
        let stepper = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &materials, &[0, 1], budget()).unwrap();
        let old = [50.0, 50.0, 50.0, 50.0, 600.0, 400.0, 400.0, 400.0];
        assert!(matches!(stepper.advance(cx, ConductionProblem {
            mesh: &mesh, boundary: &boundary, material: &conductivity,
            element_materials: None, source: &source,
        }, None, &old, 0.03, config()), Err(EnthalpyError::Phase { vertex: 4, .. })));
    });
}
