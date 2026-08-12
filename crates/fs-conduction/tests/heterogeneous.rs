//! Heterogeneous per-element conductivity (bead s93ej.2).
//!
//! One-material assignment must reproduce the uniform path. A two-layer
//! slab with k=1 then k=2 has the series-resistance interface
//! temperature 2/3 of the way from the cold face when the layers are
//! equal thickness.

mod support;

use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::{ConductivityModel, ElementMaterials, MaterialId, MaterialTable};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::{
    ConductionError, ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, SolveConfig,
    StopRule, solve, solve_with_element_materials,
};
use support::with_cx;

fn config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 4,
        },
        stop: StopRule {
            residual_rtol: 1e-12,
            residual_atol: 1e-14,
            step_atol: 1e-12,
            max_iterations: 8,
        },
        linear: LinearConfig {
            tolerance: 1e-14,
            max_iterations: 200,
            restart: 40,
        },
        initial: InitialGuess::DirichletMean,
    }
}

#[test]
fn empty_table_and_unknown_id_refuse() {
    assert!(matches!(
        MaterialTable::new(Vec::new()),
        Err(ConductionError::MaterialAssignment { .. })
    ));
    let table = MaterialTable::new([(
        MaterialId(1),
        ConductivityModel::isotropic_declared(1.0).expect("k"),
    )])
    .expect("table");
    let err = ElementMaterials::new(table, vec![MaterialId(2)]).expect_err("unknown");
    assert!(matches!(err, ConductionError::MaterialAssignment { .. }));
}

#[test]
fn one_material_assignment_matches_the_uniform_path() {
    let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(20.0).expect("k");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(340.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
    };
    let uniform = with_cx(|cx| solve(cx, problem, config()).expect("uniform"));
    let table = MaterialTable::new([(MaterialId(1), material.clone())]).expect("table");
    let assigned =
        ElementMaterials::new(table, vec![MaterialId(1); mesh.element_count()]).expect("assign");
    let hetero = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("hetero")
    });
    let worst = uniform
        .temperature
        .iter()
        .zip(&hetero.temperature)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(
        worst < 1e-12,
        "uniform vs one-material assignment differed by {worst:e}"
    );
}

#[test]
fn two_layer_slab_hits_the_series_interface_temperature() {
    let length = 2.0;
    let (complex, positions) = box_grid([8, 2, 2], [length, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left = ConductivityModel::isotropic_declared(1.0).expect("k1");
    let right = ConductivityModel::isotropic_declared(2.0).expect("k2");
    let fallback = left.clone();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(1.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cold",
            |f| on_box_face(f.centroid[0], length),
            ThermalBc::dirichlet(0.0).expect("bc"),
        )
        .expect("cold")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let table = MaterialTable::new([(MaterialId(1), left), (MaterialId(2), right)]).expect("table");
    let ids: Vec<MaterialId> = (0..mesh.element_count())
        .map(|e| {
            let tet = mesh.complex().tets[e];
            let cx = tet
                .iter()
                .map(|&v| mesh.positions()[v as usize][0])
                .sum::<f64>()
                / 4.0;
            if cx < 1.0 {
                MaterialId(1)
            } else {
                MaterialId(2)
            }
        })
        .collect();
    let assigned = ElementMaterials::new(table, ids).expect("assign");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &fallback,
        source: &source,
    };
    let solution = with_cx(|cx| {
        solve_with_element_materials(cx, problem, &assigned, config()).expect("solve")
    });
    // Equal thicknesses, k1=1, k2=2, T(0)=1, T(2)=0. Heat flows +x:
    // q = 1 / (1/1 + 1/2) = 2/3, so T(1) = 1 − q/k1 = 1/3.
    let mut worst = 0.0f64;
    let mut counted = 0usize;
    for (v, &p) in mesh.positions().iter().enumerate() {
        if (p[0] - 1.0).abs() < 1e-9 {
            worst = worst.max((solution.temperature[v] - 1.0 / 3.0).abs());
            counted += 1;
        }
    }
    assert!(counted > 0, "no interface nodes at x=1");
    assert!(
        worst < 2e-3,
        "interface temperature off by {worst:e} from 1/3"
    );
}
