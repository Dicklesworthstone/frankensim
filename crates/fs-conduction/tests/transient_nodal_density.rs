//! Graded reference density in the ordinary temperature-storage owners.

mod support;

use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{
    HeatSourceHistory, ThermalHistoryCapacity, ThermalHistoryMethod, ThermalHistoryProblem,
    ThermalHistoryRun,
};
use fs_conduction::fixtures::box_grid;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::transient::{
    HeatCapacity, TransientConfig, VolumetricHeatCapacity, assemble_capacitance,
    assemble_nodal_capacitance, march_with_capacity,
};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig,
    ScalarField, ThermalBoundaryBuilder,
};
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use support::{with_cancelled_cx, with_cx};

fn mesh(n: usize) -> ConductionMesh {
    let (complex, positions) = box_grid([n; 3], [1.0; 3]);
    ConductionMesh::new(
        TetComplex::from_tets(positions.len(), complex.tets),
        positions,
    )
    .unwrap()
}

fn capacity(mesh: &ConductionMesh) -> Vec<VolumetricHeatCapacity> {
    // A graded solid with constant specific heat 500 J/(kg K). All densities
    // and volumes refer to this same reference geometry throughout the run.
    mesh.positions()
        .iter()
        .map(|&[x, y, z]| {
            VolumetricHeatCapacity::declared(500.0 * (800.0 + 7000.0 * x + 300.0 * y + 600.0 * z))
                .unwrap()
        })
        .collect()
}

fn linear() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-12,
        max_iterations: 4096,
        restart: 32,
    }
}

fn controls() -> StepConfig {
    StepConfig {
        linear: linear(),
        energy_tolerance_j: 2e-6,
    }
}

// Independent degree-two tetrahedron quadrature integrates c_h*T_h exactly;
// this deliberately does not use the storage matrix or its row-moment formula.
fn heat(mesh: &ConductionMesh, capacities: &[VolumetricHeatCapacity], temperature: &[f64]) -> f64 {
    let mut integral = 0.0;
    for (element, tet) in mesh.complex().tets.iter().enumerate() {
        for point in 0..4 {
            let shape: [f64; 4] = std::array::from_fn(|a| {
                if a == point {
                    0.585_410_196_624_968_5
                } else {
                    0.138_196_601_125_010_5
                }
            });
            let c: f64 = tet
                .iter()
                .zip(shape)
                .map(|(&vertex, n)| n * capacities[vertex as usize].value_j_per_m3_k())
                .sum();
            let t: f64 = tet
                .iter()
                .zip(shape)
                .map(|(&vertex, n)| n * temperature[vertex as usize])
                .sum();
            integral += mesh.element_volume(element) * 0.25 * c * t;
        }
    }
    integral
}

#[test]
fn graded_heat_content_matches_the_continuum_and_preserves_density_temperature_correlation() {
    // Exact unit-box integral of two affine fields: mean(c)*mean(T) plus
    // their spatial covariance; replacing c by a cell mean drops covariance.
    let expected =
        500.0 * (800.0 + 0.5 * (7000.0 + 300.0 + 600.0)) * (300.0 + 0.5 * (100.0 - 30.0 + 10.0))
            + 500.0 * (7000.0 * 100.0 - 300.0 * 30.0 + 600.0 * 10.0) / 12.0;
    for n in [1, 2, 4] {
        let mesh = mesh(n);
        let capacities = capacity(&mesh);
        let temperature: Vec<f64> = mesh
            .positions()
            .iter()
            .map(|&[x, y, z]| 300.0 + 100.0 * x - 30.0 * y + 10.0 * z)
            .collect();
        let matrix = with_cx(|cx| assemble_nodal_capacitance(cx, &mesh, &capacities).unwrap());
        let mut stored = vec![0.0; mesh.vertex_count()];
        matrix.spmv(&temperature, &mut stored);
        let actual: f64 = stored.iter().sum();
        assert!(
            (actual - expected).abs() < 2e-14 * expected,
            "n={n}: {actual} != {expected}"
        );
        assert!((heat(&mesh, &capacities, &temperature) - expected).abs() < 2e-14 * expected);

        let averaged: f64 = mesh
            .complex()
            .tets
            .iter()
            .enumerate()
            .map(|(element, tet)| {
                let c: f64 = tet
                    .iter()
                    .map(|&v| capacities[v as usize].value_j_per_m3_k())
                    .sum();
                let t: f64 = tet.iter().map(|&v| temperature[v as usize]).sum();
                mesh.element_volume(element) * c * t / 16.0
            })
            .sum();
        assert!(
            (averaged - expected).abs() > 1e-4 * expected,
            "the fixture must expose the lost density-temperature correlation"
        );
    }
}

#[test]
fn uniform_specific_heating_has_the_exact_ramp_for_theta_backward_euler_and_history() {
    let mesh = mesh(2);
    let capacities = capacity(&mesh);
    let rate = 0.2;
    let source = ScalarField::nodal(
        "specific heating",
        mesh.vertex_count(),
        capacities
            .iter()
            .map(|c| rate * c.value_j_per_m3_k())
            .collect(),
    )
    .unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let initial = vec![300.0; mesh.vertex_count()];
    with_cx(|cx| {
        for theta in [0.5, 0.75, 1.0] {
            let config = TransientConfig::new(theta, 0.5, linear()).unwrap();
            let solved = march_with_capacity(
                cx,
                problem,
                HeatCapacity::Nodal(&capacities),
                &config,
                &initial,
                4,
            )
            .unwrap();
            assert!(solved.temperature.iter().all(|&t| (t - 300.4).abs() < 1e-9));
        }
        let stepper = BackwardEuler::nodal(cx, &mesh, &capacities).unwrap();
        let solved = stepper
            .advance(cx, problem, None, &initial, 2.0, controls())
            .unwrap();
        assert!(solved.temperature.iter().all(|&t| (t - 300.4).abs() < 1e-9));
        let expected_energy = heat(&mesh, &capacities, &vec![0.4; mesh.vertex_count()]);
        assert!((solved.stored_energy_change_j - expected_energy).abs() < 2e-6);

        let history = ThermalBoundaryHistory::new(
            cx,
            ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder(),
            &Environment::earth_lab(),
            Vec::new(),
            0.0,
        )
        .unwrap();
        let mut run = ThermalHistoryRun::new_with_capacity(
            cx,
            ThermalHistoryProblem {
                boundary: &history,
                material: &material,
                element_materials: None,
                interfaces: None,
                source: HeatSourceHistory::Fixed(source.clone()),
            },
            ThermalHistoryCapacity::Nodal(&capacities),
            &initial,
            2,
            ThermalHistoryMethod::Linear,
        )
        .unwrap();
        run.advance_to(cx, 0.75, controls()).unwrap();
        let endpoint = run.advance_to(cx, 2.0, controls()).unwrap();
        assert!(
            endpoint
                .temperature
                .iter()
                .all(|&t| (t - 300.4).abs() < 1e-9)
        );
    });
}

#[test]
fn insulated_diffusion_preserves_reference_density_weighted_heat() {
    let mesh = mesh(2);
    let capacities = capacity(&mesh);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(5000.0).unwrap();
    let source = ScalarField::uniform("no generation", 0.0).unwrap();
    let initial: Vec<f64> = mesh
        .positions()
        .iter()
        .map(|&[x, _, _]| 300.0 + 100.0 * x)
        .collect();
    let initial_heat = heat(&mesh, &capacities, &initial);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    with_cx(|cx| {
        let stepper = BackwardEuler::nodal(cx, &mesh, &capacities).unwrap();
        let solved = stepper
            .advance(cx, problem, None, &initial, 10.0, controls())
            .unwrap();
        assert!((heat(&mesh, &capacities, &solved.temperature) - initial_heat).abs() < 2e-6);
        assert!(solved.energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(
            initial
                .iter()
                .zip(&solved.temperature)
                .any(|(old, new)| (old - new).abs() > 1.0)
        );
        for theta in [0.5, 1.0] {
            let config = TransientConfig::new(theta, 2.0, linear()).unwrap();
            let solved = march_with_capacity(
                cx,
                problem,
                HeatCapacity::Nodal(&capacities),
                &config,
                &initial,
                5,
            )
            .unwrap();
            assert!((heat(&mesh, &capacities, &solved.temperature) - initial_heat).abs() < 5e-6);
        }
    });
}

#[test]
fn uniform_nodal_storage_preserves_the_original_matrix_and_rejects_incomplete_fields() {
    let mesh = mesh(2);
    let value = VolumetricHeatCapacity::declared(2e6).unwrap();
    let capacities = vec![value; mesh.vertex_count()];
    with_cx(|cx| {
        let original = assemble_capacitance(cx, &mesh, value).unwrap();
        let nodal = assemble_nodal_capacitance(cx, &mesh, &capacities).unwrap();
        for row in 0..mesh.vertex_count() {
            assert_eq!(original.row(row), nodal.row(row));
        }
        assert!(matches!(
            BackwardEuler::nodal(cx, &mesh, &capacities[1..]),
            Err(ConductionError::FieldLength {
                field: "nodal heat capacity",
                ..
            })
        ));
    });
    with_cancelled_cx(|cx| {
        assert!(matches!(
            assemble_nodal_capacitance(cx, &mesh, &capacities),
            Err(ConductionError::Cancelled { .. })
        ));
        assert!(matches!(
            BackwardEuler::nodal(cx, &mesh, &capacities),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}
