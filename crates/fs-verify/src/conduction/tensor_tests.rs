use super::*;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionMesh, ConductivityModel, ElementMaterials, InitialGuess,
    MaterialId, MaterialTable, Nonlinearity};
use fs_conduction::fixtures::unit_cube;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 73, kernel_id: 11, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}
fn config() -> MeanSolveConfig {
    let mut primal = SolveConfig::default();
    primal.nonlinearity = Nonlinearity::FixedPoint { relaxation: 1.0, max_backtracks: 4 };
    primal.stop.residual_rtol = 1e-11;
    primal.stop.residual_atol = 1e-12;
    primal.linear.tolerance = 1e-13;
    let mut dual = primal.clone();
    dual.initial = InitialGuess::Uniform(0.0);
    MeanSolveConfig { primal, dual, flux: FluxBudget::default() }
}
fn mesh(n: usize) -> ConductionMesh {
    let (complex, positions) = unit_cube(n);
    ConductionMesh::new(complex, positions).unwrap()
}
fn ends(mesh: &ConductionMesh, right: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("left", |f| f.vertices.iter().all(|&v| mesh.positions()[v as usize][0] == 0.0), ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .region("right", |f| f.vertices.iter().all(|&v| mesh.positions()[v as usize][0] == 1.0), ThermalBc::dirichlet(right).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}
fn contains(interval: crate::interval::Iv, truth: f64) {
    assert!(interval.lo <= truth && truth <= interval.hi, "{interval:?} excludes {truth}");
}

#[test]
fn real_off_diagonal_neumann_solve_and_mean_follow_the_rotated_material_frame() {
    with_cx(|cx| {
        let mut reference: Option<MeanTemperatureSolution> = None;
        for axes in [[0,1,2], [2,0,1]] {
            let (complex, original) = unit_cube(2);
            let positions = original.iter().map(|p| axes.map(|i| p[i])).collect();
            let mesh = ConductionMesh::new(complex, positions).unwrap();
            let k: ConductivityTensor = [[2.0,1.0,0.0],[1.0,3.0,0.0],[0.0,0.0,4.0]];
            let k = axes.map(|i| axes.map(|j| k[i][j]));
            let material = ConductivityModel::constant_tensor(k).unwrap();
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .region("left", |f| f.vertices.iter().all(|&v| original[v as usize][0] == 0.0), ThermalBc::dirichlet(300.0).unwrap()).unwrap()
                .region("right", |f| f.vertices.iter().all(|&v| original[v as usize][0] == 1.0), ThermalBc::dirichlet(301.0).unwrap()).unwrap()
                .region("y-in", |f| f.vertices.iter().all(|&v| original[v as usize][1] == 0.0), ThermalBc::neumann(1.0).unwrap()).unwrap()
                .region("y-out", |f| f.vertices.iter().all(|&v| original[v as usize][1] == 1.0), ThermalBc::neumann(-1.0).unwrap()).unwrap()
                .adiabatic_remainder().finish().unwrap();
            let source = ScalarField::Uniform(0.0);
            let result = solve_with_mean_bound(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source }, config()).unwrap();
            contains(result.bound.enclosure, 300.5);
            assert!(result.bound.enclosure.hi-result.bound.enclosure.lo < 1e-6);
            for (p, &t) in original.iter().zip(&result.primal.temperature) {
                assert!((t-(300.0+p[0])).abs() < 1e-8);
            }
            if let Some(old) = &reference {
                assert!((old.bound.enclosure.hi-result.bound.enclosure.hi).abs() < 1e-7);
                for (&a,&b) in old.primal.temperature.iter().zip(&result.primal.temperature) {
                    assert!((a-b).abs() < 1e-8);
                }
            }
            reference = Some(result);
        }
    });
}

#[test]
fn heterogeneous_tensors_preserve_the_series_interface_not_the_fallback() {
    with_cx(|cx| {
        let mesh = mesh(2);
        let boundary = ends(&mesh, 301.0);
        let fallback = ConductivityModel::isotropic_declared(99.0).unwrap();
        let table = MaterialTable::new([
            (MaterialId(1), ConductivityModel::constant_tensor([[1.0,0.0,0.0],[0.0,8.0,1.0],[0.0,1.0,3.0]]).unwrap()),
            (MaterialId(2), ConductivityModel::constant_tensor([[4.0,0.0,0.0],[0.0,2.0,0.5],[0.0,0.5,1.0]]).unwrap()),
        ]).unwrap();
        let ids = mesh.complex().tets.iter().map(|tet| {
            let x = tet.iter().map(|&v| mesh.positions()[v as usize][0]).sum::<f64>()/4.0;
            if x < 0.5 { MaterialId(1) } else { MaterialId(2) }
        }).collect();
        let assigned = ElementMaterials::new(table, ids).unwrap();
        let source = ScalarField::Uniform(0.0);
        let result = solve_with_mean_bound(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &fallback, element_materials: Some(&assigned), source: &source }, config()).unwrap();
        contains(result.bound.enclosure, 300.65);
        assert!(result.bound.enclosure.hi-result.bound.enclosure.lo < 1e-6);
        for (p,&t) in mesh.positions().iter().zip(&result.primal.temperature) {
            if p[0] == 0.5 { assert!((t-300.8).abs() < 1e-8); }
        }
    });
}

#[test]
fn existing_field_path_reuses_the_primal_and_matches_the_solved_bound() {
    with_cx(|cx| {
        let mesh = mesh(2);
        let boundary = ends(&mesh, 300.0);
        let material = ConductivityModel::constant_tensor([[4.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,1.0]]).unwrap();
        let source = ScalarField::Uniform(8.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let full = solve_with_mean_bound(cx, problem, config()).unwrap();
        let candidate = full.primal.temperature.clone();
        let reused = bound_temperature_mean(cx, problem, &candidate, config().dual, FluxBudget::default()).unwrap();
        assert_eq!(reused.bound.enclosure, full.bound.enclosure);
        assert_eq!(reused.bound.candidate_mean, full.bound.candidate_mean);
        assert_eq!(reused.dual.temperature, full.dual.temperature);
        assert_eq!(candidate, full.primal.temperature);
        contains(reused.bound.enclosure, 300.0+1.0/6.0);
    });
}

#[test]
fn unfinished_primal_keeps_its_algebraic_error_in_the_mean_bound() {
    with_cx(|cx| {
        let mesh = mesh(4);
        let boundary = ends(&mesh, 300.0);
        let material = ConductivityModel::constant_tensor([[4.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,1.0]]).unwrap();
        let source = ScalarField::Uniform(8.0);
        let candidate = vec![300.0;mesh.vertex_count()];
        let result = bound_temperature_mean(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source }, &candidate, config().dual, FluxBudget::default()).unwrap();
        contains(result.bound.enclosure, 300.0+1.0/6.0);
        contains(result.bound.candidate_mean, 300.0);
        assert!(result.bound.integral.residual_correction.lo > 0.1);
        assert!(candidate.iter().all(|&t| t == 300.0));
    });
}

#[test]
fn directional_conductivity_changes_the_real_field_and_bound() {
    with_cx(|cx| {
        let mesh = mesh(4);
        let boundary = ends(&mesh, 300.0);
        let source = ScalarField::Uniform(8.0);
        let mut means = Vec::new();
        for k in [1.0,4.0] {
            let material = ConductivityModel::constant_tensor([[k,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,1.0]]).unwrap();
            let result = solve_with_mean_bound(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source }, config()).unwrap();
            contains(result.bound.enclosure, 300.0+8.0/(12.0*k));
            means.push(result.bound.candidate_mean.hi);
        }
        assert!(means[0]-means[1] > 0.4);
    });
}

#[test]
fn malformed_existing_fields_refuse_before_the_dual_solver() {
    with_cx(|cx| {
        let mesh = mesh(2);
        let boundary = ends(&mesh, 300.0);
        let material = ConductivityModel::isotropic_declared(1.0).unwrap();
        let source = ScalarField::Uniform(2.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let valid = vec![300.0;mesh.vertex_count()];
        let mut not_finite = valid.clone(); not_finite[0] = f64::NAN;
        let mut wrong_trace = valid.clone(); wrong_trace[0] = 299.0;
        for candidate in [vec![300.0;mesh.vertex_count()-1], not_finite, wrong_trace] {
            // Deliberately invalid solver settings would fail differently if
            // candidate admission were postponed until after the dual solve.
            let mut dual = config().dual;
            dual.linear.restart = 0;
            assert!(matches!(bound_temperature_mean(cx, problem, &candidate, dual, FluxBudget::default()),
                Err(ConductionBoundError::Verification(TetError::Invalid(_)))));
        }
    });
}
