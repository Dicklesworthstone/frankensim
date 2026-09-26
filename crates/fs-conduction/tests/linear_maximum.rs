//! Regional maxima through the public production assembly/analyzer seam.
//! The dense oracle deliberately does not call the production Krylov solver.
mod support;

use fs_conduction::adjoint::{
    LinearGoalAnalysisConfig, LinearGoalAnalyzer, analyze_linear_maximum,
};
use fs_conduction::assemble::{DofMap, assemble_operator, reduce};
use fs_conduction::fixtures::unit_cube;
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig, ScalarField,
    ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_solver::goal::GoalResidualLimits;
use support::{with_cancelled_cx, with_cx};

struct Fixture {
    mesh: ConductionMesh,
    boundary: ThermalBoundary,
    material: ConductivityModel,
    source: ScalarField,
}
impl Fixture {
    fn new(n: usize) -> Self {
        let (complex, positions) = unit_cube(n);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("reservoir", |_| true, ThermalBc::dirichlet(300.0).unwrap())
            .unwrap().finish().unwrap();
        let source = ScalarField::Nodal(mesh.positions().iter()
            .map(|p| 1000.0 * (1.0 + p[0] + 0.2 * p[1] + 0.04 * p[2])).collect());
        Self { mesh, boundary, material: ConductivityModel::isotropic_declared(2.0).unwrap(), source }
    }
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary,
            material: &self.material, source: &self.source, element_materials: None }
    }
    fn values(&self) -> Vec<f64> { vec![300.0; self.mesh.vertex_count()] }
    fn vertices(&self) -> Vec<usize> { (0..self.mesh.vertex_count()).collect() }
}
fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-12, max_iterations: 256, restart: 32 }
}
fn policy(stability: usize) -> LinearGoalAnalysisConfig {
    LinearGoalAnalysisConfig { residual_limits: GoalResidualLimits {
        max_rows: 4096, max_nonzeros: 65536,
    }, max_stability_iterations: stability }
}

fn dense_solve(a: &fs_sparse::Csr, rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut rows: Vec<Vec<f64>> = (0..n).map(|i| {
        let mut row: Vec<_> = (0..n).map(|j| a.get(i, j)).collect();
        row.push(rhs[i]); row
    }).collect();
    for k in 0..n {
        let pivot = (k..n).max_by(|&i, &j| rows[i][k].abs().total_cmp(&rows[j][k].abs())).unwrap();
        rows.swap(k, pivot);
        assert!(rows[k][k].abs() > 0.0);
        for i in k+1..n {
            let ratio = rows[i][k] / rows[k][k];
            for j in k+1..=n { rows[i][j] -= ratio * rows[k][j]; }
            rows[i][k] = 0.0;
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let tail: f64 = (i+1..n).map(|j| rows[i][j] * x[j]).sum();
        x[i] = (rows[i][n] - tail) / rows[i][i];
    }
    x
}

#[test]
fn moving_maximum_is_enclosed_when_the_old_hottest_node_has_no_error() {
    let f = Fixture::new(3);
    with_cx(|cx| {
        let system = assemble_operator(cx, &f.mesh, &f.boundary, &f.material, &f.source, &f.values()).unwrap();
        let dofs = DofMap::new(&f.boundary, f.mesh.vertex_count()).unwrap();
        let (a, b) = reduce(&system, &dofs);
        let exact = dofs.scatter(&dense_solve(&a, &b));
        let actual_hot = (0..exact.len()).max_by(|&i, &j| exact[i].total_cmp(&exact[j])).unwrap();
        let mut approximate = exact.clone();
        approximate[actual_hot] = 300.0;
        let old_hot = (0..exact.len()).max_by(|&i, &j| approximate[i].total_cmp(&approximate[j])).unwrap();
        assert_ne!(old_hot, actual_hot);
        let gap = exact[actual_hot] - approximate[old_hot];
        assert!(gap > 1e-4, "fixture must detect an omitted maximum relocation");
        assert_eq!(approximate[old_hot], exact[old_hot]);
        let result = analyze_linear_maximum(cx, f.problem(), None, linear(),
            &approximate, &f.vertices(), policy(256)).unwrap();
        let interval = result.interval_k().expect("scaled inverse established");
        assert!(interval[0] <= exact[actual_hot] && exact[actual_hot] <= interval[1]);
        assert!(result.algebraic_half_width_k().unwrap() >= gap);
        assert_eq!(result.linear_analysis().dual_iterations, 0, "no hottest-node dual needed");
    });
}

#[test]
fn cached_assessment_replays_and_set_permutation_does_not_change_bits() {
    let f = Fixture::new(3);
    with_cx(|cx| {
        let values = f.values();
        let analyzer = LinearGoalAnalyzer::new_for_maximum(
            cx, f.problem(), None, linear(), &values, policy(256)).unwrap();
        let mut vertices = f.vertices();
        let first = analyzer.analyze_maximum(cx, &values, &vertices).unwrap();
        vertices.reverse();
        let reversed = analyzer.analyze_maximum(cx, &values, &vertices).unwrap();
        assert_eq!(first, reversed);
        assert_eq!(first, analyzer.analyze_maximum(cx, &values, &vertices).unwrap());
        assert!(first.algebraic_half_width_k().unwrap() > 0.0);
        assert!(!first.meets_absolute_tolerance(1e-12));
        assert!(first.meets_absolute_tolerance(first.algebraic_half_width_k().unwrap()));
    });
}

#[test]
fn missing_inverse_cannot_pass_but_prescribed_only_selection_is_exact() {
    // A genuinely non-dominant anisotropic operator, not a nearly zero
    // isotropic row sum whose rounding could legitimately prove a tiny margin.
    // Zero witness iterations must not manufacture an inverse bound.
    let mut f = Fixture::new(4);
    f.material = ConductivityModel::constant_tensor([
        [3.0, 0.5, 0.25], [0.5, 2.0, 0.75], [0.25, 0.75, 1.5],
    ]).unwrap();
    with_cx(|cx| {
        let values = f.values();
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(),
            &values, &vec![0.0; values.len()], policy(0)).unwrap();
        let unbounded = analyzer.analyze_maximum(cx, &values, &f.vertices()).unwrap();
        assert!(unbounded.interval_k().is_none());
        assert!(!unbounded.meets_absolute_tolerance(f64::MAX));
        let fixed = analyzer.dofs().fixed().to_vec();
        let exact = analyzer.analyze_maximum(cx, &values, &fixed).unwrap();
        assert_eq!(exact.interval_k(), Some([300.0, 300.0]));
        assert_eq!(exact.algebraic_half_width_k(), Some(0.0));
        assert_eq!(exact.free_vertices(), 0);
    });
}

#[test]
fn invalid_selections_fields_and_tolerances_refuse() {
    let f = Fixture::new(2);
    with_cx(|cx| {
        let values = f.values();
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(),
            &values, &vec![0.0; values.len()], policy(16)).unwrap();
        for vertices in [vec![], vec![0, 0], vec![usize::MAX], vec![values.len()]] {
            assert!(analyzer.analyze_maximum(cx, &values, &vertices).is_err());
        }
        assert!(analyzer.analyze_maximum(cx, &values[1..], &[0]).is_err());
        let mut invalid = values.clone(); invalid[0] = f64::NAN;
        assert!(analyzer.analyze_maximum(cx, &invalid, &[0]).is_err());
        invalid[0] = 301.0;
        assert!(analyzer.analyze_maximum(cx, &invalid, &[0]).is_err());
        let result = analyzer.analyze_maximum(cx, &values, &f.vertices()).unwrap();
        for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(!result.meets_absolute_tolerance(tolerance));
        }
    });
}

#[test]
fn nonlinear_conductivity_is_not_relabelled_as_a_linear_maximum_bound() {
    let mut f = Fixture::new(2);
    f.material = ConductivityModel::isotropic(
        fs_conduction::material::ConductivityTable::declared_curve(vec![(290.0, 2.0), (400.0, 4.0)]).unwrap());
    with_cx(|cx| {
        assert!(analyze_linear_maximum(cx, f.problem(), None, linear(),
            &f.values(), &f.vertices(), policy(32)).is_err());
    });
}

#[test]
fn cancellation_does_not_publish_a_maximum_enclosure() {
    let f = Fixture::new(2);
    with_cancelled_cx(|cx| {
        assert!(matches!(analyze_linear_maximum(cx, f.problem(), None, linear(),
            &f.values(), &f.vertices(), policy(32)), Err(fs_conduction::ConductionError::Cancelled { .. })));
    });
}

#[test]
fn maximum_preparation_improves_loose_bounds_without_repeating_witness_work() {
    let f = Fixture::new(4);
    with_cx(|cx| {
        let values = f.values();
        let unscaled = LinearGoalAnalyzer::new_for_maximum(
            cx, f.problem(), None, linear(), &values, policy(0)).unwrap()
            .analyze_maximum(cx, &values, &f.vertices()).unwrap();
        let analyzer = LinearGoalAnalyzer::new_for_maximum(
            cx, f.problem(), None, linear(), &values, policy(256)).unwrap();
        let assessed = analyzer.analyze_maximum(cx, &values, &f.vertices()).unwrap();
        let bound = assessed.algebraic_half_width_k().expect("useful scaled bound");
        assert!(bound < 1000.0, "a near-zero row margin must not leave a huge bound: {bound}");
        if let Some(old) = unscaled.algebraic_half_width_k() { assert!(bound <= old); }
        assert!(assessed.linear_analysis().stability_iterations <= 256);
        assert_eq!(assessed.linear_analysis().dual_iterations, 0);
        assert_eq!(assessed, analyzer.analyze_maximum(cx, &values, &f.vertices()).unwrap());
    });
}
