use super::*;
use super::super::RobinLinearization;
use crate::{ConductionMesh, ConductivityModel, ConductivityTable, ElementMaterials,
    InitialGuess, MaterialId, MaterialTable, ScalarField, SolveConfig, ThermalBc,
    ThermalBoundaryBuilder};
use crate::fixtures::{box_grid, on_box_face};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_gate<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 47, kernel_id: 815, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T { with_gate(&CancelGate::new_clock_free(), f) }
fn mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([4, 2, 2], [0.1, 0.04, 0.03]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn config() -> SolveConfig {
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(320.0);
    config.linear.tolerance = 1e-11;
    config.stop.residual_rtol = 1e-12;
    config.stop.step_atol = 0.0;
    config
}
fn curve() -> ConductivityModel {
    ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(250.0, 1.0), (450.0, 21.0)]).unwrap())
}
fn source_delta(vertex: usize) -> f64 { (vertex % 7) as f64 * 20.0 - 60.0 }
fn solved(cx: &Cx<'_>, eps: f64, assigned: bool) -> RobinLinearization {
    let mesh = mesh();
    let model = curve();
    // An invalid-at-this-temperature fallback must not override the assigned law.
    let fallback = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(0.0, 1.0), (1.0, 2.0)]).unwrap());
    let materials = ElementMaterials::new(MaterialTable::new([(MaterialId(7), model.clone())]).unwrap(),
        vec![MaterialId(7); mesh.element_count()]).unwrap();
    let source = ScalarField::Nodal((0..mesh.vertex_count()).map(|v| 2000.0 + eps * source_delta(v)).collect());
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("hot", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(350.0).unwrap()).unwrap()
        .region("cooled", |f| on_box_face(f.centroid[0], 0.1), ThermalBc::robin(80.0 * (0.2 * eps).exp(), 290.0 + 0.6 * eps).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    RobinLinearization::new(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: if assigned { &fallback } else { &model },
        element_materials: assigned.then_some(&materials), source: &source }, config(), &["cooled"]).unwrap()
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*b.abs().max(1.0), "{a:e} != {b:e}");
}

#[test]
fn nonlinear_material_tangent_matches_perturbed_full_fem_and_pinned_nodes() {
    with_cx(|cx| {
        let base = solved(cx, 0.0, false);
        assert!(base.uses_nonlinear_jacobian());
        assert!(base.has_smooth_material_tangent());
        let mut direction = base.zero_direction();
        direction.references_k[0] = 0.6;
        direction.log_htc[0] = 0.2;
        // Differentiate the SAME consistent P1 source load as the perturbed
        // primal, rather than mislabeling a density vector as nodal watts.
        let mesh = mesh();
        for (e, tet) in mesh.complex().tets.iter().enumerate() {
            for (a, &vertex) in tet.iter().enumerate() {
                for (b, &other) in tet.iter().enumerate() {
                    direction.nodal_load_w[vertex as usize] += mesh.element_volume(e) / 20.0
                        * if a == b { 2.0 } else { 1.0 } * source_delta(other as usize);
                }
            }
        }
        let tangent = base.apply(cx, &direction).unwrap();
        let epsilon = 2e-4;
        let plus = solved(cx, epsilon, false);
        let minus = solved(cx, -epsilon, false);
        for (v, &actual) in tangent.temperature_k.iter().enumerate() {
            close(actual, (plus.primal().temperature[v] - minus.primal().temperature[v]) / (2.0*epsilon), 3e-5);
        }
        let heat = |s: &RobinLinearization| s.primal().report.robin_fluxes[0].heat_rate_w;
        close(tangent.heat_rates_w[0], (heat(&plus) - heat(&minus)) / (2.0*epsilon), 3e-5);
        for &v in base.dofs.fixed() { assert_eq!(tangent.temperature_k[v], 0.0); }
    });
}

#[test]
fn nonlinear_adjoint_uses_transpose_and_includes_load_and_direct_heat_terms() {
    with_cx(|cx| {
        let base = solved(cx, 0.0, false);
        let n = base.dofs.n();
        let vector: Vec<_> = (0..n).map(|i| (i % 5) as f64 - 2.0).collect();
        let mut forward = vec![0.0; n];
        let mut transpose = vec![0.0; n];
        base.matrix.spmv(&vector, &mut forward);
        base.nonlinear.as_ref().unwrap().transpose.spmv(&vector, &mut transpose);
        assert!(forward.iter().zip(&transpose).any(|(a,b)| (a-b).abs() > 1e-5), "fixture must distinguish J from J^T");
        let mut direction = base.zero_direction();
        direction.references_k[0] = -0.3;
        direction.log_htc[0] = 0.4;
        for (i, v) in direction.nodal_load_w.iter_mut().enumerate() { *v = (i % 7) as f64 * 0.003 - 0.009; }
        let weights: Vec<_> = (0..direction.nodal_load_w.len()).map(|i| (i % 9) as f64 * 0.01 - 0.04).collect();
        let tangent = base.apply(cx, &direction).unwrap();
        let gradient = base.pullback(cx, &weights, &[0.7], &[-0.2]).unwrap();
        let left = dot(&weights, &tangent.temperature_k) + 0.7*tangent.mean_wall_temperatures_k[0] - 0.2*tangent.heat_rates_w[0];
        let right = dot(&gradient.references, &direction.references_k) + dot(&gradient.log_htc, &direction.log_htc)
            + dot(&gradient.nodal_load, &direction.nodal_load_w);
        close(left, right, 2e-8);
        assert!(gradient.relative_residual < base.linear.tolerance);
        for &v in base.dofs.fixed() { assert_eq!(gradient.nodal_load[v], 0.0); }
    });
}

#[test]
fn assigned_nonlinear_law_overrides_the_unused_fallback_in_both_solves() {
    with_cx(|cx| {
        let uniform = solved(cx, 0.0, false);
        let assigned = solved(cx, 0.0, true);
        assert!(assigned.uses_nonlinear_jacobian());
        for (a,b) in uniform.primal().temperature.iter().zip(&assigned.primal().temperature) { close(*a,*b,1e-12); }
        let mut direction = uniform.zero_direction(); direction.references_k[0] = 1.0;
        let a = uniform.apply(cx, &direction).unwrap();
        let b = assigned.apply(cx, &direction).unwrap();
        for (a,b) in a.temperature_k.iter().zip(&b.temperature_k) { close(*a,*b,1e-10); }
    });
}

#[test]
fn a_material_kink_retains_the_primal_but_refuses_a_unique_derivative() {
    with_cx(|cx| {
        let mesh = mesh();
        let source = ScalarField::Uniform(0.0);
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![
            (280.0, 1.0), (300.0, 3.0), (330.0, 9.0)]).unwrap());
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .remainder("cooled", ThermalBc::robin(80.0, 300.0).unwrap()).unwrap().finish().unwrap();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source };
        let mut config = config(); config.initial = InitialGuess::Uniform(300.0);
        let base = RobinLinearization::new(cx, problem, config, &["cooled"]).unwrap();
        assert!(!base.has_smooth_material_tangent());
        assert!(base.primal().temperature.iter().all(|t| (*t-300.0).abs() < 1e-9));
        assert!(matches!(base.apply(cx, &base.zero_direction()), Err(ConductionError::Config { .. })));
        assert!(base.pullback(cx, &vec![0.0; mesh.vertex_count()], &[1.0], &[0.0]).is_err());
    });
}

#[test]
fn nonlinear_derivative_budget_zero_input_and_cancellation_are_explicit() {
    let mut base = with_cx(|cx| solved(cx, 0.0, false));
    with_cx(|cx| {
        let zero = base.apply(cx, &base.zero_direction()).unwrap();
        assert_eq!(zero.iterations, 0);
        assert!(zero.temperature_k.iter().all(|&t| t == 0.0));
        base.linear.max_iterations = 1;
        base.linear.restart = 1;
        let mut direction = base.zero_direction(); direction.references_k[0] = 1.0;
        assert!(matches!(base.apply(cx, &direction), Err(ConductionError::LinearSolveFailed { krylov_iterations: 1, .. })));
        base.linear.restart = 0;
        assert!(matches!(base.apply(cx, &direction), Err(ConductionError::Config { .. })));
    });
    let gate = CancelGate::new_clock_free(); gate.request();
    with_gate(&gate, |cx| assert!(matches!(base.apply(cx, &base.zero_direction()), Err(ConductionError::Cancelled { .. }))));
}
