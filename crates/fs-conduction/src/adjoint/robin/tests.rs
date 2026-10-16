use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use crate::{ConductionMesh, ConductivityModel, InitialGuess, Nonlinearity, ThermalBoundaryBuilder};
use crate::fixtures::{box_grid, on_box_face};

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 31, kernel_id: 814, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R { with_gate(&CancelGate::new(), f) }
fn config() -> SolveConfig {
    let mut c = SolveConfig::default();
    c.nonlinearity = Nonlinearity::FixedPoint { relaxation: 1.0, max_backtracks: 8 };
    c.initial = InitialGuess::Uniform(310.0);
    c.linear.tolerance = 1e-10;
    c.stop.residual_rtol = 1e-11;
    c.stop.step_atol = 0.0;
    c
}
fn slab(cx: &Cx<'_>, h: [f64; 2], refs: [f64; 2], fixed_left: bool) -> RobinLinearization {
    let (complex, positions) = box_grid([3, 2, 2], [0.1, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("left", |f| on_box_face(f.centroid[0], 0.0),
            if fixed_left { ThermalBc::dirichlet(refs[0]).unwrap() }
            else { ThermalBc::robin(h[0], refs[0]).unwrap() }).unwrap()
        .region("right", |f| on_box_face(f.centroid[0], 0.1), ThermalBc::robin(h[1], refs[1]).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let regions: &[&str] = if fixed_left { &["right"] } else { &["left", "right"] };
    RobinLinearization::new(cx, ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source }, config(), regions).unwrap()
}
fn close(a: f64, b: f64) { assert!((a-b).abs() < 2e-7 * b.abs().max(1.0), "{a:.14e} != {b:.14e}"); }
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }

#[test]
fn robin_reference_and_conductance_tangents_match_series_resistance() {
    with_cx(|cx| {
        let linear = slab(cx, [20.0, 50.0], [330.0, 290.0], false);
        let resistance = 1.0/20.0 + 0.1/10.0 + 1.0/50.0;
        let heat = 40.0/resistance;
        let mut d = linear.zero_direction();
        d.references_k = vec![0.7, -0.3];
        d.log_htc = vec![0.2, -0.1];
        let dq = (1.0 + heat*(0.2/20.0 - 0.1/50.0))/resistance;
        let result = linear.apply(cx, &d).unwrap();
        close(result.mean_wall_temperatures_k[0], 0.7 - dq/20.0 + heat*0.2/20.0);
        close(result.mean_wall_temperatures_k[1], -0.3 + dq/50.0 + heat*0.1/50.0);
        close(result.heat_rates_w[0], -dq);
        close(result.heat_rates_w[1], dq);
        assert!(result.relative_residual < linear.linear.tolerance);
    });
}

#[test]
fn full_field_tangent_matches_independent_perturbed_fem_solves() {
    with_cx(|cx| {
        let base = slab(cx, [20.0, 50.0], [330.0, 290.0], false);
        let mut d = base.zero_direction();
        d.references_k = vec![0.4, -0.2];
        d.log_htc = vec![-0.3, 0.6];
        let result = base.apply(cx, &d).unwrap();
        let eps = 2e-4_f64;
        let plus = slab(cx, [20.0*(-0.3*eps).exp(), 50.0*(0.6*eps).exp()], [330.0+0.4*eps, 290.0-0.2*eps], false);
        let minus = slab(cx, [20.0*(0.3*eps).exp(), 50.0*(-0.6*eps).exp()], [330.0-0.4*eps, 290.0+0.2*eps], false);
        for ((a,b), d) in plus.primal.temperature.iter().zip(&minus.primal.temperature).zip(&result.temperature_k) {
            assert!(((a-b)/(2.0*eps)-d).abs() < 2e-5, "finite difference {}, tangent {d}", (a-b)/(2.0*eps));
        }
    });
}

#[test]
fn transpose_identity_includes_nodal_loads_and_direct_heat_terms() {
    with_cx(|cx| {
        let base = slab(cx, [20.0, 50.0], [330.0, 290.0], false);
        let mut d = base.zero_direction();
        d.references_k = vec![0.3, -0.4];
        d.log_htc = vec![0.2, 0.1];
        for (i,v) in d.nodal_load_w.iter_mut().enumerate() { *v = (i as f64 % 5.0 - 2.0)*0.1; }
        let nodal: Vec<f64> = (0..d.nodal_load_w.len()).map(|i| (i as f64 % 7.0 - 3.0)*0.02).collect();
        let walls = [0.8, -0.6];
        let heats = [0.003, -0.007];
        let jvp = base.apply(cx, &d).unwrap();
        let vjp = base.pullback(cx, &nodal, &walls, &heats).unwrap();
        let left = dot(&nodal, &jvp.temperature_k) + dot(&walls, &jvp.mean_wall_temperatures_k) + dot(&heats, &jvp.heat_rates_w);
        let right = dot(&d.references_k, &vjp.references) + dot(&d.log_htc, &vjp.log_htc) + dot(&d.nodal_load_w, &vjp.nodal_load);
        close(left, right);
    });
}

#[test]
fn pinned_temperatures_are_not_reintroduced_as_a_derivative_lift() {
    with_cx(|cx| {
        let base = slab(cx, [20.0, 50.0], [330.0, 290.0], true);
        let mut d = base.zero_direction();
        d.references_k[0] = 1.0;
        d.log_htc[0] = 0.2;
        for &vertex in base.dofs.fixed() { d.nodal_load_w[vertex] = 1000.0; }
        let result = base.apply(cx, &d).unwrap();
        for &vertex in base.dofs.fixed() { assert_eq!(result.temperature_k[vertex], 0.0); }
        let resistance = 0.1/10.0 + 1.0/50.0;
        let heat = 40.0/resistance;
        let dq = (-1.0 + heat*0.2/50.0)/resistance;
        close(result.heat_rates_w[0], dq);
        let g = base.pullback(cx, &vec![1.0; base.primal.temperature.len()], &[0.0], &[0.0]).unwrap();
        for &vertex in base.dofs.fixed() { assert_eq!(g.nodal_load[vertex], 0.0); }
    });
}

#[test]
fn zero_derivatives_and_bad_vectors_are_explicit() {
    with_cx(|cx| {
        let base = slab(cx, [20.0, 50.0], [330.0, 290.0], false);
        let mut d = base.zero_direction();
        let zero = base.apply(cx, &d).unwrap();
        assert_eq!(zero.iterations, 0);
        assert!(zero.temperature_k.iter().all(|v| *v == 0.0));
        d.references_k.pop();
        assert!(base.apply(cx, &d).is_err());
        d = base.zero_direction(); d.log_htc[0] = f64::NAN;
        assert!(base.apply(cx, &d).is_err());
    });
}

#[test]
fn cancellation_refuses_before_a_tangent_or_adjoint_solve() {
    let base = with_cx(|cx| slab(cx, [20.0, 50.0], [330.0, 290.0], false));
    let gate = CancelGate::new(); gate.request();
    with_gate(&gate, |cx| {
        assert!(matches!(base.apply(cx, &base.zero_direction()), Err(ConductionError::Cancelled { .. })));
        assert!(matches!(base.pullback(cx, &vec![0.0; base.primal.temperature.len()], &[0.0;2], &[0.0;2]), Err(ConductionError::Cancelled { .. })));
    });
}
