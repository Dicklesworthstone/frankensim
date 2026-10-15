//! Native SIMP physics, shared dual-output state, and live-IR consumption.
use std::ops::ControlFlow;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, TilePool};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsError, PhysicsHessian, PhysicsModel};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use fs_topopt::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget};
use fs_topopt::elasticity::second_order::ComplianceSolveOptions;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::simp_newton::SimpNewtonModel;
use fs_topopt::simp_second_order::SimpComplianceGoal;
use fs_topopt::{DensityElasticity, DensityFilter, DesignPipeline, SimpParams, SolveBudget, SolveControl};

fn fixture(beta: f64) -> (DesignPipeline, DensityElasticity, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mesh, mut positions) = fs_feec::kuhn_cube(2);
    // Unequal true cell volumes: catch accidental unweighted volume averages.
    for p in &mut positions { p[0] = p[0] * p[0]; }
    let el = DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &|p| p[0] == 0.0);
    let mut a = vec![0.0; el.n()]; let mut b = a.clone();
    for (i, p) in positions.iter().enumerate() {
        if p[0] == 1.0 { a[3*i+2] = -0.2; b[3*i+1] = 0.13*(1.0+p[2]); }
    }
    let volumes = fs_feec::element_geometry(&mesh, &positions).vol_signed.iter().map(|v| v.abs()).collect();
    let pipeline = DesignPipeline { filter: DensityFilter::new(&mesh, &positions, 0.15),
        params: SimpParams { beta, penal: 3.0, e_min: 0.01, eta: 0.5 } };
    (pipeline, el, a, b, volumes)
}
fn budget() -> ComplianceNewtonBudget {
    ComplianceNewtonBudget { solve_reservations: 20_000, linear_iterations: 2_000_000 }
}
fn options() -> ComplianceSolveOptions {
    ComplianceSolveOptions { tolerance: 1e-12, max_iterations: 2000 }
}
fn close(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    let error = a.iter().zip(b).map(|(a,b)| (a-b).abs()).fold(0.0, f64::max);
    let scale = b.iter().map(|b| b.abs()).fold(1.0, f64::max);
    assert!(error < tolerance*scale, "{error:e} / {scale:e} >= {tolerance:e}");
}
fn with_cx(test: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let pool = TilePool::for_host(1, 37); let gate = CancelGate::new_clock_free();
    pool.arena_pool().scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 37, kernel_id: 84, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        test(&gate, &cx);
    });
}

#[test]
fn both_outputs_share_preparation_and_match_native_derivatives() {
    for beta in [0.0, 2.0, 6.0] {
        let (pipeline, el, a, b, volumes) = fixture(beta); let n = el.cells();
        let loads = [LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &b, weight: 1.3 }];
        let model = SimpNewtonModel::new(&pipeline, &el, &loads, &volumes, options(), budget()).unwrap();
        let c = model.compliance("compliance").unwrap(); let v = model.volume("volume").unwrap();
        let x: Vec<_> = (0..n).map(|i| 0.38 + 0.02*(i%9) as f64).collect();
        let d: Vec<_> = (0..n).map(|i| 0.03*((i%5) as f64-2.0)).collect();
        let cs = c.value_gradient(&x, None).unwrap(); let prepared = model.work();
        let vs = v.value_gradient(&x, None).unwrap();
        assert_eq!(prepared, model.work());
        assert_eq!(prepared.solve_reservations, loads.len()+3);
        assert_eq!(prepared.started_solves, loads.len()+3);
        assert_eq!(c.signature().residual_dims, COMPLIANCE_DIMS);
        assert_eq!(v.signature().residual_dims, Dims::NONE);
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let native = SimpComplianceGoal::new(&pipeline, &el, &x, &loads, options(), &mut control).unwrap();
        let (vv, vg) = native.volume_value_gradient(&volumes, &mut control).unwrap();
        assert_eq!(cs.value, native.value()); assert_eq!(cs.gradient, native.gradient());
        assert_eq!(vs.value, vv); assert_eq!(vs.gradient, vg);
        let ch = c.hessian_vector_product(&x, &d, None).unwrap();
        let vh = v.hessian_vector_product(&x, &d, None).unwrap();
        assert_eq!(ch, native.hessian_vector_product(&d, &mut control).unwrap());
        assert_eq!(vh, native.volume_hessian_vector_product(&volumes, &d, &mut control).unwrap());
        let spent = model.work();
        assert_eq!(spent.solve_reservations, prepared.solve_reservations + loads.len()+2 + if beta == 0.0 { 0 } else { 2 });
        assert_eq!(c.hessian_vector_product(&x, &vec![0.0;n], None).unwrap(), vec![0.0;n]);
        assert_eq!(v.clone().value_gradient(&x, None).unwrap(), vs);
        assert_eq!(model.work(), spent);
        // Independent gradient differences, not a shared derivative implementation.
        let eps = 1e-4;
        let plus: Vec<_> = x.iter().zip(&d).map(|(x,d)| x+eps*d).collect();
        let minus: Vec<_> = x.iter().zip(&d).map(|(x,d)| x-eps*d).collect();
        for (goal, hv) in [(&c, &ch), (&v, &vh)] {
            let gp = goal.value_gradient(&plus, None).unwrap().gradient;
            let gm = goal.value_gradient(&minus, None).unwrap().gradient;
            let fd: Vec<_> = gp.iter().zip(&gm).map(|(a,b)| (a-b)/(2.0*eps)).collect();
            close(hv, &fd, 3e-6);
        }
        // An accepted point's Hessian cannot silently reuse the last trial state.
        assert_eq!(c.hessian_vector_product(&x, &d, None).unwrap(), ch);
        assert_eq!(v.value_gradient(&x, None).unwrap(), vs);
        assert!(el.moduli.iter().all(|&m| m == 1.0));
    }
}

#[test]
fn volume_is_projected_material_not_simp_stiffness() {
    let (pipeline, el, a, _, volumes) = fixture(0.0); let n = el.cells();
    let model = SimpNewtonModel::new(&pipeline, &el, &[LoadCase { force: &a, weight: 1.0 }],
        &volumes, options(), budget()).unwrap();
    let v = model.volume("volume").unwrap();
    let x = vec![0.4;n]; let sample = v.value_gradient(&x, None).unwrap();
    close(&[sample.value, sample.gradient.iter().sum()], &[0.4, 1.0], 1e-11);
    let work = model.work();
    assert_eq!(v.hessian_vector_product(&x, &vec![1.0;n], None).unwrap(), vec![0.0;n]);
    assert_eq!(model.work(), work);
    let raw: Vec<_> = (0..n).map(|i| 0.3+0.03*(i%9) as f64).collect();
    let filtered = pipeline.filter.apply(&raw);
    let total: f64 = volumes.iter().sum();
    let expected: f64 = volumes.iter().zip(&filtered).map(|(w,r)| w*r/total).sum();
    close(&[v.value_gradient(&raw, None).unwrap().value], &[expected], 1e-12);
}

#[test]
fn admission_and_failed_trials_preserve_shared_cache_and_actual_work() {
    let (pipeline, el, a, _, volumes) = fixture(2.0); let n = el.cells();
    let model = SimpNewtonModel::new(&pipeline, &el, &[LoadCase { force: &a, weight: 1.0 }],
        &volumes, options(), ComplianceNewtonBudget { solve_reservations: 3, ..budget() }).unwrap();
    let c = model.compliance("c").unwrap(); let v = model.volume("v").unwrap();
    let x = vec![0.5;n];
    assert!(matches!(v.value_gradient(&x, None), Err(PhysicsError::Budget { .. })));
    assert_eq!(model.work().solve_reservations, 0); assert_eq!(model.work().started_solves, 0);
    model.set_budget(budget());
    let accepted = c.value_gradient(&x, None).unwrap(); let volume = v.value_gradient(&x, None).unwrap();
    let before = model.work();
    model.set_budget(ComplianceNewtonBudget { linear_iterations: before.linear_iterations, ..budget() });
    assert!(matches!(c.value_gradient(&vec![0.55;n], None), Err(PhysicsError::Budget { .. })));
    let failed = model.work();
    assert_eq!(failed.solve_reservations, before.solve_reservations+4);
    assert_eq!(failed.started_solves, before.started_solves+1);
    assert_eq!(failed.linear_iterations, before.linear_iterations);
    assert!(model.last_native_error().is_some());
    assert_eq!(c.value_gradient(&x, None).unwrap(), accepted);
    assert_eq!(v.value_gradient(&x, None).unwrap(), volume); assert_eq!(model.work(), failed);
    model.set_budget(budget());
    assert!(v.value_gradient(&vec![0.55;n], None).unwrap().value > volume.value);
    assert!(model.work().linear_iterations > failed.linear_iterations);
    assert_eq!(c.value_gradient(&x, None).unwrap(), accepted);
}

#[test]
fn cancellation_and_malformed_directions_spend_no_work() {
    with_cx(|gate, cx| {
        let (pipeline, el, a, _, volumes) = fixture(2.0); let n = el.cells();
        let model = SimpNewtonModel::new(&pipeline, &el, &[LoadCase { force: &a, weight: 1.0 }],
            &volumes, options(), budget()).unwrap();
        assert!(model.compliance("").is_err());
        let c = model.compliance("c").unwrap(); let v = model.volume("v").unwrap();
        for bad in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
            assert!(c.value_gradient(&vec![bad;n], Some(cx)).is_err());
        }
        assert!(v.hessian_vector_product(&vec![0.5;n], &vec![f64::NAN;n], Some(cx)).is_err());
        assert!(v.hessian_vector_product(&vec![0.5;n], &[1.0], Some(cx)).is_err());
        assert_eq!(model.work().solve_reservations, 0);
        let sample = c.value_gradient(&vec![0.5;n], Some(cx)).unwrap(); let work = model.work();
        gate.request();
        assert_eq!(c.value_gradient(&vec![0.5;n], Some(cx)), Err(PhysicsError::Cancelled));
        assert_eq!(v.hessian_vector_product(&vec![0.5;n], &vec![0.0;n], Some(cx)), Err(PhysicsError::Cancelled));
        assert_eq!(model.work(), work);
        assert_eq!(c.value_gradient(&vec![0.5;n], None).unwrap(), sample);
    });
}

#[test]
fn live_ir_uses_both_physical_nodes_with_only_one_preparation() {
    let (pipeline, el, a, b, volumes) = fixture(2.0); let n = el.cells();
    let model = SimpNewtonModel::new(&pipeline, &el,
        &[LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &b, weight: 1.3 }],
        &volumes, options(), budget()).unwrap();
    let c = model.compliance("c").unwrap(); let v = model.volume("v").unwrap();
    let mut builder = ProblemBuilder::new();
    let var = builder.var("raw densities", Manifold::Rn { dim: n as u32 }, Dims::NONE).unwrap();
    let cn = builder.pde_residual("c", var, true, COMPLIANCE_DIMS).unwrap();
    let vn = builder.pde_residual("v", var, true, Dims::NONE).unwrap();
    let unit = builder.konst(1.0, COMPLIANCE_DIMS).unwrap();
    let objective = builder.div(cn, unit).unwrap();
    builder.objective(objective, Sense::Minimize, 1.0).unwrap();
    let cap = builder.konst(0.5, Dims::NONE).unwrap(); let residual = builder.sub(vn, cap).unwrap();
    builder.constraint(residual, ConstraintKind::LeZero, "physical volume").unwrap();
    let problem = builder.finish();
    let oracle = ReverseProblem::new_with_physics(&problem,
        ReverseLimits { max_nodes: 128, max_scalar_slots: 4096 },
        &[PhysicsBinding { node: cn, model: &c }, PhysicsBinding { node: vn, model: &v }], None).unwrap();
    with_cx(|_, cx| {
        let x = vec![0.45;n]; let tape = oracle.evaluate_cancellable(&[x.clone()], cx).unwrap();
        assert_eq!(model.work().solve_reservations, 5);
        let lambda = 0.7;
        let gradient = tape.lagrangian_gradient_cancellable(&[lambda], cx).unwrap();
        let cg = c.value_gradient(&x, Some(cx)).unwrap().gradient;
        let vg = v.value_gradient(&x, Some(cx)).unwrap().gradient;
        let expected: Vec<_> = cg.iter().zip(&vg).map(|(c,v)| c+lambda*v).collect();
        close(&gradient[0], &expected, 1e-12);
        assert_eq!(model.work().solve_reservations, 5);
    });
}
