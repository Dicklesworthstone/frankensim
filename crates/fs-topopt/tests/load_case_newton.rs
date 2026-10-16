//! Independent structural load cases as exact live-IR SIMP constraints.
use std::ops::ControlFlow;

use fs_opt::reverse::physics::{PhysicsBinding, PhysicsHessian, PhysicsModel};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use fs_topopt::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget};
use fs_topopt::elasticity::second_order::ComplianceSolveOptions;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::simp_newton::SimpNewtonModel;
use fs_topopt::simp_second_order::SimpComplianceGoal;
use fs_topopt::{DensityElasticity, DensityFilter, DesignPipeline, SimpParams, SolveBudget, SolveControl};

fn fixture() -> (DesignPipeline, DensityElasticity, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mesh, mut positions) = fs_feec::kuhn_cube(2);
    for p in &mut positions { p[0] *= p[0]; }
    let elasticity = DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &|p| p[0] == 0.0);
    let mut a = vec![0.0; elasticity.n()];
    let mut b = a.clone();
    for (i, p) in positions.iter().enumerate() {
        if p[0] == 1.0 {
            a[3*i+2] = -0.2;
            b[3*i+1] = 0.13*(1.0+p[2]);
        }
    }
    let volumes = fs_feec::element_geometry(&mesh, &positions)
        .vol_signed.iter().map(|v| v.abs()).collect();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&mesh, &positions, 0.15),
        params: SimpParams { beta: 2.0, penal: 3.0, e_min: 0.01, eta: 0.5 },
    };
    (pipeline, elasticity, a, b, volumes)
}
fn options() -> ComplianceSolveOptions {
    ComplianceSolveOptions { tolerance: 1e-12, max_iterations: 4000 }
}
fn budget() -> ComplianceNewtonBudget {
    ComplianceNewtonBudget { solve_reservations: 100_000, linear_iterations: 10_000_000 }
}
fn close(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    let error = a.iter().zip(b).map(|(a,b)| (a-b).abs()).fold(0.0_f64, f64::max);
    let scale = b.iter().map(|v| v.abs()).fold(1.0_f64, f64::max);
    assert!(error <= tolerance*scale, "{error:e} / {scale:e} > {tolerance:e}");
}
fn add_scaled(out: &mut [f64], scale: f64, value: &[f64]) {
    for (o, v) in out.iter_mut().zip(value) { *o += scale*v; }
}

#[test]
fn independent_case_derivatives_sum_to_the_weighted_native_goal() {
    let (pipeline, elasticity, a, b, _) = fixture();
    let loads = [LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &b, weight: 1.3 }];
    let n = elasticity.cells();
    let rho: Vec<_> = (0..n).map(|i| 0.38 + 0.02*(i%9) as f64).collect();
    let direction: Vec<_> = (0..n).map(|i| 0.025*((i%7) as f64-3.0)).collect();
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let goal = SimpComplianceGoal::new(&pipeline, &elasticity, &rho, &loads, options(), &mut control).unwrap();

    assert_eq!(goal.load_case_count(), 2);
    let (v0, g0) = goal.load_case_value_gradient(0, &mut control).unwrap();
    let (v1, g1) = goal.load_case_value_gradient(1, &mut control).unwrap();
    close(&[goal.value()], &[0.7*v0 + 1.3*v1], 2e-11);
    let mut weighted_gradient = vec![0.0; n];
    add_scaled(&mut weighted_gradient, 0.7, &g0);
    add_scaled(&mut weighted_gradient, 1.3, &g1);
    close(goal.gradient(), &weighted_gradient, 2e-10);

    let h0 = goal.load_case_hessian_vector_product(0, &direction, &mut control).unwrap();
    let h1 = goal.load_case_hessian_vector_product(1, &direction, &mut control).unwrap();
    let aggregate = goal.hessian_vector_product(&direction, &mut control).unwrap();
    let mut weighted_hessian = vec![0.0; n];
    add_scaled(&mut weighted_hessian, 0.7, &h0);
    add_scaled(&mut weighted_hessian, 1.3, &h1);
    close(&aggregate, &weighted_hessian, 3e-9);

    // Independent first-order pipeline differences check each case, rather than
    // comparing two second-order routines sharing the same contractions.
    let eps = 5e-5;
    let plus: Vec<_> = rho.iter().zip(&direction).map(|(x,d)| x+eps*d).collect();
    let minus: Vec<_> = rho.iter().zip(&direction).map(|(x,d)| x-eps*d).collect();
    for (load, expected) in [
        (LoadCase { force: &a, weight: 1.0 }, &h0),
        (LoadCase { force: &b, weight: 1.0 }, &h1),
    ] {
        let (_, mut reference, _, _, _) = fixture();
        let gp = pipeline.multi_load_compliance_and_gradient(&mut reference, &plus, &[load]).gradient;
        let gm = pipeline.multi_load_compliance_and_gradient(&mut reference, &minus, &[load]).gradient;
        let finite_difference: Vec<_> = gp.iter().zip(&gm).map(|(p,m)| (p-m)/(2.0*eps)).collect();
        close(expected, &finite_difference, 4e-6);
    }
    assert!(goal.load_case_value_gradient(2, &mut control).is_err());
}

#[test]
fn live_case_views_share_primal_state_but_keep_independent_curvature() {
    let (pipeline, elasticity, a, b, volumes) = fixture();
    let loads = [LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &b, weight: 1.3 }];
    let model = SimpNewtonModel::new(&pipeline, &elasticity, &loads, &volumes, options(), budget()).unwrap();
    let aggregate = model.compliance("aggregate").unwrap();
    let case0 = model.load_case_compliance(0, "case-0").unwrap();
    let case1 = model.load_case_compliance(1, "case-1").unwrap();
    assert!(model.load_case_compliance(2, "bad").is_err());
    assert_eq!(case0.signature().residual_dims, COMPLIANCE_DIMS);
    let n = elasticity.cells();
    let x: Vec<_> = (0..n).map(|i| 0.4 + 0.015*(i%7) as f64).collect();
    let direction: Vec<_> = (0..n).map(|i| 0.02*((i%5) as f64-2.0)).collect();

    let aggregate_sample = aggregate.value_gradient(&x, None).unwrap();
    let prepared = model.work();
    assert_eq!(prepared.solve_reservations, loads.len()+3);
    let c0 = case0.value_gradient(&x, None).unwrap();
    let after0 = model.work();
    let c1 = case1.value_gradient(&x, None).unwrap();
    let after1 = model.work();
    // Per-case gradients reuse both primal displacements; only their individual
    // transpose-filter pullbacks are new numerical solves.
    assert_eq!(after0.solve_reservations, prepared.solve_reservations+1);
    assert_eq!(after1.solve_reservations, prepared.solve_reservations+2);
    close(&[aggregate_sample.value], &[0.7*c0.value+1.3*c1.value], 2e-11);
    let mut weighted = vec![0.0;n];
    add_scaled(&mut weighted,0.7,&c0.gradient);
    add_scaled(&mut weighted,1.3,&c1.gradient);
    close(&aggregate_sample.gradient,&weighted,2e-10);

    let before_zero = model.work();
    assert_eq!(case0.hessian_vector_product(&x,&vec![0.0;n],None).unwrap(),vec![0.0;n]);
    assert_eq!(model.work(),before_zero);
    let h0 = case0.hessian_vector_product(&x,&direction,None).unwrap();
    assert_eq!(model.work().solve_reservations, before_zero.solve_reservations+3);
    assert!(h0.iter().any(|v| v.abs()>0.0));
    // Reading another case after an HVP still does not rebuild the primal point.
    assert_eq!(case1.value_gradient(&x,None).unwrap(),c1);
}

#[test]
fn live_ir_can_impose_multiple_load_case_limits_without_force_cross_terms() {
    let (pipeline, elasticity, a, b, volumes) = fixture();
    let loads = [LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &b, weight: 1.3 }];
    let model = SimpNewtonModel::new(&pipeline, &elasticity, &loads, &volumes, options(), budget()).unwrap();
    let case0 = model.load_case_compliance(0, "case-0").unwrap();
    let case1 = model.load_case_compliance(1, "case-1").unwrap();
    let volume = model.volume("volume").unwrap();
    let n = elasticity.cells();
    let x = vec![0.47;n];
    let s0 = case0.value_gradient(&x,None).unwrap();
    let s1 = case1.value_gradient(&x,None).unwrap();
    let sv = volume.value_gradient(&x,None).unwrap();

    let mut builder = ProblemBuilder::new();
    let var = builder.var("raw density", Manifold::Rn { dim:n as u32 }, Dims::NONE).unwrap();
    let p0 = builder.pde_residual("case-0",var,true,COMPLIANCE_DIMS).unwrap();
    let p1 = builder.pde_residual("case-1",var,true,COMPLIANCE_DIMS).unwrap();
    let pv = builder.pde_residual("volume",var,true,Dims::NONE).unwrap();
    builder.objective(pv,Sense::Minimize,1.0).unwrap();
    let one = builder.konst(1.0,Dims::NONE).unwrap();
    let scale0 = builder.konst(s0.value,COMPLIANCE_DIMS).unwrap();
    let scale1 = builder.konst(s1.value,COMPLIANCE_DIMS).unwrap();
    let c0 = builder.sub(builder.div(p0,scale0).unwrap(),one).unwrap();
    let c1 = builder.sub(builder.div(p1,scale1).unwrap(),one).unwrap();
    builder.constraint(c0,ConstraintKind::LeZero,"load case 0").unwrap();
    builder.constraint(c1,ConstraintKind::LeZero,"load case 1").unwrap();
    let problem = builder.finish();
    let oracle = ReverseProblem::new_with_physics(&problem,
        ReverseLimits { max_nodes:256,max_scalar_slots:8192 },
        &[PhysicsBinding {node:p0,model:&case0},PhysicsBinding {node:p1,model:&case1},
          PhysicsBinding {node:pv,model:&volume}],None).unwrap();
    let tape = oracle.evaluate(&[x.clone()]).unwrap();
    let multipliers = [0.4,0.7];
    let gradient = tape.lagrangian_gradient().unwrap();
    let mut expected_gradient = sv.gradient.clone();
    add_scaled(&mut expected_gradient,multipliers[0]/s0.value,&s0.gradient);
    add_scaled(&mut expected_gradient,multipliers[1]/s1.value,&s1.gradient);
    close(&gradient[0],&expected_gradient,3e-10);

    let direction: Vec<_> = (0..n).map(|i|0.015*((i%7) as f64-3.0)).collect();
    let hv = tape.lagrangian_hessian_vector_product(&multipliers,&[direction.clone()],None).unwrap();
    let hv0 = case0.hessian_vector_product(&x,&direction,None).unwrap();
    let hv1 = case1.hessian_vector_product(&x,&direction,None).unwrap();
    let hvv = volume.hessian_vector_product(&x,&direction,None).unwrap();
    let mut expected_hessian = hvv;
    add_scaled(&mut expected_hessian,multipliers[0]/s0.value,&hv0);
    add_scaled(&mut expected_hessian,multipliers[1]/s1.value,&hv1);
    close(&hv[0],&expected_hessian,5e-9);
}


#[test]
fn zero_weight_aggregate_allows_volume_objective_with_only_case_constraints() {
    let (pipeline, elasticity, a, b, volumes) = fixture();
    let loads = [LoadCase { force: &a, weight: 0.0 }, LoadCase { force: &b, weight: 0.0 }];
    let model = SimpNewtonModel::new(&pipeline, &elasticity, &loads, &volumes, options(), budget()).unwrap();
    let aggregate = model.compliance("unused-aggregate").unwrap();
    let case0 = model.load_case_compliance(0, "case-0").unwrap();
    let case1 = model.load_case_compliance(1, "case-1").unwrap();
    let volume = model.volume("volume").unwrap();
    let n = elasticity.cells();
    let x = vec![0.5;n];

    let zero = aggregate.value_gradient(&x,None).unwrap();
    assert_eq!(zero.value,0.0);
    assert!(zero.gradient.iter().all(|&g| g == 0.0));
    let work = model.work();
    let c0 = case0.value_gradient(&x,None).unwrap();
    let c1 = case1.value_gradient(&x,None).unwrap();
    assert!(c0.value > 0.0 && c1.value > 0.0);
    assert_eq!(model.work().solve_reservations,work.solve_reservations+2);
    let before = model.work();
    assert_eq!(aggregate.hessian_vector_product(&x,&vec![1.0;n],None).unwrap(),vec![0.0;n]);
    assert_eq!(model.work(),before,"zero weighted aggregate must not spend tangent solves");

    let mut builder = ProblemBuilder::new();
    let var = builder.var("raw density",Manifold::Rn {dim:n as u32},Dims::NONE).unwrap();
    let v = builder.pde_residual("volume",var,true,Dims::NONE).unwrap();
    let p0 = builder.pde_residual("case-0",var,true,COMPLIANCE_DIMS).unwrap();
    let p1 = builder.pde_residual("case-1",var,true,COMPLIANCE_DIMS).unwrap();
    builder.objective(v,Sense::Minimize,1.0).unwrap();
    let one = builder.konst(1.0,Dims::NONE).unwrap();
    let s0 = builder.konst(c0.value,COMPLIANCE_DIMS).unwrap();
    let s1 = builder.konst(c1.value,COMPLIANCE_DIMS).unwrap();
    let q0 = builder.div(p0,s0).unwrap();
    let q1 = builder.div(p1,s1).unwrap();
    let r0 = builder.sub(q0,one).unwrap();
    let r1 = builder.sub(q1,one).unwrap();
    builder.constraint(r0,ConstraintKind::LeZero,"case 0 cap").unwrap();
    builder.constraint(r1,ConstraintKind::LeZero,"case 1 cap").unwrap();
    let problem = builder.finish();
    let oracle = ReverseProblem::new_with_physics(&problem,
        ReverseLimits {max_nodes:256,max_scalar_slots:8192},
        &[PhysicsBinding {node:v,model:&volume},PhysicsBinding {node:p0,model:&case0},
          PhysicsBinding {node:p1,model:&case1}],None).unwrap();
    let tape = oracle.evaluate(&[x]).unwrap();
    assert_eq!(tape.constraint_values(),&[0.0,0.0]);
    assert!(tape.objective_value().is_finite());
}
