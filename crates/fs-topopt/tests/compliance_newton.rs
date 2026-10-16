//! Real 3-D elasticity in both objective and constraint roles.
use std::ops::ControlFlow;
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions, AugmentedNewtonStop, AugmentedNewtonStudy};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, TilePool};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsError, PhysicsHessian, PhysicsModel, PhysicsSample, PhysicsSignature};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use fs_topopt::{DensityElasticity, SolveBudget, SolveControl};
use fs_topopt::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget, ComplianceNewtonModel};
use fs_topopt::elasticity::second_order::{ComplianceGoal, ComplianceSolveOptions};
use fs_topopt::pipeline::LoadCase;

fn fixture() -> (DensityElasticity, Vec<f64>, Vec<f64>) {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let el = DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &|p| p[0] == 0.0);
    let mut a = vec![0.0; el.n()]; let mut b = a.clone();
    for (i, p) in positions.iter().enumerate() {
        if p[0] == 1.0 { a[3*i+2] = -0.2; b[3*i+1] = 0.13*(1.0+p[2]); }
    }
    (el,a,b)
}
fn budget() -> ComplianceNewtonBudget {
    ComplianceNewtonBudget { solve_reservations: 20_000, linear_iterations: 2_000_000 }
}
fn optimizer_budget() -> AugmentedNewtonBudget {
    AugmentedNewtonBudget { objective_evals: 300, hessian_products: 1000, multiplier_updates: 50 }
}
fn options() -> ComplianceSolveOptions {
    ComplianceSolveOptions { tolerance: 1e-12, max_iterations: 2000 }
}
fn with_cx(test: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let pool = TilePool::for_host(1, 37);
    let gate = CancelGate::new_clock_free();
    pool.arena_pool().scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 37, kernel_id: 82, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        test(&gate, &cx);
    });
}
fn close(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    let error = a.iter().zip(b).map(|(a,b)| (a-b).abs()).fold(0.0, f64::max);
    let scale = b.iter().map(|b| b.abs()).fold(1.0, f64::max);
    assert!(error < tolerance*scale, "{error:e} / {scale:e} >= {tolerance:e}");
}

#[test]
fn log_hessian_has_coordinate_curvature_and_survives_trial_eviction() {
    let (el,a,b) = fixture(); let n = el.cells();
    let model = ComplianceNewtonModel::new(&el,"native",
        &[LoadCase { force:&a,weight:0.7 },LoadCase { force:&b,weight:1.3 }],options(),budget()).unwrap();
    let x: Vec<f64> = (0..n).map(|e| 0.04*((e%7) as f64-3.0)).collect();
    let v: Vec<f64> = (0..n).map(|e| 0.03*((e%5) as f64-2.0)).collect();
    let sample = model.value_gradient(&x,None).unwrap();
    let first_work = model.work();
    assert_eq!(sample,model.value_gradient(&x,None).unwrap()); assert_eq!(first_work,model.work());
    let hv = model.hessian_vector_product(&x,&v,None).unwrap();
    let uniform = model.hessian_vector_product(&x,&vec![1.0;n],None).unwrap();
    // C(x+t*1) = exp(-t) C(x); this fails if the diagonal g_m term is omitted.
    close(&uniform,&sample.gradient.iter().map(|g| -*g).collect::<Vec<_>>(),1e-10);
    close(&[sample.gradient.iter().sum()],&[-sample.value],1e-10);
    let eps = 1e-4;
    let plus: Vec<f64> = x.iter().zip(&v).map(|(x,v)| x+eps*v).collect();
    let minus: Vec<f64> = x.iter().zip(&v).map(|(x,v)| x-eps*v).collect();
    let gp = model.value_gradient(&plus,None).unwrap();
    let gm = model.value_gradient(&minus,None).unwrap();
    let difference: Vec<f64> = gp.gradient.iter().zip(&gm.gradient).map(|(a,b)| (a-b)/(2.0*eps)).collect();
    close(&hv,&difference,2e-6);
    assert_eq!(model.hessian_vector_product(&x,&v,None).unwrap(),hv);
    let before = model.work();
    assert_eq!(model.hessian_vector_product(&x,&vec![0.0;n],None).unwrap(),vec![0.0;n]);
    assert_eq!(model.work(),before);
    assert_eq!(model.work().solve_reservations,model.work().started_solves);
}

#[test]
fn rejected_preparation_retains_cache_and_spent_attempts_then_recovers() {
    let (el,a,_) = fixture(); let n = el.cells();
    let model = ComplianceNewtonModel::new(&el,"native",&[LoadCase { force:&a,weight:1.0 }],options(),budget()).unwrap();
    let sample = model.value_gradient(&vec![0.0;n],None).unwrap();
    let before = model.work();
    model.set_budget(ComplianceNewtonBudget { linear_iterations:before.linear_iterations,..budget() });
    assert!(matches!(model.value_gradient(&vec![0.1;n],None),Err(PhysicsError::Budget { .. })));
    let stopped = model.work();
    assert_eq!(stopped.solve_reservations,before.solve_reservations+1);
    assert_eq!(stopped.started_solves,before.started_solves+1);
    assert_eq!(stopped.linear_iterations,before.linear_iterations);
    assert!(model.last_native_error().is_some());
    assert_eq!(sample,model.value_gradient(&vec![0.0;n],None).unwrap()); assert_eq!(model.work(),stopped);
    model.set_budget(budget());
    let repaired = model.value_gradient(&vec![0.1;n],None).unwrap();
    close(&[repaired.value],&[sample.value/fs_math::det::exp(0.1)],1e-10);
    assert_eq!(model.work().solve_reservations,stopped.solve_reservations+1);
    assert!(model.work().linear_iterations>stopped.linear_iterations);
}

#[test]
fn cancellation_and_invalid_coordinates_cannot_publish_or_spend_work() {
    with_cx(|gate,cx| {
        let (el,a,_) = fixture(); let n = el.cells();
        let model = ComplianceNewtonModel::new(&el,"native",&[LoadCase { force:&a,weight:1.0 }],options(),budget()).unwrap();
        for bad in [f64::NAN,f64::INFINITY,1000.0,-1000.0] {
            assert!(model.value_gradient(&vec![bad;n],Some(cx)).is_err());
        }
        assert_eq!(model.work().solve_reservations,0);
        let sample = model.value_gradient(&vec![0.0;n],Some(cx)).unwrap();
        let before = model.work(); gate.request();
        assert_eq!(model.value_gradient(&vec![0.0;n],Some(cx)),Err(PhysicsError::Cancelled));
        assert_eq!(model.hessian_vector_product(&vec![0.0;n],&vec![1.0;n],Some(cx)),Err(PhysicsError::Cancelled));
        assert_eq!(model.work(),before);
        // Explicitly disabling cancellation does not alter the scientific cache.
        assert_eq!(model.value_gradient(&vec![0.0;n],None).unwrap(),sample);
    });
}

// Test-only uniform reduction; all derivatives still run the native 3-D model.
#[derive(Debug)]
struct Uniform<'a,'m> { native:&'a ComplianceNewtonModel<'m>, n:usize }
impl PhysicsModel for Uniform<'_, '_> {
    fn signature(&self)->PhysicsSignature<'_> {
        PhysicsSignature { study:"uniform",parameter_count:1,parameter_dims:Dims::NONE,residual_dims:COMPLIANCE_DIMS }
    }
    fn value_gradient(&self,x:&[f64],cx:Option<&Cx<'_>>)->Result<PhysicsSample,PhysicsError> {
        if x.len()!=1 { return Err(PhysicsError::InvalidInput("uniform point shape")); }
        let s=self.native.value_gradient(&vec![x[0];self.n],cx)?;
        Ok(PhysicsSample { value:s.value,gradient:vec![s.gradient.iter().sum()] })
    }
    fn second_order(&self)->Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for Uniform<'_, '_> {
    fn hessian_vector_product(&self,x:&[f64],v:&[f64],cx:Option<&Cx<'_>>)->Result<Vec<f64>,PhysicsError> {
        if x.len()!=1 || v.len()!=1 { return Err(PhysicsError::InvalidInput("uniform direction shape")); }
        Ok(vec![self.native.hessian_vector_product(&vec![x[0];self.n],&vec![v[0];self.n],cx)?.iter().sum()])
    }
}

#[test]
fn native_compliance_objective_respects_stiffness_allowance_and_split_continuation() {
    with_cx(|_,cx| {
        let (el,a,bload) = fixture();
        let loads=[LoadCase { force:&a,weight:0.7 },LoadCase { force:&bload,weight:1.3 }];
        let native=ComplianceNewtonModel::new(&el,"native",&loads,options(),budget()).unwrap();
        let model=Uniform { native:&native,n:el.cells() };
        let base=model.value_gradient(&[0.0],Some(cx)).unwrap().value;
        let mut b=ProblemBuilder::new();let var=b.var("log stiffness",Manifold::Rn {dim:1},Dims::NONE).unwrap();
        let r=b.var_ref(var).unwrap();let x=b.component(r,0).unwrap();let modulus=b.exp(x).unwrap();
        let node=b.pde_residual("uniform",var,true,COMPLIANCE_DIMS).unwrap();
        let scale=b.konst(base,COMPLIANCE_DIMS).unwrap();let objective=b.div(node,scale).unwrap();
        b.objective(objective,Sense::Minimize,1.0).unwrap();let cap=b.konst(1.5,Dims::NONE).unwrap();
        let residual=b.sub(modulus,cap).unwrap();b.constraint(residual,ConstraintKind::LeZero,"stiffness allowance").unwrap();
        let p=b.finish();let oracle=ReverseProblem::new_with_physics(&p,ReverseLimits {max_nodes:128,max_scalar_slots:4096},
            &[PhysicsBinding {node,model:&model}],None).unwrap();
        let mut study=AugmentedNewtonStudy::new(&oracle,&[0.0],AugmentedNewtonOptions::default(),Some(cx)).unwrap();
        let mut split=study.clone();let report=study.run(1e-5,150,optimizer_budget(),Some(cx)).unwrap();
        assert_eq!(report.stop,AugmentedNewtonStop::Converged,"{report:?}");assert!(report.kkt.within_tolerance(1e-5));
        let result_modulus=fs_math::det::exp(report.x[0]);assert!((result_modulus-1.5).abs()<5e-5);
        assert!((report.objective-2.0/3.0).abs()<3e-5);assert!(report.multipliers[0]>0.0);
        let mut callback=|_| ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget::default(),&mut callback);
        let fresh=ComplianceGoal::new(&el,&vec![result_modulus;el.cells()],&loads,options(),&mut control).unwrap();
        close(&[fresh.value()/base],&[report.objective],1e-11);
        split.run(1e-5,4,optimizer_budget(),Some(cx)).unwrap();
        assert_eq!(split.run(1e-5,146,optimizer_budget(),Some(cx)).unwrap(),report);
    });
}

#[test]
fn native_compliance_constraint_drives_minimum_stiffness_and_budget_repair() {
    with_cx(|_,cx| {
        let (el,a,_) = fixture();let loads=[LoadCase {force:&a,weight:1.0}];
        let native=ComplianceNewtonModel::new(&el,"native",&loads,options(),
            ComplianceNewtonBudget {solve_reservations:1,..budget()}).unwrap();
        let model=Uniform {native:&native,n:el.cells()};let base=model.value_gradient(&[0.0],Some(cx)).unwrap().value;
        let mut b=ProblemBuilder::new();let var=b.var("log stiffness",Manifold::Rn {dim:1},Dims::NONE).unwrap();
        let r=b.var_ref(var).unwrap();let x=b.component(r,0).unwrap();let modulus=b.exp(x).unwrap();
        b.objective(modulus,Sense::Minimize,1.0).unwrap();
        let node=b.pde_residual("uniform",var,true,COMPLIANCE_DIMS).unwrap();
        let scale=b.konst(base,COMPLIANCE_DIMS).unwrap();let normalized=b.div(node,scale).unwrap();
        let cap=b.konst(2.0/3.0,Dims::NONE).unwrap();let residual=b.sub(normalized,cap).unwrap();
        b.constraint(residual,ConstraintKind::LeZero,"compliance limit").unwrap();let p=b.finish();
        let oracle=ReverseProblem::new_with_physics(&p,ReverseLimits {max_nodes:128,max_scalar_slots:4096},
            &[PhysicsBinding {node,model:&model}],None).unwrap();
        let mut study=AugmentedNewtonStudy::new(&oracle,&[0.0],AugmentedNewtonOptions::default(),Some(cx)).unwrap();
        let before=study.inner().snapshot();let radius=study.inner().radius();
        assert!(study.run(1e-5,1,optimizer_budget(),Some(cx)).is_err());
        assert_eq!(study.inner().snapshot().x,before.x);assert_eq!(study.inner().radius(),radius);
        assert_eq!(native.work().solve_reservations,1);assert_eq!(study.inner().snapshot().hv_evals,1);
        native.set_budget(budget());let report=study.run(1e-5,200,optimizer_budget(),Some(cx)).unwrap();
        assert_eq!(report.stop,AugmentedNewtonStop::Converged,"{report:?}");assert!(report.kkt.within_tolerance(1e-5));
        assert!((report.objective-1.5).abs()<1e-4);
        let mut callback=|_| ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget::default(),&mut callback);
        let fresh=ComplianceGoal::new(&el,&vec![report.objective;el.cells()],&loads,options(),&mut control).unwrap();
        close(&[fresh.value()/base-2.0/3.0],&report.constraints,1e-10);
        let tape=oracle.evaluate_cancellable(&[report.x.clone()],cx).unwrap();
        let gradient=tape.lagrangian_gradient_cancellable(&report.multipliers,cx).unwrap();
        close(&[gradient[0][0].abs()],&[report.kkt.stationarity],1e-10);
        assert!(native.work().started_solves>1 && report.multiplier_updates>0);
    });
}
