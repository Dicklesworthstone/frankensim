//! Physical Newton integration: production tetrahedral conduction + live IR.
use fs_ascent::conduction_design::ConductivityNewtonModel;
use fs_ascent::{ReverseNewtonStudy,StopRule};
use fs_ascent::trust::KrylovOptions;
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::{ConductionMesh,ConductionProblem,ConductivityModel,LinearConfig,ScalarField,
    ThermalBc,ThermalBoundary,ThermalBoundaryBuilder,TEMPERATURE_DIMS};
use fs_conduction::fixtures::{unit_cube,on_box_face};
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_opt::reverse::physics::{PhysicsBinding,PhysicsError,PhysicsHessian,PhysicsModel};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{Manifold,NodeId,Problem,ProblemBuilder,ReverseProblem,Sense};
use fs_qty::Dims;

fn with_cx<R>(f: impl FnOnce(&CancelGate,&Cx<'_>)->R)->R {
    let gate=CancelGate::new();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed: 1,kernel_id: 81,tile: 0,iteration: 0},
            Budget::INFINITE,ExecMode::Deterministic);
        f(&gate,&cx)
    })
}
struct Fixture {mesh: ConductionMesh,boundary: ThermalBoundary,material: ConductivityModel,source: ScalarField}
impl Fixture {
    fn new(robin: bool)->Self {
        let (complex,positions)=unit_cube(2);let mesh=ConductionMesh::new(complex,positions).unwrap();
        let mut b=ThermalBoundaryBuilder::new(&mesh).region("fixed",|f| on_box_face(f.centroid[2],0.0),
            ThermalBc::dirichlet(300.0).unwrap()).unwrap();
        if robin {b=b.region("cooling",|f| on_box_face(f.centroid[0],1.0),ThermalBc::robin(8.0,290.0).unwrap()).unwrap();}
        Self {boundary:b.adiabatic_remainder().finish().unwrap(),mesh,
            material:ConductivityModel::constant_tensor([[3.0,0.5,0.25],[0.5,2.0,0.3],[0.25,0.3,1.5]]).unwrap(),
            source:ScalarField::Uniform(4000.0)}
    }
    fn design(&self,iterations:usize)->ConductivityDesign<'_> {
        ConductivityDesign::new(ConductionProblem {mesh:&self.mesh,boundary:&self.boundary,
            material:&self.material,source:&self.source,element_materials:None},
            LinearConfig {tolerance:1e-12,max_iterations:iterations,restart:30}).unwrap()
    }
}
fn close(a:&[f64],b:&[f64],tol:f64) {
    assert_eq!(a.len(),b.len());let scale=b.iter().map(|x|x.abs()).fold(1.0_f64,f64::max);
    assert!(a.iter().zip(b).all(|(a,b)|(a-b).abs()<tol*scale),"{a:?} != {b:?}");
}
fn problem(n:usize)->(Problem,NodeId) {
    let mut b=ProblemBuilder::new();let x=b.var("log rho",Manifold::Rn {dim:n as u32},Dims::NONE).unwrap();
    let xr=b.var_ref(x).unwrap();let reg=b.norm_sq(xr).unwrap();
    let coefficient=b.konst(0.05/n as f64,Dims::NONE).unwrap();let reg=b.mul(coefficient,reg).unwrap();
    let node=b.pde_residual("thermal",x,true,TEMPERATURE_DIMS).unwrap();
    let reference=b.konst(300.0,TEMPERATURE_DIMS).unwrap();let scale=b.konst(1000.0,TEMPERATURE_DIMS).unwrap();
    let rise=b.sub(node,reference).unwrap();let rise=b.div(rise,scale).unwrap();
    let fit=b.powi(rise,2).unwrap();let total=b.add(fit,reg).unwrap();
    b.objective(total,Sense::Minimize,1.0).unwrap();(b.finish(),node)
}
fn limits()->ReverseLimits {ReverseLimits {max_nodes:128,max_scalar_slots:4096}}

#[test]
fn log_coordinate_hessian_has_the_exact_uniform_scaling_curvature() {
    let f=Fixture::new(false);let d=f.design(2000);let n=d.parameter_count();
    let model=ConductivityNewtonModel::new(&d,"thermal",vec![1.0/d.dofs().n() as f64;d.dofs().n()],100).unwrap();
    let x=vec![0.2;n];let v=vec![1.0;n];
    with_cx(|_,cx| {
        let sample=model.value_gradient(&x,Some(cx)).unwrap();
        let hv=model.hessian_vector_product(&x,&v,Some(cx)).unwrap();
        close(&hv,&sample.gradient.iter().map(|g|-g).collect::<Vec<_>>(),2e-8);
        close(&[hv.iter().sum()],&[sample.value-300.0],2e-8);
        assert_eq!(model.reserved_solves(),4);
        assert_eq!(model.value_gradient(&x,Some(cx)).unwrap(),sample);
        assert_eq!(model.reserved_solves(),4,"cached primal must not solve again");
    });
}

#[test]
fn physical_ir_composition_matches_gradient_differences_after_trial_eviction() {
    let f=Fixture::new(true);let d=f.design(2000);let n=d.parameter_count();let (p,node)=problem(n);
    let model=ConductivityNewtonModel::new(&d,"thermal",vec![1.0/d.dofs().n() as f64;d.dofs().n()],200).unwrap();
    let oracle=ReverseProblem::new_with_physics(&p,limits(),&[PhysicsBinding {node,model:&model}],None).unwrap();
    let x:Vec<_>=(0..n).map(|i|(i%7) as f64*0.03).collect();
    let v:Vec<_>=(0..n).map(|i|((i%9) as f64-4.0)/9.0).collect();
    with_cx(|_,cx| {
        let tape=oracle.evaluate_cancellable(&[x.clone()],cx).unwrap();
        let hv=tape.objective_hessian_vector_product(&[v.clone()],Some(cx)).unwrap()[0].clone();
        let h=2e-5;
        let plus:Vec<_>=x.iter().zip(&v).map(|(x,v)|x+h*v).collect();
        let minus:Vec<_>=x.iter().zip(&v).map(|(x,v)|x-h*v).collect();
        let gp=oracle.evaluate_cancellable(&[plus],cx).unwrap().objective_gradient_cancellable(cx).unwrap();
        let gm=oracle.evaluate_cancellable(&[minus],cx).unwrap().objective_gradient_cancellable(cx).unwrap();
        let fd:Vec<_>=gp[0].iter().zip(&gm[0]).map(|(a,b)|(a-b)/(2.0*h)).collect();
        close(&hv,&fd,2e-6);
        let before=model.reserved_solves();
        assert_eq!(tape.objective_hessian_vector_product(&[v],Some(cx)).unwrap()[0],hv);
        assert_eq!(model.reserved_solves(),before+4,"old tape rebuilds its own point, then applies Hv");
    });
}

#[test]
fn exact_provider_budget_retries_without_refund_or_losing_cached_value() {
    let f=Fixture::new(true);let d=f.design(2000);let n=d.parameter_count();
    let model=ConductivityNewtonModel::new(&d,"thermal",vec![1.0;d.dofs().n()],2).unwrap();
    let x=vec![0.0;n];let y=vec![0.1;n];let v=vec![1.0;n];
    with_cx(|_,cx| {
        let sample=model.value_gradient(&x,Some(cx)).unwrap();
        assert!(matches!(model.hessian_vector_product(&x,&v,Some(cx)),Err(PhysicsError::Budget {used:2,limit:2,..})));
        assert_eq!(model.reserved_solves(),2);
        model.set_solve_limit(4);model.hessian_vector_product(&x,&v,Some(cx)).unwrap();
        assert!(matches!(model.value_gradient(&y,Some(cx)),Err(PhysicsError::Budget {used:4,limit:4,..})));
        assert_eq!(model.value_gradient(&x,Some(cx)).unwrap(),sample);
        assert!(model.hessian_vector_product(&x,&vec![0.0;n],Some(cx)).unwrap().iter().all(|&x|x==0.0));
        assert_eq!(model.reserved_solves(),4);
    });
    // A failing physical solve still consumes its reserved pair.
    let bad=f.design(1);let failing=ConductivityNewtonModel::new(&bad,"thermal",vec![1.0;bad.dofs().n()],4).unwrap();
    with_cx(|_,cx| {
        for used in [2,4] {
            assert!(matches!(failing.value_gradient(&x,Some(cx)),Err(PhysicsError::NotConverged { .. })));
            assert_eq!(failing.reserved_solves(),used);assert!(failing.last_native_error().is_some());
        }
        assert!(matches!(failing.value_gradient(&x,Some(cx)),Err(PhysicsError::Budget { .. })));
    });
}

#[test]
fn native_newton_decreases_the_physical_goal_and_replays_split_iterations() {
    let f=Fixture::new(true);let d=f.design(2000);let n=d.parameter_count();let (p,node)=problem(n);
    let model=ConductivityNewtonModel::new(&d,"thermal",vec![1.0/d.dofs().n() as f64;d.dofs().n()],4000).unwrap();
    let oracle=ReverseProblem::new_with_physics(&p,limits(),&[PhysicsBinding {node,model:&model}],None).unwrap();
    with_cx(|_,cx| {
        let mut whole=ReverseNewtonStudy::new(&oracle,&vec![0.0;n],Some(cx)).unwrap();
        whole.set_krylov_options(KrylovOptions::new(0.05,4).unwrap());let mut split=whole.clone();
        let initial=whole.snapshot().f;let rule=StopRule::GradNorm(1e-8);
        let end=whole.run(&rule,8,200,Some(cx)).unwrap();
        split.run(&rule,3,200,Some(cx)).unwrap();let replay=split.run(&rule,5,200,Some(cx)).unwrap();
        assert!(end.solution.f<initial*0.9,"physical objective should materially decrease");
        assert_eq!(end.solution.x,replay.solution.x);assert_eq!(end.solution.f.to_bits(),replay.solution.f.to_bits());
        assert_eq!(whole.history(),split.history());assert_eq!(whole.krylov_options(),split.krylov_options());
        assert!(whole.history().windows(2).all(|w|w[1]<=w[0]));
        assert!(end.solution.hv_evals>0 && model.reserved_solves()>2);
    });
}

#[test]
fn provider_exhaustion_preserves_newton_checkpoint_and_repair_resumes() {
    let f=Fixture::new(true);let d=f.design(2000);let n=d.parameter_count();let (p,node)=problem(n);
    let model=ConductivityNewtonModel::new(&d,"thermal",vec![1.0/d.dofs().n() as f64;d.dofs().n()],2).unwrap();
    let oracle=ReverseProblem::new_with_physics(&p,limits(),&[PhysicsBinding {node,model:&model}],None).unwrap();
    with_cx(|_,cx| {
        let mut study=ReverseNewtonStudy::new(&oracle,&vec![0.0;n],Some(cx)).unwrap();
        study.set_krylov_options(KrylovOptions::new(0.05,2).unwrap());
        let before=study.snapshot();let radius=study.radius();let history=study.history().to_vec();
        assert!(study.run(&StopRule::GradNorm(0.0),1,100,Some(cx)).is_err());
        assert_eq!(study.snapshot().x,before.x);assert_eq!(study.snapshot().f.to_bits(),before.f.to_bits());
        assert_eq!(study.radius(),radius);assert_eq!(study.history(),history);
        assert_eq!(study.snapshot().hv_evals,1,"attempted outer Hessian is retained");
        model.set_solve_limit(100);
        let result=study.run(&StopRule::GradNorm(0.0),2,100,Some(cx)).unwrap();
        assert!(result.solution.f<before.f);assert!(model.reserved_solves()>2);
    });
}

#[test]
fn cancellation_bad_shapes_and_exp_domain_fail_before_new_solve_reservations() {
    let f=Fixture::new(true);let d=f.design(2000);let n=d.parameter_count();
    let model=ConductivityNewtonModel::new(&d,"thermal",vec![1.0;d.dofs().n()],20).unwrap();
    with_cx(|gate,cx| {
        for x in [vec![],vec![f64::NAN;n],vec![1000.0;n],vec![-1000.0;n]] {
            assert!(model.value_gradient(&x,Some(cx)).is_err());
        }
        assert_eq!(model.reserved_solves(),0);
        let x=vec![0.0;n];let saved=model.value_gradient(&x,Some(cx)).unwrap();
        assert!(model.hessian_vector_product(&x,&[],Some(cx)).is_err());
        assert!(model.hessian_vector_product(&x,&vec![f64::INFINITY;n],Some(cx)).is_err());
        assert_eq!(model.reserved_solves(),2);
        gate.request();assert_eq!(model.value_gradient(&x,Some(cx)).unwrap_err(),PhysicsError::Cancelled);
        assert_eq!(model.hessian_vector_product(&x,&vec![1.0;n],Some(cx)).unwrap_err(),PhysicsError::Cancelled);
        assert_eq!(model.reserved_solves(),2);
        with_cx(|_,cx| assert_eq!(model.value_gradient(&x,Some(cx)).unwrap(),saved));
    });
    assert!(model.value_gradient(&vec![0.0;n],None).is_err());
}
