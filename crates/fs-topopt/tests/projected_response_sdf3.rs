//! Real prescribed-motion experiments through the linear-storage optimizer.
use std::{cell::Cell,ops::ControlFlow};
use fs_ascent::projected_al::{ProjectedAlError,ProjectedAlOptions,ProjectedAlStop};
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3,ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveBudget,SolveControl,SolveProgress};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::{ProjectedResponseStudy3,ProjectedResponseOptions3,ResponseCase3,ResponseTarget3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64{(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{
        let x=Interval::new(lo[0],hi[0]);(x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval{
        if a==HeightAxis::X{Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0)}else{Interval::new(0.0,0.0)}
    }
}
fn motion(p:[f64;3],n:[f64;3])->[f64;3]{if n[0]>0.0{[0.02,0.0,0.0]}else{[0.0,0.0,0.007*p[1]]}}
fn fixture(level:u8)->(CutDensityStudy3<AdaptiveSolveSpace3>,Vec<f64>,Vec<f64>){
    let mut poll=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    let op=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),
        &Octree3::uniform(level,4,4096).unwrap(),&Slab,&IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap();
    let force=op.body_load(&|_|[0.002,0.0,-0.003],||ControlFlow::Continue(())).unwrap();
    let observation=op.body_load(&|p|[0.0,0.0,1.0+p[0]],||ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default()),force,observation)
}
fn options()->ProjectedResponseOptions3{
    ProjectedResponseOptions3{density_floor:0.1,optimizer:ProjectedAlOptions{tolerance:1e-9,..Default::default()},..Default::default()}
}
#[test]
fn g1_384_density_variables_exceed_even_the_hard_dense_sqp_cap(){
    let (mut study,force,q)=fixture(3);assert!(study.cells()>341);
    let target=[ResponseTarget3{q:&q,target:0.0,scale:0.02,weight:0.001}];
    let cases=[ResponseCase3{force:&force,prescribed:Some(&motion),targets:&target}];
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let mut opts=options();opts.response.volume_weight=1.0;opts.volume_cap=0.6;
    let rho=vec![0.45;study.cells()];let nodes=study.operator().elasticity().nodes().to_vec();
    let mut design=ProjectedResponseStudy3::new(&mut study,&cases,&rho,opts,&mut c).unwrap();
    let baseline=design.accepted().objective;let report=design.run(1).unwrap();
    assert_eq!(report.work.iterations,1);assert!(design.accepted().objective<baseline);
    assert!(design.point().iter().all(|x|*x>=0.1&&*x<=1.0));
    assert!(design.accepted().displacements[0].iter().any(|u|*u!=0.0));
    assert_eq!(design.study().operator().elasticity().nodes(),nodes);
}
#[test]
fn g3_fit_updates_use_the_existing_response_fields_and_exact_scaling(){
    let (mut study,force,q)=fixture(1);let zero=[ResponseTarget3{q:&q,target:0.0,scale:0.02,weight:1.0}];
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let target=study.evaluate_responses(&vec![0.35;study.cells()],&[ResponseCase3{force:&force,prescribed:Some(&motion),targets:&zero}],Default::default(),&mut c).unwrap().responses[0][0];
    let targets=[ResponseTarget3{target,..zero[0]}];let cases=[ResponseCase3{force:&force,prescribed:Some(&motion),targets:&targets}];
    let rho=vec![0.5;study.cells()];let mut opts=options();opts.objective_scale=2.0;
    let mut design=ProjectedResponseStudy3::new(&mut study,&cases,&rho,opts,&mut c).unwrap();
    let baseline=design.accepted().objective;let report=design.run(12).unwrap();
    assert!(design.accepted().objective<baseline);assert_eq!(report.objective,design.accepted().objective/2.0);
    let accepted=design.accepted().clone();drop(design);
    let again=study.evaluate_responses(&accepted.rho,&cases,opts.response,&mut c).unwrap();
    assert_eq!(again.displacements,accepted.displacements);assert_eq!(again.adjoints,accepted.adjoints);
    assert_eq!(again.projected_rho,accepted.projected_rho);assert_eq!(again.gradient,accepted.gradient);
}
#[test]
fn g5_segmented_runs_preserve_dual_spectral_and_physics_state(){
    let (mut a,force,q)=fixture(1);let (mut b,_,_)=fixture(1);
    let targets=[ResponseTarget3{q:&q,target:0.003,scale:0.02,weight:1.0}];
    let cases=[ResponseCase3{force:&force,prescribed:Some(&motion),targets:&targets}];let rho=vec![0.5;a.cells()];
    let mut pa=|_|ControlFlow::Continue(());let mut ca=SolveControl::new(SolveBudget::default(),&mut pa);
    let mut pb=|_|ControlFlow::Continue(());let mut cb=SolveControl::new(SolveBudget::default(),&mut pb);
    let mut da=ProjectedResponseStudy3::new(&mut a,&cases,&rho,options(),&mut ca).unwrap();
    let mut db=ProjectedResponseStudy3::new(&mut b,&cases,&rho,options(),&mut cb).unwrap();
    let ra=da.run(6).unwrap();db.run(2).unwrap();let rb=db.run(4).unwrap();
    assert_eq!(ra.stop,rb.stop);assert_eq!(ra.multiplier,rb.multiplier);assert_eq!(ra.penalty,rb.penalty);
    assert_eq!(da.point(),db.point());assert_eq!(da.accepted().displacements,db.accepted().displacements);
    assert_eq!(da.work(),db.work());assert_eq!(da.optimizer_work(),db.optimizer_work());
    let spent=db.work();db.run(0).unwrap();assert_eq!(db.work(),spent);
}
#[test]
fn g4_cancelled_trial_or_postaccept_poll_keeps_matching_physics_and_can_resume(){
    for stage in ["response-projected-evaluate","response-projected-accepted"] {
        let (mut study,force,q)=fixture(1);let targets=[ResponseTarget3{q:&q,target:0.003,scale:0.02,weight:1.0}];
        let cases=[ResponseCase3{force:&force,prescribed:Some(&motion),targets:&targets}];let rho=vec![0.5;study.cells()];
        let enabled=Cell::new(false);let mut poll=|p:SolveProgress|if enabled.get()&&p.stage==stage{ControlFlow::Break(())}else{ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
        let mut design=ProjectedResponseStudy3::new(&mut study,&cases,&rho,options(),&mut c).unwrap();
        enabled.set(true);assert!(design.run(1).is_err());enabled.set(false);
        if stage=="response-projected-evaluate"{assert_eq!(design.point(),rho);}else{assert_eq!(design.optimizer_work().iterations,1);}
        let count=design.evaluations();design.run(2).unwrap();assert!(design.evaluations()>=count);
        let accepted=design.accepted().clone();drop(design);
        let replay=study.evaluate_responses(&accepted.rho,&cases,options().response,&mut c).unwrap();
        assert_eq!(replay.displacements,accepted.displacements);
    }
}
#[test]
fn g0_dimension_refusal_precedes_physics_and_budget_limit_does_not_certify_feasibility(){
    let (mut study,force,q)=fixture(1);let targets=[ResponseTarget3{q:&q,target:0.003,scale:0.02,weight:1.0}];
    let cases=[ResponseCase3{force:&force,prescribed:Some(&motion),targets:&targets}];let rho=vec![0.8;study.cells()];
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let mut opts=options();opts.optimizer.max_dimension=4;
    assert!(matches!(ProjectedResponseStudy3::new(&mut study,&cases,&rho,opts,&mut c),Err(ProjectedAlError::Invalid(_))));
    assert_eq!(c.work().linear_solves,0);
    let mut opts=options();opts.optimizer.max_evaluations=1;opts.volume_cap=0.3;
    let mut design=ProjectedResponseStudy3::new(&mut study,&cases,&rho,opts,&mut c).unwrap();let report=design.run(100).unwrap();
    assert_eq!(report.stop,ProjectedAlStop::EvaluationLimit);assert_eq!(report.work.evaluations,1);
    assert!(design.constraint_violation()>0.0);assert!(!report.kkt.within_tolerance(opts.optimizer.tolerance));assert_eq!(design.point(),rho);
}
