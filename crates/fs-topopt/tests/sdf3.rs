//! Real 3-D implicit-cut physics; no substituted elasticity or filter callbacks.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{CutElasticity3, ElasticityOptions3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{EvaluationStop, MultiLoadOcOptions, MultiLoadOcTermination, SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3, controlled_sdf3_optimality_criteria};

#[path = "sdf3/continuation.rs"]
mod continuation;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {p[2]-0.73}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval {
        let d=if a==HeightAxis::Z {1.0} else {0.0};Interval::new(d,d)
    }
}
fn fixture(params:SimpParams)->(CutDensityStudy3,Vec<f64>,Vec<f64>) {
    let mut poll=|_|ControlFlow::Continue(());
    let mut qc=QuadratureControl3::new(QuadratureOptions3 {depth:1,..Default::default()},&mut poll).unwrap();
    let op=CutElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),[2;3],&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|p|p[0]==0.0,ElasticityOptions3::default(),&mut qc).unwrap();
    let y=op.body_load(&|_|[0.0,-1.0,0.0],||ControlFlow::Continue(())).unwrap();
    let z=op.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(op,0.15,params),y,z)
}
fn options()->MultiLoadOcOptions {MultiLoadOcOptions {max_iterations:5,change_tolerance:0.0,..Default::default()}}
#[test]
fn full_chain_compliance_and_volume_gradients_match_each_coordinate_difference() {
    for (penal,beta) in [(1.0,1.0),(3.0,2.0),(3.0,8.0)] {
        let (mut study,y,z)=fixture(SimpParams {penal,beta,..Default::default()});
        let loads=[LoadCase {force:&y,weight:0.3},LoadCase {force:&z,weight:0.7}];
        let rho:Vec<f64>=(0..study.cells()).map(|i|0.35+0.03*i as f64).collect();
        let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
        let exact=study.evaluate(&rho,&loads,&mut c).unwrap();
        for i in 0..rho.len() {
            let mut a=rho.clone();let mut b=rho.clone();let h=1e-4;a[i]+=h;b[i]-=h;
            let plus=study.evaluate(&a,&loads,&mut c).unwrap();let minus=study.evaluate(&b,&loads,&mut c).unwrap();
            let fd=(plus.objective.compliance-minus.objective.compliance)/(2.0*h);
            let vd=(plus.volume_fraction-minus.volume_fraction)/(2.0*h);
            assert!((fd-exact.objective.gradient[i]).abs()/fd.abs()<2e-4,"p={penal} beta={beta} cell={i}");
            assert!((vd-exact.volume_gradient[i]).abs()/vd.abs()<2e-5);
        }
    }
}
#[test]
fn real_cut_domain_optimization_descends_at_fixed_volume_without_remeshing() {
    let (mut study,y,z)=fixture(SimpParams::default());let rho=vec![0.5;study.cells()];
    let keys=study.operator().cell_keys();let nodes=study.operator().nodes().to_vec();
    let loads=[LoadCase {force:&y,weight:0.3},LoadCase {force:&z,weight:0.7}];
    let mut poll=|_|ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let report=controlled_sdf3_optimality_criteria(&mut study,&loads,&rho,options(),&mut control);
    assert!(report.history.len()>1,"{report:?}");
    assert!(report.history.last().unwrap().compliance<0.9*report.history[0].compliance);
    assert_eq!(study.operator().cell_keys(),keys);assert_eq!(study.operator().nodes(),nodes);
    for row in &report.history {assert!(row.volume_fraction<=0.5+1e-8);assert_eq!(row.case_compliances.len(),2);}
    for pair in report.history.windows(2) {assert!(pair[1].compliance<=pair[0].compliance);assert!(pair[1].max_change<=0.15+1e-14);}
    let final_eval=study.evaluate(&report.rho,&loads,&mut control).unwrap();
    assert_eq!(report.projected_rho,final_eval.projected_rho);
    assert_eq!(report.displacements,final_eval.objective.displacements);
    assert_eq!(report.history.last().unwrap().compliance.to_bits(),final_eval.objective.compliance.to_bits());
}
#[test]
fn opposite_independent_loads_cannot_cancel_compliance() {
    let (mut study,_,z)=fixture(SimpParams::default());let neg:Vec<f64>=z.iter().map(|f|-f).collect();
    let rho=vec![0.5;study.cells()];let mut poll=|_|ControlFlow::Continue(());
    let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let one=study.evaluate(&rho,&[LoadCase {force:&z,weight:1.0}],&mut c).unwrap();
    let both=study.evaluate(&rho,&[LoadCase {force:&z,weight:0.3},LoadCase {force:&neg,weight:0.7}],&mut c).unwrap();
    assert!(both.objective.compliance>0.0);
    assert!((both.objective.compliance-one.objective.compliance).abs()<1e-10*one.objective.compliance);
}
#[test]
fn cancellation_inside_trial_rolls_back_to_accepted_physical_fields() {
    let (mut study,y,z)=fixture(SimpParams::default());let rho=vec![0.5;study.cells()];
    let loads=[LoadCase {force:&y,weight:0.3},LoadCase {force:&z,weight:0.7}];
    let mut count=0;
    let mut poll=|p:SolveProgress| {
        if p.stage=="sdf3-evaluation" {count+=1;}
        if count==2 && p.stage=="sdf3-elasticity" && p.solve_iterations>0 {ControlFlow::Break(())}else{ControlFlow::Continue(())}
    };
    let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let report=controlled_sdf3_optimality_criteria(&mut study,&loads,&rho,options(),&mut c);
    assert_eq!(report.termination,MultiLoadOcTermination::Cancelled);
    assert_eq!(report.history.len(),1);assert_eq!(report.rho,rho);
    let before=study.operator().scales().to_vec();
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let accepted=study.evaluate(&rho,&loads,&mut c).unwrap();
    assert_eq!(study.operator().scales(),before);assert_eq!(report.displacements,accepted.objective.displacements);
}
#[test]
fn exact_shared_iteration_exhaustion_retains_no_unfunded_design() {
    let (mut study,y,z)=fixture(SimpParams::default());let rho=vec![0.5;study.cells()];
    let scales=study.operator().scales().to_vec();let loads=[LoadCase {force:&y,weight:0.3},LoadCase {force:&z,weight:0.7}];
    let mut poll=|_|ControlFlow::Continue(());
    let mut c=SolveControl::new(SolveBudget {total_iterations:1,..Default::default()},&mut poll);
    let report=controlled_sdf3_optimality_criteria(&mut study,&loads,&rho,options(),&mut c);
    assert_eq!(report.termination,MultiLoadOcTermination::LinearBudget);assert_eq!(report.work.linear_iterations,1);
    assert!(matches!(report.evaluation_stop,Some(EvaluationStop::TotalBudget{..})));
    assert!(report.history.is_empty()&&report.displacements.is_empty()&&report.projected_rho.is_empty());
    assert_eq!(study.operator().scales(),scales);
}
#[test]
fn identical_studies_replay_fields_history_and_work() {
    let (mut study,y,z)=fixture(SimpParams::default());let rho=vec![0.5;study.cells()];
    let loads=[LoadCase {force:&y,weight:0.3},LoadCase {force:&z,weight:0.7}];
    let run=|study:&mut CutDensityStudy3| {
        let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
        controlled_sdf3_optimality_criteria(study,&loads,&rho,MultiLoadOcOptions {max_iterations:2,..options()},&mut c)
    };
    let a=run(&mut study);let b=run(&mut study);
    assert!(a.history.len()>1);assert_eq!(a.termination,b.termination);assert_eq!(a.rho,b.rho);
    assert_eq!(a.displacements,b.displacements);assert_eq!(a.projected_rho,b.projected_rho);assert_eq!(a.work,b.work);
    assert_eq!(a.history.len(),b.history.len());
    for (a,b) in a.history.iter().zip(&b.history) {assert_eq!(a.compliance.to_bits(),b.compliance.to_bits());assert_eq!(a.volume_fraction.to_bits(),b.volume_fraction.to_bits());}
}
