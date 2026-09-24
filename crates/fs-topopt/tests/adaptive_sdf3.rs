use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::elastic3::{ElasticityOptions3,adaptive::AdaptiveElasticity3};
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveControl,SolveBudget,SolveProgress,EvaluationStop,MultiLoadOcOptions,MultiLoadOcTermination};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3,controlled_sdf3_optimality_criteria,inherit_raw_densities3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {p[2]-0.73}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval {
        let d=if a==HeightAxis::Z {1.0}else{0.0};Interval::new(d,d)
    }
}
fn tree()->Octree3 {
    let t=Octree3::uniform(1,4,1000).unwrap();let mark=*t.leaves().iter().find(|c|c.index()==[0,0,1]).unwrap();
    t.refined(&[mark],||ControlFlow::Continue(())).unwrap()
}
fn fixture(t:&Octree3,p:SimpParams)->(CutDensityStudy3<AdaptiveElasticity3>,Vec<f64>) {
    let mut poll=|_|ControlFlow::Continue(());
    let mut qc=QuadratureControl3::new(QuadratureOptions3{depth:1,..Default::default()},&mut poll).unwrap();
    let op=AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),t,&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|p|p[0]==0.0,ElasticityOptions3::default(),&mut qc).unwrap();
    let force=op.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(op,0.15,p),force)
}
fn options()->MultiLoadOcOptions {MultiLoadOcOptions{max_iterations:3,change_tolerance:0.0,..Default::default()}}
#[test]
fn adaptive_full_chain_coordinate_gradients_at_three_models() {
    for (penal,beta) in [(1.0,1.0),(3.0,2.0),(3.0,8.0)] {
        let (mut study,force)=fixture(&tree(),SimpParams{penal,beta,..Default::default()});
        let opposite:Vec<f64>=force.iter().map(|f|-f).collect();
        let loads=[LoadCase{force:&force,weight:0.3},LoadCase{force:&opposite,weight:0.7}];
        let rho:Vec<f64>=(0..study.cells()).map(|i|0.35+0.02*i as f64).collect();
        let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
        let exact=study.evaluate(&rho,&loads,&mut c).unwrap();assert!(exact.objective.compliance>0.0);
        for i in 0..rho.len(){
            let mut a=rho.clone();let mut b=rho.clone();a[i]+=1e-4;b[i]-=1e-4;
            let plus=study.evaluate(&a,&loads,&mut c).unwrap();let minus=study.evaluate(&b,&loads,&mut c).unwrap();
            let fd=(plus.objective.compliance-minus.objective.compliance)/2e-4;
            let vd=(plus.volume_fraction-minus.volume_fraction)/2e-4;
            assert!((fd-exact.objective.gradient[i]).abs()/fd.abs()<5e-4,"p={penal} beta={beta} cell={i}");
            assert!((vd-exact.volume_gradient[i]).abs()/vd.abs()<5e-5);
        }
    }
}
#[test]
fn optimize_refine_transfer_and_resolve_without_reusing_old_fields() {
    let first=Octree3::uniform(1,4,1000).unwrap();let (mut study,force)=fixture(&first,SimpParams::default());
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let rho=vec![0.5;study.cells()];
    let run=controlled_sdf3_optimality_criteria(&mut study,&[LoadCase{force:&force,weight:1.0}],&rho,options(),&mut c);
    assert!(run.history.len()>1);let before=c.work();
    let mark=*first.leaves().iter().find(|l|l.index()==[0,0,1]).unwrap();
    let refined=first.refined(&[mark],||ControlFlow::Continue(())).unwrap();
    let (mut next,force)=fixture(&refined,SimpParams::default());
    let seed=inherit_raw_densities3(study.operator().leaves(),&run.rho,next.operator().leaves(),&mut c).unwrap();
    assert!(next.cells()>study.cells());let seed=next.feasible_start(&seed,0.5,1e-8,&mut c).unwrap();
    let next_run=controlled_sdf3_optimality_criteria(&mut next,&[LoadCase{force:&force,weight:1.0}],&seed,options(),&mut c);
    assert!(next_run.history.len()>1);assert!(next_run.history.last().unwrap().compliance<next_run.history[0].compliance);
    assert!(next_run.history.iter().all(|row|row.volume_fraction<=0.5+1e-8));
    assert!(next_run.work.linear_iterations>before.linear_iterations);
    assert!(next.operator().physical_nodes().len()>next.operator().nodes().len());
    let fields=next.operator().physical_displacements(&next_run.displacements[0]).unwrap();
    assert_eq!(fields.len(),3*next.operator().physical_nodes().len());
}
#[test]
fn volume_restoration_is_real_and_missing_ancestors_do_not_zero_fill() {
    let (study,_)=fixture(&tree(),SimpParams::default());let original=study.operator().scales().to_vec();
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let restored=study.feasible_start(&vec![1.0;study.cells()],0.4,1e-8,&mut c).unwrap();
    assert!(restored.iter().all(|r|*r<1.0));
    assert!(study.volume_and_gradient(&restored,&mut c).unwrap().0<=0.4+1e-8);
    assert_eq!(study.operator().scales(),original);
    assert!(matches!(inherit_raw_densities3(&[],&[],study.operator().leaves(),&mut c),Err(EvaluationStop::Breakdown{stage:"sdf3-missing-source-ancestor"})));
}
#[test]
fn adaptive_trial_cancellation_preserves_the_accepted_state() {
    let (mut study,force)=fixture(&tree(),SimpParams::default());let rho=vec![0.5;study.cells()];
    let mut evaluations=0;let mut poll=|p:SolveProgress| {
        if p.stage=="sdf3-evaluation"{evaluations+=1;}
        if evaluations==2&&p.stage=="sdf3-elasticity"&&p.solve_iterations>0 {ControlFlow::Break(())}else{ControlFlow::Continue(())}
    };
    let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let run=controlled_sdf3_optimality_criteria(&mut study,&[LoadCase{force:&force,weight:1.0}],&rho,options(),&mut c);
    assert_eq!(run.termination,MultiLoadOcTermination::Cancelled);assert_eq!(run.history.len(),1);assert_eq!(run.rho,rho);
    let scales=study.operator().scales().to_vec();
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let exact=study.evaluate(&rho,&[LoadCase{force:&force,weight:1.0}],&mut c).unwrap();
    assert_eq!(run.displacements,exact.objective.displacements);assert_eq!(study.operator().scales(),scales);
}
#[test]
fn adaptive_zero_budget_and_replay_use_the_same_driver() {
    let (mut study,force)=fixture(&tree(),SimpParams::default());let rho=vec![0.5;study.cells()];
    let run=|study:&mut CutDensityStudy3<AdaptiveElasticity3>,budget| {
        let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget{total_iterations:budget,..Default::default()},&mut poll);
        controlled_sdf3_optimality_criteria(study,&[LoadCase{force:&force,weight:1.0}],&rho,options(),&mut c)
    };
    let stop=run(&mut study,0);assert_eq!(stop.termination,MultiLoadOcTermination::LinearBudget);assert!(stop.history.is_empty());
    let a=run(&mut study,100000);let b=run(&mut study,100000);assert!(a.history.len()>1);
    assert_eq!(a.rho,b.rho);assert_eq!(a.displacements,b.displacements);assert_eq!(a.work,b.work);
}

#[test]
fn continuation_reuses_the_nonuniform_cut_space_and_checks_its_pullback() {
    use fs_topopt::{ContinuationTermination, GradientCheckOptions};
    use fs_topopt::sdf3::continuation::controlled_gradient_checked_sdf3_continuation;
    let (mut study, force) = fixture(&tree(), SimpParams::default());
    let leaves = study.operator().leaves().to_vec();
    let nodes = study.operator().nodes().to_vec();
    let rho = vec![0.5; study.cells()];
    let stages = [(1.0, 1.0), (3.0, 4.0)]
        .map(|(penal, beta)| SimpParams { penal, beta, ..Default::default() });
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_gradient_checked_sdf3_continuation(&mut study,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &stages,
        MultiLoadOcOptions { max_iterations: 1, ..options() }, GradientCheckOptions::default(), &mut control);
    assert_eq!(report.termination, ContinuationTermination::ScheduleComplete, "{report:?}");
    assert_eq!(report.stages.len(), 2);
    assert_eq!(study.operator().leaves(), leaves);
    assert_eq!(study.operator().nodes(), nodes);
    for stage in &report.stages {
        assert!(stage.gradient_check.as_ref().unwrap().passed());
        assert_eq!(stage.history.len(), 2);
        assert!(stage.history[1].compliance < stage.history[0].compliance);
        assert!(stage.history[1].volume_fraction <= 0.5 + options().volume_tolerance);
    }
}
