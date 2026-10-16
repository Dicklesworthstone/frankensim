//! Actual 3-D stress quadrature, constrained fields and complete material adjoints.
use super::*;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{CutElasticity3, ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{SolveBudget, SolveProgress};
use std::cell::Cell;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { (p[0]-0.17)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], a: HeightAxis) -> Interval {
        if a == HeightAxis::X { Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0) }
        else { Interval::new(0.0,0.0) }
    }
}
fn motion(p: [f64; 3], _: [f64; 3]) -> [f64; 3] { [0.02*(p[0]-0.17),0.0,0.0] }
fn opposite(p: [f64; 3], n: [f64; 3]) -> [f64; 3] { motion(p,n).map(|v| -2.0*v) }
fn params() -> SimpParams { SimpParams { beta: 2.0, penal: 3.0, ..Default::default() } }
fn tree(mixed: bool) -> Octree3 {
    let t = Octree3::uniform(1,4,4096).unwrap();
    if mixed { t.refined(&[*t.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap() } else { t }
}
fn build(mixed: bool) -> AdaptiveSolveSpace3 {
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut p).unwrap();
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(
        HexCell::try_new([0.0;3],[1.0;3]).unwrap(), &tree(mixed), &Slab,
        &IsotropicElastic::new(1.0,0.0,1.0).unwrap(), &|_|false, &|_,_|true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut q).unwrap();
    AdaptiveSolveSpace3::jacobi(op, 100_000_000)
}
fn same(a: &StressEvaluation3, b: &StressEvaluation3) {
    assert_eq!(a.aggregate, b.aggregate); assert_eq!(a.gradient, b.gradient);
    assert_eq!(a.displacements, b.displacements); assert_eq!(a.adjoints, b.adjoints);
    assert_eq!(a.projected_rho, b.projected_rho); assert_eq!(a.volume_gradient, b.volume_gradient);
    assert_eq!(a.case_compliances, b.case_compliances);
    assert_eq!(a.sampled_physical_max, b.sampled_physical_max);
}
fn affine<O: Sdf3Elasticity>(mut s: CutDensityStudy3<O>) {
    let force = vec![0.0; s.operator.n()]; let loads = [LoadCase { force: &force, weight: 1.0 }];
    let motions: &StressMotions3<'_> = &[Some(&motion)];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let rho = vec![0.5; s.cells()]; let before = s.operator.scales().to_vec();
    let e = s.evaluate_stress_with_prescribed_motion(&rho, &loads, motions,
        StressOptions3::default(), StressMeasure3::NormalizedAverage, &mut c).unwrap();
    // E=1, nu=0, axial strain .02: reference von Mises=.02, relaxed=.5*.02.
    assert!((e.aggregate-0.01).abs()<1e-9);
    assert!((e.sampled_physical_max-0.02*e.scales[0]).abs()<1e-9);
    let analytic = 0.02*crate::filter::heaviside_derivative(0.5,params().beta,params().eta);
    assert!((e.gradient.iter().sum::<f64>()-analytic).abs()<1e-7*analytic);
    assert_eq!(e.case_compliances, vec![0.0]);
    assert_eq!(s.operator.scales(), before);
}

#[test]
fn g1_affine_motion_has_analytic_stress_and_uniform_density_derivative() {
    affine(CutDensityStudy3::new(build(false),0.15,params()));
    affine(CutDensityStudy3::new(build(true),0.15,params()));
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut p).unwrap();
    let op = CutElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),
        [2;3], &Slab, &IsotropicElastic::new(1.0,0.0,1.0).unwrap(), &|_|false, &|_,_|true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut q).unwrap();
    affine(CutDensityStudy3::new(op,0.15,params()));
}

#[test]
fn g3_complete_mixed_load_gradient_matches_fresh_solves_in_both_measures() {
    for mixed in [false,true] {
        let op = build(mixed);
        let f = op.elasticity().body_load(&|p|[0.001*(1.0+p[1]),0.0,-0.001], || ControlFlow::Continue(())).unwrap();
        let zero = vec![0.0; f.len()];
        let loads = [LoadCase { force:&f, weight:0.3 }, LoadCase { force:&zero, weight:0.7 }];
        let motions: &StressMotions3<'_> = &[Some(&motion),Some(&opposite)];
        let mut s = CutDensityStudy3::new(op,0.15,params());
        let rho: Vec<_> = (0..s.cells()).map(|i|0.35+0.02*(i%7) as f64).collect();
        let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        for measure in [StressMeasure3::NormalizedAverage,StressMeasure3::SampledPeakBound] {
            let e = s.evaluate_stress_with_prescribed_motion(&rho,&loads,motions,Default::default(),measure,&mut c).unwrap();
            for i in [0,s.cells()/2,s.cells()-1] {
                let mut a=rho.clone(); let mut b=rho.clone(); a[i]+=1e-4; b[i]-=1e-4;
                let a=s.evaluate_stress_with_prescribed_motion(&a,&loads,motions,Default::default(),measure,&mut c).unwrap();
                let b=s.evaluate_stress_with_prescribed_motion(&b,&loads,motions,Default::default(),measure,&mut c).unwrap();
                let fd=(a.aggregate-b.aggregate)/2e-4;
                assert!((fd-e.gradient[i]).abs() < 3e-4*fd.abs().max(e.gradient[i].abs()).max(1e-8));
            }
            let external: f64=f.iter().zip(&e.displacements[0]).map(|(f,u)|f*u).sum();
            assert_eq!(external,e.case_compliances[0]); assert_eq!(e.case_compliances[1],0.0);
        }
    }
}

#[test]
fn g0_motion_case_with_zero_weight_still_enters_sampled_peak_constraint() {
    let op=build(false); let f=vec![0.0;op.n()];
    let loads=[LoadCase{force:&f,weight:1.0},LoadCase{force:&f,weight:0.0}];
    let motions: &StressMotions3<'_>=&[Some(&motion),Some(&opposite)];
    let mut s=CutDensityStudy3::new(op,0.15,params());let rho=vec![0.5;s.cells()];
    let mut preparations=0;
    let mut p=|s: SolveProgress| { if s.stage=="sdf3-preconditioner-start" {preparations+=1;} ControlFlow::Continue(()) };
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let e=s.evaluate_stress_with_prescribed_motion(&rho,&loads,motions,Default::default(),StressMeasure3::SampledPeakBound,&mut c).unwrap();
    assert!((e.case_relaxed_max[1]-2.0*e.case_relaxed_max[0]).abs()<1e-8);
    assert!(e.aggregate>=e.case_relaxed_max[1]);assert!(e.aggregate>=e.case_physical_max[1]);
    drop(c); assert_eq!(preparations,1);
}

#[test]
fn g5_none_motions_preserve_the_original_path_and_work_exactly() {
    let op=build(true); let f=op.elasticity().body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();
    let loads=[LoadCase{force:&f,weight:1.0}]; let mut s=CutDensityStudy3::new(op,0.15,params());
    let rho=vec![0.5;s.cells()];
    let mut p=|_|ControlFlow::Continue(());let mut a=SolveControl::new(SolveBudget::default(),&mut p);
    let old=s.evaluate_stress_with_measure(&rho,&loads,Default::default(),StressMeasure3::SampledPeakBound,&mut a).unwrap();
    let mut p=|_|ControlFlow::Continue(());let mut b=SolveControl::new(SolveBudget::default(),&mut p);
    let new=s.evaluate_stress_with_prescribed_motion(&rho,&loads,&[None],Default::default(),StressMeasure3::SampledPeakBound,&mut b).unwrap();
    same(&old,&new);assert_eq!(a.work(),b.work());
}

#[test]
fn g4_cancellation_and_failed_motion_preserve_incoming_material() {
    let op=build(true); let f=vec![0.0;op.n()];let loads=[LoadCase{force:&f,weight:1.0}];
    let mut s=CutDensityStudy3::new(op,0.15,params());let rho=vec![0.4;s.cells()];
    let before=s.operator.scales().to_vec();let motions: &StressMotions3<'_>=&[Some(&motion)];
    for stage in ["sdf3-stress-motion-load","sdf3-stress-motion-pullback","sdf3-stress-publish"] {
        let mut p=|s: SolveProgress| if s.stage==stage {ControlFlow::Break(())}else{ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        assert!(matches!(s.evaluate_stress_with_prescribed_motion(&rho,&loads,motions,Default::default(),StressMeasure3::NormalizedAverage,&mut c),
            Err(StressError3::Evaluation(EvaluationStop::Cancelled))));
        assert_eq!(s.operator.scales(),before);
    }
    let bad=|_,_|[f64::NAN;3];let motions: &StressMotions3<'_>=&[Some(&bad)];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    assert!(s.evaluate_stress_with_prescribed_motion(&rho,&loads,motions,Default::default(),StressMeasure3::NormalizedAverage,&mut c).is_err());
    assert_eq!(s.operator.scales(),before);
}

#[test]
fn g0_motion_family_and_support_admission_precedes_callbacks_and_solves() {
    let mut p=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut p).unwrap();
    let op=AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(), &tree(false), &Slab,
        &IsotropicElastic::new(1.0,0.0,1.0).unwrap(), &|p|p[1]==0.0, ElasticityOptions3::default(), &mut q).unwrap();
    let f=vec![0.0;op.n()];let loads=[LoadCase{force:&f,weight:1.0}];let mut s=CutDensityStudy3::new(op,0.15,params());
    let rho=vec![0.5;s.cells()];let calls=Cell::new(0);
    let law=|p,n|{calls.set(calls.get()+1);motion(p,n)};let motions: &StressMotions3<'_>=&[Some(&law)];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    assert!(s.evaluate_stress_with_prescribed_motion(&rho,&loads,motions,Default::default(),StressMeasure3::NormalizedAverage,&mut c).is_err());
    assert!(s.evaluate_stress_with_prescribed_motion(&rho,&loads,&[],Default::default(),StressMeasure3::NormalizedAverage,&mut c).is_err());
    assert_eq!(calls.get(),0);assert_eq!(c.work().linear_solves,0);
}

#[test]
fn g1_motion_stress_design_reduces_volume_and_retains_matching_physics() {
    use crate::sdf3::design::{StressDesignOptions3, StressDesignStudy3};
    let op=build(true);let f=op.elasticity().body_load(&|_|[0.0,0.0,-0.001],||ControlFlow::Continue(())).unwrap();
    let loads=[LoadCase{force:&f,weight:1.0}];let motions: &StressMotions3<'_>=&[Some(&motion)];
    let mut s=CutDensityStudy3::new(op,0.15,params());let rho=vec![0.65;s.cells()];
    let options=StressDesignOptions3 { stress_limit:1.0,density_floor:0.1,..Default::default() };
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let mut session=StressDesignStudy3::new_with_prescribed_motion(&mut s,&loads,motions,&rho,options,StressMeasure3::SampledPeakBound,&mut c).unwrap();
    let initial=session.accepted().volume_fraction;
    session.run(8).unwrap();let best=session.best_feasible().unwrap();
    assert!(best.volume_fraction<initial-1e-4);assert!(best.aggregate<=options.stress_limit*(1.0+options.optimizer.tolerance));
    let endpoint=session.accepted().clone();let scales=session.study().operator().scales().to_vec();drop(session);
    let again=s.evaluate_stress_with_prescribed_motion(&endpoint.rho,&loads,motions,options.stress,StressMeasure3::SampledPeakBound,&mut c).unwrap();
    same(&endpoint,&again);assert_eq!(s.operator.scales(),scales);
}

#[test]
fn g5_motion_stress_checkpoint_restores_the_same_search_not_frozen_liftings() {
    use crate::sdf3::design::{StressDesignOptions3, StressDesignStudy3};
    let mut a=CutDensityStudy3::new(build(false),0.15,params());
    let mut b=CutDensityStudy3::new(build(false),0.15,params());
    let f=a.operator.elasticity().body_load(&|_|[0.0,0.0,-0.001],||ControlFlow::Continue(())).unwrap();
    let loads=[LoadCase{force:&f,weight:1.0}];let motions: &StressMotions3<'_>=&[Some(&motion)];
    let rho=vec![0.65;a.cells()];let options=StressDesignOptions3{stress_limit:1.0,density_floor:0.1,..Default::default()};
    let mut pa=|_|ControlFlow::Continue(());let mut ca=SolveControl::new(SolveBudget::default(),&mut pa);
    let mut pb=|_|ControlFlow::Continue(());let mut cb=SolveControl::new(SolveBudget::default(),&mut pb);
    let mut a=StressDesignStudy3::new_with_prescribed_motion(&mut a,&loads,motions,&rho,options,StressMeasure3::SampledPeakBound,&mut ca).unwrap();
    let mut split=StressDesignStudy3::new_with_prescribed_motion(&mut b,&loads,motions,&rho,options,StressMeasure3::SampledPeakBound,&mut cb).unwrap();
    a.run(4).unwrap();split.run(2).unwrap();let checkpoint=split.checkpoint();let cost=checkpoint.restoration_cost();drop(split);
    let mut resumed=StressDesignStudy3::restore_with_prescribed_motion(&mut b,&loads,motions,checkpoint,options,StressMeasure3::SampledPeakBound,&mut cb).unwrap();
    resumed.run(2).unwrap();same(a.accepted(),resumed.accepted());
    assert_eq!(a.history(),resumed.history());same(a.best_feasible().unwrap(),resumed.best_feasible().unwrap());
    let mut expected=a.checkpoint().optimizer;expected.work.evaluations+=cost;
    assert_eq!(expected,resumed.checkpoint().optimizer);
}

#[test]
fn g4_changed_motion_cannot_restore_a_saved_stress_state() {
    use crate::sdf3::design::{StressDesignOptions3, StressDesignStudy3};
    let mut s=CutDensityStudy3::new(build(false),0.15,params());let f=vec![0.0;s.operator.n()];
    let loads=[LoadCase{force:&f,weight:1.0}];let motions: &StressMotions3<'_>=&[Some(&motion)];
    let rho=vec![0.6;s.cells()];let options=StressDesignOptions3{stress_limit:1.0,density_floor:0.1,..Default::default()};
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let session=StressDesignStudy3::new_with_prescribed_motion(&mut s,&loads,motions,&rho,options,StressMeasure3::SampledPeakBound,&mut c).unwrap();
    let cp=session.checkpoint();let scales=session.study().operator().scales().to_vec();drop(session);
    let changed: &StressMotions3<'_>=&[Some(&opposite)];
    assert!(StressDesignStudy3::restore_with_prescribed_motion(&mut s,&loads,changed,cp,options,StressMeasure3::SampledPeakBound,&mut c).is_err());
    assert_eq!(s.operator.scales(),scales);
}
