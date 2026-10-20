//! Real retained CutFEM equations and full filter/projection/material gradients.
use super::*;
use std::{cell::Cell, ops::ControlFlow};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{SimpParams, SolveBudget, SolveProgress};
use crate::sdf3::PhysicalRegion3;

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { (p[0]-0.17)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        if axis == HeightAxis::X { Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0) }
        else { Interval::new(0.0,0.0) }
    }
}
fn motion(p: [f64; 3], _: [f64; 3]) -> [f64; 3] { [0.02*(p[0]-0.17),0.0,0.0] }
fn opposite(p: [f64; 3], n: [f64; 3]) -> [f64; 3] { motion(p,n).map(|v| -v) }
fn right(_: [f64; 3], n: [f64; 3]) -> [f64; 3] { if n[0]>0.0 { [1.0,0.0,0.0] } else { [0.0;3] } }
fn fixture(mixed: bool) -> (Study, Vec<f64>, Vec<f64>) {
    let tree = Octree3::uniform(1,4,4096).unwrap();
    let tree = if mixed { tree.refined(&[*tree.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap() } else { tree };
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut p).unwrap();
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),
        &tree, &Slab, &IsotropicElastic::new(1.0,0.0,1.0).unwrap(), &|_|false, &|_,_|true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut q).unwrap();
    let f = op.body_load(&|p|[0.001*(1.0+p[1]),0.0,-0.0007], || ControlFlow::Continue(())).unwrap();
    let q = op.body_load(&|p|[1.0+p[1],0.0,0.0], || ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,
        SimpParams { beta:2.0, penal:3.0, e_min:1e-4, ..Default::default() }),f,q)
}
fn limits(q: &[f64]) -> [ResponseLimit3<'_>;2] {
    [ResponseLimit3 { response:BoundedResponse3::Reaction(&right),sense:ResponseLimitSense3::AtLeast,bound:0.001,scale:0.002 },
     ResponseLimit3 { response:BoundedResponse3::Displacement(q),sense:ResponseLimitSense3::AtMost,bound:0.01,scale:0.005 }]
}
fn raw(s: &Study) -> Vec<f64> { (0..s.cells()).map(|i|0.45+0.015*(i%7) as f64).collect() }
fn same(a: &ResponseStressEvaluation3, b: &ResponseStressEvaluation3) {
    assert_eq!(a.constraint,b.constraint); assert_eq!(a.constraint_gradient,b.constraint_gradient);
    assert_eq!(a.response.rho,b.response.rho); assert_eq!(a.response.displacements,b.response.displacements);
    assert_eq!(a.response.measurements,b.response.measurements); assert_eq!(a.stress.aggregate,b.stress.aggregate);
    assert_eq!(a.stress.displacements,b.stress.displacements); assert_eq!(a.stress.gradient,b.stress.gradient);
}

#[test]
fn g3_composition_matches_standalone_physics_and_charges_both_preparations() {
    let (mut s,f,q)=fixture(true);let rho=raw(&s);let limits=limits(&q);
    let cases=[ResponseLimitCase3 { force:&f,prescribed:Some(&motion),limits:&limits }];
    let options=ResponseStressOptions3::default();let before=s.operator.scales().to_vec();
    let mut preparations=0;
    let mut p=|s:SolveProgress| { if s.stage=="sdf3-preconditioner-start" {preparations+=1;} ControlFlow::Continue(()) };
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let joint=s.evaluate_response_stress_limits(&rho,&cases,options,&mut c).unwrap();
    let spent=c.work();drop(c);assert_eq!(preparations,2);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let r=s.evaluate_response_limits(&rho,&cases,options.response,&mut c).unwrap();
    let t=s.evaluate_stress_with_prescribed_motion(&rho,&[LoadCase{force:&f,weight:1.0}],&[Some(&motion)],
        options.stress,StressMeasure3::SampledPeakBound,&mut c).unwrap();
    assert_eq!(joint.response.displacements,r.displacements);assert_eq!(joint.stress.displacements,t.displacements);
    assert_eq!(joint.response.constraint_gradient,r.constraint_gradient);assert_eq!(joint.stress.gradient,t.gradient);
    let mut residuals:Vec<_>=r.residuals.iter().flatten().copied().collect();residuals.push(t.aggregate/options.stress_limit-1.0);
    let direct=aggregate(&residuals,options.response.sharpness).unwrap().0;
    assert!((joint.constraint-direct).abs()<1e-12*direct.abs().max(1.0));
    assert_eq!(spent,c.work());assert_eq!(s.operator.scales(),before);
}

#[test]
fn g3_mixed_force_motion_joint_gradient_matches_independent_resolves() {
    for mixed in [false,true] {
        let (mut s,f,q)=fixture(mixed);let neg:Vec<_>=f.iter().map(|v|-v).collect();let rho=raw(&s);
        let a=limits(&q);let b=[ResponseLimit3{sense:ResponseLimitSense3::AtMost,bound:-0.001,..a[0]},
            ResponseLimit3{sense:ResponseLimitSense3::AtLeast,bound:-0.01,..a[1]}];
        let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a},
            ResponseLimitCase3{force:&neg,prescribed:Some(&opposite),limits:&b}];
        let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        let seed=s.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c).unwrap();
        // Arrange two active blocks: neither component tangent is numerically suppressed.
        let a=[ResponseLimit3{bound:seed.response.measurements[0][0]+0.04*a[0].scale,..a[0]},
            ResponseLimit3{bound:seed.response.measurements[0][1]+0.1*a[1].scale,..a[1]}];
        let b=[ResponseLimit3{bound:seed.response.measurements[1][0]-0.01*b[0].scale,..b[0]},
            ResponseLimit3{bound:seed.response.measurements[1][1]-0.2*b[1].scale,..b[1]}];
        let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a},
            ResponseLimitCase3{force:&neg,prescribed:Some(&opposite),limits:&b}];
        let options=ResponseStressOptions3 { stress_limit:seed.stress.aggregate/1.03,
            response:ResponseLimitOptions3{sharpness:12.0,..Default::default()},..Default::default() };
        let base=s.evaluate_response_stress_limits(&rho,&cases,options,&mut c).unwrap();
        assert!(base.component_weights.iter().all(|w|*w>0.05));
        for i in [0,s.cells()/2,s.cells()-1] {
            let mut a=rho.clone();let mut b=rho.clone();a[i]+=1e-4;b[i]-=1e-4;
            let a=s.evaluate_response_stress_limits(&a,&cases,options,&mut c).unwrap();
            let b=s.evaluate_response_stress_limits(&b,&cases,options,&mut c).unwrap();
            let fd=(a.constraint-b.constraint)/2e-4;let g=base.constraint_gradient[i];
            assert!((fd-g).abs()<4e-4*fd.abs().max(g.abs()).max(1e-8),"cell={i}: {fd} != {g}");
        }
    }
}

#[test]
fn g0_a_failed_stress_or_response_limit_cannot_be_hidden_by_the_other_block() {
    let (mut s,f,q)=fixture(false);let rho=raw(&s);let mut a=limits(&q);
    a[0].bound=-1.0;a[1].bound=1.0;
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let seed=s.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c).unwrap();
    assert!(seed.feasible(1e-6));
    let stress=s.evaluate_response_stress_limits(&rho,&cases,ResponseStressOptions3{
        stress_limit:seed.stress.aggregate/2.0,..Default::default()},&mut c).unwrap();
    assert!(stress.response.feasible(1e-6));assert!(!stress.feasible(1e-6));assert!(stress.constraint>0.99);
    a[0].bound=seed.response.measurements[0][0]+a[0].scale;
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let response=s.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c).unwrap();
    assert!(response.stress_residual<0.0);assert!(!response.feasible(1e-6));assert!(response.worst_residual>0.99);
}

#[test]
fn g4_second_phase_failure_cancellation_and_budget_keep_the_original_material() {
    let (mut s,f,q)=fixture(true);let rho=raw(&s);let a=limits(&q);
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let before=s.operator.scales().to_vec();
    for stage in ["response-stress-between","sdf3-stress-motion-pullback","response-stress-publish"] {
        let mut p=|s:SolveProgress| if s.stage==stage {ControlFlow::Break(())} else {ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        assert!(matches!(s.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c),
            Err(ResponseStressError3::Evaluation(EvaluationStop::Cancelled))));
        assert!(c.work().linear_iterations>0);assert_eq!(s.operator.scales(),before);
    }
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let options=ResponseStressOptions3{stress:StressOptions3{max_points:1,..Default::default()},..Default::default()};
    assert!(matches!(s.evaluate_response_stress_limits(&rho,&cases,options,&mut c),Err(ResponseStressError3::Stress(StressError3::PointBudget))));
    assert!(c.work().linear_iterations>0);assert_eq!(s.operator.scales(),before);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget{total_iterations:1,..Default::default()},&mut p);
    assert!(s.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c).is_err());
    assert_eq!(c.work().linear_iterations,1);assert_eq!(s.operator.scales(),before);
}

#[test]
fn g0_full_policy_admission_precedes_callbacks_and_physics() {
    let (mut s,f,q)=fixture(false);let rho=raw(&s);let a=limits(&q);let calls=Cell::new(0);
    let law=|p,n| {calls.set(calls.get()+1);motion(p,n)};
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&law),limits:&a}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    for options in [ResponseStressOptions3{stress_limit:f64::NAN,..Default::default()},
        ResponseStressOptions3{stress:StressOptions3{max_points:0,..Default::default()},..Default::default()},
        ResponseStressOptions3{response:ResponseLimitOptions3{sharpness:0.0,..Default::default()},..Default::default()}] {
        assert!(s.evaluate_response_stress_limits(&rho,&cases,options,&mut c).is_err());
    }
    assert_eq!(calls.get(),0);assert_eq!(c.work().linear_solves,0);
}

#[test]
fn g3_both_components_honor_the_same_protected_material_regions() {
    let (s,f,q)=fixture(true);let rho=raw(&s);let a=limits(&q);
    let mut labels=vec![PhysicalRegion3::Design;s.cells()];labels[0]=PhysicalRegion3::Solid;labels[1]=PhysicalRegion3::Void;
    let mut s=s.with_physical_regions(labels).unwrap();
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let a=s.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c).unwrap();
    assert_eq!(&a.stress.projected_rho[..2],&[1.0,0.0]);
    assert_eq!(a.response.projected_rho,a.stress.projected_rho);
    let b=s.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c).unwrap();same(&a,&b);
}
