//! G0/G1/G3/G4/G5: binary response fields and actual material stresses.
use super::*;
use std::cell::Cell;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{SimpParams, SolveBudget, SolveProgress};

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { (p[0]-0.17)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0],hi[0]);
        (x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], a: HeightAxis) -> Interval {
        if a == HeightAxis::X { Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0) }
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
    [ResponseLimit3 { response:BoundedResponse3::Reaction(&right),sense:ResponseLimitSense3::AtLeast,bound:-10.0,scale:1.0 },
     ResponseLimit3 { response:BoundedResponse3::Displacement(q),sense:ResponseLimitSense3::AtMost,bound:10.0,scale:1.0 }]
}
fn raw(s: &Study) -> Vec<f64> { (0..s.cells()).map(|i|0.25+0.5*i as f64/(s.cells()-1) as f64).collect() }
fn options() -> BinaryResponseStressOptions3 {
    BinaryResponseStressOptions3 { endpoint:BinaryEndpointOptions3{void_stiffness:0.05,..Default::default()},..Default::default() }
}
fn same(a: &BinaryResponseStressEvaluation3, b: &BinaryResponseStressEvaluation3) {
    assert_eq!(a.response.solid,b.response.solid); assert_eq!(a.response.displacements,b.response.displacements);
    assert_eq!(a.response.measurements,b.response.measurements); assert_eq!(a.cases,b.cases);
    assert_eq!(a.stress_aggregate,b.stress_aggregate); assert_eq!(a.constraint,b.constraint);
}
// Independent scalar formula at the small finite stresses used by the fixtures.
fn direct_vm(s: [f64;6]) -> f64 {
    (0.5*((s[0]-s[1]).powi(2)+(s[1]-s[2]).powi(2)+(s[2]-s[0]).powi(2))
        +3.0*(s[3]*s[3]+s[4]*s[4]+s[5]*s[5])).sqrt()
}

#[test]
fn g3_binary_stress_matches_direct_physics_with_no_additional_solves() {
    for mixed in [false,true] {
        let (mut s,f,q)=fixture(mixed); let a=limits(&q); let neg:Vec<_>=f.iter().map(|v|-v).collect();
        let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a},
            ResponseLimitCase3{force:&neg,prescribed:Some(&opposite),limits:&a}];
        let rho=raw(&s); let previous=s.operator.scales().to_vec(); let mut preparations=0;
        let mut p=|p:SolveProgress| {if p.stage=="sdf3-preconditioner-start" {preparations+=1;} ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        let actual=s.evaluate_binary_response_stress_limits(&rho,&cases,options(),&mut c).unwrap();
        let spent=c.work(); drop(c); assert_eq!(preparations,1); assert_eq!(spent.linear_solves,1+cases.len());
        assert_eq!(s.operator.scales(),previous); assert!(actual.response.solid.contains(&true));
        assert!(actual.response.solid.contains(&false));
        let mut p=|_|ControlFlow::Continue(()); let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        let response=s.evaluate_binary_response_limits(&rho,&cases,options().endpoint,&mut c).unwrap();
        assert_eq!(actual.response.displacements,response.displacements); assert_eq!(spent,c.work());
        s.operator.set_scales(&response.stiffness_scales()).unwrap();
        let mut total=0usize; let mut sum=0.0; let mut maximum=0.0_f64;
        for (u,case) in response.displacements.iter().zip(&actual.cases) {
            let points=s.operator.bulk_stress(u,options().max_points,||ControlFlow::Continue(())).unwrap();
            assert_eq!(case.points,points.len()); total+=points.len();
            let mut maxima=[0.0_f64;2];
            for point in points {
                let vm=direct_vm(point.stress); maximum=maximum.max(vm); sum+=vm.powf(options().aggregation_power);
                let phase=usize::from(response.solid[point.cell]); maxima[phase]=maxima[phase].max(vm);
            }
            assert!((case.void_max-maxima[0]).abs()<1e-12);
            assert!((case.solid_max-maxima[1]).abs()<1e-12);
        }
        s.operator.set_scales(&previous).unwrap();
        let expected=sum.powf(1.0/options().aggregation_power);
        assert!((actual.stress_aggregate-expected).abs()<1e-11*expected);
        assert!((actual.sampled_physical_max-maximum).abs()<1e-12); assert_eq!(actual.point_count,total);
        assert!(actual.stress_aggregate>=actual.sampled_physical_max);
    }
}

#[test]
fn g0_joint_acceptance_cannot_hide_a_failed_stress_or_response_bound() {
    let (mut s,f,q)=fixture(false); let mut a=limits(&q); let rho=raw(&s);
    let mut p=|_|ControlFlow::Continue(()); let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let first=s.evaluate_binary_response_stress_limits(&rho,&cases,options(),&mut c).unwrap();
    let strict=BinaryResponseStressOptions3{stress_limit:0.5*first.stress_aggregate,..options()};
    let rejected=s.evaluate_binary_response_stress_limits(&rho,&cases,strict,&mut c).unwrap();
    assert!(rejected.response.feasible(0.0)); assert!(!rejected.feasible(0.0)); assert!(rejected.stress_residual>0.99);
    a[0].bound=first.response.measurements[0][0]+1.0;
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let loose=BinaryResponseStressOptions3{stress_limit:100.0*first.stress_aggregate,..options()};
    let rejected=s.evaluate_binary_response_stress_limits(&rho,&cases,loose,&mut c).unwrap();
    assert!(rejected.stress_residual<0.0); assert!(!rejected.feasible(0.0)); assert!(rejected.worst_residual>0.99);
}

#[test]
fn g1_all_weak_cells_still_carry_physical_stress_under_imposed_extension() {
    let (mut s,_,q)=fixture(false); let f=vec![0.0;s.operator.n()]; let a=limits(&q); let rho=vec![0.4;s.cells()];
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let o=BinaryResponseStressOptions3{endpoint:BinaryEndpointOptions3{threshold:1.0,..options().endpoint},..options()};
    let mut p=|_|ControlFlow::Continue(()); let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let e=s.evaluate_binary_response_stress_limits(&rho,&cases,o,&mut c).unwrap();
    assert!(e.response.solid.iter().all(|s|!*s)); assert_eq!(e.response.volume_fraction,0.0);
    assert_eq!(e.cases[0].solid_max,0.0); assert!(e.cases[0].void_max>0.0);
    let expected=o.endpoint.void_stiffness*0.02;
    assert!((e.sampled_physical_max-expected).abs()<1e-8*expected);
    let strict=BinaryResponseStressOptions3{stress_limit:0.5*e.sampled_physical_max,..o};
    assert!(!s.evaluate_binary_response_stress_limits(&rho,&cases,strict,&mut c).unwrap().feasible(0.0));
}

#[test]
fn g4_sampling_stops_preserve_material_and_charge_completed_equilibria() {
    let (mut s,f,q)=fixture(true); let a=limits(&q); let rho=raw(&s);
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let previous=s.operator.scales().to_vec();
    for stage in ["response-binary-stress-sample","response-binary-stress-measure","response-binary-stress-publish"] {
        let mut p=|p:SolveProgress| if p.stage==stage {ControlFlow::Break(())}else{ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        assert!(matches!(s.evaluate_binary_response_stress_limits(&rho,&cases,options(),&mut c),
            Err(ResponseStressError3::Evaluation(EvaluationStop::Cancelled))));
        assert!(c.work().linear_iterations>0); assert_eq!(s.operator.scales(),previous);
    }
    let mut p=|_|ControlFlow::Continue(()); let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let e=s.evaluate_binary_response_stress_limits(&rho,&cases,options(),&mut c).unwrap();
    let limited=BinaryResponseStressOptions3{max_points:e.point_count,..options()};
    let both=[cases[0],cases[0]];
    assert!(matches!(s.evaluate_binary_response_stress_limits(&rho,&both,limited,&mut c),
        Err(ResponseStressError3::Stress(StressError3::PointBudget))));
    assert_eq!(s.operator.scales(),previous);
    same(&e,&s.evaluate_binary_response_stress_limits(&rho,&cases,options(),&mut c).unwrap());
}

#[test]
fn g3_protected_regions_and_changed_thresholds_recompute_actual_stress() {
    let (s,f,q)=fixture(true); let rho=raw(&s); let a=limits(&q);
    let mut regions=vec![PhysicalRegion3::Design;s.cells()]; regions[0]=PhysicalRegion3::Solid; regions[1]=PhysicalRegion3::Void;
    let mut s=s.with_physical_regions(regions).unwrap(); let previous=s.operator.scales().to_vec();
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let mut p=|_|ControlFlow::Continue(()); let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let low=BinaryResponseStressOptions3{endpoint:BinaryEndpointOptions3{threshold:0.1,..options().endpoint},..options()};
    let high=BinaryResponseStressOptions3{endpoint:BinaryEndpointOptions3{threshold:0.9,..options().endpoint},..options()};
    let first=s.evaluate_binary_response_stress_limits(&rho,&cases,low,&mut c).unwrap();
    let changed=s.evaluate_binary_response_stress_limits(&rho,&cases,high,&mut c).unwrap();
    for e in [&first,&changed] {assert!(e.response.solid[0]); assert!(!e.response.solid[1]);}
    assert_ne!(first.response.solid,changed.response.solid); assert_ne!(first.stress_aggregate,changed.stress_aggregate);
    same(&first,&s.evaluate_binary_response_stress_limits(&rho,&cases,low,&mut c).unwrap()); assert_eq!(s.operator.scales(),previous);
}

#[test]
fn g0_invalid_stress_policy_is_refused_before_callbacks_or_solves() {
    let (mut s,f,q)=fixture(false); let rho=raw(&s); let a=limits(&q); let calls=Cell::new(0);
    let g=|p,n| {calls.set(calls.get()+1); motion(p,n)};
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&g),limits:&a}];
    let mut p=|_|ControlFlow::Continue(()); let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    for o in [BinaryResponseStressOptions3{stress_limit:f64::NAN,..options()},
        BinaryResponseStressOptions3{max_points:0,..options()},
        BinaryResponseStressOptions3{aggregation_power:1.0,..options()}] {
        assert!(s.evaluate_binary_response_stress_limits(&rho,&cases,o,&mut c).is_err());
    }
    assert_eq!(calls.get(),0); assert_eq!(c.work().linear_solves,0);
}

#[test]
fn g0_scaled_physical_norm_keeps_extreme_finite_samples_and_zero_states() {
    let mut n=StressNorm::default(); assert_eq!(n.value(8.0).unwrap(),0.0);
    n.add(3.0,2.0).unwrap(); n.add(4.0,2.0).unwrap(); assert!((n.value(2.0).unwrap()-5.0).abs()<1e-14);
    for p in [2.0,8.0,64.0] {
        let mut n=StressNorm::default();
        for v in [0.0,1e-200,1e200,1e200] {n.add(v,p).unwrap();}
        let expected=1e200_f64*2.0_f64.powf(1.0/p);
        assert!((n.value(p).unwrap()-expected).abs()<1e-14*expected);
        assert!(n.add(f64::NAN,p).is_err());
    }
}
