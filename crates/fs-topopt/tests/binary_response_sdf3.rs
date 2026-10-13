//! Real binary stiffness, unchanged reference loads, independent verification.
use std::{cell::Cell, ops::ControlFlow};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3};
use fs_topopt::sdf3::response::ResponseError3;
use fs_topopt::sdf3::response::limits::*;
use fs_topopt::sdf3::response::limits::binary::*;
use fs_solver::CgState;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { (p[0]-0.17)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x-Interval::new(0.17, 0.17))*(x-Interval::new(0.83, 0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], a: HeightAxis) -> Interval {
        if a == HeightAxis::X { Interval::new(2.0, 2.0)*Interval::new(lo[0], hi[0])-Interval::new(1.0, 1.0) }
        else { Interval::new(0.0, 0.0) }
    }
}
fn motion(p: [f64; 3], _: [f64; 3]) -> [f64; 3] { [0.02*(p[0]-0.17), 0.0, 0.0] }
fn right(_: [f64; 3], n: [f64; 3]) -> [f64; 3] { if n[0] > 0.0 { [1.0, 0.0, 0.0] } else { [0.0; 3] } }
fn fixture(mixed: bool) -> (CutDensityStudy3<AdaptiveSolveSpace3>, Vec<f64>, Vec<f64>) {
    let op = endpoint_operator(mixed);
    let force = op.elasticity().body_load(&|_| [0.002, 0.0, -0.001], || ControlFlow::Continue(())).unwrap();
    let observe = op.elasticity().body_load(&|p| [1.0+p[1], 0.0, 0.0], || ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(op, 0.15, SimpParams { penal: 3.0, beta: 2.0, ..Default::default() }), force, observe)
}
fn limits(q: &[f64]) -> [ResponseLimit3<'_>; 3] {
    [ResponseLimit3 { response: BoundedResponse3::Displacement(q), sense: ResponseLimitSense3::AtMost, bound: 0.006, scale: 0.02 },
     ResponseLimit3 { response: BoundedResponse3::Reaction(&right), sense: ResponseLimitSense3::AtLeast, bound: 0.001, scale: 0.002 },
     ResponseLimit3 { response: BoundedResponse3::Reaction(&right), sense: ResponseLimitSense3::AtMost, bound: 0.004, scale: 0.002 }]
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol*a.abs().max(b.abs()).max(1e-9), "{a} vs {b}"); }

fn endpoint_operator(mixed: bool) -> AdaptiveSolveSpace3 {
    let t = Octree3::uniform(1, 4, 4096).unwrap();
    let t = if mixed { t.refined(&[*t.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap() } else { t };
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut p).unwrap();
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &t, &Slab, &IsotropicElastic::new(1.0, 0.0, 1.0).unwrap(), &|_| false, &|_, _| true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut q).unwrap();
    AdaptiveSolveSpace3::jacobi(op, 100_000_000)
}
fn compare_fields(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().map(|v| v.abs()).fold(1e-30_f64, f64::max);
    assert!(a.iter().zip(b).all(|(a,b)| (a-b).abs() < 1e-7*scale));
}

#[test]
fn g0_binary_endpoint_matches_direct_material_solves_on_uniform_and_hanging_grids() {
    for mixed in [false, true] {
        let (mut s, f, q) = fixture(mixed); let rows = limits(&q);
        let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
        let rho: Vec<_> = (0..s.cells()).map(|i| 0.25+0.055*(i%8) as f64).collect();
        let initial = s.operator().elasticity().scales().to_vec();
        let nodes = s.operator().elasticity().nodes().to_vec();
        let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        let gray = s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c).unwrap();
        let options = BinaryEndpointOptions3 { threshold: 0.45, void_stiffness: 0.01, ..Default::default() };
        let before = c.work();
        let e = s.evaluate_binary_response_limits(&rho, &cases, options, &mut c).unwrap();
        assert_eq!(e.solid, gray.projected_rho.iter().map(|r| *r>=options.threshold).collect::<Vec<_>>());
        assert!(e.solid.iter().any(|s| *s) && e.solid.iter().any(|s| !*s));
        assert_eq!(c.work().linear_solves-before.linear_solves, 2, "one filter plus one primal, no adjoint");
        assert_eq!(s.operator().elasticity().scales(), initial);
        assert_eq!(s.operator().elasticity().nodes(), nodes);
        assert_eq!(e.source_volume_fraction, gray.volume_fraction);
        let mut op = endpoint_operator(mixed); op.set_scales(&e.stiffness_scales()).unwrap();
        let lifting = op.elasticity().prescribed_displacement_load(&motion, || ControlFlow::Continue(())).unwrap();
        let rhs: Vec<_> = f.iter().zip(lifting).map(|(f,b)| f+b).collect();
        let prepared = op.prepare(|_| ControlFlow::Continue(())).unwrap();
        let mut solve = CgState::new(&op, &prepared, &rhs);
        assert!(solve.run(&op, &prepared, 1e-12, 10_000).converged);
        assert!(op.elasticity().field_residual(&solve.x, &rhs, || ControlFlow::Continue(())).unwrap()<1e-8);
        compare_fields(&e.displacements[0], &solve.x);
        let displacement: f64 = q.iter().zip(&solve.x).map(|(q,u)| q*u).sum();
        close(e.measurements[0][0], displacement, 1e-7);
        let reaction = op.elasticity().embedded_reaction(&solve.x, Some(&motion), &right, || ControlFlow::Continue(())).unwrap().value;
        close(e.measurements[0][1], reaction, 1e-7);
        assert_eq!(e.measurements[0][1], e.measurements[0][2]);
        let volumes = op.elasticity().volumes(); let total: f64 = volumes.iter().sum();
        let material: f64 = volumes.iter().zip(&e.solid).filter(|(_,s)| **s).map(|(v,_)| v/total).sum();
        assert_eq!(e.volume_fraction, material);
    }
}

#[test]
fn g1_a_feasible_gray_design_is_not_evidence_of_a_feasible_binary_endpoint() {
    let (mut s, f, _) = fixture(false); let rho = vec![0.5; s.cells()];
    let mut rows = [ResponseLimit3 { response: BoundedResponse3::Reaction(&right),
        sense: ResponseLimitSense3::AtMost, bound: 1.0, scale: 0.002 }];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let options = BinaryEndpointOptions3 { threshold: 0.25, ..Default::default() };
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let gray = s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c).unwrap();
    let solid = s.evaluate_binary_response_limits(&rho, &cases, options, &mut c).unwrap();
    let (a,b) = (gray.measurements[0][0], solid.measurements[0][0]); assert!((a-b).abs()>1e-9);
    rows[0].bound = f64::midpoint(a,b);
    rows[0].sense = if b>a { ResponseLimitSense3::AtMost } else { ResponseLimitSense3::AtLeast };
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    assert!(s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c).unwrap().feasible(0.0));
    let endpoint = s.evaluate_binary_response_limits(&rho, &cases, options, &mut c).unwrap();
    assert!(!endpoint.feasible(0.0)); assert!(endpoint.worst_residual>0.0);
}

#[test]
fn g0_ties_and_protected_regions_are_physical_not_refiltered_raw_controls() {
    let (s, f, q) = fixture(false); let mut labels = vec![PhysicalRegion3::Design; s.cells()];
    labels[0] = PhysicalRegion3::Solid; labels[1] = PhysicalRegion3::Void;
    let mut s = s.with_physical_regions(labels).unwrap(); let rho = vec![0.5; s.cells()];
    let rows = limits(&q); let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let gray = s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c).unwrap();
    let options = BinaryEndpointOptions3 { threshold: gray.projected_rho[2], void_stiffness: 0.02, ..Default::default() };
    let a = s.evaluate_binary_response_limits(&rho, &cases, options, &mut c).unwrap();
    assert!(a.solid[0] && !a.solid[1] && a.solid[2]);
    for (i,r) in gray.projected_rho.iter().enumerate().skip(2) { assert_eq!(a.solid[i], *r>=options.threshold); }
    assert_eq!(a.stiffness_scales()[0], 1.0); assert_eq!(a.stiffness_scales()[1], 0.02);
    let b = s.evaluate_binary_response_limits(&rho, &cases, options, &mut c).unwrap();
    assert_eq!(a.solid,b.solid); assert_eq!(a.measurements,b.measurements); assert_eq!(a.displacements,b.displacements);
}

#[test]
fn g4_failed_binary_preparation_solve_or_publication_preserves_incoming_material() {
    for stage in ["sdf3-preconditioner-ready", "response-limits-primal", "response-binary-publish"] {
        let (mut s, f, q) = fixture(false); let rho = vec![0.5; s.cells()]; let rows = limits(&q);
        let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
        let before = s.operator().elasticity().scales().to_vec();
        let mut p = |s: SolveProgress| if s.stage==stage && (stage!="response-limits-primal"||s.solve_iterations>0) {
            ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        assert!(matches!(s.evaluate_binary_response_limits(&rho,&cases,Default::default(),&mut c),
            Err(ResponseError3::Evaluation(EvaluationStop::Cancelled))));
        assert_eq!(s.operator().elasticity().scales(),before);
        let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut p);
        assert!(s.evaluate_binary_response_limits(&rho,&cases,Default::default(),&mut c).is_err());
        assert_eq!(c.work().linear_iterations,1); assert_eq!(s.operator().elasticity().scales(),before);
    }
}

#[test]
fn g0_endpoint_policy_and_whole_limit_family_are_admitted_before_callbacks() {
    let (mut s, f, q) = fixture(false); let rho = vec![0.5;s.cells()]; let calls=Cell::new(0);
    let mode=|p,n| { calls.set(calls.get()+1); right(p,n) };
    let mut rows=limits(&q); rows[1].response=BoundedResponse3::Reaction(&mode);
    let mut p=|_|ControlFlow::Continue(()); let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    for (threshold,void_stiffness) in [(0.0,1e-3),(1.1,1e-3),(f64::NAN,1e-3),(0.5,0.0),(0.5,1.0),(0.5,f64::NAN)] {
        let cases=[ResponseLimitCase3 { force:&f,prescribed:Some(&motion),limits:&rows }];
        assert!(s.evaluate_binary_response_limits(&rho,&cases,BinaryEndpointOptions3 { threshold,void_stiffness,..Default::default() },&mut c).is_err());
    }
    rows[2].bound=f64::NAN;
    assert!(s.evaluate_binary_response_limits(&rho,&[ResponseLimitCase3 { force:&f,prescribed:Some(&motion),limits:&rows }],Default::default(),&mut c).is_err());
    assert_eq!(calls.get(),0); assert_eq!(c.work().linear_solves,0);
}

#[path = "binary_response_sdf3/scan.rs"]
mod scan;
