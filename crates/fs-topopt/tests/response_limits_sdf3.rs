//! Actual constrained cut-cell physics, not surrogate limit values.
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
use fs_topopt::sdf3::response::{ReactionTarget3, ResponseCase3, ResponseTarget3, ResponseError3};
use fs_topopt::sdf3::response::limits::*;

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;
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
fn fixture(mixed: bool) -> (Study, Vec<f64>, Vec<f64>) {
    let t = Octree3::uniform(1, 4, 4096).unwrap();
    let t = if mixed { t.refined(&[*t.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap() } else { t };
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut p).unwrap();
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &t, &Slab, &IsotropicElastic::new(1.0, 0.0, 1.0).unwrap(), &|_| false, &|_, _| true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut q).unwrap();
    let force = op.body_load(&|_| [0.002, 0.0, -0.001], || ControlFlow::Continue(())).unwrap();
    let observe = op.body_load(&|p| [1.0+p[1], 0.0, 0.0], || ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000), 0.15,
        SimpParams { penal: 3.0, beta: 2.0, ..Default::default() }), force, observe)
}
fn limits(q: &[f64]) -> [ResponseLimit3<'_>; 3] {
    [ResponseLimit3 { response: BoundedResponse3::Displacement(q), sense: ResponseLimitSense3::AtMost, bound: 0.006, scale: 0.02 },
     ResponseLimit3 { response: BoundedResponse3::Reaction(&right), sense: ResponseLimitSense3::AtLeast, bound: 0.001, scale: 0.002 },
     ResponseLimit3 { response: BoundedResponse3::Reaction(&right), sense: ResponseLimitSense3::AtMost, bound: 0.004, scale: 0.002 }]
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol*a.abs().max(b.abs()).max(1e-9), "{a} vs {b}"); }

#[test]
fn g3_full_constraint_and_volume_gradients_match_coordinate_resolves_on_hanging_cuts() {
    for mixed in [false, true] {
        let (mut study, f, q) = fixture(mixed); let rows = limits(&q);
        let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
        let rho: Vec<_> = (0..study.cells()).map(|i| 0.35+0.025*(i%7) as f64).collect();
        let initial = study.operator().elasticity().scales().to_vec();
        let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        let exact = study.evaluate_response_limits(&rho, &cases, Default::default(), &mut c).unwrap();
        for i in 0..rho.len() {
            let mut a = rho.clone(); let mut b = rho.clone(); a[i] += 1e-4; b[i] -= 1e-4;
            let a = study.evaluate_response_limits(&a, &cases, Default::default(), &mut c).unwrap();
            let b = study.evaluate_response_limits(&b, &cases, Default::default(), &mut c).unwrap();
            close((a.constraint-b.constraint)/2e-4, exact.constraint_gradient[i], 3e-4);
            close((a.volume_fraction-b.volume_fraction)/2e-4, exact.volume_gradient[i], 2e-5);
        }
        assert_eq!(study.operator().elasticity().scales(), initial);
    }
}

#[test]
fn g0_single_signed_limits_match_existing_squared_response_pullbacks() {
    let (mut s, f, q) = fixture(false); let rho = vec![0.5; s.cells()];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    for sense in [ResponseLimitSense3::AtMost, ResponseLimitSense3::AtLeast] {
        for is_reaction in [false, true] {
            let row = [ResponseLimit3 { response: if is_reaction { BoundedResponse3::Reaction(&right) } else { BoundedResponse3::Displacement(&q) },
                sense, bound: 0.0, scale: 0.02 }];
            let case = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &row }];
            let a = s.evaluate_response_limits(&rho, &case, Default::default(), &mut c).unwrap();
            let sign = if sense == ResponseLimitSense3::AtMost { 1.0 } else { -1.0 };
            let target = a.measurements[0][0]-sign*row[0].scale;
            let d = [ResponseTarget3 { q: &q, target, scale: row[0].scale, weight: 1.0 }];
            let r = [ReactionTarget3 { mode: &right, target, scale: row[0].scale, weight: 1.0 }];
            let cases = [ResponseCase3 { force: &f, prescribed: Some(&motion), targets: if is_reaction { &[] } else { &d } }];
            let reactions: [&[ReactionTarget3<'_>]; 1] = [if is_reaction { &r } else { &[] }];
            let b = s.evaluate_responses_with_reactions(&rho, &cases, &reactions, Default::default(), &mut c).unwrap();
            for (&x, &y) in a.constraint_gradient.iter().zip(&b.gradient) { close(x, y, 1e-8); }
            assert_eq!(a.constraint, a.worst_residual);
        }
    }
}

#[test]
fn g3_independent_experiments_share_setup_but_cannot_cancel_limit_violations() {
    let (mut s, f, q) = fixture(false); let negative: Vec<_> = f.iter().map(|v| -v).collect();
    let rows = limits(&q); let rho = vec![0.45; s.cells()];
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows },
        ResponseLimitCase3 { force: &negative, prescribed: Some(&motion), limits: &rows }];
    let setups = Cell::new(0);
    let mut p = |s: SolveProgress| { if s.stage == "sdf3-preconditioner-start" { setups.set(setups.get()+1); } ControlFlow::Continue(()) };
    let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let a = s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c).unwrap();
    assert_eq!(setups.get(), 1); assert_eq!(a.displacements.len(), 2); assert_eq!(a.adjoints.len(), 2);
    let original_max = a.residuals.iter().flatten().copied().fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(a.worst_residual, original_max); assert!(a.constraint >= original_max);
    let reversed = [cases[1], cases[0]];
    let b = s.evaluate_response_limits(&rho, &reversed, Default::default(), &mut c).unwrap();
    close(a.constraint, b.constraint, 1e-10);
    for (&x, &y) in a.constraint_gradient.iter().zip(&b.constraint_gradient) { close(x, y, 1e-9); }
}

#[test]
fn g4_cancellation_and_exact_work_exhaustion_return_no_partial_limits_or_changed_material() {
    for stage in ["response-limits-primal", "response-limits-adjoint", "response-limits-publish"] {
        let (mut s, f, q) = fixture(false); let rho = vec![0.5; s.cells()]; let rows = limits(&q);
        let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
        let initial = s.operator().elasticity().scales().to_vec();
        let mut p = |p: SolveProgress| if p.stage == stage && (stage.ends_with("publish") || p.solve_iterations > 0) { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        assert!(matches!(s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c), Err(ResponseError3::Evaluation(EvaluationStop::Cancelled))));
        assert_eq!(s.operator().elasticity().scales(), initial); assert!(c.work().linear_iterations > 0);
        let mut p = |_| ControlFlow::Continue(());
        let mut c = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut p);
        assert!(matches!(s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c), Err(ResponseError3::Evaluation(EvaluationStop::TotalBudget { .. }))));
        assert_eq!(c.work().linear_iterations, 1); assert_eq!(s.operator().elasticity().scales(), initial);
    }
}

#[test]
fn g0_whole_family_admission_precedes_prescribed_or_reaction_callbacks() {
    let (mut s, f, q) = fixture(false); let rho = vec![0.5; s.cells()]; let calls = Cell::new(0);
    let mode = |p, n| { calls.set(calls.get()+1); right(p, n) };
    let mut rows = limits(&q); rows[0].response = BoundedResponse3::Reaction(&mode); rows[2].scale = f64::NAN;
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    assert!(s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c).is_err());
    assert_eq!(calls.get(), 0); assert_eq!(c.work().linear_solves, 0);
}

#[test]
fn g5_protected_regions_and_identical_designs_replay_without_geometry_changes() {
    let (s, f, q) = fixture(false); let mut labels = vec![PhysicalRegion3::Design; s.cells()];
    labels[0] = PhysicalRegion3::Solid; labels[1] = PhysicalRegion3::Void;
    let mut s = s.with_physical_regions(labels).unwrap(); let rho = vec![0.5; s.cells()]; let rows = limits(&q);
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let nodes = s.operator().elasticity().nodes().to_vec();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let a = s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c).unwrap();
    let b = s.evaluate_response_limits(&rho, &cases, Default::default(), &mut c).unwrap();
    assert_eq!(a.constraint.to_bits(), b.constraint.to_bits()); assert_eq!(a.constraint_gradient, b.constraint_gradient);
    assert_eq!(a.displacements, b.displacements); assert_eq!(a.measurements, b.measurements);
    assert_eq!(a.projected_rho[0], 1.0); assert_eq!(a.projected_rho[1], 0.0);
    assert_eq!(s.operator().elasticity().nodes(), nodes);
}
