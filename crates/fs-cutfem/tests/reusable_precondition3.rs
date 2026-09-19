//! Real adaptive operators: retained topology, refreshed numerical factors.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{
    AdaptivePreconditionError3, AdaptiveSolveOptions3, AdaptiveSolveSpace3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::{LinearOp, two_level::TwoLevelError};
use fs_sparse::precond::Precond;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval {
        let d = if a == HeightAxis::Z { 1.0 } else { 0.0 }; Interval::new(d, d)
    }
}
fn build(tree: &Octree3) -> AdaptiveElasticity3 {
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q).unwrap()
}
fn trees() -> (Octree3, Octree3) {
    let coarse = Octree3::uniform(1, 4, 4096).unwrap();
    let fine = coarse.refined(&[*coarse.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap();
    (coarse, fine)
}
fn action(p: &impl Precond, n: usize) -> Vec<f64> {
    let r: Vec<_> = (0..n).map(|i| 0.1 + (i % 7) as f64 / 9.0).collect();
    let mut z = vec![0.0; n]; p.apply(&r, &mut z); z
}

#[test]
fn retained_space_matches_fresh_transfer_at_every_density() {
    let (ct, ft) = trees(); let coarse = build(&ct); let opts = AdaptiveSolveOptions3::default();
    let mut space = AdaptiveSolveSpace3::two_level(build(&ft), &coarse, opts, || ControlFlow::Continue(())).unwrap();
    let leaves = space.elasticity().leaves().to_vec(); let entries = space.transfer_entries();
    let mut prior = None;
    for stage in 0..3 {
        let scales: Vec<_> = (0..space.elasticity().cells()).map(|i| 0.2 + 0.1 * ((i + stage) % 7) as f64).collect();
        space.set_scales(&scales).unwrap();
        let prepared = space.prepare(|_| ControlFlow::Continue(())).unwrap();
        let transfer = AdaptiveTransfer3::new(&coarse, space.elasticity(), opts.max_transfer_terms, || ControlFlow::Continue(())).unwrap();
        let fresh = transfer.prepare_two_level(opts.two_level, opts.max_diagonal_contributions, |_| ControlFlow::Continue(())).unwrap();
        let actual = action(&prepared, space.n());
        assert_eq!(actual, action(&fresh, space.n()));
        assert_eq!(prepared.work().operator_applications, space.coarse_dofs());
        assert!(std::ptr::eq(prepared.operator(), space.elasticity()));
        if let Some(previous) = prior { assert_ne!(actual, previous, "changed stiffness must change numeric preparation"); }
        prior = Some(actual);
        assert_eq!(space.elasticity().leaves(), leaves);
        assert_eq!(space.transfer_entries(), entries);
    }
}

#[test]
fn one_preparation_serves_independent_loads_with_true_residuals() {
    let (ct, ft) = trees();
    let space = AdaptiveSolveSpace3::two_level(build(&ft), &build(&ct), AdaptiveSolveOptions3::default(), || ControlFlow::Continue(())).unwrap();
    let prepared = space.prepare(|_| ControlFlow::Continue(())).unwrap();
    let setup = prepared.work();
    for direction in [[0.0, -1.0, 0.0], [0.0, 0.0, -1.0]] {
        let rhs = space.elasticity().body_load(&|_| direction, || ControlFlow::Continue(())).unwrap();
        let mut cg = fs_solver::CgState::new(&space, &prepared, &rhs);
        assert!(cg.run(&space, &prepared, 1e-12, 10_000).converged);
        assert!(space.elasticity().field_residual(&cg.x, &rhs, || ControlFlow::Continue(())).unwrap() < 1e-9);
    }
    assert_eq!(prepared.work(), setup);
}

#[test]
fn cancellation_retains_scales_and_a_retry_rebuilds_complete_factors() {
    let (ct, ft) = trees();
    let space = AdaptiveSolveSpace3::two_level(build(&ft), &build(&ct), AdaptiveSolveOptions3::default(), || ControlFlow::Continue(())).unwrap();
    let before = space.elasticity().scales().to_vec(); let mut spent = 0;
    let failed = space.prepare(|w| {
        spent = w.operator_applications;
        if spent >= 3 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
    });
    assert!(matches!(failed, Err(AdaptivePreconditionError3::Coarse(TwoLevelError::Cancelled))));
    assert_eq!(spent, 3); assert_eq!(space.elasticity().scales(), before);
    let a = space.prepare(|_| ControlFlow::Continue(())).unwrap();
    let b = space.prepare(|_| ControlFlow::Continue(())).unwrap();
    assert_eq!(action(&a, space.n()), action(&b, space.n()));
}

#[test]
fn jacobi_matches_existing_owner_and_invalid_updates_are_atomic() {
    let (_, ft) = trees(); let mut space = AdaptiveSolveSpace3::jacobi(build(&ft), 100_000_000);
    let before = space.elasticity().scales().to_vec(); let mut bad = before.clone(); bad[0] = f64::NAN;
    assert!(space.set_scales(&bad).is_err()); assert_eq!(space.elasticity().scales(), before);
    let p = space.prepare(|_| ControlFlow::Continue(())).unwrap();
    let direct = space.elasticity().prepare_jacobi(100_000_000, || ControlFlow::Continue(())).unwrap();
    assert_eq!(action(&p, space.n()), action(&direct, space.n()));
    assert_eq!(space.coarse_dofs(), 0); assert_eq!(p.work().operator_applications, 0);
}

#[test]
fn exhausted_symbolic_or_numeric_budgets_do_not_silently_fall_back() {
    let (ct, ft) = trees(); let coarse = build(&ct);
    let mut opts = AdaptiveSolveOptions3::default(); opts.two_level.max_operator_applications = 1;
    assert!(matches!(AdaptiveSolveSpace3::two_level(build(&ft), &coarse, opts, || ControlFlow::Continue(())),
        Err(AdaptivePreconditionError3::Coarse(TwoLevelError::Budget(_)))));
    assert!(AdaptiveSolveSpace3::two_level(build(&ft), &coarse, AdaptiveSolveOptions3::default(), || ControlFlow::Break(())).is_err());
    let jacobi = AdaptiveSolveSpace3::jacobi(build(&ft), 0);
    assert!(matches!(jacobi.prepare(|_| ControlFlow::Continue(())), Err(AdaptivePreconditionError3::Physics(_))));
}
