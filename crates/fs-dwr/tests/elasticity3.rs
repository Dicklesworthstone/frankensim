use std::collections::BTreeMap;
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityError3, ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_dwr::elasticity3::{GoalError3, GoalFields3, GoalOptions3, estimate_goal3, dorfler3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval { Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73) }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 }; Interval::new(d, d)
    }
}
fn build(tree: &Octree3) -> AdaptiveElasticity3 {
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut p).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut q).unwrap()
}
fn pair() -> (AdaptiveElasticity3, AdaptiveElasticity3) {
    let t = Octree3::uniform(1, 4, 4096).unwrap();
    let mark = *t.leaves().iter().find(|c| c.index() == [0, 0, 1]).unwrap();
    let t = t.refined(&[mark], || ControlFlow::Continue(())).unwrap();
    let f = t.refined(&t.leaves().iter().copied().collect::<Vec<_>>(), || ControlFlow::Continue(())).unwrap();
    let (mut coarse, mut fine) = (build(&t), build(&f));
    coarse.set_scales(&(0..coarse.cells()).map(|i| 0.25 + 0.03 * i as f64).collect::<Vec<_>>()).unwrap();
    let scales = AdaptiveTransfer3::new(&coarse, &fine, 100_000, || ControlFlow::Continue(())).unwrap().inherited_scales();
    fine.set_scales(&scales).unwrap(); (coarse, fine)
}
fn solve(op: &AdaptiveElasticity3, body: &dyn Fn([f64; 3]) -> [f64; 3]) -> Vec<f64> {
    let b = op.body_load(body, || ControlFlow::Continue(())).unwrap();
    op.solve_controlled(&b, 1e-10, 20_000, 32, |_| ControlFlow::Continue(())).unwrap().coefficients().to_vec()
}
fn body(p: [f64; 3]) -> [f64; 3] { [p[1], -1.0, p[0] - 0.5] }
#[test]
fn compliance_dwr_and_consistency_reconstruct_the_real_two_grid_difference() {
    let (c, f) = pair(); let uc = solve(&c, &body); let uf = solve(&f, &body);
    let map = AdaptiveTransfer3::new(&c, &f, 100_000, || ControlFlow::Continue(())).unwrap();
    let report = estimate_goal3(&map, &body, &body, GoalFields3::compliance(&uc, &uf), GoalOptions3::default(), || ControlFlow::Continue(())).unwrap();
    let gap = report.fine_value - report.coarse_value;
    assert!((report.correction() - gap).abs() < 1e-7 * gap.abs().max(1.0));
    assert!(report.dwr.abs() > 1e-3);
    assert!(report.coarse_space.abs() > 1e-5); // Dropping this term breaks the identity.
    assert!((report.dwr - gap).abs() > 1e-5);
    assert!(report.identity_relative_defect < 1e-8);
    assert!(report.field_residuals.iter().all(|r| *r <= 1e-8));
    assert_eq!(report.cells.len(), c.cells());
    let mass: f64 = report.cells.values().map(|c| c.marking_mass).sum();
    assert!(mass + 1e-12 >= report.dwr.abs());
}
#[test]
fn general_goal_uses_its_own_adjoint_and_sign_reversal_preserves_marking() {
    let (c, f) = pair(); let uc = solve(&c, &body); let uf = solve(&f, &body);
    let goal = |p: [f64; 3]| [0.0, 0.0, if p[0] > 0.5 { 1.0 } else { 0.0 }];
    let zc = solve(&c, &goal); let zf = solve(&f, &goal);
    let map = AdaptiveTransfer3::new(&c, &f, 100_000, || ControlFlow::Continue(())).unwrap();
    let fields = GoalFields3 { coarse_primal: &uc, fine_primal: &uf, coarse_adjoint: &zc, fine_adjoint: &zf };
    let a = estimate_goal3(&map, &body, &goal, fields, GoalOptions3::default(), || ControlFlow::Continue(())).unwrap();
    let negc: Vec<_> = zc.iter().map(|v| -v).collect(); let negf: Vec<_> = zf.iter().map(|v| -v).collect();
    let b = estimate_goal3(&map, &body, &|p| goal(p).map(|v| -v),
        GoalFields3 { coarse_adjoint: &negc, fine_adjoint: &negf, ..fields }, GoalOptions3::default(), || ControlFlow::Continue(())).unwrap();
    assert!((a.correction() + b.correction()).abs() < 1e-8);
    assert!(a.identity_relative_defect < 1e-8 && b.identity_relative_defect < 1e-8);
    assert_eq!(a.mark(0.5, 100, || ControlFlow::Continue(())).unwrap().marked,
        b.mark(0.5, 100, || ControlFlow::Continue(())).unwrap().marked);
    // A primal is not the adjoint of this different observation.
    assert!(matches!(estimate_goal3(&map, &body, &goal, GoalFields3::compliance(&uc, &uf), GoalOptions3::default(), || ControlFlow::Continue(())), Err(GoalError3::FieldResidual { .. })));
}
#[test]
fn stale_fields_and_a_refiltered_material_distribution_are_not_enrichment_evidence() {
    let (c, mut f) = pair(); let uc = solve(&c, &body); let uf = solve(&f, &body);
    let map = AdaptiveTransfer3::new(&c, &f, 100_000, || ControlFlow::Continue(())).unwrap();
    let zero = vec![0.0; uc.len()];
    assert!(matches!(estimate_goal3(&map, &body, &body, GoalFields3::compliance(&zero, &uf), GoalOptions3::default(), || ControlFlow::Continue(())), Err(GoalError3::FieldResidual { field: "coarse-primal", .. })));
    drop(map);
    f.set_scales(&vec![1.0; f.cells()]).unwrap();
    let map = AdaptiveTransfer3::new(&c, &f, 100_000, || ControlFlow::Continue(())).unwrap();
    assert!(matches!(estimate_goal3(&map, &body, &body, GoalFields3::compliance(&uc, &uf), GoalOptions3::default(), || ControlFlow::Continue(())), Err(GoalError3::Invalid(_))));
}
#[test]
fn zero_signal_and_mark_budget_exhaustion_are_not_successful_refinement_targets() {
    let tree = Octree3::uniform(1, 2, 512).unwrap(); let keys: Vec<_> = tree.leaves().iter().copied().take(2).collect();
    let masses = BTreeMap::from([(keys[1], f64::MAX), (keys[0], f64::MAX)]);
    let limited = dorfler3(&masses, 1.0, 1, || ControlFlow::Continue(())).unwrap();
    assert_eq!(limited.marked, vec![keys[0]]); assert!(!limited.target_met);
    assert!((limited.achieved_fraction - 0.5).abs() < 1e-15);
    assert!(dorfler3(&masses, 1.0, 2, || ControlFlow::Continue(())).unwrap().target_met);
    assert!(!dorfler3(&BTreeMap::new(), 0.5, 10, || ControlFlow::Continue(())).unwrap().target_met);
    assert!(dorfler3(&masses, f64::NAN, 2, || ControlFlow::Continue(())).is_err());
    assert!(dorfler3(&BTreeMap::from([(keys[0], -1.0)]), 0.5, 2, || ControlFlow::Continue(())).is_err());
}
#[test]
fn cancellation_and_absent_enrichment_cannot_produce_a_partial_estimate() {
    let (c, f) = pair(); let uc = solve(&c, &body); let uf = solve(&f, &body);
    let map = AdaptiveTransfer3::new(&c, &f, 100_000, || ControlFlow::Continue(())).unwrap();
    let mut calls = 0;
    let result = estimate_goal3(&map, &body, &body, GoalFields3::compliance(&uc, &uf), GoalOptions3::default(), || {
        calls += 1; if calls > 40 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
    });
    assert!(matches!(result, Err(GoalError3::Physics(ElasticityError3::Cancelled))));
    let map = AdaptiveTransfer3::new(&c, &c, 100_000, || ControlFlow::Continue(())).unwrap();
    assert!(matches!(estimate_goal3(&map, &body, &body, GoalFields3::compliance(&uc, &uc), GoalOptions3::default(), || ControlFlow::Continue(())), Err(GoalError3::Invalid(_))));
}
