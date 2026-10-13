use super::*;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::ElasticityOptions3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{SolveBudget, SolveProgress};
use crate::sdf3::PhysicalRegion3;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let d = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(d, d)
    }
}
pub(super) fn build(tree: &Octree3) -> AdaptiveElasticity3 {
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q).unwrap()
}
fn y(_: [f64; 3]) -> [f64; 3] { [0.0, -1.0, 0.0] }
fn z(_: [f64; 3]) -> [f64; 3] { [0.0, 0.0, -1.0] }
pub(super) fn laws() -> [ReferenceLoad3<'static>; 2] { [ReferenceLoad3::body(&y), ReferenceLoad3::body(&z)] }
pub(super) fn params() -> SimpParams { SimpParams { beta: 2.0, penal: 2.0, e_min: 1e-3, ..Default::default() } }
pub(super) fn fixture() -> (Octree3, Octree3, CutDensityStudy3<AdaptiveElasticity3>, Vec<f64>) {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let fine = tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(), || ControlFlow::Continue(())).unwrap();
    let study = CutDensityStudy3::new(build(&tree), 0.15, params());
    let rho = vec![0.45; study.cells()];
    (tree, fine, study, rho)
}
fn caps() -> ComplianceCaps3 { ComplianceCaps3::new(vec![1e6, 2e6], 0.02, 2).unwrap() }

#[test]
fn real_eroded_enrichment_checks_each_case_and_normalizes_absolute_marking() {
    let (_, fine, mut study, rho) = fixture();
    let scales = study.operator().scales().to_vec();
    let original = study.params();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let caps = caps();
    let result = study.estimate_compliance_cap_enrichment(build(&fine), &laws(), &rho, params(),
        RobustOptions3::default(), &caps, GoalRefinementOptions3::default(), &mut control).unwrap();
    assert_eq!(result.checks.len(), 2);
    assert!(result.stable_and_feasible(1e-6, 1.0));
    let mut expected = BTreeMap::new();
    for (i, case) in result.assessment.estimate.cases.iter().enumerate() {
        assert_eq!(result.checks[i].coarse, case.coarse_value);
        assert_eq!(result.checks[i].fine, case.fine_value);
        assert_eq!(result.checks[i].cap, caps.caps()[i]);
        for (&key, cell) in &case.cells {
            *expected.entry(key).or_insert(0.0) += cell.marking_mass / caps.caps()[i];
        }
    }
    assert_eq!(result.marking_mass, expected);
    assert!(!result.mark(0.5, 2, || ControlFlow::Continue(())).unwrap().marked.is_empty());
    assert_eq!(result.work, control.work());
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(study.params().beta, original.beta);
    assert_eq!(study.params().eta, original.eta);
}

#[test]
fn one_bad_case_or_an_unstable_case_refuses_even_with_a_good_mean() {
    let (_, fine, mut study, rho) = fixture();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let mut result = study.estimate_compliance_cap_enrichment(build(&fine), &laws(), &rho, params(),
        RobustOptions3::default(), &caps(), GoalRefinementOptions3::default(), &mut control).unwrap();
    // Policy falsifier only: favorable aggregate diagnostics must not waive
    // an individual requirement. It does not purport to be a physical solve.
    result.coarse_envelope = 0.5;
    result.fine_envelope = 0.5;
    result.checks[0].coarse = 0.9 * result.checks[0].cap;
    result.checks[0].fine = result.checks[0].cap.next_up();
    result.checks[0].normalized_change = 0.0;
    assert!(!result.stable_and_feasible(0.9, 1.0));
    result.checks[0].fine = 0.9 * result.checks[0].cap;
    result.checks[0].normalized_change = 0.2;
    assert!(!result.stable_and_feasible(0.1, 0.1));
    assert!(!result.stable_and_feasible(f64::NAN, 1.0));
    assert!(!result.stable_and_feasible(0.1, f64::INFINITY));
}

#[test]
fn protected_regions_and_explicit_model_are_used_without_changing_source() {
    let (_, fine, study, rho) = fixture();
    let mut mask = vec![PhysicalRegion3::Design; study.cells()];
    mask[0] = PhysicalRegion3::Solid; mask[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(mask).unwrap();
    let scales = study.operator().scales().to_vec();
    let shifted = SimpParams { beta: 3.0, ..params() };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let result = study.estimate_compliance_cap_enrichment(build(&fine), &laws(), &rho, shifted,
        RobustOptions3::default(), &caps(), GoalRefinementOptions3::default(), &mut control).unwrap();
    assert_eq!(result.assessment.eroded.projected_rho[0], 1.0);
    assert_eq!(result.assessment.eroded.projected_rho[1], 0.0);
    assert_eq!(result.assessment.params.beta, shifted.beta);
    assert_eq!(study.params().beta, params().beta);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn same_grid_and_bad_load_count_refuse_before_load_callbacks() {
    let (tree, fine, mut study, rho) = fixture();
    let forbidden = |_: [f64; 3]| -> [f64; 3] { panic!("refused input invoked physical law") };
    let laws = [ReferenceLoad3::body(&forbidden); 2];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(study.estimate_compliance_cap_enrichment(build(&tree), &laws, &rho, params(),
        RobustOptions3::default(), &caps(), GoalRefinementOptions3::default(), &mut control).is_err());
    assert!(study.estimate_compliance_cap_enrichment(build(&fine), &laws[..1], &rho, params(),
        RobustOptions3::default(), &caps(), GoalRefinementOptions3::default(), &mut control).is_err());
    assert_eq!(control.work().linear_solves, 0);
}

#[test]
fn cancellation_and_budget_exhaustion_return_no_partially_assessed_family() {
    let (_, fine, mut study, rho) = fixture();
    let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-cap-enrichment-check" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(matches!(study.estimate_compliance_cap_enrichment(build(&fine), &laws(), &rho, params(),
        RobustOptions3::default(), &caps(), GoalRefinementOptions3::default(), &mut control),
        Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
    assert!(control.work().linear_iterations > 0);
    assert_eq!(study.operator().scales(), scales);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    assert!(matches!(study.estimate_compliance_cap_enrichment(build(&fine), &laws(), &rho, params(),
        RobustOptions3::default(), &caps(), GoalRefinementOptions3::default(), &mut control),
        Err(GoalRefinementError3::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(control.work().linear_iterations, 1);
    assert_eq!(study.operator().scales(), scales);
}

#[path = "refit_tests.rs"]
mod refit;
