use super::*;
use super::super::continuation::RobustContinuationStop3;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3, adaptive::AdaptiveElasticity3, surface::ReferenceLoad3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::sdf3::PhysicalRegion3;
use crate::{SolveBudget, SolveProgress};

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
fn y(_: [f64; 3]) -> [f64; 3] { [0.0, -1.0, 0.0] }
fn z(_: [f64; 3]) -> [f64; 3] { [0.0, 0.0, -1.0] }
fn laws() -> [GoalReferenceLoad3<'static>; 2] {
    [GoalReferenceLoad3 { load: ReferenceLoad3::body(&y), weight: 0.3 },
     GoalReferenceLoad3 { load: ReferenceLoad3::body(&z), weight: 0.7 }]
}
fn params() -> SimpParams { SimpParams { beta: 2.0, penal: 2.0, ..Default::default() } }
fn policy(steps: usize) -> RobustRefitOptions3 {
    let mut options = RobustRefitOptions3::default();
    options.continuation.steps_per_stage = steps;
    options.continuation.design.volume_limit = 0.4;
    options.continuation.design.compliance_scale = 1000.0;
    options.continuation.design.density_floor = 1e-5;
    options
}
fn fixture() -> (Octree3, Octree3, CutDensityStudy3<AdaptiveElasticity3>, Vec<f64>) {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let fine = tree.refined(&[*tree.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap();
    let study = CutDensityStudy3::new(build(&tree), 0.15, SimpParams::default());
    let rho = vec![0.6; study.cells()];
    (tree, fine, study, rho)
}

#[test]
fn selective_refit_restores_actual_dilated_volume_and_reintegrates_independent_loads() {
    let (_, fine, source, rho) = fixture();
    let scales = source.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let result = source.refit_robust_reference(build(&fine), &rho, params(), &laws(), policy(8), &mut control).unwrap();
    assert_eq!(result.report.stop, RobustContinuationStop3::ScheduleComplete, "{:?}", result.report);
    assert!(result.study.cells() > source.cells());
    assert_eq!(result.study.filter_radius, source.filter_radius);
    assert!(result.inherited_rho.iter().all(|r| *r == 0.6));
    let stage = &result.report.stages[0];
    assert!(stage.restoration_scale < 1.0 && stage.gradient_check.numerical.passed());
    assert!(stage.history.len() > 1, "the fine optimizer must take real steps");
    let selected = result.report.last.as_ref().unwrap();
    assert!(selected.evaluation.dilated_volume_fraction <= 0.4 * (1.0 + policy(8).continuation.design.optimizer.tolerance));
    let mut independent = CutDensityStudy3::new(build(&fine), 0.15, params());
    let fy = independent.operator().reference_load(ReferenceLoad3::body(&y), || ControlFlow::Continue(())).unwrap();
    let fz = independent.operator().reference_load(ReferenceLoad3::body(&z), || ControlFlow::Continue(())).unwrap();
    let fresh = independent.evaluate_robust(&selected.rho, &[LoadCase { force: &fy, weight: 0.3 },
        LoadCase { force: &fz, weight: 0.7 }], policy(8).continuation.design.robustness, &mut control).unwrap();
    assert_eq!(fresh.eroded.objective.displacements, selected.evaluation.eroded.objective.displacements);
    assert_eq!(fresh.dilated_rho, selected.evaluation.dilated_rho);
    assert_eq!(source.operator().scales(), scales);
    assert_eq!(source.params().beta, SimpParams::default().beta);
}

#[test]
fn parent_map_preserves_authored_regions_and_sub_millidensity_raw_controls() {
    let (_, fine, source, mut rho) = fixture();
    rho[0] = 2e-5;
    let mut regions = vec![PhysicalRegion3::Design; source.cells()];
    regions[0] = PhysicalRegion3::Solid; regions[1] = PhysicalRegion3::Void;
    let source = source.with_physical_regions(regions.clone()).unwrap();
    let op = build(&fine);
    let parents = AdaptiveTransfer3::new(source.operator(), &op, 2_000_000,
        || ControlFlow::Continue(())).unwrap().parents().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let result = source.refit_robust_reference(op, &rho, params(), &laws(), policy(0), &mut control).unwrap();
    assert_eq!(result.report.stop, RobustContinuationStop3::ScheduleComplete, "{:?}", result.report);
    let last = result.report.last.as_ref().unwrap();
    for (i, &parent) in parents.iter().enumerate() {
        assert_eq!(result.inherited_rho[i], rho[parent]);
        assert_eq!(result.study.physical_regions().unwrap()[i], regions[parent]);
        for field in [&last.evaluation.eroded.projected_rho, &last.evaluation.nominal_rho, &last.evaluation.dilated_rho] {
            match regions[parent] {
                PhysicalRegion3::Solid => assert_eq!(field[i], 1.0),
                PhysicalRegion3::Void => assert_eq!(field[i], 0.0),
                PhysicalRegion3::Design => {}
            }
        }
    }
}

#[test]
fn unchanged_active_space_and_oversized_load_family_refuse_before_callbacks() {
    let (tree, fine, source, rho) = fixture();
    let forbidden = |_: [f64; 3]| -> [f64; 3] { panic!("refused input invoked a load") };
    let loads = [GoalReferenceLoad3 { load: ReferenceLoad3::body(&forbidden), weight: 1.0 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(source.refit_robust_reference(build(&tree), &rho, params(), &loads, policy(0), &mut control).is_err());
    let mut options = policy(0); options.max_load_cases = 0;
    assert!(source.refit_robust_reference(build(&fine), &rho, params(), &loads, options, &mut control).is_err());
    assert_eq!(control.work().linear_solves, 0);
}

#[test]
fn cancellation_during_fine_optimization_retains_candidate_prefix_not_source_mutation() {
    let (_, fine, source, rho) = fixture();
    let scales = source.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-robust-design-accepted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let candidate = source.refit_robust_reference(build(&fine), &rho, params(), &laws(), policy(8), &mut control).unwrap();
    assert_eq!(candidate.report.stop, RobustContinuationStop3::EvaluationStopped);
    assert_eq!(candidate.report.stages[0].history.len(), 2);
    assert!(candidate.report.last.is_some());
    assert!(control.work().linear_iterations > 0);
    assert_eq!(source.operator().scales(), scales);
}

#[test]
fn transfer_and_solve_budgets_never_reset_or_publish_unfunded_results() {
    let (_, fine, source, rho) = fixture();
    let mut options = policy(0); options.max_transfer_terms = 1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(source.refit_robust_reference(build(&fine), &rho, params(), &laws(), options, &mut control).is_err());
    assert_eq!(control.work().linear_solves, 0);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    let candidate = source.refit_robust_reference(build(&fine), &rho, params(), &laws(), policy(0), &mut control).unwrap();
    assert_eq!(candidate.report.stop, RobustContinuationStop3::EvaluationStopped);
    assert!(candidate.report.last.is_none());
    assert_eq!(control.work().linear_iterations, 1);
}

#[test]
fn impossible_protected_volume_keeps_the_entire_source_available() {
    let (_, fine, source, rho) = fixture();
    let count = source.cells();
    let source = source.with_physical_regions(vec![PhysicalRegion3::Solid; count]).unwrap();
    let scales = source.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let candidate = source.refit_robust_reference(build(&fine), &rho, params(), &laws(), policy(0), &mut control).unwrap();
    assert_eq!(candidate.report.stop, RobustContinuationStop3::RestorationFailed);
    assert!(candidate.report.last.is_none());
    assert_eq!(source.operator().scales(), scales);
    assert_eq!(source.physical_regions().unwrap(), vec![PhysicalRegion3::Solid; count]);
}
