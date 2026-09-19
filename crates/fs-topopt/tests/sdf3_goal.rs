use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{EvaluationStop, MultiLoadOcOptions, SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3, controlled_sdf3_optimality_criteria};
use fs_topopt::sdf3_goal::{GoalBodyLoad3, GoalRefinementError3, GoalRefinementOptions3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval { Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73) }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval {
        let v = if a == HeightAxis::Z { 1.0 } else { 0.0 }; Interval::new(v, v)
    }
}
fn y(_: [f64; 3]) -> [f64; 3] { [0.0, -1.0, 0.0] }
fn z(_: [f64; 3]) -> [f64; 3] { [0.0, 0.0, -1.0] }
fn loads() -> [GoalBodyLoad3<'static>; 2] { [GoalBodyLoad3 { density: &y, weight: 0.3 }, GoalBodyLoad3 { density: &z, weight: 0.7 }] }
fn build(t: &Octree3) -> AdaptiveElasticity3 {
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut p).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), t, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut q).unwrap()
}
fn setup() -> (Octree3, CutDensityStudy3<AdaptiveElasticity3>, fs_topopt::MultiLoadOcReport) {
    let tree = Octree3::uniform(1, 4, 4096).unwrap(); let op = build(&tree);
    let fy = op.body_load(&y, || ControlFlow::Continue(())).unwrap(); let fz = op.body_load(&z, || ControlFlow::Continue(())).unwrap();
    let mut study = CutDensityStudy3::new(op, 0.15, SimpParams::default()); let rho = vec![0.5; study.cells()];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let report = controlled_sdf3_optimality_criteria(&mut study,
        &[LoadCase { force: &fy, weight: 0.3 }, LoadCase { force: &fz, weight: 0.7 }], &rho,
        MultiLoadOcOptions { max_iterations: 2, ..Default::default() }, &mut c);
    assert!(report.history.len() > 1, "{report:?}"); (tree, study, report)
}
fn enriched(tree: &Octree3) -> AdaptiveElasticity3 {
    build(&tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(), || ControlFlow::Continue(())).unwrap())
}
#[test]
fn accepted_design_drives_real_goal_marking_without_mutating_the_coarse_state() {
    let (tree, study, accepted) = setup(); let scales = study.operator().scales().to_vec();
    let mut p = |_| ControlFlow::Continue(()); let mut control = SolveControl::new(SolveBudget::default(), &mut p);
    let estimate = study.estimate_compliance_enrichment(enriched(&tree), &loads(), &accepted.displacements,
        GoalRefinementOptions3::default(), &mut control).unwrap();
    assert_eq!(estimate.cases.len(), 2); assert_eq!(estimate.work.linear_solves, 2);
    assert!((estimate.coarse_value - accepted.history.last().unwrap().compliance).abs() < 1e-10);
    assert!((estimate.correction - (estimate.fine_value - estimate.coarse_value)).abs() < 1e-7);
    let marked = estimate.mark(0.5, 2, || ControlFlow::Continue(())).unwrap();
    assert!(!marked.marked.is_empty() && marked.marked.len() <= 2);
    assert!(marked.achieved_fraction > 0.0);
    let refined = tree.refined(&marked.marked, || ControlFlow::Continue(())).unwrap();
    assert!(refined.leaves().len() > tree.leaves().len());
    assert_eq!(study.operator().scales(), scales);
    let again = study.estimate_compliance_enrichment(enriched(&tree), &loads(), &accepted.displacements,
        GoalRefinementOptions3::default(), &mut control).unwrap();
    assert_eq!(again.marking_mass, estimate.marking_mass);
    assert_eq!(again.correction.to_bits(), estimate.correction.to_bits());
    assert!(again.work.linear_iterations > estimate.work.linear_iterations);
}
#[test]
fn enriched_work_exhaustion_is_exact_and_does_not_publish_a_partial_load_family() {
    let (tree, study, accepted) = setup(); let scales = study.operator().scales().to_vec();
    let mut p = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut p);
    let result = study.estimate_compliance_enrichment(enriched(&tree), &loads(), &accepted.displacements, GoalRefinementOptions3::default(), &mut c);
    assert!(matches!(result, Err(GoalRefinementError3::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(c.work().linear_iterations, 1); assert_eq!(study.operator().scales(), scales);
}
#[test]
fn cancellation_inside_enriched_krylov_retains_the_accepted_design() {
    let (tree, study, accepted) = setup(); let scales = study.operator().scales().to_vec();
    let mut p = |progress: SolveProgress| {
        if progress.stage == "sdf3-goal-elasticity" && progress.solve_iterations > 0 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
    };
    let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    assert!(matches!(study.estimate_compliance_enrichment(enriched(&tree), &loads(), &accepted.displacements,
        GoalRefinementOptions3::default(), &mut c), Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
    assert!(c.work().linear_iterations > 0); assert_eq!(study.operator().scales(), scales);
}
#[test]
fn wrong_body_law_and_incomplete_load_family_refuse_before_fine_solver_work() {
    let (tree, study, accepted) = setup();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let swapped = [GoalBodyLoad3 { density: &z, weight: 0.3 }, GoalBodyLoad3 { density: &y, weight: 0.7 }];
    assert!(matches!(study.estimate_compliance_enrichment(enriched(&tree), &swapped, &accepted.displacements,
        GoalRefinementOptions3::default(), &mut c), Err(GoalRefinementError3::Estimate(_))));
    assert_eq!(c.work().linear_iterations, 0);
    assert!(matches!(study.estimate_compliance_enrichment(enriched(&tree), &loads(), &accepted.displacements[..1],
        GoalRefinementOptions3::default(), &mut c), Err(GoalRefinementError3::Invalid(_))));
    assert_eq!(c.work().linear_solves, 0);
}
