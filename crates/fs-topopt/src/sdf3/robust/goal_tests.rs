use super::*;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3, surface::ReferenceLoad3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{EvaluationStop, SolveBudget, SolveProgress};
use crate::sdf3::PhysicalRegion3;

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
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut p).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q).unwrap()
}
fn fine(tree: &Octree3) -> AdaptiveElasticity3 {
    build(&tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(), || ControlFlow::Continue(())).unwrap())
}
fn y(_: [f64; 3]) -> [f64; 3] { [0.0, -1.0, 0.0] }
fn z(_: [f64; 3]) -> [f64; 3] { [0.0, 0.0, -1.0] }
fn laws() -> [GoalReferenceLoad3<'static>; 2] {
    [GoalReferenceLoad3 { load: ReferenceLoad3::body(&y), weight: 0.3 },
     GoalReferenceLoad3 { load: ReferenceLoad3::body(&z), weight: 0.7 }]
}
fn setup(tree: &Octree3, params: SimpParams, protected: bool) -> CutDensityStudy3<AdaptiveElasticity3> {
    let study = CutDensityStudy3::new(build(tree), 0.15, params);
    if !protected { return study; }
    let mut mask = vec![PhysicalRegion3::Design; study.cells()];
    mask[0] = PhysicalRegion3::Solid; mask[1] = PhysicalRegion3::Void;
    study.with_physical_regions(mask).unwrap()
}

#[test]
fn eroded_and_zero_offset_goal_experiments_match_the_existing_shifted_model() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let params = SimpParams { beta: 3.0, penal: 2.0, ..Default::default() };
    for offset in [0.0, 0.1] {
        for protected in [false, true] {
            let mut study = setup(&tree, SimpParams::default(), protected);
            let original = study.params(); let scales = study.operator().scales().to_vec();
            let rho: Vec<_> = (0..study.cells()).map(|i| 0.35 + 0.03 * i as f64).collect();
            let mut poll = |_| ControlFlow::Continue(());
            let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
            let result = study.estimate_robust_reference_compliance_enrichment(fine(&tree), &laws(),
                &rho, params, RobustOptions3 { threshold_offset: offset }, GoalRefinementOptions3::default(), &mut control).unwrap();
            let shifted = SimpParams { eta: params.eta + offset, ..params };
            let mut reference = setup(&tree, shifted, protected);
            let fy = reference.operator().body_load(&y, || ControlFlow::Continue(())).unwrap();
            let fz = reference.operator().body_load(&z, || ControlFlow::Continue(())).unwrap();
            let eroded = reference.evaluate(&rho, &[LoadCase { force: &fy, weight: 0.3 },
                LoadCase { force: &fz, weight: 0.7 }], &mut control).unwrap();
            let expected = reference.estimate_reference_compliance_enrichment(fine(&tree), &laws(),
                &eroded.objective.displacements, GoalRefinementOptions3::default(), &mut control).unwrap();
            assert_eq!(result.eroded.objective.displacements, eroded.objective.displacements);
            assert_eq!(result.eroded.projected_rho, eroded.projected_rho);
            assert_eq!(result.estimate.fine_value.to_bits(), expected.fine_value.to_bits());
            assert_eq!(result.estimate.marking_mass, expected.marking_mass);
            assert!((result.estimate.coarse_value - result.eroded.objective.compliance).abs()
                < 1e-9 * result.eroded.objective.compliance);
            assert!((result.estimate.correction - (result.estimate.fine_value - result.estimate.coarse_value)).abs()
                < 1e-7 * result.estimate.fine_value.abs());
            assert!(!result.estimate.mark(0.5, 2, || ControlFlow::Continue(())).unwrap().marked.is_empty());
            assert_eq!(study.params().beta, original.beta);
            assert_eq!(study.params().eta, original.eta);
            assert_eq!(study.operator().scales(), scales);
            if protected { assert_eq!(result.eroded.projected_rho[0], 1.0); assert_eq!(result.eroded.projected_rho[1], 0.0); }
        }
    }
}

#[test]
fn cancelled_enriched_solve_restores_a_preexisting_nominal_solution() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mut study = setup(&tree, SimpParams::default(), false);
    let rho = vec![0.45; study.cells()];
    let force = study.operator().body_load(&y, || ControlFlow::Continue(())).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    study.evaluate(&rho, &[LoadCase { force: &force, weight: 1.0 }], &mut control).unwrap();
    let original = study.params(); let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-goal-elasticity" && p.solve_iterations > 0 {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let result = study.estimate_robust_reference_compliance_enrichment(fine(&tree), &laws(), &rho,
        SimpParams { beta: 4.0, ..Default::default() }, RobustOptions3::default(), GoalRefinementOptions3::default(), &mut control);
    assert!(matches!(result, Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
    assert!(control.work().linear_iterations > 0);
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(study.params().eta, original.eta); assert_eq!(study.params().beta, original.beta);
}

#[test]
fn exhausted_shared_budget_never_publishes_a_partial_eroded_family() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mut study = setup(&tree, SimpParams::default(), false);
    let rho = vec![0.45; study.cells()]; let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    let result = study.estimate_robust_reference_compliance_enrichment(fine(&tree), &laws(), &rho,
        study_params(), RobustOptions3::default(), GoalRefinementOptions3::default(), &mut control);
    assert!(matches!(result, Err(GoalRefinementError3::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(control.work().linear_iterations, 1);
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(study.params().beta, SimpParams::default().beta);
}
fn study_params() -> SimpParams { SimpParams { beta: 4.0, ..Default::default() } }

#[test]
fn bad_reference_family_refuses_before_any_physical_callback() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mut study = setup(&tree, SimpParams::default(), false);
    let rho = vec![0.45; study.cells()];
    let forbidden = |_: [f64; 3]| -> [f64; 3] { panic!("invalid family invoked its load law") };
    let loads = [GoalReferenceLoad3 { load: ReferenceLoad3::body(&forbidden), weight: f64::NAN }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(matches!(study.estimate_robust_reference_compliance_enrichment(fine(&tree), &loads, &rho,
        study_params(), RobustOptions3::default(), GoalRefinementOptions3::default(), &mut control),
        Err(GoalRefinementError3::Invalid(_))));
    assert_eq!(control.work().linear_solves, 0);
}
