use super::*;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3, adaptive::AdaptiveElasticity3, surface::ReferenceLoad3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
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
fn build(tree: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>) -> Result<AdaptiveElasticity3, GoalRefinementError3> {
    use fs_cutfem::elastic3::ElasticityError3;
    let mut poll = |_| cp();
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll)
        .map_err(ElasticityError3::from)?;
    Ok(AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &pinned,
        ElasticityOptions3::default(), &mut q)?)
}
fn pinned(p: [f64; 3]) -> bool { p[0] == 0.0 }
fn load(_: [f64; 3]) -> [f64; 3] { [0.0, 0.0, -1.0] }
fn laws() -> [GoalReferenceLoad3<'static>; 1] {
    [GoalReferenceLoad3 { load: ReferenceLoad3::body(&load), weight: 1.0 }]
}
fn params() -> SimpParams { SimpParams { penal: 2.0, beta: 1.0, e_min: 1e-3, ..Default::default() } }
fn fixture() -> (Octree3, CutDensityStudy3<AdaptiveElasticity3>, Vec<f64>) {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let op = build(&tree, &mut || ControlFlow::Continue(())).unwrap();
    let study = CutDensityStudy3::new(op, 0.15, params());
    let raw = vec![0.3; study.cells()];
    (tree, study, raw)
}
fn policy(refinements: usize) -> RobustAdaptiveOptions3 {
    let mut options = RobustAdaptiveOptions3 { max_refinements: refinements, max_marks: 1, ..Default::default() };
    options.refit.continuation.steps_per_stage = 3;
    options.refit.continuation.design.compliance_scale = 1000.0;
    options.refit.continuation.design.volume_limit = 0.4;
    options
}

#[test]
fn two_cycles_selectively_rebuild_and_reoptimize_matching_robust_fields() {
    let (mut tree, mut study, raw) = fixture();
    let mut expected_tree = tree.clone();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let options = policy(2);
    let report = controlled_adaptive_robust_design(&mut study, &mut tree, &laws(), &raw,
        &[params()], options, &mut control, build).unwrap();
    assert_eq!(report.stop, RobustAdaptiveStop3::RefinementLimit, "{report:?}");
    assert_eq!(report.grids.len(), 3);
    assert_eq!(report.probes.len(), 2);
    for (i, probe) in report.probes.iter().enumerate() {
        assert_eq!(probe.source_grid, i); assert!(probe.installed);
        let marking = probe.marking.as_ref().unwrap();
        assert_eq!(marking.marked.len(), 1);
        expected_tree = expected_tree.refined(&marking.marked, || ControlFlow::Continue(())).unwrap();
        assert_eq!(probe.target_background_cells, Some(expected_tree.leaves().len()));
        assert_eq!(report.grids[i+1].background_cells, expected_tree.leaves().len());
        assert!(report.grids[i+1].design.stages[0].history.len() > 1);
    }
    assert_eq!(tree.leaves(), expected_tree.leaves());
    let last = report.last().unwrap();
    assert_eq!(last.rho.len(), study.cells());
    assert!(last.evaluation.dilated_volume_fraction <= 0.4 * (1.0 + options.refit.continuation.design.optimizer.tolerance));
    let force = study.operator().reference_load(ReferenceLoad3::body(&load), || ControlFlow::Continue(())).unwrap();
    let fresh = study.evaluate_robust(&last.rho, &[crate::pipeline::LoadCase { force: &force, weight: 1.0 }],
        options.refit.continuation.design.robustness, &mut control).unwrap();
    assert_eq!(fresh.eroded.objective.displacements, last.evaluation.eroded.objective.displacements);
    assert_eq!(fresh.dilated_rho, last.evaluation.dilated_rho);
}

#[test]
fn final_publication_checkpoint_protects_the_previous_grid() {
    let (mut tree, mut study, raw) = fixture();
    let original = tree.clone(); let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-robust-adaptive-promote" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_adaptive_robust_design(&mut study, &mut tree, &laws(), &raw,
        &[params()], policy(1), &mut control, build).unwrap();
    assert_eq!(report.stop, RobustAdaptiveStop3::RefinementStopped, "{report:?}");
    assert_eq!(report.grids.len(), 1);
    assert!(report.rejected_design.as_ref().unwrap().last.is_some());
    assert!(!report.probes[0].installed);
    assert_eq!(tree.leaves(), original.leaves());
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(report.last().unwrap().rho.len(), study.cells());
}

#[test]
fn post_publication_cancellation_keeps_the_complete_new_grid() {
    let (mut tree, mut study, raw) = fixture();
    let initial_cells = study.cells();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-robust-adaptive-promoted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_adaptive_robust_design(&mut study, &mut tree, &laws(), &raw,
        &[params()], policy(2), &mut control, build).unwrap();
    assert_eq!(report.stop, RobustAdaptiveStop3::RefinementStopped, "{report:?}");
    assert_eq!(report.grids.len(), 2);
    assert!(report.probes[0].installed);
    assert!(study.cells() > initial_cells);
    assert_eq!(report.last().unwrap().rho.len(), study.cells());
    assert_eq!(report.grids[1].background_cells, tree.leaves().len());
}

#[test]
fn cumulative_physics_limit_cannot_be_renewed_by_refinement() {
    let (mut tree, mut study, raw) = fixture();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let baseline = controlled_adaptive_robust_design(&mut study, &mut tree, &laws(), &raw,
        &[params()], policy(0), &mut control, build).unwrap();
    let budget = baseline.work.linear_iterations + 1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: budget, ..Default::default() }, &mut poll);
    let report = controlled_adaptive_robust_design(&mut study, &mut tree, &laws(), &raw,
        &[params()], policy(2), &mut control, build).unwrap();
    assert_eq!(report.stop, RobustAdaptiveStop3::RefinementStopped);
    assert_eq!(report.work.linear_iterations, budget);
    assert_eq!(report.grids.len(), 1);
    assert_eq!(report.last().unwrap().rho, baseline.last().unwrap().rho);
}

#[test]
fn optional_accuracy_gate_probes_the_current_grid_even_at_zero_refinement_budget() {
    let (mut tree, mut study, raw) = fixture();
    let mut options = policy(0); options.eroded_relative_tolerance = Some(1e-12);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let first = controlled_adaptive_robust_design(&mut study, &mut tree, &laws(), &raw,
        &[params()], options, &mut control, build).unwrap();
    assert_eq!(first.probes.len(), 1, "{first:?}");
    let gap = first.probes[0].maximum_relative_gap;
    assert!(gap > 1e-12 && gap < 1.0);
    assert_eq!(first.stop, RobustAdaptiveStop3::RefinementLimit);
    options.eroded_relative_tolerance = Some((gap + 1.0) / 2.0);
    let passed = controlled_adaptive_robust_design(&mut study, &mut tree, &laws(), &raw,
        &[params()], options, &mut control, build).unwrap();
    assert_eq!(passed.stop, RobustAdaptiveStop3::ErodedComplianceTolerance, "{passed:?}");
    assert_eq!(passed.grids.len(), 1);
    assert_eq!(passed.probes[0].source_grid, 0);
    assert!(passed.probes[0].marking.is_none());
}

#[test]
fn zero_residual_signal_without_a_requested_accuracy_gate_is_not_a_pass() {
    let (mut tree, mut study, raw) = fixture();
    let zero = |_: [f64; 3]| [0.0; 3];
    let loads = [GoalReferenceLoad3 { load: ReferenceLoad3::body(&zero), weight: 1.0 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_adaptive_robust_design(&mut study, &mut tree, &loads, &raw,
        &[params()], policy(2), &mut control, build).unwrap();
    assert_eq!(report.stop, RobustAdaptiveStop3::NoRefinementSignal, "{report:?}");
    assert_eq!(report.grids.len(), 1);
    assert!(!report.probes[0].installed);
}

#[test]
fn per_case_accuracy_gate_cannot_cancel_opposing_compliance_changes() {
    use fs_dwr::elasticity3::GoalEstimate3;
    use std::collections::BTreeMap;
    let case = |coarse_value, fine_value| GoalEstimate3 {
        coarse_value, fine_value, dwr: 0.0, coarse_space: 0.0, goal_transfer: 0.0,
        algebraic: 0.0, identity_relative_defect: 0.0, field_residuals: [0.0; 4], cells: BTreeMap::new(),
    };
    // Policy-only test: a cancelled aggregate must not hide per-case changes.
    let estimate = ComplianceRefinement3 {
        cases: vec![case(10.0, 12.0), case(12.0, 10.0)], weights: vec![1.0, 1.0],
        coarse_value: 22.0, fine_value: 22.0, correction: 0.0, marking_mass: BTreeMap::new(),
        work: SolveWork::default(), preconditioner: crate::sdf3_goal::GoalPreconditioner3::Identity,
    };
    assert!(maximum_gap(&estimate) > 0.16);
    assert_eq!(relative_gap(f64::MAX, -f64::MAX), 2.0);
    assert!(relative_gap(f64::NAN, 0.0).is_infinite());
}
