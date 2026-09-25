use super::*;
use fs_cutfem::elastic3::{ElasticityError3, ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{SolveBudget, SolveProgress};
use std::cell::Cell;

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
fn build(tree: &Octree3, checkpoint: &mut dyn FnMut() -> ControlFlow<()>)
    -> Result<AdaptiveElasticity3, GoalRefinementError3>
{
    let mut poll = |_| checkpoint();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll,
    ).map_err(ElasticityError3::from)?;
    Ok(AdaptiveElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut quadrature,
    )?)
}
fn fixture() -> (CutDensityStudy3<AdaptiveElasticity3>, Octree3) {
    let tree = Octree3::uniform(1, 4, 1000).unwrap();
    let operator = build(&tree, &mut || ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(operator, 0.15, SimpParams::default()), tree)
}
fn policy() -> AdaptiveContinuationOptions3 {
    AdaptiveContinuationOptions3 {
        optimization: MultiLoadOcOptions {
            max_iterations: 1, change_tolerance: 0.0, ..Default::default()
        }, max_marks: 1, ..Default::default()
    }
}
fn schedule() -> [SimpParams; 2] {
    [(1.0, 1.0), (2.0, 2.0)].map(|(penal, beta)| SimpParams {
        penal, beta, ..Default::default()
    })
}

#[test]
fn g1_g5_checkpoints_precede_more_physics_and_do_not_change_numerical_results() {
    let (mut observed, mut tree) = fixture();
    let raw = vec![0.5; observed.cells()];
    let body = |_: [f64; 3]| [0.0, 0.0, -1.0];
    let loads = [GoalReferenceLoad3 { load: ReferenceLoad3::body(&body), weight: 1.0 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let published = Cell::new(0);
    let mut endpoints = Vec::new();
    let report = controlled_adaptive_sdf3_continuation_observed(
        &mut observed, &mut tree, &loads, &raw, &schedule(), policy(), &mut control,
        |tree, checkpoint| {
            assert_eq!(published.get(), 1, "first stage is durable before enrichment");
            build(tree, checkpoint)
        },
        |study, tree, report| {
            let n = report.continuation.stages.len();
            assert_eq!(n, published.get() + 1);
            assert_eq!(report.continuation.termination, ContinuationTermination::ScheduleComplete);
            assert!(report.refinements.iter().all(|r| r.installed));
            assert_eq!(report.refinements.len(), n - 1);
            let last = report.continuation.last.as_ref().unwrap();
            assert_eq!(last.rho.len(), study.cells());
            assert!(study.operator().leaves().iter().all(|leaf| tree.leaves().contains(leaf)));
            assert!(report.continuation.work.linear_iterations > 0);
            endpoints.push((last.rho.clone(), last.displacements.clone(), tree.leaves().to_vec()));
            published.set(n);
            Ok::<(), ()>(())
        },
    ).unwrap();
    assert_eq!(published.get(), 2);
    assert_eq!(report.continuation.work, control.work());
    let (mut ordinary, mut ordinary_tree) = fixture();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let reference = controlled_adaptive_sdf3_continuation(
        &mut ordinary, &mut ordinary_tree, &loads, &raw, &schedule(), policy(),
        &mut control, build,
    );
    assert_eq!(format!("{report:?}"), format!("{reference:?}"));
    assert_eq!(observed.operator().scales(), ordinary.operator().scales());
    assert_eq!(tree.leaves(), ordinary_tree.leaves());
    assert_eq!(endpoints[1].0, report.continuation.last.as_ref().unwrap().rho);
    assert_eq!(endpoints[1].1, report.continuation.last.as_ref().unwrap().displacements);
}

#[test]
fn g4_checkpoint_failure_returns_the_storage_error_without_starting_a_proposal() {
    let (mut study, mut tree) = fixture();
    let original = tree.leaves().to_vec();
    let raw = vec![0.5; study.cells()];
    let body = |_: [f64; 3]| [0.0, 0.0, -1.0];
    let loads = [GoalReferenceLoad3 { load: ReferenceLoad3::body(&body), weight: 1.0 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let mut callbacks = 0;
    let mut solved = None;
    let error = controlled_adaptive_sdf3_continuation_observed(
        &mut study, &mut tree, &loads, &raw, &schedule(), policy(), &mut control,
        |_, _| -> Result<AdaptiveElasticity3, GoalRefinementError3> {
            panic!("no geometry after failed persistence")
        },
        |_, _, report| {
            callbacks += 1;
            solved = report.continuation.last.clone();
            Err("checkpoint-storage-full")
        },
    ).unwrap_err();
    assert_eq!(error, "checkpoint-storage-full");
    assert_eq!(callbacks, 1);
    assert_eq!(tree.leaves(), original);
    assert_eq!(solved.unwrap().rho.len(), study.cells());
}

#[test]
fn g4_unsolved_or_precancelled_states_are_never_checkpoints() {
    for cancelled in [false, true] {
        let (mut study, mut tree) = fixture();
        let raw = vec![0.5; study.cells()];
        let body = |_: [f64; 3]| [0.0, 0.0, -1.0];
        let loads = [GoalReferenceLoad3 { load: ReferenceLoad3::body(&body), weight: 1.0 }];
        let mut poll = |_: SolveProgress| if cancelled {
            ControlFlow::Break(())
        } else { ControlFlow::Continue(()) };
        let mut control = SolveControl::new(
            SolveBudget { total_iterations: 0, ..Default::default() }, &mut poll,
        );
        let report = controlled_adaptive_sdf3_continuation_observed(
            &mut study, &mut tree, &loads, &raw, &schedule(), policy(), &mut control,
            build,
            |_, _, _| -> Result<(), ()> { panic!("no solved stage to checkpoint") },
        ).unwrap();
        assert!(report.continuation.last.is_none());
        assert_eq!(report.continuation.work.linear_iterations, 0);
        assert_ne!(report.continuation.termination, ContinuationTermination::ScheduleComplete);
    }
}
