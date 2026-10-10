//! Integration tests for multilevel adaptive continuation on real cut geometry.
use super::*;
use fs_cutfem::elastic3::{ElasticityError3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{
    AdaptiveMultilevelOptions3, AdaptivePreconditionError3,
};
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::multilevel::{MultilevelBudget, MultilevelError};
use crate::{SolveBudget, SolveProgress};

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
    -> Result<AdaptiveElasticity3, GoalRefinementError3> {
    let mut poll = |_| checkpoint();
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll,
    ).map_err(ElasticityError3::from)?;
    Ok(AdaptiveElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q,
    )?)
}
fn run(extra_level: bool, interrupt_setup: bool) -> AdaptiveContinuationReport3 {
    let mut tree = Octree3::uniform(1, 4, 1000).unwrap();
    let original = tree.leaves().clone();
    let op = build(&tree, &mut || ControlFlow::Continue(())).unwrap();
    let coarse_tree = Octree3::uniform(0, 4, 1000).unwrap();
    let coarse = build(&coarse_tree, &mut || ControlFlow::Continue(())).unwrap();
    let coarse_scales = coarse.scales().to_vec();
    let mut study = CutDensityStudy3::new(op, 0.15, SimpParams::default());
    let raw = vec![0.5; study.cells()];
    let y = |_: [f64; 3]| [0.0, -1.0, 0.0];
    let z = |_: [f64; 3]| [0.0, 0.0, -1.0];
    let loads = [
        GoalReferenceLoad3 { load: ReferenceLoad3::body(&y), weight: 0.3 },
        GoalReferenceLoad3 { load: ReferenceLoad3::body(&z), weight: 0.7 },
    ];
    let schedule = [(1.0, 1.0), (2.0, 2.0)].map(|(penal, beta)| SimpParams {
        penal, beta, ..Default::default()
    });
    let options = AdaptiveContinuationOptions3 {
        optimization: MultiLoadOcOptions {
            max_iterations: 1, change_tolerance: 0.0, ..Default::default()
        },
        enrichment: GoalRefinementOptions3 {
            preconditioner: GoalPreconditioner3::Multilevel {
                options: AdaptiveMultilevelOptions3 {
                    hierarchy: MultilevelBudget { max_coarsest_dofs: 24, ..Default::default() },
                    ..Default::default()
                },
            },
            ..Default::default()
        },
        max_marks: 1, ..Default::default()
    };
    let mut poll = |p: SolveProgress| {
        if interrupt_setup && p.work.preconditioner_galerkin_products > 0 {
            ControlFlow::Break(())
        } else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let levels = [&coarse];
    let mut published = Vec::new();
    let mut first_scales = Vec::new();
    let report = controlled_adaptive_sdf3_continuation_with_coarse_levels_observed(
        &mut study, &mut tree, if extra_level { &levels } else { &[] },
        &loads, &raw, &schedule, options, &mut control, build,
        |study, _, report| {
            let last = report.continuation.last.as_ref().unwrap();
            assert!(last.history.len() > 1);
            assert!(report.continuation.stages.last().unwrap().gradient_check.as_ref().unwrap().passed());
            if published.is_empty() { first_scales = study.operator().scales().to_vec(); }
            published.push((last.rho.clone(), last.displacements.clone()));
            Ok::<(), ()>(())
        },
    ).unwrap();
    assert_eq!(report.continuation.work, control.work());
    assert_eq!(coarse.scales(), coarse_scales);
    if !extra_level || interrupt_setup {
        assert_eq!(published.len(), 1);
        assert_eq!(tree.leaves(), &original);
        assert_eq!(study.operator().scales(), first_scales);
        let last = report.continuation.last.as_ref().unwrap();
        assert_eq!(last.rho, published[0].0);
        assert_eq!(last.displacements, published[0].1);
    } else {
        assert_eq!(published.len(), 2);
        assert_eq!(report.continuation.termination, ContinuationTermination::ScheduleComplete);
        assert!(report.refinement_error.is_none());
        assert_eq!(report.refinements.len(), 1);
        assert!(report.refinements[0].installed);
        assert!(tree.leaves().len() > original.len());
        assert!(report.continuation.work.preconditioner_galerkin_products > 0);
        assert_eq!(report.continuation.work.preconditioner_operator_applications, 0);
        assert_ne!(published[1].1[0], published[1].1[1]);
    }
    report
}

#[test]
fn correction_ladder_unblocks_refinement_and_replays_real_independent_fields() {
    let refused = run(false, false);
    assert!(matches!(refused.refinement_error,
        Some(AdaptiveContinuationError3::Goal(GoalRefinementError3::Preconditioner(
            AdaptivePreconditionError3::Hierarchy(MultilevelError::Budget(_)))))));
    let a = run(true, false);
    let b = run(true, false);
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
}

#[test]
fn interruption_in_recursive_setup_retains_the_accepted_stage_and_spent_products() {
    let report = run(true, true);
    assert_eq!(report.continuation.evaluation_stop, Some(EvaluationStop::Cancelled));
    assert!(report.continuation.work.preconditioner_galerkin_products > 0);
    assert_eq!(report.continuation.stages.len(), 1);
    assert!(report.refinements.is_empty());
}
