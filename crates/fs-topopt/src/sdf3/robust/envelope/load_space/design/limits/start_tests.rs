use super::*;
use crate::sdf3::robust::envelope::tests::fixture;
use crate::sdf3::PhysicalRegion3;
use crate::{SolveBudget, SolveProgress};

fn cap_for<O: Sdf3Elasticity>(study: &mut CutDensityStudy3<O>, basis: &[&[f64]],
    control: &mut SolveControl<'_>) -> LoadSpaceCap3 {
    let e = study.evaluate_robust_load_space(&vec![1.0; study.cells()], basis,
        RobustOptions3::default(), &LoadSpace3::new(basis.len(), 1.0, 0.02).unwrap(), control).unwrap();
    LoadSpaceCap3::new(basis.len(), 2.0*e.worst_compliance, 0.02).unwrap()
}

#[test]
fn restored_combined_load_seed_replays_and_keeps_physical_regions() {
    let (study, y, z) = fixture(2.0);
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid; regions[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(regions).unwrap();
    let basis = [&y[..], &z[..]]; let raw = vec![0.15; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let cap = cap_for(&mut study, &basis, &mut control);
    let options = LoadSpaceStartOptions3 { max_bisections: 8, ..Default::default() };
    let result = controlled_load_space_feasible_start(&mut study, &basis, &raw, &cap,
        LoadSpaceLimitedOptions3::default(), options, &mut control).unwrap();
    assert_eq!(result.stop, LoadSpaceStartStop3::Restored);
    assert_eq!(result.evaluations, 10);
    let point = result.point.unwrap();
    assert!(cap.admits(&point.evaluation, 0.0));
    assert!(point.evaluation.spectrum.smooth_maximum <= 1.0-options.margin);
    assert!(result.blend.unwrap() > 0.0 && result.blend.unwrap() < 1.0);
    let fresh = study.evaluate_robust_load_space(&point.rho, &basis,
        RobustOptions3::default(), cap.policy(), &mut control).unwrap();
    assert_eq!(fresh.worst_displacement, point.evaluation.worst_displacement);
    assert_eq!(fresh.spectrum.gram, point.evaluation.spectrum.gram);
    for field in [&fresh.physical.eroded.projected_rho, &fresh.physical.nominal_rho, &fresh.physical.dilated_rho] {
        assert_eq!(field[0], 1.0); assert_eq!(field[1], 0.0);
    }
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn cancellation_in_bisection_retains_the_checked_full_endpoint() {
    let (mut study, y, z) = fixture(2.0);
    let basis = [&y[..], &z[..]]; let raw = vec![0.15; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let cap = cap_for(&mut study, &basis, &mut control);
    let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-load-space-cap-start-bisect" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let result = controlled_load_space_feasible_start(&mut study, &basis, &raw, &cap,
        LoadSpaceLimitedOptions3::default(), LoadSpaceStartOptions3::default(), &mut control).unwrap();
    assert_eq!(result.stop, LoadSpaceStartStop3::EvaluationStopped);
    assert_eq!(result.error, Some(EvaluationStop::Cancelled));
    assert_eq!(result.evaluations, 2); assert_eq!(result.blend, Some(1.0));
    let retained = result.point.unwrap();
    assert!(cap.admits(&retained.evaluation, 0.0));
    assert_eq!(retained.rho, vec![1.0; study.cells()]);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn full_refusal_invalid_policy_and_shared_budget_do_not_authorize_a_false_seed() {
    let (mut study, y, z) = fixture(2.0);
    let basis = [&y[..], &z[..]]; let raw = vec![0.15; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let cap = cap_for(&mut study, &basis, &mut control);
    let design = LoadSpaceLimitedOptions3::default();
    let before = control.work();
    assert!(controlled_load_space_feasible_start(&mut study, &basis, &raw, &cap, design,
        LoadSpaceStartOptions3 { max_bisections: 65, ..Default::default() }, &mut control).is_err());
    assert_eq!(control.work(), before);
    let refused = controlled_load_space_feasible_start(&mut study, &basis, &raw,
        &LoadSpaceCap3::new(2, 1e-20, 0.02).unwrap(), design, LoadSpaceStartOptions3::default(), &mut control).unwrap();
    assert_eq!(refused.stop, LoadSpaceStartStop3::FullDensityRejected);
    assert!(refused.point.is_none());
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let full = controlled_load_space_feasible_start(&mut study, &basis, &raw, &cap, design,
        LoadSpaceStartOptions3 { max_bisections: 0, ..Default::default() }, &mut control).unwrap();
    assert_eq!(full.stop, LoadSpaceStartStop3::Restored);
    let budget = full.work.linear_iterations+1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: budget, ..Default::default() }, &mut poll);
    let limited = controlled_load_space_feasible_start(&mut study, &basis, &raw, &cap, design,
        LoadSpaceStartOptions3::default(), &mut control).unwrap();
    assert_eq!(limited.stop, LoadSpaceStartStop3::EvaluationStopped);
    assert_eq!(limited.work.linear_iterations, budget);
    assert!(matches!(limited.error, Some(EvaluationStop::TotalBudget { .. })));
    assert_eq!(limited.point.unwrap().rho, full.point.unwrap().rho);
}
