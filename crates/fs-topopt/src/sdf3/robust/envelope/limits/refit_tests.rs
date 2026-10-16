use super::*;
use crate::sdf3::robust::envelope::limits::refit::ComplianceCapRefitOptions3;
use crate::sdf3::robust::envelope::limits::start::ComplianceStartStop3;

fn policy() -> ComplianceCapRefitOptions3 {
    let mut options = ComplianceCapRefitOptions3::default();
    options.accepted_steps = 8;
    options.start.max_bisections = 8;
    options.design.density_floor = 1e-5;
    options
}
fn requirements(study: &mut CutDensityStudy3<AdaptiveElasticity3>, control: &mut SolveControl<'_>) -> ComplianceCaps3 {
    let forces: Vec<_> = laws().iter().map(|&law|
        study.operator().reference_load(law, || ControlFlow::Continue(())).unwrap()).collect();
    let loads: Vec<_> = forces.iter().map(|force| crate::pipeline::LoadCase { force, weight: 0.5 }).collect();
    let full = study.evaluate_robust(&vec![1.0; study.cells()], &loads, RobustOptions3::default(), control).unwrap();
    // Fixture requirements, fixed before any refit or optimization.
    ComplianceCaps3::new(full.eroded.objective.case_compliances.iter().map(|c| 3.0 * c).collect(), 0.02, 2).unwrap()
}

#[test]
fn selective_cap_refit_preserves_requirements_and_rechecks_real_fine_fields() {
    let (tree, _, source, rho) = fixture();
    let fine = tree.refined(&[*tree.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap();
    let scales = source.operator().scales().to_vec();
    let mut reference = CutDensityStudy3::new(build(&fine), 0.15, params());
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let caps = requirements(&mut reference, &mut control);
    let result = source.refit_compliance_caps(build(&fine), &rho, params(), &laws(), &caps, policy(), &mut control).unwrap();
    assert!(result.completed_and_rechecked());
    assert_eq!(result.caps, caps);
    assert!(result.study.cells() > source.cells());
    assert!(result.history.len() > 1);
    let best = result.best.as_ref().unwrap();
    assert!(caps.satisfied_by(&best.evaluation.physical.eroded.objective.case_compliances));
    assert!(best.evaluation.physical.dilated_volume_fraction <= result.history[0].dilated_volume);
    let replay = result.recheck.as_ref().unwrap().as_ref().unwrap();
    assert_eq!(best.rho, replay.rho);
    assert_eq!(best.evaluation.physical.eroded.objective.displacements,
        replay.evaluation.physical.eroded.objective.displacements);
    assert_eq!(result.work, control.work());
    assert_eq!(source.operator().scales(), scales);
    assert_eq!(source.params().beta, params().beta);
}

#[test]
fn cap_refit_inherits_tiny_controls_and_protected_regions_through_same_parents() {
    let (_, fine, source, mut rho) = fixture();
    rho[0] = 2e-5;
    let mut mask = vec![PhysicalRegion3::Design; source.cells()];
    mask[0] = PhysicalRegion3::Solid; mask[1] = PhysicalRegion3::Void;
    let source = source.with_physical_regions(mask.clone()).unwrap();
    let fine_op = build(&fine);
    let parents = AdaptiveTransfer3::new(source.operator(), &fine_op, 2_000_000,
        || ControlFlow::Continue(())).unwrap().parents().to_vec();
    let expected: Vec<_> = parents.iter().map(|&p| mask[p]).collect();
    let mut reference = CutDensityStudy3::new(build(&fine), 0.15, params())
        .with_physical_regions(expected.clone()).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let caps = requirements(&mut reference, &mut control);
    let result = source.refit_compliance_caps(fine_op, &rho, params(), &laws(), &caps,
        ComplianceCapRefitOptions3 { accepted_steps: 0, ..policy() }, &mut control).unwrap();
    assert!(result.completed_and_rechecked());
    assert_eq!(result.study.physical_regions().unwrap(), expected);
    for (i, &parent) in parents.iter().enumerate() { assert_eq!(result.inherited_rho[i], rho[parent]); }
    let physical = &result.best.as_ref().unwrap().evaluation.physical;
    for (i, region) in expected.iter().enumerate() {
        for field in [&physical.eroded.projected_rho, &physical.nominal_rho, &physical.dilated_rho] {
            match region {
                PhysicalRegion3::Solid => assert_eq!(field[i], 1.0),
                PhysicalRegion3::Void => assert_eq!(field[i], 0.0),
                PhysicalRegion3::Design => {}
            }
        }
    }
}

#[test]
fn late_refit_cancellation_preserves_a_solved_prefix_and_source_grid() {
    let (_, fine, source, rho) = fixture();
    let scales = source.operator().scales().to_vec();
    let mut reference = CutDensityStudy3::new(build(&fine), 0.15, params());
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let caps = requirements(&mut reference, &mut control);
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-compliance-limited-accepted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let result = source.refit_compliance_caps(build(&fine), &rho, params(), &laws(), &caps, policy(), &mut control).unwrap();
    assert!(result.outcome.as_ref().unwrap().is_err());
    assert_eq!(result.history.len(), 2);
    assert!(result.best.is_some());
    assert!(!result.completed_and_rechecked());
    assert!(result.recheck.is_none());
    assert_eq!(source.operator().scales(), scales);
}

#[test]
fn impossible_caps_and_budget_stops_never_become_successful_refits() {
    let (_, fine, source, rho) = fixture();
    let scales = source.operator().scales().to_vec();
    let caps = ComplianceCaps3::new(vec![1e-12; 2], 0.02, 2).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let result = source.refit_compliance_caps(build(&fine), &rho, params(), &laws(), &caps, policy(), &mut control).unwrap();
    assert_eq!(result.restoration.stop, ComplianceStartStop3::FullDensityRejected);
    assert!(result.best.is_none() && result.outcome.is_none());
    assert!(!result.completed_and_rechecked());
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    let result = source.refit_compliance_caps(build(&fine), &rho, params(), &laws(), &caps, policy(), &mut control).unwrap();
    assert_eq!(result.restoration.stop, ComplianceStartStop3::EvaluationStopped);
    assert_eq!(result.work.linear_iterations, 1);
    assert!(result.best.is_none());
    assert_eq!(source.operator().scales(), scales);
}


#[test]
fn initial_reference_solve_uses_the_same_cap_optimizer_and_fresh_fields() {
    use crate::sdf3::robust::envelope::limits::ComplianceLimitedDesign3;
    let (_, _, source, rho) = fixture();
    let caps = ComplianceCaps3::new(vec![1e6, 2e6], 0.02, 2).unwrap();
    let mut options = policy(); options.accepted_steps = 2;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let mut solved = source.solve_reference_compliance_caps(&rho, &laws(), &caps, options, &mut control).unwrap();
    assert!(solved.completed_and_rechecked());
    assert_eq!(solved.inherited_rho, rho);
    let saved = solved.best.as_ref().unwrap().clone();
    let forces: Vec<_> = laws().iter().map(|&law| solved.study.operator()
        .reference_load(law, || ControlFlow::Continue(())).unwrap()).collect();
    let borrowed: Vec<_> = forces.iter().map(Vec::as_slice).collect();
    let mut direct = ComplianceLimitedDesign3::new(&mut solved.study, &borrowed,
        &solved.restoration.point.as_ref().unwrap().rho, caps.clone(), options.design, &mut control).unwrap();
    direct.run(2).unwrap();
    assert_eq!(direct.history(), solved.history);
    assert_eq!(direct.best_feasible().unwrap().rho, saved.rho);
    assert_eq!(direct.best_feasible().unwrap().evaluation.physical.eroded.objective.displacements,
        saved.evaluation.physical.eroded.objective.displacements);
}

#[test]
fn initial_reference_solve_keeps_restoration_when_optimizer_is_cancelled() {
    let (_, _, source, rho) = fixture();
    let caps = ComplianceCaps3::new(vec![1e6, 2e6], 0.02, 2).unwrap();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-compliance-limited-accepted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let solved = source.solve_reference_compliance_caps(&rho, &laws(), &caps, policy(), &mut control).unwrap();
    assert!(!solved.completed_and_rechecked());
    assert!(solved.restoration.point.is_some() && solved.best.is_some());
    assert_eq!(solved.history.len(), 2);
    assert!(matches!(solved.outcome, Some(Err(_))));
    assert!(solved.recheck.is_none());
    assert_eq!(solved.work, control.work());
}

#[test]
fn initial_reference_admission_precedes_load_callbacks() {
    let (_, _, source, rho) = fixture();
    let forbidden = |_: [f64; 3]| -> [f64; 3] { panic!("invalid policy reached load assembly") };
    let laws = [ReferenceLoad3::body(&forbidden)];
    let caps = ComplianceCaps3::new(vec![1e6, 2e6], 0.02, 2).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(source.solve_reference_compliance_caps(&rho, &laws, &caps, policy(), &mut control).is_err());
    assert_eq!(control.work().linear_solves, 0);
}
