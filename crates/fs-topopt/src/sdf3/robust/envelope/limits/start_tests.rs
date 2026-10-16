use super::*;
use super::super::super::tests::fixture;
use crate::pipeline::LoadCase;
use crate::{SolveBudget, SolveProgress};
use std::ops::ControlFlow;

fn caps(study: &mut CutDensityStudy3, y: &[f64], z: &[f64], control: &mut SolveControl<'_>) -> ComplianceCaps3 {
    let raw = vec![0.6; study.cells()];
    let e = study.evaluate_robust(&raw, &[LoadCase { force: y, weight: 0.5 },
        LoadCase { force: z, weight: 0.5 }], ComplianceLimitedOptions3::default().robustness, control).unwrap();
    ComplianceCaps3::new(e.eroded.objective.case_compliances.iter().map(|c| 1.2*c).collect(), 0.03, 2).unwrap()
}

#[test]
fn restoration_finds_actual_robust_feasibility_without_altering_source_or_caps() {
    let (mut study, y, z) = fixture(2.0); let raw = vec![0.2; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let limits = caps(&mut study, &y, &z, &mut control); let forces = [&y[..], &z[..]];
    let options = ComplianceStartOptions3 { max_bisections: 8, ..Default::default() };
    let report = controlled_compliance_feasible_start(&mut study, &forces, &raw, &limits,
        ComplianceLimitedOptions3::default(), options, &mut control).unwrap();
    assert_eq!(report.stop, ComplianceStartStop3::Restored, "{report:?}");
    assert_eq!(report.evaluations, 10);
    let blend = report.blend.unwrap(); assert!(blend > 0.0 && blend < 1.0);
    let point = report.point.unwrap();
    assert!(limits.satisfied_by(&point.evaluation.physical.eroded.objective.case_compliances));
    assert!(point.evaluation.envelope.value <= 1.0-options.margin);
    for (&r, &start) in point.rho.iter().zip(&raw) { assert_eq!(r, start+blend*(1.0-start)); }
    let repeated = controlled_compliance_feasible_start(&mut study, &forces, &point.rho, &limits,
        ComplianceLimitedOptions3::default(), options, &mut control).unwrap();
    assert_eq!(repeated.stop, ComplianceStartStop3::IncomingSatisfied);
    assert_eq!(repeated.evaluations, 1); assert_eq!(repeated.blend, Some(0.0));
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn a_refused_full_endpoint_never_becomes_a_feasible_or_global_infeasibility_claim() {
    let (mut study, _, z) = fixture(1.0); let raw = vec![0.2; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let full = study.evaluate_robust(&vec![1.0; study.cells()],
        &[LoadCase { force: &z, weight: 1.0 }], ComplianceLimitedOptions3::default().robustness, &mut control).unwrap();
    let limits = ComplianceCaps3::new(vec![full.eroded.objective.compliance*0.5], 0.03, 1).unwrap();
    let report = controlled_compliance_feasible_start(&mut study, &[&z], &raw, &limits,
        ComplianceLimitedOptions3::default(), ComplianceStartOptions3::default(), &mut control).unwrap();
    assert_eq!(report.stop, ComplianceStartStop3::FullDensityRejected);
    assert!(report.point.is_none() && report.blend.is_none());
    assert!(report.full_density_maximum_ratio.unwrap() > 1.9); assert_eq!(report.evaluations, 2);
    // Strong smoothing can refuse a full endpoint even when its PHYSICAL caps pass.
    let limits = ComplianceCaps3::new(vec![full.eroded.objective.compliance*1.01; 2], 0.1, 2).unwrap();
    let report = controlled_compliance_feasible_start(&mut study, &[&z, &z], &raw, &limits,
        ComplianceLimitedOptions3::default(), ComplianceStartOptions3::default(), &mut control).unwrap();
    assert_eq!(report.stop, ComplianceStartStop3::FullDensityRejected);
    assert!(report.full_density_maximum_ratio.unwrap() < 1.0);
    assert!(report.full_density_envelope.unwrap() > 1.0);
}

#[test]
fn interrupted_bisection_retains_the_checked_full_endpoint_and_spent_work() {
    let (mut study, y, z) = fixture(1.0); let raw = vec![0.2; study.cells()];
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-compliance-start-bisect" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let limits = caps(&mut study, &y, &z, &mut control);
    let report = controlled_compliance_feasible_start(&mut study, &[&y, &z], &raw, &limits,
        ComplianceLimitedOptions3::default(), ComplianceStartOptions3::default(), &mut control).unwrap();
    assert_eq!(report.stop, ComplianceStartStop3::EvaluationStopped);
    assert!(matches!(report.error, Some(EvaluationStop::Cancelled))); assert_eq!(report.blend, Some(1.0));
    assert_eq!(report.evaluations, 2);
    assert!(limits.satisfied_by(&report.point.unwrap().evaluation.physical.eroded.objective.case_compliances));
    assert!(report.work.linear_iterations > 0);
}

#[test]
fn zero_bisections_and_invalid_policy_have_explicit_work_bounds() {
    let (mut study, y, z) = fixture(1.0); let raw = vec![0.2; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let limits = caps(&mut study, &y, &z, &mut control);
    let options = ComplianceStartOptions3 { max_bisections: 0, ..Default::default() };
    let report = controlled_compliance_feasible_start(&mut study, &[&y, &z], &raw, &limits,
        ComplianceLimitedOptions3::default(), options, &mut control).unwrap();
    assert_eq!(report.stop, ComplianceStartStop3::Restored); assert_eq!(report.blend, Some(1.0));
    assert_eq!(report.evaluations, 2);
    let work = control.work();
    for options in [ComplianceStartOptions3 { max_bisections: 65, ..options },
        ComplianceStartOptions3 { margin: 1.0, ..options }] {
        assert!(matches!(controlled_compliance_feasible_start(&mut study, &[&y, &z], &raw, &limits,
            ComplianceLimitedOptions3::default(), options, &mut control), Err(ProjectedAlError::Invalid(_))));
    }
    assert_eq!(control.work(), work);
}
