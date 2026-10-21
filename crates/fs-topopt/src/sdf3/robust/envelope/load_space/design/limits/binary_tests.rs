use super::*;
use crate::sdf3::robust::envelope::tests::fixture;
use crate::sdf3::PhysicalRegion3;
use crate::{EvaluationStop, SolveBudget, SolveProgress};
use std::ops::ControlFlow;

fn close(a: f64, b: f64) {
    assert!((a-b).abs() <= 1e-8*a.abs().max(b.abs()).max(1e-12), "{a} != {b}");
}
fn loose() -> LoadSpaceCap3 { LoadSpaceCap3::new(2, 1e9, 0.02).unwrap() }

#[test]
fn separately_feasible_binary_basis_cases_do_not_admit_an_unsafe_combination() {
    let (mut study, y, _) = fixture(2.0);
    let raw = vec![0.8; study.cells()]; let basis = [&y[..], &y[..]];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let options = BinaryComplianceOptions3::default();
    let first = study.evaluate_binary_load_space_cap(&raw, &basis, &loose(), options, &mut control).unwrap();
    let individual = first.physical.case_compliances[0];
    let cap = LoadSpaceCap3::new(2, 1.25*individual, 0.02).unwrap();
    let before = control.work().linear_solves;
    let checked = study.evaluate_binary_load_space_cap(&raw, &basis, &cap, options, &mut control).unwrap();
    assert_eq!(control.work().linear_solves-before, 3); // one filter, two independent equilibria
    assert!(checked.physical.feasible(1e-6));
    assert!(!checked.feasible(1e-6));
    close(checked.worst_compliance, 2.0*individual);
    close(checked.spectrum.maximum, 1.6);
    close(checked.worst_coefficients.iter().map(|q| q*q).sum(), 1.0);
    assert!(checked.witness_relative_defect < 1e-9);
}

#[test]
fn binary_worst_force_replays_its_own_material_with_protected_regions() {
    let (study, y, z) = fixture(3.0);
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid; regions[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(regions).unwrap();
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.35+0.05*i as f64).collect();
    let second: Vec<_> = y.iter().zip(&z).map(|(y,z)| 0.4*y+z).collect();
    let options = BinaryComplianceOptions3 { threshold: 0.45, void_stiffness: 0.005, ..Default::default() };
    let scales = study.operator().scales().to_vec();
    let params = study.params();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let checked = study.evaluate_binary_load_space_cap(&raw, &[&y, &second], &loose(), options, &mut control).unwrap();
    let direct = study.evaluate_binary_compliance_caps(&raw, &[&checked.worst_force],
        &ComplianceCaps3::new(vec![loose().cap()], 0.02, 1).unwrap(), options, &mut control).unwrap();
    close(direct.case_compliances[0], checked.worst_compliance);
    for (&a, &b) in direct.displacements[0].iter().zip(&checked.worst_displacement) {
        assert!((a-b).abs() <= 1e-8*checked.worst_displacement.iter().map(|x| x.abs()).fold(1.0, f64::max));
    }
    assert_eq!(direct.eroded, checked.physical.eroded);
    for mask in [&checked.physical.eroded, &checked.physical.nominal, &checked.physical.dilated] {
        assert!(mask.solid[0]); assert!(!mask.solid[1]);
    }
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(study.params().eta, params.eta); assert_eq!(study.params().beta, params.beta);
}

#[test]
fn binary_spectral_cap_is_invariant_to_rotation_of_the_same_load_set() {
    let (mut study, y, z) = fixture(2.0);
    let s = 0.5_f64.sqrt();
    let a: Vec<_> = y.iter().zip(&z).map(|(y,z)| s*(y+z)).collect();
    let b: Vec<_> = y.iter().zip(&z).map(|(y,z)| s*(z-y)).collect();
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.35+0.05*i as f64).collect();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let options = BinaryComplianceOptions3::default();
    let original = study.evaluate_binary_load_space_cap(&raw, &[&y,&z], &loose(), options, &mut control).unwrap();
    let rotated = study.evaluate_binary_load_space_cap(&raw, &[&a,&b], &loose(), options, &mut control).unwrap();
    close(original.spectrum.maximum, rotated.spectrum.maximum);
    close(original.spectrum.smooth_maximum, rotated.spectrum.smooth_maximum);
    close(original.worst_compliance, rotated.worst_compliance);
    assert_eq!(original.physical.eroded, rotated.physical.eroded);
    assert_eq!(original.physical.dilated, rotated.physical.dilated);
    assert_eq!(original.feasible(1e-6), rotated.feasible(1e-6));
}

#[test]
fn zero_load_set_and_rejected_inputs_do_not_invent_compliance() {
    let (mut study, y, z) = fixture(2.0);
    let raw = vec![0.6; study.cells()]; let zero = vec![0.0; y.len()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let options = BinaryComplianceOptions3::default();
    assert!(study.evaluate_binary_load_space_cap(&raw, &[&y], &loose(), options, &mut control).is_err());
    assert!(study.evaluate_binary_load_space_cap(&raw, &[&y,&z], &loose(),
        BinaryComplianceOptions3 { max_load_cases: 1, ..options }, &mut control).is_err());
    assert!(study.evaluate_binary_load_space_cap(&raw, &[&y,&z], &loose(),
        BinaryComplianceOptions3 { threshold: f64::NAN, ..options }, &mut control).is_err());
    assert_eq!(control.work().linear_solves, 0);
    let result = study.evaluate_binary_load_space_cap(&raw, &[&zero,&zero], &loose(), options, &mut control).unwrap();
    assert_eq!(result.worst_compliance, 0.0);
    assert_eq!(result.spectrum.maximum, 0.0);
    close(result.spectrum.smooth_maximum, loose().policy().bias_bound());
    assert!(result.feasible(0.0));
    assert!(!result.feasible(f64::NAN));
}

#[test]
fn cancellation_after_binary_solves_restores_the_source_and_publishes_nothing() {
    let (mut study, y, z) = fixture(2.0);
    let raw = vec![0.6; study.cells()]; let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-binary-load-space-spectrum" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let result = study.evaluate_binary_load_space_cap(&raw, &[&y,&z], &loose(), BinaryComplianceOptions3::default(), &mut control);
    assert!(matches!(result, Err(BinaryComplianceError3::Evaluation(EvaluationStop::Cancelled))));
    assert_eq!(control.work().linear_solves, 3);
    assert_eq!(study.operator().scales(), scales);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    let result = study.evaluate_binary_load_space_cap(&raw, &[&y,&z], &loose(), BinaryComplianceOptions3::default(), &mut control);
    assert!(matches!(result, Err(BinaryComplianceError3::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(control.work().linear_iterations, 1);
    assert_eq!(study.operator().scales(), scales);
}
