use super::*;
use crate::projected::ProjectedSettings;
use crate::volume::VolumeProjectionSettings;
use crate::{Cantilever, GridSdf, OptimizeSettings};

fn policy() -> DwrGatePolicy {
    DwrGatePolicy { max_error_fraction: 0.5, max_solves_per_attempt: 18 }
}
fn measurement(coarse: f64, enriched: f64, eta: f64) -> DwrMeasurement {
    DwrMeasurement { snapshot: 1, level: 3, area: 0.6, coarse_compliance: coarse,
        enriched_compliance: enriched, eta_signed: eta,
        coarse_relative_residual: 1e-14, enriched_relative_residual: 1e-14 }
}
fn compare(b: DwrMeasurement, c: DwrMeasurement) -> DwrGainDecision {
    compare_dwr_gains(b, c, policy(), 0.0).unwrap()
}

#[test]
fn signed_dwr_gate_resolves_gains_not_just_finer_grid_descent() {
    let baseline = measurement(10.0, 11.0, 1.0);
    let resolved = measurement(9.0, 10.0, 1.0);
    assert_eq!(compare(baseline, resolved).refusal, None);
    assert_eq!(compare(baseline, baseline).refusal, Some("coarse-decrease-not-met"));
    assert_eq!(compare(baseline, measurement(9.0, 12.0, 3.0)).refusal,
        Some("enriched-decrease-not-met"));
    // Both grids improve, but the purported 1 J win loses 0.75 J on enrichment.
    assert_eq!(compare(baseline, measurement(9.0, 10.75, 1.75)).refusal,
        Some("discretization-dominates-gain"));
    // Isolate the signed DWR term: unchanged measured gains do NOT bypass it.
    assert_eq!(compare(baseline, measurement(9.0, 10.0, 2.0)).refusal,
        Some("discretization-dominates-gain"));
    // Signs matter: equal absolute corrections with opposite signs must refuse.
    assert_eq!(compare(baseline, measurement(9.0, 10.0, -1.0)).refusal,
        Some("discretization-dominates-gain"));
}

#[test]
fn signed_dwr_gate_has_exact_allowance_and_relative_decrease_boundaries() {
    let b = measurement(10.0, 10.0, 0.0);
    let edge = measurement(9.0, 9.5, 0.5);
    let at = compare(b, edge);
    assert_eq!(at.observed_gain_change, at.allowance);
    assert_eq!(at.signed_correction_change, at.allowance);
    assert!(at.refusal.is_none());
    let above = f64::from_bits(0.5f64.to_bits() + 1);
    assert_eq!(compare(b, measurement(9.0, 9.5, above)).refusal,
        Some("discretization-dominates-gain"));
    assert_eq!(compare_dwr_gains(b, edge, policy(), 0.1).unwrap().refusal,
        Some("coarse-decrease-not-met"));
    assert_eq!(compare_dwr_gains(b, edge, policy(), 0.05).unwrap().refusal,
        Some("enriched-decrease-not-met"));
}

#[test]
fn signed_dwr_gate_is_invariant_under_exact_quadratic_load_scaling() {
    for exponent in -20..=20 {
        let scale = 2.0f64.powi(2 * exponent);
        for eta in [0.0, 0.25, 0.5, 0.75, -0.75] {
            let b = measurement(16.0, 17.0, 1.0);
            let c = measurement(14.0, 15.0 + eta, 1.0 + eta);
            let scaled = |m: DwrMeasurement| DwrMeasurement {
                coarse_compliance: m.coarse_compliance * scale,
                enriched_compliance: m.enriched_compliance * scale,
                eta_signed: m.eta_signed * scale, ..m
            };
            let original = compare(b, c);
            let transformed = compare(scaled(b), scaled(c));
            assert_eq!(original.refusal, transformed.refusal);
            assert_eq!(original.coarse_gain * scale, transformed.coarse_gain);
            assert_eq!(original.allowance * scale, transformed.allowance);
        }
    }
}

#[test]
fn signed_dwr_gate_rejects_hostile_evidence_and_overflow() {
    let valid = measurement(10.0, 11.0, 1.0);
    for bad in [DwrMeasurement { level: 6, ..valid },
        DwrMeasurement { area: f64::NAN, ..valid },
        DwrMeasurement { area: 0.0, ..valid },
        DwrMeasurement { coarse_compliance: f64::INFINITY, ..valid },
        DwrMeasurement { enriched_compliance: -1.0, ..valid },
        DwrMeasurement { eta_signed: f64::NAN, ..valid },
        DwrMeasurement { coarse_relative_residual: 1e-12, ..valid },
        DwrMeasurement { enriched_relative_residual: -1e-15, ..valid }] {
        assert!(compare_dwr_gains(valid, bad, policy(), 0.0).is_err());
        assert!(compare_dwr_gains(bad, valid, policy(), 0.0).is_err());
    }
    assert!(compare_dwr_gains(valid, DwrMeasurement { level: 4, ..valid }, policy(), 0.0).is_err());
    assert!(compare_dwr_gains(valid, valid, policy(), f64::NAN).is_err());
    assert!(compare_dwr_gains(measurement(10.0, 11.0, f64::MAX),
        measurement(9.0, 10.0, -f64::MAX), policy(), 0.0).is_err());
}

#[test]
fn signed_dwr_gate_admits_the_whole_solve_family_at_n_minus_one_n_n_plus_one() {
    assert!(DwrGatePolicy { max_solves_per_attempt: 17, ..policy() }.validate(3, 8).is_err());
    assert!(policy().validate(3, 8).is_ok());
    assert!(DwrGatePolicy { max_solves_per_attempt: 19, ..policy() }.validate(3, 8).is_ok());
    for fraction in [0.0, 1.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(DwrGatePolicy { max_error_fraction: fraction, ..policy() }.validate(3, 8).is_err());
    }
    for (level, count) in [(0, 8), (6, 8), (3, 0), (3, 65), (3, usize::MAX)] {
        assert!(policy().validate(level, count).is_err());
    }
}

fn owner() -> ProjectedOptimizer {
    let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = phi.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
    let settings = OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6,
        move_cells: 0.1, nucleation_period: 0, ..OptimizeSettings::default() };
    ProjectedOptimizer::new(phi, Cantilever { load: 1.0, band: 0.125 }, settings, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        ProjectedSettings { max_candidates: 8, poll_iters: 1, ..ProjectedSettings::default() }).unwrap()
}
fn unchanged(a: &ProjectedOptimizer, b: &ProjectedOptimizer) {
    assert_eq!(a.checkpoint().geometry().nodes().iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        b.checkpoint().geometry().nodes().iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    assert_eq!(a.checkpoint().next_iteration(), b.checkpoint().next_iteration());
    assert_eq!(a.checkpoint().ell().to_bits(), b.checkpoint().ell().to_bits());
    assert_eq!(a.current().compliance.to_bits(), b.current().compliance.to_bits());
}

#[test]
fn signed_dwr_gate_cancels_real_baseline_and_candidate_assessments_without_publication() {
    let initial = owner();
    for target in [DwrCheckStage::Baseline(ComplianceDwrStage::BeforeEstimate),
        DwrCheckStage::Baseline(ComplianceDwrStage::BeforePublish),
        DwrCheckStage::Candidate { index: 0, stage: ComplianceDwrStage::BeforeEstimate },
        DwrCheckStage::Candidate { index: 0, stage: ComplianceDwrStage::BeforePublish }] {
        let mut state = initial.clone();
        // The first coarse-admissible candidate need not have index zero.
        let stopped = state.advance_one_dwr_controlled(policy(), |stage| {
            let hit = stage == target || matches!((stage, target),
                (DwrCheckStage::Candidate { stage: a, .. }, DwrCheckStage::Candidate { stage: b, .. }) if a == b);
            if hit { ControlFlow::Break(17) } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(stopped, ControlFlow::Break(17)), "stop was not reached: {target:?}");
        unchanged(&initial, &state);
    }
}

#[test]
fn signed_dwr_gate_real_replay_binds_measurements_and_preserves_failed_state() {
    let initial = owner();
    let mut a = initial.clone();
    let mut b = initial.clone();
    let run = |state: &mut ProjectedOptimizer|
        state.advance_one_dwr_controlled(policy(), |_| ControlFlow::<()>::Continue(())).unwrap();
    let first = run(&mut a);
    let replay = run(&mut b);
    assert_eq!(format!("{first:?}"), format!("{replay:?}"));
    unchanged(&a, &b);
    let ControlFlow::Continue(DwrCheckedProgress::Searched {
        baseline, candidates, reserved_solves, progress,
    }) = first else { panic!("real search must execute") };
    assert_eq!(baseline.snapshot, initial.current().snapshot);
    assert_eq!(baseline.coarse_compliance.to_bits(), initial.current().compliance.to_bits());
    assert!(!candidates.is_empty(), "must assess an actual coarse-admissible candidate");
    assert_eq!(reserved_solves, 2 * (1 + candidates.len()));
    assert!(reserved_solves <= policy().max_solves_per_attempt);
    match progress {
        ProjectedProgress::Accepted(step) => {
            let last = candidates.last().unwrap();
            assert!(last.refusal.is_none());
            assert!(last.decision.unwrap().refusal.is_none());
            assert_eq!(last.measurement.unwrap().snapshot, step.state.snapshot);
            assert_eq!(step.state.snapshot, a.current().snapshot);
        }
        ProjectedProgress::NoDescent(_) => unchanged(&initial, &a),
        ProjectedProgress::IterationLimit => panic!("initial state is not complete"),
    }
}

#[test]
fn signed_dwr_gate_underfunding_refuses_before_the_first_callback() {
    let mut state = owner();
    let initial = state.clone();
    let error = state.advance_one_dwr_controlled(
        DwrGatePolicy { max_solves_per_attempt: 17, ..policy() },
        |_| -> ControlFlow<()> { panic!("underfunded request started work") });
    assert!(error.is_err());
    unchanged(&initial, &state);
}
