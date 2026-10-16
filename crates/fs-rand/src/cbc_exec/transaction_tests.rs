//! Refused certificate/prefix/candidate transitions must remain resumable.
//! Compare the actual scan payload, not state_root(), which omits the incumbent.
use super::*;
use crate::cbc::CbcBudget;

fn executor(dimension: usize) -> CbcExecutor {
    let admission = CbcProblem::new(8, dimension).unwrap()
        .admit_for(CbcExecutionMode::Certified, CbcBudget::UNBOUNDED).unwrap();
    let mut value = CbcExecutor::new(admission).unwrap();
    value.enable_certificates().unwrap();
    value
}

fn tile() -> CbcTileShape {
    CbcTileShape::new(1, 1).unwrap()
}

fn state(value: &CbcExecutor) -> String {
    // Per-call transition counters reset on run entry; numerical state must not.
    format!("{:?}|{:?}|{:?}|{:?}|{}", value.phase, value.products,
        value.z, value.certificates, value.work_spent)
}

fn finish_and_compare(value: &mut CbcExecutor, dimension: usize) {
    let mut expected = executor(dimension);
    let mut keep_running = || CbcControl::Continue;
    assert_eq!(expected.run(&mut keep_running, tile(), u128::MAX).unwrap(),
        CbcRunStatus::Completed);
    assert_eq!(value.run(&mut keep_running, tile(), u128::MAX).unwrap(),
        CbcRunStatus::Completed);
    assert_eq!(state(value), state(&expected),
        "resumed fields, products, certificates and charged work must match uninterrupted execution");
}

fn begin_scan(dimension: usize) -> CbcExecutor {
    let mut value = executor(dimension);
    let allowance = 8 * value.schedule.initialization_visit_units()
        + value.schedule.prefix_control_units();
    assert_eq!(value.run(&mut || CbcControl::Continue, tile(), allowance).unwrap(),
        CbcRunStatus::AllowanceExhausted(CbcBoundary::Prefix));
    assert_eq!(value.prefix(), &[1]);
    value
}

fn complete_scan(dimension: usize) -> CbcExecutor {
    let mut value = begin_scan(dimension);
    // Seven candidates; four units modulo 8, each visiting all eight points.
    let allowance = 7 * (value.schedule.candidate_control_units()
        + value.schedule.certificate_candidate_units())
        + 4 * 8 * value.schedule.candidate_visit_units();
    assert_eq!(value.run(&mut || CbcControl::Continue, tile(), allowance).unwrap(),
        CbcRunStatus::AllowanceExhausted(CbcBoundary::CandidateBlock));
    let Phase::Scan { candidate: 8, accum: None, best: Some(_), runner_up, tie_class, .. }
        = &value.phase else { panic!("stop after the real scan, before certificate publication"); };
    assert!(runner_up.is_some(), "exercise the retained runner as well as the winner");
    assert!(tie_class.len() > 1, "exercise lowest-candidate tie preservation");
    value
}

#[test]
fn certificate_allowance_refusal_preserves_incumbent_and_retries_exactly_once() {
    let mut value = complete_scan(3);
    let units = value.schedule.certificate_prefix_units(2).unwrap();
    let before = state(&value);
    let work = value.work_spent();
    for allowance in [1, units - 1, 1] {
        let (status, receipt) = value.run_with_receipt(
            &mut || CbcControl::Continue, tile(), allowance).unwrap();
        assert_eq!(status, CbcRunStatus::NeedAllowance(CbcBoundary::CandidateBlock, units));
        assert_eq!(receipt.allowance_used, 0);
        assert_eq!(receipt.allowance_remaining, allowance);
        assert_eq!(receipt.committed_transitions, 0);
        assert_eq!(state(&value), before, "a replayable refusal cannot take the winning score");
    }
    assert_eq!(value.run(&mut || CbcControl::Continue, tile(), units).unwrap(),
        CbcRunStatus::AllowanceExhausted(CbcBoundary::CandidateBlock));
    assert_eq!(value.work_spent(), work + units);
    assert_eq!(value.certificates().len(), 1, "retry publishes exactly one certificate");
    assert!(matches!(value.phase, Phase::Update { next_point: 0, .. }));
    finish_and_compare(&mut value, 3);
}

#[test]
fn certificate_debit_refusal_does_not_consume_the_winner() {
    let mut value = complete_scan(3);
    let units = value.schedule.certificate_prefix_units(2).unwrap();
    let mut charges = TileCharges::sealed(&value, tile());
    // Exercise the existing debit's invariant-refusal branch without changing
    // the numerical state. This is not a newly admitted public work schedule.
    charges.admitted_work_units = value.work_spent();
    let before = state(&value);
    let transitions = value.run_transitions;
    let mut allowance = units;
    assert!(matches!(value.finish_scan(&charges, &mut allowance),
        Err(CbcExecError::ScheduleOverrun { .. })));
    assert_eq!(state(&value), before);
    assert_eq!(value.run_transitions, transitions);
    assert_eq!(allowance, units);
    finish_and_compare(&mut value, 3);
}

#[test]
fn prefix_storage_refusal_cannot_append_or_double_charge_a_component() {
    for update in [false, true] {
        let mut value = if update { complete_scan(3) } else { executor(3) };
        if update {
            let units = value.schedule.certificate_prefix_units(2).unwrap();
            assert_eq!(value.run(&mut || CbcControl::Continue, tile(), units).unwrap(),
                CbcRunStatus::AllowanceExhausted(CbcBoundary::CandidateBlock));
        }
        let pass_units = 8 * if update { value.schedule.product_update_visit_units() }
            else { value.schedule.initialization_visit_units() };
        let prefix_units = value.schedule.prefix_control_units();
        assert_eq!(value.run(&mut || CbcControl::Continue, tile(), pass_units).unwrap(),
            CbcRunStatus::NeedAllowance(CbcBoundary::Prefix, prefix_units));
        let before = state(&value);
        let prefix_len = value.prefix().len();
        let work = value.work_spent();
        let capacity = value.admissible_candidates_per_prefix;
        // Deterministic capacity-overflow refusal at the real fallible
        // reservation; no unsafe allocator replacement or resource exhaustion.
        value.admissible_candidates_per_prefix = usize::MAX;
        for _ in 0..2 {
            assert!(matches!(value.run(&mut || CbcControl::Continue, tile(), prefix_units),
                Err(CbcExecError::Storage(CbcStorageRefusal {
                    class: CbcStorageClass::CertificateTieScratch,
                    phase: CbcPhaseKind::Scan, ..
                }))));
            assert_eq!(state(&value), before, "failed next-phase allocation must not publish a prefix");
        }
        value.admissible_candidates_per_prefix = capacity;
        assert_eq!(value.run(&mut || CbcControl::Continue, tile(), prefix_units).unwrap(),
            CbcRunStatus::AllowanceExhausted(CbcBoundary::Prefix));
        assert_eq!(value.prefix().len(), prefix_len + 1);
        assert_eq!(value.work_spent(), work + prefix_units);
        finish_and_compare(&mut value, 3);
    }
}

#[test]
fn candidate_admission_refusal_preserves_the_tile_and_accumulator() {
    let mut value = begin_scan(3);
    let units = value.schedule.candidate_control_units()
        + value.schedule.certificate_candidate_units();
    let before = state(&value);
    for _ in 0..2 {
        assert_eq!(value.run(&mut || CbcControl::Continue, tile(), 1).unwrap(),
            CbcRunStatus::NeedAllowance(CbcBoundary::CandidateBlock, units));
        assert_eq!(state(&value), before, "a rejected admission cannot consume a tile slot");
    }
    let capacity = value.score_capacity_limbs;
    value.score_capacity_limbs = usize::MAX;
    assert!(matches!(value.run(&mut || CbcControl::Continue, tile(), units),
        Err(CbcExecError::Storage(CbcStorageRefusal {
            class: CbcStorageClass::ScoreAccumulator, ..
        }))));
    assert_eq!(state(&value), before);
    value.score_capacity_limbs = capacity;
    finish_and_compare(&mut value, 3);
}

#[test]
fn terminal_prefix_needs_no_next_scan_allocation() {
    let mut value = complete_scan(2);
    let units = value.schedule.certificate_prefix_units(2).unwrap();
    assert_eq!(value.run(&mut || CbcControl::Continue, tile(), units).unwrap(),
        CbcRunStatus::AllowanceExhausted(CbcBoundary::CandidateBlock));
    value.admissible_candidates_per_prefix = usize::MAX;
    // There is no following scan in dimension two; the sentinel must not be used.
    finish_and_compare(&mut value, 2);
}
