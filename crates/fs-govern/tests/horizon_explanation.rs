//! Horizon trigger B battery (bead `frankensim-epic-addendum-xpck.5.6`):
//! boundary, mutant, honesty, and Rule-4 operator-mode tests for Proposal B
//! explanation objects.

use fs_govern::horizon_explanation::{
    evaluate_trigger_b, mint_trigger_b_receipt, ExplanationCase, ExplanationDisposition,
    OperatorMode, TriggerBReceipt, TriggerBRefusal, TriggerBVerdict,
    MAX_RECONCILIATION_FAILURE_RATE,
};

fn passing_case(id: &str) -> ExplanationCase {
    ExplanationCase {
        case_id: id.into(),
        observed_qoi: 100.0,
        attributed_channels: vec![40.0, 35.0, 25.0], // Sum = 100.0, residual = 0.0
        tolerance: 1e-4,
        honesty_gate_passed: true,
        narrative_emitted: true,
    }
}

fn failing_case(id: &str) -> ExplanationCase {
    ExplanationCase {
        case_id: id.into(),
        observed_qoi: 100.0,
        attributed_channels: vec![40.0, 30.0, 20.0], // Sum = 90.0, residual = 10.0 > tol
        tolerance: 1e-4,
        honesty_gate_passed: false, // Honestly refused!
        narrative_emitted: false,   // Narrative honestly suppressed!
    }
}

#[test]
fn rule_4_human_driven_mode_always_defers() {
    let battery = vec![passing_case("case-1"), passing_case("case-2")];
    assert_eq!(
        evaluate_trigger_b(OperatorMode::HumanDriven, &battery),
        Ok(TriggerBVerdict::Rule4Defer)
    );

    let receipt = mint_trigger_b_receipt(OperatorMode::HumanDriven, Some(&battery));
    assert_eq!(receipt.disposition, ExplanationDisposition::Rule4Defer);
    assert_eq!(receipt.verdict, TriggerBVerdict::Rule4Defer);
}

#[test]
fn agent_operator_mode_activates_when_failure_rate_is_at_or_below_10_percent() {
    // 1 failing out of 20 = 5% <= 10%: activates!
    let mut battery = Vec::new();
    for i in 1..=19 {
        battery.push(passing_case(&format!("pass-{i}")));
    }
    battery.push(failing_case("fail-20"));

    assert_eq!(
        evaluate_trigger_b(OperatorMode::AgentOperator, &battery),
        Ok(TriggerBVerdict::Activate)
    );

    let receipt = mint_trigger_b_receipt(OperatorMode::AgentOperator, Some(&battery));
    assert_eq!(receipt.disposition, ExplanationDisposition::Activate);
    assert_eq!(receipt.verdict, TriggerBVerdict::Activate);
    assert!((receipt.failure_rate - 0.05).abs() < 1e-9);
}

#[test]
fn boundary_exactly_at_10_percent_failure_rate_activates() {
    // 1 failing out of 10 = 10%: activates!
    let mut battery = Vec::new();
    for i in 1..=9 {
        battery.push(passing_case(&format!("pass-{i}")));
    }
    battery.push(failing_case("fail-10"));

    assert_eq!(
        evaluate_trigger_b(OperatorMode::AgentOperator, &battery),
        Ok(TriggerBVerdict::Activate)
    );
}

#[test]
fn boundary_above_10_percent_failure_rate_defers() {
    // 2 failing out of 10 = 20% > 10%: defers.
    let mut battery = Vec::new();
    for i in 1..=8 {
        battery.push(passing_case(&format!("pass-{i}")));
    }
    battery.push(failing_case("fail-9"));
    battery.push(failing_case("fail-10"));

    assert_eq!(
        evaluate_trigger_b(OperatorMode::AgentOperator, &battery),
        Ok(TriggerBVerdict::Defer)
    );

    let receipt = mint_trigger_b_receipt(OperatorMode::AgentOperator, Some(&battery));
    assert_eq!(receipt.disposition, ExplanationDisposition::Defer);
    assert_eq!(receipt.verdict, TriggerBVerdict::Defer);
}

#[test]
fn mutant_storytelling_over_refusal_is_refused() {
    let mut bad_case = failing_case("storyteller");
    bad_case.narrative_emitted = true; // VIOLATION: emitted narrative despite honesty failure!

    let battery = vec![passing_case("pass-1"), bad_case];
    assert!(matches!(
        evaluate_trigger_b(OperatorMode::AgentOperator, &battery),
        Err(TriggerBRefusal::NarrativeOverRefusal { .. })
    ));
}

#[test]
fn mutant_smearing_unreconciled_residual_is_refused() {
    let mut smeared_case = failing_case("smeared");
    smeared_case.honesty_gate_passed = true; // VIOLATION: passed gate despite high residual!

    let battery = vec![passing_case("pass-1"), smeared_case];
    assert!(matches!(
        evaluate_trigger_b(OperatorMode::AgentOperator, &battery),
        Err(TriggerBRefusal::UnreconciledPassedGate { .. })
    ));
}

#[test]
fn empty_battery_refuses() {
    let empty: Vec<ExplanationCase> = Vec::new();
    assert_eq!(
        evaluate_trigger_b(OperatorMode::AgentOperator, &empty),
        Err(TriggerBRefusal::EmptyBattery)
    );
}

#[test]
fn mint_receipt_returns_nodata_when_battery_absent() {
    let receipt: TriggerBReceipt = mint_trigger_b_receipt(OperatorMode::AgentOperator, None);
    assert_eq!(receipt.disposition, ExplanationDisposition::NoData);
    assert_eq!(receipt.verdict, TriggerBVerdict::Defer);
    assert!(receipt.failure_rate.is_nan());
}
