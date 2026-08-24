//! Horizon trigger D battery (bead `frankensim-epic-addendum-xpck.5.7`):
//! boundary, mutant, statistical distinguishability, and Rule-4 tests for
//! Proposal D Goodhart guard.

use fs_govern::horizon_goodhart::{
    evaluate_trigger_d, mint_trigger_d_receipt, EndpointStudy, EscalationStep,
    GoodhartDisposition, OperatorMode, ProposalDPremises, StepStatus, TriggerDReceipt,
    TriggerDRefusal, TriggerDVerdict,
};

fn all_available_steps() -> [(EscalationStep, StepStatus); 4] {
    [
        (EscalationStep::RungKPlus1, StepStatus::Available),
        (EscalationStep::CrossRepresentation, StepStatus::Available),
        (EscalationStep::DeltaPerturbation, StepStatus::Available),
        (EscalationStep::IndependentEstimator, StepStatus::Available),
    ]
}

fn distinguishable_study() -> EndpointStudy {
    EndpointStudy {
        preregistration_ref: "PREREG-2026-OPT-0012".into(),
        endpoint_sample_count: 100,
        endpoint_catches: 25, // 25% catch rate
        random_sample_count: 100,
        random_catches: 5,   // 5% catch rate
        p_value: 0.001,      // p < 0.05
    }
}

fn indistinguishable_study() -> EndpointStudy {
    EndpointStudy {
        preregistration_ref: "PREREG-2026-OPT-0013".into(),
        endpoint_sample_count: 100,
        endpoint_catches: 5,
        random_sample_count: 100,
        random_catches: 5,  // Same catch rate
        p_value: 0.50,      // Not statistically distinguishable
    }
}

#[test]
fn rule_4_human_driven_mode_always_defers() {
    let premises = ProposalDPremises {
        step_statuses: all_available_steps(),
        study: Some(distinguishable_study()),
    };
    assert_eq!(
        evaluate_trigger_d(OperatorMode::HumanDriven, &premises),
        Ok(TriggerDVerdict::Rule4Defer)
    );

    let receipt = mint_trigger_d_receipt(OperatorMode::HumanDriven, Some(&premises));
    assert_eq!(receipt.disposition, GoodhartDisposition::Rule4Defer);
    assert_eq!(receipt.verdict, TriggerDVerdict::Rule4Defer);
}

#[test]
fn agent_operator_mode_activates_when_steps_available_and_study_distinguishable() {
    let premises = ProposalDPremises {
        step_statuses: all_available_steps(),
        study: Some(distinguishable_study()),
    };
    assert_eq!(
        evaluate_trigger_d(OperatorMode::AgentOperator, &premises),
        Ok(TriggerDVerdict::Activate)
    );

    let receipt = mint_trigger_d_receipt(OperatorMode::AgentOperator, Some(&premises));
    assert_eq!(receipt.disposition, GoodhartDisposition::Activate);
    assert_eq!(receipt.verdict, TriggerDVerdict::Activate);
}

#[test]
fn indistinguishable_catch_rate_returns_budget_and_defers() {
    let premises = ProposalDPremises {
        step_statuses: all_available_steps(),
        study: Some(indistinguishable_study()),
    };
    assert_eq!(
        evaluate_trigger_d(OperatorMode::AgentOperator, &premises),
        Ok(TriggerDVerdict::IndistinguishableDefer)
    );

    let receipt = mint_trigger_d_receipt(OperatorMode::AgentOperator, Some(&premises));
    assert_eq!(receipt.disposition, GoodhartDisposition::Defer);
    assert_eq!(receipt.verdict, TriggerDVerdict::IndistinguishableDefer);
}

#[test]
fn unavailable_step_defers_provisionally() {
    let mut steps = all_available_steps();
    steps[1].1 = StepStatus::Unavailable; // Step unavailable
    let premises = ProposalDPremises {
        step_statuses: steps,
        study: Some(distinguishable_study()),
    };
    assert_eq!(
        evaluate_trigger_d(OperatorMode::AgentOperator, &premises),
        Ok(TriggerDVerdict::ProvisionalDefer)
    );
}

#[test]
fn missing_study_in_agent_mode_refuses() {
    let premises = ProposalDPremises {
        step_statuses: all_available_steps(),
        study: None,
    };
    assert_eq!(
        evaluate_trigger_d(OperatorMode::AgentOperator, &premises),
        Err(TriggerDRefusal::MissingStudy)
    );
}

#[test]
fn inadmissible_catches_greater_than_samples_refuses() {
    let mut study = distinguishable_study();
    study.endpoint_catches = 150; // 150 catches out of 100 samples!
    let premises = ProposalDPremises {
        step_statuses: all_available_steps(),
        study: Some(study),
    };
    assert!(matches!(
        evaluate_trigger_d(OperatorMode::AgentOperator, &premises),
        Err(TriggerDRefusal::InadmissibleStudyData { .. })
    ));
}

#[test]
fn mint_receipt_returns_nodata_when_premises_absent() {
    let receipt: TriggerDReceipt = mint_trigger_d_receipt(OperatorMode::AgentOperator, None);
    assert_eq!(receipt.disposition, GoodhartDisposition::NoData);
    assert_eq!(receipt.verdict, TriggerDVerdict::IndistinguishableDefer);
}
