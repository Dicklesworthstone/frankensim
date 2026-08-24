//! Validation and discrepancy test battery for RANS rung (bead `frankensim-extreal-program-f85xj.5.8.3`).

#![cfg(feature = "rans-rung")]

use fs_scenario::rans_validation::{
    RansValidationCase, RansValidationLedger, RansValidationStatus,
};

#[test]
fn canonical_matrix_contains_all_required_tiers() {
    let ledger = RansValidationLedger::evaluate_canonical_matrix();
    assert_eq!(ledger.cases.len(), 6);

    let tiers: Vec<&str> = ledger.cases.iter().map(|c| c.fidelity_tier).collect();
    assert!(tiers.contains(&"Level-A"));
    assert!(tiers.contains(&"Level-B"));
    assert!(tiers.contains(&"Level-C"));
    assert!(tiers.contains(&"LBM-Comparison"));
    assert!(tiers.contains(&"Adversarial-Stress"));

    // Validation passes with documented attributions
    assert!(ledger.validate().is_ok());
}

#[test]
fn content_hash_is_deterministic() {
    let l1 = RansValidationLedger::evaluate_canonical_matrix();
    let l2 = RansValidationLedger::evaluate_canonical_matrix();
    assert_eq!(l1.content_hash().to_hex(), l2.content_hash().to_hex());
    assert!(!l1.content_hash().to_hex().is_empty());
}

#[test]
fn unattributed_falsification_refuses_validation() {
    let mut ledger = RansValidationLedger::evaluate_canonical_matrix();
    ledger.cases.push(RansValidationCase {
        case_id: "unexplained-falsified-case",
        fidelity_tier: "Level-A",
        qoi_name: "unexpected_anomaly",
        expected_envelope: (1.0, 2.0),
        observed_value: 99.0,
        attribution: "falsified without resolution",
        status: RansValidationStatus::Falsified,
    });
    assert!(ledger.validate().is_err());
}

#[test]
fn missing_attribution_refuses_validation() {
    let mut ledger = RansValidationLedger::evaluate_canonical_matrix();
    ledger.cases[0].attribution = "";
    assert!(ledger.validate().is_err());
}
