//! Adjudication test battery for RANS rung (bead `frankensim-extreal-program-f85xj.5.8.4`).

#![cfg(feature = "rans-rung")]

use fs_scenario::rans_adjudication::RansAdjudicationReceipt;
use fs_scenario::rans_card::RansCardDraft;
use fs_scenario::rans_validation::{
    RansValidationCase, RansValidationLedger, RansValidationStatus,
};

fn sample_card() -> fs_scenario::rans_card::RansModelCard {
    RansCardDraft::launder_sharma_channel("electronics-cooling/e10-rans")
        .freeze()
        .expect("canonical card freezes")
}

#[test]
fn adjudication_succeeds_for_canonical_card_and_ledger() {
    let card = sample_card();
    let ledger = RansValidationLedger::evaluate_canonical_matrix();
    let receipt = RansAdjudicationReceipt::adjudicate(&card, &ledger)
        .expect("canonical adjudication succeeds");

    assert_eq!(receipt.verdict, "ADMITTED_WITH_CONTEXTUAL_BOUNDS");
    assert_eq!(receipt.node.node_id, "e10-low-re-rans");
    assert_eq!(receipt.node.authority_class, "Estimate"); // Cost != Authority!
    assert!(receipt.edge.admitted_contexts.contains(&"heatsink_fin_array"));
    assert!(receipt.edge.refused_contexts.contains(&"massive_unsteady_separation"));
}

#[test]
fn adjudication_refuses_corrupted_validation_ledger() {
    let card = sample_card();
    let mut ledger = RansValidationLedger::evaluate_canonical_matrix();
    ledger.cases.push(RansValidationCase {
        case_id: "unexplained-failure",
        fidelity_tier: "Level-C",
        qoi_name: "thermal_resistance",
        expected_envelope: (0.1, 0.2),
        observed_value: 10.0,
        attribution: "unresolved",
        status: RansValidationStatus::Falsified,
    });

    assert!(RansAdjudicationReceipt::adjudicate(&card, &ledger).is_err());
}

#[test]
fn adjudication_receipt_hash_is_deterministic() {
    let card = sample_card();
    let ledger = RansValidationLedger::evaluate_canonical_matrix();
    let r1 = RansAdjudicationReceipt::adjudicate(&card, &ledger).unwrap();
    let r2 = RansAdjudicationReceipt::adjudicate(&card, &ledger).unwrap();
    assert_eq!(r1.content_hash().to_hex(), r2.content_hash().to_hex());
}
