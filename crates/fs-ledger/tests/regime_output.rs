//! Exact ledger retention for final operating-envelope demotion receipts.

use fs_evidence::{Color, ValidityDomain};
use fs_ledger::{Ledger, LedgerError, REGIME_DEMOTION_RECEIPT_ARTIFACT_KIND, hash_bytes};
use fs_regime::{
    AxisViolationKind, ConsumedModelCard, EnvelopeCoverage, OutputClaimReceipt, RegimeViolation,
};

fn receipt(qoi: &str, distance: f64) -> OutputClaimReceipt {
    OutputClaimReceipt {
        qoi: qoi.to_string(),
        original_color: Color::Validated {
            regime: ValidityDomain::new().with_bound("re", 0.0, 100.0),
            dataset: "wind-tunnel/reference".to_string(),
        },
        effective_color: Color::Estimated {
            estimator: "regime-extrapolation".to_string(),
            dispersion: f64::INFINITY,
        },
        in_domain_color: None,
        out_of_domain_color: Some(Color::Estimated {
            estimator: "regime-extrapolation".to_string(),
            dispersion: f64::INFINITY,
        }),
        coverage: EnvelopeCoverage::FullyOutOfDomain,
        in_domain_points: Vec::new(),
        out_of_domain_points: vec!["cruise-hot".to_string()],
        model_cards: vec![ConsumedModelCard {
            name: "forced-convection".to_string(),
            version: "2.1.0".to_string(),
        }],
        violations: vec![RegimeViolation {
            point: "cruise-hot".to_string(),
            card: "forced-convection".to_string(),
            card_version: "2.1.0".to_string(),
            axis: "re".to_string(),
            observed: Some(125.0),
            lo: 0.0,
            hi: 100.0,
            kind: AxisViolationKind::Above,
            distance,
        }],
        override_acknowledgement: None,
    }
}

#[test]
fn demotion_receipt_round_trips_exact_bytes_and_both_identities() {
    let ledger = Ledger::open(":memory:").expect("ledger opens");
    let receipt = receipt("max-junction-temperature", 0.25);
    let canonical = receipt.to_canonical_json();

    let first = ledger
        .put_regime_demotion_receipt(&receipt)
        .expect("demotion receipt is retained");
    assert_eq!(first.receipt_id, receipt.content_id());
    assert_eq!(first.artifact.hash, hash_bytes(canonical.as_bytes()));
    assert_eq!(first.artifact.len, canonical.len() as u64);
    assert!(!first.artifact.deduped);

    let info = ledger
        .artifact_info(&first.artifact.hash)
        .expect("metadata read succeeds")
        .expect("artifact exists");
    assert_eq!(info.kind, REGIME_DEMOTION_RECEIPT_ARTIFACT_KIND);
    assert_eq!(info.meta, None);
    assert_eq!(
        ledger
            .read_regime_demotion_receipt(&first.artifact.hash, &receipt)
            .expect("typed exact read succeeds"),
        canonical.as_bytes()
    );

    let retry = ledger
        .put_regime_demotion_receipt(&receipt)
        .expect("exact response-loss retry dedupes");
    assert_eq!(retry.receipt_id, first.receipt_id);
    assert_eq!(retry.artifact.hash, first.artifact.hash);
    assert!(retry.artifact.deduped);
}

#[test]
fn read_refuses_receipt_substitution_and_wrong_artifact_kind() {
    let ledger = Ledger::open(":memory:").expect("ledger opens");
    let expected = receipt("max-junction-temperature", 0.25);
    let stored = ledger
        .put_regime_demotion_receipt(&expected)
        .expect("demotion receipt is retained");

    let substituted = receipt("max-junction-temperature", 0.5);
    assert!(matches!(
        ledger.read_regime_demotion_receipt(&stored.artifact.hash, &substituted),
        Err(LedgerError::Invalid {
            field,
            ..
        }) if field == "regime_output_receipt.expected"
    ));

    let unrelated = ledger
        .put_artifact("report", b"not a regime receipt", None)
        .expect("unrelated artifact stores");
    assert!(matches!(
        ledger.read_regime_demotion_receipt(&unrelated.hash, &expected),
        Err(LedgerError::Invalid { field, .. }) if field == "artifact.kind"
    ));
}

#[test]
fn fully_in_domain_receipt_is_not_written_to_the_demotion_stream() {
    let ledger = Ledger::open(":memory:").expect("ledger opens");
    let mut in_domain = receipt("pressure-drop", 0.0);
    in_domain.effective_color = in_domain.original_color.clone();
    in_domain.in_domain_color = Some(in_domain.original_color.clone());
    in_domain.out_of_domain_color = None;
    in_domain.coverage = EnvelopeCoverage::FullyInDomain;
    in_domain.in_domain_points = vec!["nominal".to_string()];
    in_domain.out_of_domain_points.clear();
    in_domain.violations.clear();

    assert!(matches!(
        ledger.put_regime_demotion_receipt(&in_domain),
        Err(LedgerError::Invalid { field, .. })
            if field == "regime_output_receipt.coverage"
    ));

    let forged = ledger
        .put_artifact(
            REGIME_DEMOTION_RECEIPT_ARTIFACT_KIND,
            in_domain.to_canonical_json().as_bytes(),
            None,
        )
        .expect("generic storage can model a bypass attempt");
    assert!(matches!(
        ledger.read_regime_demotion_receipt(&forged.hash, &in_domain),
        Err(LedgerError::Invalid { field, .. })
            if field == "regime_output_receipt.coverage"
    ));
}
