//! Battery for the addendum risk register (Part V, R1–R10). Covers register
//! completeness + ordering, per-risk field non-emptiness, lookup, the audit
//! (complete + instrumented counts + gap detection on a deliberately
//! incomplete slice), JSON well-formedness, and determinism.

use fs_govern::{
    InstrumentationReceipt, InstrumentationStatus, MAX_RECEIPT_AGE_DAYS, ReceiptError, Risk,
    RiskId, audit, audit_slice, receipt_identity, register, risk, to_json,
};

fn evidence(label: &[u8]) -> fs_govern::ContentHash {
    fs_blake3::hash_domain("frankensim.fs-govern.test-evidence.v1", label)
}

#[test]
fn register_has_all_ten_risks_in_order() {
    let reg = register();
    assert_eq!(reg.len(), 10);
    for (r, id) in reg.iter().zip(RiskId::ALL) {
        assert_eq!(r.id, id);
    }
    // codes are R1..R10.
    let codes: Vec<&str> = reg.iter().map(|r| r.id.code()).collect();
    assert_eq!(
        codes,
        vec!["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10"]
    );
}

#[test]
fn every_risk_has_a_metric_owner_and_mitigation() {
    for r in register() {
        assert!(!r.name.is_empty(), "{:?} name", r.id);
        assert!(!r.description.is_empty(), "{:?} description", r.id);
        assert!(!r.mitigation.is_empty(), "{:?} mitigation", r.id);
        assert!(!r.early_warning.is_empty(), "{:?} early_warning", r.id);
        assert!(!r.threshold.is_empty(), "{:?} threshold", r.id);
        assert!(!r.owner.is_empty(), "{:?} owner", r.id);
    }
}

#[test]
fn owners_are_real_addendum_bead_ids_or_governance() {
    for r in register() {
        assert!(
            r.owner.starts_with("frankensim-"),
            "{:?} owner should be a bead id, got {}",
            r.id,
            r.owner
        );
    }
}

#[test]
fn lookup_returns_the_right_risk() {
    let r3 = risk(RiskId::R3);
    assert_eq!(r3.id, RiskId::R3);
    assert_eq!(r3.name, "Stable entity identity");
    // R1 is the Proposal-9 estimator-constants risk.
    assert_eq!(risk(RiskId::R1).owner, "frankensim-epic-flywheel-lmp4.1");
}

#[test]
fn audit_separates_declaration_from_live_operation() {
    let a = audit(200);
    assert_eq!(a.total, 10);
    assert_eq!(
        a.declared, 10,
        "every risk must declare a metric and an owner"
    );
    assert!(a.declared_schema_ok(), "schema gaps: {:?}", a.schema_gaps);
    // Honest baseline: nothing is verified live, so the register is
    // OPERATIONALLY RED with all ten exact gaps listed (xpck.9 — the
    // former single ok() rendered this state as a false green).
    assert_eq!(a.verified_instrumented, 0);
    assert!(!a.operationally_managed());
    assert_eq!(a.operational_gaps.len(), 10);
    assert!(
        a.operational_gaps
            .iter()
            .all(|(_, s)| *s == InstrumentationStatus::Uninstrumented)
    );
}

#[test]
fn receipts_bind_subject_provenance_and_age() {
    let mk = |receipt| Risk {
        id: RiskId::R2,
        name: "y",
        description: "y",
        mitigation: "y",
        early_warning: "a metric",
        threshold: "y",
        owner: "frankensim-epic-x",
        receipt,
    };
    // A valid receipt for a different subject cannot be replayed to turn R2
    // green. Sealed fields make direct identity tampering unrepresentable.
    let replayed = mk(Some(
        InstrumentationReceipt::new(
            "R1",
            "grafana://kills",
            "ci/governance-audit",
            evidence(b"R1-live-feed"),
            190,
        )
        .unwrap(),
    ));
    let a = audit_slice(&[replayed], 200);
    assert!(!a.operationally_managed());
    assert!(matches!(
        a.operational_gaps[0].1,
        InstrumentationStatus::BadReceipt
    ));
    // A consistent, fresh receipt verifies.
    let good_receipt = InstrumentationReceipt::new(
        "R2",
        "grafana://kills",
        "ci/governance-audit",
        evidence(b"R2-live-feed"),
        190,
    )
    .unwrap();
    let good = mk(Some(good_receipt));
    let a = audit_slice(&[good], 200);
    assert!(a.operationally_managed(), "gaps: {:?}", a.operational_gaps);
    assert_eq!(a.verified_instrumented, 1);
    assert_eq!(good_receipt.dashboard(), "grafana://kills");
    assert_eq!(good_receipt.verifier(), "ci/governance-audit");
    assert_eq!(good_receipt.verified_day(), 190);
    assert_eq!(good_receipt.evidence_artifact(), evidence(b"R2-live-feed"));
    assert!(good_receipt.is_consistent_for("R2"));
    assert!(!good_receipt.is_consistent_for("R1"));
    let json = good_receipt.to_json();
    assert!(json.contains(&good_receipt.identity().to_hex()));
    assert!(json.contains(&good_receipt.evidence_artifact().to_hex()));
    assert!(json.contains("ci/governance-audit"));
    let escaped = InstrumentationReceipt::new(
        "R2",
        "grafana://kills\u{0001}\\\"",
        "ci/governance-audit\nsecondary",
        evidence(b"escaped-provenance"),
        190,
    )
    .unwrap()
    .to_json();
    assert!(escaped.contains("grafana://kills\\u0001\\\\\\\""));
    assert!(escaped.contains("ci/governance-audit\\nsecondary"));
    assert!(
        !escaped.contains('\u{0001}'),
        "JSON may not contain raw controls"
    );
    // The same receipt, long unverified, DEMOTES to stale (dead
    // dashboards cannot keep claiming coverage).
    let a = audit_slice(&[good], 190 + MAX_RECEIPT_AGE_DAYS + 1);
    assert!(!a.operationally_managed());
    assert!(matches!(
        a.operational_gaps[0].1,
        InstrumentationStatus::Stale { .. }
    ));
    // A receipt "verified" in the future is bad, not fresh.
    let a = audit_slice(&[good], 100);
    assert!(matches!(
        a.operational_gaps[0].1,
        InstrumentationStatus::BadReceipt
    ));
}

#[test]
fn receipt_identity_binds_every_semantic_field() {
    let artifact = evidence(b"feed-snapshot-a");
    let base = receipt_identity("R2", "grafana://kills", "ci/a", artifact, 190);
    assert_eq!(
        base,
        receipt_identity("R2", "grafana://kills", "ci/a", artifact, 190),
        "canonical identity must replay exactly"
    );
    assert_ne!(
        receipt_identity("ab", "c", "ci/a", artifact, 190),
        receipt_identity("a", "bc", "ci/a", artifact, 190),
        "length framing must separate adjacent string fields"
    );
    for changed in [
        receipt_identity("R3", "grafana://kills", "ci/a", artifact, 190),
        receipt_identity("R2", "grafana://other", "ci/a", artifact, 190),
        receipt_identity("R2", "grafana://kills", "ci/b", artifact, 190),
        receipt_identity(
            "R2",
            "grafana://kills",
            "ci/a",
            evidence(b"feed-snapshot-b"),
            190,
        ),
        receipt_identity("R2", "grafana://kills", "ci/a", artifact, 191),
    ] {
        assert_ne!(base, changed, "a semantic edit must change the identity");
    }
}

#[test]
fn receipt_constructor_rejects_missing_provenance() {
    let artifact = evidence(b"feed-snapshot");
    assert_eq!(
        InstrumentationReceipt::new(" ", "grafana://kills", "ci/a", artifact, 190,),
        Err(ReceiptError::EmptySubject)
    );
    assert_eq!(
        InstrumentationReceipt::new("R2", "\t", "ci/a", artifact, 190),
        Err(ReceiptError::EmptyDashboard)
    );
    assert_eq!(
        InstrumentationReceipt::new("R2", "grafana://kills", "", artifact, 190),
        Err(ReceiptError::EmptyVerifier)
    );
    assert_eq!(
        InstrumentationReceipt::new(
            "R2",
            "grafana://kills",
            "ci/a",
            fs_govern::ContentHash([0; 32]),
            190,
        ),
        Err(ReceiptError::EmptyEvidenceArtifact)
    );
}

#[test]
fn audit_detects_a_missing_metric_or_owner() {
    // a deliberately incomplete risk must be caught (the audit is not vacuous).
    let bad = [
        Risk {
            id: RiskId::R1,
            name: "x",
            description: "x",
            mitigation: "x",
            early_warning: "", // missing
            threshold: "x",
            owner: "", // missing
            receipt: None,
        },
        Risk {
            id: RiskId::R2,
            name: "y",
            description: "y",
            mitigation: "y",
            early_warning: "a metric",
            threshold: "y",
            owner: "frankensim-epic-x",
            receipt: None,
        },
    ];
    let a = audit_slice(&bad, 200);
    assert_eq!(a.total, 2);
    assert_eq!(a.declared, 1);
    assert_eq!(a.verified_instrumented, 0);
    assert!(!a.declared_schema_ok());
    // both a missing-metric and a missing-owner gap on R1.
    assert!(
        a.schema_gaps
            .iter()
            .any(|(id, why)| *id == RiskId::R1 && why.contains("metric"))
    );
    assert!(
        a.schema_gaps
            .iter()
            .any(|(id, why)| *id == RiskId::R1 && why.contains("owner"))
    );
}

#[test]
fn audit_scope_cannot_be_empty_or_repeat_one_risk() {
    let empty = audit_slice(&[], 200);
    assert_eq!(empty.total, 0);
    assert!(!empty.declared_schema_ok());
    assert!(!empty.operationally_managed());

    let base = risk(RiskId::R1);
    let receipt = InstrumentationReceipt::new(
        base.id.code(),
        "ledger://risk/r1",
        "fs-govern-test",
        evidence(b"R1-complete"),
        200,
    )
    .expect("complete receipt");
    let repeated = Risk {
        receipt: Some(receipt),
        ..(*base).clone()
    };
    let duplicate = audit_slice(&[repeated.clone(), repeated], 200);
    assert_eq!(duplicate.total, 2);
    assert!(!duplicate.declared_schema_ok());
    assert!(!duplicate.operationally_managed());
    assert!(
        duplicate
            .schema_gaps
            .contains(&(RiskId::R1, "duplicate risk id"))
    );
}

#[test]
fn json_is_well_formed_and_complete() {
    let j = to_json(200);
    assert!(j.starts_with('[') && j.ends_with(']'));
    // one object per risk.
    assert_eq!(j.matches("\"id\":\"R").count(), 10);
    for id in RiskId::ALL {
        assert!(
            j.contains(&format!("\"id\":\"{}\"", id.code())),
            "missing {}",
            id.code()
        );
    }
    // Owner bead ids and the UNAMBIGUOUS instrumentation status are
    // present; the former boolean "instrumented" flag is gone (xpck.9).
    assert!(j.contains("frankensim-epic-flywheel-lmp4.1"));
    assert!(j.contains("\"instrumentation\":\"uninstrumented\""));
    assert_eq!(j.matches("\"receipt\":null").count(), 10);
    assert!(!j.contains("\"instrumented\""));
    // no accidental double-commas between objects.
    assert!(!j.contains(",,"));
}

#[test]
fn register_json_and_audit_are_deterministic() {
    assert_eq!(to_json(200), to_json(200));
    assert_eq!(audit(200), audit(200));
}
