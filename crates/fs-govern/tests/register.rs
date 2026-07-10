//! Battery for the addendum risk register (Part V, R1–R10). Covers register
//! completeness + ordering, per-risk field non-emptiness, lookup, the audit
//! (complete + instrumented counts + gap detection on a deliberately
//! incomplete slice), JSON well-formedness, and determinism.

use fs_govern::{
    InstrumentationReceipt, InstrumentationStatus, MAX_RECEIPT_AGE_DAYS, Risk, RiskId, audit,
    audit_slice, receipt_fingerprint, register, risk, to_json,
};

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
fn receipts_authenticate_and_age() {
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
    // A flipped flag without evidence — a receipt whose fingerprint was
    // never computed — CANNOT turn the audit green.
    let forged = mk(Some(InstrumentationReceipt {
        dashboard: "grafana://kills",
        verified_day: 190,
        fingerprint: 0xDEAD_BEEF,
    }));
    let a = audit_slice(&[forged], 200);
    assert!(!a.operationally_managed());
    assert!(matches!(
        a.operational_gaps[0].1,
        InstrumentationStatus::BadReceipt
    ));
    // A consistent, fresh receipt verifies.
    let good = mk(Some(InstrumentationReceipt {
        dashboard: "grafana://kills",
        verified_day: 190,
        fingerprint: receipt_fingerprint("R2", "grafana://kills", 190),
    }));
    let a = audit_slice(&[good], 200);
    assert!(a.operationally_managed(), "gaps: {:?}", a.operational_gaps);
    assert_eq!(a.verified_instrumented, 1);
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
    assert!(!j.contains("\"instrumented\""));
    // no accidental double-commas between objects.
    assert!(!j.contains(",,"));
}

#[test]
fn register_json_and_audit_are_deterministic() {
    assert_eq!(to_json(200), to_json(200));
    assert_eq!(audit(200), audit(200));
}
