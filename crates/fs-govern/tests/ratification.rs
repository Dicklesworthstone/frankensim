//! Battery for the vertical ratification decision record (bead f85xj.1.4):
//! the record validates against the measured world it cites, refuses to stand
//! on missing falsifiers / unbound kill criteria / placeholder baselines, and
//! the e2e governance listing logs every field of every decision record.

use fs_govern::{
    RatificationError, VERTICAL_RATIFICATION_V1, decision_records, ratification_json,
    ratified_vertical,
};
use fs_wedge::{
    CHT_BASELINE, HistoricalEvidenceTrustOrigin, RETIRED_PLACEHOLDER_BASELINE,
    default_recommendation, verify_default_comparison_evidence,
};

#[test]
fn the_ratification_record_stands_on_the_measured_world() {
    let record = ratified_vertical().expect("the ratification record must validate");
    assert_eq!(record.id, "frankensim-vertical-ratification-v1");
    assert_eq!(record.chosen_vertical, "thermal-design-assurance");
    assert_eq!(record.runner_up, "sdf-structural-topology-assurance");

    // The decision rests on the recomputable comparison, not on prose.
    let recommendation = default_recommendation().expect("comparison recomputes");
    assert_eq!(recommendation.recommended, record.chosen_vertical);
    assert_eq!(recommendation.runner_up, record.runner_up);
    let evidence = verify_default_comparison_evidence().expect("embedded evidence replays");
    assert_ne!(
        record.scoring_model_identity_blake3,
        "0000000000000000000000000000000000000000000000000000000000000000",
        "ratification must pin a non-placeholder comparison-model identity"
    );
    assert_eq!(
        record.scoring_model_identity_blake3,
        evidence.comparison_model_identity_blake3()
    );
    assert_ne!(
        record.scoring_evidence_descriptor_identity_blake3,
        "0000000000000000000000000000000000000000000000000000000000000000",
        "ratification must pin a non-placeholder evidence-descriptor identity"
    );
    assert_eq!(
        record.scoring_evidence_descriptor_identity_blake3,
        evidence.descriptor_identity_blake3()
    );
    assert_eq!(record.scoring_evidence_policy_version, 2);
    assert_eq!(
        record.scoring_evidence_policy_version,
        evidence.policy_version()
    );
    assert_eq!(
        record.scoring_evidence_trust_origin,
        HistoricalEvidenceTrustOrigin::EmbeddedDefault
    );
    assert_eq!(
        record.scoring_evidence_trust_origin,
        evidence.trust_origin()
    );
    assert!(!record.scoring_evidence_current_decision_authority);
    assert!(!record.scoring_evidence_human_review_authority);
    assert_eq!(
        record.scoring_evidence_current_decision_authority,
        evidence.current_decision_authority()
    );
    assert_eq!(
        record.scoring_evidence_human_review_authority,
        evidence.human_review_authority()
    );
    assert!(!evidence.current_decision_authority());
    assert!(!evidence.human_review_authority());

    // The kill criterion is bound to the measured baseline record.
    assert!((record.kill_target_reduction - CHT_BASELINE.target_reduction).abs() < 1e-12);
    assert_eq!(
        record.kill_within_quarters,
        CHT_BASELINE.kill_within_quarters
    );

    // Falsifiers are present and mechanically evaluable.
    assert!(record.falsifiers.len() >= 3);
    for falsifier in record.falsifiers {
        assert!(
            falsifier.is_complete(),
            "incomplete falsifier {}",
            falsifier.id
        );
    }

    // Downstream gates cite this record.
    assert!(
        record
            .downstream_gates
            .contains(&"frankensim-extreal-program-f85xj.6.1")
    );

    // The fs-wedge mirror carries the same conclusion and record id, so the
    // two crates cannot fork the story.
    assert_eq!(fs_wedge::RATIFIED_VERTICAL, record.chosen_vertical);
    assert_eq!(fs_wedge::RATIFICATION_RECORD_ID, record.id);
}

#[test]
fn a_record_without_falsifiers_is_refused() {
    let mut record = VERTICAL_RATIFICATION_V1;
    record.falsifiers = &[];
    assert_eq!(record.validate(), Err(RatificationError::MissingFalsifiers));
}

#[test]
fn an_incomplete_falsifier_is_refused() {
    use fs_govern::Falsifier;
    let mut record = VERTICAL_RATIFICATION_V1;
    record.falsifiers = &[Falsifier {
        id: "half-stated",
        statement: "something bad happens",
        measurement: "",
        threshold: "",
    }];
    assert_eq!(
        record.validate(),
        Err(RatificationError::IncompleteFalsifier { id: "half-stated" })
    );
}

#[test]
fn an_unbound_kill_criterion_is_refused() {
    let mut record = VERTICAL_RATIFICATION_V1;
    record.kill_target_reduction = 2.5;
    assert_eq!(
        record.validate(),
        Err(RatificationError::KillCriterionUnbound {
            field: "kill_target_reduction",
        })
    );

    let mut record = VERTICAL_RATIFICATION_V1;
    record.kill_within_quarters = 4;
    assert_eq!(
        record.validate(),
        Err(RatificationError::KillCriterionUnbound {
            field: "kill_within_quarters",
        })
    );
}

#[test]
fn a_placeholder_baseline_is_refused_as_kill_denominator() {
    let refusal = VERTICAL_RATIFICATION_V1.validate_against(&RETIRED_PLACEHOLDER_BASELINE);
    assert_eq!(
        refusal,
        Err(RatificationError::BaselineNotMeasured {
            provenance: "placeholder",
        })
    );
}

#[test]
fn empty_required_fields_are_refused() {
    let mut record = VERTICAL_RATIFICATION_V1;
    record.decided_on = "";
    assert_eq!(
        record.validate(),
        Err(RatificationError::EmptyField {
            field: "decided_on"
        })
    );

    let mut record = VERTICAL_RATIFICATION_V1;
    record.review_due = "  ";
    assert_eq!(
        record.validate(),
        Err(RatificationError::EmptyField {
            field: "review_due"
        })
    );

    let mut record = VERTICAL_RATIFICATION_V1;
    record.downstream_gates = &[];
    assert_eq!(
        record.validate(),
        Err(RatificationError::EmptyField {
            field: "downstream_gates",
        })
    );
}

#[test]
fn a_drifted_scoring_table_is_refused() {
    let mut record = VERTICAL_RATIFICATION_V1;
    record.chosen_vertical = "full-electronics-cooling-cht";
    match record.validate() {
        Err(RatificationError::ScoringDrift { field, .. }) => {
            assert_eq!(field, "chosen_vertical");
        }
        other => panic!("expected ScoringDrift, got {other:?}"),
    }

    let mut record = VERTICAL_RATIFICATION_V1;
    record.scoring_inventory_revision = "0000000000000000000000000000000000000000";
    match record.validate() {
        Err(RatificationError::ScoringDrift { field, .. }) => {
            assert_eq!(field, "scoring_inventory_revision");
        }
        other => panic!("expected ScoringDrift, got {other:?}"),
    }

    let mut record = VERTICAL_RATIFICATION_V1;
    record.scoring_model_identity_blake3 =
        "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff";
    match record.validate() {
        Err(RatificationError::ScoringDrift { field, .. }) => {
            assert_eq!(field, "scoring_model_identity_blake3");
        }
        other => panic!("expected model-identity ScoringDrift, got {other:?}"),
    }

    let mut record = VERTICAL_RATIFICATION_V1;
    record.scoring_evidence_descriptor_identity_blake3 =
        "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff";
    match record.validate() {
        Err(RatificationError::ScoringDrift { field, .. }) => {
            assert_eq!(field, "scoring_evidence_descriptor_identity_blake3");
        }
        other => panic!("expected descriptor-identity ScoringDrift, got {other:?}"),
    }

    let mut record = VERTICAL_RATIFICATION_V1;
    record.scoring_evidence_policy_version = 1;
    match record.validate() {
        Err(RatificationError::ScoringDrift { field, .. }) => {
            assert_eq!(field, "scoring_evidence_policy_version");
        }
        other => panic!("expected policy-version ScoringDrift, got {other:?}"),
    }

    let mut record = VERTICAL_RATIFICATION_V1;
    record.scoring_evidence_trust_origin =
        HistoricalEvidenceTrustOrigin::CallerSuppliedProtocolConsistency;
    match record.validate() {
        Err(RatificationError::ScoringDrift { field, .. }) => {
            assert_eq!(field, "scoring_evidence_trust_origin");
        }
        other => panic!("expected trust-origin ScoringDrift, got {other:?}"),
    }

    let mut record = VERTICAL_RATIFICATION_V1;
    record.scoring_evidence_current_decision_authority = true;
    match record.validate() {
        Err(RatificationError::ScoringDrift { field, .. }) => {
            assert_eq!(field, "scoring_evidence_current_decision_authority");
        }
        other => panic!("expected decision-authority ScoringDrift, got {other:?}"),
    }

    let mut record = VERTICAL_RATIFICATION_V1;
    record.scoring_evidence_human_review_authority = true;
    match record.validate() {
        Err(RatificationError::ScoringDrift { field, .. }) => {
            assert_eq!(field, "scoring_evidence_human_review_authority");
        }
        other => panic!("expected human-review-authority ScoringDrift, got {other:?}"),
    }
}

#[test]
fn governance_e2e_lists_every_decision_record_field_by_field() {
    // The e2e governance script: list all program-level decision records and
    // verify this one's completeness, with verbose field-by-field logging.
    let records = decision_records();
    assert_eq!(records.len(), 1, "exactly one ratified decision today");
    eprintln!("RESULT\tRECORD\tFIELD\tVALUE");
    for record in records {
        let mut log = |field: &str, value: &str, ok: bool| {
            eprintln!(
                "{}\t{}\t{}\t{}",
                if ok { "PASS" } else { "FAIL" },
                record.id,
                field,
                value.replace('\n', " ")
            );
            assert!(ok, "field {field} failed for {}", record.id);
        };
        log("id", record.id, !record.id.is_empty());
        log(
            "decided_on",
            record.decided_on,
            !record.decided_on.is_empty(),
        );
        log(
            "chosen_vertical",
            record.chosen_vertical,
            !record.chosen_vertical.is_empty(),
        );
        log("runner_up", record.runner_up, !record.runner_up.is_empty());
        log(
            "minority_report",
            record.minority_report,
            !record.minority_report.is_empty(),
        );
        log(
            "scoring_inventory_revision",
            record.scoring_inventory_revision,
            record.scoring_inventory_revision.len() == 40,
        );
        log(
            "scoring_model_identity_blake3",
            record.scoring_model_identity_blake3,
            record.scoring_model_identity_blake3.len() == 64,
        );
        log(
            "scoring_evidence_descriptor_identity_blake3",
            record.scoring_evidence_descriptor_identity_blake3,
            record.scoring_evidence_descriptor_identity_blake3.len() == 64,
        );
        log(
            "scoring_evidence_policy_version",
            &record.scoring_evidence_policy_version.to_string(),
            record.scoring_evidence_policy_version == 2,
        );
        log(
            "scoring_evidence_trust_origin",
            record.scoring_evidence_trust_origin.label(),
            record.scoring_evidence_trust_origin == HistoricalEvidenceTrustOrigin::EmbeddedDefault,
        );
        log(
            "scoring_evidence_current_decision_authority",
            if record.scoring_evidence_current_decision_authority {
                "true"
            } else {
                "false"
            },
            !record.scoring_evidence_current_decision_authority,
        );
        log(
            "scoring_evidence_human_review_authority",
            if record.scoring_evidence_human_review_authority {
                "true"
            } else {
                "false"
            },
            !record.scoring_evidence_human_review_authority,
        );
        for total in record.recorded_totals {
            log(
                "recorded_total",
                &format!("{}={}", total.candidate, total.weighted_total),
                total.weighted_total <= 1000,
            );
        }
        log(
            "kill_target_reduction",
            &record.kill_target_reduction.to_string(),
            record.kill_target_reduction > 1.0,
        );
        log(
            "kill_within_quarters",
            &record.kill_within_quarters.to_string(),
            record.kill_within_quarters > 0,
        );
        for falsifier in record.falsifiers {
            log("falsifier", falsifier.id, falsifier.is_complete());
        }
        log(
            "review_due",
            record.review_due,
            !record.review_due.is_empty(),
        );
        for gate in record.downstream_gates {
            log("downstream_gate", gate, !gate.is_empty());
        }
        record.validate().expect("record validates end to end");
    }
}

#[test]
fn the_ratification_json_is_fail_closed_and_deterministic() {
    let json = ratification_json().expect("valid records render");
    assert_eq!(json, ratification_json().expect("still valid"));
    assert!(json.starts_with("{\"decision_records\":["));
    assert!(json.contains("\"chosen_vertical\":\"thermal-design-assurance\""));
    assert!(
        json.contains(
            "\"scoring_inventory_revision\":\"b3b5f2c1c809eec06cde1e40cbc916d6995469b5\""
        )
    );
    assert!(!json.contains("\"declared_scoring_inventory_revision\""));
    assert!(json.contains("\"scoring_evidence_trust_origin\":\"embedded-default\""));
    assert!(json.contains("\"scoring_evidence_policy_version\":2"));
    assert!(json.contains("\"scoring_model_identity_blake3\":\""));
    assert!(json.contains("\"scoring_evidence_descriptor_identity_blake3\":\""));
    assert!(json.contains("\"scoring_evidence_current_decision_authority\":false"));
    assert!(json.contains("\"scoring_evidence_human_review_authority\":false"));
    assert!(json.contains("\"kill_target_reduction\":3"));
    assert!(json.contains("\"id\":\"level-c-thermal-data-unobtainable\""));
    assert!(json.contains("frankensim-extreal-program-f85xj.6.1"));
    assert!(json.ends_with("]}"));
}
