//! CORE UX replication sealing battery (bead `frankensim-leapfrog-2026-program-i94v.7.5.5.2.5.1`, V.5.5.H5a).

use fs_blake3::hash_domain;
use fs_govern::core_ux_replication::{
    CoreUxReplicationError, CoreUxReplicationSealV1, ReplicationProtocolSpec,
    CORE_UX_REPLICATION_SEAL_SCHEMA_V1,
};

fn sample_spec() -> ReplicationProtocolSpec {
    let h1_hash = hash_domain("test.h1", b"h1-protocol-v1");
    let h2_hash = hash_domain("test.h2", b"h2-privacy-v1");
    let task_hash = hash_domain("test.tasks", b"task-hazard-catalog-v1");
    let data_hash = hash_domain("test.data", b"disjoint-data-root");
    let artifact_hash = hash_domain("test.artifacts", b"disjoint-artifact-root");

    ReplicationProtocolSpec {
        h1_protocol_root: h1_hash,
        h2_privacy_contract_root: h2_hash,
        product_claim_revision: "rev-2026-07-30.1".into(),
        task_hazard_catalog_root: task_hash,
        cohorts: vec![
            "first-time-user".into(),
            "domain-engineer".into(),
            "safety-audit".into(),
            "assistive-tech".into(),
        ],
        recruitment_source: "external-panel-v1".into(),
        facilitator_roster: vec!["facilitator-1".into(), "facilitator-2".into()],
        accommodations: vec!["screen-reader".into(), "high-contrast".into()],
        power_target: 0.80,
        precision_target: 0.05,
        uncertainty_policy: "bootstrap-bca-95".into(),
        multiplicity_policy: "benjamini-hochberg".into(),
        missingness_policy: "conservative-bounds".into(),
        stopping_rule: "fixed-sample-preregistered".into(),
        allowed_deviations: vec!["platform-reconnect".into()],
        analysis_principal_ids: vec!["analyst-a".into(), "analyst-b".into()],
        checker_principal_ids: vec!["checker-x".into()],
        disjoint_data_root: data_hash,
        disjoint_artifact_root: artifact_hash,
        disclosure_roster: vec!["lead-adjudicator".into(), "review-board".into()],
        privacy_retention_policy: "zero-raw-retention".into(),
        no_outcome_access_attestation: true,
    }
}

#[test]
fn core_ux_replication_sealing_succeeds_and_is_deterministic() {
    let spec1 = sample_spec();
    let spec2 = sample_spec();

    let seal1 = CoreUxReplicationSealV1::seal(spec1, 1_700_000_000).expect("seal 1 succeeds");
    let seal2 = CoreUxReplicationSealV1::seal(spec2, 1_700_000_000).expect("seal 2 succeeds");

    assert_eq!(seal1.schema_version, CORE_UX_REPLICATION_SEAL_SCHEMA_V1);
    assert_eq!(seal1.seal_digest, seal2.seal_digest);
    assert_eq!(seal1.sealed_at_timestamp_s, 1_700_000_000);
}

#[test]
fn information_barrier_authorizes_only_whitelisted_disclosure_recipients() {
    let spec = sample_spec();
    let seal = CoreUxReplicationSealV1::seal(spec, 1_700_000_000).expect("seal succeeds");

    // Authorized recipient gets capability grant
    let grant = seal
        .authorize_disclosure("lead-adjudicator")
        .expect("authorized recipient must get grant");
    assert_eq!(grant.recipient, "lead-adjudicator");
    assert_eq!(grant.seal_digest, seal.seal_digest);

    // Unauthorized recipient refuses with PrematureDisclosureAttempt
    let err = seal
        .authorize_disclosure("unauthorized-party")
        .expect_err("unauthorized recipient must be blocked by information barrier");
    assert_eq!(
        err,
        CoreUxReplicationError::PrematureDisclosureAttempt {
            target: "unauthorized-party".into()
        }
    );
}

#[test]
fn sealing_refuses_without_attestation() {
    let mut spec = sample_spec();
    spec.no_outcome_access_attestation = false;

    let err = CoreUxReplicationSealV1::seal(spec, 1_700_000_000).expect_err("must refuse without attestation");
    assert_eq!(err, CoreUxReplicationError::MissingAttestation);
}

#[test]
fn sealing_refuses_empty_cohorts() {
    let mut spec = sample_spec();
    spec.cohorts.clear();

    let err = CoreUxReplicationSealV1::seal(spec, 1_700_000_000).expect_err("must refuse empty cohorts");
    assert!(matches!(err, CoreUxReplicationError::InvalidCohorts { .. }));
}

#[test]
fn sealing_refuses_invalid_power_target() {
    let mut spec = sample_spec();
    spec.power_target = 0.0; // Invalid
    assert!(matches!(
        CoreUxReplicationSealV1::seal(spec.clone(), 1_700_000_000).unwrap_err(),
        CoreUxReplicationError::InvalidPowerTarget { .. }
    ));

    spec.power_target = 1.5; // Exceeds 1.0
    assert!(matches!(
        CoreUxReplicationSealV1::seal(spec, 1_700_000_000).unwrap_err(),
        CoreUxReplicationError::InvalidPowerTarget { .. }
    ));
}

#[test]
fn sealing_refuses_duplicate_principals() {
    let mut spec = sample_spec();
    spec.analysis_principal_ids = vec!["analyst-a".into(), "analyst-a".into()];

    let err = CoreUxReplicationSealV1::seal(spec, 1_700_000_000).expect_err("must refuse duplicate principal");
    assert_eq!(
        err,
        CoreUxReplicationError::DuplicatePrincipal {
            id: "analyst-a".into()
        }
    );
}

#[test]
fn max_ux_expert_replication_sealing_succeeds_and_is_deterministic() {
    let h1_hash = hash_domain("test.h1", b"h1-max-protocol-v1");
    let h2_hash = hash_domain("test.h2", b"h2-max-privacy-v1");
    let task_hash = hash_domain("test.tasks", b"task-max-hazard-catalog-v1");
    let data_hash = hash_domain("test.data", b"disjoint-max-data-root");
    let artifact_hash = hash_domain("test.artifacts", b"disjoint-max-artifact-root");

    let spec = fs_govern::core_ux_replication::MaxReplicationProtocolSpec {
        h1_protocol_root: h1_hash,
        h2_privacy_contract_root: h2_hash,
        product_claim_revision: "rev-2026-07-30.max.1".into(),
        domain_tcb_strata: vec!["numerical-kernel".into(), "evidence-graph".into()],
        task_hazard_catalog_root: task_hash,
        expert_cohorts: vec!["domain-expert".into(), "theorem-researcher".into()],
        expert_role_criteria: vec!["tenured-faculty".into(), "lead-verification-engineer".into()],
        recruitment_source: "expert-panel-v1".into(),
        conflict_check_roster: vec!["conflict-check-passed".into()],
        facilitator_roster: vec!["lead-expert-facilitator".into()],
        accommodations: vec!["braille-display".into()],
        power_target: 0.85,
        precision_target: 0.02,
        uncertainty_policy: "exact-finite-sample".into(),
        multiplicity_policy: "familywise-error-rate".into(),
        missingness_policy: "complete-case-analysis".into(),
        stopping_rule: "preregistered-interim-bounds".into(),
        allowed_deviations: vec![],
        non_widening_restrictions: vec!["no-theorem-widening".into()],
        analysis_principal_ids: vec!["expert-analyst-1".into()],
        checker_principal_ids: vec!["expert-checker-1".into()],
        disjoint_data_root: data_hash,
        disjoint_artifact_root: artifact_hash,
        disclosure_roster: vec!["max-adjudication-committee".into()],
        privacy_retention_policy: "zero-raw-expert-retention".into(),
        no_outcome_access_attestation: true,
    };

    let seal = fs_govern::core_ux_replication::MaxUxReplicationSealV1::seal(spec, 1_700_000_100)
        .expect("max seal succeeds");

    assert_eq!(
        seal.schema_version,
        fs_govern::core_ux_replication::MAX_UX_REPLICATION_SEAL_SCHEMA_V1
    );

    let grant = seal
        .authorize_disclosure("max-adjudication-committee")
        .expect("authorized recipient gets disclosure grant");
    assert_eq!(grant.recipient, "max-adjudication-committee");

    let err = seal
        .authorize_disclosure("unauthorized-expert")
        .expect_err("unauthorized recipient blocked");
    assert_eq!(
        err,
        CoreUxReplicationError::PrematureDisclosureAttempt {
            target: "unauthorized-expert".into()
        }
    );
}
