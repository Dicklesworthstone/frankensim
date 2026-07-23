//! G0/G3 doctrine and drift checks for certificate-regime routing.

use std::path::Path;

use fs_govern::{
    CERTIFICATE_REGIME_NO_CLAIM, CERTIFICATE_REGIME_ROUTER_BEAD, CERTIFICATE_REGIME_SCHEMA_VERSION,
    CERTIFICATE_REGIMES, CERTIFICATE_REPORT_BOUNDARY, CapabilityRef, CapabilityStatus,
    CertificateRegimeError, ClaimClass, EvidenceRegime, INTERVAL_ROLES, THERMAL_EXAMPLES,
    certificate_regime, certificate_regime_json, certificate_regime_markdown_table,
    certificate_regimes, validate_certificate_regimes,
};

const DOCTRINE_DOC: &str = include_str!("../../../docs/CERTIFICATE_REGIMES.md");
const TABLE_BEGIN: &str = "<!-- BEGIN CODE-DERIVED CERTIFICATE REGIME TABLE -->";
const TABLE_END: &str = "<!-- END CODE-DERIVED CERTIFICATE REGIME TABLE -->";

#[test]
fn closed_v1_table_is_total_ordered_and_non_widening() {
    let rows = certificate_regimes().expect("canonical doctrine validates");
    assert_eq!(rows.len(), ClaimClass::ALL.len());
    assert_eq!(rows.len(), 8);

    for (index, claim) in ClaimClass::ALL.into_iter().enumerate() {
        let row = certificate_regime(claim);
        assert_eq!(row, &rows[index]);
        assert_eq!(row.claim, claim);
        assert_eq!(row.evidence, claim.required_evidence());
        assert_eq!(row.id, format!("CR-{:02}", index + 1));
        assert!(!row.scope.trim().is_empty());
        assert!(!row.no_claim.trim().is_empty());
        assert!(!row.capabilities.is_empty());
    }

    assert_eq!(
        certificate_regime(ClaimClass::ExactLongChaoticTrajectory).evidence,
        EvidenceRegime::NoUsefulBound
    );
    assert_eq!(
        certificate_regime(ClaimClass::LongHorizonMeanLoad).evidence,
        EvidenceRegime::StatisticalObservableWithModelEvidence
    );
}

#[test]
fn schema_order_identity_and_evidence_mutants_fail_closed() {
    assert_eq!(
        validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION + 1, &CERTIFICATE_REGIMES),
        Err(CertificateRegimeError::SchemaVersion {
            found: CERTIFICATE_REGIME_SCHEMA_VERSION + 1
        })
    );
    assert_eq!(
        validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION, &CERTIFICATE_REGIMES[..7]),
        Err(CertificateRegimeError::RowCount { found: 7 })
    );

    let mut rows = CERTIFICATE_REGIMES;
    rows.swap(0, 1);
    assert_eq!(
        validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION, &rows),
        Err(CertificateRegimeError::RowOrder {
            index: 0,
            expected: ClaimClass::RootOrEventTime,
            found: ClaimClass::ShortHorizonReachability,
        })
    );

    let mut rows = CERTIFICATE_REGIMES;
    rows[0].id = "CR-00";
    assert_eq!(
        validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION, &rows),
        Err(CertificateRegimeError::RowId {
            index: 0,
            expected: "CR-01",
            found: "CR-00",
        })
    );

    let mut rows = CERTIFICATE_REGIMES;
    rows[1].claim = ClaimClass::RootOrEventTime;
    assert_eq!(
        validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION, &rows),
        Err(CertificateRegimeError::DuplicateClaim {
            claim: ClaimClass::RootOrEventTime,
        })
    );

    let mut rows = CERTIFICATE_REGIMES;
    rows[7].evidence = EvidenceRegime::IntervalRootOrTaylorEnclosure;
    assert_eq!(
        validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION, &rows),
        Err(CertificateRegimeError::EvidenceMismatch {
            claim: ClaimClass::ExactLongChaoticTrajectory,
            expected: EvidenceRegime::NoUsefulBound,
            found: EvidenceRegime::IntervalRootOrTaylorEnclosure,
        })
    );
}

#[test]
fn capability_maturity_cannot_imply_unimplemented_source_evidence() {
    static STAGED_WITH_LOCATOR: [CapabilityRef; 1] = [CapabilityRef {
        crate_name: "fs-ivl",
        capability: "validated-reachability-tube",
        status: CapabilityStatus::Staged,
        source_locator: Some("crates/fs-ivl/src/taylor.rs"),
    }];
    let mut rows = CERTIFICATE_REGIMES;
    rows[1].capabilities = &STAGED_WITH_LOCATOR;
    assert_eq!(
        validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION, &rows),
        Err(CertificateRegimeError::StagedCapabilityHasLocator {
            row: "CR-02",
            capability: "validated-reachability-tube",
        })
    );

    static LIVE_WITHOUT_LOCATOR: [CapabilityRef; 1] = [CapabilityRef {
        crate_name: "fs-ivl",
        capability: "interval-root-isolation",
        status: CapabilityStatus::Available,
        source_locator: None,
    }];
    let mut rows = CERTIFICATE_REGIMES;
    rows[0].capabilities = &LIVE_WITHOUT_LOCATOR;
    assert_eq!(
        validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION, &rows),
        Err(CertificateRegimeError::LiveCapabilityMissingLocator {
            row: "CR-01",
            capability: "interval-root-isolation",
        })
    );

    static ESCAPING_LOCATOR: [CapabilityRef; 1] = [CapabilityRef {
        crate_name: "fs-ivl",
        capability: "interval-root-isolation",
        status: CapabilityStatus::Available,
        source_locator: Some("../outside.rs"),
    }];
    let mut rows = CERTIFICATE_REGIMES;
    rows[0].capabilities = &ESCAPING_LOCATOR;
    assert_eq!(
        validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION, &rows),
        Err(CertificateRegimeError::InvalidSourceLocator {
            row: "CR-01",
            locator: "../outside.rs",
        })
    );
}

#[test]
fn named_crates_and_live_capability_locators_exist() {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let repository = manifest_dir
        .parent()
        .and_then(Path::parent)
        .expect("fs-govern lives under repository/crates");

    for row in certificate_regimes().expect("canonical doctrine validates") {
        for capability in row.capabilities {
            let crate_dir = repository.join("crates").join(capability.crate_name);
            assert!(
                crate_dir.is_dir(),
                "{} names missing workspace crate {}",
                row.id,
                capability.crate_name
            );
            match capability.source_locator {
                Some(locator) => assert!(
                    repository.join(locator).is_file(),
                    "{} capability {} names missing source locator {locator}",
                    row.id,
                    capability.capability
                ),
                None => assert_eq!(
                    capability.status,
                    CapabilityStatus::Staged,
                    "{} capability {} omitted a locator without being staged",
                    row.id,
                    capability.capability
                ),
            }
        }
    }
}

#[test]
fn machine_readable_catalog_is_deterministic_and_complete() {
    let first = certificate_regime_json().expect("canonical doctrine renders");
    let second = certificate_regime_json().expect("canonical doctrine rerenders");
    assert_eq!(first, second);
    assert!(first.starts_with("{\"schema_version\":1,\"authority\":"));
    assert!(first.contains(CERTIFICATE_REGIME_NO_CLAIM));
    assert!(first.contains(CERTIFICATE_REGIME_ROUTER_BEAD));

    for claim in ClaimClass::ALL {
        let row = certificate_regime(claim);
        assert!(first.contains(&format!("\"id\":\"{}\"", row.id)));
        assert!(first.contains(&format!("\"claim\":\"{}\"", claim.code())));
        assert!(first.contains(&format!("\"evidence\":\"{}\"", row.evidence.code())));
    }
    assert!(first.contains("\"status\":\"staged\",\"source_locator\":null"));
}

#[test]
fn documentation_table_is_exactly_code_derived() {
    let start = DOCTRINE_DOC
        .find(TABLE_BEGIN)
        .expect("doctrine table begin marker")
        + TABLE_BEGIN.len();
    let tail = &DOCTRINE_DOC[start..];
    let end = tail.find(TABLE_END).expect("doctrine table end marker");
    let documented = &tail[..end];
    let generated = certificate_regime_markdown_table().expect("table renders");
    assert_eq!(documented, format!("\n{generated}"));

    assert!(DOCTRINE_DOC.contains("Intervals still have strong roles in chaotic systems"));
    assert!(DOCTRINE_DOC.contains("Ordinary repeated interval propagation"));
    assert!(DOCTRINE_DOC.contains("Worked thermal routes"));
    assert!(DOCTRINE_DOC.contains("cannot mint this reliability claim"));
}

#[test]
fn interval_roles_and_thermal_examples_retain_their_boundaries() {
    assert_eq!(INTERVAL_ROLES.len(), 5);
    for (index, role) in INTERVAL_ROLES.iter().enumerate() {
        assert_eq!(role.id, format!("IR-{:02}", index + 1));
        assert!(!role.supports.trim().is_empty());
        assert!(!role.boundary.trim().is_empty());
    }
    assert!(
        INTERVAL_ROLES
            .iter()
            .any(|role| role.supports.contains("root isolation"))
    );
    assert!(
        INTERVAL_ROLES
            .iter()
            .any(|role| role.supports.contains("conservation"))
    );
    assert!(
        INTERVAL_ROLES
            .iter()
            .any(|role| role.boundary.contains("not a validated reachability tube"))
    );

    assert_eq!(THERMAL_EXAMPLES.len(), 2);
    assert!(THERMAL_EXAMPLES[0].claim.contains("steady-state"));
    assert!(THERMAL_EXAMPLES[0].refusal.contains("duty cycles"));
    assert!(THERMAL_EXAMPLES[1].claim.contains("probability"));
    assert!(THERMAL_EXAMPLES[1].route.contains("rare-event"));
    assert!(THERMAL_EXAMPLES[1].refusal.contains("cannot mint"));

    assert!(CERTIFICATE_REPORT_BOUNDARY.contains("valid result"));
    assert!(CERTIFICATE_REPORT_BOUNDARY.contains("never promoted"));
}
