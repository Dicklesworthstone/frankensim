//! Privacy-policy conformance tests for disclosure, correlation, licensing,
//! export control, retention, and promotion decisions.

use fs_obs::privacy::{
    CorrelationDisclosure, CorrelationPolicy, Disclosure, ExportRealm, ExternalCorrelationMethod,
    ExternalCorrelationToken, FieldPolicy, LabeledField, LicenseRealm, RetentionClass,
    RetentionDisposition, Sensitivity, ShareAudience, ShareBlock, ShareError, ShareRequest,
    evaluate_share,
};
use fs_obs::{EvidenceCompleteness, EvidenceIntegrity, OperationalSupport, PromotionEffect};

fn policy(sensitivity: Sensitivity) -> FieldPolicy {
    FieldPolicy {
        sensitivity,
        license: LicenseRealm::Redistributable,
        export: ExportRealm::Unrestricted,
        retention: RetentionClass::Permanent,
    }
}

fn field(path: &str, value: &[u8], sensitivity: Sensitivity, required: bool) -> LabeledField {
    LabeledField::new(path, value, policy(sensitivity), required).expect("valid labeled field")
}

fn request(audience: ShareAudience) -> ShareRequest {
    ShareRequest::new(audience, 100, 64, 64 * 1_024).expect("bounded share request")
}

fn disclosed_bytes(disclosure: &Disclosure) -> Option<&[u8]> {
    match disclosure {
        Disclosure::Plain(bytes) => Some(bytes),
        Disclosure::Redacted { .. } | Disclosure::Omitted { .. } => None,
    }
}

#[test]
fn constructors_reject_unbounded_and_ambiguous_policy_state() {
    assert!(ShareRequest::new(ShareAudience::Public, 0, 0, 1).is_err());
    assert!(ShareRequest::new(ShareAudience::Public, 0, 1, 0).is_err());
    assert!(LabeledField::new("", b"value".to_vec(), policy(Sensitivity::Public), false).is_err());
    assert!(
        LabeledField::new(
            "field",
            b"value".to_vec(),
            FieldPolicy {
                sensitivity: Sensitivity::Public,
                license: LicenseRealm::Restricted("bad\nrealm".into()),
                export: ExportRealm::Unrestricted,
                retention: RetentionClass::Permanent,
            },
            false,
        )
        .is_err()
    );
    assert!(
        ExternalCorrelationToken::new(
            ExternalCorrelationMethod::Keyed {
                key_id: "key-1".into(),
            },
            "not hex",
        )
        .is_err()
    );
}

#[test]
fn public_manifest_never_contains_internal_personal_secret_or_credential_bytes() {
    let raw_internal = b"internal-host";
    let raw_personal = b"person@example.test";
    let raw_secret = b"launch-code";
    let raw_credential = b"Bearer top-secret-token";
    let manifest = evaluate_share(
        vec![
            field("a.public", b"publishable", Sensitivity::Public, false),
            field("b.internal", raw_internal, Sensitivity::Internal, false),
            field("c.personal", raw_personal, Sensitivity::PersonalData, true),
            field("d.secret", raw_secret, Sensitivity::Secret, true),
            field(
                "e.credential",
                raw_credential,
                Sensitivity::Credential,
                true,
            ),
        ],
        &request(ShareAudience::Public),
    )
    .expect("fixed-marker redaction is safe");

    assert_eq!(
        disclosed_bytes(&manifest.entries[0].disclosure),
        Some(b"publishable".as_slice())
    );
    let protected: [&[u8]; 4] = [
        raw_internal.as_slice(),
        raw_personal.as_slice(),
        raw_secret.as_slice(),
        raw_credential.as_slice(),
    ];
    for (entry, raw) in manifest.entries[1..].iter().zip(protected) {
        assert!(disclosed_bytes(&entry.disclosure).is_none());
        let Disclosure::Redacted {
            marker,
            correlation,
            ..
        } = &entry.disclosure
        else {
            panic!("sensitive field must retain a redaction row");
        };
        assert!(
            !marker
                .as_bytes()
                .windows(raw.len())
                .any(|window| window == raw)
        );
        assert!(correlation.is_none());
    }
    assert_eq!(manifest.completeness, EvidenceCompleteness::Incomplete);
    assert_eq!(manifest.integrity, EvidenceIntegrity::Intact);
    assert_eq!(manifest.promotion, PromotionEffect::Demoted);
    assert_eq!(manifest.support, OperationalSupport::Degraded);
}

#[test]
fn credentials_are_redacted_even_for_privileged_local_replay() {
    let local = request(ShareAudience::Local)
        .with_privilege(true)
        .with_personal_data(true)
        .with_secrets(true);
    let manifest = evaluate_share(
        vec![
            field("personal", b"name", Sensitivity::PersonalData, false),
            field("secret", b"configuration", Sensitivity::Secret, false),
            field("credential", b"api-key", Sensitivity::Credential, true),
        ],
        &local,
    )
    .expect("privileged local projection");
    assert_eq!(
        manifest
            .entries
            .iter()
            .map(|entry| entry.path.as_str())
            .collect::<Vec<_>>(),
        ["credential", "personal", "secret"],
        "share entries remain in deterministic path order",
    );
    assert!(matches!(
        &manifest.entries[0].disclosure,
        Disclosure::Redacted {
            reason: ShareBlock::Credential,
            ..
        }
    ));
    assert!(matches!(
        &manifest.entries[1].disclosure,
        Disclosure::Plain(_)
    ));
    assert!(matches!(
        &manifest.entries[2].disclosure,
        Disclosure::Plain(_)
    ));
    assert_eq!(manifest.completeness, EvidenceCompleteness::Incomplete);
}

#[test]
fn dictionary_testable_correlation_is_refused_for_personal_and_secret_data() {
    let personal = vec![field(
        "user.email",
        b"small-dictionary@example.test",
        Sensitivity::PersonalData,
        false,
    )];
    let error = evaluate_share(
        personal,
        &request(ShareAudience::Public).with_correlation(CorrelationPolicy::UnsaltedPublic),
    )
    .expect_err("unsalted PII correlation must refuse");
    assert!(matches!(
        error,
        ShareError::UnsafeCorrelation {
            sensitivity: Sensitivity::PersonalData,
            method: "unsalted",
            ..
        }
    ));

    let salted = ExternalCorrelationToken::new(
        ExternalCorrelationMethod::Salted {
            salt_id: "public-salt".into(),
        },
        "0123456789abcdef0123456789abcdef",
    )
    .expect("fixture token");
    let error = evaluate_share(
        vec![field("secret", b"one-of-four", Sensitivity::Secret, false)],
        &request(ShareAudience::Public).with_correlation(CorrelationPolicy::External(salted)),
    )
    .expect_err("salted low-entropy secret correlation must refuse");
    assert!(matches!(
        error,
        ShareError::UnsafeCorrelation {
            sensitivity: Sensitivity::Secret,
            method: "salted_dictionary_testable",
            ..
        }
    ));
}

#[test]
fn keyed_external_token_can_correlate_redacted_personal_data() {
    let keyed = ExternalCorrelationToken::new(
        ExternalCorrelationMethod::Keyed {
            key_id: "privacy-hmac-v3".into(),
        },
        "fedcba9876543210fedcba9876543210",
    )
    .expect("fixture keyed token");
    let manifest = evaluate_share(
        vec![field(
            "user.id",
            b"alice-123",
            Sensitivity::PersonalData,
            false,
        )],
        &request(ShareAudience::Public).with_correlation(CorrelationPolicy::External(keyed)),
    )
    .expect("keyed PII projection");
    assert!(matches!(
        &manifest.entries[0].disclosure,
        Disclosure::Redacted {
            correlation: Some(CorrelationDisclosure::ExternalKeyed { key_id, token }),
            ..
        } if key_id == "privacy-hmac-v3"
            && token == "fedcba9876543210fedcba9876543210"
    ));
}

#[test]
fn restricted_license_and_export_realms_fail_closed() {
    let restricted = LabeledField::new(
        "standard.excerpt",
        b"protected standard text".to_vec(),
        FieldPolicy {
            sensitivity: Sensitivity::Public,
            license: LicenseRealm::Restricted("std-license-2026".into()),
            export: ExportRealm::Controlled("ear99-review".into()),
            retention: RetentionClass::Permanent,
        },
        true,
    )
    .expect("valid restricted field");

    let public = request(ShareAudience::Public)
        .with_license_realms(["std-license-2026".to_string()])
        .expect("valid entitlement")
        .with_export_authorizations(["ear99-review".to_string()])
        .expect("valid export capability");
    let manifest = evaluate_share(vec![restricted.clone()], &public)
        .expect("public restriction is an explicit manifest row");
    assert!(matches!(
        &manifest.entries[0].disclosure,
        Disclosure::Omitted {
            reason: ShareBlock::License
        }
    ));
    assert_eq!(manifest.support, OperationalSupport::Unsupported);

    let organization = request(ShareAudience::Organization)
        .with_license_realms(["std-license-2026".to_string()])
        .expect("valid entitlement");
    let manifest = evaluate_share(vec![restricted.clone()], &organization)
        .expect("missing export authorization becomes an explicit omission");
    assert!(matches!(
        &manifest.entries[0].disclosure,
        Disclosure::Omitted {
            reason: ShareBlock::Export
        }
    ));

    let admitted = organization
        .with_export_authorizations(["ear99-review".to_string()])
        .expect("valid export capability");
    let manifest = evaluate_share(vec![restricted], &admitted).expect("authorized org share");
    assert!(matches!(
        &manifest.entries[0].disclosure,
        Disclosure::Plain(_)
    ));
}

#[test]
fn expiry_and_legal_hold_are_orthogonal_to_share_authority() {
    let expired = LabeledField::new(
        "expired",
        b"old".to_vec(),
        FieldPolicy {
            retention: RetentionClass::Until(99),
            ..policy(Sensitivity::Public)
        },
        true,
    )
    .expect("expired fixture");
    let held = LabeledField::new(
        "held",
        b"preserve".to_vec(),
        FieldPolicy {
            retention: RetentionClass::LegalHold("matter-42".into()),
            ..policy(Sensitivity::Public)
        },
        false,
    )
    .expect("hold fixture");
    let manifest = evaluate_share(vec![held, expired], &request(ShareAudience::Public))
        .expect("retention projection");
    assert_eq!(manifest.entries[0].path, "expired");
    assert!(matches!(
        &manifest.entries[0].disclosure,
        Disclosure::Omitted {
            reason: ShareBlock::Expired
        }
    ));
    assert_eq!(manifest.entries[0].retention, RetentionDisposition::Expired);
    assert_eq!(
        manifest.entries[1].retention,
        RetentionDisposition::LegalHold("matter-42".into())
    );
    assert!(matches!(
        &manifest.entries[1].disclosure,
        Disclosure::Plain(_)
    ));
}

#[test]
fn output_order_and_binary_handling_are_deterministic() {
    let binary = [0, 0xff, b'\n', b'\r', b'"', b'\\'];
    let fields = vec![
        field("z.last", &binary, Sensitivity::Public, false),
        field("a.first", b"first", Sensitivity::Public, false),
    ];
    let left = evaluate_share(fields.clone(), &request(ShareAudience::Public))
        .expect("first deterministic projection");
    let right = evaluate_share(
        fields.into_iter().rev().collect(),
        &request(ShareAudience::Public),
    )
    .expect("reordered deterministic projection");
    assert_eq!(left, right);
    assert_eq!(left.entries[0].path, "a.first");
    assert_eq!(
        disclosed_bytes(&left.entries[1].disclosure),
        Some(binary.as_slice())
    );
}

#[test]
fn duplicate_paths_and_input_budgets_refuse_before_disclosure() {
    let duplicate = vec![
        field("same", b"one", Sensitivity::Public, false),
        field("same", b"two", Sensitivity::Public, false),
    ];
    assert!(matches!(
        evaluate_share(duplicate, &request(ShareAudience::Public)),
        Err(ShareError::DuplicatePath { .. })
    ));

    let small = ShareRequest::new(ShareAudience::Public, 0, 1, 3).expect("small bound");
    assert!(matches!(
        evaluate_share(
            vec![field("large", b"four", Sensitivity::Public, false)],
            &small,
        ),
        Err(ShareError::BudgetExceeded {
            field: "input_bytes",
            ..
        })
    ));
}

#[test]
fn privilege_downgrade_changes_projection_without_mutating_labels() {
    let fields = vec![field(
        "user",
        b"private-person",
        Sensitivity::PersonalData,
        true,
    )];
    let local = request(ShareAudience::Local)
        .with_privilege(true)
        .with_personal_data(true);
    let admitted = evaluate_share(fields.clone(), &local).expect("local privileged replay");
    assert!(matches!(
        &admitted.entries[0].disclosure,
        Disclosure::Plain(_)
    ));
    assert_eq!(admitted.completeness, EvidenceCompleteness::Complete);

    let downgraded = evaluate_share(fields, &request(ShareAudience::Public))
        .expect("public downgrade produces redacted replay");
    assert!(matches!(
        &downgraded.entries[0].disclosure,
        Disclosure::Redacted { .. }
    ));
    assert_eq!(downgraded.completeness, EvidenceCompleteness::Incomplete);
    assert_eq!(downgraded.promotion, PromotionEffect::Demoted);
}
