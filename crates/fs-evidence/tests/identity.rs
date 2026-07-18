//! G0/G3/G4 coverage for typed evidence identities.
//!
//! A retained cross-ISA known-answer vector is still required for G5.

use std::collections::BTreeMap;

use fs_blake3::identity::{
    CancellationProbe, CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field,
    FieldSpec, NoClaimState, SchemaId, SourceId, StrongIdentity, TrustState, WireType,
};
use fs_evidence::{
    Ambition, COLOR_ALGEBRA_VERSION, Certified, CertifiedF64EvidenceIdV1,
    CertifiedF64EvidenceIdentityError, CertifiedF64EvidenceReceiptV1, Color,
    ColorEvidenceCompositionOpV1, ColorEvidenceIdentityError, ColorEvidenceNodeIdV1,
    ColorEvidenceNodeIdentitySchemaV1, ColorEvidenceNodeKindV1, ColorEvidenceNodeV1,
    ColorEvidenceOperationV1, ColorEvidenceParentSemanticsV1, ColorEvidenceSourceIdV1,
    ColorEvidenceSourceV1, Evidence, IdentifiedCertifiedF64EvidenceV1, IdentifiedModelCardV1,
    IdentifiedValidityDomainV1, ModelCard, ModelCardCalibrationSourceIdV1,
    ModelCardCalibrationSourceReceiptV1, ModelCardIdV1, ModelCardIdentityError, ModelCardReceiptV1,
    ModelEvidence, NumericalKind, ProvenanceHash, SensitivitySummary, StatisticalCertificate,
    ValidityDomain, ValidityDomainIdV1, ValidityDomainIdentityError,
    compose_color_evidence_nodes_v1, identify_certified_f64_evidence_v1,
    identify_color_evidence_source_node_v1, identify_color_evidence_source_v1,
    identify_model_card_v1, identify_validity_domain_v1,
};

const LIMITS: CanonicalLimits = CanonicalLimits::new(16_384, 8_192, 32, 64, 256);

enum ForeignSourceSchema {}

impl CanonicalSchema for ForeignSourceSchema {
    const DOMAIN: &'static str = "org.frankensim.tests.foreign-color-source.v1";
    const NAME: &'static str = "foreign-color-source";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "negative child-binding fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

fn verified(lo: f64, hi: f64) -> Color {
    Color::Verified { lo, hi }
}

fn source_receipt(label: &str) -> ColorEvidenceSourceV1 {
    identify_color_evidence_source_v1(
        "org.frankensim.tests.color-source",
        1,
        label.as_bytes(),
        LIMITS,
        || false,
    )
    .expect("valid retained source")
}

fn source_node(label: &str, output: Color) -> ColorEvidenceNodeV1 {
    identify_color_evidence_source_node_v1(&source_receipt(label), output, LIMITS, || false)
        .expect("valid source node")
}

fn derived(
    operation: ColorEvidenceCompositionOpV1,
    left: &ColorEvidenceNodeV1,
    right: &ColorEvidenceNodeV1,
) -> Result<ColorEvidenceNodeV1, ColorEvidenceIdentityError> {
    compose_color_evidence_nodes_v1(operation, left, right, LIMITS, || false)
}

fn parent_reference_bytes(parent: ColorEvidenceNodeIdV1) -> [u8; 65] {
    let mut output = [0_u8; 65];
    output[0] = ColorEvidenceNodeIdV1::ROLE.tag();
    output[1..33]
        .copy_from_slice(SchemaId::<ColorEvidenceNodeIdentitySchemaV1>::for_schema().as_bytes());
    output[33..].copy_from_slice(parent.as_bytes());
    output
}

fn identified_domain(domain: ValidityDomain) -> IdentifiedValidityDomainV1 {
    identify_validity_domain_v1(domain, LIMITS, || false).expect("valid normalized domain")
}

fn certified_fixture() -> Certified<f64> {
    let model = ModelEvidence {
        cards: vec!["card-a".to_string(), "card-b".to_string()],
        assumptions: vec!["assumption-a".to_string(), "assumption-b".to_string()],
        validity: ValidityDomain::unconstrained()
            .with("Mach number", 0.2, 0.8)
            .with("α", -1.0, 1.0),
        discrepancy_rel: 0.125,
        in_domain: true,
    };
    let sensitivity = SensitivitySummary {
        d_qoi: BTreeMap::from([("Mach number".to_string(), -0.0), ("α".to_string(), 2.5)]),
    };
    Evidence::enclosed(2.0, 1.0, 3.0, ProvenanceHash(0x1111_2222_3333_4444))
        .with_statistical(StatisticalCertificate::EValue {
            e: 8.0,
            alpha: 0.05,
        })
        .with_model(model)
        .with_sensitivity(sensitivity)
        .with_adjoint(ProvenanceHash(0xaaaa_bbbb_cccc_dddd))
        .certified()
        .expect("valid certified scalar fixture")
}

fn identified_certified(certified: Certified<f64>) -> IdentifiedCertifiedF64EvidenceV1 {
    identify_certified_f64_evidence_v1(certified, LIMITS, || false)
        .expect("valid certified-f64 identity")
}

#[allow(
    clippy::too_many_lines,
    reason = "independent full-schema framing is the replay oracle for the helper"
)]
fn manual_certified_receipt(
    certified: &Certified<f64>,
    value: f64,
    in_domain: bool,
) -> CertifiedF64EvidenceReceiptV1 {
    let evidence = certified.evidence();
    let validity = identified_domain(evidence.model.validity.clone());
    let numerical_tag = match evidence.numerical.kind {
        NumericalKind::Exact => 1,
        NumericalKind::Enclosure => 2,
        NumericalKind::Estimate => 3,
        NumericalKind::NoClaim => 4,
    };
    let mut statistical_payload = [0_u8; 16];
    let (statistical_tag, statistical_len) = match evidence.statistical {
        StatisticalCertificate::None => (1, 0),
        StatisticalCertificate::EValue { e, alpha } => {
            statistical_payload[..8].copy_from_slice(&e.to_bits().to_le_bytes());
            statistical_payload[8..].copy_from_slice(&alpha.to_bits().to_le_bytes());
            (2, 16)
        }
        StatisticalCertificate::HalfWidth {
            half_width,
            confidence,
        } => {
            statistical_payload[..8].copy_from_slice(&half_width.to_bits().to_le_bytes());
            statistical_payload[8..].copy_from_slice(&confidence.to_bits().to_le_bytes());
            (3, 16)
        }
    };
    let sensitivity_rows: Vec<Vec<u8>> = evidence
        .sensitivity
        .d_qoi
        .iter()
        .map(|(parameter, derivative)| {
            let mut row = Vec::new();
            row.extend_from_slice(
                &u64::try_from(parameter.len())
                    .expect("parameter length")
                    .to_le_bytes(),
            );
            row.extend_from_slice(parameter.as_bytes());
            row.extend_from_slice(&derivative.to_bits().to_le_bytes());
            row
        })
        .collect();

    CanonicalEncoder::<CertifiedF64EvidenceIdV1, _>::new(LIMITS, || false)
        .expect("certified-f64 schema")
        .finite_f64(Field::new(0, "value"), value)
        .expect("value")
        .finite_f64(Field::new(1, "qoi"), evidence.qoi)
        .expect("qoi")
        .variant(Field::new(2, "numerical-kind"), numerical_tag, &[])
        .expect("numerical kind")
        .finite_f64(Field::new(3, "numerical-lo"), evidence.numerical.lo)
        .expect("numerical lo")
        .finite_f64(Field::new(4, "numerical-hi"), evidence.numerical.hi)
        .expect("numerical hi")
        .variant(
            Field::new(5, "statistical"),
            statistical_tag,
            &statistical_payload[..statistical_len],
        )
        .expect("statistical")
        .canonical_set(
            Field::new(6, "model-cards"),
            u64::try_from(evidence.model.cards.len()).expect("card count"),
            evidence.model.cards.iter().map(|card| card.as_bytes()),
        )
        .expect("model cards")
        .canonical_set(
            Field::new(7, "model-assumptions"),
            u64::try_from(evidence.model.assumptions.len()).expect("assumption count"),
            evidence
                .model
                .assumptions
                .iter()
                .map(|assumption| assumption.as_bytes()),
        )
        .expect("model assumptions")
        .child(Field::new(8, "model-validity"), validity.id())
        .expect("typed validity")
        .u64(
            Field::new(9, "model-discrepancy-ieee754-bits"),
            evidence.model.discrepancy_rel.to_bits(),
        )
        .expect("model discrepancy")
        .flag(Field::new(10, "model-in-domain"), in_domain)
        .expect("model in-domain")
        .ordered_bytes(
            Field::new(11, "sensitivity"),
            u64::try_from(sensitivity_rows.len()).expect("sensitivity count"),
            sensitivity_rows.iter().map(|row| row.as_slice()),
        )
        .expect("sensitivity")
        .flag(
            Field::new(12, "legacy-adjoint-correlation-present"),
            evidence.adjoint_ref.is_some(),
        )
        .expect("adjoint presence")
        .finish()
        .expect("manual certified-f64 identity")
}

const CALIBRATION_BYTES: &[u8] = b"calibration-artifact-v1\0binary";
const FNV1A64_ZERO_PREIMAGE: &[u8] = &[
    0x25, 0xe4, 0xe6, 0x90, 0x73, 0xfa, 0x7c, 0x26, 0x96, 0x1d, 0xcd, 0x31, 0x29, 0x0d, 0xe9, 0x72,
    0x17,
];
const FNV1A64_ZERO_EXTENSION_COLLISION: &[u8] = &[
    0x25, 0xe4, 0xe6, 0x90, 0x73, 0xfa, 0x7c, 0x26, 0x96, 0x1d, 0xcd, 0x31, 0x29, 0x0d, 0xe9, 0x72,
    0x17, 0x00,
];

fn model_card_fixture(calibration: Option<&[u8]>) -> ModelCard {
    let card = ModelCard::new(
        "les-α",
        "1.2.3+gpu",
        Ambition::Frontier,
        vec![
            "axis units declared elsewhere".to_string(),
            "continuum regime".to_string(),
        ],
        ValidityDomain::unconstrained()
            .with("Mach number", 0.2, 0.8)
            .with("α", -1.0, 1.0),
        vec![
            "high-angle separation".to_string(),
            "wall transition".to_string(),
        ],
        0.1,
    );
    match calibration {
        Some(bytes) => card.with_calibration(ProvenanceHash::of_bytes(bytes)),
        None => card,
    }
}

fn identified_model_card(
    card: ModelCard,
    calibration_bytes: Option<Vec<u8>>,
) -> IdentifiedModelCardV1 {
    identify_model_card_v1(card, calibration_bytes, LIMITS, || false)
        .expect("valid model-card identity")
}

fn manual_model_card_receipts(
    card: &ModelCard,
    calibration_present: bool,
    calibration_bytes: &[u8],
) -> (ModelCardCalibrationSourceReceiptV1, ModelCardReceiptV1) {
    let validity = identified_domain(card.validity.clone());
    let ambition_tag = match card.ambition {
        Ambition::Solid => 1,
        Ambition::Frontier => 2,
        Ambition::Moonshot => 3,
    };
    let calibration = CanonicalEncoder::<ModelCardCalibrationSourceIdV1, _>::new(LIMITS, || false)
        .expect("calibration source schema")
        .bytes(
            Field::new(0, "canonical-calibration-artifact"),
            calibration_bytes,
        )
        .expect("calibration bytes")
        .finish()
        .expect("calibration source receipt");
    let card_receipt = CanonicalEncoder::<ModelCardIdV1, _>::new(LIMITS, || false)
        .expect("model-card schema")
        .utf8(Field::new(0, "name"), &card.name)
        .expect("name")
        .utf8(Field::new(1, "version"), &card.version)
        .expect("version")
        .variant(Field::new(2, "ambition"), ambition_tag, &[])
        .expect("ambition")
        .canonical_set(
            Field::new(3, "assumptions"),
            u64::try_from(card.assumptions.len()).expect("assumption count"),
            card.assumptions.iter().map(|value| value.as_bytes()),
        )
        .expect("assumptions")
        .child(Field::new(4, "validity"), validity.id())
        .expect("validity")
        .canonical_set(
            Field::new(5, "known-failures"),
            u64::try_from(card.known_failures.len()).expect("known-failure count"),
            card.known_failures.iter().map(|value| value.as_bytes()),
        )
        .expect("known failures")
        .flag(Field::new(6, "calibration-present"), calibration_present)
        .expect("calibration presence")
        .child(Field::new(7, "calibration-source"), calibration.id())
        .expect("calibration source")
        .u64(
            Field::new(8, "discrepancy-rel-ieee754-bits"),
            card.discrepancy_rel.to_bits(),
        )
        .expect("discrepancy")
        .finish()
        .expect("manual model-card identity");
    (calibration, card_receipt)
}

#[test]
fn validity_domain_identity_normalizes_order_and_binds_every_bound_bit() {
    let first_domain = ValidityDomain::unconstrained()
        .with("z-axis", -3.0, 4.0)
        .with("a-axis", 0.0, 1.0);
    let replay_domain = ValidityDomain::unconstrained()
        .with("a-axis", 0.0, 1.0)
        .with("z-axis", -3.0, 4.0);
    let first = identified_domain(first_domain.clone());
    let replay = identified_domain(replay_domain);

    assert_eq!(first.id(), replay.id());
    assert_eq!(first.domain(), &first_domain);
    let audit = first.receipt().audit_record();
    assert_eq!(first.trust_state(), TrustState::Unanchored);
    assert_eq!(audit.trust(), TrustState::Unanchored);
    assert_eq!(audit.no_claim(), NoClaimState::ExternalTrustRequired);
    assert_eq!(audit.id(), first.id_bytes());

    let renamed = identified_domain(
        ValidityDomain::unconstrained()
            .with("a-axis-renamed", 0.0, 1.0)
            .with("z-axis", -3.0, 4.0),
    );
    let rebound = identified_domain(
        ValidityDomain::unconstrained()
            .with("a-axis", 0.0, 1.0_f64.next_up())
            .with("z-axis", -3.0, 4.0),
    );
    assert_ne!(first.id(), renamed.id());
    assert_ne!(first.id(), rebound.id());

    let negative_zero = identified_domain(ValidityDomain::unconstrained().with("x", -0.0, 0.0));
    let positive_zero = identified_domain(ValidityDomain::unconstrained().with("x", 0.0, 0.0));
    assert_ne!(negative_zero.id(), positive_zero.id());

    let unconstrained = identified_domain(ValidityDomain::unconstrained());
    assert!(unconstrained.domain().bounds().is_empty());
}

#[test]
fn validity_domain_identity_binds_arbitrary_utf8_axis_bytes() {
    let spaced = identified_domain(
        ValidityDomain::unconstrained()
            .with("Mach number", 0.0, 1.0)
            .with("α", -1.0, 1.0),
    );
    let renamed = identified_domain(
        ValidityDomain::unconstrained()
            .with("Mach-number", 0.0, 1.0)
            .with("alpha", -1.0, 1.0),
    );

    assert!(spaced.domain().bounds().contains_key("Mach number"));
    assert!(spaced.domain().bounds().contains_key("α"));
    assert_ne!(spaced.id(), renamed.id());
}

#[test]
fn validity_domain_identity_refuses_invalid_bounds_resources_and_cancellation() {
    let non_finite = identify_validity_domain_v1(
        ValidityDomain::unconstrained().with("x", 0.0, f64::INFINITY),
        LIMITS,
        || false,
    );
    assert_eq!(
        non_finite,
        Err(ValidityDomainIdentityError::InvalidBounds {
            axis_index: 0,
            reason: "bounds must be finite",
        })
    );

    let left = ValidityDomain::unconstrained().with("x", 0.0, 1.0);
    let right = ValidityDomain::unconstrained().with("x", 2.0, 3.0);
    let disjoint = identify_validity_domain_v1(left.intersect(&right), LIMITS, || false);
    assert_eq!(
        disjoint,
        Err(ValidityDomainIdentityError::InvalidBounds {
            axis_index: 0,
            reason: "lower bound exceeds upper bound",
        })
    );

    let three_axes = ValidityDomain::unconstrained()
        .with("a", 0.0, 1.0)
        .with("b", 0.0, 1.0)
        .with("c", 0.0, 1.0);
    let collection_tiny = CanonicalLimits::new(4096, 512, 32, 2, 64);
    let bounded = identify_validity_domain_v1(three_axes, collection_tiny, || false);
    assert!(matches!(
        bounded,
        Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CollectionItems,
                requested: 3,
                limit: 2,
            }
        ))
    ));

    let field_tiny = CanonicalLimits::new(4096, 43, 32, 8, 64);
    let field_bounded = identify_validity_domain_v1(
        ValidityDomain::unconstrained().with("axis", 0.0, 1.0),
        field_tiny,
        || false,
    );
    assert!(matches!(
        field_bounded,
        Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 44,
                limit: 43,
            }
        ))
    ));

    let aggregate_tiny = CanonicalLimits::new(4096, 73, 32, 8, 64);
    let aggregate_bounded = identify_validity_domain_v1(
        ValidityDomain::unconstrained()
            .with("a", 0.0, 1.0)
            .with("b", 0.0, 1.0),
        aggregate_tiny,
        || false,
    );
    assert!(matches!(
        aggregate_bounded,
        Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 74,
                limit: 73,
            }
        ))
    ));

    let mut sixteen_axes = ValidityDomain::unconstrained();
    for axis_index in 0..16 {
        sixteen_axes = sixteen_axes.with(format!("axis-{axis_index:02}"), 0.0, 1.0);
    }
    identified_domain(sixteen_axes.clone());
    let chunk_bounded =
        identify_validity_domain_v1(sixteen_axes.with("axis-16", 0.0, 1.0), LIMITS, || false);
    assert!(matches!(
        chunk_bounded,
        Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::StreamChunks,
                requested: 68,
                limit: 64,
            }
        ))
    ));

    let canonical_tiny = CanonicalLimits::new(49, 128, 32, 8, 32);
    let canonical_bounded = identify_validity_domain_v1(
        ValidityDomain::unconstrained()
            .with("a", 0.0, 1.0)
            .with("b", 0.0, 1.0),
        canonical_tiny,
        || false,
    );
    assert!(matches!(
        canonical_bounded,
        Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit: 49,
            }
        )) if requested > 49
    ));

    let invalid_limits = identify_validity_domain_v1(
        ValidityDomain::unconstrained(),
        CanonicalLimits::new(4096, 512, 32, 8, 0),
        || false,
    );
    assert_eq!(
        invalid_limits,
        Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    );

    #[derive(Debug)]
    struct CancelAfter {
        successful_polls: usize,
    }
    impl CancellationProbe for CancelAfter {
        fn is_cancelled(&mut self) -> bool {
            if self.successful_polls == 0 {
                true
            } else {
                self.successful_polls -= 1;
                false
            }
        }
    }
    let cancelled = identify_validity_domain_v1(
        ValidityDomain::unconstrained()
            .with("a", 0.0, 1.0)
            .with("b", 0.0, 1.0),
        LIMITS,
        CancelAfter {
            successful_polls: 2,
        },
    );
    assert!(matches!(
        cancelled,
        Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::Cancelled { .. }
        ))
    ));

    let late_domain = ValidityDomain::unconstrained()
        .with("Mach number", 0.0, 1.0)
        .with("α", -1.0, 1.0);
    let poll_count = std::cell::Cell::new(0_usize);
    identify_validity_domain_v1(late_domain.clone(), LIMITS, || {
        poll_count.set(poll_count.get() + 1);
        false
    })
    .expect("baseline poll count");
    let late_cancelled = identify_validity_domain_v1(
        late_domain,
        LIMITS,
        CancelAfter {
            successful_polls: poll_count.get() - 1,
        },
    );
    assert!(matches!(
        late_cancelled,
        Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));
}

#[test]
fn validity_domain_helper_matches_independent_canonical_rows() {
    let domain = ValidityDomain::unconstrained()
        .with("z-axis", -3.0, 4.0)
        .with("a-axis", -0.0, 1.0);
    let identified = identified_domain(domain.clone());
    let mut rows = Vec::new();
    for (axis, (lo, hi)) in domain.bounds() {
        let mut row = Vec::new();
        row.extend_from_slice(
            &u64::try_from(axis.len())
                .expect("axis length")
                .to_le_bytes(),
        );
        row.extend_from_slice(axis.as_bytes());
        row.extend_from_slice(&lo.to_bits().to_le_bytes());
        row.extend_from_slice(&hi.to_bits().to_le_bytes());
        rows.push(row);
    }
    let manual = CanonicalEncoder::<ValidityDomainIdV1, _>::new(LIMITS, || false)
        .expect("validity-domain schema")
        .ordered_bytes(
            Field::new(0, "axes"),
            u64::try_from(rows.len()).expect("row count"),
            rows.iter().map(|row| row.as_slice()),
        )
        .expect("canonical rows")
        .finish()
        .expect("manual validity identity");
    assert_eq!(identified.id(), manual.id());
    assert_eq!(
        identified.receipt().canonical_preimage(),
        manual.canonical_preimage()
    );
}

#[test]
fn raw_validity_domain_receipt_is_only_schema_shaped() {
    let malformed = CanonicalEncoder::<ValidityDomainIdV1, _>::new(LIMITS, || false)
        .expect("validity-domain schema")
        .ordered_bytes(Field::new(0, "axes"), 1, [b"malformed".as_slice()])
        .expect("schema-shaped raw row")
        .finish()
        .expect("raw receipt");
    let admitted = identified_domain(ValidityDomain::unconstrained().with("malformed", 0.0, 1.0));

    assert_ne!(malformed.id(), admitted.id());
    assert_eq!(malformed.audit_record().trust(), TrustState::Unanchored);
    assert_eq!(
        malformed.audit_record().no_claim(),
        NoClaimState::ExternalTrustRequired
    );
}

#[test]
fn model_card_identity_replays_exact_children_and_retains_inputs() {
    let card = model_card_fixture(Some(CALIBRATION_BYTES));
    let first = identified_model_card(card.clone(), Some(CALIBRATION_BYTES.to_vec()));
    let replay = identified_model_card(card.clone(), Some(CALIBRATION_BYTES.to_vec()));
    let (manual_calibration, manual) = manual_model_card_receipts(&card, true, CALIBRATION_BYTES);

    assert_eq!(first.id(), replay.id());
    assert_eq!(first.id(), manual.id());
    assert_eq!(
        first.receipt().canonical_preimage(),
        manual.canonical_preimage()
    );
    assert_eq!(first.calibration_source_id(), Some(manual_calibration.id()));
    assert_eq!(
        first
            .calibration_source_receipt()
            .expect("calibrated receipt")
            .canonical_preimage(),
        manual_calibration.canonical_preimage()
    );
    assert_eq!(
        first.validity_id(),
        identified_domain(card.validity.clone()).id()
    );
    assert_eq!(first.trust_state(), TrustState::Unanchored);
    assert_eq!(
        first.receipt().audit_record().no_claim(),
        NoClaimState::ExternalTrustRequired
    );
    assert_eq!(first.card(), &card);
    assert_eq!(first.calibration_bytes(), Some(CALIBRATION_BYTES));

    let (_, raw_wrong_presence) = manual_model_card_receipts(&card, false, CALIBRATION_BYTES);
    assert_ne!(first.id(), raw_wrong_presence.id());
    assert_eq!(
        raw_wrong_presence.audit_record().trust(),
        TrustState::Unanchored
    );

    let no_calibration_card = model_card_fixture(None);
    let no_calibration = identified_model_card(no_calibration_card.clone(), None);
    let (no_calibration_child, no_calibration_manual) =
        manual_model_card_receipts(&no_calibration_card, false, &[]);
    assert_eq!(no_calibration.id(), no_calibration_manual.id());
    assert_eq!(no_calibration.calibration_source_id(), None);

    let empty_calibration_card = model_card_fixture(Some(&[]));
    let empty_calibration = identified_model_card(empty_calibration_card.clone(), Some(Vec::new()));
    let (empty_calibration_child, empty_calibration_manual) =
        manual_model_card_receipts(&empty_calibration_card, true, &[]);
    assert_eq!(empty_calibration.id(), empty_calibration_manual.id());
    assert_eq!(no_calibration_child.id(), empty_calibration_child.id());
    assert_ne!(no_calibration.id(), empty_calibration.id());
    assert_eq!(
        empty_calibration.calibration_source_id(),
        Some(empty_calibration_child.id())
    );

    let (recovered_card, recovered_calibration) = first.into_parts();
    assert_eq!(recovered_card, card);
    assert_eq!(recovered_calibration.as_deref(), Some(CALIBRATION_BYTES));
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one mutation and refusal matrix shares the exact card/calibration baseline"
)]
fn model_card_identity_binds_every_field_and_refuses_legacy_laundering() {
    let card = model_card_fixture(Some(CALIBRATION_BYTES));
    let base = identified_model_card(card.clone(), Some(CALIBRATION_BYTES.to_vec()));
    let identify = |card: ModelCard, bytes: Option<Vec<u8>>| {
        identify_model_card_v1(card, bytes, LIMITS, || false).expect("valid card mutation")
    };

    let mut name = card.clone();
    name.name.push('!');
    let name_id = identify(name, Some(CALIBRATION_BYTES.to_vec())).id();

    let mut version = card.clone();
    version.version.push_str(".1");
    let version_id = identify(version, Some(CALIBRATION_BYTES.to_vec())).id();

    let mut solid = card.clone();
    solid.ambition = Ambition::Solid;
    let solid_id = identify(solid, Some(CALIBRATION_BYTES.to_vec())).id();
    let mut moonshot = card.clone();
    moonshot.ambition = Ambition::Moonshot;
    let moonshot_id = identify(moonshot, Some(CALIBRATION_BYTES.to_vec())).id();

    let mut assumption = card.clone();
    assumption.assumptions.push("z-last assumption".to_string());
    let assumption_id = identify(assumption, Some(CALIBRATION_BYTES.to_vec())).id();

    let mut validity = card.clone();
    validity.validity = validity.validity.with("β", -2.0, 2.0);
    let validity_id = identify(validity, Some(CALIBRATION_BYTES.to_vec())).id();

    let mut failure = card.clone();
    failure.known_failures.push("z-last failure".to_string());
    let failure_id = identify(failure, Some(CALIBRATION_BYTES.to_vec())).id();

    let mut discrepancy = card.clone();
    discrepancy.discrepancy_rel = 0.1_f64.next_up();
    let discrepancy_id = identify(discrepancy, Some(CALIBRATION_BYTES.to_vec())).id();

    let changed_calibration_bytes = b"calibration-artifact-v2\0binary".to_vec();
    let mut changed_calibration = card.clone();
    changed_calibration.calibration = Some(ProvenanceHash::of_bytes(&changed_calibration_bytes));
    let changed_calibration =
        identify(changed_calibration, Some(changed_calibration_bytes.clone()));

    for (field, mutated_id) in [
        ("name", name_id),
        ("version", version_id),
        ("ambition-solid", solid_id),
        ("ambition-moonshot", moonshot_id),
        ("assumption", assumption_id),
        ("validity", validity_id),
        ("known-failure", failure_id),
        ("discrepancy", discrepancy_id),
        ("calibration-source", changed_calibration.id()),
    ] {
        assert_ne!(base.id(), mutated_id, "{field} must move the model root");
    }
    assert_ne!(solid_id, moonshot_id);
    assert_ne!(
        base.calibration_source_id(),
        changed_calibration.calibration_source_id()
    );

    let collision_hash = ProvenanceHash::of_bytes(FNV1A64_ZERO_PREIMAGE);
    assert_eq!(collision_hash, ProvenanceHash(0));
    assert_eq!(
        collision_hash,
        ProvenanceHash::of_bytes(FNV1A64_ZERO_EXTENSION_COLLISION),
        "fixed fixture must remain a real legacy FNV-1a-64 collision"
    );
    let collision_a = identified_model_card(
        model_card_fixture(Some(FNV1A64_ZERO_PREIMAGE)),
        Some(FNV1A64_ZERO_PREIMAGE.to_vec()),
    );
    let collision_b = identified_model_card(
        model_card_fixture(Some(FNV1A64_ZERO_EXTENSION_COLLISION)),
        Some(FNV1A64_ZERO_EXTENSION_COLLISION.to_vec()),
    );
    assert_eq!(
        collision_a.card().calibration,
        collision_b.card().calibration
    );
    assert_ne!(
        collision_a.calibration_source_id(),
        collision_b.calibration_source_id(),
        "exact source-byte identity must not inherit a legacy FNV collision"
    );
    assert_ne!(
        collision_a.id(),
        collision_b.id(),
        "the typed model root must preserve the exact calibration distinction"
    );

    let mut unsorted = model_card_fixture(None);
    unsorted.assumptions = vec!["z".to_string(), "a".to_string()];
    assert!(matches!(
        identify_model_card_v1(unsorted, None, LIMITS, || false),
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::NonCanonicalSetOrder { index: 1 }
        ))
    ));

    let mut duplicate = model_card_fixture(None);
    duplicate.known_failures = vec!["same".to_string(), "same".to_string()];
    assert!(matches!(
        identify_model_card_v1(duplicate, None, LIMITS, || false),
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::DuplicateSetItem { index: 1 }
        ))
    ));

    let mut invalid_validity = model_card_fixture(None);
    invalid_validity.validity = invalid_validity.validity.with("bad", f64::NAN, 1.0);
    assert!(matches!(
        identify_model_card_v1(invalid_validity, None, LIMITS, || false),
        Err(ModelCardIdentityError::Validity(
            ValidityDomainIdentityError::InvalidBounds { .. }
        ))
    ));

    for discrepancy in [f64::NAN, -0.1, f64::NEG_INFINITY] {
        let mut invalid = model_card_fixture(None);
        invalid.discrepancy_rel = discrepancy;
        assert!(matches!(
            identify_model_card_v1(invalid, None, LIMITS, || false),
            Err(ModelCardIdentityError::InvalidDiscrepancy { bits, .. })
                if bits == discrepancy.to_bits()
        ));
    }

    let missing_bytes = model_card_fixture(Some(CALIBRATION_BYTES));
    assert!(matches!(
        identify_model_card_v1(missing_bytes, None, LIMITS, || false),
        Err(ModelCardIdentityError::CalibrationPresenceMismatch {
            declared: true,
            supplied: false,
        })
    ));
    let unexpected_bytes = model_card_fixture(None);
    assert!(matches!(
        identify_model_card_v1(
            unexpected_bytes,
            Some(CALIBRATION_BYTES.to_vec()),
            LIMITS,
            || false,
        ),
        Err(ModelCardIdentityError::CalibrationPresenceMismatch {
            declared: false,
            supplied: true,
        })
    ));
    let mismatch = model_card_fixture(Some(CALIBRATION_BYTES));
    assert!(matches!(
        identify_model_card_v1(mismatch, Some(b"wrong".to_vec()), LIMITS, || false),
        Err(ModelCardIdentityError::CalibrationCorrelationMismatch { .. })
    ));

    let mut positive_zero = model_card_fixture(None);
    positive_zero.discrepancy_rel = 0.0;
    let mut negative_zero = model_card_fixture(None);
    negative_zero.discrepancy_rel = -0.0;
    assert_ne!(
        identify(positive_zero, None).id(),
        identify(negative_zero, None).id()
    );
    let mut unbounded = model_card_fixture(None);
    unbounded.discrepancy_rel = f64::INFINITY;
    let _unbounded = identify(unbounded, None);
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one exact-limit matrix shares the same model-card child schemas and poll ledger"
)]
fn model_card_identity_enforces_exact_resources_and_cancellation() {
    let field_limits = CanonicalLimits::new(65_536, 1_024, 32, 64, 64);
    let mut exact_name = model_card_fixture(None);
    exact_name.name = "n".repeat(1_024);
    identify_model_card_v1(exact_name, None, field_limits, || false)
        .expect("exact 1024-byte name is admitted");
    let mut oversized_name = model_card_fixture(None);
    oversized_name.name = "n".repeat(1_025);
    assert!(matches!(
        identify_model_card_v1(oversized_name, None, field_limits, || false),
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 1_025,
                limit: 1_024,
            }
        ))
    ));

    let exact_calibration_bytes = vec![0x5a; 1_024];
    let exact_calibration = model_card_fixture(None)
        .with_calibration(ProvenanceHash::of_bytes(&exact_calibration_bytes));
    identify_model_card_v1(
        exact_calibration,
        Some(exact_calibration_bytes),
        field_limits,
        || false,
    )
    .expect("exact 1024-byte calibration source is admitted");
    let oversized_calibration_bytes = vec![0x5a; 1_025];
    let oversized_calibration = model_card_fixture(None)
        .with_calibration(ProvenanceHash::of_bytes(&oversized_calibration_bytes));
    assert!(matches!(
        identify_model_card_v1(
            oversized_calibration,
            Some(oversized_calibration_bytes),
            field_limits,
            || false,
        ),
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 1_025,
                limit: 1_024,
            }
        ))
    ));

    let collection_limits = CanonicalLimits::new(16_384, 8_192, 32, 1, 64);
    assert!(matches!(
        identify_model_card_v1(model_card_fixture(None), None, collection_limits, || false),
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CollectionItems,
                requested: 2,
                limit: 1,
            }
        ))
    ));

    let one_byte = [0x42];
    let minimal_calibrated = ModelCard::new(
        "m",
        "v",
        Ambition::Solid,
        Vec::new(),
        ValidityDomain::unconstrained(),
        Vec::new(),
        0.0,
    )
    .with_calibration(ProvenanceHash::of_bytes(&one_byte));
    let no_chunks = CanonicalLimits::new(16_384, 8_192, 32, 0, 64);
    assert!(matches!(
        identify_model_card_v1(
            minimal_calibrated,
            Some(one_byte.to_vec()),
            no_chunks,
            || false,
        ),
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::StreamChunks,
                requested: 1,
                limit: 0,
            }
        ))
    ));

    let card = model_card_fixture(Some(CALIBRATION_BYTES));
    let baseline = identified_model_card(card.clone(), Some(CALIBRATION_BYTES.to_vec()));
    let frame_limit = baseline
        .receipt()
        .canonical_bytes()
        .checked_sub(1)
        .expect("non-empty model-card frame");
    let frame_limits = CanonicalLimits::new(frame_limit, 8_192, 32, 64, 64);
    assert!(matches!(
        identify_model_card_v1(
            card.clone(),
            Some(CALIBRATION_BYTES.to_vec()),
            frame_limits,
            || false,
        ),
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == frame_limit
    ));

    assert!(matches!(
        identify_model_card_v1(
            card.clone(),
            Some(CALIBRATION_BYTES.to_vec()),
            CanonicalLimits::new(16_384, 8_192, 32, 64, 0),
            || false,
        ),
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    ));
    assert!(matches!(
        identify_model_card_v1(
            card.clone(),
            Some(CALIBRATION_BYTES.to_vec()),
            LIMITS,
            || true,
        ),
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes: 0 }
        ))
    ));

    #[derive(Debug)]
    struct CancelAfter {
        successful_polls: usize,
    }
    impl CancellationProbe for CancelAfter {
        fn is_cancelled(&mut self) -> bool {
            if self.successful_polls == 0 {
                true
            } else {
                self.successful_polls -= 1;
                false
            }
        }
    }
    let uncalibrated = model_card_fixture(None);
    let poll_count = std::cell::Cell::new(0_usize);
    identify_model_card_v1(uncalibrated.clone(), None, LIMITS, || {
        poll_count.set(poll_count.get() + 1);
        false
    })
    .expect("baseline model-card poll count");
    let late = identify_model_card_v1(
        uncalibrated,
        None,
        LIMITS,
        CancelAfter {
            successful_polls: poll_count.get() - 1,
        },
    );
    assert!(matches!(
        late,
        Err(ModelCardIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));
}

#[test]
fn certified_f64_identity_replays_and_excludes_legacy_correlation_values() {
    let certified = certified_fixture();
    let first = identified_certified(certified.clone());
    let replay = identified_certified(certified.clone());
    let manual = manual_certified_receipt(&certified, certified.value, certified.model.in_domain);

    assert_eq!(first.id(), replay.id());
    assert_eq!(first.id(), manual.id());
    assert_eq!(
        first.receipt().canonical_preimage(),
        manual.canonical_preimage()
    );
    assert_eq!(first.id_bytes(), first.receipt().audit_record().id());
    assert_eq!(first.trust_state(), TrustState::Unanchored);
    assert_eq!(
        first.receipt().audit_record().no_claim(),
        NoClaimState::ExternalTrustRequired
    );
    assert_eq!(
        first.validity_id(),
        identified_domain(certified.model.validity.clone()).id()
    );

    let none = Evidence::exact(4.0, ProvenanceHash(7))
        .certified()
        .expect("valid no-statistical fixture");
    let mut half_width = certified.clone().into_evidence();
    half_width.statistical = StatisticalCertificate::HalfWidth {
        half_width: 0.25,
        confidence: 0.95,
    };
    let half_width = half_width.certified().expect("valid half-width fixture");
    for variant in [none, half_width] {
        let manual = manual_certified_receipt(&variant, variant.value, variant.model.in_domain);
        let helper = identified_certified(variant);
        assert_eq!(helper.id(), manual.id());
        assert_eq!(
            helper.receipt().canonical_preimage(),
            manual.canonical_preimage()
        );
    }

    let mut provenance_a = certified.clone().into_evidence();
    provenance_a.provenance = ProvenanceHash(1);
    let provenance_a = provenance_a
        .certified()
        .expect("provenance a remains certified");
    let mut provenance_b = certified.clone().into_evidence();
    provenance_b.provenance = ProvenanceHash(2);
    let provenance_b = provenance_b
        .certified()
        .expect("provenance b remains certified");
    let provenance_a_id = identified_certified(provenance_a).id();
    let provenance_b_id = identified_certified(provenance_b).id();
    assert_eq!(provenance_a_id, provenance_b_id);

    let mut adjoint_a = certified.clone().into_evidence();
    adjoint_a.adjoint_ref = Some(ProvenanceHash(1));
    let adjoint_a = adjoint_a.certified().expect("adjoint a remains certified");
    let mut adjoint_b = certified.clone().into_evidence();
    adjoint_b.adjoint_ref = Some(ProvenanceHash(2));
    let adjoint_b = adjoint_b.certified().expect("adjoint b remains certified");
    let mut no_adjoint = certified.clone().into_evidence();
    no_adjoint.adjoint_ref = None;
    let no_adjoint = no_adjoint
        .certified()
        .expect("no adjoint remains certified");
    assert_eq!(
        identified_certified(adjoint_a).id(),
        identified_certified(adjoint_b).id()
    );
    assert_ne!(first.id(), identified_certified(no_adjoint).id());

    let recovered = first.into_certified();
    assert_eq!(recovered.provenance, ProvenanceHash(0x1111_2222_3333_4444));
    assert_eq!(
        recovered.adjoint_ref,
        Some(ProvenanceHash(0xaaaa_bbbb_cccc_dddd))
    );
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one mutation matrix compares every semantic field to a shared baseline"
)]
fn certified_f64_identity_binds_every_strong_semantic_field() {
    let certified = certified_fixture();
    let base_id = identified_certified(certified.clone()).id();

    let identify_mutation = |evidence: Evidence<f64>| {
        identified_certified(evidence.certified().expect("mutation remains certified")).id()
    };

    let mut scalar = certified.clone().into_evidence();
    scalar.value = 2.0_f64.next_up();
    scalar.qoi = scalar.value;
    let scalar_id = identify_mutation(scalar);

    let mut numerical_kind = certified.clone().into_evidence();
    numerical_kind.numerical.kind = NumericalKind::Exact;
    numerical_kind.numerical.lo = numerical_kind.qoi;
    numerical_kind.numerical.hi = numerical_kind.qoi;
    let numerical_kind_id = identify_mutation(numerical_kind);

    let mut numerical_lo = certified.clone().into_evidence();
    numerical_lo.numerical.lo = 1.0_f64.next_up();
    let numerical_lo_id = identify_mutation(numerical_lo);

    let mut numerical_hi = certified.clone().into_evidence();
    numerical_hi.numerical.hi = 3.0_f64.next_down();
    let numerical_hi_id = identify_mutation(numerical_hi);

    let mut half_width_base = certified.clone().into_evidence();
    half_width_base.statistical = StatisticalCertificate::HalfWidth {
        half_width: 0.25,
        confidence: 0.95,
    };
    let half_width_base_id = identify_mutation(half_width_base);

    let mut e_value_e = certified.clone().into_evidence();
    e_value_e.statistical = StatisticalCertificate::EValue {
        e: 8.0_f64.next_up(),
        alpha: 0.05,
    };
    let e_value_e_id = identify_mutation(e_value_e);

    let mut e_value_alpha = certified.clone().into_evidence();
    e_value_alpha.statistical = StatisticalCertificate::EValue {
        e: 8.0,
        alpha: 0.05_f64.next_up(),
    };
    let e_value_alpha_id = identify_mutation(e_value_alpha);

    let mut half_width_value = certified.clone().into_evidence();
    half_width_value.statistical = StatisticalCertificate::HalfWidth {
        half_width: 0.25_f64.next_up(),
        confidence: 0.95,
    };
    let half_width_value_id = identify_mutation(half_width_value);

    let mut half_width_confidence = certified.clone().into_evidence();
    half_width_confidence.statistical = StatisticalCertificate::HalfWidth {
        half_width: 0.25,
        confidence: 0.95_f64.next_down(),
    };
    let half_width_confidence_id = identify_mutation(half_width_confidence);

    let mut card = certified.clone().into_evidence();
    card.model.cards.push("card-c".to_string());
    let card_id = identify_mutation(card);

    let mut assumption = certified.clone().into_evidence();
    assumption
        .model
        .assumptions
        .push("assumption-c".to_string());
    let assumption_id = identify_mutation(assumption);

    let mut validity = certified.clone().into_evidence();
    validity.model.validity = validity.model.validity.with("β", -2.0, 2.0);
    let validity_id = identify_mutation(validity);

    let mut discrepancy = certified.clone().into_evidence();
    discrepancy.model.discrepancy_rel = 0.125_f64.next_up();
    let discrepancy_id = identify_mutation(discrepancy);

    let mut sensitivity = certified.clone().into_evidence();
    sensitivity
        .sensitivity
        .d_qoi
        .insert("α".to_string(), 2.5_f64.next_up());
    let sensitivity_id = identify_mutation(sensitivity);

    let mut sensitivity_name = certified.clone().into_evidence();
    let derivative = sensitivity_name
        .sensitivity
        .d_qoi
        .remove("α")
        .expect("fixture sensitivity key");
    sensitivity_name
        .sensitivity
        .d_qoi
        .insert("β".to_string(), derivative);
    let sensitivity_name_id = identify_mutation(sensitivity_name);

    for (field, mutated_id) in [
        ("value-and-qoi", scalar_id),
        ("numerical-kind-and-exact-bounds", numerical_kind_id),
        ("numerical-lo", numerical_lo_id),
        ("numerical-hi", numerical_hi_id),
        ("statistical-variant", half_width_base_id),
        ("e-value-e", e_value_e_id),
        ("e-value-alpha", e_value_alpha_id),
        ("model-card", card_id),
        ("model-assumption", assumption_id),
        ("model-validity", validity_id),
        ("model-discrepancy", discrepancy_id),
        ("sensitivity-value", sensitivity_id),
        ("sensitivity-name", sensitivity_name_id),
    ] {
        assert_ne!(base_id, mutated_id, "{field} must move the root");
    }
    assert_ne!(half_width_base_id, half_width_value_id);
    assert_ne!(half_width_base_id, half_width_confidence_id);

    let positive_zero = identified_certified(
        Evidence::exact(0.0, ProvenanceHash(1))
            .certified()
            .expect("positive zero exact"),
    );
    let negative_zero = identified_certified(
        Evidence::exact(-0.0, ProvenanceHash(2))
            .certified()
            .expect("negative zero exact"),
    );
    assert_ne!(positive_zero.id(), negative_zero.id());

    let mut positive_discrepancy = certified.clone().into_evidence();
    positive_discrepancy.model.discrepancy_rel = 0.0;
    let mut negative_discrepancy = certified.clone().into_evidence();
    negative_discrepancy.model.discrepancy_rel = -0.0;
    assert_ne!(
        identify_mutation(positive_discrepancy),
        identify_mutation(negative_discrepancy)
    );

    let mut first_nan = certified.clone().into_evidence();
    first_nan.sensitivity.d_qoi.insert(
        "nan-sensitivity".to_string(),
        f64::from_bits(0x7ff8_0000_0000_0001),
    );
    let first_nan = identify_mutation(first_nan);
    let mut second_nan = certified.clone().into_evidence();
    second_nan.sensitivity.d_qoi.insert(
        "nan-sensitivity".to_string(),
        f64::from_bits(0x7ff8_0000_0000_0002),
    );
    let second_nan = identify_mutation(second_nan);
    assert_ne!(first_nan, second_nan);

    let mut unbounded_discrepancy = certified.into_evidence();
    unbounded_discrepancy.model.discrepancy_rel = f64::INFINITY;
    assert_ne!(base_id, identify_mutation(unbounded_discrepancy));
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one refusal matrix shares exact limits and cancellation accounting"
)]
fn certified_f64_identity_refuses_noncanonical_sets_resources_and_cancellation() {
    let mut unsorted = certified_fixture().into_evidence();
    unsorted.model.cards = vec!["z-card".to_string(), "a-card".to_string()];
    let unsorted = unsorted
        .certified()
        .expect("set order is not certification");
    assert!(matches!(
        identify_certified_f64_evidence_v1(unsorted, LIMITS, || false),
        Err(CertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::NonCanonicalSetOrder { index: 1 }
        ))
    ));

    let mut duplicate = certified_fixture().into_evidence();
    duplicate.model.assumptions = vec!["same".to_string(), "same".to_string()];
    let duplicate = duplicate
        .certified()
        .expect("duplicates are not certification");
    assert!(matches!(
        identify_certified_f64_evidence_v1(duplicate, LIMITS, || false),
        Err(CertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::DuplicateSetItem { index: 1 }
        ))
    ));

    let certified = certified_fixture();
    let helper = identified_certified(certified.clone());
    let raw_false_in_domain = manual_certified_receipt(&certified, certified.value, false);
    let raw_value_mismatch = manual_certified_receipt(&certified, certified.value.next_up(), true);
    assert_ne!(helper.id(), raw_false_in_domain.id());
    assert_ne!(helper.id(), raw_value_mismatch.id());
    assert_eq!(
        raw_false_in_domain.audit_record().trust(),
        TrustState::Unanchored
    );

    let frame_limit = helper
        .receipt()
        .canonical_bytes()
        .checked_sub(1)
        .expect("non-empty certified-f64 frame");
    let frame_limits = CanonicalLimits::new(frame_limit, 8_192, 32, 64, 64);
    assert!(matches!(
        identify_certified_f64_evidence_v1(certified.clone(), frame_limits, || false),
        Err(CertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == frame_limit
    ));

    let mut oversized_card = certified.clone().into_evidence();
    oversized_card.model.cards = vec!["c".repeat(250)];
    let oversized_card = oversized_card
        .certified()
        .expect("card bytes are not certification");
    let field_limits = CanonicalLimits::new(16_384, 256, 32, 64, 64);
    assert!(matches!(
        identify_certified_f64_evidence_v1(oversized_card, field_limits, || false),
        Err(CertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 266,
                limit: 256,
            }
        ))
    ));

    let sensitivity_field_limits = CanonicalLimits::new(16_384, 512, 32, 64, 64);
    let mut exact_sensitivity_field = certified.clone().into_evidence();
    exact_sensitivity_field.sensitivity.d_qoi = BTreeMap::from([("s".repeat(480), 1.0)]);
    let exact_sensitivity_field = exact_sensitivity_field
        .certified()
        .expect("exact-limit sensitivity remains certified");
    identify_certified_f64_evidence_v1(exact_sensitivity_field, sensitivity_field_limits, || false)
        .expect("exact 512-byte sensitivity field is admitted");

    let mut oversized_sensitivity_field = certified.clone().into_evidence();
    oversized_sensitivity_field.sensitivity.d_qoi = BTreeMap::from([("s".repeat(481), 1.0)]);
    let oversized_sensitivity_field = oversized_sensitivity_field
        .certified()
        .expect("over-limit sensitivity remains certified");
    assert!(matches!(
        identify_certified_f64_evidence_v1(
            oversized_sensitivity_field,
            sensitivity_field_limits,
            || false,
        ),
        Err(CertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 513,
                limit: 512,
            }
        ))
    ));

    let mut collection_bounded = Evidence::exact(1.0, ProvenanceHash(9));
    collection_bounded.sensitivity.d_qoi =
        BTreeMap::from([("a".to_string(), 1.0), ("b".to_string(), 2.0)]);
    let collection_bounded = collection_bounded
        .certified()
        .expect("collection-bounded fixture remains certified");
    let collection_limits = CanonicalLimits::new(16_384, 8_192, 32, 1, 64);
    assert!(matches!(
        identify_certified_f64_evidence_v1(collection_bounded, collection_limits, || false),
        Err(CertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CollectionItems,
                requested: 2,
                limit: 1,
            }
        ))
    ));

    let mut chunk_bounded = certified.clone().into_evidence();
    chunk_bounded.sensitivity.d_qoi = BTreeMap::from([
        ("a".to_string(), 1.0),
        ("b".to_string(), 2.0),
        ("c".to_string(), 3.0),
    ]);
    let chunk_bounded = chunk_bounded
        .certified()
        .expect("sensitivity is not certification");
    let chunk_limits = CanonicalLimits::new(16_384, 8_192, 32, 8, 64);
    assert!(matches!(
        identify_certified_f64_evidence_v1(chunk_bounded, chunk_limits, || false),
        Err(CertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::StreamChunks,
                requested: 9,
                limit: 8,
            }
        ))
    ));

    assert!(matches!(
        identify_certified_f64_evidence_v1(certified.clone(), LIMITS, || true),
        Err(CertifiedF64EvidenceIdentityError::Validity(
            ValidityDomainIdentityError::Canonical(CanonicalError::Cancelled { .. })
        ))
    ));
    assert!(matches!(
        identify_certified_f64_evidence_v1(
            certified.clone(),
            CanonicalLimits::new(4096, 512, 32, 64, 0),
            || false,
        ),
        Err(CertifiedF64EvidenceIdentityError::Validity(
            ValidityDomainIdentityError::Canonical(CanonicalError::InvalidLimits(
                "cancellation_poll_bytes must be positive"
            ))
        ))
    ));

    #[derive(Debug)]
    struct CancelAfter {
        successful_polls: usize,
    }
    impl CancellationProbe for CancelAfter {
        fn is_cancelled(&mut self) -> bool {
            if self.successful_polls == 0 {
                true
            } else {
                self.successful_polls -= 1;
                false
            }
        }
    }
    let poll_count = std::cell::Cell::new(0_usize);
    identify_certified_f64_evidence_v1(certified.clone(), LIMITS, || {
        poll_count.set(poll_count.get() + 1);
        false
    })
    .expect("baseline certified poll count");
    let late_cancelled = identify_certified_f64_evidence_v1(
        certified,
        LIMITS,
        CancelAfter {
            successful_polls: poll_count.get() - 1,
        },
    );
    assert!(matches!(
        late_cancelled,
        Err(CertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));
}

#[test]
fn typed_receipts_are_replay_stable_and_explicitly_unanchored() {
    let first = source_node("observation-a", verified(-1.0, 2.0));
    let replay = source_node("observation-a", verified(-1.0, 2.0));

    assert_eq!(first.id(), replay.id());
    assert_eq!(
        first.receipt().canonical_preimage(),
        replay.receipt().canonical_preimage()
    );

    let audit = first.receipt().audit_record();
    assert_eq!(audit.trust(), TrustState::Unanchored);
    assert_eq!(audit.no_claim(), NoClaimState::ExternalTrustRequired);
    assert_eq!(audit.id(), first.id_bytes());
    assert_eq!(first.trust_state(), TrustState::Unanchored);
}

#[test]
fn every_available_semantic_mutation_moves_the_typed_root() {
    let source_a = source_node("observation-a", verified(0.0, 1.0));
    let source_b = source_node("observation-b", verified(1.0, 2.0));
    let source_payload_changed = source_node("observation-a+", verified(0.0, 1.0));
    let output_changed = source_node("observation-a", verified(0.0, 1.5));

    assert_ne!(source_a.id(), source_b.id());
    assert_ne!(source_a.id(), source_payload_changed.id());
    assert_ne!(source_a.id(), output_changed.id());

    let add = derived(ColorEvidenceCompositionOpV1::Add, &source_a, &source_b).expect("valid add");
    let mul =
        derived(ColorEvidenceCompositionOpV1::Mul, &source_a, &source_b).expect("valid multiply");
    assert_ne!(add.id(), mul.id());
}

#[test]
fn commutative_multisets_normalize_order_but_preserve_multiplicity() {
    let a = source_node("a", verified(0.0, 1.0));
    let b = source_node("b", verified(1.0, 2.0));

    for operation in [
        ColorEvidenceCompositionOpV1::Add,
        ColorEvidenceCompositionOpV1::Mul,
        ColorEvidenceCompositionOpV1::Hull,
    ] {
        let ab = derived(operation, &a, &b).expect("a op b");
        let ba = derived(operation, &b, &a).expect("b op a");
        assert_eq!(ab.id(), ba.id());
        assert_eq!(ab.color(), ba.color());
    }

    let ab = derived(ColorEvidenceCompositionOpV1::Add, &a, &b).expect("a+b");
    let aa = derived(ColorEvidenceCompositionOpV1::Add, &a, &a).expect("a+a");
    assert_ne!(ab.id(), aa.id(), "a+a must not collapse into a+b");
}

#[test]
fn composition_recomputes_the_current_color_algebra() {
    let a = source_node("a", verified(-1.0, 2.0));
    let b = source_node("b", verified(3.0, 4.0));

    let add = derived(ColorEvidenceCompositionOpV1::Add, &a, &b).expect("add");
    let mul = derived(ColorEvidenceCompositionOpV1::Mul, &a, &b).expect("multiply");
    let hull = derived(ColorEvidenceCompositionOpV1::Hull, &a, &b).expect("hull");

    assert_eq!(
        add.color(),
        &verified(2.0_f64.next_down(), 6.0_f64.next_up())
    );
    assert_eq!(
        mul.color(),
        &verified((-4.0_f64).next_down(), 8.0_f64.next_up())
    );
    assert_eq!(hull.color(), &verified(-1.0, 4.0));
    assert_eq!(add.operation(), ColorEvidenceOperationV1::Add);
    assert_eq!(add.kind(), ColorEvidenceNodeKindV1::Composition);
    assert_eq!(
        add.parent_semantics(),
        ColorEvidenceParentSemanticsV1::CommutativeMultiset
    );
    assert_ne!(add.id(), mul.id());
    assert_ne!(mul.id(), hull.id());
}

#[test]
fn source_domain_version_and_payload_are_all_identity_bearing() {
    let a = identify_color_evidence_source_v1("domain-a", 1, b"same", LIMITS, || false)
        .expect("domain a v1");
    let b = identify_color_evidence_source_v1("domain-b", 1, b"same", LIMITS, || false)
        .expect("domain b v1");
    let c = identify_color_evidence_source_v1("domain-a", 2, b"same", LIMITS, || false)
        .expect("domain a v2");
    let d = identify_color_evidence_source_v1("domain-a", 1, b"different", LIMITS, || false)
        .expect("different bytes");
    assert_ne!(a.id(), b.id());
    assert_ne!(a.id(), c.id());
    assert_ne!(a.id(), d.id());
    assert_eq!(a.trust_state(), TrustState::Unanchored);

    assert_eq!(
        identify_color_evidence_source_v1("", 1, b"x", LIMITS, || false),
        Err(ColorEvidenceIdentityError::EmptySourceDomain)
    );
    assert_eq!(
        identify_color_evidence_source_v1("domain-a", 0, b"x", LIMITS, || false),
        Err(ColorEvidenceIdentityError::ZeroSourceSchemaVersion)
    );
}

#[test]
fn malformed_output_and_cancellation_publish_no_identity() {
    let source = source_receipt("source");
    let malformed =
        identify_color_evidence_source_node_v1(&source, verified(2.0, 1.0), LIMITS, || false);
    assert!(matches!(
        malformed,
        Err(ColorEvidenceIdentityError::MalformedColor(_))
    ));

    #[derive(Debug)]
    struct CancelNow;
    impl CancellationProbe for CancelNow {
        fn is_cancelled(&mut self) -> bool {
            true
        }
    }

    let cancelled =
        identify_color_evidence_source_node_v1(&source, verified(0.0, 1.0), LIMITS, CancelNow);
    assert!(matches!(
        cancelled,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { .. }
        ))
    ));

    let left = source_node("cancel-left", verified(0.0, 1.0));
    let right = source_node("cancel-right", verified(1.0, 2.0));
    let cancelled_composition = compose_color_evidence_nodes_v1(
        ColorEvidenceCompositionOpV1::Add,
        &left,
        &right,
        LIMITS,
        CancelNow,
    );
    assert!(matches!(
        cancelled_composition,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { .. }
        ))
    ));

    #[derive(Debug)]
    struct CancelAfter {
        successful_polls: usize,
    }
    impl CancellationProbe for CancelAfter {
        fn is_cancelled(&mut self) -> bool {
            if self.successful_polls == 0 {
                true
            } else {
                self.successful_polls -= 1;
                false
            }
        }
    }
    let midstream_regime = ValidityDomain::unconstrained()
        .with("axis-c", 0.0, 1.0)
        .with("axis-a", 0.0, 1.0)
        .with("axis-b", 0.0, 1.0);
    let midstream_cancelled = identify_color_evidence_source_node_v1(
        &source,
        Color::Validated {
            regime: midstream_regime,
            dataset: "dataset-v1".to_string(),
        },
        LIMITS,
        CancelAfter {
            successful_polls: 2,
        },
    );
    assert!(matches!(
        midstream_cancelled,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { .. }
        ))
    ));

    let stride_limits = CanonicalLimits::new(4096, 1024, 32, 64, 16);
    let stride_cancelled = identify_color_evidence_source_node_v1(
        &source,
        Color::Validated {
            regime: ValidityDomain::unconstrained().with("x", 0.0, 1.0),
            dataset: "d".repeat(64),
        },
        stride_limits,
        CancelAfter {
            successful_polls: 6,
        },
    );
    assert_eq!(
        stride_cancelled,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes: 26 }
        ))
    );

    let mut regime = ValidityDomain::unconstrained();
    for index in 0..8 {
        regime = regime.with(format!("axis-{index}"), 0.0, 1.0);
    }
    let tiny = CanonicalLimits::new(512, 256, 32, 2, 64);
    let bounded = identify_color_evidence_source_node_v1(
        &source,
        Color::Validated {
            regime,
            dataset: "dataset-v1".to_string(),
        },
        tiny,
        || false,
    );
    assert!(matches!(
        bounded,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded { .. }
        ))
    ));

    let field_tiny = CanonicalLimits::new(4096, 33, 32, 64, 64);
    let field_bounded =
        identify_color_evidence_source_node_v1(&source, verified(0.0, 1.0), field_tiny, || false);
    assert!(matches!(
        field_bounded,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 34,
                limit: 33,
            }
        ))
    ));
}

#[test]
fn source_child_role_and_schema_are_bound_by_the_node_schema() {
    let foreign = CanonicalEncoder::<SourceId<ForeignSourceSchema>, _>::new(LIMITS, || false)
        .expect("foreign source schema")
        .u64(Field::new(0, "value"), 7)
        .expect("foreign source value")
        .finish()
        .expect("foreign source receipt");

    let refusal = CanonicalEncoder::<ColorEvidenceNodeIdV1, _>::new(LIMITS, || false)
        .expect("node schema")
        .variant(Field::new(0, "node-kind"), 1, &[])
        .expect("node kind")
        .variant(Field::new(1, "operation"), 1, &[])
        .expect("operation")
        .variant(Field::new(2, "parent-semantics"), 1, &[])
        .expect("parent semantics")
        .u64(
            Field::new(3, "color-algebra-version"),
            u64::from(COLOR_ALGEBRA_VERSION),
        )
        .expect("algebra version")
        .ordered_children(Field::new(4, "source"), 1, [foreign.id()])
        .expect_err("foreign source role/schema must refuse");
    assert!(matches!(
        refusal,
        CanonicalError::ChildBindingMismatch {
            field: "source",
            ..
        }
    ));
}

#[test]
fn source_node_helper_matches_every_color_canonical_encoding() {
    let source = source_receipt("manual-parity");
    let validated_regime = ValidityDomain::unconstrained()
        .with("z-axis", -3.0, 4.0)
        .with("a-axis", 0.0, 1.0);
    let colors = [
        verified(-2.0, 3.0),
        verified(f64::NEG_INFINITY, f64::INFINITY),
        Color::Validated {
            regime: validated_regime,
            dataset: "dataset-v1".to_string(),
        },
        Color::Estimated {
            estimator: "estimator-v1".to_string(),
            dispersion: -0.0,
        },
        Color::Estimated {
            estimator: "unbounded-estimator-v1".to_string(),
            dispersion: f64::INFINITY,
        },
    ];

    for color in colors {
        let helper =
            identify_color_evidence_source_node_v1(&source, color.clone(), LIMITS, || false)
                .expect("helper source node");
        let color_bytes = color.canonical_bytes();
        let manual = CanonicalEncoder::<ColorEvidenceNodeIdV1, _>::new(LIMITS, || false)
            .expect("node schema")
            .variant(Field::new(0, "node-kind"), 1, &[])
            .expect("node kind")
            .variant(Field::new(1, "operation"), 1, &[])
            .expect("operation")
            .variant(Field::new(2, "parent-semantics"), 1, &[])
            .expect("parent semantics")
            .u64(
                Field::new(3, "color-algebra-version"),
                u64::from(COLOR_ALGEBRA_VERSION),
            )
            .expect("algebra version")
            .ordered_children(Field::new(4, "source"), 1, [source.id()])
            .expect("typed source child")
            .bytes(Field::new(5, "output-color"), &color_bytes)
            .expect("canonical color")
            .ordered_bytes(Field::new(6, "parents"), 0, core::iter::empty::<&[u8]>())
            .expect("empty parent list")
            .finish()
            .expect("manual source node");
        assert_eq!(helper.id(), manual.id());
        assert_eq!(
            helper.receipt().canonical_preimage(),
            manual.canonical_preimage()
        );
    }
}

#[test]
fn composition_helper_matches_descriptor_bound_parent_construction() {
    let a = source_node("composition-parity-a", verified(-1.0, 2.0));
    let b = source_node("composition-parity-b", verified(3.0, 4.0));
    let helper = derived(ColorEvidenceCompositionOpV1::Add, &a, &b).expect("helper add");
    let mut parents = [a.id(), b.id()];
    parents.sort_unstable();
    let parent_rows = parents.map(parent_reference_bytes);
    let output_bytes = verified(2.0_f64.next_down(), 6.0_f64.next_up()).canonical_bytes();

    let manual = CanonicalEncoder::<ColorEvidenceNodeIdV1, _>::new(LIMITS, || false)
        .expect("node schema")
        .variant(Field::new(0, "node-kind"), 2, &[])
        .expect("composition kind")
        .variant(Field::new(1, "operation"), 2, &[])
        .expect("add operation")
        .variant(Field::new(2, "parent-semantics"), 2, &[])
        .expect("commutative multiset")
        .u64(
            Field::new(3, "color-algebra-version"),
            u64::from(COLOR_ALGEBRA_VERSION),
        )
        .expect("algebra version")
        .ordered_children(
            Field::new(4, "source"),
            0,
            core::iter::empty::<ColorEvidenceSourceIdV1>(),
        )
        .expect("no source child")
        .bytes(Field::new(5, "output-color"), &output_bytes)
        .expect("canonical color")
        .ordered_bytes(
            Field::new(6, "parents"),
            2,
            parent_rows.iter().map(|row| row.as_slice()),
        )
        .expect("descriptor-bound parents")
        .finish()
        .expect("manual composition node");

    assert_eq!(helper.id(), manual.id());
    assert_eq!(
        helper.receipt().canonical_preimage(),
        manual.canonical_preimage()
    );

    let mut wrong_role_rows = parent_rows;
    wrong_role_rows[0][0] ^= 0xff;
    let wrong_role = CanonicalEncoder::<ColorEvidenceNodeIdV1, _>::new(LIMITS, || false)
        .expect("node schema")
        .variant(Field::new(0, "node-kind"), 2, &[])
        .expect("composition kind")
        .variant(Field::new(1, "operation"), 2, &[])
        .expect("add operation")
        .variant(Field::new(2, "parent-semantics"), 2, &[])
        .expect("commutative multiset")
        .u64(
            Field::new(3, "color-algebra-version"),
            u64::from(COLOR_ALGEBRA_VERSION),
        )
        .expect("algebra version")
        .ordered_children(
            Field::new(4, "source"),
            0,
            core::iter::empty::<ColorEvidenceSourceIdV1>(),
        )
        .expect("no source child")
        .bytes(Field::new(5, "output-color"), &output_bytes)
        .expect("canonical color")
        .ordered_bytes(
            Field::new(6, "parents"),
            2,
            wrong_role_rows.iter().map(|row| row.as_slice()),
        )
        .expect("schema-shaped but semantically foreign parent rows")
        .finish()
        .expect("foreign manual node");
    assert_ne!(helper.id(), wrong_role.id());
    assert_ne!(manual.id(), wrong_role.id());
}
