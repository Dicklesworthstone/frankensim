//! G0/G3/G4 coverage for typed evidence identities.
//!
//! A retained cross-ISA known-answer vector is still required for G5.

use std::collections::BTreeMap;

use fs_blake3::identity::{
    CancellationProbe, CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field,
    FieldSpec, NoClaimState, SchemaId, SemanticId, SourceId, StrongIdentity, TrustState, WireType,
};
use fs_evidence::{
    Ambition, COLOR_ALGEBRA_VERSION, CardBoundModelEvidenceIdV1,
    CardBoundModelEvidenceIdentityError, CardBoundModelEvidenceReceiptV1, Certified,
    CertifiedF64DecisionAssessmentIdV1, CertifiedF64DecisionAssessmentIdentityError,
    CertifiedF64DecisionAssessmentReceiptV1, CertifiedF64EvidenceIdV1,
    CertifiedF64EvidenceIdentityError, CertifiedF64EvidenceReceiptV1, CertifiedF64SourceIdV1,
    CertifiedF64SourceIdentityError, CertifiedF64SourceReceiptV1, CertifiedF64SourceV1, Color,
    ColorEvidenceCompositionOpV1, ColorEvidenceIdentityError, ColorEvidenceNodeIdV1,
    ColorEvidenceNodeIdentitySchemaV1, ColorEvidenceNodeKindV1, ColorEvidenceNodeV1,
    ColorEvidenceOperationV1, ColorEvidenceParentSemanticsV1, ColorEvidenceSourceIdV1,
    ColorEvidenceSourceV1, DECISION_ASSESSMENT_ALGORITHM_VERSION_V1, DecisionStatus,
    DiscrepancyBand, DiscrepancyBandIdV1, DiscrepancyBandIdentityError, DiscrepancyBandReceiptV1,
    EscalationAdvice, Evidence, FidelityPair, FidelityPairIdV1, FidelityPairIdentityError,
    FidelityPairReceiptV1, IdentifiedCardBoundModelEvidenceV1,
    IdentifiedCertifiedF64DecisionAssessmentV1, IdentifiedCertifiedF64EvidenceV1,
    IdentifiedModelBracketV1, IdentifiedModelCardV1, IdentifiedModelEvidenceV1,
    IdentifiedSourcedCertifiedF64EvidenceV1, IdentifiedValidityDomainV1,
    MAX_CERTIFIED_F64_SOURCE_FIELD_BYTES_V1, MAX_MODEL_BRACKET_MEMBERS_V1, ModelBracket,
    ModelBracketIdV1, ModelBracketIdentityError, ModelBracketReceiptV1, ModelCard,
    ModelCardCalibrationSourceIdV1, ModelCardCalibrationSourceReceiptV1, ModelCardIdV1,
    ModelCardIdentityError, ModelCardIdentitySchemaV1, ModelCardReceiptV1, ModelEvidence,
    ModelEvidenceIdV1, ModelEvidenceIdentityError, ModelEvidenceIdentitySchemaV1,
    ModelEvidenceReceiptV1, NumericalCertificate, NumericalCertificateIdV1,
    NumericalCertificateIdentityError, NumericalCertificateReceiptV1, NumericalKind,
    ProvenanceHash, SensitivitySummary, SourcedCertifiedF64EvidenceIdV1,
    SourcedCertifiedF64EvidenceIdentityError, SourcedCertifiedF64EvidenceReceiptV1,
    StatisticalCertificate, StatisticalCertificateIdV1, StatisticalCertificateIdentityError,
    StatisticalCertificateReceiptV1, UncertaintySource, ValidityDomain, ValidityDomainIdV1,
    ValidityDomainIdentityError, compose_color_evidence_nodes_v1,
    identify_card_bound_model_evidence_v1, identify_certified_f64_decision_assessment_v1,
    identify_certified_f64_evidence_v1, identify_certified_f64_source_v1,
    identify_color_evidence_source_node_v1, identify_color_evidence_source_v1,
    identify_discrepancy_band_v1, identify_fidelity_pair_v1, identify_model_bracket_v1,
    identify_model_card_v1, identify_model_evidence_v1, identify_numerical_certificate_v1,
    identify_sourced_certified_f64_evidence_v1, identify_statistical_certificate_v1,
    identify_validity_domain_v1,
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

fn numerical_kind_tag(kind: NumericalKind) -> u32 {
    match kind {
        NumericalKind::Exact => 1,
        NumericalKind::Enclosure => 2,
        NumericalKind::Estimate => 3,
        NumericalKind::NoClaim => 4,
    }
}

fn manual_numerical_certificate_receipt(
    certificate: NumericalCertificate,
) -> NumericalCertificateReceiptV1 {
    CanonicalEncoder::<NumericalCertificateIdV1, _>::new(LIMITS, || false)
        .expect("numerical-certificate schema")
        .variant(
            Field::new(0, "kind"),
            numerical_kind_tag(certificate.kind),
            &[],
        )
        .expect("numerical kind")
        .u64(Field::new(1, "lo-ieee754-bits"), certificate.lo.to_bits())
        .expect("numerical lower bits")
        .u64(Field::new(2, "hi-ieee754-bits"), certificate.hi.to_bits())
        .expect("numerical upper bits")
        .finish()
        .expect("manual numerical-certificate identity")
}

fn manual_statistical_certificate_receipt(
    certificate: StatisticalCertificate,
) -> StatisticalCertificateReceiptV1 {
    let mut payload = [0_u8; 16];
    let (tag, payload_len) = match certificate {
        StatisticalCertificate::None => (1, 0),
        StatisticalCertificate::EValue { e, alpha } => {
            payload[..8].copy_from_slice(&e.to_bits().to_le_bytes());
            payload[8..].copy_from_slice(&alpha.to_bits().to_le_bytes());
            (2, payload.len())
        }
        StatisticalCertificate::HalfWidth {
            half_width,
            confidence,
        } => {
            payload[..8].copy_from_slice(&half_width.to_bits().to_le_bytes());
            payload[8..].copy_from_slice(&confidence.to_bits().to_le_bytes());
            (3, payload.len())
        }
    };
    CanonicalEncoder::<StatisticalCertificateIdV1, _>::new(LIMITS, || false)
        .expect("statistical-certificate schema")
        .variant(Field::new(0, "certificate"), tag, &payload[..payload_len])
        .expect("statistical certificate")
        .finish()
        .expect("manual statistical-certificate identity")
}

fn fidelity_pair_fixture() -> FidelityPair {
    FidelityPair {
        params: BTreeMap::from([("alpha".to_string(), -0.0), ("mach".to_string(), 0.7)]),
        lo_fi: 4.0,
        hi_fi: 5.0,
    }
}

fn manual_fidelity_pair_receipt(pair: &FidelityPair) -> FidelityPairReceiptV1 {
    let parameter_rows: Vec<Vec<u8>> = pair
        .params
        .iter()
        .map(|(name, value)| {
            let mut row = Vec::new();
            row.extend_from_slice(
                &u64::try_from(name.len())
                    .expect("parameter-name length")
                    .to_le_bytes(),
            );
            row.extend_from_slice(name.as_bytes());
            row.extend_from_slice(&value.to_bits().to_le_bytes());
            row
        })
        .collect();
    CanonicalEncoder::<FidelityPairIdV1, _>::new(LIMITS, || false)
        .expect("fidelity-pair schema")
        .ordered_bytes(
            Field::new(0, "parameters"),
            u64::try_from(parameter_rows.len()).expect("parameter count"),
            parameter_rows.iter().map(Vec::as_slice),
        )
        .expect("parameter point")
        .u64(
            Field::new(1, "lo-fi-qoi-ieee754-bits"),
            pair.lo_fi.to_bits(),
        )
        .expect("low-fidelity QoI")
        .u64(
            Field::new(2, "hi-fi-qoi-ieee754-bits"),
            pair.hi_fi.to_bits(),
        )
        .expect("high-fidelity QoI")
        .finish()
        .expect("manual fidelity-pair identity")
}

fn manual_discrepancy_band_receipt(band: DiscrepancyBand) -> DiscrepancyBandReceiptV1 {
    CanonicalEncoder::<DiscrepancyBandIdV1, _>::new(LIMITS, || false)
        .expect("discrepancy-band schema")
        .u64(
            Field::new(0, "mean-rel-ieee754-bits"),
            band.mean_rel.to_bits(),
        )
        .expect("mean discrepancy")
        .u64(
            Field::new(1, "max-rel-ieee754-bits"),
            band.max_rel.to_bits(),
        )
        .expect("maximum discrepancy")
        .finish()
        .expect("manual discrepancy-band identity")
}

fn model_bracket(members: &[(&str, f64)]) -> ModelBracket {
    members
        .iter()
        .try_fold(ModelBracket::new(), |bracket, &(name, qoi)| {
            bracket.try_with_member(name, qoi)
        })
        .expect("valid model-bracket fixture")
}

fn identified_model_bracket(bracket: ModelBracket) -> IdentifiedModelBracketV1 {
    identify_model_bracket_v1(bracket, LIMITS, || false).expect("valid model-bracket identity")
}

fn manual_model_bracket_receipt(members: &[(&str, f64)]) -> ModelBracketReceiptV1 {
    let mut members = members.to_vec();
    members.sort_unstable_by(|(left, _), (right, _)| left.cmp(right));
    let rows: Vec<Vec<u8>> = members
        .iter()
        .map(|(name, qoi)| {
            let mut row = Vec::new();
            row.extend_from_slice(
                &u64::try_from(name.len())
                    .expect("model-name length")
                    .to_le_bytes(),
            );
            row.extend_from_slice(name.as_bytes());
            row.extend_from_slice(&qoi.to_bits().to_le_bytes());
            row
        })
        .collect();
    CanonicalEncoder::<ModelBracketIdV1, _>::new(LIMITS, || false)
        .expect("model-bracket schema")
        .ordered_bytes(
            Field::new(0, "members"),
            u64::try_from(rows.len()).expect("model-bracket member count"),
            rows.iter().map(Vec::as_slice),
        )
        .expect("model-bracket members")
        .finish()
        .expect("manual model-bracket identity")
}

fn model_evidence_fixture() -> ModelEvidence {
    ModelEvidence {
        cards: vec!["card-a".to_string(), "card-b".to_string()],
        assumptions: vec!["assumption-a".to_string(), "assumption-b".to_string()],
        validity: ValidityDomain::unconstrained()
            .with("Mach number", 0.2, 0.8)
            .with("α", -1.0, 1.0),
        discrepancy_rel: 0.125,
        in_domain: true,
    }
}

fn identified_model_evidence(model_evidence: ModelEvidence) -> IdentifiedModelEvidenceV1 {
    identify_model_evidence_v1(model_evidence, LIMITS, || false)
        .expect("valid model-evidence identity")
}

fn manual_model_evidence_receipt(model_evidence: &ModelEvidence) -> ModelEvidenceReceiptV1 {
    let validity = identified_domain(model_evidence.validity.clone());
    CanonicalEncoder::<ModelEvidenceIdV1, _>::new(LIMITS, || false)
        .expect("model-evidence schema")
        .canonical_set(
            Field::new(0, "model-card-names"),
            u64::try_from(model_evidence.cards.len()).expect("model-card count"),
            model_evidence.cards.iter().map(|card| card.as_bytes()),
        )
        .expect("model-card names")
        .canonical_set(
            Field::new(1, "assumptions"),
            u64::try_from(model_evidence.assumptions.len()).expect("assumption count"),
            model_evidence
                .assumptions
                .iter()
                .map(|assumption| assumption.as_bytes()),
        )
        .expect("assumptions")
        .child(Field::new(2, "validity"), validity.id())
        .expect("validity")
        .u64(
            Field::new(3, "discrepancy-rel-ieee754-bits"),
            model_evidence.discrepancy_rel.to_bits(),
        )
        .expect("discrepancy")
        .flag(Field::new(4, "in-domain"), model_evidence.in_domain)
        .expect("in-domain")
        .finish()
        .expect("manual model-evidence identity")
}

fn certified_fixture() -> Certified<f64> {
    let model = model_evidence_fixture();
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

fn certified_source(
    domain: &str,
    schema_version: u32,
    canonical_bytes: &[u8],
) -> CertifiedF64SourceV1 {
    identify_certified_f64_source_v1(
        domain.to_string(),
        schema_version,
        canonical_bytes.to_vec(),
        LIMITS,
        || false,
    )
    .expect("valid certified-f64 exact source")
}

fn sourced_certified(
    value: f64,
    producer_bytes: &[u8],
    adjoint_bytes: Option<&[u8]>,
) -> IdentifiedSourcedCertifiedF64EvidenceV1 {
    let mut evidence = Evidence::exact(value, ProvenanceHash::of_bytes(producer_bytes));
    if let Some(bytes) = adjoint_bytes {
        evidence = evidence.with_adjoint(ProvenanceHash::of_bytes(bytes));
    }
    let certified = evidence.certified().expect("valid exact evidence");
    let child = identified_certified(certified);
    let producer = certified_source("org.frankensim.tests.producer.v1", 1, producer_bytes);
    let adjoint =
        adjoint_bytes.map(|bytes| certified_source("org.frankensim.tests.adjoint.v1", 1, bytes));
    identify_sourced_certified_f64_evidence_v1(child, producer, adjoint, LIMITS, || false)
        .expect("matching exact sources")
}

fn manual_certified_f64_source_receipt(
    source: &CertifiedF64SourceV1,
) -> CertifiedF64SourceReceiptV1 {
    CanonicalEncoder::<CertifiedF64SourceIdV1, _>::new(LIMITS, || false)
        .expect("source schema")
        .utf8(Field::new(0, "source-domain"), source.domain())
        .expect("source domain")
        .u64(
            Field::new(1, "source-schema-version"),
            u64::from(source.schema_version()),
        )
        .expect("source schema version")
        .bytes(Field::new(2, "canonical-source"), source.canonical_bytes())
        .expect("canonical source")
        .finish()
        .expect("manual source identity")
}

fn manual_sourced_certified_f64_receipt(
    child: &IdentifiedCertifiedF64EvidenceV1,
    producer: &CertifiedF64SourceV1,
    adjoint: Option<&CertifiedF64SourceV1>,
) -> SourcedCertifiedF64EvidenceReceiptV1 {
    let adjoint_id = adjoint.map(CertifiedF64SourceV1::id);
    CanonicalEncoder::<SourcedCertifiedF64EvidenceIdV1, _>::new(LIMITS, || false)
        .expect("sourced parent schema")
        .child(Field::new(0, "certified-f64-evidence"), child.id())
        .expect("certified child")
        .child(Field::new(1, "producer-source"), producer.id())
        .expect("producer source")
        .flag(
            Field::new(2, "adjoint-source-present"),
            adjoint_id.is_some(),
        )
        .expect("adjoint presence")
        .ordered_children(
            Field::new(3, "adjoint-source"),
            u64::from(adjoint_id.is_some()),
            adjoint_id,
        )
        .expect("optional adjoint source")
        .finish()
        .expect("manual sourced parent")
}

fn identified_decision_assessment(
    certified: Certified<f64>,
    threshold_rel: f64,
) -> IdentifiedCertifiedF64DecisionAssessmentV1 {
    identify_certified_f64_decision_assessment_v1(
        identified_certified(certified),
        threshold_rel,
        LIMITS,
        || false,
    )
    .expect("valid certified-f64 decision assessment")
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

fn manual_decision_assessment_receipt(
    certified: &IdentifiedCertifiedF64EvidenceV1,
    threshold_rel: f64,
) -> CertifiedF64DecisionAssessmentReceiptV1 {
    let evidence = certified.certified();
    let breakdown = evidence.breakdown();
    let total_rel = breakdown.total_rel();
    let status = evidence.assess(threshold_rel);
    let advice = evidence.escalation_advice(threshold_rel);
    let source_tag = |source| match source {
        UncertaintySource::ModelForm => 1_u32,
        UncertaintySource::Statistical => 2,
        UncertaintySource::Numerical => 3,
    };
    let mut status_payload = [0_u8; 4];
    let (status_tag, status_payload_len) = match status {
        DecisionStatus::DecisionGrade => (1, 0),
        DecisionStatus::NotDecisionGrade { dominant, .. } => {
            status_payload.copy_from_slice(&source_tag(dominant).to_le_bytes());
            (2, status_payload.len())
        }
    };
    let advice_tag = match advice {
        EscalationAdvice::NoneNeeded => 1,
        EscalationAdvice::RefineNumerics => 2,
        EscalationAdvice::GatherMoreSamples => 3,
        EscalationAdvice::EscalateModelFidelity => 4,
    };
    CanonicalEncoder::<CertifiedF64DecisionAssessmentIdV1, _>::new(LIMITS, || false)
        .expect("decision-assessment schema")
        .child(Field::new(0, "certified-f64-evidence"), certified.id())
        .expect("certified child")
        .u64(
            Field::new(1, "assessment-algorithm-version"),
            u64::from(DECISION_ASSESSMENT_ALGORITHM_VERSION_V1),
        )
        .expect("assessment algorithm")
        .finite_f64(Field::new(2, "threshold-rel"), threshold_rel)
        .expect("threshold")
        .u64(
            Field::new(3, "numerical-rel-ieee754-bits"),
            breakdown.numerical_rel.to_bits(),
        )
        .expect("numerical band")
        .u64(
            Field::new(4, "statistical-rel-ieee754-bits"),
            breakdown.statistical_rel.to_bits(),
        )
        .expect("statistical band")
        .u64(
            Field::new(5, "model-rel-ieee754-bits"),
            breakdown.model_rel.to_bits(),
        )
        .expect("model band")
        .u64(Field::new(6, "total-rel-ieee754-bits"), total_rel.to_bits())
        .expect("total band")
        .variant(
            Field::new(7, "status"),
            status_tag,
            &status_payload[..status_payload_len],
        )
        .expect("status")
        .variant(Field::new(8, "advice"), advice_tag, &[])
        .expect("advice")
        .finish()
        .expect("manual decision-assessment identity")
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

fn card_bound_model_card(name: &str, version: &str) -> IdentifiedModelCardV1 {
    identified_model_card(
        ModelCard::new(
            name,
            version,
            Ambition::Solid,
            Vec::new(),
            ValidityDomain::unconstrained(),
            Vec::new(),
            0.0,
        ),
        None,
    )
}

fn card_bound_model_evidence_fixture() -> IdentifiedCardBoundModelEvidenceV1 {
    identify_card_bound_model_evidence_v1(
        identified_model_evidence(model_evidence_fixture()),
        vec![
            card_bound_model_card("card-a", "1.0.0"),
            card_bound_model_card("card-b", "2.0.0"),
        ],
        LIMITS,
        || false,
    )
    .expect("valid card-bound model-evidence fixture")
}

fn manual_card_bound_model_evidence_receipt(
    model_evidence: &IdentifiedModelEvidenceV1,
    model_cards: &[IdentifiedModelCardV1],
    limits: CanonicalLimits,
) -> CardBoundModelEvidenceReceiptV1 {
    CanonicalEncoder::<CardBoundModelEvidenceIdV1, _>::new(limits, || false)
        .expect("card-bound model-evidence schema")
        .child(Field::new(0, "model-evidence"), model_evidence.id())
        .expect("model-evidence child")
        .ordered_children(
            Field::new(1, "model-card-declarations"),
            u64::try_from(model_cards.len()).expect("model-card count"),
            model_cards.iter().map(IdentifiedModelCardV1::id),
        )
        .expect("model-card declarations")
        .finish()
        .expect("manual card-bound model-evidence identity")
}

fn card_bound_model_card_field_bytes(card_count: u64) -> u64 {
    let descriptor_bytes = [
        1_u64,
        u64::from(u64::BITS / 8),
        u64::try_from(ModelCardIdentitySchemaV1::DOMAIN.len()).expect("domain bytes"),
        u64::from(u64::BITS / 8),
        u64::try_from(ModelCardIdentitySchemaV1::NAME.len()).expect("name bytes"),
        32,
        u64::from(u32::BITS / 8),
        u64::from(u64::BITS / 8),
        u64::try_from(ModelCardIdentitySchemaV1::CONTEXT.len()).expect("context bytes"),
    ]
    .into_iter()
    .sum::<u64>();
    u64::from(u64::BITS / 8) + descriptor_bytes + card_count * 32
}

#[test]
fn standalone_certificate_identities_replay_manual_frames_and_retain_inputs() {
    let numerical_certificates = [
        NumericalCertificate::exact(-0.0),
        NumericalCertificate::exact(f64::INFINITY),
        NumericalCertificate::enclosure(f64::NEG_INFINITY, 4.0),
        NumericalCertificate::estimate(-2.0, f64::INFINITY),
        NumericalCertificate::no_claim(),
    ];
    for certificate in numerical_certificates {
        let first = identify_numerical_certificate_v1(certificate, LIMITS, || false)
            .expect("admitted numerical certificate");
        let replay = identify_numerical_certificate_v1(certificate, LIMITS, || false)
            .expect("replayed numerical certificate");
        let manual = manual_numerical_certificate_receipt(certificate);
        assert_eq!(first.id(), replay.id());
        assert_eq!(first.id(), manual.id());
        assert_eq!(
            first.receipt().canonical_preimage(),
            manual.canonical_preimage()
        );
        assert_eq!(first.certificate(), &certificate);
        assert_eq!(first.id_bytes(), first.receipt().audit_record().id());
        assert_eq!(first.trust_state(), TrustState::Unanchored);
        assert_eq!(
            first.receipt().audit_record().no_claim(),
            NoClaimState::ExternalTrustRequired
        );
        assert_eq!(first.into_certificate(), certificate);
    }

    let statistical_certificates = [
        StatisticalCertificate::None,
        StatisticalCertificate::EValue {
            e: -0.0,
            alpha: 0.05,
        },
        StatisticalCertificate::HalfWidth {
            half_width: -0.0,
            confidence: 0.95,
        },
    ];
    for certificate in statistical_certificates {
        let first = identify_statistical_certificate_v1(certificate, LIMITS, || false)
            .expect("admitted statistical certificate");
        let replay = identify_statistical_certificate_v1(certificate, LIMITS, || false)
            .expect("replayed statistical certificate");
        let manual = manual_statistical_certificate_receipt(certificate);
        assert_eq!(first.id(), replay.id());
        assert_eq!(first.id(), manual.id());
        assert_eq!(
            first.receipt().canonical_preimage(),
            manual.canonical_preimage()
        );
        assert_eq!(first.certificate(), &certificate);
        assert_eq!(first.id_bytes(), first.receipt().audit_record().id());
        assert_eq!(first.trust_state(), TrustState::Unanchored);
        assert_eq!(
            first.receipt().audit_record().no_claim(),
            NoClaimState::ExternalTrustRequired
        );
        assert_eq!(first.into_certificate(), certificate);
    }
}

#[test]
fn numerical_certificate_identity_binds_bits_and_refuses_malformed_shapes() {
    let enclosure = NumericalCertificate::enclosure(1.0, 2.0);
    let base =
        identify_numerical_certificate_v1(enclosure, LIMITS, || false).expect("baseline enclosure");
    let estimate = identify_numerical_certificate_v1(
        NumericalCertificate {
            kind: NumericalKind::Estimate,
            ..enclosure
        },
        LIMITS,
        || false,
    )
    .expect("same bounds with estimate tag");
    let lower = identify_numerical_certificate_v1(
        NumericalCertificate::enclosure(1.0_f64.next_down(), 2.0),
        LIMITS,
        || false,
    )
    .expect("mutated lower bound");
    let upper = identify_numerical_certificate_v1(
        NumericalCertificate::enclosure(1.0, 2.0_f64.next_up()),
        LIMITS,
        || false,
    )
    .expect("mutated upper bound");
    assert_ne!(base.id(), estimate.id());
    assert_ne!(base.id(), lower.id());
    assert_ne!(base.id(), upper.id());

    let positive_zero =
        identify_numerical_certificate_v1(NumericalCertificate::exact(0.0), LIMITS, || false)
            .expect("positive-zero exact");
    let negative_zero =
        identify_numerical_certificate_v1(NumericalCertificate::exact(-0.0), LIMITS, || false)
            .expect("negative-zero exact");
    assert_ne!(positive_zero.id(), negative_zero.id());

    let normalized = identify_numerical_certificate_v1(
        NumericalCertificate::enclosure(2.0, 1.0),
        LIMITS,
        || false,
    )
    .expect("constructor normalizes before identity admission");
    assert_eq!(normalized.id(), base.id());

    for certificate in [
        NumericalCertificate::exact(f64::NEG_INFINITY),
        NumericalCertificate::enclosure(f64::NEG_INFINITY, f64::INFINITY),
        NumericalCertificate::estimate(f64::NEG_INFINITY, f64::INFINITY),
    ] {
        identify_numerical_certificate_v1(certificate, LIMITS, || false)
            .expect("ordered infinite structural state is admitted");
    }

    let nan_bits = 0x7ff8_0000_0000_0042;
    let nan = f64::from_bits(nan_bits);
    let malformed = [
        (
            NumericalCertificate {
                kind: NumericalKind::Estimate,
                lo: nan,
                hi: 1.0,
            },
            "bounds must not be NaN",
        ),
        (
            NumericalCertificate {
                kind: NumericalKind::Enclosure,
                lo: 2.0,
                hi: 1.0,
            },
            "lower bound must not exceed upper bound",
        ),
        (
            NumericalCertificate {
                kind: NumericalKind::Exact,
                lo: 1.0,
                hi: 1.0_f64.next_up(),
            },
            "Exact bounds must be bit-identical",
        ),
        (
            NumericalCertificate {
                kind: NumericalKind::Exact,
                lo: -0.0,
                hi: 0.0,
            },
            "Exact bounds must be bit-identical",
        ),
        (
            NumericalCertificate {
                kind: NumericalKind::NoClaim,
                lo: -1.0,
                hi: 1.0,
            },
            "NoClaim must use the canonical negative-infinity to positive-infinity bounds",
        ),
    ];
    for (certificate, expected_reason) in malformed {
        assert!(matches!(
            identify_numerical_certificate_v1(certificate, LIMITS, || false),
            Err(NumericalCertificateIdentityError::InvalidShape {
                kind,
                lo_bits,
                hi_bits,
                reason,
            }) if kind == certificate.kind
                && lo_bits == certificate.lo.to_bits()
                && hi_bits == certificate.hi.to_bits()
                && reason == expected_reason
        ));
    }
}

#[test]
fn statistical_certificate_identity_binds_bits_and_refuses_invalid_parameters() {
    let none = identify_statistical_certificate_v1(StatisticalCertificate::None, LIMITS, || false)
        .expect("none state");
    let e_value = StatisticalCertificate::EValue {
        e: 8.0,
        alpha: 0.05,
    };
    let e_base =
        identify_statistical_certificate_v1(e_value, LIMITS, || false).expect("baseline e-value");
    let e_mutated = identify_statistical_certificate_v1(
        StatisticalCertificate::EValue {
            e: 8.0_f64.next_up(),
            alpha: 0.05,
        },
        LIMITS,
        || false,
    )
    .expect("mutated e-value");
    let alpha_mutated = identify_statistical_certificate_v1(
        StatisticalCertificate::EValue {
            e: 8.0,
            alpha: 0.05_f64.next_up(),
        },
        LIMITS,
        || false,
    )
    .expect("mutated alpha");
    let half_base = identify_statistical_certificate_v1(
        StatisticalCertificate::HalfWidth {
            half_width: 0.25,
            confidence: 0.95,
        },
        LIMITS,
        || false,
    )
    .expect("baseline half-width");
    let half_mutated = identify_statistical_certificate_v1(
        StatisticalCertificate::HalfWidth {
            half_width: 0.25_f64.next_up(),
            confidence: 0.95,
        },
        LIMITS,
        || false,
    )
    .expect("mutated half-width");
    let confidence_mutated = identify_statistical_certificate_v1(
        StatisticalCertificate::HalfWidth {
            half_width: 0.25,
            confidence: 0.95_f64.next_down(),
        },
        LIMITS,
        || false,
    )
    .expect("mutated confidence");
    for mutated in [
        e_base.id(),
        e_mutated.id(),
        alpha_mutated.id(),
        half_base.id(),
        half_mutated.id(),
        confidence_mutated.id(),
    ] {
        assert_ne!(none.id(), mutated);
    }
    assert_ne!(e_base.id(), e_mutated.id());
    assert_ne!(e_base.id(), alpha_mutated.id());
    assert_ne!(half_base.id(), half_mutated.id());
    assert_ne!(half_base.id(), confidence_mutated.id());

    let positive_e = identify_statistical_certificate_v1(
        StatisticalCertificate::EValue {
            e: 0.0,
            alpha: 0.05,
        },
        LIMITS,
        || false,
    )
    .expect("positive-zero e-value");
    let negative_e = identify_statistical_certificate_v1(
        StatisticalCertificate::EValue {
            e: -0.0,
            alpha: 0.05,
        },
        LIMITS,
        || false,
    )
    .expect("negative-zero e-value");
    assert_ne!(positive_e.id(), negative_e.id());

    let positive_width = identify_statistical_certificate_v1(
        StatisticalCertificate::HalfWidth {
            half_width: 0.0,
            confidence: 0.95,
        },
        LIMITS,
        || false,
    )
    .expect("positive-zero width");
    let negative_width = identify_statistical_certificate_v1(
        StatisticalCertificate::HalfWidth {
            half_width: -0.0,
            confidence: 0.95,
        },
        LIMITS,
        || false,
    )
    .expect("negative-zero width");
    assert_ne!(positive_width.id(), negative_width.id());

    let invalid = [
        (
            StatisticalCertificate::EValue {
                e: -1.0,
                alpha: 0.05,
            },
            "e",
            -1.0_f64,
        ),
        (
            StatisticalCertificate::EValue {
                e: f64::INFINITY,
                alpha: 0.05,
            },
            "e",
            f64::INFINITY,
        ),
        (
            StatisticalCertificate::EValue {
                e: f64::NAN,
                alpha: 0.05,
            },
            "e",
            f64::NAN,
        ),
        (
            StatisticalCertificate::EValue { e: 1.0, alpha: 0.0 },
            "alpha",
            0.0,
        ),
        (
            StatisticalCertificate::EValue {
                e: 1.0,
                alpha: -0.0,
            },
            "alpha",
            -0.0,
        ),
        (
            StatisticalCertificate::EValue { e: 1.0, alpha: 1.0 },
            "alpha",
            1.0,
        ),
        (
            StatisticalCertificate::EValue {
                e: 1.0,
                alpha: f64::INFINITY,
            },
            "alpha",
            f64::INFINITY,
        ),
        (
            StatisticalCertificate::EValue {
                e: 1.0,
                alpha: f64::NAN,
            },
            "alpha",
            f64::NAN,
        ),
        (
            StatisticalCertificate::HalfWidth {
                half_width: -1.0,
                confidence: 0.95,
            },
            "half_width",
            -1.0,
        ),
        (
            StatisticalCertificate::HalfWidth {
                half_width: f64::NAN,
                confidence: 0.95,
            },
            "half_width",
            f64::NAN,
        ),
        (
            StatisticalCertificate::HalfWidth {
                half_width: f64::INFINITY,
                confidence: 0.95,
            },
            "half_width",
            f64::INFINITY,
        ),
        (
            StatisticalCertificate::HalfWidth {
                half_width: 1.0,
                confidence: 0.0,
            },
            "confidence",
            0.0,
        ),
        (
            StatisticalCertificate::HalfWidth {
                half_width: 1.0,
                confidence: -0.0,
            },
            "confidence",
            -0.0,
        ),
        (
            StatisticalCertificate::HalfWidth {
                half_width: 1.0,
                confidence: f64::INFINITY,
            },
            "confidence",
            f64::INFINITY,
        ),
        (
            StatisticalCertificate::HalfWidth {
                half_width: 1.0,
                confidence: 1.0,
            },
            "confidence",
            1.0,
        ),
        (
            StatisticalCertificate::HalfWidth {
                half_width: 1.0,
                confidence: f64::NAN,
            },
            "confidence",
            f64::NAN,
        ),
    ];
    for (certificate, expected_field, expected_value) in invalid {
        assert!(matches!(
            identify_statistical_certificate_v1(certificate, LIMITS, || false),
            Err(StatisticalCertificateIdentityError::InvalidParameter {
                field,
                bits,
                ..
            }) if field == expected_field && bits == expected_value.to_bits()
        ));
    }
}

#[test]
fn standalone_certificate_identities_enforce_resources_cancellation_and_no_transitivity() {
    let numerical = NumericalCertificate::enclosure(-1.0, 2.0);
    let numerical_baseline =
        identify_numerical_certificate_v1(numerical, LIMITS, || false).expect("baseline numerical");
    let numerical_frame = numerical_baseline.receipt().canonical_bytes();
    identify_numerical_certificate_v1(
        numerical,
        CanonicalLimits::new(numerical_frame, 8_192, 3, 64, 64),
        || false,
    )
    .expect("exact numerical frame limit");
    assert!(matches!(
        identify_numerical_certificate_v1(
            numerical,
            CanonicalLimits::new(numerical_frame - 1, 8_192, 3, 64, 64),
            || false,
        ),
        Err(NumericalCertificateIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == numerical_frame - 1
    ));
    assert!(matches!(
        identify_numerical_certificate_v1(
            numerical,
            CanonicalLimits::new(16_384, 8_192, 2, 64, 64),
            || false,
        ),
        Err(NumericalCertificateIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::Fields,
                requested: 3,
                limit: 2,
            }
        ))
    ));

    let statistical = StatisticalCertificate::EValue {
        e: 4.0,
        alpha: 0.05,
    };
    let statistical_baseline = identify_statistical_certificate_v1(statistical, LIMITS, || false)
        .expect("baseline statistical");
    let statistical_frame = statistical_baseline.receipt().canonical_bytes();
    identify_statistical_certificate_v1(
        statistical,
        CanonicalLimits::new(statistical_frame, 8_192, 1, 64, 64),
        || false,
    )
    .expect("exact statistical frame limit");
    assert!(matches!(
        identify_statistical_certificate_v1(
            statistical,
            CanonicalLimits::new(statistical_frame - 1, 8_192, 1, 64, 64),
            || false,
        ),
        Err(StatisticalCertificateIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == statistical_frame - 1
    ));
    assert!(matches!(
        identify_statistical_certificate_v1(
            statistical,
            CanonicalLimits::new(16_384, 8_192, 0, 64, 64),
            || false,
        ),
        Err(StatisticalCertificateIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::Fields,
                requested: 1,
                limit: 0,
            }
        ))
    ));

    assert_eq!(
        identify_numerical_certificate_v1(
            numerical,
            CanonicalLimits::new(16_384, 8_192, 3, 64, 0),
            || false,
        ),
        Err(NumericalCertificateIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    );
    assert_eq!(
        identify_statistical_certificate_v1(
            statistical,
            CanonicalLimits::new(16_384, 8_192, 1, 64, 0),
            || false,
        ),
        Err(StatisticalCertificateIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    );
    assert!(matches!(
        identify_numerical_certificate_v1(numerical, LIMITS, || true),
        Err(NumericalCertificateIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes: 0 }
        ))
    ));
    assert!(matches!(
        identify_statistical_certificate_v1(statistical, LIMITS, || true),
        Err(StatisticalCertificateIdentityError::Canonical(
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
    let numerical_polls = std::cell::Cell::new(0_usize);
    identify_numerical_certificate_v1(numerical, LIMITS, || {
        numerical_polls.set(numerical_polls.get() + 1);
        false
    })
    .expect("baseline numerical poll count");
    assert!(matches!(
        identify_numerical_certificate_v1(
            numerical,
            LIMITS,
            CancelAfter {
                successful_polls: numerical_polls.get() - 1,
            },
        ),
        Err(NumericalCertificateIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));

    let statistical_polls = std::cell::Cell::new(0_usize);
    identify_statistical_certificate_v1(statistical, LIMITS, || {
        statistical_polls.set(statistical_polls.get() + 1);
        false
    })
    .expect("baseline statistical poll count");
    assert!(matches!(
        identify_statistical_certificate_v1(
            statistical,
            LIMITS,
            CancelAfter {
                successful_polls: statistical_polls.get() - 1,
            },
        ),
        Err(StatisticalCertificateIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));

    let mut first_evidence = Evidence::exact(1.0, ProvenanceHash(1));
    first_evidence.numerical = numerical;
    first_evidence.statistical = statistical;
    let mut second_evidence = Evidence::exact(9.0, ProvenanceHash(2));
    second_evidence.numerical = numerical;
    second_evidence.statistical = statistical;
    assert_ne!(first_evidence.qoi.to_bits(), second_evidence.qoi.to_bits());
    assert_eq!(
        identify_numerical_certificate_v1(first_evidence.numerical, LIMITS, || false)
            .expect("first structural projection")
            .id(),
        identify_numerical_certificate_v1(second_evidence.numerical, LIMITS, || false)
            .expect("second structural projection")
            .id()
    );
    assert_eq!(
        identify_statistical_certificate_v1(first_evidence.statistical, LIMITS, || false)
            .expect("first statistical projection")
            .id(),
        identify_statistical_certificate_v1(second_evidence.statistical, LIMITS, || false)
            .expect("second statistical projection")
            .id()
    );
}

#[test]
fn fidelity_pair_and_discrepancy_band_identities_replay_and_retain_inputs() {
    let pair = fidelity_pair_fixture();
    let first =
        identify_fidelity_pair_v1(pair.clone(), LIMITS, || false).expect("admitted fidelity pair");
    let replay =
        identify_fidelity_pair_v1(pair.clone(), LIMITS, || false).expect("replayed fidelity pair");
    let manual = manual_fidelity_pair_receipt(&pair);
    assert_eq!(first.id(), replay.id());
    assert_eq!(first.id(), manual.id());
    assert_eq!(
        first.receipt().canonical_preimage(),
        manual.canonical_preimage()
    );
    assert_eq!(first.pair(), &pair);
    assert_eq!(first.id_bytes(), first.receipt().audit_record().id());
    assert_eq!(first.trust_state(), TrustState::Unanchored);
    assert_eq!(
        first.receipt().audit_record().no_claim(),
        NoClaimState::ExternalTrustRequired
    );

    let mut reverse_insertion = BTreeMap::new();
    reverse_insertion.insert("mach".to_string(), 0.7);
    reverse_insertion.insert("alpha".to_string(), -0.0);
    let reordered = FidelityPair {
        params: reverse_insertion,
        lo_fi: pair.lo_fi,
        hi_fi: pair.hi_fi,
    };
    assert_eq!(
        first.id(),
        identify_fidelity_pair_v1(reordered, LIMITS, || false)
            .expect("BTreeMap-normalized pair")
            .id()
    );
    assert_eq!(first.into_pair(), pair);

    for band in [
        DiscrepancyBand {
            mean_rel: 0.1,
            max_rel: 0.25,
        },
        DiscrepancyBand {
            mean_rel: f64::INFINITY,
            max_rel: f64::INFINITY,
        },
    ] {
        let first = identify_discrepancy_band_v1(band, LIMITS, || false)
            .expect("admitted discrepancy band");
        let replay = identify_discrepancy_band_v1(band, LIMITS, || false)
            .expect("replayed discrepancy band");
        let manual = manual_discrepancy_band_receipt(band);
        assert_eq!(first.id(), replay.id());
        assert_eq!(first.id(), manual.id());
        assert_eq!(
            first.receipt().canonical_preimage(),
            manual.canonical_preimage()
        );
        assert_eq!(first.band(), &band);
        assert_eq!(first.id_bytes(), first.receipt().audit_record().id());
        assert_eq!(first.trust_state(), TrustState::Unanchored);
        assert_eq!(
            first.receipt().audit_record().no_claim(),
            NoClaimState::ExternalTrustRequired
        );
        assert_eq!(first.into_band(), band);
    }
}

#[test]
fn fidelity_pair_identity_binds_every_field_and_refuses_malformed_observations() {
    let pair = fidelity_pair_fixture();
    let base =
        identify_fidelity_pair_v1(pair.clone(), LIMITS, || false).expect("baseline fidelity pair");

    let mut coordinate = pair.clone();
    coordinate
        .params
        .insert("mach".to_string(), 0.7_f64.next_up());
    let mut renamed = pair.clone();
    let alpha = renamed
        .params
        .remove("alpha")
        .expect("fixture alpha coordinate");
    renamed.params.insert("beta".to_string(), alpha);
    let mut added = pair.clone();
    added.params.insert("reynolds".to_string(), 1_000.0);
    let mut lo_fi = pair.clone();
    lo_fi.lo_fi = 4.0_f64.next_up();
    let mut hi_fi = pair.clone();
    hi_fi.hi_fi = 5.0_f64.next_down();
    for (field, mutated) in [
        ("coordinate", coordinate),
        ("parameter-name", renamed),
        ("parameter-count", added),
        ("lo-fi-qoi", lo_fi),
        ("hi-fi-qoi", hi_fi),
    ] {
        assert_ne!(
            base.id(),
            identify_fidelity_pair_v1(mutated, LIMITS, || false)
                .expect("valid field mutation")
                .id(),
            "{field} must move the root"
        );
    }

    let mut positive_coordinate = pair.clone();
    positive_coordinate.params.insert("alpha".to_string(), 0.0);
    assert_ne!(
        base.id(),
        identify_fidelity_pair_v1(positive_coordinate, LIMITS, || false)
            .expect("positive-zero coordinate")
            .id()
    );
    let positive_qoi = FidelityPair {
        params: BTreeMap::from([("x".to_string(), 1.0)]),
        lo_fi: 0.0,
        hi_fi: 1.0,
    };
    let negative_qoi = FidelityPair {
        lo_fi: -0.0,
        ..positive_qoi.clone()
    };
    assert_ne!(
        identify_fidelity_pair_v1(positive_qoi, LIMITS, || false)
            .expect("positive-zero QoI")
            .id(),
        identify_fidelity_pair_v1(negative_qoi, LIMITS, || false)
            .expect("negative-zero QoI")
            .id()
    );

    assert_eq!(
        identify_fidelity_pair_v1(
            FidelityPair {
                params: BTreeMap::new(),
                lo_fi: 1.0,
                hi_fi: 2.0,
            },
            LIMITS,
            || false,
        ),
        Err(FidelityPairIdentityError::EmptyParameterPoint)
    );
    assert!(matches!(
        identify_fidelity_pair_v1(
            FidelityPair {
                params: BTreeMap::from([("alpha beta".to_string(), 1.0)]),
                lo_fi: 1.0,
                hi_fi: 2.0,
            },
            LIMITS,
            || false,
        ),
        Err(FidelityPairIdentityError::InvalidParameterName {
            parameter_index: 0,
            reason: "invalid-character",
        })
    ));
    let coordinate_nan = f64::from_bits(0x7ff8_0000_0000_0042);
    assert!(matches!(
        identify_fidelity_pair_v1(
            FidelityPair {
                params: BTreeMap::from([("x".to_string(), coordinate_nan)]),
                lo_fi: 1.0,
                hi_fi: 2.0,
            },
            LIMITS,
            || false,
        ),
        Err(FidelityPairIdentityError::NonFiniteParameter {
            parameter_index: 0,
            bits,
        }) if bits == coordinate_nan.to_bits()
    ));
    for (field, lo_fi, hi_fi, bits) in [
        ("lo_fi", f64::INFINITY, 2.0, f64::INFINITY.to_bits()),
        ("hi_fi", 1.0, f64::NAN, f64::NAN.to_bits()),
    ] {
        assert!(matches!(
            identify_fidelity_pair_v1(
                FidelityPair {
                    params: BTreeMap::from([("x".to_string(), 0.0)]),
                    lo_fi,
                    hi_fi,
                },
                LIMITS,
                || false,
            ),
            Err(FidelityPairIdentityError::NonFiniteQoi {
                field: actual_field,
                bits: actual_bits,
            }) if actual_field == field && actual_bits == bits
        ));
    }
}

#[test]
fn discrepancy_band_identity_binds_bits_and_refuses_malformed_bands() {
    let base_band = DiscrepancyBand {
        mean_rel: 0.1,
        max_rel: 0.25,
    };
    let base = identify_discrepancy_band_v1(base_band, LIMITS, || false)
        .expect("baseline discrepancy band");
    let mean = identify_discrepancy_band_v1(
        DiscrepancyBand {
            mean_rel: 0.1_f64.next_up(),
            ..base_band
        },
        LIMITS,
        || false,
    )
    .expect("mutated mean");
    let maximum = identify_discrepancy_band_v1(
        DiscrepancyBand {
            max_rel: 0.25_f64.next_up(),
            ..base_band
        },
        LIMITS,
        || false,
    )
    .expect("mutated maximum");
    assert_ne!(base.id(), mean.id());
    assert_ne!(base.id(), maximum.id());

    let positive_zero = identify_discrepancy_band_v1(
        DiscrepancyBand {
            mean_rel: 0.0,
            max_rel: 0.0,
        },
        LIMITS,
        || false,
    )
    .expect("positive-zero band");
    let negative_zero = identify_discrepancy_band_v1(
        DiscrepancyBand {
            mean_rel: -0.0,
            max_rel: 0.0,
        },
        LIMITS,
        || false,
    )
    .expect("negative-zero mean");
    assert_ne!(positive_zero.id(), negative_zero.id());
    identify_discrepancy_band_v1(
        DiscrepancyBand {
            mean_rel: 1.0,
            max_rel: f64::INFINITY,
        },
        LIMITS,
        || false,
    )
    .expect("positive infinity is explicit unbounded state");

    let invalid = [
        (
            DiscrepancyBand {
                mean_rel: f64::NAN,
                max_rel: 1.0,
            },
            "mean_rel",
            f64::NAN,
            "value must not be NaN",
        ),
        (
            DiscrepancyBand {
                mean_rel: 0.0,
                max_rel: f64::NAN,
            },
            "max_rel",
            f64::NAN,
            "value must not be NaN",
        ),
        (
            DiscrepancyBand {
                mean_rel: -1.0,
                max_rel: 1.0,
            },
            "mean_rel",
            -1.0,
            "value must be non-negative",
        ),
        (
            DiscrepancyBand {
                mean_rel: 0.0,
                max_rel: -1.0,
            },
            "max_rel",
            -1.0,
            "value must be non-negative",
        ),
        (
            DiscrepancyBand {
                mean_rel: 0.5,
                max_rel: 0.25,
            },
            "mean_rel",
            0.5,
            "mean_rel must not exceed max_rel",
        ),
    ];
    for (band, expected_field, expected_value, expected_reason) in invalid {
        assert!(matches!(
            identify_discrepancy_band_v1(band, LIMITS, || false),
            Err(DiscrepancyBandIdentityError::InvalidBand {
                field,
                bits,
                reason,
            }) if field == expected_field
                && bits == expected_value.to_bits()
                && reason == expected_reason
        ));
    }
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one resource matrix shares exact frame, field, and cancellation accounting"
)]
fn fidelity_pair_and_discrepancy_band_identities_enforce_resources_and_cancellation() {
    let pair = fidelity_pair_fixture();
    let pair_baseline =
        identify_fidelity_pair_v1(pair.clone(), LIMITS, || false).expect("baseline fidelity pair");
    let pair_frame = pair_baseline.receipt().canonical_bytes();
    identify_fidelity_pair_v1(
        pair.clone(),
        CanonicalLimits::new(pair_frame, 8_192, 3, 64, 64),
        || false,
    )
    .expect("exact fidelity-pair frame limit");
    assert!(matches!(
        identify_fidelity_pair_v1(
            pair.clone(),
            CanonicalLimits::new(pair_frame - 1, 8_192, 3, 64, 64),
            || false,
        ),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == pair_frame - 1
    ));
    assert!(matches!(
        identify_fidelity_pair_v1(
            pair.clone(),
            CanonicalLimits::new(16_384, 8_192, 2, 64, 64),
            || false,
        ),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::Fields,
                requested: 3,
                limit: 2,
            }
        ))
    ));

    let longest_name = "a".repeat(256);
    let longest_pair = FidelityPair {
        params: BTreeMap::from([(longest_name, 1.0)]),
        lo_fi: 1.0,
        hi_fi: 2.0,
    };
    identify_fidelity_pair_v1(
        longest_pair.clone(),
        CanonicalLimits::new(16_384, 288, 3, 3, 64),
        || false,
    )
    .expect("exact 288-byte parameter field");
    assert!(matches!(
        identify_fidelity_pair_v1(
            longest_pair,
            CanonicalLimits::new(16_384, 287, 3, 3, 64),
            || false,
        ),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 288,
                limit: 287,
            }
        ))
    ));

    assert!(matches!(
        identify_fidelity_pair_v1(
            pair.clone(),
            CanonicalLimits::new(16_384, 8_192, 3, 1, 64),
            || false,
        ),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CollectionItems,
                requested: 2,
                limit: 1,
            }
        ))
    ));
    let one_parameter = FidelityPair {
        params: BTreeMap::from([("x".to_string(), 1.0)]),
        lo_fi: 1.0,
        hi_fi: 2.0,
    };
    assert!(matches!(
        identify_fidelity_pair_v1(
            one_parameter,
            CanonicalLimits::new(16_384, 8_192, 3, 2, 64),
            || false,
        ),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::StreamChunks,
                requested: 3,
                limit: 2,
            }
        ))
    ));
    let oversized_params = (0..1_025)
        .map(|index| (format!("p{index:04}"), index as f64))
        .collect();
    assert!(matches!(
        identify_fidelity_pair_v1(
            FidelityPair {
                params: oversized_params,
                lo_fi: 1.0,
                hi_fi: 2.0,
            },
            CanonicalLimits::new(1 << 20, 1 << 20, 3, 4_000, 64),
            || false,
        ),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CollectionItems,
                requested: 1_025,
                limit: 1_024,
            }
        ))
    ));

    let band = DiscrepancyBand {
        mean_rel: 0.1,
        max_rel: 0.25,
    };
    let band_baseline =
        identify_discrepancy_band_v1(band, LIMITS, || false).expect("baseline discrepancy band");
    let band_frame = band_baseline.receipt().canonical_bytes();
    identify_discrepancy_band_v1(
        band,
        CanonicalLimits::new(band_frame, 8_192, 2, 64, 64),
        || false,
    )
    .expect("exact discrepancy-band frame limit");
    assert!(matches!(
        identify_discrepancy_band_v1(
            band,
            CanonicalLimits::new(band_frame - 1, 8_192, 2, 64, 64),
            || false,
        ),
        Err(DiscrepancyBandIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == band_frame - 1
    ));
    assert!(matches!(
        identify_discrepancy_band_v1(band, CanonicalLimits::new(16_384, 8_192, 1, 64, 64), || {
            false
        },),
        Err(DiscrepancyBandIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::Fields,
                requested: 2,
                limit: 1,
            }
        ))
    ));

    assert_eq!(
        identify_fidelity_pair_v1(
            pair.clone(),
            CanonicalLimits::new(16_384, 8_192, 3, 64, 0),
            || false,
        ),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    );
    assert_eq!(
        identify_discrepancy_band_v1(band, CanonicalLimits::new(16_384, 8_192, 2, 64, 0), || {
            false
        },),
        Err(DiscrepancyBandIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    );
    assert!(matches!(
        identify_fidelity_pair_v1(pair.clone(), LIMITS, || true),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes: 0 }
        ))
    ));
    assert!(matches!(
        identify_discrepancy_band_v1(band, LIMITS, || true),
        Err(DiscrepancyBandIdentityError::Canonical(
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
    assert!(matches!(
        identify_fidelity_pair_v1(
            pair.clone(),
            LIMITS,
            CancelAfter {
                successful_polls: 1,
            },
        ),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes: 0 }
        ))
    ));
    let pair_polls = std::cell::Cell::new(0_usize);
    identify_fidelity_pair_v1(pair.clone(), LIMITS, || {
        pair_polls.set(pair_polls.get() + 1);
        false
    })
    .expect("baseline fidelity-pair poll count");
    assert!(matches!(
        identify_fidelity_pair_v1(
            pair,
            LIMITS,
            CancelAfter {
                successful_polls: pair_polls.get() - 1,
            },
        ),
        Err(FidelityPairIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));
    let band_polls = std::cell::Cell::new(0_usize);
    identify_discrepancy_band_v1(band, LIMITS, || {
        band_polls.set(band_polls.get() + 1);
        false
    })
    .expect("baseline discrepancy-band poll count");
    assert!(matches!(
        identify_discrepancy_band_v1(
            band,
            LIMITS,
            CancelAfter {
                successful_polls: band_polls.get() - 1,
            },
        ),
        Err(DiscrepancyBandIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));
}

#[test]
fn model_bracket_identity_normalizes_construction_order_and_binds_every_member() {
    let members = [("model-z", 3.0), ("model-a", 1.0), ("model-m", 2.0)];
    let first_bracket = model_bracket(&members);
    let reordered_bracket = model_bracket(&[("model-m", 2.0), ("model-z", 3.0), ("model-a", 1.0)]);
    assert_eq!(first_bracket, reordered_bracket);

    let manual = manual_model_bracket_receipt(&members);
    let first = identified_model_bracket(first_bracket.clone());
    let replay = identified_model_bracket(first_bracket);
    let reordered = identified_model_bracket(reordered_bracket);
    assert_eq!(first.id(), replay.id());
    assert_eq!(first.id(), reordered.id());
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
    let evidence = first
        .bracket()
        .evidence(ProvenanceHash(0x1234))
        .expect("identified bracket remains usable");
    assert_eq!(
        evidence.model.cards,
        vec![
            "model-a".to_string(),
            "model-m".to_string(),
            "model-z".to_string()
        ]
    );
    let different_legacy_provenance = first
        .bracket()
        .evidence(ProvenanceHash(0x5678))
        .expect("same bracket with different legacy correlation");
    assert_ne!(evidence.provenance, different_legacy_provenance.provenance);
    assert_ne!(evidence, different_legacy_provenance);

    let card_v1 = ModelCard::new(
        "model-a",
        "1.0.0",
        Ambition::Solid,
        vec!["assumption-a".to_string()],
        ValidityDomain::unconstrained(),
        Vec::new(),
        0.1,
    );
    let card_v2 = ModelCard::new(
        "model-a",
        "2.0.0",
        Ambition::Moonshot,
        vec!["different-assumption".to_string()],
        ValidityDomain::unconstrained().with("x", 0.0, 1.0),
        vec!["different-failure".to_string()],
        0.9,
    );
    assert_ne!(card_v1, card_v2);
    let card_v1_bracket = identified_model_bracket(model_bracket(&[
        (card_v1.name.as_str(), 1.0),
        ("model-b", 2.0),
    ]));
    let card_v2_bracket = identified_model_bracket(model_bracket(&[
        (card_v2.name.as_str(), 1.0),
        ("model-b", 2.0),
    ]));
    assert_eq!(card_v1_bracket.id(), card_v2_bracket.id());

    let renamed = identified_model_bracket(model_bracket(&[
        ("model-z", 3.0),
        ("model-b", 1.0),
        ("model-m", 2.0),
    ]));
    let changed_qoi = identified_model_bracket(model_bracket(&[
        ("model-z", 3.0),
        ("model-a", 1.0_f64.next_up()),
        ("model-m", 2.0),
    ]));
    let added = identified_model_bracket(model_bracket(&[
        ("model-z", 3.0),
        ("model-a", 1.0),
        ("model-m", 2.0),
        ("model-x", 4.0),
    ]));
    assert_ne!(first.id(), renamed.id());
    assert_ne!(first.id(), changed_qoi.id());
    assert_ne!(first.id(), added.id());

    let positive_zero =
        identified_model_bracket(model_bracket(&[("model-a", 0.0), ("model-b", 1.0)]));
    let negative_zero =
        identified_model_bracket(model_bracket(&[("model-a", -0.0), ("model-b", 1.0)]));
    assert_ne!(positive_zero.id(), negative_zero.id());

    assert_eq!(replay.into_bracket(), model_bracket(&members));
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one model-bracket matrix shares exact resource and cancellation accounting"
)]
fn model_bracket_identity_enforces_admission_resources_and_cancellation() {
    assert_eq!(
        identify_model_bracket_v1(ModelBracket::new(), LIMITS, || false),
        Err(ModelBracketIdentityError::TooFewMembers {
            count: 0,
            minimum: 2,
        })
    );
    assert_eq!(
        identify_model_bracket_v1(model_bracket(&[("model-a", 1.0)]), LIMITS, || false),
        Err(ModelBracketIdentityError::TooFewMembers {
            count: 1,
            minimum: 2,
        })
    );

    let bracket = model_bracket(&[("model-a", 1.0), ("model-b", 2.0)]);
    let baseline = identified_model_bracket(bracket.clone());
    let frame = baseline.receipt().canonical_bytes();
    identify_model_bracket_v1(
        bracket.clone(),
        CanonicalLimits::new(frame, 8_192, 1, 6, 64),
        || false,
    )
    .expect("exact model-bracket frame, field, and chunk limits");
    assert!(matches!(
        identify_model_bracket_v1(
            bracket.clone(),
            CanonicalLimits::new(frame - 1, 8_192, 1, 6, 64),
            || false,
        ),
        Err(ModelBracketIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == frame - 1
    ));
    assert!(matches!(
        identify_model_bracket_v1(
            bracket.clone(),
            CanonicalLimits::new(16_384, 8_192, 0, 6, 64),
            || false,
        ),
        Err(ModelBracketIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::Fields,
                requested: 1,
                limit: 0,
            }
        ))
    ));
    assert!(matches!(
        identify_model_bracket_v1(
            model_bracket(&[("model-a", 1.0), ("model-b", 2.0), ("model-c", 3.0)]),
            CanonicalLimits::new(16_384, 8_192, 1, 2, 64),
            || false,
        ),
        Err(ModelBracketIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CollectionItems,
                requested: 3,
                limit: 2,
            }
        ))
    ));

    let maximum = usize::try_from(MAX_MODEL_BRACKET_MEMBERS_V1)
        .expect("model-bracket hard maximum fits usize");
    let maximal_bracket = (0..maximum)
        .try_fold(ModelBracket::new(), |bracket, index| {
            bracket.try_with_member(format!("model-{index:04}"), index as f64)
        })
        .expect("exact hard-maximum model bracket");
    identify_model_bracket_v1(
        maximal_bracket,
        CanonicalLimits::new(1 << 20, 1 << 20, 1, MAX_MODEL_BRACKET_MEMBERS_V1 * 3, 256),
        || false,
    )
    .expect("exact 1,024-member hard boundary");
    assert!(matches!(
        identify_model_bracket_v1(
            bracket.clone(),
            CanonicalLimits::new(16_384, 8_192, 1, 5, 64),
            || false,
        ),
        Err(ModelBracketIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::StreamChunks,
                requested: 6,
                limit: 5,
            }
        ))
    ));

    let long_a = "a".repeat(256);
    let long_b = "b".repeat(256);
    let long_bracket = model_bracket(&[(long_a.as_str(), 1.0), (long_b.as_str(), 2.0)]);
    identify_model_bracket_v1(
        long_bracket.clone(),
        CanonicalLimits::new(16_384, 568, 1, 6, 64),
        || false,
    )
    .expect("exact 568-byte model-bracket member field");
    assert!(matches!(
        identify_model_bracket_v1(
            long_bracket,
            CanonicalLimits::new(16_384, 567, 1, 6, 64),
            || false,
        ),
        Err(ModelBracketIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 568,
                limit: 567,
            }
        ))
    ));

    assert_eq!(
        identify_model_bracket_v1(
            bracket.clone(),
            CanonicalLimits::new(16_384, 8_192, 1, 6, 0),
            || false,
        ),
        Err(ModelBracketIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    );
    assert!(matches!(
        identify_model_bracket_v1(bracket.clone(), LIMITS, || true),
        Err(ModelBracketIdentityError::Canonical(
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
    assert!(matches!(
        identify_model_bracket_v1(
            bracket.clone(),
            LIMITS,
            CancelAfter {
                successful_polls: 1,
            },
        ),
        Err(ModelBracketIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes: 0 }
        ))
    ));
    let polls = std::cell::Cell::new(0_usize);
    identify_model_bracket_v1(bracket.clone(), LIMITS, || {
        polls.set(polls.get() + 1);
        false
    })
    .expect("baseline model-bracket poll count");
    assert!(matches!(
        identify_model_bracket_v1(
            bracket,
            LIMITS,
            CancelAfter {
                successful_polls: polls.get() - 1,
            },
        ),
        Err(ModelBracketIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));
}

fn decision_model_evidence(discrepancy_rel: f64) -> ModelEvidence {
    ModelEvidence {
        cards: vec!["model-a".to_string()],
        assumptions: vec!["local-model-assumption".to_string()],
        validity: ValidityDomain::unconstrained(),
        discrepancy_rel,
        in_domain: true,
    }
}

#[test]
fn certified_f64_decision_assessment_replays_and_covers_every_advice_lane() {
    let grade = Evidence::exact(1.0, ProvenanceHash(1))
        .certified()
        .expect("exact grade fixture");
    let numerical = Evidence::enclosed(1.0, 0.8, 1.2, ProvenanceHash(2))
        .certified()
        .expect("numerical fixture");
    let statistical = Evidence::exact(1.0, ProvenanceHash(3))
        .with_statistical(StatisticalCertificate::HalfWidth {
            half_width: 0.2,
            confidence: 0.95,
        })
        .certified()
        .expect("statistical fixture");
    let model = Evidence::exact(1.0, ProvenanceHash(4))
        .with_model(decision_model_evidence(0.2))
        .certified()
        .expect("model fixture");

    let cases = [
        (grade.clone(), 0.0, EscalationAdvice::NoneNeeded, None),
        (
            numerical.clone(),
            0.1,
            EscalationAdvice::RefineNumerics,
            Some(UncertaintySource::Numerical),
        ),
        (
            statistical.clone(),
            0.1,
            EscalationAdvice::GatherMoreSamples,
            Some(UncertaintySource::Statistical),
        ),
        (
            model.clone(),
            0.1,
            EscalationAdvice::EscalateModelFidelity,
            Some(UncertaintySource::ModelForm),
        ),
    ];
    for (certified, threshold, expected_advice, expected_dominant) in cases {
        let child = identified_certified(certified);
        let manual = manual_decision_assessment_receipt(&child, threshold);
        let first =
            identify_certified_f64_decision_assessment_v1(child.clone(), threshold, LIMITS, || {
                false
            })
            .expect("decision assessment");
        let replay =
            identify_certified_f64_decision_assessment_v1(child, threshold, LIMITS, || false)
                .expect("replayed decision assessment");
        assert_eq!(first.id(), replay.id());
        assert_eq!(first.id(), manual.id());
        assert_eq!(
            first.receipt().canonical_preimage(),
            manual.canonical_preimage()
        );
        assert_eq!(
            first.certified_evidence_id(),
            first.certified_evidence().id()
        );
        assert_eq!(first.threshold_rel().to_bits(), threshold.to_bits());
        assert_eq!(
            first.total_rel().to_bits(),
            first.breakdown().total_rel().to_bits()
        );
        assert_eq!(first.advice(), expected_advice);
        match (first.status(), expected_dominant) {
            (DecisionStatus::DecisionGrade, None) => {}
            (DecisionStatus::NotDecisionGrade { dominant, detail }, Some(expected_dominant)) => {
                assert_eq!(*dominant, expected_dominant);
                assert!(!detail.is_empty());
            }
            (status, dominant) => panic!("unexpected status {status:?} for {dominant:?}"),
        }
        assert_eq!(first.id_bytes(), first.receipt().audit_record().id());
        assert_eq!(first.trust_state(), TrustState::Unanchored);
        assert_eq!(
            first.receipt().audit_record().no_claim(),
            NoClaimState::ExternalTrustRequired
        );
        let (recovered_child, recovered_threshold) = first.into_inputs();
        assert_eq!(recovered_child.id(), replay.certified_evidence_id());
        assert_eq!(recovered_threshold.to_bits(), threshold.to_bits());
    }

    let tie = Evidence::enclosed(1.0, 0.8, 1.2, ProvenanceHash(5))
        .with_statistical(StatisticalCertificate::HalfWidth {
            half_width: 0.2,
            confidence: 0.95,
        })
        .with_model(decision_model_evidence(0.2))
        .certified()
        .expect("tie fixture");
    let tie = identified_decision_assessment(tie, 0.1);
    assert_eq!(tie.advice(), EscalationAdvice::EscalateModelFidelity);
    assert!(matches!(
        tie.status(),
        DecisionStatus::NotDecisionGrade {
            dominant: UncertaintySource::ModelForm,
            ..
        }
    ));

    let statistical_numerical_tie = Evidence::enclosed(1.0, 0.75, 1.25, ProvenanceHash(8))
        .with_statistical(StatisticalCertificate::HalfWidth {
            half_width: 0.25,
            confidence: 0.95,
        })
        .certified()
        .expect("statistical-numerical tie fixture");
    let statistical_numerical_tie = identified_decision_assessment(statistical_numerical_tie, 0.1);
    assert_eq!(
        statistical_numerical_tie.advice(),
        EscalationAdvice::GatherMoreSamples
    );
    assert!(matches!(
        statistical_numerical_tie.status(),
        DecisionStatus::NotDecisionGrade {
            dominant: UncertaintySource::Statistical,
            ..
        }
    ));

    let unbounded_numerical = Evidence::enclosed(1.0, -f64::MAX, f64::MAX, ProvenanceHash(9))
        .certified()
        .expect("unbounded relative numerical fixture");
    let unbounded_numerical = identified_decision_assessment(unbounded_numerical, 1.0);
    assert!(unbounded_numerical.breakdown().numerical_rel.is_infinite());
    assert!(unbounded_numerical.total_rel().is_infinite());
    assert_eq!(
        unbounded_numerical.advice(),
        EscalationAdvice::RefineNumerics
    );

    let unbounded_statistical = Evidence::exact(f64::MIN_POSITIVE, ProvenanceHash(6))
        .with_statistical(StatisticalCertificate::HalfWidth {
            half_width: f64::MAX,
            confidence: 0.95,
        })
        .certified()
        .expect("unbounded relative statistical fixture");
    let unbounded_statistical = identified_decision_assessment(unbounded_statistical, 1.0);
    assert!(
        unbounded_statistical
            .breakdown()
            .statistical_rel
            .is_infinite()
    );
    assert!(unbounded_statistical.total_rel().is_infinite());
    assert_eq!(
        unbounded_statistical.advice(),
        EscalationAdvice::GatherMoreSamples
    );

    let unbounded_model = Evidence::exact(1.0, ProvenanceHash(7))
        .with_model(decision_model_evidence(f64::INFINITY))
        .certified()
        .expect("unbounded model fixture");
    let unbounded_model = identified_decision_assessment(unbounded_model, 1.0);
    assert!(unbounded_model.breakdown().model_rel.is_infinite());
    assert_eq!(
        unbounded_model.advice(),
        EscalationAdvice::EscalateModelFidelity
    );
}

#[test]
fn certified_f64_decision_assessment_binds_child_threshold_and_derived_state() {
    let grade = Evidence::exact(1.0, ProvenanceHash(10))
        .certified()
        .expect("grade fixture");
    let positive_zero = identified_decision_assessment(grade.clone(), 0.0);
    let negative_zero = identified_decision_assessment(grade.clone(), -0.0);
    assert!(matches!(
        positive_zero.status(),
        DecisionStatus::DecisionGrade
    ));
    assert!(matches!(
        negative_zero.status(),
        DecisionStatus::DecisionGrade
    ));
    assert_ne!(positive_zero.id(), negative_zero.id());

    let looser = identified_decision_assessment(grade.clone(), 0.1);
    assert_ne!(positive_zero.id(), looser.id());

    let numerical = Evidence::enclosed(1.0, 0.9, 1.1, ProvenanceHash(11))
        .certified()
        .expect("numerical child fixture");
    let numerical = identified_decision_assessment(numerical, 0.0);
    assert_ne!(
        positive_zero.certified_evidence_id(),
        numerical.certified_evidence_id()
    );
    assert_ne!(positive_zero.id(), numerical.id());

    let provenance_a = Evidence::exact(1.0, ProvenanceHash(100))
        .certified()
        .expect("provenance-a fixture");
    let provenance_b = Evidence::exact(1.0, ProvenanceHash(200))
        .certified()
        .expect("provenance-b fixture");
    let provenance_a = identified_decision_assessment(provenance_a, 0.0);
    let provenance_b = identified_decision_assessment(provenance_b, 0.0);
    assert_eq!(
        provenance_a.certified_evidence_id(),
        provenance_b.certified_evidence_id()
    );
    assert_eq!(provenance_a.id(), provenance_b.id());
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one refusal matrix shares exact frame and cancellation accounting"
)]
fn certified_f64_decision_assessment_refuses_thresholds_resources_and_cancellation() {
    let certified = Evidence::exact(1.0, ProvenanceHash(20))
        .certified()
        .expect("decision resource fixture");
    for (threshold, reason) in [
        (f64::NAN, "threshold must be finite"),
        (f64::INFINITY, "threshold must be finite"),
        (f64::NEG_INFINITY, "threshold must be finite"),
        (-0.1, "threshold must be non-negative"),
    ] {
        let child = identified_certified(certified.clone());
        assert!(matches!(
            identify_certified_f64_decision_assessment_v1(
                child,
                threshold,
                LIMITS,
                || false,
            ),
            Err(CertifiedF64DecisionAssessmentIdentityError::InvalidThreshold {
                bits,
                reason: actual_reason,
            }) if bits == threshold.to_bits() && actual_reason == reason
        ));
    }

    let child = identified_certified(certified.clone());
    let baseline =
        identify_certified_f64_decision_assessment_v1(child.clone(), 0.0, LIMITS, || false)
            .expect("baseline decision assessment");
    let frame = baseline.receipt().canonical_bytes();
    identify_certified_f64_decision_assessment_v1(
        child.clone(),
        0.0,
        CanonicalLimits::new(frame, 8_192, 13, 64, 64),
        || false,
    )
    .expect("exact decision-assessment frame limit");
    assert!(matches!(
        identify_certified_f64_decision_assessment_v1(
            child.clone(),
            0.0,
            CanonicalLimits::new(frame - 1, 8_192, 13, 64, 64),
            || false,
        ),
        Err(CertifiedF64DecisionAssessmentIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == frame - 1
    ));
    assert!(matches!(
        identify_certified_f64_decision_assessment_v1(
            child.clone(),
            0.0,
            CanonicalLimits::new(16_384, 8_192, 8, 64, 64),
            || false,
        ),
        Err(CertifiedF64DecisionAssessmentIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::Fields,
                requested: 9,
                limit: 8,
            }
        ))
    ));
    assert!(matches!(
        identify_certified_f64_decision_assessment_v1(
            child.clone(),
            0.0,
            CanonicalLimits::new(16_384, 8_192, 13, 64, 0),
            || false,
        ),
        Err(CertifiedF64DecisionAssessmentIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    ));
    assert!(matches!(
        identify_certified_f64_decision_assessment_v1(child.clone(), 0.0, LIMITS, || true,),
        Err(CertifiedF64DecisionAssessmentIdentityError::Canonical(
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
    let polls = std::cell::Cell::new(0_usize);
    identify_certified_f64_decision_assessment_v1(child.clone(), 0.0, LIMITS, || {
        polls.set(polls.get() + 1);
        false
    })
    .expect("baseline decision-assessment poll count");
    assert!(matches!(
        identify_certified_f64_decision_assessment_v1(
            child,
            0.0,
            LIMITS,
            CancelAfter {
                successful_polls: polls.get() - 1,
            },
        ),
        Err(CertifiedF64DecisionAssessmentIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));
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
fn model_evidence_identity_replays_manual_frame_and_retains_input() {
    let model_evidence = model_evidence_fixture();
    let first = identified_model_evidence(model_evidence.clone());
    let replay = identified_model_evidence(model_evidence.clone());
    let manual = manual_model_evidence_receipt(&model_evidence);

    assert_eq!(first.id(), replay.id());
    assert_eq!(first.id(), manual.id());
    assert_eq!(
        first.receipt().canonical_preimage(),
        manual.canonical_preimage()
    );
    assert_eq!(
        first.validity_id(),
        identified_domain(model_evidence.validity.clone()).id()
    );
    assert_eq!(first.id_bytes(), first.receipt().audit_record().id());
    assert_eq!(first.trust_state(), TrustState::Unanchored);
    assert_eq!(
        first.receipt().audit_record().no_claim(),
        NoClaimState::ExternalTrustRequired
    );
    assert_eq!(first.model_evidence(), &model_evidence);

    let none = ModelEvidence::none();
    let identified_none = identified_model_evidence(none.clone());
    assert_eq!(
        identified_none.id(),
        manual_model_evidence_receipt(&none).id()
    );
    assert_ne!(first.id(), identified_none.id());

    assert_eq!(first.into_model_evidence(), model_evidence);
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one mutation and refusal matrix shares an exact model-evidence baseline"
)]
fn model_evidence_identity_binds_every_field_and_refuses_malformed_claims() {
    let model_evidence = model_evidence_fixture();
    let base = identified_model_evidence(model_evidence.clone());
    let identify = |model_evidence: ModelEvidence| {
        identify_model_evidence_v1(model_evidence, LIMITS, || false)
            .expect("valid model-evidence mutation")
    };

    let mut card = model_evidence.clone();
    card.cards.push("card-c".to_string());
    let card_id = identify(card).id();

    let mut assumption = model_evidence.clone();
    assumption.assumptions.push("assumption-c".to_string());
    let assumption_id = identify(assumption).id();

    let mut validity = model_evidence.clone();
    validity.validity = validity.validity.with("β", -2.0, 2.0);
    let validity_id = identify(validity).id();

    let mut discrepancy = model_evidence.clone();
    discrepancy.discrepancy_rel = 0.125_f64.next_up();
    let discrepancy_id = identify(discrepancy).id();

    let mut out_of_domain = model_evidence.clone();
    out_of_domain.in_domain = false;
    let out_of_domain_id = identify(out_of_domain).id();

    for (field, mutated_id) in [
        ("model-card", card_id),
        ("assumption", assumption_id),
        ("validity", validity_id),
        ("discrepancy", discrepancy_id),
        ("in-domain", out_of_domain_id),
    ] {
        assert_ne!(base.id(), mutated_id, "{field} must move the model root");
    }

    let mut unsorted = ModelEvidence::none();
    unsorted.cards = vec!["z".to_string(), "a".to_string()];
    assert!(matches!(
        identify_model_evidence_v1(unsorted, LIMITS, || false),
        Err(ModelEvidenceIdentityError::Canonical(
            CanonicalError::NonCanonicalSetOrder { index: 1 }
        ))
    ));

    let mut duplicate = ModelEvidence::none();
    duplicate.assumptions = vec!["same".to_string(), "same".to_string()];
    assert!(matches!(
        identify_model_evidence_v1(duplicate, LIMITS, || false),
        Err(ModelEvidenceIdentityError::Canonical(
            CanonicalError::DuplicateSetItem { index: 1 }
        ))
    ));

    let mut invalid_validity = ModelEvidence::none();
    invalid_validity.validity = invalid_validity.validity.with("bad", f64::NAN, 1.0);
    assert!(matches!(
        identify_model_evidence_v1(invalid_validity, LIMITS, || false),
        Err(ModelEvidenceIdentityError::Validity(
            ValidityDomainIdentityError::InvalidBounds { .. }
        ))
    ));

    for discrepancy_rel in [f64::NAN, -0.1, f64::NEG_INFINITY] {
        let mut invalid = ModelEvidence::none();
        invalid.discrepancy_rel = discrepancy_rel;
        assert!(matches!(
            identify_model_evidence_v1(invalid, LIMITS, || false),
            Err(ModelEvidenceIdentityError::InvalidDiscrepancy { bits, .. })
                if bits == discrepancy_rel.to_bits()
        ));
    }

    let mut positive_zero = ModelEvidence::none();
    positive_zero.discrepancy_rel = 0.0;
    let mut negative_zero = ModelEvidence::none();
    negative_zero.discrepancy_rel = -0.0;
    assert_ne!(identify(positive_zero).id(), identify(negative_zero).id());

    let mut unbounded = ModelEvidence::none();
    unbounded.discrepancy_rel = f64::INFINITY;
    let _unbounded = identify(unbounded);

    let none_id = identified_model_evidence(ModelEvidence::none()).id();
    let mut uncarded_assumption = ModelEvidence::none();
    uncarded_assumption.assumptions = vec!["declared assumption".to_string()];
    let mut uncarded_validity = ModelEvidence::none();
    uncarded_validity.validity = uncarded_validity.validity.with("Re", 1.0, 2.0);
    let mut uncarded_discrepancy = ModelEvidence::none();
    uncarded_discrepancy.discrepancy_rel = 0.25;
    let mut uncarded_out_of_domain = ModelEvidence::none();
    uncarded_out_of_domain.in_domain = false;
    for (field, diagnostic_id) in [
        ("assumption", identify(uncarded_assumption).id()),
        ("validity", identify(uncarded_validity).id()),
        ("discrepancy", identify(uncarded_discrepancy).id()),
        ("in-domain", identify(uncarded_out_of_domain).id()),
    ] {
        assert_ne!(
            none_id, diagnostic_id,
            "empty cards must not erase {field} diagnostic state"
        );
    }

    let original_card = model_card_fixture(None);
    let mut revised_card = original_card.clone();
    revised_card.version.push_str(".revision");
    revised_card.ambition = Ambition::Moonshot;
    revised_card
        .known_failures
        .push("new failure outside this projection".to_string());
    revised_card.calibration = Some(ProvenanceHash(0xfeed_face));
    let point = BTreeMap::from([("Mach number".to_string(), 0.5), ("α".to_string(), 0.0)]);
    let original_projection = ModelEvidence::from_card(&original_card, &point);
    let revised_projection = ModelEvidence::from_card(&revised_card, &point);
    assert_eq!(original_projection, revised_projection);
    assert_eq!(
        identify(original_projection).id(),
        identify(revised_projection).id(),
        "card fields absent from ModelEvidence must remain explicit no-claims"
    );
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one exact-limit matrix shares the same model-evidence schema and poll ledger"
)]
fn model_evidence_identity_enforces_exact_resources_and_cancellation() {
    let field_limits = CanonicalLimits::new(65_536, 1_024, 32, 64, 64);
    let mut exact_card = ModelEvidence::none();
    exact_card.cards = vec!["c".repeat(1_008)];
    identify_model_evidence_v1(exact_card, field_limits, || false)
        .expect("exact 1024-byte canonical-set payload is admitted");
    let mut oversized_card = ModelEvidence::none();
    oversized_card.cards = vec!["c".repeat(1_009)];
    assert!(matches!(
        identify_model_evidence_v1(oversized_card, field_limits, || false),
        Err(ModelEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 1_025,
                limit: 1_024,
            }
        ))
    ));

    let mut two_cards = ModelEvidence::none();
    two_cards.cards = vec!["a".to_string(), "b".to_string()];
    let collection_limits = CanonicalLimits::new(16_384, 8_192, 32, 1, 64);
    assert!(matches!(
        identify_model_evidence_v1(two_cards, collection_limits, || false),
        Err(ModelEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CollectionItems,
                requested: 2,
                limit: 1,
            }
        ))
    ));

    let model_evidence = model_evidence_fixture();
    let baseline = identified_model_evidence(model_evidence.clone());
    let frame_limit = baseline
        .receipt()
        .canonical_bytes()
        .checked_sub(1)
        .expect("non-empty model-evidence frame");
    let frame_limits = CanonicalLimits::new(frame_limit, 8_192, 32, 64, 64);
    assert!(matches!(
        identify_model_evidence_v1(model_evidence.clone(), frame_limits, || false),
        Err(ModelEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == frame_limit
    ));

    assert!(matches!(
        identify_model_evidence_v1(
            model_evidence.clone(),
            CanonicalLimits::new(16_384, 8_192, 32, 64, 0),
            || false,
        ),
        Err(ModelEvidenceIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    ));
    assert!(matches!(
        identify_model_evidence_v1(model_evidence.clone(), LIMITS, || true),
        Err(ModelEvidenceIdentityError::Canonical(
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
    let poll_count = std::cell::Cell::new(0_usize);
    identify_model_evidence_v1(model_evidence.clone(), LIMITS, || {
        poll_count.set(poll_count.get() + 1);
        false
    })
    .expect("baseline model-evidence poll count");
    let late = identify_model_evidence_v1(
        model_evidence,
        LIMITS,
        CancelAfter {
            successful_polls: poll_count.get() - 1,
        },
    );
    assert!(matches!(
        late,
        Err(ModelEvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));
}

#[test]
fn card_bound_model_evidence_replays_manual_frame_and_retains_children() {
    let first = card_bound_model_evidence_fixture();
    let replay = card_bound_model_evidence_fixture();
    let manual = manual_card_bound_model_evidence_receipt(
        first.model_evidence(),
        first.model_cards(),
        LIMITS,
    );

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
        first.model_evidence().model_evidence(),
        &model_evidence_fixture()
    );
    assert_eq!(
        first
            .model_cards()
            .iter()
            .map(|card| card.card().name.as_str())
            .collect::<Vec<_>>(),
        ["card-a", "card-b"]
    );

    let none_child = identified_model_evidence(ModelEvidence::none());
    let none = identify_card_bound_model_evidence_v1(none_child, Vec::new(), LIMITS, || false)
        .expect("empty evidence binds to an empty card vector");
    assert_eq!(
        none.id(),
        manual_card_bound_model_evidence_receipt(
            none.model_evidence(),
            none.model_cards(),
            LIMITS,
        )
        .id()
    );

    let mut diagnostic = ModelEvidence::none();
    diagnostic.assumptions = vec!["diagnostic without a card".to_string()];
    diagnostic.validity = diagnostic.validity.with("Re", 1.0, 2.0);
    diagnostic.discrepancy_rel = 0.25;
    diagnostic.in_domain = false;
    let diagnostic = identify_card_bound_model_evidence_v1(
        identified_model_evidence(diagnostic),
        Vec::new(),
        LIMITS,
        || false,
    )
    .expect("card-free diagnostic model evidence remains bindable");
    assert_ne!(none.id(), diagnostic.id());

    let expected_evidence_id = first.model_evidence().id();
    let expected_card_ids = first
        .model_cards()
        .iter()
        .map(IdentifiedModelCardV1::id)
        .collect::<Vec<_>>();
    let (model_evidence, model_cards) = first.into_parts();
    assert_eq!(model_evidence.id(), expected_evidence_id);
    assert_eq!(
        model_cards
            .iter()
            .map(IdentifiedModelCardV1::id)
            .collect::<Vec<_>>(),
        expected_card_ids
    );
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one adversarial matrix shares the exact positional card-binding baseline"
)]
fn card_bound_model_evidence_binds_card_content_and_refuses_name_mismatch() {
    let base = card_bound_model_evidence_fixture();
    let changed_card = identified_model_card(
        ModelCard::new(
            "card-a",
            "9.9.9",
            Ambition::Moonshot,
            vec!["contradicts the evidence assumptions".to_string()],
            ValidityDomain::unconstrained().with("outside-evidence-axis", 100.0, 200.0),
            vec!["known incompatible regime".to_string()],
            0.875,
        ),
        None,
    );
    let contradictory = identify_card_bound_model_evidence_v1(
        identified_model_evidence(model_evidence_fixture()),
        vec![changed_card, card_bound_model_card("card-b", "2.0.0")],
        LIMITS,
        || false,
    )
    .expect("the join checks names, not semantic agreement");
    assert_eq!(
        base.model_evidence().id(),
        contradictory.model_evidence().id()
    );
    assert_ne!(
        base.model_cards()[0].id(),
        contradictory.model_cards()[0].id()
    );
    assert_ne!(base.id(), contradictory.id());

    let bind_names = |names: &[&str]| {
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(model_evidence_fixture()),
            names
                .iter()
                .map(|name| card_bound_model_card(name, "test"))
                .collect(),
            LIMITS,
            || false,
        )
    };
    assert_eq!(
        bind_names(&["card-a"]).expect_err("one missing declaration must refuse"),
        CardBoundModelEvidenceIdentityError::CardCountMismatch {
            declared: 2,
            supplied: 1,
        }
    );
    assert_eq!(
        bind_names(&["card-a", "card-b", "card-c"]).expect_err("one extra declaration must refuse"),
        CardBoundModelEvidenceIdentityError::CardCountMismatch {
            declared: 2,
            supplied: 3,
        }
    );

    for (names, expected_index, declared_bytes, supplied_bytes, first_difference) in [
        (&["card-b", "card-a"][..], 0, 6, 6, 5),
        (&["card-a", "card-a"][..], 1, 6, 6, 5),
        (&["card-a-long", "card-b"][..], 0, 6, 11, 6),
        (&["Card-a", "card-b"][..], 0, 6, 6, 0),
        (&["card-a ", "card-b"][..], 0, 6, 7, 6),
    ] {
        assert_eq!(
            bind_names(names).expect_err("byte-inexact positional name must refuse"),
            CardBoundModelEvidenceIdentityError::CardNameMismatch {
                card_index: expected_index,
                declared_bytes,
                supplied_bytes,
                first_different_byte: first_difference,
            }
        );
    }

    let unicode_evidence = ModelEvidence {
        cards: vec!["café".to_string()],
        ..ModelEvidence::none()
    };
    assert_eq!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(unicode_evidence),
            vec![card_bound_model_card("cafe\u{301}", "test")],
            LIMITS,
            || false,
        )
        .expect_err("Unicode-normalization differences must remain visible"),
        CardBoundModelEvidenceIdentityError::CardNameMismatch {
            card_index: 0,
            declared_bytes: 5,
            supplied_bytes: 6,
            first_different_byte: 3,
        }
    );

    let long_prefix = "p".repeat(300);
    let late_declared = format!("{long_prefix}x");
    let late_supplied = format!("{long_prefix}y");
    let late_mismatch_evidence = ModelEvidence {
        cards: vec![late_declared],
        ..ModelEvidence::none()
    };
    assert_eq!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(late_mismatch_evidence),
            vec![card_bound_model_card(&late_supplied, "test")],
            LIMITS,
            || false,
        )
        .expect_err("a mismatch after one complete poll stride must retain its absolute offset"),
        CardBoundModelEvidenceIdentityError::CardNameMismatch {
            card_index: 0,
            declared_bytes: 301,
            supplied_bytes: 301,
            first_different_byte: 300,
        }
    );

    let long_prefix_evidence = ModelEvidence {
        cards: vec![format!("{long_prefix}x")],
        ..ModelEvidence::none()
    };
    assert_eq!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(long_prefix_evidence),
            vec![card_bound_model_card(&long_prefix, "test")],
            LIMITS,
            || false,
        )
        .expect_err("a strict-prefix mismatch after one poll stride must retain its offset"),
        CardBoundModelEvidenceIdentityError::CardNameMismatch {
            card_index: 0,
            declared_bytes: 301,
            supplied_bytes: 300,
            first_different_byte: 300,
        }
    );

    let model_evidence = identified_model_evidence(model_evidence_fixture());
    let reversed_cards = vec![
        card_bound_model_card("card-b", "2.0.0"),
        card_bound_model_card("card-a", "1.0.0"),
    ];
    let raw_mismatched =
        manual_card_bound_model_evidence_receipt(&model_evidence, &reversed_cards, LIMITS);
    assert_eq!(
        raw_mismatched.audit_record().trust(),
        TrustState::Unanchored
    );
    assert_eq!(
        raw_mismatched.audit_record().no_claim(),
        NoClaimState::ExternalTrustRequired
    );
    assert_ne!(raw_mismatched.id(), base.id());
    assert!(matches!(
        identify_card_bound_model_evidence_v1(
            model_evidence.clone(),
            reversed_cards,
            LIMITS,
            || false,
        ),
        Err(CardBoundModelEvidenceIdentityError::CardNameMismatch { card_index: 0, .. })
    ));

    let wrong_role =
        SourceId::<ModelEvidenceIdentitySchemaV1>::parse_slice(model_evidence.id().as_bytes())
            .expect("same digest bytes under an intentionally wrong source role");
    let role_refusal = CanonicalEncoder::<CardBoundModelEvidenceIdV1, _>::new(LIMITS, || false)
        .expect("card-bound schema")
        .child(Field::new(0, "model-evidence"), wrong_role)
        .expect_err("the model-evidence child requires the semantic role");
    assert!(matches!(
        role_refusal,
        CanonicalError::ChildBindingMismatch {
            field: "model-evidence",
            what: "child role",
        }
    ));
    let wrong_schema =
        SemanticId::<ModelCardIdentitySchemaV1>::parse_slice(base.model_cards()[0].id().as_bytes())
            .expect("same digest bytes under an intentionally wrong semantic schema");
    let schema_refusal = CanonicalEncoder::<CardBoundModelEvidenceIdV1, _>::new(LIMITS, || false)
        .expect("card-bound schema")
        .child(Field::new(0, "model-evidence"), wrong_schema)
        .expect_err("the model-evidence child requires its exact schema");
    assert!(matches!(
        schema_refusal,
        CanonicalError::ChildBindingMismatch {
            field: "model-evidence",
            what: "child schema domain",
        }
    ));
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one exact-limit matrix shares the card-bound schema and poll ledger"
)]
fn card_bound_model_evidence_enforces_exact_resources_and_cancellation() {
    let baseline = card_bound_model_evidence_fixture();
    let exact_field_bytes = card_bound_model_card_field_bytes(2);
    let exact_field_limits = CanonicalLimits::new(16_384, exact_field_bytes, 32, 2, 256);
    identify_card_bound_model_evidence_v1(
        identified_model_evidence(model_evidence_fixture()),
        vec![
            card_bound_model_card("card-a", "1.0.0"),
            card_bound_model_card("card-b", "2.0.0"),
        ],
        exact_field_limits,
        || false,
    )
    .expect("exact ordered-card field and collection budgets");
    assert!(matches!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(model_evidence_fixture()),
            vec![
                card_bound_model_card("card-a", "1.0.0"),
                card_bound_model_card("card-b", "2.0.0"),
            ],
            CanonicalLimits::new(16_384, exact_field_bytes - 1, 32, 2, 256),
            || false,
        ),
        Err(CardBoundModelEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested,
                limit,
            }
        )) if requested == exact_field_bytes && limit == exact_field_bytes - 1
    ));
    assert!(matches!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(model_evidence_fixture()),
            vec![
                card_bound_model_card("card-a", "1.0.0"),
                card_bound_model_card("card-b", "2.0.0"),
            ],
            CanonicalLimits::new(16_384, 8_192, 32, 1, 256),
            || false,
        ),
        Err(CardBoundModelEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CollectionItems,
                requested: 2,
                limit: 1,
            }
        ))
    ));
    assert!(matches!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(model_evidence_fixture()),
            vec![
                card_bound_model_card("card-a", "1.0.0"),
                card_bound_model_card("card-b", "2.0.0"),
            ],
            CanonicalLimits::new(16_384, 8_192, 8, 2, 256),
            || false,
        ),
        Err(CardBoundModelEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::Fields,
                requested,
                limit: 8,
            }
        )) if requested > 8
    ));

    let frame = baseline.receipt().canonical_bytes();
    identify_card_bound_model_evidence_v1(
        identified_model_evidence(model_evidence_fixture()),
        vec![
            card_bound_model_card("card-a", "1.0.0"),
            card_bound_model_card("card-b", "2.0.0"),
        ],
        CanonicalLimits::new(frame, 8_192, 32, 2, 256),
        || false,
    )
    .expect("exact complete-frame budget");
    assert!(matches!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(model_evidence_fixture()),
            vec![
                card_bound_model_card("card-a", "1.0.0"),
                card_bound_model_card("card-b", "2.0.0"),
            ],
            CanonicalLimits::new(frame - 1, 8_192, 32, 2, 256),
            || false,
        ),
        Err(CardBoundModelEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == frame - 1
    ));
    assert_eq!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(model_evidence_fixture()),
            vec![
                card_bound_model_card("card-a", "1.0.0"),
                card_bound_model_card("card-b", "2.0.0"),
            ],
            CanonicalLimits::new(16_384, 8_192, 32, 2, 0),
            || false,
        )
        .expect_err("zero cancellation stride must refuse"),
        CardBoundModelEvidenceIdentityError::Canonical(CanonicalError::InvalidLimits(
            "cancellation_poll_bytes must be positive"
        ))
    );
    assert!(matches!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(model_evidence_fixture()),
            vec![
                card_bound_model_card("card-a", "1.0.0"),
                card_bound_model_card("card-b", "2.0.0"),
            ],
            LIMITS,
            || true,
        ),
        Err(CardBoundModelEvidenceIdentityError::Canonical(
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

    let long_name = "l".repeat(600);
    let long_evidence = ModelEvidence {
        cards: vec![long_name.clone()],
        ..ModelEvidence::none()
    };
    assert_eq!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(long_evidence),
            vec![card_bound_model_card(&long_name, "1.0.0")],
            CanonicalLimits::new(16_384, 8_192, 32, 1, 256),
            CancelAfter {
                successful_polls: 2,
            },
        )
        .expect_err("long-name comparison must observe cancellation"),
        CardBoundModelEvidenceIdentityError::CardNameComparisonCancelled {
            card_index: 0,
            compared_bytes: 256,
        }
    );
    assert_eq!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(model_evidence_fixture()),
            vec![
                card_bound_model_card("card-a", "1.0.0"),
                card_bound_model_card("card-b", "2.0.0"),
            ],
            LIMITS,
            CancelAfter {
                successful_polls: 3,
            },
        )
        .expect_err("cancellation between positional card comparisons must refuse"),
        CardBoundModelEvidenceIdentityError::CardNameComparisonCancelled {
            card_index: 1,
            compared_bytes: 0,
        }
    );

    let polls = std::cell::Cell::new(0_usize);
    identify_card_bound_model_evidence_v1(
        identified_model_evidence(model_evidence_fixture()),
        vec![
            card_bound_model_card("card-a", "1.0.0"),
            card_bound_model_card("card-b", "2.0.0"),
        ],
        LIMITS,
        || {
            polls.set(polls.get() + 1);
            false
        },
    )
    .expect("baseline card-bound poll count");
    assert!(polls.get() > 1);
    assert!(matches!(
        identify_card_bound_model_evidence_v1(
            identified_model_evidence(model_evidence_fixture()),
            vec![
                card_bound_model_card("card-a", "1.0.0"),
                card_bound_model_card("card-b", "2.0.0"),
            ],
            LIMITS,
            CancelAfter {
                successful_polls: polls.get() - 1,
            },
        ),
        Err(CardBoundModelEvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { absorbed_bytes }
        )) if absorbed_bytes > 0
    ));
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
fn certified_f64_source_identity_binds_descriptor_bytes_and_retains_inputs() {
    let source = certified_source("org.frankensim.tests.solver-output.v7", 7, b"payload\0v1");
    let replay = certified_source("org.frankensim.tests.solver-output.v7", 7, b"payload\0v1");
    let manual = manual_certified_f64_source_receipt(&source);

    assert_eq!(source.id(), replay.id());
    assert_eq!(source.id(), manual.id());
    assert_eq!(
        source.receipt().canonical_preimage(),
        manual.canonical_preimage()
    );
    assert_ne!(
        source.id(),
        certified_source("org.frankensim.tests.solver-output.v8", 7, b"payload\0v1").id()
    );
    assert_ne!(
        source.id(),
        certified_source("org.frankensim.tests.solver-output.v7", 8, b"payload\0v1").id()
    );
    assert_ne!(
        source.id(),
        certified_source("org.frankensim.tests.solver-output.v7", 7, b"payload\0v2").id()
    );
    assert_eq!(source.trust_state(), TrustState::Unanchored);
    assert_eq!(
        source.receipt().audit_record().no_claim(),
        NoClaimState::ExternalTrustRequired
    );

    let (domain, version, bytes) = source.into_parts();
    assert_eq!(domain, "org.frankensim.tests.solver-output.v7");
    assert_eq!(version, 7);
    assert_eq!(bytes.as_slice(), b"payload\0v1");

    let empty = certified_source("org.frankensim.tests.empty.v1", 1, b"");
    assert!(empty.canonical_bytes().is_empty());
    assert!(matches!(
        identify_certified_f64_source_v1(String::new(), 1, Vec::new(), LIMITS, || false,),
        Err(CertifiedF64SourceIdentityError::EmptySourceDomain)
    ));
    assert!(matches!(
        identify_certified_f64_source_v1(
            "org.frankensim.tests.source".to_string(),
            0,
            Vec::new(),
            LIMITS,
            || false,
        ),
        Err(CertifiedF64SourceIdentityError::ZeroSourceSchemaVersion)
    ));
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one G0/G3 matrix keeps source crosswalk, collision, and mutation semantics together"
)]
fn sourced_certified_f64_identity_crosswalks_exact_sources_and_separates_fnv_collisions() {
    let parent = sourced_certified(3.0, b"producer-v1", Some(b"adjoint-v1"));
    let manual = manual_sourced_certified_f64_receipt(
        parent.certified_evidence(),
        parent.producer_source(),
        parent.adjoint_source(),
    );
    assert_eq!(parent.id(), manual.id());
    assert_eq!(
        parent.receipt().canonical_preimage(),
        manual.canonical_preimage()
    );
    assert_eq!(parent.producer_source().canonical_bytes(), b"producer-v1");
    assert_eq!(
        parent
            .adjoint_source()
            .expect("present adjoint source")
            .canonical_bytes(),
        b"adjoint-v1"
    );
    assert_eq!(parent.trust_state(), TrustState::Unanchored);
    assert_eq!(
        parent.receipt().audit_record().no_claim(),
        NoClaimState::ExternalTrustRequired
    );

    let base_id = parent.id();
    let producer_bytes = b"producer-v1";
    let adjoint_bytes = b"adjoint-v1";
    let make_child = |value| {
        identified_certified(
            Evidence::exact(value, ProvenanceHash::of_bytes(producer_bytes))
                .with_adjoint(ProvenanceHash::of_bytes(adjoint_bytes))
                .certified()
                .expect("valid sourced fixture"),
        )
    };
    let identify = |child: IdentifiedCertifiedF64EvidenceV1,
                    producer: CertifiedF64SourceV1,
                    adjoint: Option<CertifiedF64SourceV1>| {
        identify_sourced_certified_f64_evidence_v1(child, producer, adjoint, LIMITS, || false)
            .expect("matching source mutation")
    };

    let producer_domain_changed = identify(
        make_child(3.0),
        certified_source("org.frankensim.tests.producer-v2", 1, producer_bytes),
        Some(certified_source(
            "org.frankensim.tests.adjoint.v1",
            1,
            adjoint_bytes,
        )),
    );
    let producer_version_changed = identify(
        make_child(3.0),
        certified_source("org.frankensim.tests.producer.v1", 2, producer_bytes),
        Some(certified_source(
            "org.frankensim.tests.adjoint.v1",
            1,
            adjoint_bytes,
        )),
    );
    let adjoint_domain_changed = identify(
        make_child(3.0),
        certified_source("org.frankensim.tests.producer.v1", 1, producer_bytes),
        Some(certified_source(
            "org.frankensim.tests.adjoint-v2",
            1,
            adjoint_bytes,
        )),
    );
    let child_changed = identify(
        make_child(4.0),
        certified_source("org.frankensim.tests.producer.v1", 1, producer_bytes),
        Some(certified_source(
            "org.frankensim.tests.adjoint.v1",
            1,
            adjoint_bytes,
        )),
    );
    for (field, mutation) in [
        ("producer-domain", producer_domain_changed.id()),
        ("producer-schema-version", producer_version_changed.id()),
        ("adjoint-domain", adjoint_domain_changed.id()),
        ("certified-child", child_changed.id()),
    ] {
        assert_ne!(base_id, mutation, "{field} must move the parent root");
    }

    let collision_hash = ProvenanceHash::of_bytes(FNV1A64_ZERO_PREIMAGE);
    assert_eq!(
        collision_hash,
        ProvenanceHash::of_bytes(FNV1A64_ZERO_EXTENSION_COLLISION),
        "fixed fixture must remain a real legacy FNV-1a-64 collision"
    );
    let collision_child = || {
        identified_certified(
            Evidence::exact(9.0, collision_hash)
                .certified()
                .expect("collision fixture remains certifiable"),
        )
    };
    let first_source = certified_source(
        "org.frankensim.tests.collision.v1",
        1,
        FNV1A64_ZERO_PREIMAGE,
    );
    let second_source = certified_source(
        "org.frankensim.tests.collision.v1",
        1,
        FNV1A64_ZERO_EXTENSION_COLLISION,
    );
    let first_source_id = first_source.id();
    let second_source_id = second_source.id();
    let first = identify(collision_child(), first_source, None);
    let second = identify(collision_child(), second_source, None);
    assert_ne!(first_source_id, second_source_id);
    assert_ne!(first.id(), second.id());

    let empty_adjoint = sourced_certified(1.0, b"same-producer", Some(b""));
    let absent_adjoint = sourced_certified(1.0, b"same-producer", None);
    let absent_manual = manual_sourced_certified_f64_receipt(
        absent_adjoint.certified_evidence(),
        absent_adjoint.producer_source(),
        None,
    );
    let same_child_present_empty_manual = manual_sourced_certified_f64_receipt(
        absent_adjoint.certified_evidence(),
        absent_adjoint.producer_source(),
        empty_adjoint.adjoint_source(),
    );
    assert_eq!(absent_adjoint.id(), absent_manual.id());
    assert_ne!(absent_manual.id(), same_child_present_empty_manual.id());
    assert_ne!(empty_adjoint.id(), absent_adjoint.id());
    assert!(
        empty_adjoint
            .adjoint_source()
            .expect("present empty artifact")
            .canonical_bytes()
            .is_empty()
    );

    let (child, producer, adjoint) = parent.into_parts();
    assert_eq!(child.certified().value, 3.0);
    assert_eq!(producer.canonical_bytes(), b"producer-v1");
    assert_eq!(
        adjoint.expect("recovered adjoint source").canonical_bytes(),
        b"adjoint-v1"
    );
}

#[test]
fn sourced_certified_f64_identity_refuses_shape_and_correlation_mismatch() {
    let producer_bytes = b"producer";
    let adjoint_bytes = b"adjoint";
    let child_with_adjoint = || {
        identified_certified(
            Evidence::exact(1.0, ProvenanceHash::of_bytes(producer_bytes))
                .with_adjoint(ProvenanceHash::of_bytes(adjoint_bytes))
                .certified()
                .expect("valid mismatch fixture"),
        )
    };
    let producer = || certified_source("org.frankensim.tests.producer.v1", 1, producer_bytes);
    let adjoint = || certified_source("org.frankensim.tests.adjoint.v1", 1, adjoint_bytes);

    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child_with_adjoint(),
            producer(),
            None,
            LIMITS,
            || false,
        ),
        Err(
            SourcedCertifiedF64EvidenceIdentityError::AdjointSourcePresenceMismatch {
                declared: true,
                supplied: false,
            }
        )
    ));

    let child_without_adjoint = identified_certified(
        Evidence::exact(1.0, ProvenanceHash::of_bytes(producer_bytes))
            .certified()
            .expect("valid no-adjoint fixture"),
    );
    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child_without_adjoint,
            producer(),
            Some(adjoint()),
            LIMITS,
            || false,
        ),
        Err(
            SourcedCertifiedF64EvidenceIdentityError::AdjointSourcePresenceMismatch {
                declared: false,
                supplied: true,
            }
        )
    ));

    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child_with_adjoint(),
            certified_source("org.frankensim.tests.producer.v1", 1, b"wrong"),
            Some(adjoint()),
            LIMITS,
            || false,
        ),
        Err(SourcedCertifiedF64EvidenceIdentityError::SourceCorrelationMismatch {
            source: "producer",
            declared,
            computed,
        }) if declared != computed
    ));

    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child_with_adjoint(),
            producer(),
            Some(certified_source(
                "org.frankensim.tests.adjoint.v1",
                1,
                b"wrong",
            )),
            LIMITS,
            || false,
        ),
        Err(SourcedCertifiedF64EvidenceIdentityError::SourceCorrelationMismatch {
            source: "adjoint",
            declared,
            computed,
        }) if declared != computed
    ));
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one G4 refusal matrix shares exact source and parent budgets"
)]
fn sourced_certified_f64_identity_enforces_resources_and_cancellation() {
    let large_limits = CanonicalLimits::new(2 << 20, 2 << 20, 32, 64, 256);
    let exact_hard_limit = vec![0x5a; MAX_CERTIFIED_F64_SOURCE_FIELD_BYTES_V1 as usize];
    identify_certified_f64_source_v1(
        "org.frankensim.tests.large.v1".to_string(),
        1,
        exact_hard_limit,
        large_limits,
        || false,
    )
    .expect("exact hard source-byte limit is admitted");
    let over_hard_limit = vec![0x5a; MAX_CERTIFIED_F64_SOURCE_FIELD_BYTES_V1 as usize + 1];
    assert!(matches!(
        identify_certified_f64_source_v1(
            "org.frankensim.tests.large.v1".to_string(),
            1,
            over_hard_limit,
            large_limits,
            || false,
        ),
        Err(CertifiedF64SourceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested,
                limit,
            }
        )) if requested == MAX_CERTIFIED_F64_SOURCE_FIELD_BYTES_V1 + 1
            && limit == MAX_CERTIFIED_F64_SOURCE_FIELD_BYTES_V1
    ));

    let field_limits = CanonicalLimits::new(1_024, 64, 32, 64, 8);
    identify_certified_f64_source_v1("d".to_string(), 1, vec![0x11; 64], field_limits, || false)
        .expect("exact caller field limit is admitted");
    assert!(matches!(
        identify_certified_f64_source_v1("d".to_string(), 1, vec![0x11; 65], field_limits, || {
            false
        },),
        Err(CertifiedF64SourceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 65,
                limit: 64,
            }
        ))
    ));
    assert!(matches!(
        identify_certified_f64_source_v1(
            "d".to_string(),
            1,
            vec![0x11; 8],
            CanonicalLimits::new(1_024, 64, 32, 64, 0),
            || false,
        ),
        Err(CertifiedF64SourceIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    ));
    assert!(matches!(
        identify_certified_f64_source_v1("d".to_string(), 1, vec![0x11; 8], field_limits, || true,),
        Err(CertifiedF64SourceIdentityError::Canonical(
            CanonicalError::Cancelled { .. }
        ))
    ));

    let make_parent_inputs = |bytes: &[u8]| {
        let child = identified_certified(
            Evidence::exact(1.0, ProvenanceHash::of_bytes(bytes))
                .certified()
                .expect("valid resource fixture"),
        );
        let source = certified_source("org.frankensim.tests.producer.v1", 1, bytes);
        (child, source)
    };
    let (child, source) = make_parent_inputs(b"frame");
    let baseline =
        identify_sourced_certified_f64_evidence_v1(child, source, None, LIMITS, || false)
            .expect("baseline parent frame");
    let frame_limit = baseline
        .receipt()
        .canonical_bytes()
        .checked_sub(1)
        .expect("non-empty parent frame");
    let (child, source) = make_parent_inputs(b"frame");
    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child,
            source,
            None,
            CanonicalLimits::new(frame_limit, 8_192, 32, 64, 256),
            || false,
        ),
        Err(SourcedCertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CanonicalBytes,
                requested,
                limit,
            }
        )) if requested > limit && limit == frame_limit
    ));

    let (child, source) = make_parent_inputs(b"cancel-entry");
    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(child, source, None, LIMITS, || true,),
        Err(SourcedCertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { .. }
        ))
    ));
    let (child, source) = make_parent_inputs(b"zero-stride");
    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child,
            source,
            None,
            CanonicalLimits::new(16_384, 8_192, 32, 64, 0),
            || false,
        ),
        Err(SourcedCertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive")
        ))
    ));

    let (child, source) = make_parent_inputs(b"no-adjoint-collection");
    identify_sourced_certified_f64_evidence_v1(
        child,
        source,
        None,
        CanonicalLimits::new(16_384, 8_192, 32, 0, 256),
        || false,
    )
    .expect("zero ordered adjoint children fit a zero collection budget");
    let producer_bytes = b"producer-collection";
    let adjoint_bytes = b"adjoint-collection";
    let make_adjoint_inputs = || {
        let child = identified_certified(
            Evidence::exact(1.0, ProvenanceHash::of_bytes(producer_bytes))
                .with_adjoint(ProvenanceHash::of_bytes(adjoint_bytes))
                .certified()
                .expect("valid adjoint collection fixture"),
        );
        let producer = certified_source("org.frankensim.tests.producer.v1", 1, producer_bytes);
        let adjoint = certified_source("org.frankensim.tests.adjoint.v1", 1, adjoint_bytes);
        (child, producer, adjoint)
    };
    let (child, producer, adjoint) = make_adjoint_inputs();
    identify_sourced_certified_f64_evidence_v1(
        child,
        producer,
        Some(adjoint),
        CanonicalLimits::new(16_384, 8_192, 32, 1, 256),
        || false,
    )
    .expect("one ordered adjoint child fits an exact collection budget");
    let (child, producer, adjoint) = make_adjoint_inputs();
    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child,
            producer,
            Some(adjoint),
            CanonicalLimits::new(16_384, 8_192, 32, 0, 256),
            || false,
        ),
        Err(SourcedCertifiedF64EvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::CollectionItems,
                requested: 1,
                limit: 0,
            }
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

    let crosswalk_bytes = vec![0x33; 600];
    let (child, source) = make_parent_inputs(&crosswalk_bytes);
    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child,
            source,
            None,
            LIMITS,
            CancelAfter {
                successful_polls: 2,
            },
        ),
        Err(
            SourcedCertifiedF64EvidenceIdentityError::SourceCrosswalkCancelled {
                source: "producer",
                processed_bytes: 256,
            }
        )
    ));

    let adjoint_crosswalk_bytes = vec![0x44; 600];
    let producer_bytes = b"adjoint-cancel-producer";
    let child = identified_certified(
        Evidence::exact(1.0, ProvenanceHash::of_bytes(producer_bytes))
            .with_adjoint(ProvenanceHash::of_bytes(&adjoint_crosswalk_bytes))
            .certified()
            .expect("valid adjoint cancellation fixture"),
    );
    let producer = certified_source("org.frankensim.tests.producer.v1", 1, producer_bytes);
    let adjoint = certified_source(
        "org.frankensim.tests.adjoint.v1",
        1,
        &adjoint_crosswalk_bytes,
    );
    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child,
            producer,
            Some(adjoint),
            LIMITS,
            CancelAfter {
                successful_polls: 4,
            },
        ),
        Err(
            SourcedCertifiedF64EvidenceIdentityError::SourceCrosswalkCancelled {
                source: "adjoint",
                processed_bytes: 256,
            }
        )
    ));

    let poll_count = std::cell::Cell::new(0_usize);
    let (child, source) = make_parent_inputs(b"late-frame");
    identify_sourced_certified_f64_evidence_v1(child, source, None, LIMITS, || {
        poll_count.set(poll_count.get() + 1);
        false
    })
    .expect("baseline sourced-parent poll count");
    let (child, source) = make_parent_inputs(b"late-frame");
    assert!(matches!(
        identify_sourced_certified_f64_evidence_v1(
            child,
            source,
            None,
            LIMITS,
            CancelAfter {
                successful_polls: poll_count.get() - 1,
            },
        ),
        Err(SourcedCertifiedF64EvidenceIdentityError::Canonical(
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
