//! I10.1 G0/G3 conformance for the authority-separated multi-case schema.
//!
//! Tests use deterministic JSON diagnostics so the central batch verifier can
//! retain exact refusal/identity context.  No test treats a hash as laboratory
//! authentication or an identifiability theorem.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_evidence::vv::*;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, MATDB_SCHEMA_VERSION,
    MaterialCard, MaterialStateId, Provenance,
};
use fs_material::identifiability::*;
use fs_qty::{Dims, QuantitySpec};

const STRESS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const DIMENSIONLESS: Dims = Dims([0; 6]);
const TEST_HASH_DOMAIN: &str = "org.frankensim.fs-material.identifiability-authority-test.v1";

fn escape_json(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            character if character.is_control() => {
                write!(&mut escaped, "\\u{:04x}", character as u32)
                    .expect("writing JSON escape to String cannot fail");
            }
            character => escaped.push(character),
        }
    }
    escaped
}

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-material/identifiability-authority\",\
         \"case\":\"{}\",\"verdict\":\"{}\",\"detail\":\"{}\"}}",
        escape_json(case),
        escape_json(verdict),
        escape_json(detail),
    );
}

fn hash(label: &str) -> ContentHash {
    hash_domain(TEST_HASH_DOMAIN, label.as_bytes())
}

fn case_physics_hash(domain: &str, label: &str) -> ContentHash {
    hash_domain(domain, label.as_bytes())
}

fn case_physics_source(value: &str, kind: SourceKind, domain: &str) -> SourceRef {
    SourceRef::try_new(
        source_key(value),
        kind,
        case_physics_hash(domain, value),
        domain,
        CASE_PHYSICS_SOURCE_CONTRACT_VERSION,
    )
    .expect("case-physics source fixture")
}

fn artifact(value: &str) -> ArtifactId {
    ArtifactId::try_new(value).expect("fixture artifact token")
}

fn qoi(value: &str) -> QoiId {
    QoiId::try_new(value).expect("fixture QoI token")
}

fn unit(value: &str) -> UnitId {
    UnitId::try_new(value).expect("fixture unit token")
}

fn axis(value: &str) -> AxisId {
    AxisId::try_new(value).expect("fixture axis token")
}

fn role(value: &str) -> ParameterRoleId {
    ParameterRoleId::try_new(value).expect("fixture parameter role")
}

fn case_id(value: &str) -> CaseId {
    CaseId::try_new(value).expect("fixture case id")
}

fn source_key(value: &str) -> SourceKey {
    SourceKey::try_new(value).expect("fixture source key")
}

fn external_trust(label: &str, subject: &SourceRef) -> AuthorityDisposition {
    AuthorityDisposition::ExternalTrustReceipt {
        trust_receipt: TrustReceiptRef::try_new(
            source(
                "fixture-trust-receipt",
                SourceKind::EvidenceReceipt,
                hash(label),
            ),
            subject.clone(),
            TrustAuthentication::Unauthenticated,
        )
        .expect("typed external trust receipt fixture"),
    }
}

fn channel(value: &str) -> ObservationChannelId {
    ObservationChannelId::try_new(value).expect("fixture channel")
}

fn header(id: &str, capability: &str) -> ArtifactHeader {
    header_with_units(id, capability, &["Pa"])
}

fn header_with_units(id: &str, capability: &str, units: &[&str]) -> ArtifactHeader {
    ArtifactHeader::try_new(
        artifact(id),
        units.iter().map(|value| unit(value)).collect(),
        SeedDeclaration::Fixed(0x1d3_171f),
        DeclaredBudget::Limit(1.0e-9),
        DeclaredBudget::Limit(30_000),
        DeclaredBudget::Limit(32 << 20),
        vec![("fixture".to_string(), "1".to_string())],
        vec![capability.to_string()],
    )
    .expect("Five Explicits fixture")
}

fn read_u32_le(bytes: &[u8], at: &mut usize, field: &str) -> u32 {
    let end = *at + 4;
    let encoded: [u8; 4] = bytes
        .get(*at..end)
        .unwrap_or_else(|| panic!("missing {field} at byte {}", *at))
        .try_into()
        .expect("four-byte canonical u32");
    *at = end;
    u32::from_le_bytes(encoded)
}

fn read_u64_le(bytes: &[u8], at: &mut usize, field: &str) -> u64 {
    let end = *at + 8;
    let encoded: [u8; 8] = bytes
        .get(*at..end)
        .unwrap_or_else(|| panic!("missing {field} at byte {}", *at))
        .try_into()
        .expect("eight-byte canonical u64");
    *at = end;
    u64::from_le_bytes(encoded)
}

fn read_text<'a>(bytes: &'a [u8], at: &mut usize, field: &str) -> &'a str {
    let len = usize::try_from(read_u32_le(bytes, at, field)).expect("u32 fits usize");
    let end = *at + len;
    let value = std::str::from_utf8(
        bytes
            .get(*at..end)
            .unwrap_or_else(|| panic!("truncated {field} at byte {}", *at)),
    )
    .unwrap_or_else(|error| panic!("non-UTF-8 {field}: {error}"));
    *at = end;
    value
}

/// Assert the exact identity-mode header grammar and return the first byte
/// after it. This is intentionally independent of the production decoder: it
/// pins the projection marker, field order, collection framing, and numeric
/// endianness that the identity declarations promise.
fn assert_identity_header_layout(bytes: &[u8], magic: &[u8], header: &ArtifactHeader) -> usize {
    assert!(bytes.starts_with(magic), "identity wire magic moved");
    let mut at = magic.len();
    assert_eq!(
        read_u32_le(bytes, &mut at, "schema version"),
        IDENTIFIABILITY_AUTHORITY_SCHEMA_VERSION,
    );
    assert_eq!(
        bytes.get(at),
        Some(&0),
        "identity header marker must be zero"
    );
    at += 1;

    assert_eq!(
        read_u32_le(bytes, &mut at, "header unit count"),
        u32::try_from(header.units().len()).expect("bounded unit count"),
    );
    for expected in header.units() {
        assert_eq!(read_text(bytes, &mut at, "header unit"), expected.as_str());
    }

    match header.seed() {
        SeedDeclaration::Fixed(seed) => {
            assert_eq!(bytes.get(at), Some(&0), "fixed-seed tag moved");
            at += 1;
            assert_eq!(read_u64_le(bytes, &mut at, "seed"), *seed);
        }
        SeedDeclaration::NotApplicable { .. } => {
            panic!("wire-layout fixture unexpectedly uses a non-numeric seed")
        }
    }
    match header.accuracy() {
        DeclaredBudget::Limit(value) => {
            assert_eq!(bytes.get(at), Some(&0), "accuracy-limit tag moved");
            at += 1;
            assert_eq!(read_u64_le(bytes, &mut at, "accuracy"), value.to_bits());
        }
        DeclaredBudget::NotApplicable { .. } => {
            panic!("wire-layout fixture unexpectedly omits accuracy")
        }
    }
    for (field, budget) in [
        ("time budget", header.time_ms()),
        ("memory budget", header.memory_bytes()),
    ] {
        match budget {
            DeclaredBudget::Limit(value) => {
                assert_eq!(bytes.get(at), Some(&0), "{field} limit tag moved");
                at += 1;
                assert_eq!(read_u64_le(bytes, &mut at, field), *value);
            }
            DeclaredBudget::NotApplicable { .. } => {
                panic!("wire-layout fixture unexpectedly omits {field}")
            }
        }
    }

    assert_eq!(
        read_u32_le(bytes, &mut at, "header version count"),
        u32::try_from(header.versions().len()).expect("bounded version count"),
    );
    for (component, version) in header.versions() {
        assert_eq!(read_text(bytes, &mut at, "version component"), component);
        assert_eq!(read_text(bytes, &mut at, "version value"), version);
    }
    assert_eq!(
        read_u32_le(bytes, &mut at, "header capability count"),
        u32::try_from(header.capabilities().len()).expect("bounded capability count"),
    );
    for capability in header.capabilities() {
        assert_eq!(read_text(bytes, &mut at, "capability"), capability);
    }
    at
}

fn project_exact_header_to_identity(
    exact: &[u8],
    magic: &[u8],
    header: &ArtifactHeader,
) -> Vec<u8> {
    assert!(exact.starts_with(magic), "exact transport magic moved");
    let marker_at = magic.len() + 4;
    assert_eq!(
        exact.get(marker_at),
        Some(&1),
        "exact header marker must be one"
    );
    let mut id_len_at = marker_at + 1;
    let id_len = usize::try_from(read_u32_le(exact, &mut id_len_at, "artifact id length"))
        .expect("artifact id length fits usize");
    assert_eq!(id_len, header.id().as_str().len());
    let id_end = id_len_at + id_len;
    assert_eq!(
        exact.get(id_len_at..id_end),
        Some(header.id().as_str().as_bytes()),
        "exact artifact label moved",
    );

    let mut projected = Vec::with_capacity(exact.len() - 4 - id_len);
    projected.extend_from_slice(&exact[..marker_at]);
    projected.push(0);
    projected.extend_from_slice(&exact[id_end..]);
    projected
}

fn context() -> ContextOfUse {
    ContextOfUse::try_new(
        header_with_units("context-1", "fixture.context", &["Pa", "K"]),
        "Choose a constitutive calibration that predicts stress response.",
        vec![
            QoiSpec::try_new(
                qoi("stress"),
                "axial stress",
                unit("Pa"),
                AcceptanceCriterion::ClosedRange {
                    lo: -2.0e9,
                    hi: 2.0e9,
                },
            )
            .expect("stress QoI"),
            QoiSpec::try_new(
                qoi("tangent"),
                "algorithmic tangent",
                unit("Pa"),
                AcceptanceCriterion::ClosedRange {
                    lo: -2.0e9,
                    hi: 2.0e9,
                },
            )
            .expect("tangent QoI"),
        ],
        ApplicabilityDomain::try_new(
            vec![
                NumericDomainAxis::try_new(axis("temperature"), unit("K"), 250.0, 450.0)
                    .expect("temperature applicability"),
            ],
            Vec::new(),
        )
        .expect("applicability domain"),
        ApplicabilityPolicy::Demote,
    )
    .expect("context fixture")
}

fn model_cards() -> (MaterialCard, ConstitutiveModelCard) {
    let mut parameters = BTreeMap::new();
    parameters.insert(
        "yield_stress".to_string(),
        LawParameter {
            value: 276.0e6,
            dims: STRESS,
        },
    );
    parameters.insert(
        "hardening_modulus".to_string(),
        LawParameter {
            value: 1.2e9,
            dims: STRESS,
        },
    );
    let model = ConstitutiveModelCard {
        law: LawId("j2-identifiability-authority-fixture".to_string()),
        law_version: 3,
        parameters,
        state_schema_version: 2,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("temperature", 250.0, 450.0),
        sources: vec![hash("model-source")],
        provenance: Provenance {
            source: "authority fixture".to_string(),
            license: "test-only".to_string(),
            artifact: Some(hash("model-provenance")),
        },
    };
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "AA6061".to_string(),
            phase: "wrought".to_string(),
            process: "T6".to_string(),
            revision: 0,
        },
        ClaimSet::new(),
        vec![model.clone()],
    )
    .expect("material card fixture");
    (material, model)
}

fn source(value: &str, kind: SourceKind, content: ContentHash) -> SourceRef {
    let (domain, version) = match kind {
        SourceKind::ContextOfUse
        | SourceKind::ExperimentArtifact
        | SourceKind::CalibrationSplit => (VV_ARTIFACT_SOURCE_DOMAIN, VV_SCHEMA_VERSION),
        SourceKind::MaterialCard => (MATERIAL_CARD_SOURCE_DOMAIN, MATDB_SCHEMA_VERSION),
        SourceKind::ConstitutiveModelCard => {
            (CONSTITUTIVE_MODEL_CARD_SOURCE_DOMAIN, MATDB_SCHEMA_VERSION)
        }
        SourceKind::Parser => (TEST_HASH_DOMAIN, 2),
        SourceKind::ObservationOperator => (TEST_HASH_DOMAIN, 4),
        _ => (TEST_HASH_DOMAIN, 1),
    };
    SourceRef::try_new(source_key(value), kind, content, domain, version)
        .expect("source reference fixture")
}

fn frame(case: &str) -> FrameBinding {
    let label = format!("frame-transform-{case}");
    FrameBinding::try_new(
        artifact(&format!("frame-{case}")),
        case_physics_hash(FRAME_TRANSFORM_SOURCE_DOMAIN, &label),
        "right-handed-cartesian",
    )
    .expect("frame fixture")
}

fn protocol(case: &str) -> ProtocolBinding {
    let load = format!("load-path-{case}");
    let environment = format!("environment-path-{case}");
    let time = format!("time-grid-{case}");
    ProtocolBinding::try_new(
        artifact(&format!("protocol-{case}")),
        7,
        2,
        3,
        case_physics_hash(LOAD_PATH_SOURCE_DOMAIN, &load),
        case_physics_hash(ENVIRONMENT_PATH_SOURCE_DOMAIN, &environment),
        case_physics_hash(TIME_GRID_SOURCE_DOMAIN, &time),
        artifact(&format!("clock-{case}")),
    )
    .expect("protocol fixture")
}

fn specimen(case: &str, frame: FrameBinding) -> SpecimenBinding {
    let geometry = format!("geometry-{case}");
    let process = format!("process-{case}");
    let preparation = format!("preparation-{case}");
    SpecimenBinding::try_new(
        artifact(&format!("specimen-{case}")),
        case_physics_hash(SPECIMEN_GEOMETRY_SOURCE_DOMAIN, &geometry),
        case_physics_hash(SPECIMEN_PROCESS_SOURCE_DOMAIN, &process),
        case_physics_hash(SPECIMEN_PREPARATION_SOURCE_DOMAIN, &preparation),
        frame,
    )
    .expect("specimen fixture")
}

fn parameter(
    name: &str,
    treatment: ParameterTreatment,
    coverage: InfluenceCoverage,
    prior_version: u32,
    domain_override: Option<ParameterDomain>,
    prior_absent: bool,
    scope: ParameterScopeBinding,
) -> StudyParameter {
    let domain = domain_override.unwrap_or_else(|| {
        if name == "yield_stress" {
            ParameterDomain::try_new(1.0e6, 1.0e9).expect("yield domain")
        } else {
            ParameterDomain::try_new(1.0e7, 5.0e9).expect("hardening domain")
        }
    });
    StudyParameter::try_new(
        role(name),
        QuantitySpec::dimensional(STRESS),
        domain,
        if name == "yield_stress" {
            ParameterPurpose::Estimand
        } else {
            ParameterPurpose::Nuisance
        },
        treatment,
        ParameterOwnerBinding::ConstitutiveModel,
        scope,
        if prior_absent {
            PriorPolicy::Absent {
                reason: "no probability measure has been declared for this estimand".to_string(),
            }
        } else {
            PriorPolicy::Distribution(ParameterPrior::Uniform {
                version: prior_version,
                domain,
            })
        },
        coverage,
    )
    .expect("study parameter fixture")
}

#[derive(Debug, Clone, Copy, Default)]
struct ProblemOptions {
    reverse_cases: bool,
    missing_hardening_influence: bool,
    dangling_operator: bool,
    dense_with_bounded_marginal: bool,
    overlapping_gauges: bool,
    bad_constraint_units: bool,
    derived_cycle: bool,
    retrospective_reuse: bool,
    declared_sharing: bool,
    bad_observation_endpoint: bool,
    self_correlation: bool,
    alternate_graph_domain: bool,
    context_contract_mutation: u8,
    observation_contract_mutation: u8,
    blind_prospective_case: bool,
    second_case_complementary: bool,
    claim_strata_in_problem: bool,
    parameter_prior_version: u32,
    valid_constraint: bool,
    ordered_constraint_case: u8,
    yield_log_scale: bool,
    one_gauge: bool,
    external_noise: bool,
    alternate_sharing_justification: bool,
    yield_prior_absent: bool,
    yield_case_a_only: bool,
    case_physics_mutation: u8,
    modeled_discrepancy_case: u8,
    composite_influence_chain: bool,
    gauge_case: u8,
    yield_influence_case_b: bool,
    joint_prior_choice: u8,
    declared_gauge_composition: bool,
    independent_gauge_composition: bool,
}

struct ProblemFixture {
    context: ContextOfUse,
    material: MaterialCard,
    model: ConstitutiveModelCard,
    graph: ContentHash,
    document: Result<IdentifiabilityProblemDocument, IdentifiabilityError>,
}

#[derive(Debug, Clone, Copy)]
enum ProblemRoot {
    Context,
    Material,
    Model,
    Graph,
}

fn make_case(
    name: &str,
    purpose: CasePurpose,
    qoi_name: &str,
    channel_name: &str,
    bounded_noise: bool,
    retrospective: bool,
    experiment_key: &str,
    observation_contract_mutation: u8,
) -> StudyCaseDocument {
    let frame = frame(name);
    let protocol = protocol(name);
    let (observation_clock, observation_protocol_version, observation_refinement_version) =
        match observation_contract_mutation {
            0 => (protocol.clock().clone(), 7, 3),
            1 => (protocol.clock().clone(), 8, 3),
            2 => (protocol.clock().clone(), 7, 4),
            3 => (artifact(&format!("wrong-clock-{name}")), 7, 3),
            other => panic!("unsupported observation-contract mutation {other}"),
        };
    let rows = if retrospective {
        ObservationRows::Retrospective(BTreeSet::from([ObservationId::try_new(format!(
            "row-{name}"
        ))
        .expect("row fixture")]))
    } else {
        ObservationRows::Prospective
    };
    let observation = StudyObservation::try_new(
        channel(channel_name),
        qoi(qoi_name),
        unit("Pa"),
        QuantitySpec::dimensional(STRESS),
        source_key("unit-pa"),
        frame.clone(),
        format!("node-{name}"),
        "stress-output",
        source_key(if name == "a" {
            "operator-a"
        } else {
            "operator-b"
        }),
        source_key(if name == "a" {
            "aggregation-a"
        } else {
            "aggregation-b"
        }),
        source_key(if name == "a" { "sensor-a" } else { "sensor-b" }),
        artifact(&format!("instrument-{name}")),
        artifact(&format!("acquisition-{name}")),
        observation_clock,
        4,
        if bounded_noise {
            MarginalNoiseSpec::Bounded { half_width: 1.0 }
        } else {
            MarginalNoiseSpec::Gaussian {
                standard_deviation: 2.0e5,
            }
        },
        MissingnessAssumption::Unknown {
            reason: "missingness has not yet been characterized".to_string(),
        },
        None,
        observation_protocol_version,
        observation_refinement_version,
        rows,
    )
    .expect("observation fixture");
    let data = if retrospective {
        CaseDataDeclaration::Retrospective {
            experiment: source_key(experiment_key),
            split: source_key(if name == "a" { "split-a" } else { "split-b" }),
            parser: source_key("parser"),
            preprocessing: source_key("preprocessing"),
            parser_version: 2,
            split_grouping: artifact("split-by-specimen"),
        }
    } else {
        CaseDataDeclaration::Prospective
    };
    StudyCaseDocument::try_new(
        case_id(name),
        purpose,
        InitialStateBinding::Zero { schema_version: 2 },
        specimen(name, frame),
        protocol,
        CasePhysicsSources::new(
            source_key(&format!("frame-transform-{name}")),
            source_key(&format!("geometry-{name}")),
            source_key(&format!("process-{name}")),
            source_key(&format!("preparation-{name}")),
            source_key(&format!("load-path-{name}")),
            source_key(&format!("environment-path-{name}")),
            source_key(&format!("time-grid-{name}")),
            None,
        ),
        source_key(if name == "a" {
            "forward-a"
        } else {
            "forward-b"
        }),
        data,
        vec![observation],
        vec![(
            channel(channel_name),
            StudyDiscrepancy::Uncharacterized {
                reason: "no discrepancy family is admitted for this channel".to_string(),
            },
        )],
        Vec::new(),
    )
    .expect("case fixture")
}

fn retrospective_artifacts(
    origin: ExperimentOrigin,
) -> (ExperimentArtifact, CalibrationSplit, CalibrationSplit) {
    let source_bytes_hash = hash("retrospective-origin-source-bytes");
    let instrument_a = artifact("instrument-a");
    let instrument_b = artifact("instrument-b");
    let clock_a = artifact("clock-a");
    let clock_b = artifact("clock-b");
    let rows = [
        ("row-a", "stress", instrument_a.clone(), clock_a.clone()),
        ("row-b", "tangent", instrument_b.clone(), clock_b.clone()),
        (
            "row-validation",
            "stress",
            instrument_a.clone(),
            clock_a.clone(),
        ),
        (
            "row-blind",
            "tangent",
            instrument_b.clone(),
            clock_b.clone(),
        ),
    ];
    let manifest = ObservationManifest::try_new(
        rows.iter()
            .map(|(row, qoi_id, instrument, clock)| {
                (
                    ObservationId::try_new(*row).expect("retrospective row id"),
                    ObservationManifestRow::try_new(
                        ObservationSourceRef::try_new(
                            source_bytes_hash,
                            TEST_HASH_DOMAIN,
                            1,
                            hash(&format!("locator-{row}")),
                            hash(&format!("extraction-{row}")),
                        )
                        .expect("retrospective row source"),
                        qoi(qoi_id),
                        instrument.clone(),
                        artifact(match *row {
                            "row-a" => "acquisition-a",
                            "row-b" => "acquisition-b",
                            "row-validation" => "acquisition-validation",
                            "row-blind" => "acquisition-blind",
                            _ => unreachable!("bounded retrospective row fixture"),
                        }),
                        clock.clone(),
                    )
                    .expect("retrospective manifest row"),
                )
            })
            .collect(),
    )
    .expect("retrospective manifest");
    let experiment = ExperimentArtifact::try_new(
        header("experiment-shared", "fixture.experiment"),
        artifact("dataset-shared"),
        origin,
        vec![qoi("stress"), qoi("tangent")],
        manifest,
        vec![
            InstrumentCalibration::new(instrument_a, hash("calibration-a"), true),
            InstrumentCalibration::new(instrument_b, hash("calibration-b"), true),
        ],
        ClockSynchronization::synchronized(
            vec![clock_a, clock_b],
            "fixture-synchronized-clocks",
            1.0e-9,
            hash("clock-synchronization"),
        )
        .expect("retrospective clock topology"),
        RepeatabilitySummary::try_new(
            3,
            CovarianceMatrix::try_new(2, vec![1.0, 0.0, 1.0]).expect("retrospective covariance"),
        )
        .expect("retrospective repeatability"),
        DataAuthenticity::new(source_bytes_hash, hash("retrospective-custody"), true),
    )
    .expect("retrospective experiment");
    let experiment_hash = experiment.content_hash().expect("experiment content hash");
    let split = |id: &str| {
        CalibrationSplit::try_new(
            header(id, "fixture.split"),
            ArtifactRef::new(
                ArtifactKind::ExperimentArtifact,
                experiment.id().clone(),
                experiment_hash,
            ),
            hash(&format!("preregistration-{id}")),
            vec![
                ObservationId::try_new("row-a").expect("row-a"),
                ObservationId::try_new("row-b").expect("row-b"),
            ],
            vec![ObservationId::try_new("row-validation").expect("validation row")],
            vec![(
                ObservationId::try_new("row-blind").expect("blind row"),
                hash("locator-row-blind"),
            )],
        )
        .expect("retrospective split")
    };
    (experiment, split("split-a"), split("split-b"))
}

#[derive(Debug, Clone, Copy)]
enum DiscrepancyOriginFixture {
    Physical,
    DeclaredSynthetic {
        declared_producer: &'static str,
        stale_forward_binding: bool,
        production_binding_key: &'static str,
    },
    Uncharacterized,
}

struct RetrospectiveOriginFixture {
    problem: ProblemFixture,
    experiment: ExperimentArtifact,
    split_a: CalibrationSplit,
    split_b: CalibrationSplit,
}

fn retrospective_origin_fixture(
    origin: ExperimentOrigin,
    purpose: CasePurpose,
    discrepancy: DiscrepancyOriginFixture,
) -> RetrospectiveOriginFixture {
    retrospective_origin_fixture_with_options(
        origin,
        purpose,
        discrepancy,
        ProblemOptions::default(),
    )
}

fn retrospective_origin_fixture_with_options(
    origin: ExperimentOrigin,
    purpose: CasePurpose,
    discrepancy: DiscrepancyOriginFixture,
    options: ProblemOptions,
) -> RetrospectiveOriginFixture {
    let (experiment, split_a, split_b) = retrospective_artifacts(origin);
    let base = problem_fixture(ProblemOptions {
        retrospective_reuse: true,
        declared_sharing: true,
        ..options
    });
    let document = base.document.expect("retrospective structural document");
    let mut sources = document.sources().values().cloned().collect::<Vec<_>>();
    for source in &mut sources {
        *source = match source.key().as_str() {
            "experiment-shared" => SourceRef::experiment(source.key().clone(), &experiment)
                .expect("typed experiment source"),
            "split-a" => SourceRef::calibration_split(source.key().clone(), &split_a)
                .expect("typed split-a source"),
            "split-b" => SourceRef::calibration_split(source.key().clone(), &split_b)
                .expect("typed split-b source"),
            _ => source.clone(),
        };
    }
    let assumption_key = source_key("discrepancy-origin-assumption");
    let basis = match discrepancy {
        DiscrepancyOriginFixture::Physical => {
            sources.push(source(
                assumption_key.as_str(),
                SourceKind::Assumption,
                hash(assumption_key.as_str()),
            ));
            Some(DiscrepancyInapplicability::PhysicalApplicability {
                assumption: assumption_key.clone(),
            })
        }
        DiscrepancyOriginFixture::DeclaredSynthetic {
            declared_producer,
            stale_forward_binding,
            production_binding_key,
        } => {
            let producer = artifact(declared_producer);
            let generator_key = source_key("forward-a");
            let generator_index = sources
                .iter()
                .position(|source| source.key() == &generator_key)
                .expect("forward source fixture");
            let binding_forward = sources[generator_index].clone();
            let binding_preimage =
                forward_model_production_binding_preimage(&producer, &binding_forward)
                    .expect("producer binding preimage");
            if stale_forward_binding {
                sources[generator_index] = source(
                    "forward-a",
                    SourceKind::ForwardModel,
                    hash("different-forward-a-content"),
                );
            }
            let binding_key = source_key(production_binding_key);
            sources.extend([
                source(
                    assumption_key.as_str(),
                    SourceKind::Assumption,
                    hash(assumption_key.as_str()),
                ),
                SourceRef::try_new(
                    binding_key.clone(),
                    SourceKind::ForwardModelProductionBinding,
                    hash_domain(FORWARD_MODEL_PRODUCTION_BINDING_DOMAIN, &binding_preimage),
                    FORWARD_MODEL_PRODUCTION_BINDING_DOMAIN,
                    FORWARD_MODEL_PRODUCTION_BINDING_VERSION,
                )
                .expect("producer binding source"),
            ]);
            Some(DiscrepancyInapplicability::DeclaredSyntheticSelfModel {
                generator: generator_key,
                producer,
                production_binding: binding_key,
                assumption: assumption_key.clone(),
            })
        }
        DiscrepancyOriginFixture::Uncharacterized => None,
    };
    let cases = document
        .cases()
        .values()
        .map(|case| {
            if case.id().as_str() != "a" {
                return case.clone();
            }
            StudyCaseDocument::try_new(
                case.id().clone(),
                purpose.clone(),
                case.initial_state(),
                case.specimen().clone(),
                case.protocol().clone(),
                case.physics_sources().clone(),
                case.forward_model().clone(),
                case.data().clone(),
                case.observations().values().cloned().collect(),
                case.discrepancies()
                    .keys()
                    .cloned()
                    .map(|channel| {
                        let discrepancy = basis.clone().map_or_else(
                            || StudyDiscrepancy::Uncharacterized {
                                reason: "fixture leaves discrepancy applicability open".to_string(),
                            },
                            |basis| StudyDiscrepancy::NotApplicable { basis },
                        );
                        (channel, discrepancy)
                    })
                    .collect(),
                case.observation_sharing().to_vec(),
            )
            .expect("rebuilt origin case")
        })
        .collect();
    let document = IdentifiabilityProblemDocument::try_new(
        document.context_source().clone(),
        document.material_source().clone(),
        document.model_source().clone(),
        document.graph_source().clone(),
        document.joint_prior().cloned(),
        sources,
        document.parameters().values().cloned().collect(),
        document.constraints().values().cloned().collect(),
        document.admissible_domain().clone(),
        cases,
        document.influences().values().cloned().collect(),
        document.gauges().values().cloned().collect(),
        document.gauge_compositions().values().cloned().collect(),
        document.joint_noise().clone(),
        document.data_reuse().clone(),
    );
    RetrospectiveOriginFixture {
        problem: ProblemFixture {
            context: base.context,
            material: base.material,
            model: base.model,
            graph: base.graph,
            document,
        },
        experiment,
        split_a,
        split_b,
    }
}

fn admit_retrospective_origin_fixture(
    fixture: &RetrospectiveOriginFixture,
) -> Result<AdmittedIdentifiabilityProblem, IdentifiabilityError> {
    let document = fixture
        .problem
        .document
        .clone()
        .expect("origin document admits structurally");
    let opaque = opaque_resolutions(&document);
    AdmittedIdentifiabilityProblem::resolve_and_admit(
        document,
        ProblemSourceBundle::new(
            &fixture.problem.context,
            &fixture.problem.material,
            &fixture.problem.model,
            BTreeMap::from([
                (
                    case_id("a"),
                    CaseSourceBundle::new(&fixture.experiment, &fixture.split_a),
                ),
                (
                    case_id("b"),
                    CaseSourceBundle::new(&fixture.experiment, &fixture.split_b),
                ),
            ]),
            opaque,
        ),
    )
}

fn gauge_validity_fixture(
    members: &BTreeSet<ParameterRoleId>,
    local_obstruction: bool,
) -> GaugeValidityScope {
    let local = if local_obstruction {
        members.clone()
    } else {
        BTreeSet::new()
    };
    let per_case = ["a", "b"]
        .into_iter()
        .map(|name| {
            (
                case_id(name),
                GaugeExtentSupport::try_new(local.clone(), members.clone())
                    .expect("gauge extent support fixture"),
            )
        })
        .collect();
    GaugeValidityScope::try_new(BTreeMap::from([(
        GaugeApplicabilityAxes::new(
            GaugeInformationRegime::StructuralExactModel,
            GaugeScalarDomain::Real,
            GaugeLocus::WholeDomain,
            GaugeQuantifierScope::AtRealization {
                realization: source_key("fixture-gauge-realization"),
            },
        ),
        GaugeCellDomain::try_new(per_case).expect("gauge cell fixture"),
    )]))
    .expect("gauge validity fixture")
}

fn gauge_fixture(
    id: &str,
    action: &str,
    members: BTreeSet<ParameterRoleId>,
    algebra: GaugeAlgebra,
    orbit_geometry: GaugeOrbitGeometry,
    local_obstruction: bool,
) -> GaugeDeclaration {
    let validity = gauge_validity_fixture(&members, local_obstruction);
    GaugeDeclaration::try_new(
        GaugeClassId::try_new(id).expect("gauge id"),
        members,
        source_key(action),
        algebra,
        orbit_geometry,
        GaugeStatus::Assumed {
            assumption: source_key("fixture-gauge-assumption"),
        },
        validity,
    )
    .expect("gauge fixture")
}

fn problem_fixture(options: ProblemOptions) -> ProblemFixture {
    let context = context();
    let (material, model) = model_cards();
    let graph = hash("constitutive-graph");
    let context_hash = context.content_hash().expect("context hashes");
    let context_source = match options.context_contract_mutation {
        0 => source("context", SourceKind::ContextOfUse, context_hash),
        1 => source("context", SourceKind::ContextOfUse, hash("wrong-context")),
        2 => SourceRef::try_new(
            source_key("context"),
            SourceKind::ContextOfUse,
            context_hash,
            TEST_HASH_DOMAIN,
            VV_SCHEMA_VERSION,
        )
        .expect("wrong-domain context reference"),
        3 => SourceRef::try_new(
            source_key("context"),
            SourceKind::ContextOfUse,
            context_hash,
            VV_ARTIFACT_SOURCE_DOMAIN,
            VV_SCHEMA_VERSION + 1,
        )
        .expect("wrong-version context reference"),
        other => panic!("unsupported context contract mutation {other}"),
    };
    let graph_source = if options.alternate_graph_domain {
        SourceRef::try_new(
            source_key("graph"),
            SourceKind::ConstitutiveGraph,
            graph,
            "org.frankensim.test.alternate-graph-domain.v1",
            1,
        )
        .expect("alternate graph-domain reference")
    } else {
        source("graph", SourceKind::ConstitutiveGraph, graph)
    };
    let mut sources = vec![
        context_source,
        source(
            "material",
            SourceKind::MaterialCard,
            material.content_hash(),
        ),
        source(
            "model",
            SourceKind::ConstitutiveModelCard,
            model.content_hash(),
        ),
        graph_source,
        source("forward-a", SourceKind::ForwardModel, hash("forward-a")),
        source("forward-b", SourceKind::ForwardModel, hash("forward-b")),
        source(
            "operator-a",
            SourceKind::ObservationOperator,
            hash("operator-a"),
        ),
        source(
            "operator-b",
            SourceKind::ObservationOperator,
            hash("operator-b"),
        ),
        source(
            "aggregation-a",
            SourceKind::ObservationOperator,
            hash("aggregation-a"),
        ),
        source(
            "aggregation-b",
            SourceKind::ObservationOperator,
            hash("aggregation-b"),
        ),
        source("sensor-a", SourceKind::Metrology, hash("sensor-a")),
        source("sensor-b", SourceKind::Metrology, hash("sensor-b")),
        source("unit-pa", SourceKind::UnitDefinition, hash("unit-pa")),
    ];
    match options.joint_prior_choice {
        0 => {}
        1 => sources.push(source(
            "joint-prior-measure-a",
            SourceKind::ProbabilityMeasure,
            hash("joint-prior-measure-a"),
        )),
        2 => sources.push(source(
            "joint-prior-measure-b",
            SourceKind::ProbabilityMeasure,
            hash("joint-prior-measure-b"),
        )),
        other => panic!("unsupported joint-prior fixture choice {other}"),
    }
    for case in ["a", "b"] {
        sources.extend([
            case_physics_source(
                &format!("frame-transform-{case}"),
                SourceKind::Geometry,
                FRAME_TRANSFORM_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("geometry-{case}"),
                SourceKind::Geometry,
                SPECIMEN_GEOMETRY_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("process-{case}"),
                SourceKind::Process,
                SPECIMEN_PROCESS_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("preparation-{case}"),
                SourceKind::Process,
                SPECIMEN_PREPARATION_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("load-path-{case}"),
                SourceKind::Protocol,
                LOAD_PATH_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("environment-path-{case}"),
                SourceKind::Protocol,
                ENVIRONMENT_PATH_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("time-grid-{case}"),
                SourceKind::Protocol,
                TIME_GRID_SOURCE_DOMAIN,
            ),
        ]);
    }
    match options.case_physics_mutation {
        0 => {}
        1 => sources.retain(|source| source.key().as_str() != "geometry-a"),
        mutation @ 2..=5 => {
            let target = sources
                .iter_mut()
                .find(|source| source.key().as_str() == "geometry-a")
                .expect("case-physics mutation target");
            let (kind, expected_hash, domain, version) = match mutation {
                2 => (
                    SourceKind::Process,
                    target.expected_hash(),
                    SPECIMEN_GEOMETRY_SOURCE_DOMAIN,
                    CASE_PHYSICS_SOURCE_CONTRACT_VERSION,
                ),
                3 => (
                    SourceKind::Geometry,
                    case_physics_hash(SPECIMEN_GEOMETRY_SOURCE_DOMAIN, "wrong-geometry-a"),
                    SPECIMEN_GEOMETRY_SOURCE_DOMAIN,
                    CASE_PHYSICS_SOURCE_CONTRACT_VERSION,
                ),
                4 => (
                    SourceKind::Geometry,
                    target.expected_hash(),
                    TEST_HASH_DOMAIN,
                    CASE_PHYSICS_SOURCE_CONTRACT_VERSION,
                ),
                5 => (
                    SourceKind::Geometry,
                    target.expected_hash(),
                    SPECIMEN_GEOMETRY_SOURCE_DOMAIN,
                    CASE_PHYSICS_SOURCE_CONTRACT_VERSION + 1,
                ),
                _ => unreachable!(),
            };
            *target = SourceRef::try_new(
                source_key("geometry-a"),
                kind,
                expected_hash,
                domain,
                version,
            )
            .expect("mutated case-physics source");
        }
        other => panic!("unsupported case-physics mutation {other}"),
    }
    if options.dangling_operator {
        sources.retain(|source| source.key().as_str() != "operator-b");
    }
    if options.retrospective_reuse {
        sources.extend([
            source(
                "experiment-shared",
                SourceKind::ExperimentArtifact,
                hash("experiment-shared"),
            ),
            source("split-a", SourceKind::CalibrationSplit, hash("split-a")),
            source("split-b", SourceKind::CalibrationSplit, hash("split-b")),
            source("parser", SourceKind::Parser, hash("parser")),
            source(
                "preprocessing",
                SourceKind::Preprocessing,
                hash("preprocessing"),
            ),
        ]);
        if options.declared_sharing {
            sources.push(source(
                "joint-likelihood",
                SourceKind::Likelihood,
                hash("joint-likelihood"),
            ));
        }
    }
    let mut cases = vec![
        make_case(
            "a",
            if options.blind_prospective_case {
                CasePurpose::BlindFalsification
            } else {
                CasePurpose::Calibration
            },
            "stress",
            "stress",
            options.dense_with_bounded_marginal,
            options.retrospective_reuse,
            "experiment-shared",
            options.observation_contract_mutation,
        ),
        make_case(
            "b",
            if options.second_case_complementary {
                CasePurpose::Complementary {
                    reason: "case b supplies a complementary excitation".to_string(),
                }
            } else {
                CasePurpose::SymmetryBreaking
            },
            "tangent",
            "tangent",
            false,
            options.retrospective_reuse,
            "experiment-shared",
            0,
        ),
    ];
    if options.reverse_cases {
        cases.reverse();
    }
    let ordered_domains = match options.ordered_constraint_case {
        0 => (None, None),
        1 => (
            Some(ParameterDomain::try_new(10.0, 20.0).expect("ordered left domain")),
            Some(ParameterDomain::try_new(0.0, 1.0).expect("ordered right domain")),
        ),
        2 | 4 => (
            Some(ParameterDomain::try_new(10.0, 20.0).expect("ordered left domain")),
            Some(ParameterDomain::try_new(20.0, 30.0).expect("ordered right domain")),
        ),
        3 => (
            Some(ParameterDomain::try_new(10.0, 20.0).expect("ordered left domain")),
            Some(ParameterDomain::try_new(0.0, 10.0).expect("ordered right domain")),
        ),
        other => panic!("unsupported ordered-constraint fixture {other}"),
    };
    let mut parameters = vec![
        parameter(
            "yield_stress",
            ParameterTreatment::Estimated,
            InfluenceCoverage::Declared,
            options.parameter_prior_version.max(1),
            ordered_domains.0,
            options.yield_prior_absent,
            if options.yield_case_a_only {
                ParameterScopeBinding::Cases(BTreeSet::from([case_id("a")]))
            } else {
                ParameterScopeBinding::Global
            },
        ),
        parameter(
            "hardening_modulus",
            ParameterTreatment::Marginalized,
            InfluenceCoverage::Declared,
            options.parameter_prior_version.max(1),
            ordered_domains.1,
            false,
            ParameterScopeBinding::Global,
        ),
    ];
    if options.modeled_discrepancy_case != 0 {
        sources.extend([
            source(
                "fixture-discrepancy-family",
                SourceKind::Discrepancy,
                hash("fixture-discrepancy-family"),
            ),
            source(
                "fixture-discrepancy-guard",
                SourceKind::Constraint,
                hash("fixture-discrepancy-guard"),
            ),
        ]);
        let discrepancy_domain =
            ParameterDomain::try_new(-1.0e6, 1.0e6).expect("discrepancy domain");
        parameters.push(
            StudyParameter::try_new(
                role("discrepancy_bias"),
                QuantitySpec::dimensional(STRESS),
                discrepancy_domain,
                ParameterPurpose::Nuisance,
                ParameterTreatment::Profiled,
                ParameterOwnerBinding::Discrepancy {
                    family: source_key("fixture-discrepancy-family"),
                },
                ParameterScopeBinding::Cases(BTreeSet::from([case_id("a")])),
                PriorPolicy::Distribution(ParameterPrior::Uniform {
                    version: 1,
                    domain: discrepancy_domain,
                }),
                InfluenceCoverage::IntentionallyAbsent {
                    reason: "this fixture isolates discrepancy ownership and case scope"
                        .to_string(),
                },
            )
            .expect("modeled discrepancy parameter"),
        );

        let target_case = match options.modeled_discrepancy_case {
            1 => "a",
            2 => "b",
            other => panic!("unsupported modeled-discrepancy fixture {other}"),
        };
        let position = cases
            .iter()
            .position(|case| case.id() == &case_id(target_case))
            .expect("modeled-discrepancy target case");
        let case = cases.remove(position);
        let target_channel = if target_case == "a" {
            channel("stress")
        } else {
            channel("tangent")
        };
        let rebuilt = StudyCaseDocument::try_new(
            case.id().clone(),
            case.purpose().clone(),
            case.initial_state(),
            case.specimen().clone(),
            case.protocol().clone(),
            case.physics_sources().clone(),
            case.forward_model().clone(),
            case.data().clone(),
            case.observations().values().cloned().collect(),
            vec![(
                target_channel,
                StudyDiscrepancy::Modeled {
                    family: source_key("fixture-discrepancy-family"),
                    parameters: BTreeSet::from([role("discrepancy_bias")]),
                    support: source_key(&format!("geometry-{target_case}")),
                    confounding_guard: source_key("fixture-discrepancy-guard"),
                },
            )],
            case.observation_sharing().to_vec(),
        )
        .expect("rebuilt modeled-discrepancy case");
        cases.insert(position, rebuilt);
    }
    if options.derived_cycle {
        sources.push(source(
            "derived-definition",
            SourceKind::Constraint,
            hash("derived-definition"),
        ));
        for (name, parent) in [("derived-a", "derived-b"), ("derived-b", "derived-a")] {
            parameters.push(
                StudyParameter::try_new(
                    role(name),
                    QuantitySpec::dimensional(DIMENSIONLESS),
                    ParameterDomain::try_new(0.0, 1.0).expect("derived domain"),
                    ParameterPurpose::Hyperparameter,
                    ParameterTreatment::Derived {
                        definition: source_key("derived-definition"),
                        parents: BTreeSet::from([role(parent)]),
                    },
                    ParameterOwnerBinding::Population {
                        hierarchy: source_key("derived-definition"),
                    },
                    ParameterScopeBinding::Global,
                    PriorPolicy::NotApplicable {
                        reason: "derived values do not own independent priors".to_string(),
                    },
                    InfluenceCoverage::Declared,
                )
                .expect("derived parameter fixture"),
            );
        }
    }
    let mut influences = vec![InfluenceDeclaration::new(
        InfluenceId::try_new("yield-to-stress").expect("influence id"),
        role("yield_stress"),
        if options.yield_influence_case_b {
            DistributionFunctional::Location {
                observation: ObservationKey::new(case_id("b"), channel("tangent")),
            }
        } else if options.yield_log_scale {
            DistributionFunctional::LogScale {
                observation: ObservationKey::new(case_id("a"), channel("stress")),
            }
        } else {
            DistributionFunctional::Location {
                observation: ObservationKey::new(case_id("a"), channel("stress")),
            }
        },
        InfluenceRepresentation::Direct,
    )];
    if !options.missing_hardening_influence {
        let tangent = ObservationKey::new(
            if options.bad_observation_endpoint {
                case_id("missing")
            } else {
                case_id("b")
            },
            channel("tangent"),
        );
        influences.push(InfluenceDeclaration::new(
            InfluenceId::try_new("hardening-to-tangent").expect("influence id"),
            role("hardening_modulus"),
            if options.self_correlation {
                DistributionFunctional::Correlation {
                    left: tangent.clone(),
                    right: tangent,
                }
            } else {
                DistributionFunctional::Location {
                    observation: tangent,
                }
            },
            InfluenceRepresentation::Direct,
        ));
    }
    if options.composite_influence_chain {
        sources.push(source(
            "composite-influence-operator",
            SourceKind::InfluenceComposition,
            hash("composite-influence-operator"),
        ));
        influences.extend([
            InfluenceDeclaration::new(
                InfluenceId::try_new("composite-middle").expect("influence id"),
                role("yield_stress"),
                DistributionFunctional::Location {
                    observation: ObservationKey::new(case_id("a"), channel("stress")),
                },
                InfluenceRepresentation::Composite {
                    operator: source_key("composite-influence-operator"),
                    inputs: BTreeSet::from([
                        InfluenceId::try_new("hardening-to-tangent").expect("influence id")
                    ]),
                },
            ),
            InfluenceDeclaration::new(
                InfluenceId::try_new("composite-top").expect("influence id"),
                role("yield_stress"),
                DistributionFunctional::Location {
                    observation: ObservationKey::new(case_id("a"), channel("stress")),
                },
                InfluenceRepresentation::Composite {
                    operator: source_key("composite-influence-operator"),
                    inputs: BTreeSet::from([
                        InfluenceId::try_new("composite-middle").expect("influence id")
                    ]),
                },
            ),
        ]);
    }
    let mut constraints = Vec::new();
    if options.valid_constraint {
        constraints.push(JointConstraint::new(
            ConstraintId::try_new("stress-balance").expect("constraint id"),
            JointConstraintKind::Affine {
                terms: vec![
                    AffineConstraintTerm::try_new(
                        role("yield_stress"),
                        1.0,
                        QuantitySpec::dimensional(DIMENSIONLESS),
                    )
                    .expect("yield term"),
                    AffineConstraintTerm::try_new(
                        role("hardening_modulus"),
                        -1.0,
                        QuantitySpec::dimensional(DIMENSIONLESS),
                    )
                    .expect("hardening term"),
                ],
                relation: ConstraintRelation::LessOrEqual,
                rhs_si: 0.0,
                residual_quantity: QuantitySpec::dimensional(STRESS),
            },
        ));
    }
    if options.ordered_constraint_case != 0 {
        constraints.push(JointConstraint::new(
            ConstraintId::try_new("ordered-stresses").expect("constraint id"),
            JointConstraintKind::Ordered {
                members: vec![role("yield_stress"), role("hardening_modulus")],
                strict: matches!(options.ordered_constraint_case, 3 | 4),
            },
        ));
    }
    if options.bad_constraint_units {
        constraints.push(JointConstraint::new(
            ConstraintId::try_new("bad-units").expect("constraint id"),
            JointConstraintKind::Affine {
                terms: vec![
                    AffineConstraintTerm::try_new(
                        role("yield_stress"),
                        1.0,
                        QuantitySpec::dimensional(DIMENSIONLESS),
                    )
                    .expect("term"),
                    AffineConstraintTerm::try_new(
                        role("hardening_modulus"),
                        -1.0,
                        QuantitySpec::dimensional(DIMENSIONLESS),
                    )
                    .expect("term"),
                ],
                relation: ConstraintRelation::Equal,
                rhs_si: 0.0,
                residual_quantity: QuantitySpec::dimensional(DIMENSIONLESS),
            },
        ));
    }
    let mut gauges = Vec::new();
    let mut gauge_compositions = Vec::new();
    if options.gauge_case != 0
        || options.one_gauge
        || options.claim_strata_in_problem
        || options.overlapping_gauges
    {
        sources.extend([
            source(
                "fixture-gauge-assumption",
                SourceKind::Assumption,
                hash("fixture-gauge-assumption"),
            ),
            source(
                "fixture-gauge-realization",
                SourceKind::QuantifierRealization,
                hash("fixture-gauge-realization"),
            ),
        ]);
    }
    if options.gauge_case != 0 {
        let gauge_action = if options.gauge_case == 11 {
            "fixture-unary-scaling-action"
        } else {
            "fixture-gauge-action"
        };
        sources.push(source(
            gauge_action,
            SourceKind::GaugeAction,
            hash(gauge_action),
        ));
        let members = if options.gauge_case == 11 {
            BTreeSet::from([role("yield_stress")])
        } else {
            BTreeSet::from([role("yield_stress"), role("hardening_modulus")])
        };
        let (algebra, orbit_geometry, local_obstruction) = match options.gauge_case {
            1 => (
                GaugeAlgebra::Discrete {
                    size: GaugeDiscreteSize::Finite { order: 2 },
                },
                GaugeOrbitGeometry::Regular {
                    principal: RegularGaugeOrbit::new(
                        GaugeContinuousDimension::Finite { dimension: 0 },
                        GaugeDiscreteOrbitCardinality::Finite { cardinality: 2 },
                    ),
                    stabilizer_profile: None,
                },
                false,
            ),
            2 | 3 | 7 => (
                GaugeAlgebra::Mixed {
                    continuous_group_dimension: GaugeContinuousDimension::Finite { dimension: 1 },
                    component_group: GaugeDiscreteSize::Finite { order: 2 },
                },
                GaugeOrbitGeometry::Regular {
                    principal: RegularGaugeOrbit::new(
                        GaugeContinuousDimension::Finite { dimension: 1 },
                        GaugeDiscreteOrbitCardinality::Finite { cardinality: 2 },
                    ),
                    stabilizer_profile: None,
                },
                true,
            ),
            9 => {
                sources.push(source(
                    "fixture-gauge-strata",
                    SourceKind::GaugeOrbitTypeProfile,
                    hash("fixture-gauge-strata"),
                ));
                (
                    GaugeAlgebra::Continuous {
                        group_dimension: GaugeContinuousDimension::Finite { dimension: 1 },
                    },
                    GaugeOrbitGeometry::Stratified {
                        principal: RegularGaugeOrbit::new(
                            GaugeContinuousDimension::Finite { dimension: 1 },
                            GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 },
                        ),
                        orbit_type_stabilizer_profile: source_key("fixture-gauge-strata"),
                    },
                    true,
                )
            }
            mode @ (4..=6 | 8 | 10 | 11) => {
                let dimension = if mode == 10 { 2 } else { 1 };
                (
                    GaugeAlgebra::Continuous {
                        group_dimension: GaugeContinuousDimension::Finite { dimension },
                    },
                    GaugeOrbitGeometry::Regular {
                        principal: RegularGaugeOrbit::new(
                            GaugeContinuousDimension::Finite { dimension },
                            GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 },
                        ),
                        stabilizer_profile: None,
                    },
                    true,
                )
            }
            other => panic!("unsupported gauge fixture {other}"),
        };
        gauges.push(gauge_fixture(
            "fixture-gauge",
            gauge_action,
            members,
            algebra,
            orbit_geometry,
            local_obstruction,
        ));
    }
    if options.one_gauge {
        sources.push(source(
            "single-gauge-action",
            SourceKind::GaugeAction,
            hash("single-gauge-action"),
        ));
        gauges.push(gauge_fixture(
            "single-gauge",
            "single-gauge-action",
            BTreeSet::from([role("yield_stress"), role("hardening_modulus")]),
            GaugeAlgebra::Continuous {
                group_dimension: GaugeContinuousDimension::Finite { dimension: 1 },
            },
            GaugeOrbitGeometry::Regular {
                principal: RegularGaugeOrbit::new(
                    GaugeContinuousDimension::Finite { dimension: 1 },
                    GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 },
                ),
                stabilizer_profile: None,
            },
            true,
        ));
    }
    if options.claim_strata_in_problem {
        sources.extend([
            source(
                "claim-domain-action",
                SourceKind::GaugeAction,
                hash("claim-domain-action"),
            ),
            source(
                "claim-strata",
                SourceKind::GaugeOrbitTypeProfile,
                hash("claim-strata"),
            ),
        ]);
        gauges.push(gauge_fixture(
            "claim-domain-gauge",
            "claim-domain-action",
            BTreeSet::from([role("yield_stress"), role("hardening_modulus")]),
            GaugeAlgebra::Continuous {
                group_dimension: GaugeContinuousDimension::Finite { dimension: 1 },
            },
            GaugeOrbitGeometry::Stratified {
                principal: RegularGaugeOrbit::new(
                    GaugeContinuousDimension::Finite { dimension: 1 },
                    GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 },
                ),
                orbit_type_stabilizer_profile: source_key("claim-strata"),
            },
            true,
        ));
    }
    if options.overlapping_gauges {
        for index in 0..2 {
            let action = format!("gauge-action-{index}");
            sources.push(source(&action, SourceKind::GaugeAction, hash(&action)));
            gauges.push(gauge_fixture(
                &format!("gauge-{index}"),
                &action,
                BTreeSet::from([role("yield_stress"), role("hardening_modulus")]),
                GaugeAlgebra::Continuous {
                    group_dimension: GaugeContinuousDimension::Finite { dimension: 1 },
                },
                GaugeOrbitGeometry::Regular {
                    principal: RegularGaugeOrbit::new(
                        GaugeContinuousDimension::Finite { dimension: 1 },
                        GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 },
                    ),
                    stabilizer_profile: None,
                },
                true,
            ));
        }
    }
    if options.declared_gauge_composition {
        assert!(
            options.overlapping_gauges,
            "a composition fixture needs the two overlapping member gauges"
        );
        sources.push(source(
            "gauge-composition-law",
            SourceKind::GaugeComposition,
            hash("gauge-composition-law"),
        ));
        let effective_dimension = if options.independent_gauge_composition {
            2
        } else {
            1
        };
        gauge_compositions.push(
            GaugeCompositionDeclaration::try_new(
                GaugeCompositionId::try_new("gauge-system").expect("composition id"),
                BTreeSet::from([
                    GaugeClassId::try_new("gauge-0").expect("gauge id"),
                    GaugeClassId::try_new("gauge-1").expect("gauge id"),
                ]),
                if options.independent_gauge_composition {
                    GaugeCompositionKind::IndependentProduct
                } else {
                    GaugeCompositionKind::Generated
                },
                source_key("gauge-composition-law"),
                GaugeAlgebra::Continuous {
                    group_dimension: GaugeContinuousDimension::Finite {
                        dimension: effective_dimension,
                    },
                },
                GaugeOrbitGeometry::Regular {
                    principal: RegularGaugeOrbit::new(
                        GaugeContinuousDimension::Finite {
                            dimension: effective_dimension,
                        },
                        GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 },
                    ),
                    stabilizer_profile: None,
                },
                GaugeStatus::Assumed {
                    assumption: source_key("fixture-gauge-assumption"),
                },
                gauge_validity_fixture(
                    &BTreeSet::from([role("yield_stress"), role("hardening_modulus")]),
                    true,
                ),
            )
            .expect("declared gauge composition fixture"),
        );
    }
    let joint_noise = if options.dense_with_bounded_marginal {
        sources.push(source(
            "correlation-model",
            SourceKind::Likelihood,
            hash("correlation-model"),
        ));
        JointNoiseModel::DenseCorrelation {
            order: vec![
                ObservationKey::new(case_id("a"), channel("stress")),
                ObservationKey::new(case_id("b"), channel("tangent")),
            ],
            correlation: CovarianceMatrix::try_new(2, vec![1.0, 0.0, 1.0])
                .expect("correlation matrix"),
            model: source_key("correlation-model"),
        }
    } else if options.declared_sharing {
        JointNoiseModel::ExternalKernel {
            model: source_key("joint-likelihood"),
        }
    } else if options.external_noise {
        sources.push(source(
            "external-noise",
            SourceKind::Likelihood,
            hash("external-noise"),
        ));
        JointNoiseModel::ExternalKernel {
            model: source_key("external-noise"),
        }
    } else {
        sources.push(source(
            "independent-noise-assumption",
            SourceKind::Assumption,
            hash("independent-noise-assumption"),
        ));
        JointNoiseModel::Independent {
            assumption: source_key("independent-noise-assumption"),
        }
    };
    let data_reuse = if options.retrospective_reuse && options.declared_sharing {
        DataReusePolicy::Shared {
            groups: vec![
                DataSharingGroup::try_new(
                    BTreeSet::from([case_id("a"), case_id("b")]),
                    source_key("joint-likelihood"),
                    if options.alternate_sharing_justification {
                        "the exact same shared campaign is justified by an alternate audit note"
                    } else {
                        "the cases intentionally derive complementary channels from one raw campaign"
                    },
                )
                .expect("sharing group"),
            ],
        }
    } else {
        DataReusePolicy::Disjoint
    };
    let admissible_domain = AdmissibleDomainWitness::try_new(
        parameters
            .iter()
            .map(|parameter| (parameter.role().clone(), parameter.domain().bounds().0))
            .collect(),
        (options.gauge_case == 10).then(|| {
            OpaqueDomainMembershipClaim::new(source_key("fixture-domain-membership-certificate"))
        }),
    )
    .expect("constructive admissible-domain witness");
    if options.gauge_case == 10 {
        let parameter_map = parameters
            .iter()
            .cloned()
            .map(|parameter| (parameter.role().clone(), parameter))
            .collect();
        let constraint_map = constraints
            .iter()
            .cloned()
            .map(|constraint| (constraint.id().clone(), constraint))
            .collect();
        let source_map = sources
            .iter()
            .cloned()
            .map(|source| (source.key().clone(), source))
            .collect();
        let certificate_preimage = admissible_domain_membership_certificate_preimage(
            &admissible_domain,
            &parameter_map,
            &constraint_map,
            &source_map,
        )
        .expect("canonical membership-certificate preimage");
        sources.push(
            SourceRef::try_new(
                source_key("fixture-domain-membership-certificate"),
                SourceKind::AdmissibleDomainCertificate,
                hash_domain(
                    ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_DOMAIN,
                    &certificate_preimage,
                ),
                ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_DOMAIN,
                ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_VERSION,
            )
            .expect("witness-bound membership certificate"),
        );
    }
    let joint_prior = match options.joint_prior_choice {
        0 => None,
        1 => Some(source_key("joint-prior-measure-a")),
        2 => Some(source_key("joint-prior-measure-b")),
        _ => unreachable!("joint-prior fixture choice checked above"),
    };
    let document = IdentifiabilityProblemDocument::try_new(
        source_key("context"),
        source_key("material"),
        source_key("model"),
        source_key("graph"),
        joint_prior,
        sources,
        parameters,
        constraints,
        admissible_domain,
        cases,
        influences,
        gauges,
        gauge_compositions,
        joint_noise,
        data_reuse,
    );
    ProblemFixture {
        context,
        material,
        model,
        graph,
        document,
    }
}

fn rekey_problem_root(
    fixture: ProblemFixture,
    root: ProblemRoot,
    replacement: &str,
) -> ProblemFixture {
    let ProblemFixture {
        context,
        material,
        model,
        graph,
        document,
    } = fixture;
    let document = document.expect("problem before root rekey");
    let old_key = match root {
        ProblemRoot::Context => document.context_source(),
        ProblemRoot::Material => document.material_source(),
        ProblemRoot::Model => document.model_source(),
        ProblemRoot::Graph => document.graph_source(),
    }
    .clone();
    let replacement = source_key(replacement);
    let sources = document
        .sources()
        .values()
        .map(|source| {
            if source.key() == &old_key {
                SourceRef::try_new(
                    replacement.clone(),
                    source.kind(),
                    source.expected_hash(),
                    source.content_hash_domain(),
                    source.contract_version(),
                )
                .expect("rekeyed source")
            } else {
                source.clone()
            }
        })
        .collect();
    let context_source = if matches!(root, ProblemRoot::Context) {
        replacement.clone()
    } else {
        document.context_source().clone()
    };
    let material_source = if matches!(root, ProblemRoot::Material) {
        replacement.clone()
    } else {
        document.material_source().clone()
    };
    let model_source = if matches!(root, ProblemRoot::Model) {
        replacement.clone()
    } else {
        document.model_source().clone()
    };
    let graph_source = if matches!(root, ProblemRoot::Graph) {
        replacement
    } else {
        document.graph_source().clone()
    };
    let document = IdentifiabilityProblemDocument::try_new(
        context_source,
        material_source,
        model_source,
        graph_source,
        document.joint_prior().cloned(),
        sources,
        document.parameters().values().cloned().collect(),
        document.constraints().values().cloned().collect(),
        document.admissible_domain().clone(),
        document.cases().values().cloned().collect(),
        document.influences().values().cloned().collect(),
        document.gauges().values().cloned().collect(),
        document.gauge_compositions().values().cloned().collect(),
        document.joint_noise().clone(),
        document.data_reuse().clone(),
    );
    ProblemFixture {
        context,
        material,
        model,
        graph,
        document,
    }
}

fn unresolved_problem_identity(document: &IdentifiabilityProblemDocument) -> ContentHash {
    hash_domain(
        IDENTIFIABILITY_PROBLEM_IDENTITY_DOMAIN,
        &document.canonical_bytes().expect("canonical problem bytes"),
    )
}

fn opaque_source_preimage(
    document: &IdentifiabilityProblemDocument,
    source: &SourceRef,
) -> Vec<u8> {
    if source.kind() == SourceKind::AdmissibleDomainCertificate {
        document
            .admissible_domain()
            .opaque_membership_claim()
            .and_then(OpaqueDomainMembershipClaim::witness_binding)
            .expect("bound admissible-domain membership certificate")
            .as_bytes()
            .to_vec()
    } else if source.kind() == SourceKind::ForwardModelProductionBinding {
        document
            .cases()
            .values()
            .flat_map(|case| case.discrepancies().values())
            .find_map(|discrepancy| match discrepancy {
                StudyDiscrepancy::NotApplicable {
                    basis:
                        DiscrepancyInapplicability::DeclaredSyntheticSelfModel {
                            generator,
                            producer,
                            production_binding,
                            ..
                        },
                } if production_binding == source.key() => Some(
                    forward_model_production_binding_preimage(
                        producer,
                        &document.sources()[generator],
                    )
                    .expect("exact production-binding preimage"),
                ),
                _ => None,
            })
            .expect("reachable production-binding source")
    } else if source.kind() == SourceKind::ConstitutiveGraph {
        b"constitutive-graph".to_vec()
    } else {
        source.key().as_str().as_bytes().to_vec()
    }
}

fn opaque_resolutions(document: &IdentifiabilityProblemDocument) -> SourceResolutionSet {
    let entries = document
        .sources()
        .values()
        .filter(|source| {
            !matches!(
                source.kind(),
                SourceKind::ContextOfUse
                    | SourceKind::MaterialCard
                    | SourceKind::ConstitutiveModelCard
                    | SourceKind::ExperimentArtifact
                    | SourceKind::CalibrationSplit
            )
        })
        .map(|source| {
            SourceResolution::verify(
                source,
                &opaque_source_preimage(document, source),
                AuthorityDisposition::ContentVerified,
            )
            .expect("opaque resolution fixture")
        })
        .collect();
    SourceResolutionSet::try_new(entries).expect("resolution set fixture")
}

fn admit_fixture(fixture: ProblemFixture) -> AdmittedIdentifiabilityProblem {
    admit_fixture_with_authority(fixture, false)
}

fn admit_fixture_with_authority(
    fixture: ProblemFixture,
    authenticated: bool,
) -> AdmittedIdentifiabilityProblem {
    admit_fixture_with_authority_and_order(fixture, authenticated, false)
}

fn admit_fixture_with_authority_and_order(
    fixture: ProblemFixture,
    authenticated: bool,
    reverse_resolutions: bool,
) -> AdmittedIdentifiabilityProblem {
    let ProblemFixture {
        context,
        material,
        model,
        graph: _,
        document,
    } = fixture;
    let document = document.expect("problem structurally admits");
    let mut entries = if authenticated {
        document
            .sources()
            .values()
            .filter(|source| {
                !matches!(
                    source.kind(),
                    SourceKind::ContextOfUse
                        | SourceKind::MaterialCard
                        | SourceKind::ConstitutiveModelCard
                        | SourceKind::ExperimentArtifact
                        | SourceKind::CalibrationSplit
                )
            })
            .map(|source| {
                SourceResolution::verify(
                    source,
                    &opaque_source_preimage(&document, source),
                    external_trust(&format!("trust-{}", source.key()), source),
                )
                .expect("authenticated resolution")
            })
            .collect::<Vec<_>>()
    } else {
        opaque_resolutions(&document)
            .entries()
            .values()
            .cloned()
            .collect::<Vec<_>>()
    };
    if reverse_resolutions {
        entries.reverse();
    }
    let opaque = SourceResolutionSet::try_new(entries).expect("canonical resolution set");
    AdmittedIdentifiabilityProblem::resolve_and_admit(
        document,
        ProblemSourceBundle::new(&context, &material, &model, BTreeMap::new(), opaque),
    )
    .expect("source-resolved problem admits")
}

fn admit_fixture_with_single_external_authority(
    fixture: ProblemFixture,
    trusted_key: &str,
) -> AdmittedIdentifiabilityProblem {
    let ProblemFixture {
        context,
        material,
        model,
        graph: _,
        document,
    } = fixture;
    let document = document.expect("problem structurally admits");
    let mut changed = 0;
    let entries = document
        .sources()
        .values()
        .filter(|source| {
            !matches!(
                source.kind(),
                SourceKind::ContextOfUse
                    | SourceKind::MaterialCard
                    | SourceKind::ConstitutiveModelCard
                    | SourceKind::ExperimentArtifact
                    | SourceKind::CalibrationSplit
            )
        })
        .map(|source| {
            let authority = if source.key().as_str() == trusted_key {
                changed += 1;
                external_trust(&format!("trust-{trusted_key}"), source)
            } else {
                AuthorityDisposition::ContentVerified
            };
            SourceResolution::verify(
                source,
                &opaque_source_preimage(&document, source),
                authority,
            )
            .expect("single-authority resolution")
        })
        .collect::<Vec<_>>();
    assert_eq!(
        changed, 1,
        "exactly one requested source authority must move"
    );
    let opaque = SourceResolutionSet::try_new(entries).expect("single-authority resolution set");
    AdmittedIdentifiabilityProblem::resolve_and_admit(
        document,
        ProblemSourceBundle::new(&context, &material, &model, BTreeMap::new(), opaque),
    )
    .expect("single-authority problem admits")
}

fn execution_header(seed: u64) -> ArtifactHeader {
    ArtifactHeader::try_new(
        artifact("execution-1"),
        vec![unit("Pa")],
        SeedDeclaration::Fixed(seed),
        DeclaredBudget::Limit(1.0e-9),
        DeclaredBudget::Limit(30_000),
        DeclaredBudget::Limit(32 << 20),
        vec![("fixture".to_string(), "1".to_string())],
        vec!["identifiability.execute".to_string()],
    )
    .expect("execution header")
}

fn execution_header_with_semantics(
    units: &[&str],
    seed: u64,
    accuracy: f64,
    time_ms: u64,
    memory_bytes: u64,
    version: &str,
    extra_capability: bool,
) -> ArtifactHeader {
    let mut capabilities = vec!["identifiability.execute".to_string()];
    if extra_capability {
        capabilities.push("identifiability.symbolic".to_string());
    }
    ArtifactHeader::try_new(
        artifact("execution-1"),
        units.iter().map(|value| unit(value)).collect(),
        SeedDeclaration::Fixed(seed),
        DeclaredBudget::Limit(accuracy),
        DeclaredBudget::Limit(time_ms),
        DeclaredBudget::Limit(memory_bytes),
        vec![("fixture".to_string(), version.to_string())],
        capabilities,
    )
    .expect("execution semantic header")
}

fn assessment_header_with_semantics(
    units: &[&str],
    seed: u64,
    accuracy: f64,
    time_ms: u64,
    memory_bytes: u64,
    version: &str,
    extra_capability: bool,
) -> ArtifactHeader {
    let mut capabilities = vec!["identifiability.assess".to_string()];
    if extra_capability {
        capabilities.push("identifiability.explain".to_string());
    }
    ArtifactHeader::try_new(
        artifact("assessment-1"),
        units.iter().map(|value| unit(value)).collect(),
        SeedDeclaration::Fixed(seed),
        DeclaredBudget::Limit(accuracy),
        DeclaredBudget::Limit(time_ms),
        DeclaredBudget::Limit(memory_bytes),
        vec![("fixture".to_string(), version.to_string())],
        capabilities,
    )
    .expect("assessment semantic header")
}

fn coordinate(name: &str, affine: bool) -> ParameterCoordinate {
    let domain = if name == "yield_stress" {
        ParameterDomain::try_new(1.0e6, 1.0e9).expect("domain")
    } else {
        ParameterDomain::try_new(1.0e7, 5.0e9).expect("domain")
    };
    let (quantity, coordinate_domain, transform, suffix) = if affine {
        let (lo, hi) = domain.bounds();
        (
            QuantitySpec::dimensional(DIMENSIONLESS),
            ParameterDomain::try_new(lo / 1.0e6, hi / 1.0e6).expect("scaled domain"),
            CoordinateTransform::Affine {
                scale: 1.0e6,
                scale_quantity: QuantitySpec::dimensional(STRESS),
                offset: 0.0,
            },
            "mpa",
        )
    } else {
        (
            QuantitySpec::dimensional(STRESS),
            domain,
            CoordinateTransform::Identity,
            "si",
        )
    };
    ParameterCoordinate::try_new(
        CoordinateId::try_new(format!("{name}-{suffix}")).expect("coordinate id"),
        quantity,
        coordinate_domain,
        transform,
    )
    .expect("coordinate fixture")
}

fn resolved_sources(sources: &[&SourceRef], external_trust: bool) -> SourceResolutionSet {
    SourceResolutionSet::try_new(
        sources
            .iter()
            .map(|source| {
                let authority = if external_trust {
                    external_trust(&format!("execution-trust-{}", source.key()), source)
                } else {
                    AuthorityDisposition::ContentVerified
                };
                SourceResolution::verify(source, source.key().as_str().as_bytes(), authority)
                    .expect("content-verified source fixture")
            })
            .collect(),
    )
    .expect("closed content-verified source set")
}

fn default_quantifier_domain(problem: &AdmittedIdentifiabilityProblem) -> SourceRef {
    problem
        .document()
        .sources()
        .get(&source_key("claim-domain"))
        .cloned()
        .unwrap_or_else(|| {
            source(
                "claim-domain",
                SourceKind::QuantifierDomain,
                hash("claim-domain"),
            )
        })
}

fn fixture_claim_quantifier(problem: &AdmittedIdentifiabilityProblem) -> ClaimQuantifier {
    problem
        .document()
        .sources()
        .get(&source_key("fixture-gauge-realization"))
        .map_or_else(
            || ClaimQuantifier::ForAll {
                domain: default_quantifier_domain(problem),
            },
            |realization| ClaimQuantifier::AtRealization {
                realization: realization.clone(),
            },
        )
}

fn default_claim(problem: &AdmittedIdentifiabilityProblem) -> TypedIdentifiabilityClaim {
    let claimed_role = role("yield_stress");
    let fiber = if let Some(strata) = problem
        .document()
        .sources()
        .get(&source_key("claim-strata"))
    {
        FiberStructure::Stratified {
            strata: strata.clone(),
        }
    } else if problem
        .document()
        .gauges()
        .values()
        .any(|gauge| gauge.members().contains(&claimed_role))
    {
        FiberStructure::PositiveDimensional {
            lower_bound: FiberDimensionLowerBound::Finite {
                minimum_dimension: 1,
            },
        }
    } else {
        FiberStructure::Unique
    };
    TypedIdentifiabilityClaim::new(
        ClaimId::try_new("yield-structural-global").expect("claim id"),
        InformationRegime::StructuralExactModel,
        IdentifiabilityExtent::Global,
        fiber,
        fixture_claim_quantifier(problem),
        ScalarDomain::Real,
        ClaimSubject::Parameter(claimed_role),
        ClaimScope::WholeCampaign,
    )
}

fn structural_claim(
    problem: &AdmittedIdentifiabilityProblem,
    id: &str,
    fiber: FiberStructure,
    subject: ClaimSubject,
    scope: ClaimScope,
) -> TypedIdentifiabilityClaim {
    TypedIdentifiabilityClaim::new(
        ClaimId::try_new(id).expect("claim id"),
        InformationRegime::StructuralExactModel,
        IdentifiabilityExtent::Global,
        fiber,
        fixture_claim_quantifier(problem),
        ScalarDomain::Real,
        subject,
        scope,
    )
}

fn default_error_policy() -> DimensionlessErrorPolicy {
    DimensionlessErrorPolicy::try_new(
        source(
            "claim-error-metric",
            SourceKind::DimensionlessErrorMetric,
            hash("claim-error-metric"),
        ),
        source(
            "claim-nondimensionalization",
            SourceKind::Nondimensionalization,
            hash("claim-nondimensionalization"),
        ),
        1.0e-8,
    )
    .expect("dimensionless claim error policy")
}

fn default_claim_request(problem: &AdmittedIdentifiabilityProblem) -> ClaimRequest {
    ClaimRequest::new(default_claim(problem), default_error_policy())
}

fn request_for_claim(claim: TypedIdentifiabilityClaim) -> ClaimRequest {
    ClaimRequest::new(claim, default_error_policy())
}

fn claim_sources(claim: &TypedIdentifiabilityClaim) -> Vec<SourceRef> {
    let mut sources = Vec::new();
    if let InformationRegime::PosteriorUnderDeclaredPrior { joint_prior } = claim.information() {
        sources.push(joint_prior.clone());
    }
    match claim.scalar_domain() {
        ScalarDomain::Real => {}
        ScalarDomain::Complex { extension } => sources.push(extension.clone()),
        ScalarDomain::MixedDiscreteContinuous { stratification } => {
            sources.push(stratification.clone());
        }
    }
    if let FiberStructure::Stratified { strata } = claim.fiber() {
        sources.push(strata.clone());
    }
    if let ClaimSubject::DerivedFunctional { definition, .. } = claim.subject() {
        sources.push(definition.clone());
    }
    sources.push(match claim.quantifier() {
        ClaimQuantifier::AtRealization { realization } => realization.clone(),
        ClaimQuantifier::AlmostEverywhere { measure }
        | ClaimQuantifier::ProbabilityAtLeast { measure, .. } => measure.clone(),
        ClaimQuantifier::ForAll { domain } => domain.clone(),
    });
    sources
}

fn resolve_owned_sources(
    sources: impl IntoIterator<Item = SourceRef>,
    external_trust: bool,
    authority_override: Option<(&SourceKey, AuthorityDisposition)>,
) -> SourceResolutionSet {
    let mut unique = BTreeMap::<SourceKey, SourceRef>::new();
    for source in sources {
        if let Some(prior) = unique.insert(source.key().clone(), source.clone()) {
            assert_eq!(
                prior, source,
                "fixture source aliases must have exactly equal semantics",
            );
        }
    }
    SourceResolutionSet::try_new(
        unique
            .values()
            .map(|source| {
                let authority = authority_override
                    .as_ref()
                    .filter(|(key, _)| source.key() == *key)
                    .map(|(_, authority)| authority.clone())
                    .unwrap_or_else(|| {
                        if external_trust {
                            external_trust(&format!("execution-trust-{}", source.key()), source)
                        } else {
                            AuthorityDisposition::ContentVerified
                        }
                    });
                SourceResolution::verify(source, source.key().as_str().as_bytes(), authority)
                    .expect("content-verified owned source fixture")
            })
            .collect(),
    )
    .expect("deduplicated source authority")
}

fn execution(
    problem: &AdmittedIdentifiabilityProblem,
    affine: bool,
    seed: u64,
    tolerance: f64,
    wrong_action: bool,
) -> Result<IdentifiabilityExecutionPlan, IdentifiabilityError> {
    execution_with_claim_requests_and_authority(
        problem,
        affine,
        seed,
        tolerance,
        wrong_action,
        vec![default_claim_request(problem)],
        Vec::new(),
        false,
    )
}

fn execution_for_claim(
    problem: &AdmittedIdentifiabilityProblem,
    affine: bool,
    seed: u64,
    tolerance: f64,
    wrong_action: bool,
    claim: TypedIdentifiabilityClaim,
) -> Result<IdentifiabilityExecutionPlan, IdentifiabilityError> {
    execution_with_claim_requests_and_authority(
        problem,
        affine,
        seed,
        tolerance,
        wrong_action,
        vec![request_for_claim(claim)],
        Vec::new(),
        false,
    )
}

fn execution_for_claim_with_gauge_reductions(
    problem: &AdmittedIdentifiabilityProblem,
    claim: TypedIdentifiabilityClaim,
    gauge_reductions: Vec<GaugeReductionBinding>,
) -> Result<IdentifiabilityExecutionPlan, IdentifiabilityError> {
    execution_with_claim_requests_and_authority(
        problem,
        false,
        17,
        1.0e-10,
        false,
        vec![request_for_claim(claim)],
        gauge_reductions,
        false,
    )
}

fn gauge_reduction_authority_sources(binding: &GaugeReductionBinding) -> Vec<SourceRef> {
    fn quotient_sources(quotient: &GaugeQuotientPlan, sources: &mut Vec<SourceRef>) {
        match quotient {
            GaugeQuotientPlan::RegularAtlas {
                quotient_map,
                local_section_atlas,
                coverage,
            } => sources.extend([
                quotient_map.clone(),
                local_section_atlas.clone(),
                coverage.clone(),
            ]),
            GaugeQuotientPlan::SingularOrGeneralized {
                quotient_map,
                quotient_profile,
                local_models,
            } => {
                sources.extend([quotient_map.clone(), quotient_profile.clone()]);
                sources.extend(local_models.iter().cloned());
            }
            GaugeQuotientPlan::InvariantMap {
                invariants,
                completeness_profile,
            } => sources.extend([invariants.clone(), completeness_profile.clone()]),
            GaugeQuotientPlan::GroupoidOrStack {
                presentation,
                quotient_profile,
            } => sources.extend([presentation.clone(), quotient_profile.clone()]),
        }
    }

    fn slice_sources(slice: &GaugeSlicePlan, sources: &mut Vec<SourceRef>) {
        sources.extend([slice.constraint().clone(), slice.coverage().clone()]);
        match slice.expected_codimension() {
            GaugeSliceCodimension::FixedFinite { .. } => {}
            GaugeSliceCodimension::FixedInfinite {
                codimension_model,
                compatibility,
            } => sources.extend([codimension_model.clone(), compatibility.clone()]),
            GaugeSliceCodimension::Stratified { profile } => sources.push(profile.clone()),
        }
    }

    fn continuous_sources(reduction: &ContinuousGaugeReductionPlan, sources: &mut Vec<SourceRef>) {
        match reduction {
            ContinuousGaugeReductionPlan::Quotient { quotient } => {
                quotient_sources(quotient, sources);
            }
            ContinuousGaugeReductionPlan::Slice { slice } => slice_sources(slice, sources),
        }
    }

    let mut sources = Vec::new();
    match binding.plan() {
        GaugeReductionPlan::Unreduced { .. } => {}
        GaugeReductionPlan::Quotient { quotient } => quotient_sources(quotient, &mut sources),
        GaugeReductionPlan::Slice { slice } => slice_sources(slice, &mut sources),
        GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
            reduction,
            normal_subgroup,
            factor_extension,
            residual_quotient_action,
            compatibility,
        } => {
            continuous_sources(reduction, &mut sources);
            sources.extend([
                normal_subgroup.clone(),
                factor_extension.clone(),
                residual_quotient_action.clone(),
                compatibility.clone(),
            ]);
        }
    }
    if let GaugeReductionStage::After {
        composition_law,
        relation,
        ..
    } = binding.stage()
    {
        sources.push(composition_law.clone());
        match relation {
            GaugeReductionStageRelation::NormalSubgroupTower {
                normality,
                induced_residual_action,
            } => sources.extend([normality.clone(), induced_residual_action.clone()]),
            GaugeReductionStageRelation::SemidirectOrGenerated {
                extension,
                induced_residual_action,
            } => sources.extend([extension.clone(), induced_residual_action.clone()]),
            GaugeReductionStageRelation::TransverseSlices { transversality } => {
                sources.push(transversality.clone());
            }
            GaugeReductionStageRelation::GaugeForGauge {
                reducibility,
                induced_residual_action,
            } => sources.extend([reducibility.clone(), induced_residual_action.clone()]),
        }
    }
    if let GaugeMeasureSemantics::Pushforward {
        source_measure,
        reduced_measure,
        transport,
        jacobian_or_disintegration,
    } = binding.measure()
    {
        sources.extend([
            source_measure.clone(),
            reduced_measure.clone(),
            transport.clone(),
            jacobian_or_disintegration.clone(),
        ]);
    }
    sources
}

#[allow(clippy::too_many_arguments)]
fn execution_with_claim_requests_and_authority(
    problem: &AdmittedIdentifiabilityProblem,
    affine: bool,
    seed: u64,
    tolerance: f64,
    wrong_action: bool,
    claim_requests: Vec<ClaimRequest>,
    gauge_reductions: Vec<GaugeReductionBinding>,
    external_trust: bool,
) -> Result<IdentifiabilityExecutionPlan, IdentifiabilityError> {
    let analyzer = source("analyzer", SourceKind::Analyzer, hash("analyzer"));
    let build = source("build", SourceKind::Build, hash("build"));
    let derivatives = source(
        "derivatives",
        SourceKind::DerivativeProvider,
        hash("derivatives"),
    );
    let quadrature = source("quadrature", SourceKind::Analyzer, hash("quadrature"));
    let measure_transport = source(
        "hardening-measure-transport",
        SourceKind::MeasureTransport,
        hash("hardening-measure-transport"),
    );
    let initialization = source(
        "initialization",
        SourceKind::Assumption,
        hash("initialization"),
    );
    let stopping = source("stopping", SourceKind::Assumption, hash("stopping"));
    let determinism = source("determinism", SourceKind::Assumption, hash("determinism"));
    let numerical_nondimensionalization = source(
        "numerical-nondimensionalization",
        SourceKind::Nondimensionalization,
        hash("numerical-nondimensionalization"),
    );
    let mut authority_sources = vec![
        analyzer.clone(),
        build.clone(),
        derivatives.clone(),
        initialization.clone(),
        stopping.clone(),
        determinism.clone(),
        numerical_nondimensionalization.clone(),
    ];
    if !wrong_action {
        authority_sources.push(quadrature.clone());
        authority_sources.push(measure_transport.clone());
    }
    for request in &claim_requests {
        authority_sources.extend(claim_sources(request.claim()));
        authority_sources.push(request.error_policy().metric().clone());
        authority_sources.push(request.error_policy().nondimensionalization().clone());
    }
    for reduction in &gauge_reductions {
        authority_sources.extend(gauge_reduction_authority_sources(reduction));
    }
    let source_authority = resolve_owned_sources(authority_sources, external_trust, None);
    IdentifiabilityExecutionPlan::try_new(
        execution_header(seed),
        problem,
        analyzer,
        build,
        Some(derivatives),
        claim_requests,
        vec![
            (
                role("yield_stress"),
                ParameterExecutionAction::Optimize {
                    coordinate: coordinate("yield_stress", affine),
                },
            ),
            (
                role("hardening_modulus"),
                if wrong_action {
                    ParameterExecutionAction::Optimize {
                        coordinate: coordinate("hardening_modulus", affine),
                    }
                } else {
                    ParameterExecutionAction::Marginalize {
                        coordinate: coordinate("hardening_modulus", affine),
                        integrator: quadrature,
                        measure_transport,
                    }
                },
            ),
        ],
        gauge_reductions,
        IdentifiabilityNumericalPolicy::try_new(
            tolerance,
            0.0,
            1.0e12,
            ArithmeticPolicy::CertifiedInterval,
            numerical_nondimensionalization,
        )?,
        initialization,
        stopping,
        determinism,
        source_authority,
    )
}

#[derive(Clone)]
struct ExecutionParts {
    header: ArtifactHeader,
    analyzer: SourceRef,
    build: SourceRef,
    derivative_provider: Option<SourceRef>,
    claim_requests: Vec<ClaimRequest>,
    actions: Vec<(ParameterRoleId, ParameterExecutionAction)>,
    gauge_reductions: Vec<GaugeReductionBinding>,
    numerical: IdentifiabilityNumericalPolicy,
    initialization: SourceRef,
    stopping: SourceRef,
    determinism: SourceRef,
}

impl ExecutionParts {
    fn from_plan(plan: &IdentifiabilityExecutionPlan) -> Self {
        Self {
            header: plan.header().clone(),
            analyzer: plan.analyzer().clone(),
            build: plan.build().clone(),
            derivative_provider: plan.derivative_provider().cloned(),
            claim_requests: plan.claim_requests().values().cloned().collect(),
            actions: plan
                .actions()
                .iter()
                .map(|(role, action)| (role.clone(), action.clone()))
                .collect(),
            gauge_reductions: plan.gauge_reductions().values().cloned().collect(),
            numerical: plan.numerical_policy().clone(),
            initialization: plan.initialization().clone(),
            stopping: plan.stopping().clone(),
            determinism: plan.determinism_contract().clone(),
        }
    }

    fn build(
        self,
        problem: &AdmittedIdentifiabilityProblem,
        external_trust: bool,
    ) -> IdentifiabilityExecutionPlan {
        let source_authority = {
            let mut sources = vec![
                self.analyzer.clone(),
                self.build.clone(),
                self.initialization.clone(),
                self.stopping.clone(),
                self.determinism.clone(),
                self.numerical.nondimensionalization().clone(),
            ];
            if let Some(provider) = &self.derivative_provider {
                sources.push(provider.clone());
            }
            for (_, action) in &self.actions {
                if let ParameterExecutionAction::Marginalize {
                    integrator,
                    measure_transport,
                    ..
                } = action
                {
                    sources.push(integrator.clone());
                    sources.push(measure_transport.clone());
                }
            }
            for request in &self.claim_requests {
                sources.extend(claim_sources(request.claim()));
                sources.push(request.error_policy().metric().clone());
                sources.push(request.error_policy().nondimensionalization().clone());
            }
            for reduction in &self.gauge_reductions {
                sources.extend(gauge_reduction_authority_sources(reduction));
            }
            resolve_owned_sources(sources, external_trust, None)
        };
        IdentifiabilityExecutionPlan::try_new(
            self.header,
            problem,
            self.analyzer,
            self.build,
            self.derivative_provider,
            self.claim_requests,
            self.actions,
            self.gauge_reductions,
            self.numerical,
            self.initialization,
            self.stopping,
            self.determinism,
            source_authority,
        )
        .expect("rebuilt execution plan")
    }
}

#[derive(Clone)]
struct AssessmentParts {
    header: ArtifactHeader,
    claims: Vec<TypedIdentifiabilityClaim>,
    evidence: Vec<(ClaimId, ClaimAssessment)>,
    source_authority: SourceResolutionSet,
}

impl AssessmentParts {
    fn from_assessment(assessment: &IdentifiabilityAssessment) -> Self {
        Self {
            header: assessment.header().clone(),
            claims: assessment.claims().values().cloned().collect(),
            evidence: assessment
                .evidence()
                .iter()
                .map(|(id, value)| (id.clone(), value.clone()))
                .collect(),
            source_authority: assessment.source_authority().clone(),
        }
    }

    fn build(
        self,
        problem: &AdmittedIdentifiabilityProblem,
        execution: &IdentifiabilityExecutionPlan,
    ) -> IdentifiabilityAssessment {
        IdentifiabilityAssessment::try_new(
            self.header,
            problem,
            execution,
            self.claims,
            self.evidence,
            self.source_authority,
        )
        .expect("rebuilt assessment")
    }
}

fn assessment(
    problem: &AdmittedIdentifiabilityProblem,
    execution: &IdentifiabilityExecutionPlan,
    receipt_label: &str,
) -> IdentifiabilityAssessment {
    assessment_result(problem, execution, receipt_label).expect("assessment fixture")
}

fn assessment_result(
    problem: &AdmittedIdentifiabilityProblem,
    execution: &IdentifiabilityExecutionPlan,
    receipt_label: &str,
) -> Result<IdentifiabilityAssessment, IdentifiabilityError> {
    assessment_result_with_claim_source_authority(
        problem,
        execution,
        receipt_label,
        AuthorityDisposition::ContentVerified,
    )
}

fn assessment_result_with_claim_source_authority(
    problem: &AdmittedIdentifiabilityProblem,
    execution: &IdentifiabilityExecutionPlan,
    receipt_label: &str,
    claim_source_authority: AuthorityDisposition,
) -> Result<IdentifiabilityAssessment, IdentifiabilityError> {
    let claim = default_claim(problem);
    let claim_id = claim.id().clone();
    let authority_key = match claim.fiber() {
        FiberStructure::Stratified { strata } => strata.key().clone(),
        _ => match claim.quantifier() {
            ClaimQuantifier::ForAll { domain } => domain.key().clone(),
            _ => unreachable!("default claim uses universal quantification"),
        },
    };
    let request =
        execution
            .claim_requests()
            .get(&claim_id)
            .ok_or(IdentifiabilityError::SourceMismatch {
                field: "assessment/execution exact claim preregistration",
            })?;
    let method = execution.analyzer().clone();
    let receipt = source(
        receipt_label,
        SourceKind::EvidenceReceipt,
        hash(receipt_label),
    );
    let mut assessment_sources = claim_sources(&claim);
    assessment_sources.extend([
        method.clone(),
        receipt.clone(),
        request.error_policy().metric().clone(),
        request.error_policy().nondimensionalization().clone(),
    ]);
    let source_authority = resolve_owned_sources(
        assessment_sources,
        false,
        Some((&authority_key, claim_source_authority)),
    );
    IdentifiabilityAssessment::try_new(
        header("assessment-1", "identifiability.assess"),
        problem,
        execution,
        vec![claim],
        vec![(
            claim_id,
            ClaimAssessment::ClaimedEstablished {
                method,
                receipt,
                metric: request.error_policy().metric().clone(),
                nondimensionalization: request.error_policy().nondimensionalization().clone(),
                certified_error_bound: 5.0e-9,
                gauge_resolutions: BTreeMap::new(),
            },
        )],
        source_authority,
    )
}

fn assessment_with_claim(
    problem: &AdmittedIdentifiabilityProblem,
    execution: &IdentifiabilityExecutionPlan,
    claim: TypedIdentifiabilityClaim,
    mut claim_sources: Vec<SourceRef>,
) -> Result<IdentifiabilityAssessment, IdentifiabilityError> {
    let claim_id = claim.id().clone();
    let request =
        execution
            .claim_requests()
            .get(&claim_id)
            .ok_or(IdentifiabilityError::SourceMismatch {
                field: "assessment/execution exact claim preregistration",
            })?;
    let method = execution.analyzer().clone();
    let receipt = source(
        "custom-claim-receipt",
        SourceKind::EvidenceReceipt,
        hash("custom-claim-receipt"),
    );
    claim_sources.extend([
        method.clone(),
        receipt.clone(),
        request.error_policy().metric().clone(),
        request.error_policy().nondimensionalization().clone(),
    ]);
    IdentifiabilityAssessment::try_new(
        header("assessment-custom-claim", "identifiability.assess"),
        problem,
        execution,
        vec![claim],
        vec![(
            claim_id,
            ClaimAssessment::ClaimedEstablished {
                method,
                receipt,
                metric: request.error_policy().metric().clone(),
                nondimensionalization: request.error_policy().nondimensionalization().clone(),
                certified_error_bound: 5.0e-9,
                gauge_resolutions: BTreeMap::new(),
            },
        )],
        resolve_owned_sources(claim_sources, false, None),
    )
}

fn two_claims() -> Vec<TypedIdentifiabilityClaim> {
    let domain_left = source(
        "claim-domain-left",
        SourceKind::QuantifierDomain,
        hash("claim-domain-left"),
    );
    let domain_right = source(
        "claim-domain-right",
        SourceKind::QuantifierDomain,
        hash("claim-domain-right"),
    );
    let left_id = ClaimId::try_new("claim-left").expect("left claim id");
    let right_id = ClaimId::try_new("claim-right").expect("right claim id");
    vec![
        TypedIdentifiabilityClaim::new(
            left_id.clone(),
            InformationRegime::StructuralExactModel,
            IdentifiabilityExtent::Global,
            FiberStructure::Unique,
            ClaimQuantifier::ForAll {
                domain: domain_left.clone(),
            },
            ScalarDomain::Real,
            ClaimSubject::Parameter(role("yield_stress")),
            ClaimScope::WholeCampaign,
        ),
        TypedIdentifiabilityClaim::new(
            right_id.clone(),
            InformationRegime::StructuralExactModel,
            IdentifiabilityExtent::Global,
            FiberStructure::Unique,
            ClaimQuantifier::ForAll {
                domain: domain_right.clone(),
            },
            ScalarDomain::Real,
            ClaimSubject::Parameter(role("hardening_modulus")),
            ClaimScope::WholeCampaign,
        ),
    ]
}

fn two_claim_assessment(
    problem: &AdmittedIdentifiabilityProblem,
    execution: &IdentifiabilityExecutionPlan,
) -> IdentifiabilityAssessment {
    let method = execution.analyzer().clone();
    let claims = two_claims();
    let receipt_left = source(
        "claim-receipt-left",
        SourceKind::EvidenceReceipt,
        hash("claim-receipt-left"),
    );
    let receipt_right = source(
        "claim-receipt-right",
        SourceKind::EvidenceReceipt,
        hash("claim-receipt-right"),
    );
    let left_id = ClaimId::try_new("claim-left").expect("left claim id");
    let right_id = ClaimId::try_new("claim-right").expect("right claim id");
    let left_policy = execution
        .claim_requests()
        .get(&left_id)
        .expect("left claim preregistered")
        .error_policy();
    let right_policy = execution
        .claim_requests()
        .get(&right_id)
        .expect("right claim preregistered")
        .error_policy();
    let evidence = vec![
        (
            left_id,
            ClaimAssessment::ClaimedEstablished {
                method: method.clone(),
                receipt: receipt_left.clone(),
                metric: left_policy.metric().clone(),
                nondimensionalization: left_policy.nondimensionalization().clone(),
                certified_error_bound: 5.0e-9,
                gauge_resolutions: BTreeMap::new(),
            },
        ),
        (
            right_id,
            ClaimAssessment::ClaimedEstablished {
                method: method.clone(),
                receipt: receipt_right.clone(),
                metric: right_policy.metric().clone(),
                nondimensionalization: right_policy.nondimensionalization().clone(),
                certified_error_bound: 5.0e-9,
                gauge_resolutions: BTreeMap::new(),
            },
        ),
    ];
    let mut assessment_sources = claims.iter().flat_map(claim_sources).collect::<Vec<_>>();
    assessment_sources.extend([
        method.clone(),
        receipt_left.clone(),
        receipt_right.clone(),
        left_policy.metric().clone(),
        left_policy.nondimensionalization().clone(),
        right_policy.metric().clone(),
        right_policy.nondimensionalization().clone(),
    ]);
    let source_authority = resolve_owned_sources(assessment_sources, false, None);
    IdentifiabilityAssessment::try_new(
        header("assessment-2", "identifiability.assess"),
        problem,
        execution,
        claims,
        evidence,
        source_authority,
    )
    .expect("two-claim assessment")
}

#[test]
fn assessment_requires_the_exact_preregistered_proposition() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let plan = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    assert_eq!(
        plan.requested_axes(),
        BTreeSet::from([RequestedClaimAxis::Structural, RequestedClaimAxis::Global,]),
        "coarse axes remain a derived planner projection",
    );
    let baseline = default_claim(&problem);
    let substitutions = [
        (
            "quantifier",
            TypedIdentifiabilityClaim::new(
                baseline.id().clone(),
                baseline.information().clone(),
                baseline.extent(),
                baseline.fiber().clone(),
                ClaimQuantifier::ForAll {
                    domain: source(
                        "substituted-quantifier-domain",
                        SourceKind::QuantifierDomain,
                        hash("substituted-quantifier-domain"),
                    ),
                },
                baseline.scalar_domain().clone(),
                baseline.subject().clone(),
                baseline.scope().clone(),
            ),
        ),
        (
            "subject",
            TypedIdentifiabilityClaim::new(
                baseline.id().clone(),
                baseline.information().clone(),
                baseline.extent(),
                baseline.fiber().clone(),
                baseline.quantifier().clone(),
                baseline.scalar_domain().clone(),
                ClaimSubject::Parameter(role("hardening_modulus")),
                baseline.scope().clone(),
            ),
        ),
        (
            "scope",
            TypedIdentifiabilityClaim::new(
                baseline.id().clone(),
                baseline.information().clone(),
                baseline.extent(),
                baseline.fiber().clone(),
                baseline.quantifier().clone(),
                baseline.scalar_domain().clone(),
                baseline.subject().clone(),
                ClaimScope::Cases(BTreeSet::from([case_id("a")])),
            ),
        ),
        (
            "fiber",
            TypedIdentifiabilityClaim::new(
                baseline.id().clone(),
                baseline.information().clone(),
                baseline.extent(),
                FiberStructure::FiniteToOne {
                    maximum_cardinality: Some(FiberCardinalityBound::UniformU64(2)),
                },
                baseline.quantifier().clone(),
                baseline.scalar_domain().clone(),
                baseline.subject().clone(),
                baseline.scope().clone(),
            ),
        ),
    ];
    for (field, substituted) in substitutions {
        let substituted_plan =
            execution_for_claim(&problem, false, 17, 1.0e-10, false, substituted.clone())
                .unwrap_or_else(|error| panic!("valid {field} substitution refused: {error}"));
        assert_eq!(
            substituted_plan.requested_axes(),
            plan.requested_axes(),
            "{field} substitution must retain identical coarse planner axes",
        );
        let substituted_sources = claim_sources(&substituted);
        assert!(matches!(
            assessment_with_claim(&problem, &plan, substituted, substituted_sources),
            Err(IdentifiabilityError::SourceMismatch {
                field: "assessment/execution exact claim preregistration",
            })
        ));
    }
    assessment_result(&problem, &plan, "exact-proposition-receipt")
        .expect("the exact preregistered proposition admits");
    log(
        "exact-claim-preregistration",
        "pass",
        "matching coarse axes cannot authorize a post-hoc proposition substitution",
    );
}

#[test]
fn claimed_evidence_is_bound_to_metric_scale_and_preregistered_error_ceiling() {
    let metric = source(
        "claim-error-metric",
        SourceKind::DimensionlessErrorMetric,
        hash("claim-error-metric"),
    );
    let nondimensionalization = source(
        "claim-nondimensionalization",
        SourceKind::Nondimensionalization,
        hash("claim-nondimensionalization"),
    );
    assert!(matches!(
        DimensionlessErrorPolicy::try_new(
            source(
                "dimensional-error-metric",
                SourceKind::Assumption,
                hash("dimensional-error-metric"),
            ),
            nondimensionalization.clone(),
            1.0e-8,
        ),
        Err(IdentifiabilityError::InvalidText {
            field: "claim error metric",
            ..
        })
    ));
    assert!(matches!(
        DimensionlessErrorPolicy::try_new(
            metric,
            source(
                "implicit-claim-scale",
                SourceKind::Assumption,
                hash("implicit-claim-scale"),
            ),
            1.0e-8,
        ),
        Err(IdentifiabilityError::InvalidText {
            field: "claim nondimensionalization",
            ..
        })
    ));
    assert!(matches!(
        IdentifiabilityNumericalPolicy::try_new(
            1.0e-10,
            0.0,
            1.0e12,
            ArithmeticPolicy::CertifiedInterval,
            source(
                "dimensional-rank-policy",
                SourceKind::Assumption,
                hash("dimensional-rank-policy"),
            ),
        ),
        Err(IdentifiabilityError::InvalidText {
            field: "numerical nondimensionalization",
            ..
        })
    ));
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let execution = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    let baseline = assessment(&problem, &execution, "metric-policy-receipt");
    let claim_id = baseline
        .claims()
        .keys()
        .next()
        .expect("baseline claim")
        .clone();
    let (method, receipt, metric, nondimensionalization) = match &baseline.evidence()[&claim_id] {
        ClaimAssessment::ClaimedEstablished {
            method,
            receipt,
            metric,
            nondimensionalization,
            ..
        } => (
            method.clone(),
            receipt.clone(),
            metric.clone(),
            nondimensionalization.clone(),
        ),
        _ => panic!("baseline evidence unexpectedly changed variant"),
    };
    let assess = |evidence| {
        IdentifiabilityAssessment::try_new(
            baseline.header().clone(),
            &problem,
            &execution,
            baseline.claims().values().cloned().collect(),
            vec![(claim_id.clone(), evidence)],
            baseline.source_authority().clone(),
        )
    };

    assert!(matches!(
        assess(ClaimAssessment::ClaimedEstablished {
            method: method.clone(),
            receipt: receipt.clone(),
            metric: source(
                "post-hoc-error-metric",
                SourceKind::DimensionlessErrorMetric,
                hash("post-hoc-error-metric"),
            ),
            nondimensionalization: nondimensionalization.clone(),
            certified_error_bound: 5.0e-9,
            gauge_resolutions: BTreeMap::new(),
        }),
        Err(IdentifiabilityError::SourceMismatch {
            field: "claim evidence/error policy",
        })
    ));
    assess(ClaimAssessment::ClaimedRefuted {
        method: method.clone(),
        receipt: receipt.clone(),
        metric: metric.clone(),
        nondimensionalization: nondimensionalization.clone(),
        certified_error_bound: 5.0e-9,
    })
    .expect("refuting evidence uses the same preregistered dimensionless policy");
    assert!(matches!(
        assess(ClaimAssessment::ClaimedEstablished {
            method: method.clone(),
            receipt: receipt.clone(),
            metric: metric.clone(),
            nondimensionalization: source(
                "post-hoc-nondimensionalization",
                SourceKind::Nondimensionalization,
                hash("post-hoc-nondimensionalization"),
            ),
            certified_error_bound: 5.0e-9,
            gauge_resolutions: BTreeMap::new(),
        }),
        Err(IdentifiabilityError::SourceMismatch {
            field: "claim evidence/error policy",
        })
    ));
    assert!(matches!(
        assess(ClaimAssessment::ClaimedEstablished {
            method,
            receipt,
            metric,
            nondimensionalization,
            certified_error_bound: 2.0e-8,
            gauge_resolutions: BTreeMap::new(),
        }),
        Err(IdentifiabilityError::InvalidNumeric {
            field: "certified claim error",
            ..
        })
    ));
    log(
        "claim-error-policy-binding",
        "pass",
        "claimed evidence cannot substitute its dimensionless metric, scaling policy, or preregistered error ceiling",
    );
}

#[test]
fn claim_quantifiers_require_semantically_typed_sources() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let claim_for = |quantifier| {
        TypedIdentifiabilityClaim::new(
            ClaimId::try_new("quantifier-kind-contract").expect("claim id"),
            InformationRegime::StructuralExactModel,
            IdentifiabilityExtent::Global,
            FiberStructure::Unique,
            quantifier,
            ScalarDomain::Real,
            ClaimSubject::Parameter(role("yield_stress")),
            ClaimScope::WholeCampaign,
        )
    };

    for quantifier in [
        ClaimQuantifier::AtRealization {
            realization: source(
                "typed-realization",
                SourceKind::QuantifierRealization,
                hash("typed-realization"),
            ),
        },
        ClaimQuantifier::AlmostEverywhere {
            measure: source(
                "reference-measure",
                SourceKind::ReferenceMeasure,
                hash("reference-measure"),
            ),
        },
        ClaimQuantifier::AlmostEverywhere {
            measure: source(
                "probability-reference-measure",
                SourceKind::ProbabilityMeasure,
                hash("probability-reference-measure"),
            ),
        },
        ClaimQuantifier::ForAll {
            domain: source(
                "universal-domain",
                SourceKind::QuantifierDomain,
                hash("universal-domain"),
            ),
        },
        ClaimQuantifier::ProbabilityAtLeast {
            probability: 0.95,
            measure: source(
                "probability-measure",
                SourceKind::ProbabilityMeasure,
                hash("probability-measure"),
            ),
        },
    ] {
        let claim = claim_for(quantifier);
        let execution = execution_for_claim(&problem, false, 17, 1.0e-10, false, claim.clone())
            .expect("typed quantifier preregisters");
        assessment_with_claim(&problem, &execution, claim.clone(), claim_sources(&claim))
            .expect("typed quantifier assesses");
    }

    for (field, quantifier) in [
        (
            "claim realization source",
            ClaimQuantifier::AtRealization {
                realization: source(
                    "bad-realization-prior",
                    SourceKind::ProbabilityMeasure,
                    hash("bad-realization-prior"),
                ),
            },
        ),
        (
            "claim almost-everywhere measure",
            ClaimQuantifier::AlmostEverywhere {
                measure: source(
                    "bad-ae-analyzer",
                    SourceKind::QuantifierRealization,
                    hash("bad-ae-analyzer"),
                ),
            },
        ),
        (
            "claim universal domain",
            ClaimQuantifier::ForAll {
                domain: source(
                    "bad-domain-build",
                    SourceKind::ReferenceMeasure,
                    hash("bad-domain-build"),
                ),
            },
        ),
        (
            "claim probability measure",
            ClaimQuantifier::ProbabilityAtLeast {
                probability: 0.95,
                measure: source(
                    "bad-probability-manifold",
                    SourceKind::QuantifierDomain,
                    hash("bad-probability-manifold"),
                ),
            },
        ),
    ] {
        let claim = claim_for(quantifier);
        assert!(matches!(
            execution_for_claim(&problem, false, 17, 1.0e-10, false, claim),
            Err(IdentifiabilityError::InvalidText { field: actual, .. }) if actual == field
        ));
    }
    log(
        "claim-quantifier-source-kinds",
        "pass",
        "realizations, universal domains, almost-everywhere measures, and probability measures enforce distinct source-kind contracts",
    );
}

#[test]
fn continuous_uniform_prior_refuses_atomic_or_out_of_domain_support() {
    let physical = ParameterDomain::try_new(0.0, 2.0).expect("physical parameter domain");
    let parameter_with_prior = |name: &str, prior_domain: ParameterDomain| {
        StudyParameter::try_new(
            role(name),
            QuantitySpec::dimensional(STRESS),
            physical,
            ParameterPurpose::Estimand,
            ParameterTreatment::Estimated,
            ParameterOwnerBinding::ConstitutiveModel,
            ParameterScopeBinding::Global,
            PriorPolicy::Distribution(ParameterPrior::Uniform {
                domain: prior_domain,
                version: 1,
            }),
            InfluenceCoverage::Declared,
        )
    };

    let singleton = ParameterDomain::try_new(1.0, 1.0).expect("singleton support");
    assert!(matches!(
        parameter_with_prior("singleton-uniform", singleton),
        Err(IdentifiabilityError::InvalidNumeric {
            field: "uniform prior support",
            ..
        })
    ));
    let signed_zero = ParameterDomain::try_new(-0.0, 0.0).expect("signed-zero support");
    assert!(matches!(
        parameter_with_prior("signed-zero-uniform", signed_zero),
        Err(IdentifiabilityError::InvalidNumeric {
            field: "uniform prior support",
            ..
        })
    ));
    let outside = ParameterDomain::try_new(-1.0, 1.0).expect("out-of-domain support");
    assert!(matches!(
        parameter_with_prior("outside-uniform", outside),
        Err(IdentifiabilityError::InvalidNumeric {
            field: "uniform prior support",
            ..
        })
    ));
    let positive = ParameterDomain::try_new(0.5, 1.5).expect("positive-width support");
    parameter_with_prior("positive-width-uniform", positive)
        .expect("positive-width contained uniform support");
    log(
        "continuous-uniform-support",
        "pass",
        "continuous Uniform priors require positive-width contained support; atomic mass awaits explicit discrete/Dirac semantics",
    );
}

#[test]
fn posterior_claims_bind_one_exact_joint_prior_measure() {
    let fixture = retrospective_origin_fixture_with_options(
        ExperimentOrigin::Physical {
            apparatus_id: artifact("posterior-apparatus"),
            facility_id: artifact("posterior-facility"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::Uncharacterized,
        ProblemOptions {
            joint_prior_choice: 1,
            ..ProblemOptions::default()
        },
    );
    let problem = admit_retrospective_origin_fixture(&fixture).expect("posterior problem admits");
    let posterior_claim = |joint_prior| {
        TypedIdentifiabilityClaim::new(
            ClaimId::try_new("posterior-joint-prior-contract").expect("claim id"),
            InformationRegime::PosteriorUnderDeclaredPrior { joint_prior },
            IdentifiabilityExtent::Global,
            FiberStructure::Unique,
            ClaimQuantifier::ForAll {
                domain: source(
                    "posterior-domain",
                    SourceKind::QuantifierDomain,
                    hash("posterior-domain"),
                ),
            },
            ScalarDomain::Real,
            ClaimSubject::Parameter(role("yield_stress")),
            ClaimScope::WholeCampaign,
        )
    };
    let claim =
        posterior_claim(problem.document().sources()[&source_key("joint-prior-measure-a")].clone());
    let execution = execution_for_claim(&problem, false, 17, 1.0e-10, false, claim.clone())
        .expect("posterior claim preregisters");
    let assessment = IdentifiabilityAssessment::try_new(
        header("posterior-assessment", "identifiability.assess"),
        &problem,
        &execution,
        vec![claim.clone()],
        vec![(
            claim.id().clone(),
            ClaimAssessment::ClaimedInconclusive {
                method: None,
                receipt: None,
                reason: "the exact joint prior is bound, but no decisive posterior theorem receipt is claimed"
                    .to_string(),
            },
        )],
        resolve_owned_sources(claim_sources(&claim), false, None),
    )
    .expect("posterior claim records an honest inconclusive assessment");
    assert!(
        assessment
            .source_authority()
            .entries()
            .get(&source_key("joint-prior-measure-a"))
            .is_some(),
        "the joint prior must remain in the assessment authority envelope",
    );

    let wrong_kind = posterior_claim(source(
        "bare-prior-family",
        SourceKind::Prior,
        hash("bare-prior-family"),
    ));
    assert!(matches!(
        execution_for_claim(&problem, false, 17, 1.0e-10, false, wrong_kind),
        Err(IdentifiabilityError::InvalidText {
            field: "claim joint-prior measure",
            ..
        })
    ));
    let mismatched_measure = posterior_claim(source(
        "joint-prior-measure-b",
        SourceKind::ProbabilityMeasure,
        hash("joint-prior-measure-b"),
    ));
    assert!(matches!(
        execution_for_claim(&problem, false, 17, 1.0e-10, false, mismatched_measure,),
        Err(IdentifiabilityError::SourceMismatch {
            field: "claim joint prior/problem joint prior",
        })
    ));
    log(
        "posterior-joint-prior",
        "pass",
        "posterior semantics carry one exact authorized probability measure rather than an implicit product-of-marginals assumption",
    );
}

#[test]
fn claim_product_rejects_cross_axis_semantic_contradictions() {
    let domain = source(
        "compatibility-domain",
        SourceKind::QuantifierDomain,
        hash("compatibility-domain"),
    );
    let make_claim = |information, fiber, scope| {
        TypedIdentifiabilityClaim::new(
            ClaimId::try_new("compatibility-claim").expect("claim id"),
            information,
            IdentifiabilityExtent::Global,
            fiber,
            ClaimQuantifier::ForAll {
                domain: domain.clone(),
            },
            ScalarDomain::Real,
            ClaimSubject::Parameter(role("yield_stress")),
            scope,
        )
    };

    let prospective = admit_fixture(problem_fixture(ProblemOptions::default()));
    let claim = make_claim(
        InformationRegime::NoisyFiniteData,
        FiberStructure::Unique,
        ClaimScope::WholeCampaign,
    );
    assert!(matches!(
        execution_for_claim(&prospective, false, 17, 1.0e-10, false, claim,),
        Err(IdentifiabilityError::InvalidText {
            field: "finite-data claim scope",
            ..
        })
    ));

    let missing_prior_fixture = retrospective_origin_fixture_with_options(
        ExperimentOrigin::Physical {
            apparatus_id: artifact("missing-prior-apparatus"),
            facility_id: artifact("missing-prior-facility"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::Uncharacterized,
        ProblemOptions {
            yield_prior_absent: true,
            joint_prior_choice: 1,
            ..ProblemOptions::default()
        },
    );
    let missing_prior = admit_retrospective_origin_fixture(&missing_prior_fixture)
        .expect("retrospective missing-parameter-prior problem");
    let joint_prior =
        missing_prior.document().sources()[&source_key("joint-prior-measure-a")].clone();
    let claim = make_claim(
        InformationRegime::PosteriorUnderDeclaredPrior { joint_prior },
        FiberStructure::Unique,
        ClaimScope::WholeCampaign,
    );
    assert!(matches!(
        execution_for_claim(&missing_prior, false, 17, 1.0e-10, false, claim,),
        Err(IdentifiabilityError::InvalidText {
            field: "posterior claim prior",
            ..
        })
    ));

    let scoped = admit_fixture(problem_fixture(ProblemOptions {
        yield_case_a_only: true,
        ..ProblemOptions::default()
    }));
    let claim = make_claim(
        InformationRegime::StructuralExactModel,
        FiberStructure::Unique,
        ClaimScope::Cases(BTreeSet::from([case_id("b")])),
    );
    assert!(matches!(
        execution_for_claim(&scoped, false, 17, 1.0e-10, false, claim,),
        Err(IdentifiabilityError::InvalidText {
            field: "claim parameter/case scope",
            ..
        })
    ));

    let claim = make_claim(
        InformationRegime::StructuralExactModel,
        FiberStructure::FiniteToOne {
            maximum_cardinality: Some(FiberCardinalityBound::UniformU64(1)),
        },
        ClaimScope::WholeCampaign,
    );
    assert!(matches!(
        execution_for_claim(&prospective, false, 17, 1.0e-10, false, claim,),
        Err(IdentifiabilityError::InvalidNumeric {
            field: "finite-to-one cardinality",
            ..
        })
    ));
    log(
        "claim-product-compatibility",
        "pass",
        "data regime, prior policy, parameter scope, and fiber cardinality fail closed as one product",
    );
}

#[test]
fn influence_endpoints_must_lie_inside_parameter_applicability() {
    let result = problem_fixture(ProblemOptions {
        yield_case_a_only: true,
        yield_influence_case_b: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::InvalidText {
            field: "influence parameter/case scope",
            ..
        })
    ));
    log(
        "influence-parameter-case-scope",
        "pass",
        "a declared influence endpoint cannot escape the exact applicability of its physical parameter",
    );
}

#[test]
fn influence_claim_scope_closes_over_the_transitive_composite_dag() {
    let problem = admit_fixture(problem_fixture(ProblemOptions {
        composite_influence_chain: true,
        ..ProblemOptions::default()
    }));
    let case_a_only = structural_claim(
        &problem,
        "transitive-influence-scope",
        FiberStructure::Unique,
        ClaimSubject::Influence(InfluenceId::try_new("composite-top").expect("top influence id")),
        ClaimScope::Cases(BTreeSet::from([case_id("a")])),
    );
    assert!(matches!(
        execution_for_claim(&problem, false, 17, 1.0e-10, false, case_a_only),
        Err(IdentifiabilityError::InvalidText {
            field: "claim influence/case scope",
            ..
        })
    ));

    let closed_scope = structural_claim(
        &problem,
        "transitive-influence-scope",
        FiberStructure::Unique,
        ClaimSubject::Influence(InfluenceId::try_new("composite-top").expect("top influence id")),
        ClaimScope::Cases(BTreeSet::from([case_id("a"), case_id("b")])),
    );
    let execution = execution_for_claim(&problem, false, 17, 1.0e-10, false, closed_scope.clone())
        .expect("claim scope containing every transitive endpoint");
    assessment_with_claim(
        &problem,
        &execution,
        closed_scope.clone(),
        claim_sources(&closed_scope),
    )
    .expect("transitively closed influence claim assesses");
    log(
        "influence-transitive-claim-scope",
        "pass",
        "ClaimSubject::Influence closes over every recursive composite input endpoint and parameter",
    );
}

#[test]
fn source_bound_derived_complex_and_local_set_valued_claims_round_trip() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let domain = source(
        "derived-claim-domain",
        SourceKind::QuantifierDomain,
        hash("derived-claim-domain"),
    );
    let definition = source(
        "derived-functional",
        SourceKind::DerivedFunctional,
        hash("derived-functional"),
    );
    let extension = source(
        "complex-extension",
        SourceKind::AlgebraicExtension,
        hash("complex-extension"),
    );
    let claim = TypedIdentifiabilityClaim::new(
        ClaimId::try_new("derived-complex-positive-dimensional-fiber").expect("claim id"),
        InformationRegime::StructuralExactModel,
        IdentifiabilityExtent::Local,
        FiberStructure::PositiveDimensional {
            lower_bound: FiberDimensionLowerBound::Finite {
                minimum_dimension: 1,
            },
        },
        ClaimQuantifier::ForAll {
            domain: domain.clone(),
        },
        ScalarDomain::Complex {
            extension: extension.clone(),
        },
        ClaimSubject::DerivedFunctional {
            definition: definition.clone(),
            parameters: BTreeSet::from([role("yield_stress")]),
        },
        ClaimScope::Cases(BTreeSet::from([case_id("a")])),
    );
    let execution = execution_for_claim(&problem, false, 17, 1.0e-10, false, claim.clone())
        .expect("source-bound claim execution");
    let assessment = assessment_with_claim(
        &problem,
        &execution,
        claim,
        vec![domain, definition, extension],
    )
    .expect("source-bound local positive-dimensional claim");
    let bytes = assessment.canonical_bytes().expect("assessment bytes");
    let replay = IdentifiabilityAssessment::from_canonical_bytes(
        &bytes,
        &problem,
        &execution,
        assessment.source_authority(),
    )
    .expect("assessment replay");
    assert_eq!(assessment, replay);
    assert_eq!(
        assessment.id().expect("assessment id"),
        replay.id().expect("replay id")
    );
    log(
        "claim-taxonomy-round-trip",
        "pass",
        "local positive-dimensional derived functionals over a source-bound complex extension are canonical",
    );
}

#[test]
fn problem_roundtrip_is_canonical_but_remains_unresolved() {
    let document = problem_fixture(ProblemOptions::default())
        .document
        .expect("valid problem");
    let bytes = document.canonical_bytes().expect("problem encodes");
    let decoded = IdentifiabilityProblemDocument::from_canonical_bytes(&bytes)
        .expect("unresolved problem decodes");
    assert_eq!(decoded, document);
    assert_eq!(decoded.cases().len(), 2);
    log(
        "problem-roundtrip-unresolved",
        "pass",
        "decode returned only an unresolved multi-case document",
    );
}

#[test]
fn case_and_registry_input_order_are_nonsemantic() {
    let left = problem_fixture(ProblemOptions::default())
        .document
        .expect("left problem");
    let right = problem_fixture(ProblemOptions {
        reverse_cases: true,
        ..ProblemOptions::default()
    })
    .document
    .expect("right problem");
    assert_eq!(
        left.canonical_bytes().expect("left bytes"),
        right.canonical_bytes().expect("right bytes")
    );
    log(
        "case-order",
        "pass",
        "canonical maps erase caller insertion order without erasing cases",
    );
}

#[test]
fn source_resolution_input_order_is_nonsemantic() {
    let forward = admit_fixture_with_authority_and_order(
        problem_fixture(ProblemOptions::default()),
        false,
        false,
    );
    let reverse = admit_fixture_with_authority_and_order(
        problem_fixture(ProblemOptions::default()),
        false,
        true,
    );
    assert_eq!(forward.id(), reverse.id());
    assert_eq!(forward.source_admission_id(), reverse.source_admission_id());
    assert_eq!(
        forward
            .source_admission_canonical_bytes()
            .expect("forward authority bytes"),
        reverse
            .source_admission_canonical_bytes()
            .expect("reverse authority bytes"),
    );
    log(
        "source-resolution-order",
        "pass",
        "source authority canonicalization erases caller insertion order",
    );
}

#[test]
fn multi_case_campaign_qualifies_local_channel_names() {
    let problem = problem_fixture(ProblemOptions::default())
        .document
        .expect("problem");
    assert!(
        problem
            .cases()
            .contains_key(&CaseId::try_new("a").expect("case"))
    );
    assert!(
        problem
            .cases()
            .contains_key(&CaseId::try_new("b").expect("case"))
    );
    assert_ne!(
        ObservationKey::new(case_id("a"), channel("stress")),
        ObservationKey::new(case_id("b"), channel("stress"))
    );
    log(
        "composite-observation-key",
        "pass",
        "case qualification prevents cross-protocol channel aliasing",
    );
}

#[test]
fn lot_field_and_hierarchical_scopes_retain_explicit_case_support() {
    let supported_cases = BTreeSet::from([case_id("a")]);
    for scope in [
        ParameterScopeBinding::MaterialLot {
            lot: artifact("lot-a"),
            cases: supported_cases.clone(),
        },
        ParameterScopeBinding::Field {
            support: source_key("field-support-a"),
            cases: supported_cases.clone(),
        },
        ParameterScopeBinding::Hierarchical {
            population: artifact("population-a"),
            level: 2,
            hierarchy: source_key("hierarchy-a"),
            cases: supported_cases.clone(),
        },
    ] {
        let parameter = parameter(
            "yield_stress",
            ParameterTreatment::Estimated,
            InfluenceCoverage::Declared,
            1,
            None,
            false,
            scope.clone(),
        );
        assert_eq!(parameter.scope(), &scope);
        let cases = match parameter.scope() {
            ParameterScopeBinding::MaterialLot { cases, .. }
            | ParameterScopeBinding::Field { cases, .. }
            | ParameterScopeBinding::Hierarchical { cases, .. } => cases,
            _ => panic!("fixture scope unexpectedly changed variant"),
        };
        assert_eq!(cases, &supported_cases);
    }
    log(
        "explicit-parameter-case-support",
        "pass",
        "lot, field, and hierarchical parameters retain exact case applicability instead of silently becoming campaign-global",
    );
}

#[test]
fn dangling_composite_observation_endpoint_refuses() {
    let result = problem_fixture(ProblemOptions {
        bad_observation_endpoint: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::UnknownReference {
            field: "composite observation key",
            ..
        })
    ));
    log(
        "dangling-observation",
        "pass",
        "unknown case/channel refused",
    );
}

#[test]
fn disconnected_free_parameter_refuses_without_false_theorem() {
    let result = problem_fixture(ProblemOptions {
        missing_hardening_influence: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::DisconnectedEstimatedParameter { .. })
    ));
    log(
        "disconnected-parameter",
        "pass",
        "schema connectivity refused without claiming nonzero sensitivity",
    );
}

#[test]
fn dangling_source_key_refuses_before_authority_admission() {
    let result = problem_fixture(ProblemOptions {
        dangling_operator: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::UnknownReference {
            field: "observation operator",
            ..
        })
    ));
    log("dangling-source", "pass", "source registry closure refused");
}

#[test]
fn case_physics_sources_require_exact_role_hash_domain_and_version() {
    for mutation in 1..=5 {
        let error = problem_fixture(ProblemOptions {
            case_physics_mutation: mutation,
            ..ProblemOptions::default()
        })
        .document
        .expect_err("mutated case-physics source must refuse");
        assert!(
            matches!(
                error,
                IdentifiabilityError::UnknownReference {
                    field: "case specimen-geometry source",
                    ..
                } | IdentifiabilityError::SourceMismatch {
                    field: "case specimen-geometry source"
                }
            ),
            "mutation {mutation} returned the wrong diagnostic: {error}"
        );
    }
    let admitted = admit_fixture(problem_fixture(ProblemOptions::default()));
    let geometry = admitted
        .document()
        .sources()
        .get(&source_key("geometry-a"))
        .expect("admitted geometry source");
    let resolution = admitted
        .source_resolutions()
        .get(geometry.key())
        .expect("admitted geometry resolution");
    assert!(matches!(
        resolution.verification(),
        SourceVerification::HashPreimage { byte_len }
            if *byte_len == u64::try_from("geometry-a".len()).expect("fixture length")
    ));
    log(
        "case-physics-source-closure",
        "pass",
        "every embedded physics digest is role/domain/version bound and hash-preimage verified before ProblemId",
    );
}

#[test]
fn observation_contract_mismatches_report_the_exact_field() {
    for (mutation, expected) in [
        (1, "case observation protocol version"),
        (2, "case observation refinement version"),
    ] {
        assert!(matches!(
            problem_fixture(ProblemOptions {
                observation_contract_mutation: mutation,
                ..ProblemOptions::default()
            })
            .document,
            Err(IdentifiabilityError::VersionMismatch { field, .. }) if field == expected
        ));
    }
    assert!(matches!(
        problem_fixture(ProblemOptions {
            observation_contract_mutation: 3,
            ..ProblemOptions::default()
        })
        .document,
        Err(IdentifiabilityError::InvalidText {
            field: "case observation protocol clock",
            detail,
        }) if detail.contains("wrong-clock-a") && detail.contains("clock-a")
    ));
    log(
        "observation-contract-diagnostics",
        "pass",
        "protocol, refinement, and clock mismatches report distinct fields and exact clock identities",
    );
}

#[test]
fn blind_falsification_cannot_bypass_release_with_prospective_data() {
    assert!(matches!(
        problem_fixture(ProblemOptions {
            blind_prospective_case: true,
            ..ProblemOptions::default()
        })
        .document,
        Err(IdentifiabilityError::InvalidText {
            field: "blind-falsification case data",
            ..
        })
    ));
    log(
        "blind-prospective-bypass",
        "pass",
        "blind falsification requires retrospective blind rows before release authority can be considered",
    );
}

#[test]
fn derived_parameter_cycles_refuse() {
    let result = problem_fixture(ProblemOptions {
        derived_cycle: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::InvalidNumeric {
            field: "derived parameter graph",
            ..
        })
    ));
    log(
        "derived-cycle",
        "pass",
        "derived parameter DAG is fail-closed",
    );
}

#[test]
fn joint_constraint_units_are_checked_term_by_term() {
    let result = problem_fixture(ProblemOptions {
        bad_constraint_units: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::InvalidNumeric {
            field: "affine constraint units",
            ..
        })
    ));
    log(
        "constraint-units",
        "pass",
        "coefficient times parameter must equal residual dimensions",
    );
}

#[test]
fn ordered_constraints_require_an_interval_witness() {
    for fixture in [1, 3] {
        assert!(matches!(
            problem_fixture(ProblemOptions {
                ordered_constraint_case: fixture,
                ..ProblemOptions::default()
            })
            .document,
            Err(IdentifiabilityError::InvalidNumeric {
                field: "ordered constraint feasibility",
                ..
            })
        ));
    }
    for fixture in [2, 4] {
        problem_fixture(ProblemOptions {
            ordered_constraint_case: fixture,
            ..ProblemOptions::default()
        })
        .document
        .unwrap_or_else(|error| panic!("feasible ordered fixture {fixture} refused: {error}"));
    }
    log(
        "ordered-constraint-feasibility",
        "pass",
        "non-strict boundary witnesses admit while impossible and strict-boundary chains refuse",
    );
}

#[test]
fn modeled_discrepancy_parameters_obey_exact_case_scope() {
    problem_fixture(ProblemOptions {
        modeled_discrepancy_case: 1,
        ..ProblemOptions::default()
    })
    .document
    .expect("case-a discrepancy uses a case-a parameter");
    assert!(matches!(
        problem_fixture(ProblemOptions {
            modeled_discrepancy_case: 2,
            ..ProblemOptions::default()
        })
        .document,
        Err(IdentifiabilityError::InvalidText {
            field: "discrepancy parameter/case scope",
            ..
        })
    ));
    log(
        "modeled-discrepancy-case-scope",
        "pass",
        "matching owner families cannot authorize a discrepancy parameter outside its exact case applicability",
    );
}

#[test]
fn discrepancy_inapplicability_is_closed_against_typed_experiment_origin() {
    let physical = retrospective_origin_fixture(
        ExperimentOrigin::Physical {
            apparatus_id: artifact("origin-apparatus"),
            facility_id: artifact("origin-facility"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::Physical,
    );
    admit_retrospective_origin_fixture(&physical)
        .expect("physical applicability basis matches physical experiment origin");

    let synthetic_for_physical = retrospective_origin_fixture(
        ExperimentOrigin::SyntheticHighFidelity {
            producer: artifact("forward-producer-a"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::Physical,
    );
    assert!(matches!(
        admit_retrospective_origin_fixture(&synthetic_for_physical),
        Err(IdentifiabilityError::SourceMismatch {
            field: "physical discrepancy/experiment origin",
        })
    ));

    let declared_synthetic = retrospective_origin_fixture(
        ExperimentOrigin::SyntheticHighFidelity {
            producer: artifact("forward-producer-a"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::DeclaredSynthetic {
            declared_producer: "forward-producer-a",
            stale_forward_binding: false,
            production_binding_key: "forward-a-production-binding",
        },
    );
    admit_retrospective_origin_fixture(&declared_synthetic).expect(
        "declared-synthetic basis binds typed synthetic origin, producer, and full forward SourceRef",
    );
    let alternate_production_binding = retrospective_origin_fixture(
        ExperimentOrigin::SyntheticHighFidelity {
            producer: artifact("forward-producer-a"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::DeclaredSynthetic {
            declared_producer: "forward-producer-a",
            stale_forward_binding: false,
            production_binding_key: "forward-a-production-binding-alternate",
        },
    );
    assert_ne!(
        unresolved_problem_identity(
            declared_synthetic
                .problem
                .document
                .as_ref()
                .expect("declared-synthetic document"),
        ),
        unresolved_problem_identity(
            alternate_production_binding
                .problem
                .document
                .as_ref()
                .expect("alternate production-binding document"),
        ),
        "the exact production-binding SourceRef key and discrepancy pointer must move ProblemId",
    );
    admit_retrospective_origin_fixture(&alternate_production_binding)
        .expect("alternate exact production-binding key remains structurally admissible");

    let physical_for_synthetic = retrospective_origin_fixture(
        ExperimentOrigin::Physical {
            apparatus_id: artifact("origin-apparatus"),
            facility_id: artifact("origin-facility"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::DeclaredSynthetic {
            declared_producer: "forward-producer-a",
            stale_forward_binding: false,
            production_binding_key: "forward-a-production-binding",
        },
    );
    assert!(matches!(
        admit_retrospective_origin_fixture(&physical_for_synthetic),
        Err(IdentifiabilityError::SourceMismatch {
            field: "declared-synthetic discrepancy/experiment origin",
        })
    ));

    let wrong_producer = retrospective_origin_fixture(
        ExperimentOrigin::SyntheticHighFidelity {
            producer: artifact("different-producer"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::DeclaredSynthetic {
            declared_producer: "forward-producer-a",
            stale_forward_binding: false,
            production_binding_key: "forward-a-production-binding",
        },
    );
    assert!(matches!(
        admit_retrospective_origin_fixture(&wrong_producer),
        Err(IdentifiabilityError::SourceMismatch {
            field: "declared-synthetic producer/forward-model binding",
        })
    ));

    let independent_implementation = retrospective_origin_fixture(
        ExperimentOrigin::SecondImplementation {
            producer: artifact("forward-producer-a"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::DeclaredSynthetic {
            declared_producer: "forward-producer-a",
            stale_forward_binding: false,
            production_binding_key: "forward-a-production-binding",
        },
    );
    assert!(matches!(
        admit_retrospective_origin_fixture(&independent_implementation),
        Err(IdentifiabilityError::SourceMismatch {
            field: "declared-synthetic discrepancy/experiment origin",
        })
    ));

    let stale_forward_binding = retrospective_origin_fixture(
        ExperimentOrigin::SyntheticHighFidelity {
            producer: artifact("forward-producer-a"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::DeclaredSynthetic {
            declared_producer: "forward-producer-a",
            stale_forward_binding: true,
            production_binding_key: "forward-a-production-binding",
        },
    );
    assert!(matches!(
        stale_forward_binding.problem.document,
        Err(IdentifiabilityError::SourceMismatch {
            field: "declared-synthetic producer/forward-model production binding",
        })
    ));
    log(
        "discrepancy-origin-closure",
        "pass",
        "physical and declared-synthetic inapplicability are checked against typed origin, producer identity, and an independently resolved full-SourceRef production binding",
    );
}

#[test]
fn validation_only_cases_require_physical_experiments() {
    for origin in [
        ExperimentOrigin::SyntheticHighFidelity {
            producer: artifact("validation-synthetic-producer"),
        },
        ExperimentOrigin::SecondImplementation {
            producer: artifact("validation-second-implementation"),
        },
    ] {
        let fixture = retrospective_origin_fixture(
            origin,
            CasePurpose::ValidationOnly,
            DiscrepancyOriginFixture::Uncharacterized,
        );
        assert!(matches!(
            admit_retrospective_origin_fixture(&fixture),
            Err(IdentifiabilityError::SourceMismatch {
                field: "validation-only case/physical experiment origin",
            })
        ));
    }
    log(
        "validation-only-origin",
        "pass",
        "synthetic and second-implementation artifacts cannot inherit the physics-validation meaning of ValidationOnly",
    );
}

#[test]
fn overlapping_assumed_gauges_require_one_exact_composition_hyperedge() {
    let result = problem_fixture(ProblemOptions {
        overlapping_gauges: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::InvalidText {
            field: "assumed gauge composition",
            ..
        })
    ));
    log(
        "overlapping-gauges",
        "pass",
        "simultaneously active assumed gauges refuse without one exact product/generated composition declaration",
    );
}

#[test]
fn admissible_domain_certificate_binds_exact_witness_preimage() {
    let certified = problem_fixture(ProblemOptions {
        gauge_case: 10,
        ..ProblemOptions::default()
    })
    .document
    .expect("witness-bound opaque domain certificate");
    let stale_membership_source = IdentifiabilityProblemDocument::try_new(
        certified.context_source().clone(),
        certified.material_source().clone(),
        certified.model_source().clone(),
        certified.graph_source().clone(),
        certified.joint_prior().cloned(),
        certified
            .sources()
            .values()
            .map(|source| {
                if source.kind() == SourceKind::AdmissibleDomainCertificate {
                    SourceRef::try_new(
                        source.key().clone(),
                        source.kind(),
                        hash("stale-membership-certificate-content"),
                        source.content_hash_domain(),
                        source.contract_version(),
                    )
                    .expect("stale membership source fixture")
                } else {
                    source.clone()
                }
            })
            .collect(),
        certified.parameters().values().cloned().collect(),
        certified.constraints().values().cloned().collect(),
        certified.admissible_domain().clone(),
        certified.cases().values().cloned().collect(),
        certified.influences().values().cloned().collect(),
        certified.gauges().values().cloned().collect(),
        certified.gauge_compositions().values().cloned().collect(),
        certified.joint_noise().clone(),
        certified.data_reuse().clone(),
    );
    assert!(matches!(
        stale_membership_source,
        Err(IdentifiabilityError::SourceMismatch {
            field: "admissible-domain membership certificate content"
        })
    ));
    log(
        "admissible-domain-certificate-binding",
        "pass",
        "opaque domain membership authority is bound to the exact witness/parameter/constraint/source preimage",
    );
}

fn gauge_reduction_binding(
    id: &str,
    claim: &TypedIdentifiabilityClaim,
    plan: GaugeReductionPlan,
    measure: GaugeMeasureSemantics,
) -> GaugeReductionBinding {
    GaugeReductionBinding::try_new(
        GaugeReductionId::try_new(id).expect("gauge reduction id"),
        GaugeActionReference::Single(GaugeClassId::try_new("fixture-gauge").expect("gauge id")),
        BTreeSet::from([claim.id().clone()]),
        plan,
        GaugeReductionStage::Root,
        measure,
    )
    .expect("gauge reduction binding fixture")
}

fn regular_quotient(prefix: &str) -> GaugeQuotientPlan {
    GaugeQuotientPlan::RegularAtlas {
        quotient_map: source(
            &format!("{prefix}-quotient-map"),
            SourceKind::GaugeQuotientMap,
            hash(&format!("{prefix}-quotient-map")),
        ),
        local_section_atlas: source(
            &format!("{prefix}-section-atlas"),
            SourceKind::GaugeSection,
            hash(&format!("{prefix}-section-atlas")),
        ),
        coverage: source(
            &format!("{prefix}-coverage"),
            SourceKind::GaugeSection,
            hash(&format!("{prefix}-coverage")),
        ),
    }
}

fn gauge_pushforward(prefix: &str) -> GaugeMeasureSemantics {
    GaugeMeasureSemantics::Pushforward {
        source_measure: source(
            &format!("{prefix}-source-measure"),
            SourceKind::ProbabilityMeasure,
            hash(&format!("{prefix}-source-measure")),
        ),
        reduced_measure: source(
            &format!("{prefix}-reduced-measure"),
            SourceKind::ProbabilityMeasure,
            hash(&format!("{prefix}-reduced-measure")),
        ),
        transport: source(
            &format!("{prefix}-pushforward"),
            SourceKind::GaugeMeasureTransport,
            hash(&format!("{prefix}-pushforward")),
        ),
        jacobian_or_disintegration: source(
            &format!("{prefix}-jacobian"),
            SourceKind::GaugeMeasureTransport,
            hash(&format!("{prefix}-jacobian")),
        ),
    }
}

fn gauge_slice_plan(
    prefix: &str,
    support: BTreeSet<ParameterRoleId>,
    expected_codimension: GaugeSliceCodimension,
    coverage_kind: SourceKind,
) -> GaugeSlicePlan {
    GaugeSlicePlan::try_new(
        support,
        source(
            &format!("{prefix}-constraint"),
            SourceKind::Constraint,
            hash(&format!("{prefix}-constraint")),
        ),
        expected_codimension,
        source(
            &format!("{prefix}-coverage"),
            coverage_kind,
            hash(&format!("{prefix}-coverage")),
        ),
    )
    .expect("gauge slice plan fixture")
}

#[test]
fn gauge_slices_are_execution_plans_with_exact_support_codimension_and_coverage() {
    let problem = admit_fixture(problem_fixture(ProblemOptions {
        gauge_case: 4,
        ..ProblemOptions::default()
    }));
    let claim = structural_claim(
        &problem,
        "fixed-slice-quotient",
        FiberStructure::OrbitQuotientUnique {
            action: GaugeActionReference::Single(
                GaugeClassId::try_new("fixture-gauge").expect("gauge id"),
            ),
        },
        ClaimSubject::Parameter(role("yield_stress")),
        ClaimScope::WholeCampaign,
    );
    let valid_slice = gauge_slice_plan(
        "fixed-slice",
        BTreeSet::from([role("yield_stress"), role("hardening_modulus")]),
        GaugeSliceCodimension::FixedFinite { codimension: 1 },
        SourceKind::GaugeSection,
    );
    let valid = execution_for_claim_with_gauge_reductions(
        &problem,
        claim.clone(),
        vec![gauge_reduction_binding(
            "fixed-slice-reduction",
            &claim,
            GaugeReductionPlan::Slice { slice: valid_slice },
            gauge_pushforward("fixed-slice"),
        )],
    )
    .expect("transverse fixed-codimension execution slice");
    assert_eq!(
        IdentifiabilityExecutionPlan::from_canonical_bytes(
            &valid.canonical_bytes().expect("slice execution bytes"),
            &problem,
            valid.source_authority(),
        )
        .expect("slice execution replay"),
        valid,
    );

    let invalid_support = gauge_slice_plan(
        "invalid-support",
        BTreeSet::from([role("not-in-gauge-carrier")]),
        GaugeSliceCodimension::FixedFinite { codimension: 1 },
        SourceKind::GaugeSection,
    );
    assert!(matches!(
        execution_for_claim_with_gauge_reductions(
            &problem,
            claim.clone(),
            vec![gauge_reduction_binding(
                "invalid-support-reduction",
                &claim,
                GaugeReductionPlan::Slice {
                    slice: invalid_support,
                },
                gauge_pushforward("invalid-support"),
            )],
        ),
        Err(IdentifiabilityError::InvalidText {
            field: "execution gauge slice support",
            ..
        })
    ));
    let invalid_codimension = gauge_slice_plan(
        "invalid-codimension",
        BTreeSet::from([role("yield_stress")]),
        GaugeSliceCodimension::FixedFinite { codimension: 2 },
        SourceKind::GaugeSection,
    );
    assert!(matches!(
        execution_for_claim_with_gauge_reductions(
            &problem,
            claim.clone(),
            vec![gauge_reduction_binding(
                "invalid-codimension-reduction",
                &claim,
                GaugeReductionPlan::Slice {
                    slice: invalid_codimension,
                },
                gauge_pushforward("invalid-codimension"),
            )],
        ),
        Err(IdentifiabilityError::InvalidText {
            field: "execution gauge slice codimension",
            ..
        })
    ));
    let invalid_coverage = gauge_slice_plan(
        "invalid-coverage",
        BTreeSet::from([role("yield_stress")]),
        GaugeSliceCodimension::FixedFinite { codimension: 1 },
        SourceKind::Assumption,
    );
    assert!(matches!(
        execution_for_claim_with_gauge_reductions(
            &problem,
            claim.clone(),
            vec![gauge_reduction_binding(
                "invalid-coverage-reduction",
                &claim,
                GaugeReductionPlan::Slice {
                    slice: invalid_coverage,
                },
                gauge_pushforward("invalid-coverage"),
            )],
        ),
        Err(IdentifiabilityError::InvalidText {
            field: "execution gauge slice sources",
            ..
        })
    ));

    let stratified_problem = admit_fixture(problem_fixture(ProblemOptions {
        gauge_case: 9,
        ..ProblemOptions::default()
    }));
    let stratified_claim = structural_claim(
        &stratified_problem,
        "stratified-slice-quotient",
        FiberStructure::OrbitQuotientUnique {
            action: GaugeActionReference::Single(
                GaugeClassId::try_new("fixture-gauge").expect("gauge id"),
            ),
        },
        ClaimSubject::Parameter(role("yield_stress")),
        ClaimScope::WholeCampaign,
    );
    let orbit_profile =
        stratified_problem.document().sources()[&source_key("fixture-gauge-strata")].clone();
    let stratified_slice = gauge_slice_plan(
        "stratified-slice",
        BTreeSet::from([role("yield_stress"), role("hardening_modulus")]),
        GaugeSliceCodimension::Stratified {
            profile: orbit_profile,
        },
        SourceKind::GaugeSection,
    );
    execution_for_claim_with_gauge_reductions(
        &stratified_problem,
        stratified_claim.clone(),
        vec![gauge_reduction_binding(
            "stratified-slice-reduction",
            &stratified_claim,
            GaugeReductionPlan::Slice {
                slice: stratified_slice,
            },
            gauge_pushforward("stratified-slice"),
        )],
    )
    .expect("stratified execution slice binds exact orbit profile");

    log(
        "execution-gauge-slice-binding",
        "pass",
        "execution-time slices bind exact action support, effective-orbit codimension, section coverage, stratified profiles, measure transport, authority, and canonical replay",
    );
}

#[test]
fn orbit_fibers_distinguish_pure_discrete_mixed_residual_and_full_quotients() {
    let pure_discrete = admit_fixture(problem_fixture(ProblemOptions {
        gauge_case: 1,
        ..ProblemOptions::default()
    }));
    let discrete_claim = structural_claim(
        &pure_discrete,
        "pure-discrete-orbit",
        FiberStructure::DiscreteOrbit {
            action: GaugeActionReference::Single(
                GaugeClassId::try_new("fixture-gauge").expect("gauge id"),
            ),
        },
        ClaimSubject::Parameter(role("yield_stress")),
        ClaimScope::WholeCampaign,
    );
    assert!(matches!(
        execution_for_claim(
            &pure_discrete,
            false,
            17,
            1.0e-10,
            false,
            discrete_claim.clone(),
        ),
        Err(IdentifiabilityError::Cardinality {
            field: "gauge reduction coverage",
            ..
        })
    ));
    execution_for_claim_with_gauge_reductions(
        &pure_discrete,
        discrete_claim.clone(),
        vec![gauge_reduction_binding(
            "pure-discrete-unreduced",
            &discrete_claim,
            GaugeReductionPlan::Unreduced {
                reason: "the analyzer evaluates the exact declared discrete orbit without quotienting it"
                    .to_string(),
            },
            GaugeMeasureSemantics::NotApplicable {
                reason: "an unreduced structural discrete-orbit proposition performs no measure transport"
                    .to_string(),
            },
        )],
    )
    .expect("pure discrete orbit has an exact explicit unreduced plan");

    let mixed_retained = admit_fixture(problem_fixture(ProblemOptions {
        gauge_case: 2,
        ..ProblemOptions::default()
    }));
    let wrong_discrete_claim = structural_claim(
        &mixed_retained,
        "mixed-is-not-pure-discrete",
        FiberStructure::DiscreteOrbit {
            action: GaugeActionReference::Single(
                GaugeClassId::try_new("fixture-gauge").expect("gauge id"),
            ),
        },
        ClaimSubject::Parameter(role("yield_stress")),
        ClaimScope::WholeCampaign,
    );
    assert!(matches!(
        execution_for_claim_with_gauge_reductions(
            &mixed_retained,
            wrong_discrete_claim.clone(),
            vec![gauge_reduction_binding(
                "wrong-pure-discrete-plan",
                &wrong_discrete_claim,
                GaugeReductionPlan::Unreduced {
                    reason: "fixture retains the mixed action".to_string(),
                },
                GaugeMeasureSemantics::NotApplicable {
                    reason: "unreduced structural analysis has no quotient measure".to_string(),
                },
            )],
        ),
        Err(IdentifiabilityError::InvalidNumeric {
            field: "discrete-orbit fiber",
            ..
        })
    ));
    let mixed_claim = structural_claim(
        &mixed_retained,
        "declared-mixed-orbit",
        FiberStructure::MixedOrbit {
            action: GaugeActionReference::Single(
                GaugeClassId::try_new("fixture-gauge").expect("gauge id"),
            ),
        },
        ClaimSubject::Parameter(role("yield_stress")),
        ClaimScope::WholeCampaign,
    );
    let unreduced_execution = execution_for_claim_with_gauge_reductions(
        &mixed_retained,
        mixed_claim.clone(),
        vec![gauge_reduction_binding(
            "mixed-unreduced",
            &mixed_claim,
            GaugeReductionPlan::Unreduced {
                reason: "the mixed physical orbit is analyzed without computational quotienting"
                    .to_string(),
            },
            GaugeMeasureSemantics::NotApplicable {
                reason: "the unreduced structural proposition changes no probability measure"
                    .to_string(),
            },
        )],
    )
    .expect("mixed orbit admits only after explicit unreduced coverage");
    let mixed_residual_plan = GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
        reduction: ContinuousGaugeReductionPlan::Quotient {
            quotient: regular_quotient("mixed-residual"),
        },
        normal_subgroup: source(
            "mixed-residual-normal-subgroup",
            SourceKind::GaugeSubgroupCertificate,
            hash("mixed-residual-normal-subgroup"),
        ),
        factor_extension: source(
            "mixed-residual-factor-extension",
            SourceKind::GaugeReductionLaw,
            hash("mixed-residual-factor-extension"),
        ),
        residual_quotient_action: source(
            "mixed-residual-action",
            SourceKind::GaugeResidualAction,
            hash("mixed-residual-action"),
        ),
        compatibility: source(
            "mixed-residual-compatibility",
            SourceKind::GaugeReductionLaw,
            hash("mixed-residual-compatibility"),
        ),
    };
    let mixed_execution = execution_for_claim_with_gauge_reductions(
        &mixed_retained,
        mixed_claim.clone(),
        vec![gauge_reduction_binding(
            "mixed-residual-reduction",
            &mixed_claim,
            mixed_residual_plan,
            gauge_pushforward("mixed-residual"),
        )],
    )
    .expect("continuous normal subgroup is reduced with an explicit residual discrete action");
    assert_ne!(
        unreduced_execution.id().expect("unreduced execution id"),
        mixed_execution.id().expect("residual execution id"),
        "unreduced and continuously reduced plans are distinct execution semantics",
    );
    assert_eq!(
        IdentifiabilityExecutionPlan::from_canonical_bytes(
            &mixed_execution
                .canonical_bytes()
                .expect("mixed-orbit execution bytes"),
            &mixed_retained,
            mixed_execution.source_authority(),
        )
        .expect("mixed-orbit execution replay"),
        mixed_execution,
    );
    let quotient_claim = structural_claim(
        &mixed_retained,
        "fully-quotiented-mixed-orbit",
        FiberStructure::OrbitQuotientUnique {
            action: GaugeActionReference::Single(
                GaugeClassId::try_new("fixture-gauge").expect("gauge id"),
            ),
        },
        ClaimSubject::Parameter(role("yield_stress")),
        ClaimScope::WholeCampaign,
    );
    execution_for_claim_with_gauge_reductions(
        &mixed_retained,
        quotient_claim.clone(),
        vec![gauge_reduction_binding(
            "mixed-full-quotient",
            &quotient_claim,
            GaugeReductionPlan::Quotient {
                quotient: regular_quotient("mixed-full"),
            },
            gauge_pushforward("mixed-full"),
        )],
    )
    .expect("full mixed quotient has exact atlas and measure semantics");
    log(
        "orbit-fiber-component-semantics",
        "pass",
        "pure discrete, unreduced mixed, continuously reduced mixed residual, and fully quotiented mixed semantics remain explicit and canonically distinct",
    );
}

#[test]
fn self_correlation_is_not_an_identifiability_route() {
    let result = problem_fixture(ProblemOptions {
        self_correlation: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::InvalidNumeric {
            field: "correlation functional",
            ..
        })
    ));
    log(
        "self-correlation",
        "pass",
        "constant self-correlation cannot masquerade as sensitivity",
    );
}

#[test]
fn derivative_units_are_derived_from_functional_and_parameter() {
    let problem = problem_fixture(ProblemOptions::default())
        .document
        .expect("problem");
    let quantity = problem
        .influence_derivative_quantity(&InfluenceId::try_new("yield-to-stress").expect("influence"))
        .expect("derived quantity");
    assert_eq!(quantity, QuantitySpec::dimensional(DIMENSIONLESS));
    log(
        "derived-influence-units",
        "pass",
        "caller cannot inject contradictory derivative dimensions",
    );
}

#[test]
fn dense_correlation_refuses_marginals_without_finite_standard_deviation() {
    let result = problem_fixture(ProblemOptions {
        dense_with_bounded_marginal: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::Covariance { .. })
    ));
    log(
        "dense-correlation-marginals",
        "pass",
        "bounded noise was not silently converted into Gaussian scale",
    );
}

#[test]
fn accidental_raw_experiment_reuse_refuses_under_disjoint_policy() {
    let result = problem_fixture(ProblemOptions {
        retrospective_reuse: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(matches!(
        result,
        Err(IdentifiabilityError::InvalidText {
            field: "data reuse policy",
            ..
        })
    ));
    log(
        "accidental-data-reuse",
        "pass",
        "same experiment cannot be double-counted under Disjoint",
    );
}

#[test]
fn declared_raw_reuse_requires_joint_likelihood_and_justification() {
    let result = problem_fixture(ProblemOptions {
        retrospective_reuse: true,
        declared_sharing: true,
        ..ProblemOptions::default()
    })
    .document;
    assert!(result.is_ok());
    log(
        "declared-data-reuse",
        "pass",
        "shared campaign admitted only through explicit group and likelihood",
    );
}

#[test]
fn wrong_concrete_source_hash_refuses_problem_identity() {
    for mutation in 1..=3 {
        let fixture = problem_fixture(ProblemOptions {
            context_contract_mutation: mutation,
            ..ProblemOptions::default()
        });
        let ProblemFixture {
            context,
            material,
            model,
            document,
            ..
        } = fixture;
        let document = document.expect("structural document");
        let opaque = opaque_resolutions(&document);
        let result = AdmittedIdentifiabilityProblem::resolve_and_admit(
            document,
            ProblemSourceBundle::new(&context, &material, &model, BTreeMap::new(), opaque),
        );
        assert!(matches!(
            result,
            Err(IdentifiabilityError::SourceMismatch { .. })
        ));
    }
    log(
        "wrong-concrete-source",
        "pass",
        "typed context hash, digest domain, and contract version are resolver-derived",
    );
}

#[test]
fn opaque_resolution_cannot_replay_across_hash_domains() {
    let good = problem_fixture(ProblemOptions::default());
    let good_document = good.document.expect("good document");
    let resolutions = opaque_resolutions(&good_document);
    let bad = problem_fixture(ProblemOptions {
        alternate_graph_domain: true,
        ..ProblemOptions::default()
    });
    let bad_document = bad.document.expect("alternate-domain document");
    let result = AdmittedIdentifiabilityProblem::resolve_and_admit(
        bad_document,
        ProblemSourceBundle::new(
            &bad.context,
            &bad.material,
            &bad.model,
            BTreeMap::new(),
            resolutions,
        ),
    );
    assert!(matches!(
        result,
        Err(IdentifiabilityError::SourceMismatch {
            field: "opaque source resolution"
        })
    ));
    log(
        "cross-domain-resolution-replay",
        "pass",
        "a digest verified under one domain cannot authorize an equal digest under another",
    );
}

#[test]
fn unverified_opaque_source_cannot_mint_problem_id() {
    let fixture = problem_fixture(ProblemOptions::default());
    let ProblemFixture {
        context,
        material,
        model,
        graph: _,
        document,
    } = fixture;
    let document = document.expect("document");
    let entries = document
        .sources()
        .values()
        .filter(|source| {
            !matches!(
                source.kind(),
                SourceKind::ContextOfUse
                    | SourceKind::MaterialCard
                    | SourceKind::ConstitutiveModelCard
                    | SourceKind::ExperimentArtifact
                    | SourceKind::CalibrationSplit
            )
        })
        .map(|source| {
            if source.key().as_str() == "forward-a" {
                SourceResolution::unresolved(source, "resolver has not fetched this artifact")
                    .expect("unresolved diagnostic")
            } else {
                let preimage = if source.key().as_str() == "graph" {
                    b"constitutive-graph".as_slice()
                } else {
                    source.key().as_str().as_bytes()
                };
                SourceResolution::verify(source, preimage, AuthorityDisposition::ContentVerified)
                    .expect("resolution")
            }
        })
        .collect();
    let opaque = SourceResolutionSet::try_new(entries).expect("resolution set");
    let result = AdmittedIdentifiabilityProblem::resolve_and_admit(
        document,
        ProblemSourceBundle::new(&context, &material, &model, BTreeMap::new(), opaque),
    );
    assert!(matches!(
        result,
        Err(IdentifiabilityError::InvalidText {
            field: "source authority",
            ..
        })
    ));
    log(
        "unverified-source",
        "pass",
        "content reference alone did not grant source authority",
    );
}

#[test]
fn external_trust_receipts_bind_exact_source_and_typed_subject_artifact() {
    let expected = source(
        "trust-subject-expected",
        SourceKind::Assumption,
        hash("trust-subject-expected"),
    );
    let different = source(
        "trust-subject-different",
        SourceKind::Assumption,
        hash("trust-subject-different"),
    );
    assert!(matches!(
        SourceResolution::verify(
            &expected,
            b"trust-subject-expected",
            external_trust("cross-source-trust-receipt", &different),
        ),
        Err(IdentifiabilityError::SourceMismatch {
            field: "trust receipt subject/source resolution",
        })
    ));

    let fixture = retrospective_origin_fixture(
        ExperimentOrigin::Physical {
            apparatus_id: artifact("trust-artifact-apparatus"),
            facility_id: artifact("trust-artifact-facility"),
        },
        CasePurpose::Calibration,
        DiscrepancyOriginFixture::Uncharacterized,
    );
    let document = fixture
        .problem
        .document
        .clone()
        .expect("trust subject artifact document");
    let split_key = source_key("split-a");
    let split_source = document.sources()[&split_key].clone();
    let malformed_blind_namespace = SourceRef::try_new(
        source_key("malformed-blind-receipt"),
        SourceKind::EvidenceReceipt,
        hash("malformed-blind-receipt"),
        BLIND_RELEASE_TRUST_RECEIPT_DOMAIN,
        BLIND_RELEASE_TRUST_RECEIPT_VERSION,
    )
    .expect("syntactically valid blind namespace source");
    assert!(matches!(
        TrustReceiptRef::try_new(
            malformed_blind_namespace,
            split_source.clone(),
            TrustAuthentication::Unauthenticated,
        ),
        Err(IdentifiabilityError::InvalidText {
            field: "trust receipt subject artifact",
            ..
        })
    ));
    let stale_blind_namespace = SourceRef::try_new(
        source_key("stale-blind-receipt"),
        SourceKind::EvidenceReceipt,
        hash("stale-blind-receipt"),
        BLIND_RELEASE_TRUST_RECEIPT_DOMAIN,
        BLIND_RELEASE_TRUST_RECEIPT_VERSION + 1,
    )
    .expect("syntactically valid stale blind namespace source");
    assert!(matches!(
        TrustReceiptRef::try_new(
            stale_blind_namespace,
            split_source.clone(),
            TrustAuthentication::Unauthenticated,
        ),
        Err(IdentifiabilityError::VersionMismatch {
            field: "blind-release trust receipt",
            ..
        })
    ));
    assert!(matches!(
        TrustReceiptRef::blind_release(
            &expected,
            artifact("not-a-calibration-split"),
            hash("wrong-kind-blind-receipt"),
        ),
        Err(IdentifiabilityError::InvalidText {
            field: "trust receipt subject artifact",
            ..
        })
    ));
    let generic_split_authority = external_trust("generic-non-blind-split-trust", &split_source);
    let generic_bundle = ProblemSourceBundle::new(
        &fixture.problem.context,
        &fixture.problem.material,
        &fixture.problem.model,
        BTreeMap::from([
            (
                case_id("a"),
                CaseSourceBundle::new(&fixture.experiment, &fixture.split_a),
            ),
            (
                case_id("b"),
                CaseSourceBundle::new(&fixture.experiment, &fixture.split_b),
            ),
        ]),
        opaque_resolutions(&document),
    )
    .with_concrete_authority(vec![(split_key.clone(), generic_split_authority)])
    .expect("generic non-blind split authority envelope");
    let generic_admission = AdmittedIdentifiabilityProblem::resolve_and_admit(
        document.clone(),
        generic_bundle,
    )
    .expect(
        "generic exact-SourceRef trust may authorize a non-blind split without artifact metadata",
    );

    let correct_blind_authority = AuthorityDisposition::ExternalTrustReceipt {
        trust_receipt: TrustReceiptRef::blind_release(
            &split_source,
            fixture.split_a.id().clone(),
            hash("correct-subject-artifact-receipt"),
        )
        .expect("exact blind-release subject artifact"),
    };
    let correct_blind_bundle = ProblemSourceBundle::new(
        &fixture.problem.context,
        &fixture.problem.material,
        &fixture.problem.model,
        BTreeMap::from([
            (
                case_id("a"),
                CaseSourceBundle::new(&fixture.experiment, &fixture.split_a),
            ),
            (
                case_id("b"),
                CaseSourceBundle::new(&fixture.experiment, &fixture.split_b),
            ),
        ]),
        opaque_resolutions(&document),
    )
    .with_concrete_authority(vec![(split_key.clone(), correct_blind_authority)])
    .expect("bounded correct blind authority envelope");
    let blind_admission =
        AdmittedIdentifiabilityProblem::resolve_and_admit(document.clone(), correct_blind_bundle)
            .expect("exact blind-release subject artifact admits");
    assert_eq!(generic_admission.id(), blind_admission.id());
    assert_ne!(
        generic_admission.source_admission_id(),
        blind_admission.source_admission_id(),
        "typed blind-release authority must move only the authority envelope",
    );

    let policy_receipt = source(
        "issuer-policy-trust-receipt",
        SourceKind::EvidenceReceipt,
        hash("issuer-policy-trust-receipt"),
    );
    let trust_policy = source(
        "declared-trust-policy",
        SourceKind::Assumption,
        hash("declared-trust-policy"),
    );
    assert!(matches!(
        TrustReceiptRef::try_new(
            policy_receipt.clone(),
            split_source.clone(),
            TrustAuthentication::IssuerPolicy {
                issuer: artifact("policy-issuer"),
                trust_policy: source(
                    "wrong-kind-trust-policy",
                    SourceKind::EvidenceReceipt,
                    hash("wrong-kind-trust-policy"),
                ),
            },
        ),
        Err(IdentifiabilityError::InvalidText {
            field: "trust receipt policy",
            ..
        })
    ));
    let admit_with_trust_receipt = |trust_receipt| {
        let bundle = ProblemSourceBundle::new(
            &fixture.problem.context,
            &fixture.problem.material,
            &fixture.problem.model,
            BTreeMap::from([
                (
                    case_id("a"),
                    CaseSourceBundle::new(&fixture.experiment, &fixture.split_a),
                ),
                (
                    case_id("b"),
                    CaseSourceBundle::new(&fixture.experiment, &fixture.split_b),
                ),
            ]),
            opaque_resolutions(&document),
        )
        .with_concrete_authority(vec![(
            split_key.clone(),
            AuthorityDisposition::ExternalTrustReceipt { trust_receipt },
        )])
        .expect("bounded issuer-policy authority envelope");
        AdmittedIdentifiabilityProblem::resolve_and_admit(document.clone(), bundle)
            .expect("issuer-policy declaration admits structurally")
    };
    let unauthenticated_declaration = admit_with_trust_receipt(
        TrustReceiptRef::try_new(
            policy_receipt.clone(),
            split_source.clone(),
            TrustAuthentication::Unauthenticated,
        )
        .expect("unauthenticated generic trust declaration"),
    );
    let issuer_policy_declaration = admit_with_trust_receipt(
        TrustReceiptRef::try_new(
            policy_receipt,
            split_source.clone(),
            TrustAuthentication::IssuerPolicy {
                issuer: artifact("policy-issuer"),
                trust_policy,
            },
        )
        .expect("issuer-policy-bound trust declaration"),
    );
    assert_eq!(
        unauthenticated_declaration.id(),
        issuer_policy_declaration.id(),
    );
    assert_ne!(
        unauthenticated_declaration.source_admission_id(),
        issuer_policy_declaration.source_admission_id(),
        "declared issuer/policy semantics must move SourceAdmissionId without claiming issuer verification",
    );

    let wrong_subject_artifact = AuthorityDisposition::ExternalTrustReceipt {
        trust_receipt: TrustReceiptRef::blind_release(
            &document.sources()[&split_key],
            artifact("not-the-resolved-split"),
            hash("wrong-subject-artifact-receipt"),
        )
        .expect("structurally typed blind-release receipt"),
    };
    let opaque = opaque_resolutions(&document);
    let bundle = ProblemSourceBundle::new(
        &fixture.problem.context,
        &fixture.problem.material,
        &fixture.problem.model,
        BTreeMap::from([
            (
                case_id("a"),
                CaseSourceBundle::new(&fixture.experiment, &fixture.split_a),
            ),
            (
                case_id("b"),
                CaseSourceBundle::new(&fixture.experiment, &fixture.split_b),
            ),
        ]),
        opaque,
    )
    .with_concrete_authority(vec![(split_key, wrong_subject_artifact)])
    .expect("bounded concrete authority envelope");
    assert!(matches!(
        AdmittedIdentifiabilityProblem::resolve_and_admit(document, bundle),
        Err(IdentifiabilityError::SourceMismatch {
            field: "trust receipt subject artifact/source resolution",
        })
    ));
    log(
        "typed-trust-subject-binding",
        "pass",
        "trust receipts cannot replay across SourceRefs, and blind-release subject artifacts must equal the exact resolved split ID",
    );
}

#[test]
fn source_bundle_public_ingress_enforces_aggregate_cardinality_before_admission_work() {
    let fixture = problem_fixture(ProblemOptions::default());
    let authority_overflow = (0..=MAX_IDENTIFIABILITY_ITEMS)
        .map(|index| {
            (
                source_key(&format!("concrete-authority-{index:04}")),
                AuthorityDisposition::ContentVerified,
            )
        })
        .collect();
    assert!(matches!(
        ProblemSourceBundle::new(
            &fixture.context,
            &fixture.material,
            &fixture.model,
            BTreeMap::new(),
            SourceResolutionSet::default(),
        )
        .with_concrete_authority(authority_overflow),
        Err(IdentifiabilityError::Cardinality {
            field: "concrete source authority",
            ..
        })
    ));

    let (experiment, split, _) = retrospective_artifacts(ExperimentOrigin::Physical {
        apparatus_id: artifact("cardinality-apparatus"),
        facility_id: artifact("cardinality-facility"),
    });
    let case_overflow = (0..=MAX_IDENTIFIABILITY_ITEMS)
        .map(|index| {
            (
                case_id(&format!("concrete-case-{index:04}")),
                CaseSourceBundle::new(&experiment, &split),
            )
        })
        .collect();
    let document = fixture.document.expect("cardinality document");
    assert!(matches!(
        AdmittedIdentifiabilityProblem::resolve_and_admit(
            document,
            ProblemSourceBundle::new(
                &fixture.context,
                &fixture.material,
                &fixture.model,
                case_overflow,
                SourceResolutionSet::default(),
            ),
        ),
        Err(IdentifiabilityError::Cardinality {
            field: "retrospective case source bundles",
            ..
        })
    ));
    log(
        "source-bundle-cardinality",
        "pass",
        "max-plus-one case bundles and concrete authority entries refuse before hashing, resolution, or derived-lineage work",
    );
}

#[test]
fn execution_and_assessment_collection_caps_match_the_canonical_decoder() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let baseline = execution(&problem, false, 17, 1.0e-10, false).expect("baseline execution");
    let attempt_execution = |parts: ExecutionParts| {
        IdentifiabilityExecutionPlan::try_new(
            parts.header,
            &problem,
            parts.analyzer,
            parts.build,
            parts.derivative_provider,
            parts.claim_requests,
            parts.actions,
            parts.gauge_reductions,
            parts.numerical,
            parts.initialization,
            parts.stopping,
            parts.determinism,
            baseline.source_authority().clone(),
        )
    };

    let mut oversized_actions = ExecutionParts::from_plan(&baseline);
    oversized_actions.actions = vec![
        (
            role("overflow-action"),
            ParameterExecutionAction::Conditioned
        );
        MAX_IDENTIFIABILITY_ITEMS + 1
    ];
    assert!(matches!(
        attempt_execution(oversized_actions),
        Err(IdentifiabilityError::Cardinality {
            field: "execution parameter actions",
            ..
        })
    ));

    let claim = default_claim(&problem);
    let reduction = gauge_reduction_binding(
        "overflow-reduction",
        &claim,
        GaugeReductionPlan::Unreduced {
            reason: "cardinality sentinel deliberately leaves the action unreduced".to_string(),
        },
        GaugeMeasureSemantics::NotApplicable {
            reason: "no quotient or slice is applied".to_string(),
        },
    );
    let mut oversized_reductions = ExecutionParts::from_plan(&baseline);
    oversized_reductions.gauge_reductions = vec![reduction; MAX_IDENTIFIABILITY_ITEMS + 1];
    assert!(matches!(
        attempt_execution(oversized_reductions),
        Err(IdentifiabilityError::Cardinality {
            field: "execution gauge reductions",
            ..
        })
    ));

    let baseline_assessment = assessment(&problem, &baseline, "collection-cap-receipt");
    let parts = AssessmentParts::from_assessment(&baseline_assessment);
    let evidence_row = parts.evidence[0].clone();
    assert!(matches!(
        IdentifiabilityAssessment::try_new(
            parts.header.clone(),
            &problem,
            &baseline,
            parts.claims.clone(),
            vec![evidence_row; MAX_IDENTIFIABILITY_ITEMS + 1],
            parts.source_authority.clone(),
        ),
        Err(IdentifiabilityError::Cardinality {
            field: "claim assessments",
            ..
        })
    ));

    let (claim_id, method, receipt, metric, nondimensionalization, certified_error_bound) = {
        let (claim_id, conclusion) = &parts.evidence[0];
        let ClaimAssessment::ClaimedEstablished {
            method,
            receipt,
            metric,
            nondimensionalization,
            certified_error_bound,
            ..
        } = conclusion
        else {
            panic!("baseline assessment must be claimed established");
        };
        (
            claim_id.clone(),
            method.clone(),
            receipt.clone(),
            metric.clone(),
            nondimensionalization.clone(),
            *certified_error_bound,
        )
    };
    let oversized_gauge_resolutions = (0..=MAX_IDENTIFIABILITY_ITEMS)
        .map(|index| {
            let action = GaugeActionReference::Single(
                GaugeClassId::try_new(format!("overflow-gauge-{index:04}"))
                    .expect("bounded gauge id"),
            );
            (
                action.clone(),
                GaugeResolutionEvidence::new(
                    action,
                    GaugeResolutionDisposition::CandidateRefuted,
                    method.clone(),
                    receipt.clone(),
                ),
            )
        })
        .collect();
    assert!(matches!(
        IdentifiabilityAssessment::try_new(
            parts.header,
            &problem,
            &baseline,
            parts.claims,
            vec![(
                claim_id,
                ClaimAssessment::ClaimedEstablished {
                    method,
                    receipt,
                    metric,
                    nondimensionalization,
                    certified_error_bound,
                    gauge_resolutions: oversized_gauge_resolutions,
                },
            )],
            parts.source_authority,
        ),
        Err(IdentifiabilityError::Cardinality {
            field: "positive-claim gauge resolutions",
            ..
        })
    ));
    log(
        "canonical-collection-caps",
        "pass",
        "execution actions/reductions and assessment evidence/gauge-resolution maps refuse MAX+1 at public ingress, matching their decoder count bounds",
    );
}

#[test]
fn problem_and_source_admission_identities_separate_question_from_trust_envelope() {
    let content_only =
        admit_fixture_with_authority(problem_fixture(ProblemOptions::default()), false);
    let authenticated =
        admit_fixture_with_authority(problem_fixture(ProblemOptions::default()), true);
    assert_eq!(content_only.id(), authenticated.id());
    assert_ne!(
        content_only.source_admission_id(),
        authenticated.source_admission_id()
    );
    log(
        "problem-vs-authority-identity",
        "pass",
        "trust receipt moves authority envelope without rewriting physical question",
    );
}

#[test]
fn identifiability_problem_identity_bindings_have_exact_mutation_evidence() {
    let baseline_document = problem_fixture(ProblemOptions::default())
        .document
        .expect("baseline problem");
    let baseline = unresolved_problem_identity(&baseline_document);
    let mut witness_values = baseline_document
        .admissible_domain()
        .values()
        .iter()
        .map(|(role, value)| (role.clone(), *value))
        .collect::<Vec<_>>();
    let (_, yield_witness) = witness_values
        .iter_mut()
        .find(|(role, _)| role == &role("yield_stress"))
        .expect("yield witness");
    *yield_witness = 2.0e6;
    let witness_variant = IdentifiabilityProblemDocument::try_new(
        baseline_document.context_source().clone(),
        baseline_document.material_source().clone(),
        baseline_document.model_source().clone(),
        baseline_document.graph_source().clone(),
        baseline_document.joint_prior().cloned(),
        baseline_document.sources().values().cloned().collect(),
        baseline_document.parameters().values().cloned().collect(),
        baseline_document.constraints().values().cloned().collect(),
        AdmissibleDomainWitness::try_new(witness_values, None)
            .expect("mutated admissible-domain witness"),
        baseline_document.cases().values().cloned().collect(),
        baseline_document.influences().values().cloned().collect(),
        baseline_document.gauges().values().cloned().collect(),
        baseline_document
            .gauge_compositions()
            .values()
            .cloned()
            .collect(),
        baseline_document.joint_noise().clone(),
        baseline_document.data_reuse().clone(),
    )
    .expect("independent admissible-domain mutation");
    let variants = [
        (
            "context_source",
            rekey_problem_root(
                problem_fixture(ProblemOptions::default()),
                ProblemRoot::Context,
                "context-rekeyed",
            )
            .document
            .expect("context-root mutation"),
        ),
        (
            "material_source",
            rekey_problem_root(
                problem_fixture(ProblemOptions::default()),
                ProblemRoot::Material,
                "material-rekeyed",
            )
            .document
            .expect("material-root mutation"),
        ),
        (
            "model_source",
            rekey_problem_root(
                problem_fixture(ProblemOptions::default()),
                ProblemRoot::Model,
                "model-rekeyed",
            )
            .document
            .expect("model-root mutation"),
        ),
        (
            "graph_source",
            rekey_problem_root(
                problem_fixture(ProblemOptions::default()),
                ProblemRoot::Graph,
                "graph-rekeyed",
            )
            .document
            .expect("graph-root mutation"),
        ),
        (
            "sources",
            problem_fixture(ProblemOptions {
                alternate_graph_domain: true,
                ..ProblemOptions::default()
            })
            .document
            .expect("source-registry mutation"),
        ),
        (
            "parameters",
            problem_fixture(ProblemOptions {
                parameter_prior_version: 2,
                ..ProblemOptions::default()
            })
            .document
            .expect("parameter-registry mutation"),
        ),
        (
            "constraints",
            problem_fixture(ProblemOptions {
                valid_constraint: true,
                ..ProblemOptions::default()
            })
            .document
            .expect("constraint-registry mutation"),
        ),
        ("admissible_domain", witness_variant),
        (
            "cases",
            problem_fixture(ProblemOptions {
                second_case_complementary: true,
                ..ProblemOptions::default()
            })
            .document
            .expect("case-registry mutation"),
        ),
        (
            "influences",
            problem_fixture(ProblemOptions {
                yield_log_scale: true,
                ..ProblemOptions::default()
            })
            .document
            .expect("influence-registry mutation"),
        ),
        (
            "gauges",
            problem_fixture(ProblemOptions {
                one_gauge: true,
                ..ProblemOptions::default()
            })
            .document
            .expect("gauge-registry mutation"),
        ),
        (
            "joint_noise",
            problem_fixture(ProblemOptions {
                external_noise: true,
                ..ProblemOptions::default()
            })
            .document
            .expect("joint-noise mutation"),
        ),
    ];
    for (field, variant) in variants {
        assert_ne!(
            baseline,
            unresolved_problem_identity(&variant),
            "problem semantic field {field} did not move identity",
        );
    }

    let joint_prior_a = problem_fixture(ProblemOptions {
        joint_prior_choice: 1,
        ..ProblemOptions::default()
    })
    .document
    .expect("joint-prior A problem");
    let joint_prior_b = problem_fixture(ProblemOptions {
        joint_prior_choice: 2,
        ..ProblemOptions::default()
    })
    .document
    .expect("joint-prior B problem");
    assert_eq!(
        joint_prior_a.sources().len(),
        joint_prior_b.sources().len(),
        "joint-prior alternatives must have equal-sized closed source registries",
    );
    assert!(
        joint_prior_a
            .sources()
            .contains_key(&source_key("joint-prior-measure-a"))
            && !joint_prior_a
                .sources()
                .contains_key(&source_key("joint-prior-measure-b"))
            && joint_prior_b
                .sources()
                .contains_key(&source_key("joint-prior-measure-b"))
            && !joint_prior_b
                .sources()
                .contains_key(&source_key("joint-prior-measure-a")),
        "source reachability must retain exactly the selected joint measure",
    );
    assert_ne!(
        unresolved_problem_identity(&joint_prior_a),
        unresolved_problem_identity(&joint_prior_b),
        "joint_prior semantic field did not move identity",
    );

    let generated_composition = problem_fixture(ProblemOptions {
        overlapping_gauges: true,
        declared_gauge_composition: true,
        ..ProblemOptions::default()
    })
    .document
    .expect("generated gauge composition");
    let independent_composition = problem_fixture(ProblemOptions {
        overlapping_gauges: true,
        declared_gauge_composition: true,
        independent_gauge_composition: true,
        ..ProblemOptions::default()
    })
    .document
    .expect("independent-product gauge composition");
    assert_eq!(
        generated_composition.sources(),
        independent_composition.sources()
    );
    assert_eq!(
        generated_composition.gauges(),
        independent_composition.gauges()
    );
    for document in [&generated_composition, &independent_composition] {
        assert_eq!(
            IdentifiabilityProblemDocument::from_canonical_bytes(
                &document.canonical_bytes().expect("composition transport")
            )
            .expect("composition replay"),
            *document,
        );
    }
    assert_ne!(
        unresolved_problem_identity(&generated_composition),
        unresolved_problem_identity(&independent_composition),
        "gauge_compositions semantic field did not move identity",
    );

    let shared_left = problem_fixture(ProblemOptions {
        retrospective_reuse: true,
        declared_sharing: true,
        ..ProblemOptions::default()
    })
    .document
    .expect("shared-data baseline");
    let shared_right = problem_fixture(ProblemOptions {
        retrospective_reuse: true,
        declared_sharing: true,
        alternate_sharing_justification: true,
        ..ProblemOptions::default()
    })
    .document
    .expect("shared-data mutation");
    assert_ne!(
        unresolved_problem_identity(&shared_left),
        unresolved_problem_identity(&shared_right),
        "data_reuse semantic field did not move identity",
    );
    assert_eq!(baseline_document.schema_version(), 3);
    log(
        "problem-identity-semantic-fields",
        "pass",
        "every direct problem field has independent or validity-coupled mutation evidence",
    );
}

#[test]
fn identifiability_source_admission_identity_bindings_have_exact_mutation_evidence() {
    let content_only =
        admit_fixture_with_authority(problem_fixture(ProblemOptions::default()), false);
    let external_trust = admit_fixture_with_single_external_authority(
        problem_fixture(ProblemOptions::default()),
        "forward-a",
    );
    assert_eq!(content_only.id(), external_trust.id());
    assert_ne!(
        content_only.source_admission_id(),
        external_trust.source_admission_id(),
    );
    assert_ne!(
        content_only
            .source_admission_canonical_bytes()
            .expect("content-only source-admission bytes"),
        external_trust
            .source_admission_canonical_bytes()
            .expect("external-trust source-admission bytes"),
    );
    let authority_deltas = content_only
        .source_resolutions()
        .iter()
        .filter_map(|(key, resolution)| {
            (external_trust.source_resolutions().get(key) != Some(resolution)).then_some(key)
        })
        .collect::<Vec<_>>();
    assert_eq!(
        authority_deltas
            .iter()
            .map(|key| key.as_str())
            .collect::<Vec<_>>(),
        vec!["forward-a"],
        "one exact authority disposition must move the source-admission identity",
    );
    let different_problem = admit_fixture(problem_fixture(ProblemOptions {
        second_case_complementary: true,
        ..ProblemOptions::default()
    }));
    assert_ne!(content_only.id(), different_problem.id());
    assert_eq!(
        content_only.source_resolutions(),
        different_problem.source_resolutions(),
        "problem-only mutation must retain the exact resolution registry",
    );
    assert_ne!(
        content_only.source_admission_id(),
        different_problem.source_admission_id(),
        "problem_id must move SourceAdmissionId independently of authority disposition",
    );
    assert_eq!(
        content_only.source_admission_id().digest().as_bytes().len(),
        32
    );
    log(
        "source-admission-identity-semantic-fields",
        "pass",
        "problem id and exact resolution authority independently move SourceAdmissionId",
    );
}

#[test]
fn source_admission_id_is_stable_across_execution_variants() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let source_admission_id = problem.source_admission_id();
    let si = execution(&problem, false, 17, 1.0e-10, false).expect("SI execution");
    let affine = execution(&problem, true, 18, 1.0e-8, false).expect("affine execution");
    assert_ne!(si.id().expect("SI id"), affine.id().expect("affine id"));
    assert_eq!(problem.source_admission_id(), source_admission_id);
    log(
        "source-admission-execution-noninterference",
        "pass",
        "execution coordinates, seeds, and tolerances cannot rewrite source admission",
    );
}

#[test]
fn coordinates_do_not_move_problem_identity() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let identity = problem.id();
    let si = execution(&problem, false, 17, 1.0e-10, false).expect("SI plan");
    let affine = execution(&problem, true, 17, 1.0e-10, false).expect("affine plan");
    assert_eq!(problem.id(), identity);
    assert_ne!(si.id().expect("SI id"), affine.id().expect("affine id"));
    log(
        "coordinate-noninterference",
        "pass",
        "coordinates move execution only, never ProblemId",
    );
}

#[test]
fn seed_and_tolerance_move_execution_identity() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let baseline = execution(&problem, false, 17, 1.0e-10, false).expect("baseline");
    let seed = execution(&problem, false, 18, 1.0e-10, false).expect("seed variant");
    let tolerance = execution(&problem, false, 17, 1.0e-8, false).expect("tol variant");
    assert_ne!(baseline.id().expect("id"), seed.id().expect("id"));
    assert_ne!(baseline.id().expect("id"), tolerance.id().expect("id"));
    log(
        "execution-semantic-fields",
        "pass",
        "Five Explicits and numerical policy are execution semantics",
    );
}

#[test]
fn execution_source_authority_moves_execution_identity() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let content_only = execution_with_claim_requests_and_authority(
        &problem,
        false,
        17,
        1.0e-10,
        false,
        vec![default_claim_request(&problem)],
        Vec::new(),
        false,
    )
    .expect("content-verified execution");
    let external_trust = execution_with_claim_requests_and_authority(
        &problem,
        false,
        17,
        1.0e-10,
        false,
        vec![default_claim_request(&problem)],
        Vec::new(),
        true,
    )
    .expect("externally trusted execution");
    assert_ne!(
        content_only.id().expect("content-only id"),
        external_trust.id().expect("external-trust id"),
    );
    log(
        "execution-source-authority-identity",
        "pass",
        "execution authority is transitive identity state, not an unverified annotation",
    );
}

#[test]
fn identifiability_execution_identity_bindings_have_exact_mutation_evidence() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let baseline = execution(&problem, false, 17, 1.0e-10, false).expect("baseline");
    let baseline_id = baseline.id().expect("baseline id");
    let baseline_parts = ExecutionParts::from_plan(&baseline);
    let assert_moves = |field: &str, parts: ExecutionParts, external_trust: bool| {
        let variant = parts.build(&problem, external_trust);
        assert_ne!(
            baseline_id,
            variant.id().expect("variant execution id"),
            "execution semantic field {field} did not move identity",
        );
    };

    for (field, header) in [
        (
            "header.units",
            execution_header_with_semantics(&["K", "Pa"], 17, 1.0e-9, 30_000, 32 << 20, "1", false),
        ),
        (
            "header.seed",
            execution_header_with_semantics(&["Pa"], 18, 1.0e-9, 30_000, 32 << 20, "1", false),
        ),
        (
            "header.accuracy",
            execution_header_with_semantics(&["Pa"], 17, 1.0e-8, 30_000, 32 << 20, "1", false),
        ),
        (
            "header.time_ms",
            execution_header_with_semantics(&["Pa"], 17, 1.0e-9, 30_001, 32 << 20, "1", false),
        ),
        (
            "header.memory_bytes",
            execution_header_with_semantics(
                &["Pa"],
                17,
                1.0e-9,
                30_000,
                (32 << 20) + 1,
                "1",
                false,
            ),
        ),
        (
            "header.versions",
            execution_header_with_semantics(&["Pa"], 17, 1.0e-9, 30_000, 32 << 20, "2", false),
        ),
        (
            "header.capabilities",
            execution_header_with_semantics(&["Pa"], 17, 1.0e-9, 30_000, 32 << 20, "1", true),
        ),
    ] {
        let mut parts = baseline_parts.clone();
        parts.header = header;
        assert_moves(field, parts, false);
    }

    let mut parts = baseline_parts.clone();
    parts.analyzer = source("analyzer-v2", SourceKind::Analyzer, hash("analyzer-v2"));
    assert_moves("analyzer", parts, false);
    let mut parts = baseline_parts.clone();
    parts.build = source("build-v2", SourceKind::Build, hash("build-v2"));
    assert_moves("build", parts, false);
    let mut parts = baseline_parts.clone();
    parts.derivative_provider = None;
    assert_moves("derivative_provider", parts, false);
    let mut parts = baseline_parts.clone();
    let baseline_request = parts
        .claim_requests
        .first()
        .expect("baseline claim request")
        .clone();
    parts.claim_requests[0] = ClaimRequest::new(
        baseline_request.claim().clone(),
        DimensionlessErrorPolicy::try_new(
            baseline_request.error_policy().metric().clone(),
            baseline_request
                .error_policy()
                .nondimensionalization()
                .clone(),
            9.0e-9,
        )
        .expect("claim error-bound mutation"),
    );
    assert_moves("claim_requests", parts, false);
    let mut parts = baseline_parts.clone();
    for (role, action) in &mut parts.actions {
        if role.as_str() == "yield_stress" {
            *action = ParameterExecutionAction::Optimize {
                coordinate: coordinate("yield_stress", true),
            };
        }
    }
    assert_moves("actions", parts, false);
    let mut parts = baseline_parts.clone();
    let numerical_nondimensionalization = parts.numerical.nondimensionalization().clone();
    parts.numerical = IdentifiabilityNumericalPolicy::try_new(
        1.0e-8,
        0.0,
        1.0e12,
        ArithmeticPolicy::CertifiedInterval,
        numerical_nondimensionalization,
    )
    .expect("numerical mutation");
    assert_moves("numerical", parts, false);
    let mut parts = baseline_parts.clone();
    parts.initialization = source(
        "initialization-v2",
        SourceKind::Assumption,
        hash("initialization-v2"),
    );
    assert_moves("initialization", parts, false);
    let mut parts = baseline_parts.clone();
    parts.stopping = source("stopping-v2", SourceKind::Assumption, hash("stopping-v2"));
    assert_moves("stopping", parts, false);
    let mut parts = baseline_parts.clone();
    parts.determinism = source(
        "determinism-v2",
        SourceKind::Assumption,
        hash("determinism-v2"),
    );
    assert_moves("determinism_contract", parts, false);
    assert_moves("source_authority", baseline_parts.clone(), true);

    let gauge_problem = admit_fixture(problem_fixture(ProblemOptions {
        gauge_case: 4,
        ..ProblemOptions::default()
    }));
    let action =
        GaugeActionReference::Single(GaugeClassId::try_new("fixture-gauge").expect("gauge id"));
    let gauge_claim = structural_claim(
        &gauge_problem,
        "identity-unreduced-gauge",
        FiberStructure::OrbitQuotientUnique {
            action: action.clone(),
        },
        ClaimSubject::Parameter(role("yield_stress")),
        ClaimScope::WholeCampaign,
    );
    let unreduced = |reason: &str| {
        GaugeReductionBinding::try_new(
            GaugeReductionId::try_new("identity-unreduced-plan").expect("reduction id"),
            action.clone(),
            BTreeSet::from([gauge_claim.id().clone()]),
            GaugeReductionPlan::Unreduced {
                reason: reason.to_string(),
            },
            GaugeReductionStage::Root,
            GaugeMeasureSemantics::NotApplicable {
                reason: "no quotient transport is claimed for an unreduced action".to_string(),
            },
        )
        .expect("unreduced gauge identity fixture")
    };
    let unreduced_left = execution_for_claim_with_gauge_reductions(
        &gauge_problem,
        gauge_claim.clone(),
        vec![unreduced(
            "the action remains explicit in the original coordinates",
        )],
    )
    .expect("left unreduced execution");
    let unreduced_right = execution_for_claim_with_gauge_reductions(
        &gauge_problem,
        gauge_claim.clone(),
        vec![unreduced(
            "the action is intentionally retained for a different declared reason",
        )],
    )
    .expect("right unreduced execution");
    assert_eq!(
        unreduced_left.source_authority(),
        unreduced_right.source_authority(),
        "reason-only reduction mutation must not alter source authority",
    );
    assert_ne!(
        unreduced_left.id().expect("left unreduced id"),
        unreduced_right.id().expect("right unreduced id"),
        "gauge_reductions semantic field did not move identity",
    );

    let physical_variant = admit_fixture(problem_fixture(ProblemOptions {
        second_case_complementary: true,
        ..ProblemOptions::default()
    }));
    let physical_execution = baseline_parts.clone().build(&physical_variant, false);
    assert_ne!(baseline.problem_id(), physical_execution.problem_id());
    assert_ne!(
        baseline_id,
        physical_execution.id().expect("physical variant id")
    );

    let authority_variant =
        admit_fixture_with_authority(problem_fixture(ProblemOptions::default()), true);
    let authority_execution = baseline_parts.build(&authority_variant, false);
    assert_eq!(baseline.problem_id(), authority_execution.problem_id());
    assert_ne!(
        baseline.source_admission_id(),
        authority_execution.source_admission_id(),
    );
    assert_ne!(
        baseline_id,
        authority_execution.id().expect("authority variant id")
    );
    assert_eq!(baseline.schema_version(), 3);
    log(
        "execution-identity-semantic-fields",
        "pass",
        "every execution field and header projection has direct mutation evidence",
    );
}

#[test]
fn execution_action_input_order_is_nonsemantic() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let baseline = execution(&problem, false, 17, 1.0e-10, false).expect("baseline");
    let mut reversed_actions = baseline
        .actions()
        .iter()
        .map(|(role, action)| (role.clone(), action.clone()))
        .collect::<Vec<_>>();
    reversed_actions.reverse();
    let rebuilt = IdentifiabilityExecutionPlan::try_new(
        baseline.header().clone(),
        &problem,
        baseline.analyzer().clone(),
        baseline.build().clone(),
        baseline.derivative_provider().cloned(),
        baseline.claim_requests().values().cloned().collect(),
        reversed_actions,
        baseline.gauge_reductions().values().cloned().collect(),
        baseline.numerical_policy().clone(),
        baseline.initialization().clone(),
        baseline.stopping().clone(),
        baseline.determinism_contract().clone(),
        baseline.source_authority().clone(),
    )
    .expect("reordered execution");
    assert_eq!(
        baseline.id().expect("baseline id"),
        rebuilt.id().expect("rebuilt id")
    );
    assert_eq!(
        baseline.canonical_bytes().expect("baseline transport"),
        rebuilt.canonical_bytes().expect("rebuilt transport"),
    );
    log(
        "execution-action-order",
        "pass",
        "parameter actions are keyed semantics rather than caller sequence semantics",
    );
}

#[test]
fn execution_action_must_match_physical_treatment() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let result = execution(&problem, false, 17, 1.0e-10, true);
    assert!(matches!(
        result,
        Err(IdentifiabilityError::InvalidText {
            field: "execution parameter treatment",
            ..
        })
    ));
    log(
        "treatment-action-coverage",
        "pass",
        "marginalized parameter cannot silently become optimized",
    );
}

#[test]
fn execution_roundtrip_requires_the_exact_admitted_problem() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let plan = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    let bytes = plan.canonical_bytes().expect("execution bytes");
    assert!(matches!(
        IdentifiabilityExecutionPlan::from_canonical_bytes(
            &bytes,
            &problem,
            &SourceResolutionSet::default(),
        ),
        Err(IdentifiabilityError::SourceMismatch {
            field: "execution source-resolution replay",
        })
    ));
    let decoded = IdentifiabilityExecutionPlan::from_canonical_bytes(
        &bytes,
        &problem,
        plan.source_authority(),
    )
    .expect("execution roundtrip");
    assert_eq!(decoded, plan);
    log(
        "execution-roundtrip",
        "pass",
        "transport revalidates ProblemId and SourceAdmissionId",
    );
}

#[test]
fn artifact_labels_do_not_move_execution_or_assessment_identity() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let plan = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    let execution_id = plan.id().expect("execution id");
    let mut execution_bytes = plan.canonical_bytes().expect("execution transport");
    let execution_label = b"execution-1";
    let execution_label_at = execution_bytes
        .windows(execution_label.len())
        .position(|window| window == execution_label)
        .expect("execution artifact label in exact transport");
    execution_bytes[execution_label_at..execution_label_at + execution_label.len()]
        .copy_from_slice(b"execution-2");
    let relabeled_execution = IdentifiabilityExecutionPlan::from_canonical_bytes(
        &execution_bytes,
        &problem,
        plan.source_authority(),
    )
    .expect("relabeled execution transport");
    assert_eq!(
        relabeled_execution.id().expect("relabeled id"),
        execution_id
    );
    assert_ne!(
        relabeled_execution
            .canonical_bytes()
            .expect("exact transport"),
        plan.canonical_bytes().expect("baseline exact transport"),
    );

    let assessment = assessment(&problem, &plan, "receipt");
    let assessment_id = assessment.id().expect("assessment id");
    let mut assessment_bytes = assessment.canonical_bytes().expect("assessment transport");
    let assessment_label = b"assessment-1";
    let assessment_label_at = assessment_bytes
        .windows(assessment_label.len())
        .position(|window| window == assessment_label)
        .expect("assessment artifact label in exact transport");
    assessment_bytes[assessment_label_at..assessment_label_at + assessment_label.len()]
        .copy_from_slice(b"assessment-2");
    let relabeled_assessment = IdentifiabilityAssessment::from_canonical_bytes(
        &assessment_bytes,
        &problem,
        &plan,
        assessment.source_authority(),
    )
    .expect("relabeled assessment transport");
    assert_eq!(
        relabeled_assessment.id().expect("relabeled id"),
        assessment_id,
    );
    log(
        "artifact-label-nonsemantic",
        "pass",
        "exact transport retains ledger labels while scientific identities exclude them",
    );
}

#[test]
fn evidence_changes_assessment_not_problem_or_execution() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let plan = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    let problem_id = problem.id();
    let execution_id = plan.id().expect("execution id");
    let left = assessment(&problem, &plan, "receipt-left");
    let right = assessment(&problem, &plan, "receipt-right");
    assert_ne!(left.id().expect("left id"), right.id().expect("right id"));
    assert_eq!(problem.id(), problem_id);
    assert_eq!(plan.id().expect("execution id"), execution_id);
    log(
        "assessment-noninterference",
        "pass",
        "evidence cannot rewrite problem or execution preimages",
    );
}

#[test]
fn identifiability_assessment_identity_bindings_have_exact_mutation_evidence() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let plan = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    let baseline = assessment(&problem, &plan, "receipt-left");
    let baseline_id = baseline.id().expect("baseline assessment id");
    let baseline_parts = AssessmentParts::from_assessment(&baseline);
    let assert_moves = |field: &str, parts: AssessmentParts| {
        let variant = parts.build(&problem, &plan);
        assert_eq!(baseline.problem_id(), variant.problem_id());
        assert_eq!(baseline.execution_id(), variant.execution_id());
        assert_ne!(
            baseline_id,
            variant.id().expect("variant assessment id"),
            "assessment semantic field {field} did not move identity",
        );
    };

    const ASSESSMENT_SEED: u64 = 0x1d3_171f;
    for (field, header) in [
        (
            "header.units",
            assessment_header_with_semantics(
                &["K", "Pa"],
                ASSESSMENT_SEED,
                1.0e-9,
                30_000,
                32 << 20,
                "1",
                false,
            ),
        ),
        (
            "header.seed",
            assessment_header_with_semantics(
                &["Pa"],
                ASSESSMENT_SEED + 1,
                1.0e-9,
                30_000,
                32 << 20,
                "1",
                false,
            ),
        ),
        (
            "header.accuracy",
            assessment_header_with_semantics(
                &["Pa"],
                ASSESSMENT_SEED,
                1.0e-8,
                30_000,
                32 << 20,
                "1",
                false,
            ),
        ),
        (
            "header.time_ms",
            assessment_header_with_semantics(
                &["Pa"],
                ASSESSMENT_SEED,
                1.0e-9,
                30_001,
                32 << 20,
                "1",
                false,
            ),
        ),
        (
            "header.memory_bytes",
            assessment_header_with_semantics(
                &["Pa"],
                ASSESSMENT_SEED,
                1.0e-9,
                30_000,
                (32 << 20) + 1,
                "1",
                false,
            ),
        ),
        (
            "header.versions",
            assessment_header_with_semantics(
                &["Pa"],
                ASSESSMENT_SEED,
                1.0e-9,
                30_000,
                32 << 20,
                "2",
                false,
            ),
        ),
        (
            "header.capabilities",
            assessment_header_with_semantics(
                &["Pa"],
                ASSESSMENT_SEED,
                1.0e-9,
                30_000,
                32 << 20,
                "1",
                true,
            ),
        ),
    ] {
        let mut parts = baseline_parts.clone();
        parts.header = header;
        assert_eq!(parts.claims, baseline_parts.claims);
        assert_eq!(parts.evidence, baseline_parts.evidence);
        assert_eq!(parts.source_authority, baseline_parts.source_authority);
        assert_moves(field, parts);
    }

    let baseline_claim = baseline.claims().values().next().expect("baseline claim");
    let claim_variant = TypedIdentifiabilityClaim::new(
        baseline_claim.id().clone(),
        baseline_claim.information().clone(),
        baseline_claim.extent(),
        FiberStructure::FiniteToOne {
            maximum_cardinality: Some(FiberCardinalityBound::UniformU64(2)),
        },
        baseline_claim.quantifier().clone(),
        baseline_claim.scalar_domain().clone(),
        baseline_claim.subject().clone(),
        baseline_claim.scope().clone(),
    );
    let variant_plan =
        execution_for_claim(&problem, false, 17, 1.0e-10, false, claim_variant.clone())
            .expect("claim-variant execution");
    let mut claim_parts = baseline_parts.clone();
    *claim_parts
        .claims
        .iter_mut()
        .find(|claim| claim.id() == baseline_claim.id())
        .expect("claim mutation target") = claim_variant;
    assert_eq!(claim_parts.evidence, baseline_parts.evidence);
    assert_eq!(
        claim_parts.source_authority,
        baseline_parts.source_authority
    );
    assert_eq!(claim_parts.header, baseline_parts.header);
    let claim_variant_assessment = claim_parts.build(&problem, &variant_plan);
    assert_ne!(
        baseline.execution_id(),
        claim_variant_assessment.execution_id()
    );
    assert_ne!(
        baseline_id,
        claim_variant_assessment
            .id()
            .expect("claim-variant assessment id"),
        "validity-coupled claim and execution projections must move assessment identity",
    );

    let mut parts = baseline_parts.clone();
    let (_, conclusion) = parts
        .evidence
        .iter_mut()
        .find(|(id, _)| id == baseline_claim.id())
        .expect("evidence mutation target");
    let (method, receipt, metric, nondimensionalization) = match conclusion.clone() {
        ClaimAssessment::ClaimedEstablished {
            method,
            receipt,
            metric,
            nondimensionalization,
            ..
        } => (method, receipt, metric, nondimensionalization),
        _ => panic!("baseline evidence unexpectedly changed variant"),
    };
    *conclusion = ClaimAssessment::ClaimedEstablished {
        method,
        receipt,
        metric,
        nondimensionalization,
        certified_error_bound: 7.5e-9,
        gauge_resolutions: BTreeMap::new(),
    };
    assert_eq!(parts.claims, baseline_parts.claims);
    assert_eq!(parts.source_authority, baseline_parts.source_authority);
    assert_eq!(parts.header, baseline_parts.header);
    assert_moves("evidence", parts);

    let gauge_problem = admit_fixture(problem_fixture(ProblemOptions {
        one_gauge: true,
        ..ProblemOptions::default()
    }));
    let derived_definition = source(
        "gauge-invariant-functional",
        SourceKind::DerivedFunctional,
        hash("gauge-invariant-functional"),
    );
    let gauge_claim = structural_claim(
        &gauge_problem,
        "gauge-resolution-identity",
        FiberStructure::Unique,
        ClaimSubject::DerivedFunctional {
            definition: derived_definition,
            parameters: BTreeSet::from([role("yield_stress")]),
        },
        ClaimScope::WholeCampaign,
    );
    let gauge_execution = execution_for_claim(
        &gauge_problem,
        false,
        17,
        1.0e-10,
        false,
        gauge_claim.clone(),
    )
    .expect("gauge-resolution execution");
    let build_gauge_assessment = |disposition| {
        let request = &gauge_execution.claim_requests()[gauge_claim.id()];
        let method = gauge_execution.analyzer().clone();
        let receipt = source(
            "gauge-resolution-receipt",
            SourceKind::EvidenceReceipt,
            hash("gauge-resolution-receipt"),
        );
        let action =
            GaugeActionReference::Single(GaugeClassId::try_new("single-gauge").expect("gauge id"));
        let mut assessment_sources = claim_sources(&gauge_claim);
        assessment_sources.extend([
            method.clone(),
            receipt.clone(),
            request.error_policy().metric().clone(),
            request.error_policy().nondimensionalization().clone(),
        ]);
        IdentifiabilityAssessment::try_new(
            header("gauge-resolution-assessment", "identifiability.assess"),
            &gauge_problem,
            &gauge_execution,
            vec![gauge_claim.clone()],
            vec![(
                gauge_claim.id().clone(),
                ClaimAssessment::ClaimedEstablished {
                    method: method.clone(),
                    receipt: receipt.clone(),
                    metric: request.error_policy().metric().clone(),
                    nondimensionalization: request.error_policy().nondimensionalization().clone(),
                    certified_error_bound: 5.0e-9,
                    gauge_resolutions: BTreeMap::from([(
                        action.clone(),
                        GaugeResolutionEvidence::new(action, disposition, method, receipt),
                    )]),
                },
            )],
            resolve_owned_sources(assessment_sources, false, None),
        )
        .expect("positive nonempty gauge-resolution assessment")
    };
    let no_projection = build_gauge_assessment(GaugeResolutionDisposition::NoProjectionOnSubject);
    let descends = build_gauge_assessment(GaugeResolutionDisposition::SubjectDescendsToQuotient);
    assert_eq!(no_projection.problem_id(), descends.problem_id());
    assert_eq!(no_projection.execution_id(), descends.execution_id());
    assert_eq!(no_projection.claims(), descends.claims());
    assert_eq!(
        no_projection.source_authority(),
        descends.source_authority(),
    );
    assert_ne!(
        no_projection.id().expect("no-projection assessment id"),
        descends.id().expect("quotient-descent assessment id"),
        "nonempty gauge_resolutions semantic field did not move identity",
    );

    let receipt = match baseline
        .evidence()
        .values()
        .next()
        .expect("baseline evidence")
    {
        ClaimAssessment::ClaimedEstablished { receipt, .. } => receipt.clone(),
        _ => panic!("baseline evidence unexpectedly changed variant"),
    };
    let mut authority_entries = baseline
        .source_authority()
        .entries()
        .values()
        .cloned()
        .collect::<Vec<_>>();
    let replacement = SourceResolution::verify(
        &receipt,
        b"receipt-left",
        external_trust("assessment-receipt-external-trust", &receipt),
    )
    .expect("assessment-exclusive trusted receipt resolution");
    let mut replacement_count = 0;
    for resolution in &mut authority_entries {
        if resolution.key() == receipt.key() {
            *resolution = replacement.clone();
            replacement_count += 1;
        }
    }
    assert_eq!(
        replacement_count, 1,
        "one assessment-exclusive authority moved"
    );
    let mut parts = baseline_parts.clone();
    parts.source_authority =
        SourceResolutionSet::try_new(authority_entries).expect("authority mutation");
    assert_eq!(parts.claims, baseline_parts.claims);
    assert_eq!(parts.evidence, baseline_parts.evidence);
    assert_eq!(parts.header, baseline_parts.header);
    assert_moves("source_authority", parts);

    let execution_variant =
        execution(&problem, false, 18, 1.0e-10, false).expect("execution-id variant");
    let execution_assessment = assessment(&problem, &execution_variant, "receipt-left");
    assert_eq!(baseline.problem_id(), execution_assessment.problem_id());
    assert_ne!(baseline.execution_id(), execution_assessment.execution_id());
    assert_eq!(baseline.header(), execution_assessment.header());
    assert_eq!(baseline.claims(), execution_assessment.claims());
    assert_eq!(baseline.evidence(), execution_assessment.evidence());
    assert_eq!(
        baseline.source_authority(),
        execution_assessment.source_authority(),
    );
    assert_ne!(
        baseline_id,
        execution_assessment
            .id()
            .expect("execution variant assessment id"),
    );

    let problem_variant = admit_fixture(problem_fixture(ProblemOptions {
        second_case_complementary: true,
        ..ProblemOptions::default()
    }));
    let variant_plan =
        execution(&problem_variant, false, 17, 1.0e-10, false).expect("problem variant plan");
    let problem_assessment = assessment(&problem_variant, &variant_plan, "receipt-left");
    assert_ne!(baseline.problem_id(), problem_assessment.problem_id());
    assert_ne!(baseline.execution_id(), problem_assessment.execution_id());
    assert_eq!(baseline.header(), problem_assessment.header());
    assert_eq!(baseline.claims(), problem_assessment.claims());
    assert_eq!(baseline.evidence(), problem_assessment.evidence());
    assert_eq!(
        baseline.source_authority(),
        problem_assessment.source_authority(),
    );
    assert_ne!(
        baseline_id,
        problem_assessment
            .id()
            .expect("problem variant assessment id"),
    );

    let authority_problem = admit_fixture_with_single_external_authority(
        problem_fixture(ProblemOptions::default()),
        "forward-a",
    );
    assert_eq!(problem.id(), authority_problem.id());
    assert_ne!(
        problem.source_admission_id(),
        authority_problem.source_admission_id(),
    );
    let authority_plan =
        execution(&authority_problem, false, 17, 1.0e-10, false).expect("authority plan");
    let authority_assessment = assessment(&authority_problem, &authority_plan, "receipt-left");
    assert_ne!(baseline.execution_id(), authority_assessment.execution_id());
    assert_eq!(baseline.header(), authority_assessment.header());
    assert_eq!(baseline.claims(), authority_assessment.claims());
    assert_eq!(baseline.evidence(), authority_assessment.evidence());
    assert_eq!(
        baseline.source_authority(),
        authority_assessment.source_authority(),
    );
    assert_ne!(
        baseline_id,
        authority_assessment
            .id()
            .expect("authority-envelope assessment id"),
    );

    const ASSESSMENT_MAGIC: &[u8] = b"fs-material-identifiability-assessment\0";
    let mut stale = baseline.canonical_bytes().expect("assessment transport");
    stale[ASSESSMENT_MAGIC.len()..ASSESSMENT_MAGIC.len() + 4].copy_from_slice(&1_u32.to_le_bytes());
    assert!(matches!(
        IdentifiabilityAssessment::from_canonical_bytes(
            &stale,
            &problem,
            &plan,
            baseline.source_authority(),
        ),
        Err(IdentifiabilityError::UnsupportedSchemaVersion { .. })
    ));
    log(
        "assessment-identity-semantic-fields",
        "pass",
        "every direct field and each validity-coupled parent projection has mutation evidence",
    );
}

#[test]
fn assessment_authority_must_agree_with_problem_and_execution_on_transitive_overlap() {
    let problem = admit_fixture(problem_fixture(ProblemOptions {
        claim_strata_in_problem: true,
        ..ProblemOptions::default()
    }));
    let execution = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    assert!(matches!(
        execution
            .claim_requests()
            .values()
            .next()
            .expect("default claim request")
            .claim()
            .fiber(),
        FiberStructure::Stratified { .. }
    ));
    let overlap_key = source_key("claim-strata");
    assert_eq!(
        execution.source_authority().entries().get(&overlap_key),
        problem.source_resolutions().get(&overlap_key),
        "preregistering the claim must carry forward the problem's admitted strata authority",
    );
    assessment_result(&problem, &execution, "matching-domain-authority")
        .expect("matching problem/assessment authority admits");
    let result = assessment_result_with_claim_source_authority(
        &problem,
        &execution,
        "conflicting-domain-authority",
        external_trust(
            "conflicting-domain-trust",
            &problem.document().sources()[&overlap_key],
        ),
    );
    assert!(matches!(
        result,
        Err(IdentifiabilityError::SourceMismatch {
            field: "assessment/execution source authority",
        })
    ));
    log(
        "assessment-problem-execution-authority-overlap",
        "pass",
        "assessment cannot relabel claim strata whose problem authority was carried into exact execution preregistration",
    );
}

#[test]
fn assessment_authority_must_agree_with_execution_on_transitive_overlap() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let execution = execution_with_claim_requests_and_authority(
        &problem,
        false,
        17,
        1.0e-10,
        false,
        vec![default_claim_request(&problem)],
        Vec::new(),
        true,
    )
    .expect("externally trusted execution");
    let result = assessment_result(&problem, &execution, "execution-overlap-receipt");
    assert!(matches!(
        result,
        Err(IdentifiabilityError::SourceMismatch {
            field: "assessment/execution source authority",
        })
    ));
    log(
        "assessment-execution-authority-overlap",
        "pass",
        "assessment cannot relabel authority for its execution analyzer source",
    );
}

#[test]
fn assessment_input_order_is_nonsemantic() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let execution = execution_with_claim_requests_and_authority(
        &problem,
        false,
        17,
        1.0e-10,
        false,
        two_claims().into_iter().map(request_for_claim).collect(),
        Vec::new(),
        false,
    )
    .expect("two-claim execution");
    let baseline = two_claim_assessment(&problem, &execution);
    let mut claims = baseline.claims().values().cloned().collect::<Vec<_>>();
    let mut evidence = baseline
        .evidence()
        .iter()
        .map(|(id, value)| (id.clone(), value.clone()))
        .collect::<Vec<_>>();
    let mut resolutions = baseline
        .source_authority()
        .entries()
        .values()
        .cloned()
        .collect::<Vec<_>>();
    assert!(claims.len() >= 2, "claim order test needs multiple claims");
    assert!(
        evidence.len() >= 2,
        "evidence order test needs multiple entries"
    );
    assert!(
        resolutions.len() >= 2,
        "authority order test needs multiple resolutions"
    );
    let baseline_claims = claims.clone();
    let baseline_evidence = evidence.clone();
    let baseline_resolutions = resolutions.clone();
    let assert_same = |field: &str,
                       claims: Vec<TypedIdentifiabilityClaim>,
                       evidence: Vec<(ClaimId, ClaimAssessment)>,
                       resolutions: Vec<SourceResolution>| {
        let rebuilt = IdentifiabilityAssessment::try_new(
            baseline.header().clone(),
            &problem,
            &execution,
            claims,
            evidence,
            SourceResolutionSet::try_new(resolutions).expect("reordered assessment authority"),
        )
        .unwrap_or_else(|error| panic!("reordered {field} refused: {error}"));
        assert_eq!(
            baseline.id().expect("baseline assessment id"),
            rebuilt.id().expect("rebuilt assessment id"),
            "caller order for {field} moved assessment identity",
        );
        assert_eq!(
            baseline.canonical_bytes().expect("baseline transport"),
            rebuilt.canonical_bytes().expect("rebuilt transport"),
            "caller order for {field} moved canonical transport",
        );
    };

    claims.reverse();
    assert_ne!(
        claims, baseline_claims,
        "claim reversal must be non-vacuous"
    );
    assert_same(
        "claims",
        claims,
        baseline_evidence.clone(),
        baseline_resolutions.clone(),
    );
    evidence.reverse();
    assert_ne!(
        evidence, baseline_evidence,
        "evidence reversal must be non-vacuous"
    );
    assert_same(
        "evidence",
        baseline_claims.clone(),
        evidence,
        baseline_resolutions.clone(),
    );
    resolutions.reverse();
    assert_ne!(
        resolutions, baseline_resolutions,
        "authority reversal must be non-vacuous"
    );
    assert_same(
        "source authority",
        baseline_claims,
        baseline_evidence,
        resolutions,
    );
    log(
        "assessment-input-order",
        "pass",
        "claims, evidence, and authority are canonical keyed collections",
    );
}

#[test]
fn assessment_roundtrip_preserves_product_typed_claim() {
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let plan = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    let assessment = assessment(&problem, &plan, "receipt");
    let bytes = assessment.canonical_bytes().expect("assessment bytes");
    let decoded = IdentifiabilityAssessment::from_canonical_bytes(
        &bytes,
        &problem,
        &plan,
        assessment.source_authority(),
    )
    .expect("assessment roundtrip");
    assert_eq!(decoded, assessment);
    log(
        "assessment-roundtrip",
        "pass",
        "regime, extent, quantifier, scalar domain, subject, and scope retained",
    );
}

#[test]
fn identity_domains_and_wire_magics_are_stage_separated() {
    const PROBLEM_MAGIC: &[u8] = b"fs-material-identifiability-problem\0";
    const SOURCE_ADMISSION_MAGIC: &[u8] = b"fs-material-identifiability-source-admission\0";
    const EXECUTION_MAGIC: &[u8] = b"fs-material-identifiability-execution\0";
    const ASSESSMENT_MAGIC: &[u8] = b"fs-material-identifiability-assessment\0";

    let domains = [
        IDENTIFIABILITY_PROBLEM_IDENTITY_DOMAIN,
        IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_DOMAIN,
        IDENTIFIABILITY_EXECUTION_IDENTITY_DOMAIN,
        IDENTIFIABILITY_ASSESSMENT_IDENTITY_DOMAIN,
    ];
    assert_eq!(domains.into_iter().collect::<BTreeSet<_>>().len(), 4);
    assert_eq!(
        domains
            .into_iter()
            .map(|domain| hash_domain(domain, b"same-stage-preimage"))
            .collect::<BTreeSet<_>>()
            .len(),
        4,
    );

    let document = problem_fixture(ProblemOptions::default())
        .document
        .expect("problem");
    assert!(
        document
            .canonical_bytes()
            .expect("problem transport")
            .starts_with(PROBLEM_MAGIC)
    );
    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    assert!(
        problem
            .source_admission_canonical_bytes()
            .expect("source admission transport")
            .starts_with(SOURCE_ADMISSION_MAGIC)
    );
    let execution = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    assert!(
        execution
            .canonical_bytes()
            .expect("execution transport")
            .starts_with(EXECUTION_MAGIC)
    );
    let assessment = assessment(&problem, &execution, "domain-receipt");
    assert!(
        assessment
            .canonical_bytes()
            .expect("assessment transport")
            .starts_with(ASSESSMENT_MAGIC)
    );
    log(
        "identity-domain-and-magic-separation",
        "pass",
        "all four authority stages have distinct domains and exact wire magics",
    );
}

#[test]
fn identifiability_identity_preimages_have_exact_wire_layout() {
    const PROBLEM_MAGIC: &[u8] = b"fs-material-identifiability-problem\0";
    const SOURCE_ADMISSION_MAGIC: &[u8] = b"fs-material-identifiability-source-admission\0";
    const EXECUTION_MAGIC: &[u8] = b"fs-material-identifiability-execution\0";
    const ASSESSMENT_MAGIC: &[u8] = b"fs-material-identifiability-assessment\0";

    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let problem_bytes = problem
        .document()
        .canonical_bytes()
        .expect("problem identity preimage");
    assert!(problem_bytes.starts_with(PROBLEM_MAGIC));
    let mut problem_at = PROBLEM_MAGIC.len();
    assert_eq!(
        read_u32_le(&problem_bytes, &mut problem_at, "problem version"),
        IDENTIFIABILITY_PROBLEM_IDENTITY_VERSION,
    );
    for expected in [
        problem.document().context_source(),
        problem.document().material_source(),
        problem.document().model_source(),
        problem.document().graph_source(),
    ] {
        assert_eq!(
            read_text(&problem_bytes, &mut problem_at, "problem root source"),
            expected.as_str(),
            "problem root-binding order moved",
        );
    }
    assert_eq!(
        problem_bytes.get(problem_at),
        Some(&0),
        "absent joint-prior tag moved",
    );
    problem_at += 1;
    assert_eq!(
        read_u32_le(&problem_bytes, &mut problem_at, "source registry count"),
        u32::try_from(problem.document().sources().len()).expect("bounded source count"),
    );
    let known_parameter_bound = 1.0e6_f64.to_bits().to_le_bytes();
    assert!(
        problem_bytes
            .windows(known_parameter_bound.len())
            .any(|window| window == known_parameter_bound),
        "problem numeric fields must use canonical f64 little-endian encoding",
    );
    assert_eq!(
        hash_domain(IDENTIFIABILITY_PROBLEM_IDENTITY_DOMAIN, &problem_bytes),
        problem.id().digest(),
    );

    let posterior_document = problem_fixture(ProblemOptions {
        joint_prior_choice: 1,
        ..ProblemOptions::default()
    })
    .document
    .expect("posterior wire-layout problem");
    let posterior_bytes = posterior_document
        .canonical_bytes()
        .expect("posterior problem preimage");
    let mut posterior_at = PROBLEM_MAGIC.len();
    assert_eq!(
        read_u32_le(&posterior_bytes, &mut posterior_at, "posterior version"),
        IDENTIFIABILITY_PROBLEM_IDENTITY_VERSION,
    );
    for _ in 0..4 {
        let _ = read_text(&posterior_bytes, &mut posterior_at, "posterior root source");
    }
    assert_eq!(
        posterior_bytes.get(posterior_at),
        Some(&1),
        "present joint-prior tag moved",
    );
    posterior_at += 1;
    assert_eq!(
        read_text(&posterior_bytes, &mut posterior_at, "joint prior source"),
        "joint-prior-measure-a",
    );

    let source_admission_bytes = problem
        .source_admission_canonical_bytes()
        .expect("source-admission identity preimage");
    assert!(source_admission_bytes.starts_with(SOURCE_ADMISSION_MAGIC));
    let mut admission_at = SOURCE_ADMISSION_MAGIC.len();
    assert_eq!(
        read_u32_le(
            &source_admission_bytes,
            &mut admission_at,
            "source-admission version",
        ),
        IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_VERSION,
    );
    let problem_id_end = admission_at + 32;
    assert_eq!(
        source_admission_bytes.get(admission_at..problem_id_end),
        Some(problem.id().digest().as_bytes().as_slice()),
        "source admission must place ProblemId before the resolution registry",
    );
    admission_at = problem_id_end;
    assert_eq!(
        read_u32_le(
            &source_admission_bytes,
            &mut admission_at,
            "source-admission resolution count",
        ),
        u32::try_from(problem.source_resolutions().len()).expect("bounded resolution count"),
    );
    assert_eq!(
        hash_domain(
            IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_DOMAIN,
            &source_admission_bytes,
        ),
        problem.source_admission_id().digest(),
    );

    let execution = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    let exact_execution = execution
        .canonical_bytes()
        .expect("exact execution transport");
    let execution_preimage =
        project_exact_header_to_identity(&exact_execution, EXECUTION_MAGIC, execution.header());
    let execution_body_at =
        assert_identity_header_layout(&execution_preimage, EXECUTION_MAGIC, execution.header());
    assert_eq!(
        execution_preimage.get(execution_body_at..execution_body_at + 32),
        Some(execution.problem_id().digest().as_bytes().as_slice()),
        "execution identity must place ProblemId immediately after the projected header",
    );
    assert_eq!(
        execution_preimage.get(execution_body_at + 32..execution_body_at + 64),
        Some(
            execution
                .source_admission_id()
                .digest()
                .as_bytes()
                .as_slice()
        ),
        "execution identity must place SourceAdmissionId after ProblemId",
    );
    assert_eq!(
        hash_domain(
            IDENTIFIABILITY_EXECUTION_IDENTITY_DOMAIN,
            &execution_preimage,
        ),
        execution.id().expect("execution id").digest(),
    );
    let mut malformed_execution = exact_execution;
    malformed_execution[EXECUTION_MAGIC.len() + 4] = 0;
    assert!(matches!(
        IdentifiabilityExecutionPlan::from_canonical_bytes(
            &malformed_execution,
            &problem,
            execution.source_authority(),
        ),
        Err(IdentifiabilityError::Canonical { .. })
    ));

    let assessment = assessment(&problem, &execution, "layout-receipt");
    let exact_assessment = assessment
        .canonical_bytes()
        .expect("exact assessment transport");
    let assessment_preimage =
        project_exact_header_to_identity(&exact_assessment, ASSESSMENT_MAGIC, assessment.header());
    let assessment_body_at =
        assert_identity_header_layout(&assessment_preimage, ASSESSMENT_MAGIC, assessment.header());
    assert_eq!(
        assessment_preimage.get(assessment_body_at..assessment_body_at + 32),
        Some(assessment.problem_id().digest().as_bytes().as_slice()),
        "assessment identity must place ProblemId immediately after the projected header",
    );
    assert_eq!(
        assessment_preimage.get(assessment_body_at + 32..assessment_body_at + 64),
        Some(assessment.execution_id().digest().as_bytes().as_slice()),
        "assessment identity must place ExecutionId after ProblemId",
    );
    let mut claims_at = assessment_body_at + 64;
    assert_eq!(
        read_u32_le(&assessment_preimage, &mut claims_at, "claim count"),
        u32::try_from(assessment.claims().len()).expect("bounded claim count"),
    );
    assert_eq!(
        hash_domain(
            IDENTIFIABILITY_ASSESSMENT_IDENTITY_DOMAIN,
            &assessment_preimage,
        ),
        assessment.id().expect("assessment id").digest(),
    );
    let mut malformed_assessment = exact_assessment;
    malformed_assessment[ASSESSMENT_MAGIC.len() + 4] = 0;
    assert!(matches!(
        IdentifiabilityAssessment::from_canonical_bytes(
            &malformed_assessment,
            &problem,
            &execution,
            assessment.source_authority(),
        ),
        Err(IdentifiabilityError::Canonical { .. })
    ));
    log(
        "identity-wire-layout",
        "pass",
        "domains, version/count framing, numeric endianness, parent order, and header projection are exact",
    );
}

#[test]
fn trailing_bytes_and_stale_versions_refuse() {
    const MAGIC: &[u8] = b"fs-material-identifiability-problem\0";
    let problem = problem_fixture(ProblemOptions::default())
        .document
        .expect("problem");
    let mut trailing = problem.canonical_bytes().expect("bytes");
    trailing.push(0);
    assert!(IdentifiabilityProblemDocument::from_canonical_bytes(&trailing).is_err());
    let mut stale = problem.canonical_bytes().expect("bytes");
    stale[MAGIC.len()..MAGIC.len() + 4].copy_from_slice(&1_u32.to_le_bytes());
    assert!(matches!(
        IdentifiabilityProblemDocument::from_canonical_bytes(&stale),
        Err(IdentifiabilityError::UnsupportedSchemaVersion { .. })
    ));
    log(
        "canonical-adversaries",
        "pass",
        "trailing bytes and stale/future schema versions fail closed",
    );
}

#[test]
fn identifiability_identity_versions_and_transports_fail_closed() {
    const PROBLEM_MAGIC: &[u8] = b"fs-material-identifiability-problem\0";
    const EXECUTION_MAGIC: &[u8] = b"fs-material-identifiability-execution\0";
    const ASSESSMENT_MAGIC: &[u8] = b"fs-material-identifiability-assessment\0";

    let unresolved = problem_fixture(ProblemOptions::default())
        .document
        .expect("problem document");
    let mut problem_bytes = unresolved.canonical_bytes().expect("problem bytes");
    let mut bad_problem_magic = problem_bytes.clone();
    bad_problem_magic[0] ^= 0x01;
    assert!(IdentifiabilityProblemDocument::from_canonical_bytes(&bad_problem_magic).is_err());
    problem_bytes[PROBLEM_MAGIC.len()..PROBLEM_MAGIC.len() + 4]
        .copy_from_slice(&1_u32.to_le_bytes());
    assert!(matches!(
        IdentifiabilityProblemDocument::from_canonical_bytes(&problem_bytes),
        Err(IdentifiabilityError::UnsupportedSchemaVersion { .. })
    ));

    let problem = admit_fixture(problem_fixture(ProblemOptions::default()));
    let execution = execution(&problem, false, 17, 1.0e-10, false).expect("execution");
    let mut execution_bytes = execution.canonical_bytes().expect("execution bytes");
    let mut bad_execution_magic = execution_bytes.clone();
    bad_execution_magic[0] ^= 0x01;
    assert!(
        IdentifiabilityExecutionPlan::from_canonical_bytes(
            &bad_execution_magic,
            &problem,
            execution.source_authority(),
        )
        .is_err()
    );
    execution_bytes[EXECUTION_MAGIC.len()..EXECUTION_MAGIC.len() + 4]
        .copy_from_slice(&1_u32.to_le_bytes());
    assert!(matches!(
        IdentifiabilityExecutionPlan::from_canonical_bytes(
            &execution_bytes,
            &problem,
            execution.source_authority(),
        ),
        Err(IdentifiabilityError::UnsupportedSchemaVersion { .. })
    ));

    let assessment = assessment(&problem, &execution, "version-guard-receipt");
    let mut assessment_bytes = assessment.canonical_bytes().expect("assessment bytes");
    let mut bad_assessment_magic = assessment_bytes.clone();
    bad_assessment_magic[0] ^= 0x01;
    assert!(
        IdentifiabilityAssessment::from_canonical_bytes(
            &bad_assessment_magic,
            &problem,
            &execution,
            assessment.source_authority(),
        )
        .is_err()
    );
    assessment_bytes[ASSESSMENT_MAGIC.len()..ASSESSMENT_MAGIC.len() + 4]
        .copy_from_slice(&1_u32.to_le_bytes());
    assert!(matches!(
        IdentifiabilityAssessment::from_canonical_bytes(
            &assessment_bytes,
            &problem,
            &execution,
            assessment.source_authority(),
        ),
        Err(IdentifiabilityError::UnsupportedSchemaVersion { .. })
    ));
    assert!(check_source_admission_identity_version(1).is_err());
    log(
        "identity-stage-version-transports",
        "pass",
        "problem, source-admission, execution, and assessment versions fail closed independently",
    );
}

#[test]
fn source_ref_semantics_version_and_hash_are_mandatory() {
    assert!(
        SourceRef::try_new(
            source_key("zero"),
            SourceKind::Assumption,
            ContentHash([0; 32]),
            "fixture",
            1,
        )
        .is_err()
    );
    assert!(
        SourceRef::try_new(
            source_key("version-zero"),
            SourceKind::Assumption,
            hash("nonzero"),
            "fixture",
            0,
        )
        .is_err()
    );
    log(
        "source-ref-bounds",
        "pass",
        "zero identity and unpublished source semantics refused",
    );
}

#[test]
fn gauge_algebra_and_orbits_retain_ambitious_continuous_discrete_and_stratified_space() {
    let members = BTreeSet::from([role("yield_stress"), role("hardening_modulus")]);
    for (index, (algebra, orbit_geometry)) in [
        (
            GaugeAlgebra::Continuous {
                group_dimension: GaugeContinuousDimension::Finite { dimension: 1 },
            },
            GaugeOrbitGeometry::Regular {
                principal: RegularGaugeOrbit::new(
                    GaugeContinuousDimension::Finite { dimension: 1 },
                    GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 },
                ),
                stabilizer_profile: None,
            },
        ),
        (
            GaugeAlgebra::Discrete {
                size: GaugeDiscreteSize::Finite { order: 2 },
            },
            GaugeOrbitGeometry::Regular {
                principal: RegularGaugeOrbit::new(
                    GaugeContinuousDimension::Finite { dimension: 0 },
                    GaugeDiscreteOrbitCardinality::Finite { cardinality: 2 },
                ),
                stabilizer_profile: None,
            },
        ),
        (
            GaugeAlgebra::Mixed {
                continuous_group_dimension: GaugeContinuousDimension::Finite { dimension: 1 },
                component_group: GaugeDiscreteSize::Finite { order: 2 },
            },
            GaugeOrbitGeometry::Regular {
                principal: RegularGaugeOrbit::new(
                    GaugeContinuousDimension::Finite { dimension: 1 },
                    GaugeDiscreteOrbitCardinality::Finite { cardinality: 2 },
                ),
                stabilizer_profile: None,
            },
        ),
        (
            GaugeAlgebra::Continuous {
                group_dimension: GaugeContinuousDimension::Finite { dimension: 1 },
            },
            GaugeOrbitGeometry::Stratified {
                principal: RegularGaugeOrbit::new(
                    GaugeContinuousDimension::Finite { dimension: 1 },
                    GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 },
                ),
                orbit_type_stabilizer_profile: source_key("strata"),
            },
        ),
    ]
    .into_iter()
    .enumerate()
    {
        assert!(
            GaugeDeclaration::try_new(
                GaugeClassId::try_new(format!("kind-{index}")).expect("gauge id"),
                members.clone(),
                source_key("action"),
                algebra,
                orbit_geometry,
                GaugeStatus::Candidate {
                    rationale: source_key("candidate-rationale"),
                },
                gauge_validity_fixture(&members, true),
            )
            .is_ok()
        );
    }
    log(
        "gauge-kind-space",
        "pass",
        "schema preserves theorem-ready continuous/discrete/mixed/stratified gauges",
    );
}

#[test]
fn identity_version_guard_is_exact() {
    for checker in [
        check_authority_schema_version as fn(u32) -> Result<(), IdentifiabilityError>,
        check_problem_identity_version,
        check_source_admission_identity_version,
        check_execution_identity_version,
        check_assessment_identity_version,
    ] {
        assert!(checker(3).is_ok());
        assert!(checker(0).is_err());
        assert!(checker(1).is_err());
        assert!(checker(2).is_err());
        assert!(checker(4).is_err());
    }
    log(
        "authority-version-guard",
        "pass",
        "umbrella and all four identity-stage versions fail closed independently",
    );
}
