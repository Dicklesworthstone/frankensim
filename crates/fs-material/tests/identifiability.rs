#![cfg(any())]

//! Historical single-case prototype tests retained for design archaeology.
//! The authority-separated, multi-case v3 conformance suite lives in
//! `identifiability_authority.rs`; this obsolete prototype cannot mint current
//! identities and is deliberately excluded from compilation.
//! The body intentionally preserves removed constants and obsolete APIs, is
//! neither API-current nor compilable, and must not be treated as a re-enable
//! checklist. Any future migration needs a new semantic review against the
//! active contracts rather than mechanical uncommenting.
//!
//! I10.1 G0/G3 conformance for the law/experiment identifiability schema.
//! Each case emits one deterministic JSON line so a batch-verification run can
//! retain useful diagnostics without depending on a test runner's prose.

use std::collections::BTreeMap;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_evidence::vv::{
    ArtifactHeader, ArtifactId, ArtifactKind, ArtifactRef, CalibrationSplit, ClockSynchronization,
    CovarianceMatrix, DataAuthenticity, DeclaredBudget, ExperimentArtifact, ExperimentOrigin,
    InstrumentCalibration, ObservationId, ObservationManifest, ObservationManifestRow,
    ObservationSourceRef, QoiId, RepeatabilitySummary, SeedDeclaration, UnitId, VV_SCHEMA_VERSION,
};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, MATDB_SCHEMA_VERSION,
    MaterialCard, MaterialStateId, Provenance,
};
use fs_material::identifiability::*;
use fs_qty::{Dims, QuantitySpec};

const STRESS_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const DIMENSIONLESS: Dims = Dims([0; 6]);
const MAGIC: &[u8] = b"fs-material-identifiability-study\0";

fn log_case(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-material/identifiability\",\"case\":\"{case}\",\
         \"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn hash(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.fs-material.identifiability-test.v1",
        label.as_bytes(),
    )
}

fn artifact_id(value: &str) -> ArtifactId {
    ArtifactId::try_new(value).expect("fixture artifact id is valid")
}

fn qoi_id(value: &str) -> QoiId {
    QoiId::try_new(value).expect("fixture QoI id is valid")
}

fn row_id(value: &str) -> ObservationId {
    ObservationId::try_new(value).expect("fixture observation-row id is valid")
}

fn basic_header(id: &str) -> ArtifactHeader {
    ArtifactHeader::try_new(
        artifact_id(id),
        vec![UnitId::try_new("unitless").expect("fixture unit")],
        SeedDeclaration::Fixed(17),
        DeclaredBudget::Limit(1.0e-9),
        DeclaredBudget::Limit(10_000),
        DeclaredBudget::Limit(1 << 20),
        vec![("fixture".to_string(), "1".to_string())],
        vec!["fixture".to_string()],
    )
    .expect("basic Five-Explicits fixture header")
}

fn study_header(id: &str, law: u32, state: u32, protocol: u32, refinement: u32) -> ArtifactHeader {
    ArtifactHeader::try_new(
        artifact_id(id),
        vec![
            UnitId::try_new("Pa").expect("fixture pressure unit"),
            UnitId::try_new("s0").expect("equal-width fixture time unit"),
        ],
        SeedDeclaration::Fixed(0x1d3_171f_1ab1_17),
        DeclaredBudget::Limit(1.0e-8),
        DeclaredBudget::Limit(30_000),
        DeclaredBudget::Limit(16 << 20),
        vec![
            (
                "fs-material-identifiability".to_string(),
                IDENTIFIABILITY_SCHEMA_VERSION.to_string(),
            ),
            ("fs-evidence-vv".to_string(), VV_SCHEMA_VERSION.to_string()),
            ("fs-matdb".to_string(), MATDB_SCHEMA_VERSION.to_string()),
            ("constitutive-law".to_string(), law.to_string()),
            ("constitutive-state".to_string(), state.to_string()),
            ("experiment-protocol".to_string(), protocol.to_string()),
            ("refinement-policy".to_string(), refinement.to_string()),
        ],
        vec!["identifiability.study".to_string()],
    )
    .expect("study header is valid")
}

fn model_cards() -> (MaterialCard, ConstitutiveModelCard) {
    let mut parameters = BTreeMap::new();
    parameters.insert(
        "yield_stress".to_string(),
        LawParameter {
            value: 276.0e6,
            dims: STRESS_DIMS,
        },
    );
    parameters.insert(
        "hardening_modulus".to_string(),
        LawParameter {
            value: 1.2e9,
            dims: STRESS_DIMS,
        },
    );
    let model = ConstitutiveModelCard {
        law: LawId("j2-identifiability-fixture".to_string()),
        law_version: 3,
        parameters,
        state_schema_version: 2,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("temperature", 250.0, 450.0),
        sources: vec![hash("model-source")],
        provenance: Provenance {
            source: "fixture calibration report".to_string(),
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
    .expect("fixture material card is valid");
    (material, model)
}

fn experiment_and_split() -> (ExperimentArtifact, CalibrationSplit) {
    let source_bytes = hash("source-bytes");
    let blind_source = hash("source-row-blind-1");
    let instrument = artifact_id("instrument-1");
    let clock = artifact_id("clock-1");
    let typed_source = |locator_hash: ContentHash, extraction_receipt: &str| {
        ObservationSourceRef::try_new(
            source_bytes,
            "org.frankensim.fs-material.identifiability-test.row-locator.v1",
            1,
            locator_hash,
            hash(extraction_receipt),
        )
        .expect("typed observation source")
    };
    let manifest = ObservationManifest::try_new(vec![
        (
            row_id("cal-1"),
            ObservationManifestRow::try_new(
                typed_source(hash("source-row-cal-1"), "extraction-receipt-cal-1"),
                qoi_id("stress"),
                instrument.clone(),
                artifact_id("acquisition-channel-cal-1"),
                clock.clone(),
            )
            .expect("typed calibration row"),
        ),
        (
            row_id("cal-2"),
            ObservationManifestRow::try_new(
                typed_source(hash("source-row-cal-2"), "extraction-receipt-cal-2"),
                qoi_id("tangent"),
                instrument.clone(),
                artifact_id("acquisition-channel-cal-2"),
                clock.clone(),
            )
            .expect("typed calibration row"),
        ),
        (
            row_id("val-1"),
            ObservationManifestRow::try_new(
                typed_source(hash("source-row-val-1"), "extraction-receipt-val-1"),
                qoi_id("stress"),
                instrument.clone(),
                artifact_id("acquisition-channel-val-1"),
                clock.clone(),
            )
            .expect("typed validation row"),
        ),
        (
            row_id("blind-1"),
            ObservationManifestRow::try_new(
                typed_source(blind_source, "extraction-receipt-blind-1"),
                qoi_id("tangent"),
                instrument.clone(),
                artifact_id("acquisition-channel-blind-1"),
                clock.clone(),
            )
            .expect("typed blind row"),
        ),
    ])
    .expect("injective fixture manifest");
    let experiment = ExperimentArtifact::try_new(
        basic_header("experiment-1"),
        artifact_id("dataset-1"),
        ExperimentOrigin::Physical {
            apparatus_id: artifact_id("load-frame-1"),
            facility_id: artifact_id("lab-1"),
        },
        vec![qoi_id("stress"), qoi_id("tangent")],
        manifest,
        vec![InstrumentCalibration::new(
            instrument,
            hash("instrument-calibration"),
            true,
        )],
        ClockSynchronization::SingleClock { clock_id: clock },
        RepeatabilitySummary::try_new(
            4,
            CovarianceMatrix::try_new(2, vec![4.0, 0.5, 9.0])
                .expect("positive-definite repeatability covariance"),
        )
        .expect("repeatability fixture"),
        DataAuthenticity::new(source_bytes, hash("custody-receipt"), true),
    )
    .expect("experiment fixture admits");
    let experiment_ref = ArtifactRef::new(
        ArtifactKind::ExperimentArtifact,
        experiment.id().clone(),
        experiment.content_hash().expect("experiment hashes"),
    );
    let split = CalibrationSplit::try_new(
        basic_header("split-1"),
        experiment_ref,
        hash("split-preregistration"),
        vec![row_id("cal-1"), row_id("cal-2")],
        vec![row_id("val-1")],
        vec![(row_id("blind-1"), blind_source)],
    )
    .expect("split fixture admits");
    (experiment, split)
}

#[derive(Debug, Clone, Copy)]
enum CoordinateFixture {
    Identity,
    Affine,
}

fn parameter(
    role: &str,
    class: ParameterClass,
    observability: ParameterObservability,
    coordinate_fixture: CoordinateFixture,
    prior_version: u32,
) -> ParameterSpec {
    let domain = if role == "yield_stress" {
        ParameterDomain::try_new(1.0e6, 1.0e9).expect("yield domain")
    } else {
        ParameterDomain::try_new(1.0e7, 5.0e9).expect("hardening domain")
    };
    let (coordinate_domain, coordinate_quantity, transform, coordinate_name) =
        match coordinate_fixture {
            CoordinateFixture::Identity => (
                domain,
                QuantitySpec::dimensional(STRESS_DIMS),
                CoordinateTransform::Identity,
                format!("{role}-si"),
            ),
            CoordinateFixture::Affine => {
                let (lo, hi) = domain.bounds();
                (
                    ParameterDomain::try_new(lo / 1.0e6, hi / 1.0e6)
                        .expect("scaled coordinate domain"),
                    QuantitySpec::dimensional(STRESS_DIMS),
                    CoordinateTransform::Affine {
                        scale: 1.0e6,
                        offset: 0.0,
                    },
                    format!("{role}-mpa"),
                )
            }
        };
    ParameterSpec::try_new(
        ParameterRoleId::try_new(role).expect("parameter role"),
        QuantitySpec::dimensional(STRESS_DIMS),
        domain,
        ParameterPrior::Uniform {
            domain,
            version: prior_version,
        },
        ParameterCoordinate::try_new(
            CoordinateId::try_new(coordinate_name).expect("coordinate id"),
            coordinate_quantity,
            coordinate_domain,
            transform,
        )
        .expect("coordinate fixture"),
        ParameterOwner::ConstitutiveModel,
        ParameterScope::Global,
        class,
        observability,
    )
    .expect("parameter fixture")
}

fn sensor(channel: &str, clock: &str) -> SensorBinding {
    SensorBinding::try_new(
        artifact_id("instrument-1"),
        artifact_id(&format!("instrument-1-{channel}")),
        hash(&format!("sensor-model-{channel}")),
        5,
        hash("instrument-calibration"),
        hash(&format!("transfer-{channel}")),
        hash(&format!("filter-{channel}")),
        hash(&format!("support-{channel}")),
        artifact_id(clock),
        2_500,
        hash(&format!("anti-aliasing-{channel}")),
    )
    .expect("sensor fixture")
}

fn observation(
    channel: &str,
    qoi: &str,
    row: &str,
    frame: &FrameBinding,
    clock: &str,
    protocol_version: u32,
    refinement_version: u32,
) -> ObservationSpec {
    ObservationSpec::try_new(
        ObservationChannelId::try_new(channel).expect("observation channel"),
        qoi_id(qoi),
        QuantitySpec::dimensional(STRESS_DIMS),
        frame.clone(),
        format!("node-{channel}"),
        "stress-output",
        hash(&format!("operator-{channel}")),
        4,
        hash(&format!("aggregation-{channel}")),
        sensor(channel, clock),
        NoiseModel::Gaussian {
            standard_deviation: if channel == "stress" { 2.0e5 } else { 3.0e5 },
        },
        MissingnessModel::Complete {
            evidence: hash(&format!("complete-{channel}")),
        },
        None,
        protocol_version,
        refinement_version,
        vec![row_id(row)],
    )
    .expect("observation fixture")
}

fn not_assessed(reason: &str) -> EvidenceStatus {
    EvidenceStatus::NotAssessed {
        reason: reason.to_string(),
    }
}

fn evidence() -> IdentifiabilityEvidence {
    IdentifiabilityEvidence::try_new(
        not_assessed("structural checker has not run"),
        not_assessed("local checker has not run"),
        not_assessed("generic checker has not run"),
        not_assessed("global checker has not run"),
        not_assessed("practical checker has not run"),
    )
    .expect("explicit five-axis no-claim")
}

#[derive(Debug, Clone)]
struct StudyFixture {
    header: ArtifactHeader,
    context: ArtifactRef,
    model: MaterialModelBinding,
    initial_state: InitialStateBinding,
    specimen: SpecimenBinding,
    protocol: ProtocolBinding,
    data: DataLineage,
    parameters: Vec<ParameterSpec>,
    gauges: Vec<GaugeClass>,
    observations: Vec<ObservationSpec>,
    paths: Vec<ObservationPath>,
    noise_dependence: NoiseDependence,
    discrepancies: Vec<DiscrepancySpec>,
    evidence: IdentifiabilityEvidence,
}

impl StudyFixture {
    fn spec(self) -> Result<IdentifiabilityStudySpec, IdentifiabilityError> {
        IdentifiabilityStudySpec::try_new(
            self.header,
            self.context,
            self.model,
            self.initial_state,
            self.specimen,
            self.protocol,
            self.data,
            self.parameters,
            self.gauges,
            self.observations,
            self.paths,
            self.noise_dependence,
            self.discrepancies,
            self.evidence,
        )
    }

    fn admit(self) -> Result<AdmittedIdentifiabilityStudy, IdentifiabilityError> {
        AdmittedIdentifiabilityStudy::admit(self.spec()?)
    }
}

fn fixture() -> StudyFixture {
    let (material, model_card) = model_cards();
    let model = MaterialModelBinding::from_cards(&material, &model_card, hash("graph-v1"))
        .expect("material/model binding");
    let (experiment, split) = experiment_and_split();
    let data = DataLineage::from_vv(
        &experiment,
        &split,
        hash("parser"),
        2,
        hash("preprocessing"),
        artifact_id("split-by-specimen"),
    )
    .expect("data lineage fixture");
    let frame = FrameBinding::try_new(
        artifact_id("specimen-frame"),
        hash("specimen-frame-transform"),
        "right-handed-cartesian",
    )
    .expect("frame fixture");
    let specimen = SpecimenBinding::try_new(
        artifact_id("specimen-1"),
        hash("specimen-geometry"),
        hash("specimen-process"),
        hash("specimen-preparation"),
        frame.clone(),
    )
    .expect("specimen fixture");
    let protocol = ProtocolBinding::try_new(
        artifact_id("protocol-1"),
        7,
        2,
        3,
        hash("load-path"),
        hash("environment-path"),
        hash("time-grid"),
        artifact_id("clock-1"),
    )
    .expect("protocol fixture");
    let yield_parameter = parameter(
        "yield_stress",
        ParameterClass::Target,
        ParameterObservability::Candidate,
        CoordinateFixture::Identity,
        2,
    );
    let hardening_parameter = parameter(
        "hardening_modulus",
        ParameterClass::Nuisance {
            calibration: data.split().clone(),
        },
        ParameterObservability::Candidate,
        CoordinateFixture::Identity,
        2,
    );
    let stress = observation("stress", "stress", "cal-1", &frame, "clock-1", 7, 3);
    let tangent = observation("tangent", "tangent", "cal-2", &frame, "clock-1", 7, 3);
    let paths = vec![
        ObservationPath::try_new(
            ParameterRoleId::try_new("yield_stress").expect("role"),
            ObservationChannelId::try_new("stress").expect("channel"),
            InfluenceMechanism::Mean,
            hash("path-yield-stress"),
            QuantitySpec::dimensional(DIMENSIONLESS),
            InfluenceStatus::DeclaredConnectivity,
        )
        .expect("yield path"),
        ObservationPath::try_new(
            ParameterRoleId::try_new("hardening_modulus").expect("role"),
            ObservationChannelId::try_new("tangent").expect("channel"),
            InfluenceMechanism::StateMediated,
            hash("path-hardening-tangent"),
            QuantitySpec::dimensional(DIMENSIONLESS),
            InfluenceStatus::NumericallyWitnessed {
                receipt: hash("hardening-path-witness"),
            },
        )
        .expect("hardening path"),
    ];
    let noise_dependence = NoiseDependence::try_new(
        vec![
            ObservationChannelId::try_new("tangent").expect("channel"),
            ObservationChannelId::try_new("stress").expect("channel"),
        ],
        CovarianceMatrix::try_new(2, vec![1.0, 0.25, 1.0]).expect("correlation fixture"),
        hash("correlation-evidence"),
    )
    .expect("noise dependence fixture");
    StudyFixture {
        header: study_header("study-1", 3, 2, 7, 3),
        context: ArtifactRef::new(
            ArtifactKind::ContextOfUse,
            artifact_id("context-1"),
            hash("context-1"),
        ),
        model,
        initial_state: InitialStateBinding::Zero { schema_version: 2 },
        specimen,
        protocol,
        data,
        parameters: vec![hardening_parameter, yield_parameter],
        gauges: Vec::new(),
        observations: vec![tangent, stress],
        paths,
        noise_dependence,
        discrepancies: vec![
            DiscrepancySpec::try_new(
                ObservationChannelId::try_new("tangent").expect("channel"),
                DiscrepancyModel::NoModel {
                    reason: "no tangent discrepancy family has been admitted".to_string(),
                },
            )
            .expect("tangent discrepancy"),
            DiscrepancySpec::try_new(
                ObservationChannelId::try_new("stress").expect("channel"),
                DiscrepancyModel::NoModel {
                    reason: "no stress discrepancy family has been admitted".to_string(),
                },
            )
            .expect("stress discrepancy"),
        ],
        evidence: evidence(),
    }
}

fn unit_count_offset(bytes: &[u8]) -> usize {
    let mut at = MAGIC.len() + 4 + 1 + 1;
    let id_len = u32::from_le_bytes(bytes[at..at + 4].try_into().expect("id length")) as usize;
    at += 4 + id_len;
    at
}

fn swap_two_equal_width_header_units(bytes: &mut [u8]) {
    let mut at = unit_count_offset(bytes) + 4;
    let first_len =
        u32::from_le_bytes(bytes[at..at + 4].try_into().expect("first unit length")) as usize;
    at += 4;
    let first = at;
    at += first_len;
    let second_len =
        u32::from_le_bytes(bytes[at..at + 4].try_into().expect("second unit length")) as usize;
    at += 4;
    assert_eq!(first_len, second_len, "fixture units must be equal width");
    for offset in 0..first_len {
        bytes.swap(first + offset, at + offset);
    }
}

#[test]
fn canonical_roundtrip_retains_exact_and_physical_preimages() {
    let admitted = fixture().admit().expect("canonical fixture admits");
    let bytes = admitted.exact_receipt().canonical_bytes().to_vec();
    let decoded = IdentifiabilityStudySpec::from_canonical_bytes(&bytes)
        .expect("exact canonical bytes decode");
    let replay = AdmittedIdentifiabilityStudy::admit(decoded).expect("decoded study re-admits");

    assert_eq!(replay.exact_receipt(), admitted.exact_receipt());
    assert_eq!(replay.physical_receipt(), admitted.physical_receipt());
    assert_eq!(admitted.exact_receipt().schema_version(), 2);
    assert_eq!(admitted.exact_receipt().item_count(), 8);
    assert_ne!(
        admitted.exact_receipt().id().digest(),
        admitted.physical_receipt().id().digest(),
        "identity domains remain type- and digest-separated"
    );
    log_case(
        "i10-g0-canonical-roundtrip",
        "pass",
        "decode re-admission preserves both complete identity receipts",
    );
}

#[test]
fn insertion_order_is_canonical_and_covariance_order_is_quotiented() {
    let canonical = fixture().admit().expect("baseline admits");
    let mut permuted = fixture();
    permuted.parameters.reverse();
    permuted.observations.reverse();
    permuted.paths.reverse();
    permuted.discrepancies.reverse();
    permuted.noise_dependence = NoiseDependence::try_new(
        vec![
            ObservationChannelId::try_new("stress").expect("channel"),
            ObservationChannelId::try_new("tangent").expect("channel"),
        ],
        CovarianceMatrix::try_new(2, vec![1.0, 0.25, 1.0]).expect("permuted correlation"),
        hash("correlation-evidence"),
    )
    .expect("permuted covariance order is valid");
    let permuted = permuted.admit().expect("permuted fixture admits");

    assert_eq!(canonical.exact_receipt(), permuted.exact_receipt());
    assert_eq!(canonical.physical_receipt(), permuted.physical_receipt());
    log_case(
        "i10-g0-order-canonicalization",
        "pass",
        "map set path and covariance order permutations mint identical receipts",
    );
}

#[test]
fn affine_reparameterization_moves_replay_identity_but_not_physical_identity() {
    let baseline = fixture().admit().expect("baseline admits");
    let mut reparameterized = fixture();
    reparameterized.parameters[1] = parameter(
        "yield_stress",
        ParameterClass::Target,
        ParameterObservability::Candidate,
        CoordinateFixture::Affine,
        2,
    );
    let reparameterized = reparameterized.admit().expect("affine chart admits");

    assert_ne!(
        baseline.exact_receipt().id(),
        reparameterized.exact_receipt().id()
    );
    assert_eq!(
        baseline.physical_receipt().id(),
        reparameterized.physical_receipt().id()
    );
    assert_ne!(
        baseline.exact_receipt().canonical_bytes(),
        reparameterized.exact_receipt().canonical_bytes()
    );
    assert_eq!(
        baseline.physical_receipt().canonical_bytes(),
        reparameterized.physical_receipt().canonical_bytes()
    );
    log_case(
        "i10-g3-coordinate-quotient",
        "pass",
        "validated affine chart changes exact replay but preserves physical problem identity",
    );
}

#[test]
fn physical_prior_semantics_are_not_erased_by_coordinate_quotient() {
    let baseline = fixture().admit().expect("baseline admits");
    let mut changed = fixture();
    changed.parameters[1] = parameter(
        "yield_stress",
        ParameterClass::Target,
        ParameterObservability::Candidate,
        CoordinateFixture::Identity,
        3,
    );
    let changed = changed.admit().expect("new prior version admits");

    assert_ne!(baseline.exact_receipt().id(), changed.exact_receipt().id());
    assert_ne!(
        baseline.physical_receipt().id(),
        changed.physical_receipt().id()
    );
    log_case(
        "i10-g3-prior-identity",
        "pass",
        "prior semantics version moves exact and physical identities",
    );
}

#[test]
fn disconnected_estimated_parameters_refuse_until_explicitly_unidentifiable() {
    let mut disconnected = fixture();
    disconnected
        .paths
        .retain(|path| path.parameter().as_str() == "hardening_modulus");
    assert!(matches!(
        disconnected.clone().spec(),
        Err(IdentifiabilityError::DisconnectedEstimatedParameter { parameter })
            if parameter.as_str() == "yield_stress"
    ));

    disconnected.parameters[1] = parameter(
        "yield_stress",
        ParameterClass::Target,
        ParameterObservability::ExplicitlyUnidentifiable {
            reason: "the selected protocol supplies no yield-sensitive observable".to_string(),
            witness: hash("yield-unidentifiable-witness"),
        },
        CoordinateFixture::Identity,
        2,
    );
    disconnected
        .admit()
        .expect("honest explicit-unidentifiable declaration admits");
    log_case(
        "i10-g0-disconnected-parameter",
        "pass",
        "estimated disconnection refuses and explicit witnessed no-claim admits",
    );
}

#[test]
fn nuisance_calibration_and_blind_data_leakage_fail_closed() {
    let mut wrong_nuisance = fixture();
    wrong_nuisance.parameters[0] = parameter(
        "hardening_modulus",
        ParameterClass::Nuisance {
            calibration: ArtifactRef::new(
                ArtifactKind::CalibrationSplit,
                artifact_id("split-1"),
                hash("wrong-split"),
            ),
        },
        ParameterObservability::Candidate,
        CoordinateFixture::Identity,
        2,
    );
    assert!(matches!(
        wrong_nuisance.spec(),
        Err(IdentifiabilityError::NuisanceCalibration { parameter })
            if parameter.as_str() == "hardening_modulus"
    ));

    let mut leaked = fixture();
    let frame = leaked.specimen.frame().clone();
    leaked.observations[1] = observation("stress", "stress", "val-1", &frame, "clock-1", 7, 3);
    assert!(matches!(
        leaked.spec(),
        Err(IdentifiabilityError::Vv { .. })
    ));
    log_case(
        "i10-g0-nuisance-and-blind-lineage",
        "pass",
        "wrong split and validation-row consumption are both refused before admission",
    );
}

#[test]
fn state_protocol_refinement_and_clock_versions_refuse_independently() {
    let mut wrong_state = fixture();
    wrong_state.initial_state = InitialStateBinding::Zero { schema_version: 9 };
    assert!(matches!(
        wrong_state.spec(),
        Err(IdentifiabilityError::VersionMismatch {
            field: "initial state schema",
            ..
        })
    ));

    let mut wrong_protocol = fixture();
    let frame = wrong_protocol.specimen.frame().clone();
    wrong_protocol.observations[1] =
        observation("stress", "stress", "cal-1", &frame, "clock-1", 8, 3);
    assert!(matches!(
        wrong_protocol.spec(),
        Err(IdentifiabilityError::VersionMismatch {
            field: "observation protocol",
            ..
        })
    ));

    let mut wrong_refinement = fixture();
    let frame = wrong_refinement.specimen.frame().clone();
    wrong_refinement.observations[1] =
        observation("stress", "stress", "cal-1", &frame, "clock-1", 7, 4);
    assert!(matches!(
        wrong_refinement.spec(),
        Err(IdentifiabilityError::VersionMismatch {
            field: "observation refinement",
            ..
        })
    ));

    let mut wrong_clock = fixture();
    let frame = wrong_clock.specimen.frame().clone();
    wrong_clock.observations[1] =
        observation("stress", "stress", "cal-1", &frame, "clock-rogue", 7, 3);
    assert!(matches!(
        wrong_clock.spec(),
        Err(IdentifiabilityError::UnknownReference {
            field: "observation clock",
            ..
        })
    ));
    log_case(
        "i10-g0-version-and-clock-closure",
        "pass",
        "state protocol refinement and clock mismatches reach distinct typed refusals",
    );
}

#[test]
fn covariance_channel_closure_and_normalization_fail_closed() {
    let bad_diagonal = NoiseDependence::try_new(
        vec![ObservationChannelId::try_new("stress").expect("channel")],
        CovarianceMatrix::try_new(1, vec![2.0]).expect("positive scalar matrix"),
        hash("correlation-evidence"),
    );
    assert!(matches!(
        bad_diagonal,
        Err(IdentifiabilityError::Covariance { .. })
    ));

    let mut missing_channel = fixture();
    missing_channel.noise_dependence = NoiseDependence::try_new(
        vec![ObservationChannelId::try_new("stress").expect("channel")],
        CovarianceMatrix::try_new(1, vec![1.0]).expect("unit scalar correlation"),
        hash("correlation-evidence"),
    )
    .expect("locally valid one-channel dependence");
    assert!(matches!(
        missing_channel.spec(),
        Err(IdentifiabilityError::Covariance { .. })
    ));
    log_case(
        "i10-g0-unit-safe-correlation",
        "pass",
        "nonunit diagonal and incomplete channel order both refuse",
    );
}

#[test]
fn gauge_membership_and_fixed_parameter_rules_fail_closed() {
    let valid_gauge = GaugeClass::try_new(
        GaugeClassId::try_new("product-confounding").expect("gauge id"),
        vec![
            ParameterRoleId::try_new("yield_stress").expect("role"),
            ParameterRoleId::try_new("hardening_modulus").expect("role"),
        ],
        1,
        hash("gauge-action"),
        hash("gauge-quotient"),
        hash("gauge-slice"),
        hash("gauge-strata"),
        not_assessed("gauge action is declared but not proved"),
    )
    .expect("structurally valid gauge declaration");
    let mut gauged = fixture();
    gauged.gauges.push(valid_gauge);
    gauged.admit().expect("known nonfixed gauge members admit");

    let mut dangling = fixture();
    dangling.gauges.push(
        GaugeClass::try_new(
            GaugeClassId::try_new("dangling-gauge").expect("gauge id"),
            vec![
                ParameterRoleId::try_new("yield_stress").expect("role"),
                ParameterRoleId::try_new("not-a-parameter").expect("role"),
            ],
            1,
            hash("dangling-action"),
            hash("dangling-quotient"),
            hash("dangling-slice"),
            hash("dangling-strata"),
            not_assessed("gauge declaration is only a fixture"),
        )
        .expect("gauge constructor cannot yet see study endpoints"),
    );
    assert!(matches!(
        dangling.spec(),
        Err(IdentifiabilityError::UnknownReference {
            field: "gauge member",
            ..
        })
    ));
    log_case(
        "i10-g0-gauge-closure",
        "pass",
        "valid gauge declaration admits while dangling membership refuses",
    );
}

#[test]
fn discrepancy_absence_zero_and_modeled_states_remain_distinct() {
    let baseline = fixture().admit().expect("baseline admits");
    let mut missing = fixture();
    missing.discrepancies.pop();
    assert!(matches!(
        missing.spec(),
        Err(IdentifiabilityError::Cardinality {
            field: "discrepancy rows",
            ..
        })
    ));

    let mut zero = fixture();
    zero.discrepancies[1] = DiscrepancySpec::try_new(
        ObservationChannelId::try_new("stress").expect("channel"),
        DiscrepancyModel::Zero {
            evidence: hash("zero-discrepancy-proof"),
        },
    )
    .expect("evidence-backed zero discrepancy");
    let zero = zero.admit().expect("zero-discrepancy study admits");
    assert_ne!(baseline.exact_receipt().id(), zero.exact_receipt().id());
    assert_ne!(
        baseline.physical_receipt().id(),
        zero.physical_receipt().id()
    );
    log_case(
        "i10-g0-discrepancy-tristate",
        "pass",
        "missing refuses and explicit zero differs from explicit no-model",
    );
}

#[test]
fn five_identifiability_axes_are_orthogonal_identity_fields() {
    let supported = EvidenceStatus::Supported {
        method: "differential-algebra".to_string(),
        receipt: hash("structural-receipt"),
    };
    let mut structural = fixture();
    structural.evidence = IdentifiabilityEvidence::try_new(
        supported.clone(),
        not_assessed("local checker has not run"),
        not_assessed("generic checker has not run"),
        not_assessed("global checker has not run"),
        not_assessed("practical checker has not run"),
    )
    .expect("structural evidence fixture");
    let structural = structural
        .admit()
        .expect("structural-evidence study admits");

    let mut local = fixture();
    local.evidence = IdentifiabilityEvidence::try_new(
        not_assessed("structural checker has not run"),
        supported,
        not_assessed("generic checker has not run"),
        not_assessed("global checker has not run"),
        not_assessed("practical checker has not run"),
    )
    .expect("local evidence fixture");
    let local = local.admit().expect("local-evidence study admits");

    assert!(matches!(
        structural.spec().evidence().structural(),
        EvidenceStatus::Supported { .. }
    ));
    assert!(matches!(
        local.spec().evidence().local(),
        EvidenceStatus::Supported { .. }
    ));
    assert_ne!(structural.exact_receipt().id(), local.exact_receipt().id());
    log_case(
        "i10-g0-evidence-axis-orthogonality",
        "pass",
        "structural and local support are separate fields and separate identities",
    );
}

#[test]
fn canonical_decoder_refuses_version_trailing_truncation_bombs_and_reordering() {
    let admitted = fixture().admit().expect("baseline admits");
    let canonical = admitted.exact_receipt().canonical_bytes();

    let mut future = canonical.to_vec();
    future[MAGIC.len()..MAGIC.len() + 4]
        .copy_from_slice(&(IDENTIFIABILITY_SCHEMA_VERSION + 1).to_le_bytes());
    assert!(matches!(
        IdentifiabilityStudySpec::from_canonical_bytes(&future),
        Err(IdentifiabilityError::UnsupportedSchemaVersion { .. })
    ));

    let mut trailing = canonical.to_vec();
    trailing.push(0);
    assert!(matches!(
        IdentifiabilityStudySpec::from_canonical_bytes(&trailing),
        Err(IdentifiabilityError::Canonical { .. })
    ));
    assert!(matches!(
        IdentifiabilityStudySpec::from_canonical_bytes(&canonical[..canonical.len() - 1]),
        Err(IdentifiabilityError::Canonical { .. })
    ));

    let mut count_bomb = canonical.to_vec();
    let count_at = unit_count_offset(&count_bomb);
    count_bomb[count_at..count_at + 4].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(matches!(
        IdentifiabilityStudySpec::from_canonical_bytes(&count_bomb),
        Err(IdentifiabilityError::Canonical { .. })
    ));

    let mut reordered = canonical.to_vec();
    swap_two_equal_width_header_units(&mut reordered);
    assert!(matches!(
        IdentifiabilityStudySpec::from_canonical_bytes(&reordered),
        Err(IdentifiabilityError::Canonical { .. })
    ));
    log_case(
        "i10-g0-adversarial-codec",
        "pass",
        "future trailing truncated oversized-count and noncanonical-order inputs all refuse",
    );
}

#[test]
fn numeric_admission_canonicalizes_signed_zero_and_refuses_invalid_transforms() {
    let domain = ParameterDomain::try_new(-0.0, 1.0).expect("signed zero is canonicalizable");
    assert_eq!(domain.bounds().0.to_bits(), 0.0f64.to_bits());
    assert!(ParameterDomain::try_new(f64::NAN, 1.0).is_err());
    assert!(ParameterDomain::try_new(2.0, 1.0).is_err());
    assert!(
        ParameterCoordinate::try_new(
            CoordinateId::try_new("bad-affine").expect("coordinate id"),
            QuantitySpec::dimensional(STRESS_DIMS),
            ParameterDomain::try_new(0.0, 1.0).expect("coordinate domain"),
            CoordinateTransform::Affine {
                scale: 0.0,
                offset: 0.0,
            },
        )
        .is_err()
    );
    assert!(
        ParameterCoordinate::try_new(
            CoordinateId::try_new("bad-log").expect("coordinate id"),
            QuantitySpec::dimensional(DIMENSIONLESS),
            ParameterDomain::try_new(0.0, 1.0).expect("coordinate domain"),
            CoordinateTransform::LogPositive { reference: 0.0 },
        )
        .is_err()
    );
    log_case(
        "i10-g0-numeric-canonicalization",
        "pass",
        "signed zero canonicalizes while NaN reversed domains and singular transforms refuse",
    );
}
