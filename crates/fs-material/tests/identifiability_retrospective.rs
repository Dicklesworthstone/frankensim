//! I10.1 G0/G3 tests for retrospective identifiability admission.
//!
//! These tests deliberately cross the real `fs-evidence` artifact boundary:
//! every admitted case carries a canonical `ExperimentArtifact` and matching
//! `CalibrationSplit`.  Fixed JSON logs make the eventual batch-verification
//! evidence useful without treating structural admission as a theorem about
//! scientific identifiability or laboratory authenticity.

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
const TEST_HASH_DOMAIN: &str = "org.frankensim.fs-material.identifiability-retrospective-test.v1";
const TEST_ROW_LOCATOR_DOMAIN: &str = "org.frankensim.fixture.row-locator.v1";
const TEST_ALTERNATIVE_ROW_LOCATOR_DOMAIN: &str =
    "org.frankensim.fixture.alternative-row-locator.v2";

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

fn log(case: &str, verdict: &str, expected: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-material/identifiability-retrospective\",\
         \"case\":\"{}\",\"verdict\":\"{}\",\"expected\":\"{}\",\
         \"detail\":\"{}\"}}",
        escape_json(case),
        escape_json(verdict),
        escape_json(expected),
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
    ArtifactId::try_new(value).expect("fixture artifact id")
}

fn qoi(value: &str) -> QoiId {
    QoiId::try_new(value).expect("fixture QoI id")
}

fn unit(value: &str) -> UnitId {
    UnitId::try_new(value).expect("fixture unit id")
}

fn axis(value: &str) -> AxisId {
    AxisId::try_new(value).expect("fixture axis id")
}

fn observation(value: &str) -> ObservationId {
    ObservationId::try_new(value).expect("fixture observation id")
}

fn case_id(value: &str) -> CaseId {
    CaseId::try_new(value).expect("fixture case id")
}

fn channel(value: &str) -> ObservationChannelId {
    ObservationChannelId::try_new(value).expect("fixture observation channel")
}

fn role(value: &str) -> ParameterRoleId {
    ParameterRoleId::try_new(value).expect("fixture parameter role")
}

fn source_key(value: &str) -> SourceKey {
    SourceKey::try_new(value).expect("fixture source key")
}

fn header(id: &str, units: &[&str], capability: &str) -> ArtifactHeader {
    ArtifactHeader::try_new(
        artifact(id),
        units.iter().copied().map(unit).collect(),
        SeedDeclaration::Fixed(0x171f_10_1),
        DeclaredBudget::Limit(1.0e-9),
        DeclaredBudget::Limit(30_000),
        DeclaredBudget::Limit(32 << 20),
        vec![(
            "fixture".to_string(),
            "identifiability-retrospective-v3".to_string(),
        )],
        vec![capability.to_string()],
    )
    .expect("Five Explicits fixture")
}

fn context() -> ContextOfUse {
    ContextOfUse::try_new(
        header("retrospective-context", &["Pa", "K"], "fixture.context"),
        "Calibrate a constitutive parameter without crossing preregistered evidence partitions.",
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
                    .expect("temperature applicability axis"),
            ],
            Vec::new(),
        )
        .expect("applicability domain"),
        ApplicabilityPolicy::Demote,
    )
    .expect("context fixture")
}

fn model_cards() -> (MaterialCard, ConstitutiveModelCard) {
    let model = ConstitutiveModelCard {
        law: LawId("retrospective-identifiability-fixture".to_string()),
        law_version: 1,
        parameters: BTreeMap::from([(
            "yield_stress".to_string(),
            LawParameter {
                value: 276.0e6,
                dims: STRESS,
            },
        )]),
        state_schema_version: 2,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("temperature", 250.0, 450.0),
        sources: vec![hash("model-source")],
        provenance: Provenance {
            source: "retrospective admission fixture".to_string(),
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

fn physical_data(
    label: &str,
    experiment_qois: &[&str],
    rows: &[(&str, &str)],
    calibration: &[&str],
    validation: &[&str],
    blind: &[&str],
    source_bytes_label: &str,
) -> (ExperimentArtifact, CalibrationSplit) {
    physical_data_with_metrology(
        label,
        label,
        experiment_qois,
        rows,
        calibration,
        validation,
        blind,
        source_bytes_label,
    )
}

#[allow(clippy::too_many_arguments)]
fn physical_data_with_metrology(
    label: &str,
    metrology_label: &str,
    experiment_qois: &[&str],
    rows: &[(&str, &str)],
    calibration: &[&str],
    validation: &[&str],
    blind: &[&str],
    source_bytes_label: &str,
) -> (ExperimentArtifact, CalibrationSplit) {
    let source_by_row = rows
        .iter()
        .map(|(id, source)| ((*id).to_string(), hash(source)))
        .collect::<BTreeMap<_, _>>();
    let qois = experiment_qois.iter().copied().map(qoi).collect::<Vec<_>>();
    assert!(
        !qois.is_empty() && rows.len() >= qois.len(),
        "fixture needs at least one row for every declared QoI",
    );
    let instrument = artifact(&format!("instrument-{metrology_label}"));
    let clock = artifact(&format!("clock-{metrology_label}"));
    let dataset_source_bytes_hash = hash(source_bytes_label);
    let manifest = ObservationManifest::try_new(
        rows.iter()
            .enumerate()
            .map(|(index, (id, source))| {
                (
                    observation(id),
                    ObservationManifestRow::try_new(
                        ObservationSourceRef::try_new(
                            dataset_source_bytes_hash,
                            TEST_ROW_LOCATOR_DOMAIN,
                            1,
                            hash(source),
                            hash(&format!("extraction-{metrology_label}-{id}")),
                        )
                        .expect("typed row-source fixture"),
                        qois[index % qois.len()].clone(),
                        instrument.clone(),
                        artifact(&format!("acquisition-channel-{id}")),
                        clock.clone(),
                    )
                    .expect("typed retrospective manifest row"),
                )
            })
            .collect(),
    )
    .expect("injective experiment manifest");
    let experiment = ExperimentArtifact::try_new(
        header(
            &format!("experiment-{label}"),
            &["Pa"],
            "fixture.experiment",
        ),
        artifact(&format!("dataset-{label}")),
        ExperimentOrigin::Physical {
            apparatus_id: artifact(&format!("apparatus-{label}")),
            facility_id: artifact("facility-retrospective-fixture"),
        },
        qois,
        manifest,
        vec![InstrumentCalibration::new(
            instrument,
            hash(&format!("sensor-{metrology_label}")),
            true,
        )],
        ClockSynchronization::SingleClock { clock_id: clock },
        RepeatabilitySummary::try_new(
            3,
            CovarianceMatrix::try_new(
                experiment_qois.len(),
                (0..experiment_qois.len())
                    .flat_map(|row| {
                        (0..=row).map(move |column| if row == column { 0.25 } else { 0.0 })
                    })
                    .collect(),
            )
            .expect("positive-semidefinite repeatability covariance"),
        )
        .expect("repeatability fixture"),
        DataAuthenticity::new(
            dataset_source_bytes_hash,
            hash(&format!("custody-{label}")),
            true,
        ),
    )
    .expect("experiment fixture");
    let experiment_hash = experiment.content_hash().expect("experiment hashes");
    let split = CalibrationSplit::try_new(
        header(&format!("split-{label}"), &["unitless"], "fixture.split"),
        ArtifactRef::new(
            ArtifactKind::ExperimentArtifact,
            experiment.id().clone(),
            experiment_hash,
        ),
        hash(&format!("preregistration-{label}")),
        calibration.iter().copied().map(observation).collect(),
        validation.iter().copied().map(observation).collect(),
        blind
            .iter()
            .map(|id| {
                (
                    observation(id),
                    *source_by_row
                        .get(*id)
                        .expect("blind row must exist in manifest"),
                )
            })
            .collect(),
    )
    .expect("calibration split fixture");
    (experiment, split)
}

fn blind_release_for(split: &CalibrationSplit, authority_label: &str) -> BlindReleaseReceipt {
    BlindReleaseReceipt::new(
        ArtifactRef::new(
            ArtifactKind::CalibrationSplit,
            split.id().clone(),
            split.content_hash().expect("split release hash"),
        ),
        split.blind_commitment(),
        hash(authority_label),
    )
    .expect("blind release fixture")
}

#[derive(Clone)]
struct RetrospectiveCaseFixture {
    id: &'static str,
    purpose: CasePurpose,
    observation_qoi: &'static str,
    observation_row: &'static str,
    experiment_key: &'static str,
    split_key: &'static str,
    experiment: ExperimentArtifact,
    split: CalibrationSplit,
    observation_instrument: Option<ArtifactId>,
    observation_acquisition_channel: Option<ArtifactId>,
    observation_clock: Option<ArtifactId>,
    observation_sensor: Option<SourceKey>,
    duplicate_row_channel: bool,
    declare_duplicate_row_sharing: bool,
    blind_release: Option<BlindReleaseReceipt>,
}

impl RetrospectiveCaseFixture {
    fn with_observation_instrument(mut self, instrument: &str) -> Self {
        self.observation_instrument = Some(artifact(instrument));
        self
    }

    fn with_observation_acquisition_channel(mut self, acquisition_channel: &str) -> Self {
        self.observation_acquisition_channel = Some(artifact(acquisition_channel));
        self
    }

    fn with_observation_clock(mut self, clock: &str) -> Self {
        self.observation_clock = Some(artifact(clock));
        self
    }

    fn with_observation_sensor(mut self, sensor: &str) -> Self {
        self.observation_sensor = Some(source_key(sensor));
        self
    }

    fn with_duplicate_row_channel(mut self) -> Self {
        self.duplicate_row_channel = true;
        self
    }

    fn with_declared_duplicate_row_sharing(mut self) -> Self {
        self.duplicate_row_channel = true;
        self.declare_duplicate_row_sharing = true;
        self
    }

    fn without_blind_release(mut self) -> Self {
        self.blind_release = None;
        self
    }

    fn with_blind_release_authority(mut self, authority_label: &str) -> Self {
        self.blind_release = Some(blind_release_for(&self.split, authority_label));
        self
    }

    fn with_blind_release(mut self, release: BlindReleaseReceipt) -> Self {
        self.blind_release = Some(release);
        self
    }
}

fn case_fixture(
    id: &'static str,
    purpose: CasePurpose,
    observation_qoi: &'static str,
    observation_row: &'static str,
    experiment_key: &'static str,
    split_key: &'static str,
    data: (ExperimentArtifact, CalibrationSplit),
) -> RetrospectiveCaseFixture {
    let blind_release = matches!(&purpose, CasePurpose::BlindFalsification)
        .then(|| blind_release_for(&data.1, &format!("blind-release-authority-{id}")));
    RetrospectiveCaseFixture {
        id,
        purpose,
        observation_qoi,
        observation_row,
        experiment_key,
        split_key,
        experiment: data.0,
        split: data.1,
        observation_instrument: None,
        observation_acquisition_channel: None,
        observation_clock: None,
        observation_sensor: None,
        duplicate_row_channel: false,
        declare_duplicate_row_sharing: false,
        blind_release,
    }
}

fn case_sensor_source(case: &RetrospectiveCaseFixture) -> SourceKey {
    if let Some(sensor) = &case.observation_sensor {
        return sensor.clone();
    }
    let instrument = case
        .experiment
        .instruments()
        .first()
        .expect("retrospective fixture instrument")
        .instrument_id()
        .as_str();
    let suffix = instrument
        .strip_prefix("instrument-")
        .expect("fixture instrument naming contract");
    source_key(&format!("sensor-{suffix}"))
}

fn study_case(case: &RetrospectiveCaseFixture) -> StudyCaseDocument {
    let frame_transform = format!("frame-transform-{}", case.id);
    let frame = FrameBinding::try_new(
        artifact(&format!("frame-{}", case.id)),
        case_physics_hash(FRAME_TRANSFORM_SOURCE_DOMAIN, &frame_transform),
        "right-handed-cartesian",
    )
    .expect("frame fixture");
    let experiment_clock = match case.experiment.clocks() {
        ClockSynchronization::SingleClock { clock_id } => clock_id.clone(),
        ClockSynchronization::Synchronized { clock_ids, .. } => clock_ids
            .first()
            .expect("synchronized fixture clock")
            .clone(),
    };
    let instrument = case
        .experiment
        .instruments()
        .first()
        .expect("retrospective fixture instrument")
        .instrument_id()
        .clone();
    let observation_instrument = case
        .observation_instrument
        .clone()
        .unwrap_or_else(|| instrument.clone());
    let observation_acquisition_channel = case
        .observation_acquisition_channel
        .clone()
        .unwrap_or_else(|| artifact(&format!("acquisition-channel-{}", case.observation_row)));
    let observation_clock = case
        .observation_clock
        .clone()
        .unwrap_or_else(|| experiment_clock.clone());
    let load_path = format!("load-path-{}", case.id);
    let environment_path = format!("environment-path-{}", case.id);
    let time_grid = format!("time-grid-{}", case.id);
    let protocol = ProtocolBinding::try_new(
        artifact(&format!("protocol-{}", case.id)),
        7,
        2,
        3,
        case_physics_hash(LOAD_PATH_SOURCE_DOMAIN, &load_path),
        case_physics_hash(ENVIRONMENT_PATH_SOURCE_DOMAIN, &environment_path),
        case_physics_hash(TIME_GRID_SOURCE_DOMAIN, &time_grid),
        observation_clock.clone(),
    )
    .expect("protocol fixture");
    let observation_channel = channel(&format!("signal-{}", case.id));
    let observation = StudyObservation::try_new(
        observation_channel.clone(),
        qoi(case.observation_qoi),
        unit("Pa"),
        QuantitySpec::dimensional(STRESS),
        frame.clone(),
        format!("node-{}", case.id),
        "stress-output",
        source_key(&format!("operator-{}", case.id)),
        source_key(&format!("aggregation-{}", case.id)),
        case_sensor_source(case),
        observation_instrument,
        observation_acquisition_channel,
        observation_clock,
        4,
        MarginalNoiseSpec::Gaussian {
            standard_deviation: 2.0e5,
        },
        MissingnessAssumption::Unknown {
            reason: "missingness has not yet been characterized".to_string(),
        },
        None,
        7,
        3,
        ObservationRows::Retrospective(BTreeSet::from([observation(case.observation_row)])),
    )
    .expect("retrospective observation fixture");
    let mut observations = vec![observation];
    if case.duplicate_row_channel {
        let original = &observations[0];
        observations.push(
            StudyObservation::try_new(
                channel(&format!("signal-{}-duplicate", case.id)),
                original.qoi().clone(),
                original.unit().clone(),
                original.quantity(),
                original.frame().clone(),
                original.graph_node().to_string(),
                original.graph_port().to_string(),
                original.operator().clone(),
                original.aggregation().clone(),
                original.sensor().clone(),
                original.instrument().clone(),
                original.acquisition_channel().clone(),
                original.clock().clone(),
                original.operator_version(),
                original.noise().clone(),
                original.missingness().clone(),
                original.saturation(),
                original.protocol_version(),
                original.refinement_version(),
                original.rows().clone(),
            )
            .expect("duplicate-row observation fixture"),
        );
    }
    let discrepancies = observations
        .iter()
        .map(|observation| {
            (
                observation.id().clone(),
                StudyDiscrepancy::Uncharacterized {
                    reason: "no discrepancy model is admitted for this test channel".to_string(),
                },
            )
        })
        .collect();
    let observation_sharing = if case.declare_duplicate_row_sharing {
        vec![
            ObservationSharingGroup::try_new(
                BTreeSet::from([
                    observation_channel,
                    channel(&format!("signal-{}-duplicate", case.id)),
                ]),
                BTreeSet::from([observation(case.observation_row)]),
                source_key("joint-likelihood"),
                "the two channels intentionally share one row under the exact joint likelihood",
            )
            .expect("within-case sharing group"),
        ]
    } else {
        Vec::new()
    };
    let geometry = format!("geometry-{}", case.id);
    let process = format!("process-{}", case.id);
    let preparation = format!("preparation-{}", case.id);
    StudyCaseDocument::try_new(
        case_id(case.id),
        case.purpose.clone(),
        InitialStateBinding::Zero { schema_version: 2 },
        SpecimenBinding::try_new(
            artifact(&format!("specimen-{}", case.id)),
            case_physics_hash(SPECIMEN_GEOMETRY_SOURCE_DOMAIN, &geometry),
            case_physics_hash(SPECIMEN_PROCESS_SOURCE_DOMAIN, &process),
            case_physics_hash(SPECIMEN_PREPARATION_SOURCE_DOMAIN, &preparation),
            frame,
        )
        .expect("specimen fixture"),
        protocol,
        CasePhysicsSources::new(
            source_key(&frame_transform),
            source_key(&geometry),
            source_key(&process),
            source_key(&preparation),
            source_key(&load_path),
            source_key(&environment_path),
            source_key(&time_grid),
            None,
        ),
        source_key(&format!("forward-{}", case.id)),
        CaseDataDeclaration::Retrospective {
            experiment: source_key(case.experiment_key),
            split: source_key(case.split_key),
            parser: source_key("parser"),
            preprocessing: source_key("preprocessing"),
            parser_version: 2,
            split_grouping: artifact("split-by-specimen"),
        },
        observations,
        discrepancies,
        observation_sharing,
    )
    .expect("study case fixture")
}

struct ProblemFixture {
    context: ContextOfUse,
    material: MaterialCard,
    model: ConstitutiveModelCard,
    document: IdentifiabilityProblemDocument,
    cases: Vec<RetrospectiveCaseFixture>,
}

fn problem_fixture(
    cases: Vec<RetrospectiveCaseFixture>,
    data_reuse: DataReusePolicy,
) -> ProblemFixture {
    try_problem_fixture_with_global_likelihood(cases, data_reuse, "joint-likelihood")
        .expect("retrospective problem is structurally valid")
}

fn try_problem_fixture_with_global_likelihood(
    cases: Vec<RetrospectiveCaseFixture>,
    data_reuse: DataReusePolicy,
    global_likelihood: &str,
) -> Result<ProblemFixture, IdentifiabilityError> {
    let within_case_sharing = cases.iter().any(|case| case.declare_duplicate_row_sharing);
    let any_sharing = matches!(&data_reuse, DataReusePolicy::Shared { .. }) || within_case_sharing;
    let context = context();
    let (material, model) = model_cards();
    let mut sources = vec![
        source(
            "context",
            SourceKind::ContextOfUse,
            context.content_hash().expect("context hashes"),
        ),
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
        source("graph", SourceKind::ConstitutiveGraph, hash("graph")),
        source("parser", SourceKind::Parser, hash("parser")),
        source(
            "preprocessing",
            SourceKind::Preprocessing,
            hash("preprocessing"),
        ),
    ];
    let mut registered_sensor_sources = BTreeSet::new();
    for case in &cases {
        sources.extend([
            source(
                &format!("forward-{}", case.id),
                SourceKind::ForwardModel,
                hash(&format!("forward-{}", case.id)),
            ),
            source(
                &format!("operator-{}", case.id),
                SourceKind::ObservationOperator,
                hash(&format!("operator-{}", case.id)),
            ),
            source(
                &format!("aggregation-{}", case.id),
                SourceKind::ObservationOperator,
                hash(&format!("aggregation-{}", case.id)),
            ),
            case_physics_source(
                &format!("frame-transform-{}", case.id),
                SourceKind::Geometry,
                FRAME_TRANSFORM_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("geometry-{}", case.id),
                SourceKind::Geometry,
                SPECIMEN_GEOMETRY_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("process-{}", case.id),
                SourceKind::Process,
                SPECIMEN_PROCESS_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("preparation-{}", case.id),
                SourceKind::Process,
                SPECIMEN_PREPARATION_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("load-path-{}", case.id),
                SourceKind::Protocol,
                LOAD_PATH_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("environment-path-{}", case.id),
                SourceKind::Protocol,
                ENVIRONMENT_PATH_SOURCE_DOMAIN,
            ),
            case_physics_source(
                &format!("time-grid-{}", case.id),
                SourceKind::Protocol,
                TIME_GRID_SOURCE_DOMAIN,
            ),
        ]);
        for candidate in [
            source(
                case.experiment_key,
                SourceKind::ExperimentArtifact,
                case.experiment
                    .content_hash()
                    .expect("experiment source hashes"),
            ),
            source(
                case.split_key,
                SourceKind::CalibrationSplit,
                case.split.content_hash().expect("split source hashes"),
            ),
        ] {
            if let Some(existing) = sources
                .iter()
                .find(|existing| existing.key() == candidate.key())
            {
                assert_eq!(
                    existing, &candidate,
                    "a shared concrete source key must retain exactly one source reference"
                );
            } else {
                sources.push(candidate);
            }
        }
        let sensor = case_sensor_source(case);
        if registered_sensor_sources.insert(sensor.clone()) {
            let certificate_hash = case
                .experiment
                .instruments()
                .first()
                .expect("fixture experiment instrument")
                .certificate_hash();
            let expected_hash = if case.observation_sensor.is_some() {
                hash(sensor.as_str())
            } else {
                certificate_hash
            };
            sources.push(source(
                sensor.as_str(),
                SourceKind::Metrology,
                expected_hash,
            ));
        }
    }
    if any_sharing {
        let mut likelihoods = BTreeSet::from([source_key(global_likelihood)]);
        if within_case_sharing {
            likelihoods.insert(source_key("joint-likelihood"));
        }
        if let DataReusePolicy::Shared { groups } = &data_reuse {
            likelihoods.extend(groups.iter().map(|group| group.joint_likelihood().clone()));
        }
        sources.extend(likelihoods.into_iter().map(|likelihood| {
            source(
                likelihood.as_str(),
                SourceKind::Likelihood,
                hash(likelihood.as_str()),
            )
        }));
    }
    let parameter_domain =
        ParameterDomain::try_new(1.0e6, 1.0e9).expect("yield-stress parameter domain");
    let parameters = vec![
        StudyParameter::try_new(
            role("yield_stress"),
            QuantitySpec::dimensional(STRESS),
            parameter_domain,
            ParameterPurpose::Estimand,
            ParameterTreatment::Estimated,
            ParameterOwnerBinding::ConstitutiveModel,
            ParameterScopeBinding::Global,
            PriorPolicy::Distribution(ParameterPrior::Uniform {
                version: 1,
                domain: parameter_domain,
            }),
            InfluenceCoverage::Declared,
        )
        .expect("study parameter fixture"),
    ];
    let influences = cases
        .iter()
        .map(|case| {
            InfluenceDeclaration::new(
                InfluenceId::try_new(format!("yield-to-observation-{}", case.id))
                    .expect("influence id"),
                role("yield_stress"),
                DistributionFunctional::Location {
                    observation: ObservationKey::new(
                        case_id(case.id),
                        channel(&format!("signal-{}", case.id)),
                    ),
                },
                InfluenceRepresentation::Direct,
            )
        })
        .collect();
    let joint_noise = if any_sharing {
        // This fixture conservatively declines an independence assumption for
        // reused provenance. The sharing group's likelihood is also the
        // explicit cross-case noise kernel; provenance reuse alone is not a
        // theorem of stochastic dependence.
        JointNoiseModel::ExternalKernel {
            model: source_key(global_likelihood),
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
    let admissible_domain = AdmissibleDomainWitness::try_new(
        parameters
            .iter()
            .map(|parameter| (parameter.role().clone(), parameter.domain().bounds().0))
            .collect(),
        None,
    )
    .expect("constructive admissible-domain witness");
    let document = IdentifiabilityProblemDocument::try_new(
        source_key("context"),
        source_key("material"),
        source_key("model"),
        source_key("graph"),
        sources,
        parameters,
        Vec::new(),
        admissible_domain,
        cases.iter().map(study_case).collect(),
        influences,
        Vec::new(),
        joint_noise,
        data_reuse,
    )?;
    Ok(ProblemFixture {
        context,
        material,
        model,
        document,
        cases,
    })
}

fn opaque_resolutions(document: &IdentifiabilityProblemDocument) -> SourceResolutionSet {
    SourceResolutionSet::try_new(
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
                    source.key().as_str().as_bytes(),
                    AuthorityDisposition::ContentVerified,
                )
                .expect("opaque source resolution")
            })
            .collect(),
    )
    .expect("closed opaque source resolution set")
}

#[derive(Debug, Clone, Copy)]
enum BundleMode {
    Exact,
    Missing,
    Extra,
}

fn admit(
    fixture: ProblemFixture,
    bundle_mode: BundleMode,
) -> Result<AdmittedIdentifiabilityProblem, IdentifiabilityError> {
    admit_with_concrete_authority(fixture, bundle_mode, Vec::new())
}

fn admit_with_concrete_authority(
    fixture: ProblemFixture,
    bundle_mode: BundleMode,
    concrete_authority: Vec<(SourceKey, AuthorityDisposition)>,
) -> Result<AdmittedIdentifiabilityProblem, IdentifiabilityError> {
    let ProblemFixture {
        context,
        material,
        model,
        document,
        cases,
    } = fixture;
    let opaque = opaque_resolutions(&document);
    let mut bundles = BTreeMap::new();
    if !matches!(bundle_mode, BundleMode::Missing) {
        for case in &cases {
            let mut bundle = CaseSourceBundle::new(&case.experiment, &case.split);
            if let Some(release) = &case.blind_release {
                bundle = bundle.with_blind_release(release);
            }
            bundles.insert(case_id(case.id), bundle);
        }
    }
    if matches!(bundle_mode, BundleMode::Extra) {
        let source = cases.first().expect("extra-bundle fixture needs one case");
        bundles.insert(
            case_id("unknown-case"),
            CaseSourceBundle::new(&source.experiment, &source.split),
        );
    }
    let bundle = ProblemSourceBundle::new(&context, &material, &model, bundles, opaque)
        .with_concrete_authority(concrete_authority)?;
    AdmittedIdentifiabilityProblem::resolve_and_admit(document, bundle)
}

fn ordinary_data(label: &str) -> (ExperimentArtifact, CalibrationSplit) {
    let calibration = format!("cal-{label}");
    let validation = format!("val-{label}");
    let blind = format!("blind-{label}");
    let calibration_source = format!("source-cal-{label}");
    let validation_source = format!("source-val-{label}");
    let blind_source = format!("source-blind-{label}");
    let source_bytes = format!("source-bytes-{label}");
    physical_data(
        label,
        &["stress"],
        &[
            (calibration.as_str(), calibration_source.as_str()),
            (validation.as_str(), validation_source.as_str()),
            (blind.as_str(), blind_source.as_str()),
        ],
        &[calibration.as_str()],
        &[validation.as_str()],
        &[blind.as_str()],
        &source_bytes,
    )
}

/// A two-QoI, two-instrument, two-clock campaign whose primary and alternative
/// rows are both globally admissible.  Tests can therefore cross-wire one
/// row's interpretation without relying on an obviously unknown identifier.
fn cross_wire_data(
    swap_primary_and_alternative_sources: bool,
) -> (ExperimentArtifact, CalibrationSplit) {
    let dataset_source_bytes_hash = hash("source-bytes-cross-wire");
    let primary_instrument = artifact("instrument-cross-wire-a");
    let alternative_instrument = artifact("instrument-cross-wire-z");
    let primary_clock = artifact("clock-cross-wire-a");
    let alternative_clock = artifact("clock-cross-wire-z");
    let primary_source = ObservationSourceRef::try_new(
        dataset_source_bytes_hash,
        TEST_ROW_LOCATOR_DOMAIN,
        1,
        hash("locator-cross-wire-primary"),
        hash("extraction-cross-wire-primary"),
    )
    .expect("primary typed row source");
    let alternative_source = ObservationSourceRef::try_new(
        dataset_source_bytes_hash,
        TEST_ALTERNATIVE_ROW_LOCATOR_DOMAIN,
        2,
        hash("locator-cross-wire-alternative"),
        hash("extraction-cross-wire-alternative"),
    )
    .expect("alternative typed row source");
    let (primary_source, alternative_source) = if swap_primary_and_alternative_sources {
        (alternative_source, primary_source)
    } else {
        (primary_source, alternative_source)
    };
    let manifest = ObservationManifest::try_new(vec![
        (
            observation("cal-cross-primary"),
            ObservationManifestRow::try_new(
                primary_source,
                qoi("stress"),
                primary_instrument.clone(),
                artifact("acquisition-channel-cross-primary"),
                primary_clock.clone(),
            )
            .expect("primary manifest row"),
        ),
        (
            observation("cal-cross-alternative"),
            ObservationManifestRow::try_new(
                alternative_source,
                qoi("tangent"),
                alternative_instrument.clone(),
                artifact("acquisition-channel-cross-alternative"),
                alternative_clock.clone(),
            )
            .expect("alternative manifest row"),
        ),
        (
            observation("val-cross-wire"),
            ObservationManifestRow::try_new(
                ObservationSourceRef::try_new(
                    dataset_source_bytes_hash,
                    TEST_ROW_LOCATOR_DOMAIN,
                    1,
                    hash("locator-cross-wire-validation"),
                    hash("extraction-cross-wire-validation"),
                )
                .expect("validation typed row source"),
                qoi("stress"),
                primary_instrument.clone(),
                artifact("acquisition-channel-cross-validation"),
                primary_clock.clone(),
            )
            .expect("validation manifest row"),
        ),
        (
            observation("blind-cross-wire"),
            ObservationManifestRow::try_new(
                ObservationSourceRef::try_new(
                    dataset_source_bytes_hash,
                    TEST_ROW_LOCATOR_DOMAIN,
                    1,
                    hash("locator-cross-wire-blind"),
                    hash("extraction-cross-wire-blind"),
                )
                .expect("blind typed row source"),
                qoi("stress"),
                primary_instrument.clone(),
                artifact("acquisition-channel-cross-blind"),
                primary_clock.clone(),
            )
            .expect("blind manifest row"),
        ),
    ])
    .expect("injective cross-wire manifest");
    let experiment = ExperimentArtifact::try_new(
        header("experiment-cross-wire", &["Pa"], "fixture.experiment"),
        artifact("dataset-cross-wire"),
        ExperimentOrigin::Physical {
            apparatus_id: artifact("apparatus-cross-wire"),
            facility_id: artifact("facility-retrospective-fixture"),
        },
        vec![qoi("stress"), qoi("tangent")],
        manifest,
        vec![
            InstrumentCalibration::new(primary_instrument, hash("sensor-cross-wire-a"), true),
            InstrumentCalibration::new(alternative_instrument, hash("sensor-cross-wire-z"), true),
        ],
        ClockSynchronization::synchronized(
            vec![primary_clock, alternative_clock],
            "fixture cross-clock synchronization",
            1.0e-9,
            hash("cross-clock-synchronization-evidence"),
        )
        .expect("cross-wire clock topology"),
        RepeatabilitySummary::try_new(
            3,
            CovarianceMatrix::try_new(2, vec![0.25, 0.0, 0.25]).expect("cross-wire covariance"),
        )
        .expect("cross-wire repeatability"),
        DataAuthenticity::new(dataset_source_bytes_hash, hash("custody-cross-wire"), true),
    )
    .expect("cross-wire experiment");
    let split = CalibrationSplit::try_new(
        header("split-cross-wire", &["unitless"], "fixture.split"),
        experiment_reference(&experiment),
        hash("preregistration-cross-wire"),
        vec![
            observation("cal-cross-primary"),
            observation("cal-cross-alternative"),
        ],
        vec![observation("val-cross-wire")],
        vec![(
            observation("blind-cross-wire"),
            hash("locator-cross-wire-blind"),
        )],
    )
    .expect("cross-wire split");
    (experiment, split)
}

fn replacement_split(
    label: &str,
    experiment_reference: ArtifactRef,
    calibration: &[&str],
    validation: &[&str],
    blind: &[(&str, ContentHash)],
) -> CalibrationSplit {
    CalibrationSplit::try_new(
        header(
            &format!("replacement-split-{label}"),
            &["unitless"],
            "fixture.split",
        ),
        experiment_reference,
        hash(&format!("replacement-preregistration-{label}")),
        calibration.iter().copied().map(observation).collect(),
        validation.iter().copied().map(observation).collect(),
        blind
            .iter()
            .map(|(id, source)| (observation(id), *source))
            .collect(),
    )
    .expect("replacement split is structurally valid")
}

fn experiment_reference(experiment: &ExperimentArtifact) -> ArtifactRef {
    ArtifactRef::new(
        ArtifactKind::ExperimentArtifact,
        experiment.id().clone(),
        experiment.content_hash().expect("experiment content hash"),
    )
}

fn assert_manifest_cross_wire_refuses(name: &str, case: RetrospectiveCaseFixture) {
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("a globally valid value from another manifest row must not be cross-wired");
    assert!(matches!(
        &error,
        IdentifiabilityError::SourceMismatch {
            field: "observation/manifest row binding",
        }
    ));
    log(
        name,
        "pass",
        "each consumed row retains its own exact scientific and metrology interpretation",
        &error.to_string(),
    );
}

#[test]
fn calibration_case_admits_only_its_calibration_partition() {
    let fixture = problem_fixture(
        vec![case_fixture(
            "a",
            CasePurpose::Calibration,
            "stress",
            "cal-a",
            "experiment-a",
            "split-a",
            ordinary_data("a"),
        )],
        DataReusePolicy::Disjoint,
    );
    let admitted = admit(fixture, BundleMode::Exact).expect("calibration case admits");
    assert_eq!(admitted.data().len(), 1);
    log(
        "calibration-partition-admission",
        "pass",
        "one source-resolved calibration lineage",
        &format!(
            "admitted_cases={}; problem_id={:?}",
            admitted.data().len(),
            admitted.id().digest(),
        ),
    );
}

#[test]
fn public_lineage_accessors_do_not_release_validation_or_blind_rows() {
    let fixture = problem_fixture(
        vec![case_fixture(
            "a",
            CasePurpose::Calibration,
            "stress",
            "cal-a",
            "experiment-a",
            "split-a",
            ordinary_data("a"),
        )],
        DataReusePolicy::Disjoint,
    );
    let admitted = admit(fixture, BundleMode::Exact).expect("calibration case admits");
    let lineage = &admitted.data()[&case_id("a")];
    let calibration = observation("cal-a");
    let validation = observation("val-a");
    let blind = observation("blind-a");
    let unknown = observation("not-in-the-manifest");

    assert_eq!(
        lineage.calibration_ids(),
        &BTreeSet::from([calibration.clone()])
    );
    assert_eq!(lineage.partition_counts(), (1, 1, 1));
    assert!(
        lineage
            .row_binding(&calibration)
            .map(ObservationManifestRow::locator_hash)
            .is_some()
    );
    for protected in [&validation, &blind, &unknown] {
        assert_eq!(lineage.row_binding(protected), None);
    }
    log(
        "public-lineage-capability-boundary",
        "pass",
        "only calibration-authorized row identities and bindings are exposed",
        "validation, blind, and unknown row accessors returned None; aggregate counts remain visible",
    );
}

#[test]
fn data_lineage_debug_redacts_row_and_source_capabilities() {
    fn push_hash_fragments(fragments: &mut Vec<String>, hash: ContentHash) {
        fragments.push(hash.to_string());
        fragments.push(format!("{hash:?}"));
    }

    let (experiment, split) = ordinary_data("a");
    let fixture = problem_fixture(
        vec![case_fixture(
            "a",
            CasePurpose::Calibration,
            "stress",
            "cal-a",
            "experiment-a",
            "split-a",
            (experiment.clone(), split.clone()),
        )],
        DataReusePolicy::Disjoint,
    );
    let admitted = admit(fixture, BundleMode::Exact).expect("calibration case admits");
    let lineage = &admitted.data()[&case_id("a")];
    let calibration = observation("cal-a");
    let debug = format!("{lineage:?}");

    let mut protected_fragments = vec![
        experiment.id().as_str().to_string(),
        experiment.dataset_id().as_str().to_string(),
        split.id().as_str().to_string(),
        split.experiment().id().as_str().to_string(),
        lineage.experiment().id().as_str().to_string(),
        lineage.split().id().as_str().to_string(),
        lineage.split_grouping().as_str().to_string(),
    ];
    for hash in [
        lineage.experiment().hash(),
        lineage.split().hash(),
        lineage.raw_manifest(),
        lineage.source_bytes(),
        lineage.custody_receipt(),
        lineage.preregistration(),
        lineage.parser(),
        lineage.preprocessing(),
        lineage.blind_commitment(),
        lineage
            .row_binding(&calibration)
            .expect("calibration locator remains available through its capability")
            .locator_hash(),
    ] {
        push_hash_fragments(&mut protected_fragments, hash);
    }
    for (row_id, row) in experiment.manifest().rows() {
        protected_fragments.extend([
            row_id.as_str().to_string(),
            row.qoi().as_str().to_string(),
            row.instrument().as_str().to_string(),
            row.acquisition_channel().as_str().to_string(),
            row.clock().as_str().to_string(),
            row.source_ref().locator_domain().to_string(),
        ]);
        for hash in [
            row.source_ref().dataset_source_bytes_hash(),
            row.source_ref().locator_hash(),
            row.source_ref().extraction_receipt_hash(),
        ] {
            push_hash_fragments(&mut protected_fragments, hash);
        }
    }
    for instrument in experiment.instruments() {
        protected_fragments.push(instrument.instrument_id().as_str().to_string());
        push_hash_fragments(&mut protected_fragments, instrument.certificate_hash());
    }
    for hash in [
        experiment.authenticity().source_bytes_hash(),
        experiment.authenticity().custody_receipt_hash(),
        split.experiment().hash(),
        split.preregistration_hash(),
        split.blind_commitment(),
    ] {
        push_hash_fragments(&mut protected_fragments, hash);
    }
    protected_fragments.extend(
        split
            .calibration_ids()
            .iter()
            .chain(split.validation_ids())
            .map(|id| id.as_str().to_string()),
    );
    for (id, source) in split.blind_sources() {
        protected_fragments.push(id.as_str().to_string());
        push_hash_fragments(&mut protected_fragments, *source);
    }
    protected_fragments.sort();
    protected_fragments.dedup();

    for protected_identity in &protected_fragments {
        assert!(
            !debug.contains(protected_identity),
            "Debug leaked `{protected_identity}` from source, custody, metrology, partition, or row authority: {debug}",
        );
    }
    assert!(debug.len() < 1_024, "Debug must remain bounded: {debug}");
    assert!(debug.contains("experiment_binding: \"<redacted>\""));
    assert!(debug.contains("split_binding: \"<redacted>\""));
    assert!(debug.contains("blind_commitment: \"<redacted>\""));
    assert!(debug.contains("row_bindings: \"<redacted>\""));
    assert!(debug.contains("source_and_custody_hashes: \"<redacted>\""));
    log(
        "lineage-debug-capability-boundary",
        "pass",
        "diagnostics retain aggregate counts without leaking row/source/provenance capabilities",
        &debug,
    );
}

#[test]
fn retrospective_transport_is_v3_fixed_point_and_refuses_pre_v3_bytes() {
    assert_eq!(VV_SCHEMA_VERSION, 3, "this fixture exercises V&V schema v3");
    let (experiment, _) = ordinary_data("transport-v3");
    let bytes = experiment
        .canonical_bytes()
        .expect("canonical v3 experiment transport");
    assert_eq!(&bytes[4..8], &VV_SCHEMA_VERSION.to_le_bytes());
    let round_trip =
        ExperimentArtifact::from_canonical_bytes(&bytes).expect("current v3 transport decodes");
    assert_eq!(round_trip, experiment);
    assert_eq!(
        round_trip
            .canonical_bytes()
            .expect("round-trip canonical bytes"),
        bytes,
    );
    for stale_schema in [1_u32, 2] {
        let mut stale = bytes.clone();
        stale[4..8].copy_from_slice(&stale_schema.to_le_bytes());
        let error = ExperimentArtifact::from_canonical_bytes(&stale)
            .expect_err("pre-v3 observation-manifest transport must refuse");
        assert_eq!(error.rule_name(), "vv-canonical-identity");
        assert_eq!(error.offset(), 4);
        assert!(error.detail().contains("unsupported V&V schema version"));
    }
    log(
        "retrospective-v3-transport",
        "pass",
        "schema-v3 fixed point and fail-closed v1/v2 transport refusal",
        &format!("canonical_bytes={}", bytes.len()),
    );
}

#[test]
fn validation_only_case_admits_only_its_validation_partition() {
    let fixture = problem_fixture(
        vec![case_fixture(
            "a",
            CasePurpose::ValidationOnly,
            "stress",
            "val-a",
            "experiment-a",
            "split-a",
            ordinary_data("a"),
        )],
        DataReusePolicy::Disjoint,
    );
    let admitted = admit(fixture, BundleMode::Exact).expect("validation-only case admits");
    assert_eq!(admitted.data().len(), 1);
    log(
        "validation-partition-admission",
        "pass",
        "one source-resolved validation lineage",
        &format!("admitted_cases={}", admitted.data().len()),
    );
}

#[test]
fn blind_falsification_case_admits_only_its_blind_partition() {
    let fixture = problem_fixture(
        vec![case_fixture(
            "a",
            CasePurpose::BlindFalsification,
            "stress",
            "blind-a",
            "experiment-a",
            "split-a",
            ordinary_data("a"),
        )],
        DataReusePolicy::Disjoint,
    );
    let admitted = admit(fixture, BundleMode::Exact).expect("blind-falsification case admits");
    assert_eq!(admitted.data().len(), 1);
    log(
        "blind-partition-admission",
        "pass",
        "one source-resolved blind lineage",
        &format!("admitted_cases={}", admitted.data().len()),
    );
}

#[test]
fn blind_falsification_without_release_refuses() {
    let case = case_fixture(
        "a",
        CasePurpose::BlindFalsification,
        "stress",
        "blind-a",
        "experiment-a",
        "split-a",
        ordinary_data("a"),
    )
    .without_blind_release();
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("sealed blind rows require a release");
    assert!(matches!(
        &error,
        IdentifiabilityError::InvalidText {
            field: "blind release",
            ..
        }
    ));
    log(
        "blind-release-required",
        "pass",
        "blind holdout remains sealed without authority",
        &error.to_string(),
    );
}

#[test]
fn blind_release_bound_to_another_split_refuses() {
    let data = ordinary_data("a");
    let wrong_release = BlindReleaseReceipt::new(
        ArtifactRef::new(
            ArtifactKind::CalibrationSplit,
            artifact("another-split"),
            hash("another-split"),
        ),
        data.1.blind_commitment(),
        hash("wrong-split-release-authority"),
    )
    .expect("structurally valid wrong-split receipt");
    let case = case_fixture(
        "a",
        CasePurpose::BlindFalsification,
        "stress",
        "blind-a",
        "experiment-a",
        "split-a",
        data,
    )
    .with_blind_release(wrong_release);
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("release for another split must refuse");
    assert!(matches!(&error, IdentifiabilityError::Vv { .. }));
    log(
        "blind-release-split-binding",
        "pass",
        "release must bind the exact split id and content hash",
        &error.to_string(),
    );
}

#[test]
fn blind_release_with_wrong_split_id_refuses() {
    let data = ordinary_data("a");
    let wrong_release = BlindReleaseReceipt::new(
        ArtifactRef::new(
            ArtifactKind::CalibrationSplit,
            artifact("another-split"),
            data.1.content_hash().expect("real split hash"),
        ),
        data.1.blind_commitment(),
        hash("wrong-split-id-release-authority"),
    )
    .expect("structurally valid wrong-id receipt");
    let case = case_fixture(
        "a",
        CasePurpose::BlindFalsification,
        "stress",
        "blind-a",
        "experiment-a",
        "split-a",
        data,
    )
    .with_blind_release(wrong_release);
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("right hash under another split id must refuse");
    assert!(matches!(&error, IdentifiabilityError::Vv { .. }));
    log(
        "blind-release-split-id-binding",
        "pass",
        "release binds split id independently of the split hash",
        &error.to_string(),
    );
}

#[test]
fn blind_release_with_wrong_split_hash_refuses() {
    let data = ordinary_data("a");
    let wrong_release = BlindReleaseReceipt::new(
        ArtifactRef::new(
            ArtifactKind::CalibrationSplit,
            data.1.id().clone(),
            hash("another-split-hash"),
        ),
        data.1.blind_commitment(),
        hash("wrong-split-hash-release-authority"),
    )
    .expect("structurally valid wrong-hash receipt");
    let case = case_fixture(
        "a",
        CasePurpose::BlindFalsification,
        "stress",
        "blind-a",
        "experiment-a",
        "split-a",
        data,
    )
    .with_blind_release(wrong_release);
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("right split id under another hash must refuse");
    assert!(matches!(&error, IdentifiabilityError::Vv { .. }));
    log(
        "blind-release-split-hash-binding",
        "pass",
        "release binds split content hash independently of the split id",
        &error.to_string(),
    );
}

#[test]
fn blind_release_with_wrong_commitment_refuses() {
    let data = ordinary_data("a");
    let wrong_release = BlindReleaseReceipt::new(
        ArtifactRef::new(
            ArtifactKind::CalibrationSplit,
            data.1.id().clone(),
            data.1.content_hash().expect("split hash"),
        ),
        hash("wrong-blind-commitment"),
        hash("wrong-commitment-release-authority"),
    )
    .expect("structurally valid wrong-commitment receipt");
    let case = case_fixture(
        "a",
        CasePurpose::BlindFalsification,
        "stress",
        "blind-a",
        "experiment-a",
        "split-a",
        data,
    )
    .with_blind_release(wrong_release);
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("release for another commitment must refuse");
    assert!(matches!(&error, IdentifiabilityError::Vv { .. }));
    log(
        "blind-release-commitment-binding",
        "pass",
        "release cannot be replayed against a different sealed row commitment",
        &error.to_string(),
    );
}

#[test]
fn nonblind_case_rejects_surplus_blind_release() {
    let case = case_fixture(
        "a",
        CasePurpose::Calibration,
        "stress",
        "cal-a",
        "experiment-a",
        "split-a",
        ordinary_data("a"),
    )
    .with_blind_release_authority("surplus-release-authority");
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("non-blind case must not consume blind authority");
    assert!(matches!(
        &error,
        IdentifiabilityError::InvalidText {
            field: "blind release",
            ..
        }
    ));
    log(
        "blind-release-surplus",
        "pass",
        "blind-release authority is purpose-scoped and cannot be laundered into calibration",
        &error.to_string(),
    );
}

#[test]
fn blind_release_authority_moves_admission_not_problem_identity() {
    let make_fixture = |authority: &'static str| {
        problem_fixture(
            vec![
                case_fixture(
                    "a",
                    CasePurpose::BlindFalsification,
                    "stress",
                    "blind-a",
                    "experiment-a",
                    "split-a",
                    ordinary_data("a"),
                )
                .with_blind_release_authority(authority),
            ],
            DataReusePolicy::Disjoint,
        )
    };
    let left = admit(make_fixture("blind-authority-left"), BundleMode::Exact)
        .expect("left release admits");
    let right = admit(make_fixture("blind-authority-right"), BundleMode::Exact)
        .expect("right release admits");
    assert_eq!(left.id(), right.id());
    assert_ne!(left.source_admission_id(), right.source_admission_id());
    assert_ne!(
        left.source_admission_canonical_bytes()
            .expect("left source admission"),
        right
            .source_admission_canonical_bytes()
            .expect("right source admission"),
    );
    log(
        "blind-release-authority-identity",
        "pass",
        "release authority moves SourceAdmissionId while leaving the physical question stable",
        "two exact authority receipts retained in distinct admission preimages",
    );
}

#[test]
fn blind_release_authority_must_agree_with_explicit_concrete_authority() {
    let make_fixture = || {
        problem_fixture(
            vec![case_fixture(
                "a",
                CasePurpose::BlindFalsification,
                "stress",
                "blind-a",
                "experiment-a",
                "split-a",
                ordinary_data("a"),
            )],
            DataReusePolicy::Disjoint,
        )
    };
    let matching = AuthorityDisposition::ExternalTrustReceipt {
        trust_receipt: hash("blind-release-authority-a"),
    };
    admit_with_concrete_authority(
        make_fixture(),
        BundleMode::Exact,
        vec![(source_key("split-a"), matching)],
    )
    .expect("matching explicit split authority admits");

    let conflicting = AuthorityDisposition::ExternalTrustReceipt {
        trust_receipt: hash("different-explicit-split-authority"),
    };
    let error = admit_with_concrete_authority(
        make_fixture(),
        BundleMode::Exact,
        vec![(source_key("split-a"), conflicting)],
    )
    .expect_err("conflicting explicit split authority must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::SourceMismatch {
            field: "blind release/concrete source authority",
        }
    ));
    log(
        "blind-release-explicit-authority-agreement",
        "pass",
        "release-derived and caller-declared split authority must agree exactly",
        &error.to_string(),
    );
}

#[test]
fn shared_split_with_matching_blind_release_admits() {
    let shared = ordinary_data("shared-blind-positive");
    let release = blind_release_for(&shared.1, "shared-blind-authority");
    let cases = vec![
        case_fixture(
            "a",
            CasePurpose::BlindFalsification,
            "stress",
            "blind-shared-blind-positive",
            "experiment-shared",
            "split-shared",
            shared.clone(),
        )
        .with_blind_release(release.clone()),
        case_fixture(
            "b",
            CasePurpose::BlindFalsification,
            "stress",
            "blind-shared-blind-positive",
            "experiment-shared",
            "split-shared",
            shared,
        )
        .with_blind_release(release),
    ];
    let reuse = DataReusePolicy::Shared {
        groups: vec![
            DataSharingGroup::try_new(
                BTreeSet::from([case_id("a"), case_id("b")]),
                source_key("joint-likelihood"),
                "both cases intentionally consume one released sealed campaign",
            )
            .expect("sharing group"),
        ],
    };
    let admitted = admit(problem_fixture(cases, reuse), BundleMode::Exact)
        .expect("one exact shared release admits deterministically");
    assert_eq!(admitted.data().len(), 2);
    log(
        "shared-split-release-agreement",
        "pass",
        "one exact release authority can safely authorize a shared split",
        &format!("admitted_cases={}", admitted.data().len()),
    );
}

#[test]
fn shared_split_with_conflicting_blind_releases_refuses() {
    let shared = ordinary_data("shared-blind");
    let cases = vec![
        case_fixture(
            "a",
            CasePurpose::BlindFalsification,
            "stress",
            "blind-shared-blind",
            "experiment-shared",
            "split-shared",
            shared.clone(),
        )
        .with_blind_release_authority("shared-authority-left"),
        case_fixture(
            "b",
            CasePurpose::BlindFalsification,
            "stress",
            "blind-shared-blind",
            "experiment-shared",
            "split-shared",
            shared,
        )
        .with_blind_release_authority("shared-authority-right"),
    ];
    let reuse = DataReusePolicy::Shared {
        groups: vec![
            DataSharingGroup::try_new(
                BTreeSet::from([case_id("a"), case_id("b")]),
                source_key("joint-likelihood"),
                "both cases intentionally consume the same sealed campaign",
            )
            .expect("sharing group"),
        ],
    };
    let error = admit(problem_fixture(cases, reuse), BundleMode::Exact)
        .expect_err("one split key cannot carry contradictory release authority");
    assert!(matches!(
        &error,
        IdentifiabilityError::SourceMismatch {
            field: "shared split blind release",
        }
    ));
    log(
        "shared-split-release-conflict",
        "pass",
        "shared split authority must be exact and order-independent",
        &error.to_string(),
    );
}

#[test]
fn split_bound_to_another_experiment_refuses() {
    let (experiment, _) = ordinary_data("a");
    let split = replacement_split(
        "wrong-experiment",
        ArtifactRef::new(
            ArtifactKind::ExperimentArtifact,
            artifact("another-experiment"),
            hash("another-experiment"),
        ),
        &["cal-a"],
        &["val-a"],
        &[("blind-a", hash("source-blind-a"))],
    );
    let error = admit(
        problem_fixture(
            vec![case_fixture(
                "a",
                CasePurpose::Calibration,
                "stress",
                "cal-a",
                "experiment-a",
                "split-a",
                (experiment, split),
            )],
            DataReusePolicy::Disjoint,
        ),
        BundleMode::Exact,
    )
    .expect_err("split bound to another experiment must refuse");
    assert!(matches!(&error, IdentifiabilityError::Vv { .. }));
    log(
        "split-experiment-binding",
        "pass",
        "CalibrationSplit must bind the exact admitted ExperimentArtifact",
        &error.to_string(),
    );
}

#[test]
fn split_partition_union_different_from_manifest_refuses() {
    let (experiment, _) = ordinary_data("a");
    let split = replacement_split(
        "partition-union",
        experiment_reference(&experiment),
        &["cal-a"],
        &["val-not-in-manifest"],
        &[("blind-a", hash("source-blind-a"))],
    );
    let error = admit(
        problem_fixture(
            vec![case_fixture(
                "a",
                CasePurpose::Calibration,
                "stress",
                "cal-a",
                "experiment-a",
                "split-a",
                (experiment, split),
            )],
            DataReusePolicy::Disjoint,
        ),
        BundleMode::Exact,
    )
    .expect_err("split partition union different from the manifest must refuse");
    assert!(matches!(&error, IdentifiabilityError::Vv { .. }));
    log(
        "split-manifest-partition-union",
        "pass",
        "calibration, validation, and blind IDs must exactly cover the manifest",
        &error.to_string(),
    );
}

#[test]
fn split_blind_source_different_from_manifest_refuses() {
    let (experiment, _) = ordinary_data("a");
    let split = replacement_split(
        "blind-source",
        experiment_reference(&experiment),
        &["cal-a"],
        &["val-a"],
        &[("blind-a", hash("wrong-blind-source"))],
    );
    let error = admit(
        problem_fixture(
            vec![case_fixture(
                "a",
                CasePurpose::Calibration,
                "stress",
                "cal-a",
                "experiment-a",
                "split-a",
                (experiment, split),
            )],
            DataReusePolicy::Disjoint,
        ),
        BundleMode::Exact,
    )
    .expect_err("blind row rebound to another immutable source must refuse");
    assert!(matches!(&error, IdentifiabilityError::Vv { .. }));
    log(
        "blind-row-source-binding",
        "pass",
        "blind row identity and immutable source identity remain jointly sealed",
        &error.to_string(),
    );
}

#[test]
fn observation_row_absent_from_manifest_refuses() {
    let error = admit(
        problem_fixture(
            vec![case_fixture(
                "a",
                CasePurpose::Calibration,
                "stress",
                "cal-not-in-manifest",
                "experiment-a",
                "split-a",
                ordinary_data("a"),
            )],
            DataReusePolicy::Disjoint,
        ),
        BundleMode::Exact,
    )
    .expect_err("observation row outside the admitted manifest must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::UnknownReference {
            field: "observation raw row",
            ..
        }
    ));
    log(
        "observation-manifest-row-closure",
        "pass",
        "each retrospective observation row must occur in the admitted manifest",
        &error.to_string(),
    );
}

fn assert_case_partition_refuses(
    name: &str,
    purpose: CasePurpose,
    row: &'static str,
    expected_partition: &str,
) {
    let fixture = problem_fixture(
        vec![case_fixture(
            "a",
            purpose,
            "stress",
            row,
            "experiment-a",
            "split-a",
            ordinary_data("a"),
        )],
        DataReusePolicy::Disjoint,
    );
    let error =
        admit(fixture, BundleMode::Exact).expect_err("case-purpose partition leakage must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::InvalidText {
            field: "case-purpose data partition",
            ..
        }
    ));
    log(
        name,
        "pass",
        &format!("only {expected_partition} rows are authorized"),
        &error.to_string(),
    );
}

#[test]
fn calibration_case_refuses_validation_partition() {
    assert_case_partition_refuses(
        "calibration-consuming-validation",
        CasePurpose::Calibration,
        "val-a",
        "calibration",
    );
}

#[test]
fn validation_case_refuses_blind_partition() {
    assert_case_partition_refuses(
        "validation-consuming-blind",
        CasePurpose::ValidationOnly,
        "blind-a",
        "validation",
    );
}

#[test]
fn blind_case_refuses_calibration_partition() {
    assert_case_partition_refuses(
        "blind-consuming-calibration",
        CasePurpose::BlindFalsification,
        "cal-a",
        "blind holdout",
    );
}

#[test]
fn observation_qoi_absent_from_experiment_refuses() {
    let fixture = problem_fixture(
        vec![case_fixture(
            "a",
            CasePurpose::Calibration,
            "tangent",
            "cal-a",
            "experiment-a",
            "split-a",
            ordinary_data("a"),
        )],
        DataReusePolicy::Disjoint,
    );
    let error = admit(fixture, BundleMode::Exact)
        .expect_err("context QoI without experiment QoI must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::UnknownReference {
            field: "experiment observation QoI",
            ..
        }
    ));
    log(
        "experiment-qoi-closure",
        "pass",
        "observation QoI must occur in the admitted experiment",
        &error.to_string(),
    );
}

#[test]
fn observation_instrument_absent_from_experiment_refuses() {
    let case = case_fixture(
        "a",
        CasePurpose::Calibration,
        "stress",
        "cal-a",
        "experiment-a",
        "split-a",
        ordinary_data("a"),
    )
    .with_observation_instrument("instrument-not-in-experiment");
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("observation instrument outside the experiment roster must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::UnknownReference {
            field: "experiment observation instrument",
            ..
        }
    ));
    log(
        "experiment-instrument-closure",
        "pass",
        "observation instrument must occur in the admitted experiment roster",
        &error.to_string(),
    );
}

#[test]
fn observation_sensor_not_bound_to_instrument_calibration_refuses() {
    let case = case_fixture(
        "a",
        CasePurpose::Calibration,
        "stress",
        "cal-a",
        "experiment-a",
        "split-a",
        ordinary_data("a"),
    )
    .with_observation_sensor("sensor-wrong-calibration");
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("metrology source not bound to instrument certificate must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::SourceMismatch {
            field: "observation sensor/instrument calibration",
        }
    ));
    log(
        "instrument-calibration-source-binding",
        "pass",
        "observation metrology source must equal the admitted instrument certificate",
        &error.to_string(),
    );
}

#[test]
fn observation_clock_absent_from_experiment_refuses() {
    let case = case_fixture(
        "a",
        CasePurpose::Calibration,
        "stress",
        "cal-a",
        "experiment-a",
        "split-a",
        ordinary_data("a"),
    )
    .with_observation_clock("clock-not-in-experiment");
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("observation clock outside the experiment topology must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::UnknownReference {
            field: "experiment observation clock",
            ..
        }
    ));
    log(
        "experiment-clock-closure",
        "pass",
        "observation and protocol clock must occur in the admitted experiment topology",
        &error.to_string(),
    );
}

#[test]
fn globally_valid_alternative_qoi_cannot_be_cross_wired_to_primary_row() {
    assert_manifest_cross_wire_refuses(
        "manifest-row-qoi-cross-wire",
        case_fixture(
            "a",
            CasePurpose::Calibration,
            "tangent",
            "cal-cross-primary",
            "experiment-a",
            "split-a",
            cross_wire_data(false),
        ),
    );
}

#[test]
fn globally_valid_alternative_instrument_cannot_be_cross_wired_to_primary_row() {
    assert_manifest_cross_wire_refuses(
        "manifest-row-instrument-cross-wire",
        case_fixture(
            "a",
            CasePurpose::Calibration,
            "stress",
            "cal-cross-primary",
            "experiment-a",
            "split-a",
            cross_wire_data(false),
        )
        .with_observation_instrument("instrument-cross-wire-z")
        .with_observation_sensor("sensor-cross-wire-z"),
    );
}

#[test]
fn globally_valid_alternative_acquisition_channel_cannot_be_cross_wired_to_primary_row() {
    assert_manifest_cross_wire_refuses(
        "manifest-row-acquisition-channel-cross-wire",
        case_fixture(
            "a",
            CasePurpose::Calibration,
            "stress",
            "cal-cross-primary",
            "experiment-a",
            "split-a",
            cross_wire_data(false),
        )
        .with_observation_acquisition_channel("acquisition-channel-cross-alternative"),
    );
}

#[test]
fn globally_valid_alternative_clock_cannot_be_cross_wired_to_primary_row() {
    assert_manifest_cross_wire_refuses(
        "manifest-row-clock-cross-wire",
        case_fixture(
            "a",
            CasePurpose::Calibration,
            "stress",
            "cal-cross-primary",
            "experiment-a",
            "split-a",
            cross_wire_data(false),
        )
        .with_observation_clock("clock-cross-wire-z"),
    );
}

#[test]
fn swapping_globally_valid_typed_row_sources_changes_exact_artifact_identity() {
    let (baseline, _) = cross_wire_data(false);
    let (swapped, _) = cross_wire_data(true);
    let primary = observation("cal-cross-primary");
    let alternative = observation("cal-cross-alternative");
    assert_eq!(
        baseline
            .manifest()
            .source_ref_of(&primary)
            .expect("baseline primary source"),
        swapped
            .manifest()
            .source_ref_of(&alternative)
            .expect("swapped alternative source"),
    );
    assert_eq!(
        baseline
            .manifest()
            .source_ref_of(&alternative)
            .expect("baseline alternative source"),
        swapped
            .manifest()
            .source_ref_of(&primary)
            .expect("swapped primary source"),
    );
    assert_ne!(
        baseline.manifest().canonical_hash(),
        swapped.manifest().canonical_hash()
    );
    assert_ne!(
        baseline.content_hash().expect("baseline experiment hash"),
        swapped.content_hash().expect("swapped experiment hash"),
    );
    log(
        "manifest-row-source-cross-wire-identity",
        "pass",
        "dataset, locator contract, locator, and extraction receipt remain row-bound",
        "swapping two otherwise valid typed sources changes manifest and artifact identity",
    );
}

#[test]
fn one_raw_row_reused_across_channels_refuses_independent_noise() {
    let case = case_fixture(
        "a",
        CasePurpose::Calibration,
        "stress",
        "cal-a",
        "experiment-a",
        "split-a",
        ordinary_data("a"),
    )
    .with_duplicate_row_channel();
    let error = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect_err("one raw row reused across channels needs explicit dependence");
    assert!(matches!(
        &error,
        IdentifiabilityError::InvalidText {
            field: "observation-sharing coverage",
            ..
        }
    ));
    log(
        "within-case-row-reuse-needs-dependence",
        "pass",
        "one immutable row cannot feed two channels without an exact sharing likelihood",
        &error.to_string(),
    );
}

#[test]
fn one_raw_row_reused_across_channels_admits_with_exact_joint_likelihood() {
    let case = case_fixture(
        "a",
        CasePurpose::Calibration,
        "stress",
        "cal-a",
        "experiment-a",
        "split-a",
        ordinary_data("a"),
    )
    .with_declared_duplicate_row_sharing();
    let admitted = admit(
        problem_fixture(vec![case], DataReusePolicy::Disjoint),
        BundleMode::Exact,
    )
    .expect("exact row-sharing declaration and global joint likelihood admit");
    assert_eq!(
        admitted.document().cases()[&case_id("a")]
            .observation_sharing()
            .len(),
        1,
    );
    log(
        "within-case-row-reuse-explicit-likelihood",
        "pass",
        "repeated raw rows remain representable only through one exact source-bound factor",
        "admitted",
    );
}

#[test]
fn repeated_row_likelihood_must_match_the_global_joint_model() {
    let case = case_fixture(
        "a",
        CasePurpose::Calibration,
        "stress",
        "cal-a",
        "experiment-a",
        "split-a",
        ordinary_data("a"),
    )
    .with_declared_duplicate_row_sharing();
    let error = match try_problem_fixture_with_global_likelihood(
        vec![case],
        DataReusePolicy::Disjoint,
        "different-global-likelihood",
    ) {
        Ok(_) => panic!("a local repeated-row factor cannot differ from the global joint model"),
        Err(error) => error,
    };
    assert!(matches!(&error, IdentifiabilityError::Covariance { .. }));
    log(
        "within-case-row-reuse-global-likelihood-mismatch",
        "pass",
        "every repeated-row factor names the exact global joint-likelihood source",
        &error.to_string(),
    );
}

#[test]
fn distinct_two_case_campaign_admits_under_disjoint_policy() {
    let data_a = ordinary_data("a");
    let data_b = ordinary_data("b");
    assert_ne!(
        data_a.0.content_hash().expect("experiment a hash"),
        data_b.0.content_hash().expect("experiment b hash"),
        "disjoint baseline must use distinct experiment content",
    );
    let fixture = problem_fixture(
        vec![
            case_fixture(
                "a",
                CasePurpose::Calibration,
                "stress",
                "cal-a",
                "experiment-a",
                "split-a",
                data_a,
            ),
            case_fixture(
                "b",
                CasePurpose::Complementary {
                    reason: "independent loading path complements case a".to_string(),
                },
                "stress",
                "cal-b",
                "experiment-b",
                "split-b",
                data_b,
            ),
        ],
        DataReusePolicy::Disjoint,
    );
    let admitted = admit(fixture, BundleMode::Exact).expect("disjoint two-case campaign admits");
    assert_eq!(admitted.data().len(), 2);
    log(
        "disjoint-two-case-baseline",
        "pass",
        "two distinct source-resolved experiments admit under Disjoint",
        &format!("admitted_cases={}", admitted.data().len()),
    );
}

#[test]
fn exact_source_bytes_reuse_refuses_under_disjoint_policy() {
    let data_a = physical_data(
        "a",
        &["stress"],
        &[
            ("cal-a", "source-cal-a"),
            ("val-a", "source-val-a"),
            ("blind-a", "source-blind-a"),
        ],
        &["cal-a"],
        &["val-a"],
        &["blind-a"],
        "shared-raw-source-bytes",
    );
    let data_b = physical_data(
        "b",
        &["stress"],
        &[
            ("cal-b", "source-cal-b"),
            ("val-b", "source-val-b"),
            ("blind-b", "source-blind-b"),
        ],
        &["cal-b"],
        &["val-b"],
        &["blind-b"],
        "shared-raw-source-bytes",
    );
    assert_eq!(
        data_a.0.authenticity().source_bytes_hash(),
        data_b.0.authenticity().source_bytes_hash(),
    );
    assert_ne!(
        data_a.0.manifest().canonical_hash(),
        data_b.0.manifest().canonical_hash(),
    );
    let error = admit(
        problem_fixture(
            vec![
                case_fixture(
                    "a",
                    CasePurpose::Calibration,
                    "stress",
                    "cal-a",
                    "experiment-a",
                    "split-a",
                    data_a,
                ),
                case_fixture(
                    "b",
                    CasePurpose::Calibration,
                    "stress",
                    "cal-b",
                    "experiment-b",
                    "split-b",
                    data_b,
                ),
            ],
            DataReusePolicy::Disjoint,
        ),
        BundleMode::Exact,
    )
    .expect_err("identical raw-source bytes under Disjoint must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::InvalidText {
            field: "data reuse policy",
            ..
        }
    ));
    log(
        "disjoint-source-bytes-alias",
        "pass",
        "distinct manifests cannot hide one exact raw-source byte stream",
        &error.to_string(),
    );
}

#[test]
fn exact_manifest_reuse_refuses_under_disjoint_policy() {
    let rows = [
        ("cal-shared", "source-cal-shared"),
        ("val-shared", "source-val-shared"),
        ("blind-shared", "source-blind-shared"),
    ];
    let data_a = physical_data_with_metrology(
        "a",
        "shared-manifest",
        &["stress"],
        &rows,
        &["cal-shared"],
        &["val-shared"],
        &["blind-shared"],
        "source-bytes-shared-manifest",
    );
    let data_b = physical_data_with_metrology(
        "b",
        "shared-manifest",
        &["stress"],
        &rows,
        &["cal-shared"],
        &["val-shared"],
        &["blind-shared"],
        "source-bytes-shared-manifest",
    );
    assert_eq!(
        data_a.0.manifest().canonical_hash(),
        data_b.0.manifest().canonical_hash(),
    );
    assert_eq!(
        data_a.0.authenticity().source_bytes_hash(),
        data_b.0.authenticity().source_bytes_hash(),
        "a schema-v3 manifest binds its exact dataset byte stream",
    );
    let error = admit(
        problem_fixture(
            vec![
                case_fixture(
                    "a",
                    CasePurpose::Calibration,
                    "stress",
                    "cal-shared",
                    "experiment-a",
                    "split-a",
                    data_a,
                ),
                case_fixture(
                    "b",
                    CasePurpose::Calibration,
                    "stress",
                    "cal-shared",
                    "experiment-b",
                    "split-b",
                    data_b,
                ),
            ],
            DataReusePolicy::Disjoint,
        ),
        BundleMode::Exact,
    )
    .expect_err("identical manifests under Disjoint must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::InvalidText {
            field: "data reuse policy",
            ..
        }
    ));
    log(
        "disjoint-manifest-alias",
        "pass",
        "distinct authenticity wrappers cannot hide one exact observation manifest",
        &error.to_string(),
    );
}

#[test]
fn equal_bare_locator_hash_in_distinct_dataset_scopes_admits_under_disjoint_policy() {
    let data_a = physical_data(
        "a",
        &["stress"],
        &[
            ("cal-a", "shared-immutable-row-source"),
            ("val-a", "source-val-a"),
            ("blind-a", "source-blind-a"),
        ],
        &["cal-a"],
        &["val-a"],
        &["blind-a"],
        "source-bytes-a",
    );
    let data_b = physical_data(
        "b",
        &["stress"],
        &[
            ("cal-b", "shared-immutable-row-source"),
            ("val-b", "source-val-b"),
            ("blind-b", "source-blind-b"),
        ],
        &["cal-b"],
        &["val-b"],
        &["blind-b"],
        "source-bytes-b",
    );
    assert_ne!(
        data_a.0.authenticity().source_bytes_hash(),
        data_b.0.authenticity().source_bytes_hash(),
        "fixture must place the equal bare locator hash in different dataset scopes",
    );
    assert_ne!(
        data_a.0.manifest().canonical_hash(),
        data_b.0.manifest().canonical_hash(),
        "fixture must isolate row-source aliasing from manifest identity",
    );
    let source_a = data_a
        .0
        .manifest()
        .source_ref_of(&observation("cal-a"))
        .expect("case-a typed row source");
    let source_b = data_b
        .0
        .manifest()
        .source_ref_of(&observation("cal-b"))
        .expect("case-b typed row source");
    assert_eq!(
        source_a.locator_hash(),
        source_b.locator_hash(),
        "fixture deliberately collides the bare locator hash",
    );
    assert_ne!(
        source_a.locator_identity(),
        source_b.locator_identity(),
        "different dataset-byte identities must keep equal bare locator hashes in distinct receipt-independent locator scopes",
    );
    assert_ne!(
        source_a, source_b,
        "the complete typed row sources retain their distinct dataset and provenance fields",
    );
    let fixture = problem_fixture(
        vec![
            case_fixture(
                "a",
                CasePurpose::Calibration,
                "stress",
                "cal-a",
                "experiment-a",
                "split-a",
                data_a,
            ),
            case_fixture(
                "b",
                CasePurpose::Calibration,
                "stress",
                "cal-b",
                "experiment-b",
                "split-b",
                data_b,
            ),
        ],
        DataReusePolicy::Disjoint,
    );
    let admitted = admit(fixture, BundleMode::Exact)
        .expect("a bare locator hash alone is not cross-case raw-row identity");
    assert_eq!(admitted.data().len(), 2);
    log(
        "disjoint-locator-hash-coincidence",
        "pass",
        "dataset-scoped locator identity defines reuse independently of extraction-receipt relabeling; a bare locator hash does not",
        &format!("admitted_cases={}", admitted.data().len()),
    );
}

#[test]
fn declared_sharing_group_admits_one_joint_raw_campaign() {
    let shared = physical_data(
        "shared",
        &["stress"],
        &[
            ("cal-a", "source-cal-a"),
            ("cal-b", "source-cal-b"),
            ("val-shared", "source-val-shared"),
            ("blind-shared", "source-blind-shared"),
        ],
        &["cal-a", "cal-b"],
        &["val-shared"],
        &["blind-shared"],
        "source-bytes-shared",
    );
    let fixture = problem_fixture(
        vec![
            case_fixture(
                "a",
                CasePurpose::Calibration,
                "stress",
                "cal-a",
                "experiment-a",
                "split-a",
                shared.clone(),
            ),
            case_fixture(
                "b",
                CasePurpose::SymmetryBreaking,
                "stress",
                "cal-b",
                "experiment-b",
                "split-b",
                shared,
            ),
        ],
        DataReusePolicy::Shared {
            groups: vec![
                DataSharingGroup::try_new(
                    BTreeSet::from([case_id("a"), case_id("b")]),
                    source_key("joint-likelihood"),
                    "both cases intentionally consume complementary channels from one raw campaign",
                )
                .expect("sharing group fixture"),
            ],
        },
    );
    let admitted = admit(fixture, BundleMode::Exact).expect("declared sharing admits");
    assert_eq!(admitted.data().len(), 2);
    log(
        "declared-joint-sharing",
        "pass",
        "one shared experiment plus an explicit joint likelihood",
        &format!("admitted_cases={}", admitted.data().len()),
    );
}

#[test]
fn cross_case_sharing_likelihood_must_match_the_global_joint_model() {
    let shared = physical_data(
        "shared-mismatched-likelihood",
        &["stress"],
        &[
            ("cal-a", "source-cal-a"),
            ("cal-b", "source-cal-b"),
            ("val-shared", "source-val-shared"),
            ("blind-shared", "source-blind-shared"),
        ],
        &["cal-a", "cal-b"],
        &["val-shared"],
        &["blind-shared"],
        "source-bytes-shared-mismatched-likelihood",
    );
    let cases = vec![
        case_fixture(
            "a",
            CasePurpose::Calibration,
            "stress",
            "cal-a",
            "experiment-a",
            "split-a",
            shared.clone(),
        ),
        case_fixture(
            "b",
            CasePurpose::SymmetryBreaking,
            "stress",
            "cal-b",
            "experiment-b",
            "split-b",
            shared,
        ),
    ];
    let reuse = DataReusePolicy::Shared {
        groups: vec![
            DataSharingGroup::try_new(
                BTreeSet::from([case_id("a"), case_id("b")]),
                source_key("case-sharing-likelihood"),
                "the exact shared campaign deliberately names a non-global factor",
            )
            .expect("cross-case sharing group"),
        ],
    };
    let error = match try_problem_fixture_with_global_likelihood(cases, reuse, "joint-likelihood") {
        Ok(_) => panic!("a cross-case sharing factor cannot differ from the global joint model"),
        Err(error) => error,
    };
    assert!(matches!(&error, IdentifiabilityError::Covariance { .. }));
    log(
        "cross-case-reuse-global-likelihood-mismatch",
        "pass",
        "every cross-case sharing factor names the exact global joint-likelihood source",
        &error.to_string(),
    );
}

#[test]
fn declared_sharing_group_without_actual_overlap_refuses() {
    let fixture = problem_fixture(
        vec![
            case_fixture(
                "a",
                CasePurpose::Calibration,
                "stress",
                "cal-a",
                "experiment-a",
                "split-a",
                ordinary_data("a"),
            ),
            case_fixture(
                "b",
                CasePurpose::Complementary {
                    reason: "independent campaign used as a false sharing declaration".to_string(),
                },
                "stress",
                "cal-b",
                "experiment-b",
                "split-b",
                ordinary_data("b"),
            ),
        ],
        DataReusePolicy::Shared {
            groups: vec![
                DataSharingGroup::try_new(
                    BTreeSet::from([case_id("a"), case_id("b")]),
                    source_key("joint-likelihood"),
                    "this deliberately false declaration must not create sharing authority",
                )
                .expect("sharing group fixture"),
            ],
        },
    );
    let error = admit(fixture, BundleMode::Exact)
        .expect_err("a sharing declaration without admitted overlap must refuse");
    assert!(matches!(
        &error,
        IdentifiabilityError::InvalidText {
            field: "data sharing group",
            ..
        }
    ));
    log(
        "sharing-declaration-needs-overlap",
        "pass",
        "declaration alone cannot manufacture raw-data sharing authority",
        &error.to_string(),
    );
}

#[test]
fn missing_retrospective_case_bundle_refuses() {
    let missing_fixture = problem_fixture(
        vec![case_fixture(
            "a",
            CasePurpose::Calibration,
            "stress",
            "cal-a",
            "experiment-a",
            "split-a",
            ordinary_data("a"),
        )],
        DataReusePolicy::Disjoint,
    );
    let missing = admit(missing_fixture, BundleMode::Missing)
        .expect_err("missing retrospective case bundle must refuse");
    assert!(matches!(
        &missing,
        IdentifiabilityError::UnknownReference {
            field: "retrospective case source bundle",
            ..
        }
    ));
    log(
        "missing-case-bundle",
        "pass",
        "every retrospective case has one concrete bundle",
        &missing.to_string(),
    );
}

#[test]
fn unknown_retrospective_case_bundle_refuses() {
    let extra_fixture = problem_fixture(
        vec![case_fixture(
            "a",
            CasePurpose::Calibration,
            "stress",
            "cal-a",
            "experiment-a",
            "split-a",
            ordinary_data("a"),
        )],
        DataReusePolicy::Disjoint,
    );
    let extra = admit(extra_fixture, BundleMode::Extra)
        .expect_err("extra retrospective case bundle must refuse");
    assert!(matches!(
        &extra,
        IdentifiabilityError::Cardinality {
            field: "case source bundles",
            ..
        }
    ));
    log(
        "extra-case-bundle",
        "pass",
        "no unknown case bundle is accepted",
        &extra.to_string(),
    );
}
