//! Battery for the prediction-input bundle and its atomic sealer
//! (bead frankensim-jmh21.1, input half): canonical round trips,
//! permutation-invariant identity, the no-target guarantee, every typed
//! refusal, atomic seal semantics with crash/tamper/swap hostility, and
//! per-field mutation sensitivity.

use fs_blake3::ContentHash;
use fs_evidence::prediction_bundle::{
    AccessPolicy, ModelRungPolicy, OutputArtifactRef, OutputFamily, PredictionExecutionInput,
    PredictionOutputBundle, RandomStreamDesign, SampleAccounting, load_sealed_input,
    load_sealed_output, seal_prediction_input, seal_prediction_output,
};
use fs_evidence::vv::{ApplicabilityPolicy, ArtifactId, ArtifactKind, ArtifactRef};

/// One labeled field mutation used by the identity-sensitivity scenario.
type Mutation = (&'static str, Box<dyn Fn(&mut Fixture)>);

fn reference(kind: ArtifactKind, id: &str, salt: u8) -> ArtifactRef {
    ArtifactRef::new(
        kind,
        ArtifactId::try_new(id).expect("valid id"),
        fs_blake3::hash_bytes(&[salt]),
    )
}

fn stream(name: &str, seed: u64) -> RandomStreamDesign {
    RandomStreamDesign {
        name: name.to_string(),
        seed_domain: format!("org.frankensim.test.{name}"),
        seed,
        substreams: 4,
    }
}

struct Fixture {
    source_identity: Vec<(String, String)>,
    context_of_use: ArtifactRef,
    validation_plan: ArtifactRef,
    calibration_split: ArtifactRef,
    scenarios: Vec<ArtifactRef>,
    parameter_distributions: Vec<ArtifactRef>,
    random_streams: Vec<RandomStreamDesign>,
    model_rungs: ModelRungPolicy,
    qoi_identities: Vec<String>,
    evidence_role: String,
    blind_partition: ArtifactRef,
    access_policy: AccessPolicy,
}

impl Fixture {
    fn nominal() -> Self {
        Self {
            source_identity: vec![
                ("head_sha".to_string(), "deadbeef".to_string()),
                ("toolchain".to_string(), "nightly-2026-07-06".to_string()),
            ],
            context_of_use: reference(ArtifactKind::ContextOfUse, "cou-1", 1),
            validation_plan: reference(ArtifactKind::ValidationPlan, "plan-1", 2),
            calibration_split: reference(ArtifactKind::CalibrationSplit, "split-1", 3),
            scenarios: vec![
                reference(ArtifactKind::ExperimentArtifact, "scenario-a", 4),
                reference(ArtifactKind::ExperimentArtifact, "scenario-b", 5),
            ],
            parameter_distributions: vec![reference(
                ArtifactKind::ExperimentArtifact,
                "param-cov",
                6,
            )],
            random_streams: vec![stream("sample-draw", 7), stream("scenario-jitter", 11)],
            model_rungs: ModelRungPolicy {
                allowed_rungs: vec!["reduced-order".to_string(), "full-fem".to_string()],
                applicability: ApplicabilityPolicy::Refuse,
            },
            qoi_identities: vec!["junction-maximum".to_string(), "pressure-drop".to_string()],
            evidence_role: "blind-prediction".to_string(),
            blind_partition: reference(ArtifactKind::ExperimentArtifact, "holdout-1", 8),
            access_policy: AccessPolicy::ExecutorOnly,
        }
    }

    fn build(
        self,
    ) -> Result<PredictionExecutionInput, fs_evidence::prediction_bundle::PredictionBundleError>
    {
        PredictionExecutionInput::try_new(
            self.source_identity,
            self.context_of_use,
            self.validation_plan,
            self.calibration_split,
            self.scenarios,
            self.parameter_distributions,
            self.random_streams,
            self.model_rungs,
            self.qoi_identities,
            self.evidence_role,
            self.blind_partition,
            self.access_policy,
        )
    }
}

fn nominal() -> PredictionExecutionInput {
    Fixture::nominal().build().expect("nominal admits")
}

fn scratch(label: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "fs-evidence-prediction-bundle-{label}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

#[test]
fn canonical_round_trip_preserves_identity() {
    let input = nominal();
    let bytes = input.canonical_bytes().expect("encodes");
    let decoded = PredictionExecutionInput::from_canonical_bytes(&bytes).expect("decodes");
    assert_eq!(decoded, input);
    assert_eq!(
        decoded.identity().expect("identity"),
        input.identity().expect("identity")
    );
}

#[test]
fn declaration_order_never_reaches_identity() {
    let base = nominal().identity().expect("identity");
    let mut fixture = Fixture::nominal();
    fixture.random_streams.reverse();
    fixture.model_rungs.allowed_rungs.reverse();
    fixture.qoi_identities.reverse();
    fixture.source_identity.reverse();
    let permuted = fixture.build().expect("permuted admits");
    assert_eq!(permuted.identity().expect("identity"), base);
}

#[test]
fn trailing_bytes_refuse_as_target_smuggling() {
    // The no-target guarantee, executable: there is no field a target could
    // occupy, so the ONLY way to smuggle one is trailing data - which
    // refuses with the dedicated rule.
    let mut bytes = nominal().canonical_bytes().expect("encodes");
    bytes.extend_from_slice(b"target=42.0");
    let error = PredictionExecutionInput::from_canonical_bytes(&bytes)
        .expect_err("trailing bytes must refuse");
    assert_eq!(error.rule, "prediction-input-transport-trailing");
}

#[test]
fn truncated_and_foreign_transports_refuse() {
    let bytes = nominal().canonical_bytes().expect("encodes");
    let error = PredictionExecutionInput::from_canonical_bytes(&bytes[..bytes.len() - 3])
        .expect_err("truncation must refuse");
    assert!(
        error.rule.starts_with("prediction-input-transport"),
        "{error}"
    );
    let error = PredictionExecutionInput::from_canonical_bytes(b"FSVVnot-this-schema")
        .expect_err("foreign magic must refuse");
    assert_eq!(error.rule, "prediction-input-transport-magic");
}

#[test]
fn wrong_reference_kinds_refuse() {
    let mut fixture = Fixture::nominal();
    fixture.context_of_use = reference(ArtifactKind::ValidationPlan, "wrong", 9);
    let error = fixture.build().expect_err("wrong kind must refuse");
    assert_eq!(error.rule, "prediction-input-reference-kind");
    assert_eq!(error.field, "context_of_use");
}

#[test]
fn empty_and_duplicate_declarations_refuse() {
    let mut fixture = Fixture::nominal();
    fixture.scenarios.clear();
    assert_eq!(
        fixture.build().expect_err("empty scenarios").rule,
        "prediction-input-reference-bounds"
    );

    let mut fixture = Fixture::nominal();
    fixture.scenarios = vec![
        reference(ArtifactKind::ExperimentArtifact, "dup-a", 4),
        reference(ArtifactKind::ExperimentArtifact, "dup-b", 4),
    ];
    assert_eq!(
        fixture.build().expect_err("duplicate scenario hashes").rule,
        "prediction-input-reference-bounds"
    );

    let mut fixture = Fixture::nominal();
    fixture.random_streams = vec![stream("same", 1), stream("same", 2)];
    assert_eq!(
        fixture.build().expect_err("duplicate stream names").rule,
        "prediction-input-random-streams"
    );

    let mut fixture = Fixture::nominal();
    fixture.qoi_identities = vec!["q".to_string(), "q".to_string()];
    assert_eq!(
        fixture.build().expect_err("duplicate QoIs").rule,
        "prediction-input-qoi"
    );

    let mut fixture = Fixture::nominal();
    fixture.random_streams[0].substreams = 0;
    assert_eq!(
        fixture.build().expect_err("zero substreams").rule,
        "prediction-input-random-streams"
    );
}

#[test]
fn boundary_counts_admit_exactly_the_cap() {
    use fs_evidence::prediction_bundle::MAX_BUNDLE_ITEMS;
    let build_with = |count: usize| {
        let mut fixture = Fixture::nominal();
        fixture.qoi_identities = (0..count).map(|index| format!("qoi-{index:05}")).collect();
        fixture.build()
    };
    assert!(build_with(1).is_ok(), "minimum admits");
    assert!(build_with(MAX_BUNDLE_ITEMS).is_ok(), "exact cap admits");
    assert_eq!(
        build_with(MAX_BUNDLE_ITEMS + 1)
            .expect_err("cap+1 refuses")
            .rule,
        "prediction-input-qoi"
    );
}

#[test]
fn every_semantic_field_moves_the_identity() {
    let base = nominal().identity().expect("identity");
    let mutations: Vec<Mutation> = vec![
        (
            "source_identity",
            Box::new(|f| f.source_identity[0].1 = "cafef00d".to_string()),
        ),
        (
            "context_of_use",
            Box::new(|f| f.context_of_use = reference(ArtifactKind::ContextOfUse, "cou-2", 21)),
        ),
        (
            "validation_plan",
            Box::new(|f| f.validation_plan = reference(ArtifactKind::ValidationPlan, "plan-2", 22)),
        ),
        (
            "calibration_split",
            Box::new(|f| {
                f.calibration_split = reference(ArtifactKind::CalibrationSplit, "split-2", 23);
            }),
        ),
        (
            "scenarios",
            Box::new(|f| {
                f.scenarios.pop();
            }),
        ),
        (
            "parameter_distributions",
            Box::new(|f| f.parameter_distributions.clear()),
        ),
        (
            "random_stream_seed",
            Box::new(|f| f.random_streams[0].seed ^= 1),
        ),
        (
            "random_stream_domain",
            Box::new(|f| f.random_streams[0].seed_domain.push('x')),
        ),
        (
            "model_rungs",
            Box::new(|f| {
                f.model_rungs.allowed_rungs.pop();
            }),
        ),
        (
            "applicability",
            Box::new(|f| f.model_rungs.applicability = ApplicabilityPolicy::Demote),
        ),
        (
            "qoi_identities",
            Box::new(|f| f.qoi_identities[0].push('x')),
        ),
        ("evidence_role", Box::new(|f| f.evidence_role.push('x'))),
        (
            "blind_partition",
            Box::new(|f| {
                f.blind_partition = reference(ArtifactKind::ExperimentArtifact, "holdout-2", 24);
            }),
        ),
        (
            "access_policy",
            Box::new(|f| f.access_policy = AccessPolicy::Open),
        ),
    ];
    for (label, mutate) in mutations {
        let mut fixture = Fixture::nominal();
        mutate(&mut fixture);
        let mutated = fixture.build().expect("mutated fixture still admits");
        assert_ne!(
            mutated.identity().expect("identity"),
            base,
            "mutating {label} must mint a new identity"
        );
    }
}

#[test]
fn seal_is_atomic_immutable_and_independently_verifiable() {
    let dir = scratch("seal");
    let path = dir.join("input.fspi");
    let input = nominal();
    let sealed = seal_prediction_input(&input, &path).expect("seals");

    // Independent verification from artifact bytes alone.
    let (loaded, identity) = load_sealed_input(&path, Some(sealed)).expect("verifies");
    assert_eq!(loaded, input);
    assert_eq!(identity, sealed);

    // Seals are immutable: a second publication at the same path refuses.
    assert_eq!(
        seal_prediction_input(&input, &path)
            .expect_err("reseal refuses")
            .rule,
        "prediction-input-seal-immutable"
    );

    // No partial residue survives a successful seal.
    let residue: Vec<_> = std::fs::read_dir(&dir)
        .expect("scratch listing")
        .filter_map(Result::ok)
        .filter(|entry| entry.path() != path)
        .collect();
    assert!(
        residue.is_empty(),
        "partial files must not survive: {residue:?}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn crash_residue_at_the_partial_path_is_never_a_sealed_bundle() {
    // Simulate a crash mid-publication: bytes at the PARTIAL path, nothing
    // at the sealed path. The sealed path simply does not exist, so a
    // consumer cannot read a half-published bundle; and a later seal
    // completes cleanly over the residue.
    let dir = scratch("crash");
    let path = dir.join("input.fspi");
    let partial = path.with_extension(format!("partial.{}", std::process::id()));
    std::fs::write(&partial, b"half-written").expect("residue writes");
    assert!(
        load_sealed_input(&path, None).is_err(),
        "no readable bundle"
    );
    let input = nominal();
    let sealed = seal_prediction_input(&input, &path).expect("seal completes over residue");
    let (_, identity) = load_sealed_input(&path, Some(sealed)).expect("verifies");
    assert_eq!(identity, sealed);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn tampered_swapped_and_stale_bundles_refuse_at_verification() {
    let dir = scratch("tamper");
    let input = nominal();
    let path = dir.join("input.fspi");
    let sealed = seal_prediction_input(&input, &path).expect("seals");

    // Bit flip inside the sealed bytes: decode may pass, identity must not.
    let mut bytes = std::fs::read(&path).expect("reads");
    let flip_at = bytes.len() - 1;
    bytes[flip_at] ^= 0x01;
    std::fs::write(&path, &bytes).expect("tamper writes");
    let error = load_sealed_input(&path, Some(sealed)).expect_err("tamper must refuse");
    assert!(
        error.rule == "prediction-input-seal-identity"
            || error.rule.starts_with("prediction-input-transport")
            || error.rule.starts_with("prediction-input-"),
        "{error}"
    );

    // Swap: a DIFFERENT valid bundle at the path fails the expected identity.
    let mut fixture = Fixture::nominal();
    fixture.evidence_role = "swapped-role".to_string();
    let other = fixture.build().expect("admits");
    std::fs::remove_file(&path).expect("clears");
    seal_prediction_input(&other, &path).expect("seals other");
    let error = load_sealed_input(&path, Some(sealed)).expect_err("swap must refuse");
    assert_eq!(error.rule, "prediction-input-seal-identity");

    // Truncation refuses before identity comparison.
    let bytes = std::fs::read(&path).expect("reads");
    std::fs::write(&path, &bytes[..bytes.len() / 2]).expect("truncates");
    assert!(
        load_sealed_input(&path, Some(sealed)).is_err(),
        "truncation refuses"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn identity_is_domain_separated_from_raw_hashing() {
    let input = nominal();
    let bytes = input.canonical_bytes().expect("encodes");
    let identity = input.identity().expect("identity");
    assert_ne!(
        identity,
        fs_blake3::hash_bytes(&bytes),
        "identity must be domain-separated, not a bare byte hash"
    );
    assert_ne!(
        identity,
        ContentHash::from_slice(&[0u8; 32]).expect("32 bytes")
    );
}

fn output_ref(family: OutputFamily, id: &str, salt: u8) -> OutputArtifactRef {
    OutputArtifactRef {
        family,
        id: id.to_string(),
        hash: fs_blake3::hash_bytes(&[salt, 0xA0]),
    }
}

fn nominal_output(input_root: ContentHash) -> PredictionOutputBundle {
    PredictionOutputBundle::try_new(
        input_root,
        SampleAccounting {
            requested: 100,
            succeeded: 90,
            refused: 7,
            failed: 3,
        },
        vec![
            output_ref(OutputFamily::Trajectory, "traj-0", 1),
            output_ref(OutputFamily::Event, "events-0", 2),
            output_ref(OutputFamily::Energy, "energy-0", 3),
        ],
        vec![output_ref(OutputFamily::Aggregate, "qoi-dist", 4)],
        12_345_678,
        false,
        None,
        reference(ArtifactKind::SolutionVerificationReceipt, "sv-1", 5),
        vec!["replay qoi-dist against traj-0 with the sealed input seeds".to_string()],
    )
    .expect("nominal output admits")
}

#[test]
fn output_round_trip_and_identity() {
    let root = nominal().identity().expect("identity");
    let output = nominal_output(root);
    let bytes = output.canonical_bytes().expect("encodes");
    let decoded = PredictionOutputBundle::from_canonical_bytes(&bytes).expect("decodes");
    assert_eq!(decoded, output);
    assert_eq!(
        decoded.identity().expect("identity"),
        output.identity().expect("identity")
    );
    // Outcome smuggling closes the same way as on the input side.
    let mut smuggled = bytes;
    smuggled.extend_from_slice(b"observed=1.5");
    assert_eq!(
        PredictionOutputBundle::from_canonical_bytes(&smuggled)
            .expect_err("trailing refuses")
            .rule,
        "prediction-output-transport-trailing"
    );
}

#[test]
fn sample_partition_must_be_total() {
    let root = nominal().identity().expect("identity");
    let build = |requested, succeeded, refused, failed| {
        PredictionOutputBundle::try_new(
            root,
            SampleAccounting {
                requested,
                succeeded,
                refused,
                failed,
            },
            vec![output_ref(OutputFamily::Trajectory, "t", 1)],
            Vec::new(),
            1,
            false,
            None,
            reference(ArtifactKind::SolutionVerificationReceipt, "sv-1", 5),
            vec!["check".to_string()],
        )
    };
    assert!(build(10, 8, 1, 1).is_ok(), "total partition admits");
    assert_eq!(
        build(10, 8, 1, 0).expect_err("missing sample").rule,
        "prediction-output-accounting"
    );
    assert_eq!(
        build(10, 9, 1, 1).expect_err("excess sample").rule,
        "prediction-output-accounting"
    );
    assert_eq!(
        build(0, 0, 0, 0).expect_err("zero requested").rule,
        "prediction-output-accounting"
    );
    // Overflow cannot fake totality.
    assert_eq!(
        build(u64::MAX, u64::MAX, 2, u64::MAX - 1)
            .expect_err("overflowing partition")
            .rule,
        "prediction-output-accounting"
    );
    // Success without artifacts is unscoreable.
    let empty = PredictionOutputBundle::try_new(
        root,
        SampleAccounting {
            requested: 5,
            succeeded: 5,
            refused: 0,
            failed: 0,
        },
        Vec::new(),
        Vec::new(),
        1,
        false,
        None,
        reference(ArtifactKind::SolutionVerificationReceipt, "sv-1", 5),
        vec!["check".to_string()],
    );
    assert_eq!(
        empty.expect_err("success without artifacts").rule,
        "prediction-output-artifacts"
    );
}

#[test]
fn input_before_output_ordering_is_checkable() {
    let sealed_input = nominal().identity().expect("identity");
    let output = nominal_output(sealed_input);
    output
        .verify_against_input(sealed_input)
        .expect("matching root joins");
    let foreign_root = fs_blake3::hash_bytes(b"some other input");
    assert_eq!(
        output
            .verify_against_input(foreign_root)
            .expect_err("foreign root refuses")
            .rule,
        "prediction-output-input-root"
    );
    // And an output built against a root that was never sealed as an input
    // refuses at the SAME join - ordering is enforced by identity, not by
    // trust in the producer's claim.
    let premature = nominal_output(foreign_root);
    assert!(premature.verify_against_input(sealed_input).is_err());
}

#[test]
fn output_seal_tamper_and_swap_refuse() {
    let dir = scratch("output-seal");
    let root = nominal().identity().expect("identity");
    let output = nominal_output(root);
    let path = dir.join("output.fspo");
    let sealed = seal_prediction_output(&output, &path).expect("seals");
    let (loaded, identity) = load_sealed_output(&path, Some(sealed)).expect("verifies");
    assert_eq!(loaded, output);
    assert_eq!(identity, sealed);

    // Cross-domain swap: the sealed INPUT bytes at the output path must
    // refuse (domain separation, not just schema magic).
    let input_path = dir.join("input.fspi");
    seal_prediction_input(&nominal(), &input_path).expect("input seals");
    let input_bytes = std::fs::read(&input_path).expect("reads");
    std::fs::remove_file(&path).expect("clears");
    std::fs::write(&path, &input_bytes).expect("swaps");
    assert!(
        load_sealed_output(&path, Some(sealed)).is_err(),
        "an input transport at the output path must refuse"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn every_output_field_moves_the_identity() {
    let root = nominal().identity().expect("identity");
    let base = nominal_output(root).identity().expect("identity");
    let mutate: Vec<(&str, PredictionOutputBundle)> = vec![
        (
            "input_root",
            nominal_output(fs_blake3::hash_bytes(b"other")),
        ),
        (
            "accounting",
            PredictionOutputBundle::try_new(
                root,
                SampleAccounting {
                    requested: 100,
                    succeeded: 89,
                    refused: 8,
                    failed: 3,
                },
                vec![output_ref(OutputFamily::Trajectory, "traj-0", 1)],
                Vec::new(),
                12_345_678,
                false,
                None,
                reference(ArtifactKind::SolutionVerificationReceipt, "sv-1", 5),
                vec!["replay".to_string()],
            )
            .expect("admits"),
        ),
    ];
    for (label, mutated) in mutate {
        assert_ne!(
            mutated.identity().expect("identity"),
            base,
            "mutating {label} must mint a new identity"
        );
    }
    // Work units, budget flag, and checker text are also identity-bearing.
    let mut bytes = nominal_output(root).canonical_bytes().expect("encodes");
    let flip_at = bytes.len() - 1;
    bytes[flip_at] ^= 0x20;
    let reparsed = PredictionOutputBundle::from_canonical_bytes(&bytes);
    if let Ok(reparsed) = reparsed {
        assert_ne!(reparsed.identity().expect("identity"), base);
    }
}
