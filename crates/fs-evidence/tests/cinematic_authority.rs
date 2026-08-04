//! G0/G3 authority and negative-twin tests for Euler cinematic artifacts.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::cinematic::{
    CINEMATIC_AUTHORITY_SCHEMA_VERSION, CinematicArtifactKind, CinematicAudioMaster,
    CinematicAuthorityClass, CinematicAuthorityError, CinematicAuthorityInput,
    CinematicAuthorityRecord, CinematicClock, CinematicClockDomain, CinematicDeliverableContract,
    CinematicDeliverableError, CinematicDisplayPreview, CinematicImageMaster,
    CinematicMuxedDerivative, CinematicNoClaim, CinematicTerm, CinematicTransformDisposition,
    CinematicUnitContract, DeclaredAcousticCalibrationReceipt, MAX_CINEMATIC_LABEL_BYTES,
    MAX_CINEMATIC_NO_CLAIMS, SoundAuthority, required_no_claims,
};

fn identity(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.test.cinematic-authority.v1",
        label.as_bytes(),
    )
}

fn calibration() -> DeclaredAcousticCalibrationReceipt {
    DeclaredAcousticCalibrationReceipt::try_new(
        identity("calibration-dataset"),
        identity("calibration-method"),
        identity("calibration-validity-domain"),
        7,
    )
    .expect("valid test calibration declaration")
}

fn input_for(
    label: &str,
    artifact_kind: CinematicArtifactKind,
    authority_class: CinematicAuthorityClass,
    transform_disposition: CinematicTransformDisposition,
) -> CinematicAuthorityInput {
    let (unit_contract, clock) = match authority_class {
        CinematicAuthorityClass::SimulatedState => (
            CinematicUnitContract::SiMechanicsRadians,
            CinematicClock::try_new(CinematicClockDomain::Simulation, 1_000, 1, 0, 2_000)
                .expect("valid simulation clock"),
        ),
        CinematicAuthorityClass::MonteCarloRender => (
            CinematicUnitContract::SpectralRadianceSi,
            CinematicClock::try_new(CinematicClockDomain::Video, 24, 1, 0, 48)
                .expect("valid video clock"),
        ),
        CinematicAuthorityClass::VisualizationDerivative => (
            CinematicUnitContract::DisplayEncoded,
            CinematicClock::try_new(CinematicClockDomain::Video, 24_000, 1_001, 0, 48)
                .expect("valid rational video clock"),
        ),
        CinematicAuthorityClass::Sound(_) => (
            CinematicUnitContract::DigitalAudioFullScale,
            CinematicClock::try_new(CinematicClockDomain::Audio, 48_000, 1, 0, 96_000)
                .expect("valid audio clock"),
        ),
    };
    CinematicAuthorityInput {
        schema_version: CINEMATIC_AUTHORITY_SCHEMA_VERSION,
        artifact_kind,
        authority_class,
        artifact_identity: identity(&format!("artifact-{label}")),
        source_identity: identity(&format!("source-{label}")),
        transform_identity: identity(&format!("transform-{label}")),
        transform_name: format!("cinematic/{label}"),
        configuration_identity: identity(&format!("configuration-{label}")),
        configuration_version: 3,
        unit_contract,
        clock,
        transform_disposition,
        no_claims: required_no_claims(authority_class),
        acoustic_calibration: match authority_class {
            CinematicAuthorityClass::Sound(SoundAuthority::Calibrated) => Some(calibration()),
            _ => None,
        },
    }
}

fn admitted(input: CinematicAuthorityInput) -> CinematicAuthorityRecord {
    CinematicAuthorityRecord::try_new(input).expect("test record should be admitted")
}

fn assert_error_code(
    result: Result<CinematicAuthorityRecord, CinematicAuthorityError>,
    expected: &'static str,
    context: &str,
) {
    let error = result.expect_err(context);
    assert_eq!(error.code(), expected, "{context}: {error}");
}

fn disposition_offset(bytes: &[u8]) -> usize {
    const TRANSFORM_NAME_LENGTH_OFFSET: usize = 171;
    let name_len = usize::from(u16::from_le_bytes([
        bytes[TRANSFORM_NAME_LENGTH_OFFSET],
        bytes[TRANSFORM_NAME_LENGTH_OFFSET + 1],
    ]));
    TRANSFORM_NAME_LENGTH_OFFSET + 2 + name_len
}

fn calibration_presence_offset(bytes: &[u8]) -> usize {
    let disposition = disposition_offset(bytes);
    let label_len_offset = disposition + 1;
    let label_len = usize::from(u16::from_le_bytes([
        bytes[label_len_offset],
        bytes[label_len_offset + 1],
    ]));
    label_len_offset + 2 + label_len
}

#[test]
fn canonical_deliverable_freezes_4k_av_masters_and_exact_timeline() {
    let contract = CinematicDeliverableContract::euler_disc_v1();
    assert_eq!(contract.width_pixels(), 3_840);
    assert_eq!(contract.height_pixels(), 2_160);
    assert_eq!(contract.frames_per_second_numerator(), 24);
    assert_eq!(contract.frames_per_second_denominator(), 1);
    assert_eq!(contract.minimum_frame_count(), 192);
    assert_eq!(contract.maximum_frame_count(), 288);
    assert_eq!(contract.image_master(), CinematicImageMaster::OpenExrFloat);
    assert_eq!(
        contract.display_preview(),
        CinematicDisplayPreview::DisplayReferred
    );
    assert_eq!(contract.audio_sample_rate_hz(), 48_000);
    assert_eq!(contract.audio_channels(), 2);
    assert_eq!(contract.audio_master(), CinematicAudioMaster::WaveFloat32);
    assert!(contract.sequence_manifest_required());
    assert_eq!(
        contract.muxed_derivative(),
        CinematicMuxedDerivative::OptionalNonAuthoritative
    );

    contract
        .validate_timeline(192, 384_000)
        .expect("8-second inclusive boundary must synchronize");
    contract
        .validate_timeline(288, 576_000)
        .expect("12-second inclusive boundary must synchronize");
    assert_eq!(
        contract.validate_timeline(191, 382_000),
        Err(CinematicDeliverableError::FrameCountOutOfRange {
            got: 191,
            minimum: 192,
            maximum: 288,
        })
    );
    assert_eq!(
        contract.validate_timeline(240, 479_999),
        Err(CinematicDeliverableError::AudioVideoClockMismatch {
            video_frame_count: 240,
            audio_sample_frames: 479_999,
            expected_audio_sample_frames: 480_000,
        })
    );

    let manifest = contract.to_manifest_json();
    for required in [
        "euler-disc-cinematic-v1",
        "\"width_pixels\":3840",
        "\"height_pixels\":2160",
        "openexr-float",
        "display-referred-preview",
        "wave-float32",
        "optional-non-authoritative",
    ] {
        assert!(
            manifest.contains(required),
            "deliverable manifest omitted {required}"
        );
    }
}

#[test]
fn frozen_terminology_keeps_aesthetic_and_scientific_claims_separate() {
    let mut codes = Vec::new();
    for term in CinematicTerm::ALL {
        assert!(!term.code().is_empty());
        assert!(!term.definition().is_empty());
        codes.push(term.code());
    }
    codes.sort_unstable();
    codes.dedup();
    assert_eq!(codes.len(), CinematicTerm::ALL.len());
    assert!(
        CinematicTerm::RenderConvergence
            .definition()
            .contains("nothing about physical-model validity")
    );
    assert!(
        CinematicTerm::VisualApproval
            .definition()
            .contains("cannot promote scientific")
    );
    assert!(
        CinematicTerm::MediaEncoding
            .definition()
            .contains("without increasing their authority")
    );
}

#[test]
fn all_authority_classes_are_distinct_admitted_canonical_records() {
    let cases = [
        (
            "state",
            CinematicArtifactKind::SimulationState,
            CinematicAuthorityClass::SimulatedState,
            CinematicTransformDisposition::ModelState,
        ),
        (
            "raw-render",
            CinematicArtifactKind::RenderEstimate,
            CinematicAuthorityClass::MonteCarloRender,
            CinematicTransformDisposition::MonteCarloEstimator,
        ),
        (
            "visualization",
            CinematicArtifactKind::Visualization,
            CinematicAuthorityClass::VisualizationDerivative,
            CinematicTransformDisposition::BiasedVisualization("oidn-aces-v1".into()),
        ),
        (
            "artistic-sound",
            CinematicArtifactKind::Audio,
            CinematicAuthorityClass::Sound(SoundAuthority::Artistic),
            CinematicTransformDisposition::SoundSynthesis("designed-chirp-v1".into()),
        ),
        (
            "informed-sound",
            CinematicArtifactKind::Audio,
            CinematicAuthorityClass::Sound(SoundAuthority::PhysicallyInformed),
            CinematicTransformDisposition::SoundSynthesis("state-driven-modal-v1".into()),
        ),
        (
            "calibrated-sound",
            CinematicArtifactKind::Audio,
            CinematicAuthorityClass::Sound(SoundAuthority::Calibrated),
            CinematicTransformDisposition::SoundSynthesis("calibrated-modal-v1".into()),
        ),
    ];

    let expected_cases = cases.len();
    let mut record_identities = Vec::new();
    for (label, kind, class, disposition) in cases {
        eprintln!("cinematic-authority phase=admit case={label} class={class:?}");
        let record = admitted(input_for(label, kind, class, disposition));
        let bytes = record.canonical_bytes();
        let decoded = CinematicAuthorityRecord::from_canonical_bytes(&bytes)
            .expect("canonical record must decode");
        assert_eq!(decoded, record, "round trip changed {label}");
        assert_eq!(
            decoded.canonical_bytes(),
            bytes,
            "re-encoding changed {label}"
        );
        assert_eq!(decoded.identity(), record.identity());
        record_identities.push(record.identity());
    }
    record_identities.sort_unstable();
    record_identities.dedup();
    assert_eq!(
        record_identities.len(),
        expected_cases,
        "authority cases aliased"
    );
}

#[test]
fn legal_pipeline_derivations_preserve_sources_and_refuse_promotion() {
    let state = admitted(input_for(
        "pipeline-state",
        CinematicArtifactKind::SimulationState,
        CinematicAuthorityClass::SimulatedState,
        CinematicTransformDisposition::ModelState,
    ));
    let state_before = state.clone();

    let mut render_input = input_for(
        "pipeline-render",
        CinematicArtifactKind::RenderEstimate,
        CinematicAuthorityClass::MonteCarloRender,
        CinematicTransformDisposition::MonteCarloEstimator,
    );
    render_input.source_identity = state.artifact_identity();
    let render = CinematicAuthorityRecord::derive(&state, render_input)
        .expect("simulation state may feed a raw MC render");

    let mut visualization_input = input_for(
        "pipeline-visualization",
        CinematicArtifactKind::MediaDerivative,
        CinematicAuthorityClass::VisualizationDerivative,
        CinematicTransformDisposition::BiasedVisualization("denoise-tone-encode-v1".into()),
    );
    visualization_input.source_identity = render.artifact_identity();
    let visualization = CinematicAuthorityRecord::derive(&render, visualization_input)
        .expect("raw render may feed a disclosed visualization derivative");

    let mut promoted_render = input_for(
        "illegal-promoted-render",
        CinematicArtifactKind::RenderEstimate,
        CinematicAuthorityClass::MonteCarloRender,
        CinematicTransformDisposition::MonteCarloEstimator,
    );
    promoted_render.source_identity = visualization.artifact_identity();
    assert_error_code(
        CinematicAuthorityRecord::derive(&visualization, promoted_render),
        "cinematic-authority-illegal-promotion",
        "biased visualization must not become raw MC evidence",
    );

    let mut wrong_source = input_for(
        "wrong-source",
        CinematicArtifactKind::Visualization,
        CinematicAuthorityClass::VisualizationDerivative,
        CinematicTransformDisposition::BiasedVisualization("display-map-v1".into()),
    );
    wrong_source.source_identity = identity("unrelated-parent");
    assert_error_code(
        CinematicAuthorityRecord::derive(&render, wrong_source),
        "cinematic-authority-source-parent-mismatch",
        "derived record must bind its actual parent",
    );
    assert_eq!(
        state, state_before,
        "derivation mutated the admitted parent"
    );
}

#[test]
fn acoustic_calibration_is_explicit_and_cannot_be_laundered() {
    let mut informed_input = input_for(
        "acoustic-informed",
        CinematicArtifactKind::Audio,
        CinematicAuthorityClass::Sound(SoundAuthority::PhysicallyInformed),
        CinematicTransformDisposition::SoundSynthesis("state-driven-modal-v1".into()),
    );
    let state = admitted(input_for(
        "acoustic-state",
        CinematicArtifactKind::SimulationState,
        CinematicAuthorityClass::SimulatedState,
        CinematicTransformDisposition::ModelState,
    ));
    informed_input.source_identity = state.artifact_identity();
    let informed = CinematicAuthorityRecord::derive(&state, informed_input)
        .expect("state may drive physically-informed synthesis");

    let mut calibrated_input = input_for(
        "acoustic-calibrated",
        CinematicArtifactKind::Audio,
        CinematicAuthorityClass::Sound(SoundAuthority::Calibrated),
        CinematicTransformDisposition::SoundSynthesis("calibrated-modal-v1".into()),
    );
    calibrated_input.source_identity = informed.artifact_identity();
    let calibrated = CinematicAuthorityRecord::derive(&informed, calibrated_input)
        .expect("informed sound plus declared calibration may enter calibrated tier");
    let receipt = calibrated
        .acoustic_calibration()
        .expect("calibrated tier retains the declaration");
    let manifest = calibrated.to_manifest_json();
    assert!(manifest.contains(&receipt.dataset_identity().to_string()));
    assert!(manifest.contains(&receipt.method_identity().to_string()));
    assert!(manifest.contains(&receipt.validity_identity().to_string()));
    assert!(manifest.contains("\"version\":7"));

    let mut missing_receipt = input_for(
        "acoustic-missing-receipt",
        CinematicArtifactKind::Audio,
        CinematicAuthorityClass::Sound(SoundAuthority::Calibrated),
        CinematicTransformDisposition::SoundSynthesis("calibrated-modal-v1".into()),
    );
    missing_receipt.source_identity = informed.artifact_identity();
    missing_receipt.acoustic_calibration = None;
    assert_error_code(
        CinematicAuthorityRecord::derive(&informed, missing_receipt),
        "cinematic-authority-missing-calibration",
        "physically-informed modal audio cannot be relabeled calibrated without a receipt",
    );

    let artistic = admitted(input_for(
        "acoustic-artistic",
        CinematicArtifactKind::Audio,
        CinematicAuthorityClass::Sound(SoundAuthority::Artistic),
        CinematicTransformDisposition::SoundSynthesis("designed-chirp-v1".into()),
    ));
    let mut artistic_promotion = input_for(
        "artistic-promotion",
        CinematicArtifactKind::Audio,
        CinematicAuthorityClass::Sound(SoundAuthority::Calibrated),
        CinematicTransformDisposition::SoundSynthesis("calibrated-modal-v1".into()),
    );
    artistic_promotion.source_identity = artistic.artifact_identity();
    assert_error_code(
        CinematicAuthorityRecord::derive(&artistic, artistic_promotion),
        "cinematic-authority-illegal-promotion",
        "artistic audio cannot be relabeled calibrated",
    );
}

#[test]
fn machine_manifest_and_human_disclosure_share_every_no_claim() {
    let record = admitted(input_for(
        "disclosure-twin",
        CinematicArtifactKind::Visualization,
        CinematicAuthorityClass::VisualizationDerivative,
        CinematicTransformDisposition::BiasedVisualization("denoise-display-v1".into()),
    ));
    let machine = record.to_manifest_json();
    let human = record.human_disclosure();
    assert!(machine.contains("\"source_identity\""));
    assert!(machine.contains("\"unit_contract\""));
    assert!(machine.contains("\"clock\""));
    assert!(machine.contains("\"configuration_identity\""));
    assert!(human.contains(&record.source_identity().to_string()));

    for claim in record.no_claims() {
        assert!(
            machine.contains(claim.code()),
            "machine manifest omitted {}",
            claim.code()
        );
        assert!(
            machine.contains(claim.statement()),
            "machine manifest omitted human statement for {}",
            claim.code()
        );
        assert!(
            human.contains(claim.statement()),
            "human disclosure omitted {}",
            claim.code()
        );
    }
}

#[test]
fn structural_admission_refuses_ambiguous_or_incomplete_records() {
    let base = input_for(
        "structural-refusal",
        CinematicArtifactKind::SimulationState,
        CinematicAuthorityClass::SimulatedState,
        CinematicTransformDisposition::ModelState,
    );

    let mut bad = base.clone();
    bad.schema_version += 1;
    assert_error_code(
        CinematicAuthorityRecord::try_new(bad),
        "cinematic-authority-unsupported-schema",
        "unknown schema",
    );

    let mut bad = base.clone();
    bad.source_identity = ContentHash([0; 32]);
    assert_error_code(
        CinematicAuthorityRecord::try_new(bad),
        "cinematic-authority-missing-identity",
        "zero source identity",
    );

    let mut bad = base.clone();
    bad.configuration_version = 0;
    assert_error_code(
        CinematicAuthorityRecord::try_new(bad),
        "cinematic-authority-invalid-config-version",
        "zero configuration version",
    );

    let mut bad = base.clone();
    bad.transform_name = "Uppercase Is Not A Machine Label".into();
    assert_error_code(
        CinematicAuthorityRecord::try_new(bad),
        "cinematic-authority-invalid-transform-label",
        "invalid transform label grammar",
    );

    let mut bad = base.clone();
    bad.transform_name = "a".repeat(MAX_CINEMATIC_LABEL_BYTES + 1);
    assert_error_code(
        CinematicAuthorityRecord::try_new(bad),
        "cinematic-authority-invalid-transform-label",
        "oversized transform label",
    );

    let mut bad = base.clone();
    bad.artifact_kind = CinematicArtifactKind::Audio;
    assert_error_code(
        CinematicAuthorityRecord::try_new(bad),
        "cinematic-authority-incompatible-artifact-class",
        "audio payload with simulation authority",
    );

    let mut bad = base.clone();
    bad.transform_disposition = CinematicTransformDisposition::MonteCarloEstimator;
    assert_error_code(
        CinematicAuthorityRecord::try_new(bad),
        "cinematic-authority-incompatible-transform-disposition",
        "raw estimator disposition on model state",
    );

    let mut bad = base.clone();
    let missing = bad.no_claims.pop().expect("required disclosures exist");
    let error = CinematicAuthorityRecord::try_new(bad).expect_err("missing no-claim must refuse");
    assert_eq!(error, CinematicAuthorityError::MissingNoClaim(missing));

    let mut bad = base.clone();
    bad.no_claims.push(bad.no_claims[0]);
    assert_error_code(
        CinematicAuthorityRecord::try_new(bad),
        "cinematic-authority-duplicate-no-claim",
        "duplicate no-claim",
    );

    let mut bad = base;
    bad.no_claims = vec![CinematicNoClaim::PhysicalStoppingTime; MAX_CINEMATIC_NO_CLAIMS + 1];
    assert_error_code(
        CinematicAuthorityRecord::try_new(bad),
        "cinematic-authority-too-many-no-claims",
        "bounded no-claim vector",
    );
}

#[test]
fn clocks_refuse_undefined_rates_ranges_and_noncanonical_timeless_values() {
    assert_eq!(
        CinematicClock::try_new(CinematicClockDomain::Video, 0, 1, 0, 1),
        Err(CinematicAuthorityError::InvalidClockRate)
    );
    assert_eq!(
        CinematicClock::try_new(CinematicClockDomain::Video, 24, 1, 2, 1),
        Err(CinematicAuthorityError::InvalidClockRange)
    );
    assert_eq!(
        CinematicClock::try_new(CinematicClockDomain::Timeless, 1, 1, 0, 1),
        Err(CinematicAuthorityError::InvalidTimelessClock)
    );
    assert_eq!(
        CinematicClock::timeless().domain(),
        CinematicClockDomain::Timeless
    );
}

#[test]
fn canonical_decoder_refuses_truncation_trailing_data_and_tag_mutations() {
    let record = admitted(input_for(
        "codec-mutations",
        CinematicArtifactKind::Visualization,
        CinematicAuthorityClass::VisualizationDerivative,
        CinematicTransformDisposition::BiasedVisualization("denoise-display-v1".into()),
    ));
    let bytes = record.canonical_bytes();

    for prefix_len in 0..bytes.len() {
        let result = CinematicAuthorityRecord::from_canonical_bytes(&bytes[..prefix_len]);
        assert!(
            result.is_err(),
            "truncated prefix {prefix_len} was admitted"
        );
    }
    let mut mutated = bytes.clone();
    mutated.push(0);
    assert_error_code(
        CinematicAuthorityRecord::from_canonical_bytes(&mutated),
        "cinematic-authority-trailing-bytes",
        "trailing bytes",
    );

    let mutation_cases = [
        (8, 2, "cinematic-authority-unsupported-schema", "schema"),
        (
            10,
            99,
            "cinematic-authority-unknown-artifact-kind",
            "artifact kind",
        ),
        (
            11,
            99,
            "cinematic-authority-unknown-class",
            "authority class",
        ),
        (
            12,
            1,
            "cinematic-authority-unexpected-sound-class",
            "sound tier",
        ),
        (
            145,
            99,
            "cinematic-authority-unknown-unit-contract",
            "unit contract",
        ),
        (
            146,
            99,
            "cinematic-authority-unknown-clock-domain",
            "clock domain",
        ),
    ];
    for (offset, value, code, context) in mutation_cases {
        let mut mutated = bytes.clone();
        mutated[offset] = value;
        assert_error_code(
            CinematicAuthorityRecord::from_canonical_bytes(&mutated),
            code,
            context,
        );
    }

    let mut mutated = bytes.clone();
    let offset = disposition_offset(&mutated);
    mutated[offset] = 99;
    assert_error_code(
        CinematicAuthorityRecord::from_canonical_bytes(&mutated),
        "cinematic-authority-unknown-disposition",
        "disposition tag",
    );

    let mut mutated = bytes.clone();
    let offset = calibration_presence_offset(&mutated);
    mutated[offset] = 99;
    assert_error_code(
        CinematicAuthorityRecord::from_canonical_bytes(&mutated),
        "cinematic-authority-unknown-calibration-presence",
        "calibration option tag",
    );

    let mut mutated = bytes;
    let last = mutated.len() - 1;
    mutated[last] = 99;
    assert_error_code(
        CinematicAuthorityRecord::from_canonical_bytes(&mutated),
        "cinematic-authority-unknown-no-claim",
        "no-claim tag",
    );
}

#[test]
fn canonical_mutations_cannot_strip_source_or_required_disclosure() {
    let record = admitted(input_for(
        "negative-twin",
        CinematicArtifactKind::Visualization,
        CinematicAuthorityClass::VisualizationDerivative,
        CinematicTransformDisposition::BiasedVisualization("denoise-display-v1".into()),
    ));
    let bytes = record.canonical_bytes();

    let mut zero_source = bytes.clone();
    zero_source[45..77].fill(0);
    assert_error_code(
        CinematicAuthorityRecord::from_canonical_bytes(&zero_source),
        "cinematic-authority-missing-identity",
        "negative twin with stripped source identity",
    );

    let mut missing_disclosure = bytes;
    let count_offset = missing_disclosure.len() - record.no_claims().len() - 2;
    let reduced_count = u16::try_from(record.no_claims().len() - 1).expect("small claim count");
    missing_disclosure[count_offset..count_offset + 2]
        .copy_from_slice(&reduced_count.to_le_bytes());
    missing_disclosure.pop();
    assert_error_code(
        CinematicAuthorityRecord::from_canonical_bytes(&missing_disclosure),
        "cinematic-authority-missing-no-claim",
        "negative twin with stripped mandatory disclosure",
    );
}
