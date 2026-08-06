//! Focused Gauntlet battery for the frame-sequence artifact contract.
//!
//! G0 cases exercise input validation, exact resource boundaries, lineage,
//! and transactional registration. G3/G5 cases exercise order independence,
//! relocation, snapshot/resume equivalence, and byte-identical finalization.

// Each scenario keeps its full before/refusal/after transaction narrative in
// one test so a failure names the violated contract rather than a helper step.
#![allow(clippy::too_many_lines)]

use std::collections::BTreeMap;

use fs_blake3::DomainHasher;
use fs_img::{
    Channel, ContentHash, ExpectedFrameArtifact, FRAME_SEQUENCE_MANIFEST_VERSION,
    FrameArtifactDescriptor, FrameArtifactFileState, FrameArtifactFormat, FrameArtifactKey,
    FrameArtifactObservation, FrameArtifactRole, FrameChannel, FrameChannelType,
    FrameSamplingStats, FrameSequenceContext, FrameSequenceError, FrameSequenceLimits,
    FrameSequenceManifest, FrameSequenceState, PixelType, PngColor, RegistrationOutcome, read_exr,
    read_png, write_exr, write_png8,
};

fn content_hash(label: &[u8]) -> ContentHash {
    FrameArtifactFileState::from_bytes(label)
        .expect("small fixture hash")
        .content_hash()
}

fn snapshot_identity(bytes: &[u8]) -> ContentHash {
    let domain =
        format!("org.frankensim.fs-img.frame-sequence-manifest.v{FRAME_SEQUENCE_MANIFEST_VERSION}");
    let mut hasher = DomainHasher::new(&domain);
    hasher.update(bytes);
    hasher.finalize()
}

fn context() -> FrameSequenceContext {
    FrameSequenceContext::try_new(
        content_hash(b"sequence-battery-shot"),
        content_hash(b"sequence-battery-trajectory"),
        content_hash(b"sequence-battery-render-config"),
        content_hash(b"sequence-battery-scene"),
        content_hash(b"sequence-battery-build"),
        content_hash(b"sequence-battery-profile"),
    )
    .expect("fixture identities are non-placeholder")
}

fn limits(max_output_bytes: u64) -> FrameSequenceLimits {
    FrameSequenceLimits::try_new(64, 16, 512, 1 << 20, max_output_bytes)
        .expect("fixture limits are nonzero")
}

fn channel(name: &str, sample_type: FrameChannelType) -> FrameChannel {
    FrameChannel::try_new(name, sample_type).expect("fixture channel")
}

fn exr_channels() -> Vec<FrameChannel> {
    // Deliberately noncanonical input order. The descriptor must agree with
    // the EXR writer's alphabetical stored order, independent of callers.
    vec![
        channel("R", FrameChannelType::Float32),
        channel("B", FrameChannelType::Float32),
        channel("G", FrameChannelType::Float32),
    ]
}

fn png8_channels() -> Vec<FrameChannel> {
    // Deliberately noncanonical input order. PNG has a fixed packed order.
    vec![
        channel("B", FrameChannelType::Uint8),
        channel("R", FrameChannelType::Uint8),
        channel("G", FrameChannelType::Uint8),
    ]
}

#[allow(clippy::too_many_arguments)]
fn descriptor(
    frame_index: u64,
    segment_index: u32,
    role: FrameArtifactRole,
    frame_time_s: f64,
    format: FrameArtifactFormat,
    width: u32,
    height: u32,
    channels: Vec<FrameChannel>,
    sampling: FrameSamplingStats,
) -> FrameArtifactDescriptor {
    FrameArtifactDescriptor::try_new(
        frame_index,
        segment_index,
        role,
        frame_time_s,
        format,
        width,
        height,
        channels,
        sampling,
    )
    .expect("fixture descriptor")
}

fn raw_descriptor(
    frame_index: u64,
    segment_index: u32,
    frame_time_s: f64,
) -> FrameArtifactDescriptor {
    descriptor(
        frame_index,
        segment_index,
        FrameArtifactRole::RawMaster,
        frame_time_s,
        FrameArtifactFormat::OpenExr,
        2,
        1,
        exr_channels(),
        FrameSamplingStats::Uniform { spp: 8 },
    )
}

#[derive(Clone, Debug)]
struct FixtureArtifact {
    descriptor: FrameArtifactDescriptor,
    bytes: Vec<u8>,
    source: Option<FrameArtifactKey>,
    source_content_hash: Option<ContentHash>,
}

impl FixtureArtifact {
    fn key(&self) -> FrameArtifactKey {
        self.descriptor.key()
    }

    fn file_state(&self) -> FrameArtifactFileState {
        FrameArtifactFileState::from_bytes(&self.bytes).expect("small fixture file state")
    }

    fn expected(&self) -> ExpectedFrameArtifact {
        ExpectedFrameArtifact::try_new(
            self.descriptor.clone(),
            u64::try_from(self.bytes.len()).expect("fixture length fits u64"),
            self.source,
        )
        .expect("fixture lineage and reservation")
    }
}

fn tiny_mixed_sequence() -> (FrameSequenceContext, Vec<FixtureArtifact>) {
    let context = context();
    let mut artifacts = Vec::new();

    for frame_index in 0_u64..2 {
        let segment_index = u32::try_from(frame_index).expect("two fixture frames");
        let frame_time_s = frame_index as f64 / 24.0;
        let base = frame_index as f32 * 0.125;

        let raw = raw_descriptor(frame_index, segment_index, frame_time_s);
        let raw_bytes = write_exr(
            2,
            1,
            &[
                Channel {
                    name: "R".to_owned(),
                    ty: PixelType::Float,
                    data: vec![base + 0.10, base + 0.20],
                },
                Channel {
                    name: "B".to_owned(),
                    ty: PixelType::Float,
                    data: vec![base + 0.30, base + 0.40],
                },
                Channel {
                    name: "G".to_owned(),
                    ty: PixelType::Float,
                    data: vec![base + 0.50, base + 0.60],
                },
            ],
        )
        .expect("encode raw EXR fixture");
        read_exr(&raw_bytes).expect("raw EXR fixture must decode");
        let raw_key = raw.key();
        let raw_hash = FrameArtifactFileState::from_bytes(&raw_bytes)
            .expect("raw file state")
            .content_hash();
        artifacts.push(FixtureArtifact {
            descriptor: raw,
            bytes: raw_bytes,
            source: None,
            source_content_hash: None,
        });

        let denoised = descriptor(
            frame_index,
            segment_index,
            FrameArtifactRole::DenoisedIntermediate,
            frame_time_s,
            FrameArtifactFormat::OpenExr,
            2,
            1,
            exr_channels(),
            FrameSamplingStats::Uniform { spp: 8 },
        );
        let denoised_bytes = write_exr(
            2,
            1,
            &[
                Channel {
                    name: "G".to_owned(),
                    ty: PixelType::Float,
                    data: vec![base + 0.45, base + 0.55],
                },
                Channel {
                    name: "R".to_owned(),
                    ty: PixelType::Float,
                    data: vec![base + 0.15, base + 0.25],
                },
                Channel {
                    name: "B".to_owned(),
                    ty: PixelType::Float,
                    data: vec![base + 0.35, base + 0.45],
                },
            ],
        )
        .expect("encode derived EXR fixture");
        read_exr(&denoised_bytes).expect("derived EXR fixture must decode");
        let denoised_key = denoised.key();
        let denoised_hash = FrameArtifactFileState::from_bytes(&denoised_bytes)
            .expect("denoised file state")
            .content_hash();
        artifacts.push(FixtureArtifact {
            descriptor: denoised,
            bytes: denoised_bytes,
            source: Some(raw_key),
            source_content_hash: Some(raw_hash),
        });

        let preview = descriptor(
            frame_index,
            segment_index,
            FrameArtifactRole::DisplayPreview,
            frame_time_s,
            FrameArtifactFormat::Png8,
            2,
            1,
            png8_channels(),
            FrameSamplingStats::Uniform { spp: 8 },
        );
        let preview_bytes = write_png8(
            2,
            1,
            PngColor::Rgb,
            &[
                20 + u8::try_from(frame_index).expect("fixture frame"),
                40,
                60,
                80,
                100,
                120,
            ],
        )
        .expect("encode PNG fixture");
        read_png(&preview_bytes).expect("PNG fixture must decode");
        artifacts.push(FixtureArtifact {
            descriptor: preview,
            bytes: preview_bytes,
            source: Some(denoised_key),
            source_content_hash: Some(denoised_hash),
        });
    }

    (context, artifacts)
}

fn total_reserved(artifacts: &[FixtureArtifact]) -> u64 {
    artifacts.iter().fold(0_u64, |total, artifact| {
        total
            .checked_add(u64::try_from(artifact.bytes.len()).expect("fixture length fits u64"))
            .expect("fixture total fits u64")
    })
}

fn manifest_for(
    context: FrameSequenceContext,
    artifacts: &[FixtureArtifact],
    reverse_expectations: bool,
) -> FrameSequenceManifest {
    let mut expected: Vec<_> = artifacts.iter().map(FixtureArtifact::expected).collect();
    if reverse_expectations {
        expected.reverse();
    }
    let reserved = total_reserved(artifacts);
    FrameSequenceManifest::try_new(context, expected, limits(reserved), reserved)
        .expect("admit tiny fixture sequence")
}

fn path_for(manifest: &FrameSequenceManifest, key: FrameArtifactKey) -> String {
    manifest
        .entries()
        .iter()
        .find(|entry| entry.descriptor().key() == key)
        .unwrap_or_else(|| panic!("missing fixture key {key:?}"))
        .relative_path()
        .to_owned()
}

fn planned_path(
    context: FrameSequenceContext,
    expected: Vec<ExpectedFrameArtifact>,
    key: FrameArtifactKey,
) -> String {
    let reserved = expected.iter().fold(0_u64, |total, artifact| {
        total
            .checked_add(artifact.max_bytes())
            .expect("small planned reservation")
    });
    let manifest = FrameSequenceManifest::try_new(context, expected, limits(reserved), reserved)
        .expect("admit path-identity fixture");
    path_for(&manifest, key)
}

fn register_indices(
    manifest: &mut FrameSequenceManifest,
    artifacts: &[FixtureArtifact],
    indices: impl IntoIterator<Item = usize>,
) {
    for index in indices {
        let artifact = &artifacts[index];
        let path = path_for(manifest, artifact.key());
        let outcome = manifest
            .register_artifact_bytes(
                &path,
                artifact.descriptor.clone(),
                manifest.context().profile_id(),
                &artifact.bytes,
                artifact.source_content_hash,
            )
            .unwrap_or_else(|error| panic!("register {path}: {error}"));
        assert_eq!(outcome, RegistrationOutcome::Recorded, "register {path}");
    }
}

fn file_states(
    manifest: &FrameSequenceManifest,
    artifacts: &[FixtureArtifact],
) -> BTreeMap<String, FrameArtifactFileState> {
    artifacts
        .iter()
        .map(|artifact| (path_for(manifest, artifact.key()), artifact.file_state()))
        .collect()
}

#[test]
fn g0_zero_invalid_canonical_paths_and_segment_distinction() {
    assert!(matches!(
        FrameSequenceLimits::try_new(0, 1, 1, 1, 1),
        Err(FrameSequenceError::InvalidLimit {
            field: "max_artifacts"
        })
    ));
    assert!(matches!(
        FrameSequenceContext::try_new(
            ContentHash([0; 32]),
            content_hash(b"trajectory"),
            content_hash(b"render"),
            content_hash(b"scene"),
            content_hash(b"build"),
            content_hash(b"profile"),
        ),
        Err(FrameSequenceError::PlaceholderIdentity { field: "shot_id" })
    ));
    assert!(matches!(
        FrameSequenceManifest::try_new(context(), Vec::new(), limits(1), 1),
        Err(FrameSequenceError::EmptySequence)
    ));

    let bad_time = FrameArtifactDescriptor::try_new(
        7,
        0,
        FrameArtifactRole::RawMaster,
        f64::NAN,
        FrameArtifactFormat::OpenExr,
        1,
        1,
        vec![channel("R", FrameChannelType::Float32)],
        FrameSamplingStats::Uniform { spp: 1 },
    );
    assert!(matches!(
        bad_time,
        Err(FrameSequenceError::InvalidFrameTime { frame_index: 7 })
    ));
    assert!(matches!(
        FrameArtifactDescriptor::try_new(
            7,
            0,
            FrameArtifactRole::RawMaster,
            0.0,
            FrameArtifactFormat::OpenExr,
            0,
            1,
            vec![channel("R", FrameChannelType::Float32)],
            FrameSamplingStats::Uniform { spp: 1 },
        ),
        Err(FrameSequenceError::InvalidDimensions { frame_index: 7 })
    ));
    assert!(matches!(
        FrameArtifactDescriptor::try_new(
            7,
            0,
            FrameArtifactRole::DisplayPreview,
            0.0,
            FrameArtifactFormat::Png8,
            1,
            1,
            vec![
                channel("X", FrameChannelType::Uint8),
                channel("Y", FrameChannelType::Uint8),
                channel("Z", FrameChannelType::Uint8),
            ],
            FrameSamplingStats::Uniform { spp: 1 },
        ),
        Err(FrameSequenceError::InvalidChannelSet { frame_index: 7 })
    ));

    let positive_zero = raw_descriptor(7, 0, 0.0);
    let negative_zero = raw_descriptor(7, 0, -0.0);
    assert_eq!(
        positive_zero, negative_zero,
        "signed zero must canonicalize"
    );
    assert_eq!(negative_zero.frame_time_bits(), 0);
    assert!(matches!(
        ExpectedFrameArtifact::try_new(positive_zero.clone(), 0, None),
        Err(FrameSequenceError::InvalidArtifactLimit { .. })
    ));

    let segment_one = raw_descriptor(7, 1, 0.0);
    let expected_a = ExpectedFrameArtifact::try_new(positive_zero.clone(), 1, None).unwrap();
    let expected_b = ExpectedFrameArtifact::try_new(segment_one.clone(), 1, None).unwrap();
    assert!(matches!(
        FrameSequenceManifest::try_new(
            context(),
            vec![expected_a.clone(), expected_a.clone()],
            limits(2),
            2,
        ),
        Err(FrameSequenceError::DuplicateExpectedArtifact)
    ));
    let manifest =
        FrameSequenceManifest::try_new(context(), vec![expected_b, expected_a], limits(2), 2)
            .expect("segments are distinct expected artifacts");
    let path_zero = path_for(&manifest, positive_zero.key());
    let path_one = path_for(&manifest, segment_one.key());
    assert_ne!(path_zero, path_one, "segments must never alias one file");
    assert!(path_zero.contains("-segment-0000000000-"), "{path_zero}");
    assert!(path_one.contains("-segment-0000000001-"), "{path_one}");
    for path in [&path_zero, &path_one] {
        assert!(
            !path.starts_with('/'),
            "path must be root-independent: {path}"
        );
        assert!(!path.contains("//"), "path must be normalized: {path}");
        assert!(
            path.split('/').all(|part| !matches!(part, "" | "." | "..")),
            "path contains an unsafe component: {path}"
        );
    }

    let base_context = context();
    let changed_context = FrameSequenceContext::try_new(
        base_context.shot_id(),
        base_context.trajectory_id(),
        base_context.render_config_id(),
        base_context.scene_id(),
        content_hash(b"different-build-identity"),
        base_context.profile_id(),
    )
    .unwrap();
    let base_path = planned_path(
        base_context,
        vec![ExpectedFrameArtifact::try_new(positive_zero.clone(), 1, None).unwrap()],
        positive_zero.key(),
    );
    let changed_context_path = planned_path(
        changed_context,
        vec![ExpectedFrameArtifact::try_new(positive_zero.clone(), 1, None).unwrap()],
        positive_zero.key(),
    );
    assert_ne!(
        base_path, changed_context_path,
        "the full sequence context must participate in path identity"
    );

    let changed_descriptor = descriptor(
        7,
        0,
        FrameArtifactRole::RawMaster,
        0.0,
        FrameArtifactFormat::OpenExr,
        3,
        1,
        exr_channels(),
        FrameSamplingStats::Uniform { spp: 8 },
    );
    let changed_descriptor_path = planned_path(
        base_context,
        vec![ExpectedFrameArtifact::try_new(changed_descriptor.clone(), 1, None).unwrap()],
        changed_descriptor.key(),
    );
    assert_ne!(
        base_path, changed_descriptor_path,
        "descriptor semantics must participate in path identity"
    );

    let source_raw = raw_descriptor(13, 0, 1.0);
    let source_denoised = descriptor(
        13,
        0,
        FrameArtifactRole::DenoisedIntermediate,
        1.0,
        FrameArtifactFormat::OpenExr,
        2,
        1,
        exr_channels(),
        FrameSamplingStats::Uniform { spp: 8 },
    );
    let overlay = descriptor(
        13,
        0,
        FrameArtifactRole::ScientificOverlay,
        1.0,
        FrameArtifactFormat::OpenExr,
        2,
        1,
        exr_channels(),
        FrameSamplingStats::Uniform { spp: 8 },
    );
    let sourced_from_raw = planned_path(
        base_context,
        vec![
            ExpectedFrameArtifact::try_new(source_raw.clone(), 1, None).unwrap(),
            ExpectedFrameArtifact::try_new(overlay.clone(), 1, Some(source_raw.key())).unwrap(),
        ],
        overlay.key(),
    );
    let sourced_from_denoised = planned_path(
        base_context,
        vec![
            ExpectedFrameArtifact::try_new(source_raw.clone(), 1, None).unwrap(),
            ExpectedFrameArtifact::try_new(source_denoised.clone(), 1, Some(source_raw.key()))
                .unwrap(),
            ExpectedFrameArtifact::try_new(overlay.clone(), 1, Some(source_denoised.key()))
                .unwrap(),
        ],
        overlay.key(),
    );
    assert_ne!(
        sourced_from_raw, sourced_from_denoised,
        "the same derived descriptor must not alias across valid source roles"
    );
}

#[test]
fn g5_real_mixed_sequence_is_order_independent_and_relocatable() {
    let (context, artifacts) = tiny_mixed_sequence();
    assert_eq!(artifacts.len(), 6, "two frames times three artifact roles");
    assert_eq!(
        artifacts[0]
            .descriptor
            .channels()
            .iter()
            .map(FrameChannel::name)
            .collect::<Vec<_>>(),
        ["B", "G", "R"],
        "EXR descriptor order must match the writer's canonical channel order"
    );
    assert_eq!(
        artifacts[2]
            .descriptor
            .channels()
            .iter()
            .map(FrameChannel::name)
            .collect::<Vec<_>>(),
        ["R", "G", "B"],
        "PNG descriptors must normalize to packed channel order"
    );

    let mut forward = manifest_for(context, &artifacts, false);
    let mut reverse = manifest_for(context, &artifacts, true);
    register_indices(&mut forward, &artifacts, 0..artifacts.len());
    register_indices(&mut reverse, &artifacts, (0..artifacts.len()).rev());
    let states = file_states(&forward, &artifacts);
    let forward_seal = forward
        .finalize_with(|| true, |path| states.get(path).copied())
        .expect("finalize forward registration");
    let reverse_states = file_states(&reverse, &artifacts);
    let reverse_seal = reverse
        .finalize_with(|| true, |path| reverse_states.get(path).copied())
        .expect("finalize reverse registration");

    assert_eq!(forward_seal.bytes(), reverse_seal.bytes());
    assert_eq!(forward_seal.identity(), reverse_seal.identity());
    assert_eq!(forward_seal.artifact_count(), 6);
    assert_eq!(forward_seal.output_bytes(), total_reserved(&artifacts));

    for root in ["/render-node-a/job-19", "/archive-volume/relocated/job-19"] {
        let relocated: BTreeMap<_, _> = states
            .iter()
            .map(|(relative, state)| (format!("{root}/{relative}"), *state))
            .collect();
        forward
            .audit_with(
                || true,
                |relative| relocated.get(&format!("{root}/{relative}")).copied(),
            )
            .unwrap_or_else(|error| panic!("audit after relocation to {root}: {error}"));
    }
    assert_eq!(
        forward.audit_with(|| false, |_| None),
        Err(FrameSequenceError::Cancelled),
        "finalized audit must expose cancellation before observation"
    );
}

#[test]
fn g0_descriptor_mutations_and_duplicate_conflicts_are_transactional() {
    let context = context();
    let raw = raw_descriptor(11, 3, 0.5);
    let raw_bytes = b"raw-source-bytes";
    let raw_hash = content_hash(raw_bytes);
    let overlay = descriptor(
        11,
        3,
        FrameArtifactRole::ScientificOverlay,
        0.5,
        FrameArtifactFormat::OpenExr,
        2,
        1,
        exr_channels(),
        FrameSamplingStats::Uniform { spp: 8 },
    );
    let expected = vec![
        ExpectedFrameArtifact::try_new(raw.clone(), 64, None).unwrap(),
        ExpectedFrameArtifact::try_new(overlay.clone(), 64, Some(raw.key())).unwrap(),
    ];
    let mut manifest = FrameSequenceManifest::try_new(context, expected, limits(128), 128).unwrap();
    let overlay_path = path_for(&manifest, overlay.key());
    let payload = b"correct-overlay";

    let mutations = [
        (
            "artifact key",
            descriptor(
                12,
                3,
                FrameArtifactRole::ScientificOverlay,
                0.5,
                FrameArtifactFormat::OpenExr,
                2,
                1,
                exr_channels(),
                FrameSamplingStats::Uniform { spp: 8 },
            ),
        ),
        (
            "artifact key",
            descriptor(
                11,
                4,
                FrameArtifactRole::ScientificOverlay,
                0.5,
                FrameArtifactFormat::OpenExr,
                2,
                1,
                exr_channels(),
                FrameSamplingStats::Uniform { spp: 8 },
            ),
        ),
        (
            "frame time",
            descriptor(
                11,
                3,
                FrameArtifactRole::ScientificOverlay,
                0.75,
                FrameArtifactFormat::OpenExr,
                2,
                1,
                exr_channels(),
                FrameSamplingStats::Uniform { spp: 8 },
            ),
        ),
        (
            "format",
            descriptor(
                11,
                3,
                FrameArtifactRole::ScientificOverlay,
                0.5,
                FrameArtifactFormat::Png8,
                2,
                1,
                png8_channels(),
                FrameSamplingStats::Uniform { spp: 8 },
            ),
        ),
        (
            "dimensions",
            descriptor(
                11,
                3,
                FrameArtifactRole::ScientificOverlay,
                0.5,
                FrameArtifactFormat::OpenExr,
                3,
                1,
                exr_channels(),
                FrameSamplingStats::Uniform { spp: 8 },
            ),
        ),
        (
            "channels",
            descriptor(
                11,
                3,
                FrameArtifactRole::ScientificOverlay,
                0.5,
                FrameArtifactFormat::OpenExr,
                2,
                1,
                vec![channel("depth.Z", FrameChannelType::Float32)],
                FrameSamplingStats::Uniform { spp: 8 },
            ),
        ),
        (
            "sampling statistics",
            descriptor(
                11,
                3,
                FrameArtifactRole::ScientificOverlay,
                0.5,
                FrameArtifactFormat::OpenExr,
                2,
                1,
                exr_channels(),
                FrameSamplingStats::Uniform { spp: 16 },
            ),
        ),
    ];

    for (expected_field, mutation) in mutations {
        let before = manifest.clone();
        let error = manifest
            .register_artifact_bytes(
                &overlay_path,
                mutation,
                context.profile_id(),
                payload,
                Some(raw_hash),
            )
            .expect_err("mutated descriptor must refuse");
        assert!(
            matches!(
                error,
                FrameSequenceError::DescriptorMismatch { field, .. } if field == expected_field
            ),
            "expected descriptor field {expected_field:?}, got {error:?}"
        );
        assert_eq!(manifest, before, "{expected_field} refusal mutated state");
    }

    let before_profile_refusal = manifest.clone();
    let error = manifest
        .register_artifact_bytes(
            &overlay_path,
            overlay.clone(),
            content_hash(b"wrong-profile"),
            payload,
            Some(raw_hash),
        )
        .expect_err("wrong profile must refuse");
    assert!(matches!(
        error,
        FrameSequenceError::DescriptorMismatch {
            field: "profile identity",
            ..
        }
    ));
    assert_eq!(manifest, before_profile_refusal);

    let before_unexpected = manifest.clone();
    let mut unexpected_hash_polls = 0_u32;
    let unexpected = manifest
        .register_artifact_bytes_with_poll(
            "raw/../caller-selected.exr",
            overlay.clone(),
            context.profile_id(),
            payload,
            Some(raw_hash),
            || {
                unexpected_hash_polls += 1;
                true
            },
        )
        .expect_err("caller-selected or traversal-like names must not enter the inventory");
    assert!(matches!(
        unexpected,
        FrameSequenceError::UnexpectedArtifact { .. }
    ));
    assert_eq!(manifest, before_unexpected);
    assert_eq!(
        unexpected_hash_polls, 0,
        "unknown paths must refuse before payload hashing"
    );

    assert_eq!(
        manifest
            .register_artifact_bytes(
                &overlay_path,
                overlay.clone(),
                context.profile_id(),
                payload,
                Some(raw_hash),
            )
            .unwrap(),
        RegistrationOutcome::Recorded
    );
    let after_record = manifest.clone();
    assert_eq!(
        manifest
            .register_artifact_bytes(
                &overlay_path,
                overlay.clone(),
                context.profile_id(),
                payload,
                Some(raw_hash),
            )
            .unwrap(),
        RegistrationOutcome::AlreadyRecorded
    );
    assert_eq!(manifest, after_record, "exact retry must be a no-op");

    let conflict = manifest
        .register_artifact_bytes(
            &overlay_path,
            overlay,
            context.profile_id(),
            b"conflicting-overlay",
            Some(raw_hash),
        )
        .expect_err("conflicting retry must refuse");
    assert!(matches!(
        conflict,
        FrameSequenceError::ConflictingDuplicate { .. }
    ));
    assert_eq!(manifest, after_record, "conflict must not alter completion");
}

#[test]
fn g0_storage_admission_and_artifact_limits_are_exact_boundaries() {
    const MAX_BYTES: u64 = 16;
    let context = context();
    let raw = raw_descriptor(0, 0, 0.0);
    let expected = || ExpectedFrameArtifact::try_new(raw.clone(), MAX_BYTES, None).unwrap();

    let mut exact =
        FrameSequenceManifest::try_new(context, vec![expected()], limits(MAX_BYTES), MAX_BYTES)
            .expect("exact output and available-space boundaries must admit");
    assert_eq!(exact.remaining_reserved_bytes().unwrap(), MAX_BYTES);

    let finalized_manifest_bytes = exact
        .finalized_manifest_bytes()
        .expect("compute eventual completed-manifest size");
    let exact_manifest_limits =
        FrameSequenceLimits::try_new(64, 16, 512, finalized_manifest_bytes, MAX_BYTES).unwrap();
    FrameSequenceManifest::try_new(context, vec![expected()], exact_manifest_limits, MAX_BYTES)
        .expect("exact finalized-manifest byte ceiling must admit");
    let one_short_manifest_limits =
        FrameSequenceLimits::try_new(64, 16, 512, finalized_manifest_bytes - 1, MAX_BYTES).unwrap();
    let manifest_error = FrameSequenceManifest::try_new(
        context,
        vec![expected()],
        one_short_manifest_limits,
        MAX_BYTES,
    )
    .expect_err("one byte below eventual completed-manifest size must refuse");
    assert!(matches!(
        manifest_error,
        FrameSequenceError::ResourceLimit {
            resource: "manifest bytes",
            requested,
            limit,
        } if requested == finalized_manifest_bytes && limit == finalized_manifest_bytes - 1
    ));

    let output_error =
        FrameSequenceManifest::try_new(context, vec![expected()], limits(MAX_BYTES - 1), MAX_BYTES)
            .expect_err("reservation one byte over output ceiling must refuse");
    assert!(matches!(
        output_error,
        FrameSequenceError::ResourceLimit {
            resource: "reserved output bytes",
            requested: MAX_BYTES,
            limit,
        } if limit == MAX_BYTES - 1
    ));

    let available_error =
        FrameSequenceManifest::try_new(context, vec![expected()], limits(MAX_BYTES), MAX_BYTES - 1)
            .expect_err("reservation one byte over available space must refuse");
    assert!(matches!(
        available_error,
        FrameSequenceError::ResourceLimit {
            resource: "available output bytes",
            requested: MAX_BYTES,
            limit,
        } if limit == MAX_BYTES - 1
    ));

    let path = path_for(&exact, raw.key());
    assert_eq!(
        exact
            .register_artifact_bytes(
                &path,
                raw.clone(),
                context.profile_id(),
                &[0x5a; MAX_BYTES as usize],
                None,
            )
            .unwrap(),
        RegistrationOutcome::Recorded
    );
    assert_eq!(exact.completed_bytes(), MAX_BYTES);
    assert_eq!(exact.remaining_reserved_bytes().unwrap(), 0);

    let mut oversized =
        FrameSequenceManifest::try_new(context, vec![expected()], limits(MAX_BYTES), MAX_BYTES)
            .unwrap();
    let before = oversized.clone();
    let mut oversize_hash_polls = 0_u32;
    let error = oversized
        .register_artifact_bytes_with_poll(
            &path_for(&oversized, raw.key()),
            raw,
            context.profile_id(),
            &[0xa5; MAX_BYTES as usize + 1],
            None,
            || {
                oversize_hash_polls += 1;
                true
            },
        )
        .expect_err("artifact one byte over its reservation must refuse");
    assert!(matches!(
        error,
        FrameSequenceError::ResourceLimit {
            resource: "artifact bytes",
            requested,
            limit: MAX_BYTES,
        } if requested == MAX_BYTES + 1
    ));
    assert_eq!(oversized, before, "oversize refusal must be transactional");
    assert_eq!(
        oversize_hash_polls, 0,
        "oversized payloads must refuse before hashing"
    );
}

#[test]
fn g5_cancel_snapshot_decode_resume_matches_uninterrupted() {
    let (context, artifacts) = tiny_mixed_sequence();
    let reserved = total_reserved(&artifacts);
    let mut resumed = manifest_for(context, &artifacts, false);

    let first = &artifacts[0];
    let first_path = path_for(&resumed, first.key());
    let before_hash_cancel = resumed.clone();
    let hash_cancel = resumed
        .register_artifact_bytes_with_poll(
            &first_path,
            first.descriptor.clone(),
            context.profile_id(),
            &first.bytes,
            first.source_content_hash,
            || false,
        )
        .expect_err("artifact hashing cancellation must refuse registration");
    assert_eq!(hash_cancel, FrameSequenceError::Cancelled);
    assert_eq!(
        resumed, before_hash_cancel,
        "cancelled artifact hashing mutated completion state"
    );
    register_indices(&mut resumed, &artifacts, 0..3);

    let before_cancel = resumed.clone();
    let mut polls = 0_u32;
    let cancelled = resumed
        .snapshot_with_poll(|| {
            polls += 1;
            polls < 3
        })
        .expect_err("snapshot cancellation must be observed");
    assert_eq!(cancelled, FrameSequenceError::Cancelled);
    assert_eq!(resumed, before_cancel, "cancelled snapshot mutated state");

    let snapshot = resumed.snapshot().expect("create resumable snapshot");
    assert_eq!(snapshot.state(), FrameSequenceState::Incomplete);
    assert_eq!(snapshot.identity(), snapshot_identity(snapshot.bytes()));
    let remaining = resumed.remaining_reserved_bytes().unwrap();
    assert!(remaining < reserved && remaining > 0);

    let wrong_pin = FrameSequenceManifest::decode_snapshot(
        snapshot.bytes(),
        content_hash(b"wrong-independent-pin"),
        limits(reserved),
        remaining,
    )
    .expect_err("wrong identity pin must refuse");
    assert!(matches!(
        wrong_pin,
        FrameSequenceError::IdentityMismatch { .. }
    ));

    for prefix_len in 0..snapshot.bytes().len() {
        let truncated = &snapshot.bytes()[..prefix_len];
        let truncation_error = FrameSequenceManifest::decode_snapshot(
            truncated,
            snapshot_identity(truncated),
            limits(reserved),
            remaining,
        )
        .unwrap_err();
        assert_eq!(
            truncation_error,
            FrameSequenceError::Truncated,
            "correctly pinned prefix {prefix_len} of {} bytes must refuse as truncated",
            snapshot.bytes().len()
        );
    }

    let mut trailing = snapshot.bytes().to_vec();
    trailing.push(0xa5);
    let trailing_error = FrameSequenceManifest::decode_snapshot(
        &trailing,
        snapshot_identity(&trailing),
        limits(reserved),
        remaining,
    )
    .expect_err("correctly pinned trailing data must refuse");
    assert_eq!(trailing_error, FrameSequenceError::TrailingBytes);

    let mut unsupported_version = snapshot.bytes().to_vec();
    let next_version = FRAME_SEQUENCE_MANIFEST_VERSION
        .checked_add(1)
        .expect("fixture schema version can increment");
    unsupported_version[8..10].copy_from_slice(&next_version.to_le_bytes());
    let version_error = FrameSequenceManifest::decode_snapshot(
        &unsupported_version,
        snapshot_identity(&unsupported_version),
        limits(reserved),
        remaining,
    )
    .expect_err("correctly pinned unsupported grammar must refuse");
    assert_eq!(
        version_error,
        FrameSequenceError::UnsupportedVersion {
            version: next_version
        }
    );

    let mut decoded = FrameSequenceManifest::decode_snapshot(
        snapshot.bytes(),
        snapshot.identity(),
        limits(reserved),
        remaining,
    )
    .expect("decode and re-admit resumable state");
    assert_eq!(decoded, resumed);
    register_indices(&mut decoded, &artifacts, 3..artifacts.len());

    let states = file_states(&decoded, &artifacts);
    let before_finalize_cancel = decoded.clone();
    let finalize_cancel = decoded
        .finalize_with(|| false, |_| None)
        .expect_err("finalization cancellation must not publish a seal");
    assert_eq!(finalize_cancel, FrameSequenceError::Cancelled);
    assert_eq!(decoded, before_finalize_cancel);
    let resumed_seal = decoded
        .finalize_with(|| true, |path| states.get(path).copied())
        .expect("finalize resumed sequence");

    let mut uninterrupted = manifest_for(context, &artifacts, true);
    register_indices(&mut uninterrupted, &artifacts, 0..artifacts.len());
    let uninterrupted_states = file_states(&uninterrupted, &artifacts);
    let uninterrupted_seal = uninterrupted
        .finalize_with(|| true, |path| uninterrupted_states.get(path).copied())
        .expect("finalize uninterrupted sequence");
    assert_eq!(resumed_seal.bytes(), uninterrupted_seal.bytes());
    assert_eq!(resumed_seal.identity(), uninterrupted_seal.identity());

    let decoded_final = FrameSequenceManifest::decode_snapshot(
        resumed_seal.bytes(),
        resumed_seal.identity(),
        limits(reserved),
        0,
    )
    .expect("finalized snapshots have no remaining storage reservation");
    assert_eq!(decoded_final.state(), FrameSequenceState::Finalized);
}

#[test]
fn g0_missing_stale_source_lineage_and_adaptive_partitions_refuse() {
    let adaptive = FrameSamplingStats::Adaptive {
        min_spp: 2,
        max_spp: 8,
        total_samples: 20,
        converged_pixels: 3,
        maximum_sample_pixels: 1,
    };
    descriptor(
        20,
        0,
        FrameArtifactRole::RawMaster,
        1.0,
        FrameArtifactFormat::OpenExr,
        2,
        2,
        vec![channel("R", FrameChannelType::Float32)],
        adaptive,
    );
    for invalid in [
        FrameSamplingStats::Adaptive {
            min_spp: 2,
            max_spp: 8,
            total_samples: 20,
            converged_pixels: 2,
            maximum_sample_pixels: 1,
        },
        FrameSamplingStats::Adaptive {
            min_spp: 2,
            max_spp: 8,
            total_samples: 13,
            converged_pixels: 3,
            maximum_sample_pixels: 1,
        },
        FrameSamplingStats::Adaptive {
            min_spp: 8,
            max_spp: 8,
            total_samples: 32,
            converged_pixels: 1,
            maximum_sample_pixels: 3,
        },
    ] {
        assert!(matches!(
            FrameArtifactDescriptor::try_new(
                20,
                0,
                FrameArtifactRole::RawMaster,
                1.0,
                FrameArtifactFormat::OpenExr,
                2,
                2,
                vec![channel("R", FrameChannelType::Float32)],
                invalid,
            ),
            Err(FrameSequenceError::InvalidSampling)
        ));
    }

    let context = context();
    let raw = raw_descriptor(21, 5, 1.25);
    let derived = descriptor(
        21,
        5,
        FrameArtifactRole::DenoisedIntermediate,
        1.25,
        FrameArtifactFormat::OpenExr,
        2,
        1,
        exr_channels(),
        FrameSamplingStats::Uniform { spp: 8 },
    );
    assert!(matches!(
        ExpectedFrameArtifact::try_new(
            derived.clone(),
            64,
            Some(FrameArtifactKey::new(21, 4, FrameArtifactRole::RawMaster)),
        ),
        Err(FrameSequenceError::InvalidSource { .. })
    ));
    let derived_expected =
        ExpectedFrameArtifact::try_new(derived.clone(), 64, Some(raw.key())).unwrap();
    let missing_source =
        FrameSequenceManifest::try_new(context, vec![derived_expected.clone()], limits(64), 64)
            .expect_err("derived row without its declared source row must refuse");
    assert!(matches!(
        missing_source,
        FrameSequenceError::InvalidSource { .. }
    ));

    let expected = vec![
        ExpectedFrameArtifact::try_new(raw.clone(), 64, None).unwrap(),
        derived_expected,
    ];
    let raw_bytes = b"lineage-raw";
    let derived_bytes = b"lineage-derived";
    let raw_hash = content_hash(raw_bytes);

    let mut placeholders =
        FrameSequenceManifest::try_new(context, expected.clone(), limits(128), 128).unwrap();
    let placeholder_raw_path = path_for(&placeholders, raw.key());
    let placeholder_derived_path = path_for(&placeholders, derived.key());
    let before_artifact_placeholder = placeholders.clone();
    let artifact_placeholder = placeholders
        .register_artifact(
            &placeholder_raw_path,
            &FrameArtifactObservation::new(
                raw.clone(),
                context.profile_id(),
                FrameArtifactFileState::new(ContentHash([0; 32]), 1),
                None,
            ),
        )
        .expect_err("zero artifact content identity must refuse");
    assert!(matches!(
        artifact_placeholder,
        FrameSequenceError::PlaceholderIdentity {
            field: "artifact content hash"
        }
    ));
    assert_eq!(placeholders, before_artifact_placeholder);
    let before_source_placeholder = placeholders.clone();
    let source_placeholder = placeholders
        .register_artifact(
            &placeholder_derived_path,
            &FrameArtifactObservation::new(
                derived.clone(),
                context.profile_id(),
                FrameArtifactFileState::from_bytes(derived_bytes).unwrap(),
                Some(ContentHash([0; 32])),
            ),
        )
        .expect_err("zero source content identity must refuse");
    assert!(matches!(
        source_placeholder,
        FrameSequenceError::PlaceholderIdentity {
            field: "source content hash"
        }
    ));
    assert_eq!(placeholders, before_source_placeholder);

    let matching_late_source_bytes = b"matching-late-source";
    let matching_late_source_hash = content_hash(matching_late_source_bytes);
    let mut derived_first =
        FrameSequenceManifest::try_new(context, expected.clone(), limits(128), 128).unwrap();
    let derived_first_raw_path = path_for(&derived_first, raw.key());
    let derived_first_path = path_for(&derived_first, derived.key());
    derived_first
        .register_artifact_bytes(
            &derived_first_path,
            derived.clone(),
            context.profile_id(),
            derived_bytes,
            Some(matching_late_source_hash),
        )
        .expect("a derived row may arrive before its declared source bytes");
    let before_mismatching_late_source = derived_first.clone();
    let mismatching_late_source = derived_first
        .register_artifact_bytes(
            &derived_first_raw_path,
            raw.clone(),
            context.profile_id(),
            b"mismatching-late-source",
            None,
        )
        .expect_err("later source bytes must agree with an existing dependent claim");
    assert!(matches!(
        mismatching_late_source,
        FrameSequenceError::SourceHashMismatch { .. }
    ));
    assert_eq!(derived_first, before_mismatching_late_source);
    assert_eq!(
        derived_first
            .register_artifact_bytes(
                &derived_first_raw_path,
                raw.clone(),
                context.profile_id(),
                matching_late_source_bytes,
                None,
            )
            .unwrap(),
        RegistrationOutcome::Recorded
    );
    assert_eq!(derived_first.completed_artifacts(), 2);

    let mut incomplete =
        FrameSequenceManifest::try_new(context, expected.clone(), limits(128), 128).unwrap();
    let missing = incomplete
        .finalize_with(|| true, |_| None)
        .expect_err("unregistered expected output must refuse finalization");
    assert!(matches!(
        missing,
        FrameSequenceError::MissingArtifact { .. }
    ));
    assert_eq!(incomplete.state(), FrameSequenceState::Incomplete);

    let raw_path = path_for(&incomplete, raw.key());
    let derived_path = path_for(&incomplete, derived.key());
    incomplete
        .register_artifact_bytes(
            &raw_path,
            raw.clone(),
            context.profile_id(),
            raw_bytes,
            None,
        )
        .unwrap();

    let before_missing_lineage = incomplete.clone();
    let missing_lineage = incomplete
        .register_artifact_bytes(
            &derived_path,
            derived.clone(),
            context.profile_id(),
            derived_bytes,
            None,
        )
        .expect_err("derived completion without source hash must refuse");
    assert!(matches!(
        missing_lineage,
        FrameSequenceError::DescriptorMismatch {
            field: "source content hash",
            ..
        }
    ));
    assert_eq!(incomplete, before_missing_lineage);

    let before_stale_source = incomplete.clone();
    let mut stale_source_hash_polls = 0_u32;
    let stale_source = incomplete
        .register_artifact_bytes_with_poll(
            &derived_path,
            derived.clone(),
            context.profile_id(),
            derived_bytes,
            Some(content_hash(b"stale-source")),
            || {
                stale_source_hash_polls += 1;
                true
            },
        )
        .expect_err("known stale source identity must refuse transactionally");
    assert!(matches!(
        stale_source,
        FrameSequenceError::SourceHashMismatch {
            expected,
            ..
        } if expected == raw_hash
    ));
    assert_eq!(incomplete, before_stale_source);
    assert_eq!(
        stale_source_hash_polls, 0,
        "known stale source hashes must refuse before hashing derived bytes"
    );

    let raw_state = FrameArtifactFileState::from_bytes(raw_bytes).unwrap();
    let derived_state = FrameArtifactFileState::from_bytes(derived_bytes).unwrap();

    let mut stale_file =
        FrameSequenceManifest::try_new(context, expected, limits(128), 128).unwrap();
    stale_file
        .register_artifact_bytes(&raw_path, raw, context.profile_id(), raw_bytes, None)
        .unwrap();
    stale_file
        .register_artifact_bytes(
            &derived_path,
            derived,
            context.profile_id(),
            derived_bytes,
            Some(raw_hash),
        )
        .unwrap();
    let changed_derived = FrameArtifactFileState::from_bytes(b"changed-on-disk").unwrap();
    let stale = stale_file
        .finalize_with(
            || true,
            |path| match path {
                path if path == raw_path => Some(raw_state),
                path if path == derived_path => Some(changed_derived),
                _ => None,
            },
        )
        .expect_err("fresh observation must catch a changed file");
    assert!(matches!(
        stale,
        FrameSequenceError::StaleArtifact {
            expected,
            actual,
            ..
        } if expected == derived_state && actual == changed_derived
    ));
    assert_eq!(stale_file.state(), FrameSequenceState::Incomplete);
}
