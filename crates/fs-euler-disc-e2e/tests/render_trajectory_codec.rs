//! G0/G3 integration checks for the canonical Euler render-trajectory codec.
//!
//! These tests exercise the public durable boundary rather than private wire
//! helpers: deterministic identity, exact binary64 retention, bounded refusal,
//! chunk-boundary behavior, corruption/truncation handling, and successful use
//! of a decoded artifact by the timeline and control consumers.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::{
    DeclaredDiscontinuityKind, DeclaredTimelineDiscontinuity, DerivedEulerQois,
    EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
    EulerControlStream, EulerRenderTrajectoryArtifact, EventEvaluationSide, RenderBaseFrame,
    RenderBaseModeState, RenderChannelAvailability, RenderContactBranch, RenderContactGeometry,
    RenderContactTransition, RenderMassProperties, RenderSampleDisposition, RenderSupportFeature,
    RenderTrajectory, RenderTrajectoryAuthority, RenderTrajectoryCodecBudget,
    RenderTrajectoryCodecError, RenderTrajectoryMetadata, RenderTrajectorySampleInput,
    RenderUnitSystem, RenderWorldFrame, TimelineEvent, TimelineResampler, TimelineSampleSource,
    coupled_runner::{ChannelOwnership, ChannelWrench, ContactTransitionKind},
    render_motion_bridge::EulerRenderMotionBridge,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const RAW_STEP_BITS: u64 = 0x3fc0_0000_0000_0001;
const WIRE_HEADER_LEN: usize = 116;

fn wire_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap())
}

fn wire_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

fn first_chunk_offset(bytes: &[u8]) -> usize {
    WIRE_HEADER_LEN + wire_u32(bytes, 36) as usize + wire_u64(bytes, 44) as usize
}

fn next_chunk_offset(bytes: &[u8], chunk_offset: usize) -> usize {
    chunk_offset + 56 + wire_u64(bytes, chunk_offset + 16) as usize
}

fn with_cx<R>(cancelled: bool, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x434f_4445_435f_5431,
                kernel_id: 0x4555_4c45,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn identity(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.test.render-trajectory-codec.v1",
        label.as_bytes(),
    )
}

fn mass() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap()
}

fn state(time_s: f64) -> RigidBodyState {
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.35).unwrap();
    RigidBodyState::new(
        Pose::new(Vec3::new(time_s, -0.5 * time_s, 1.0), orientation).unwrap(),
        Vec3::new(2.0, -1.0, 0.0),
        Vec3::new(0.25, 0.5, 1.0),
    )
    .unwrap()
}

fn wrench(force: Vec3, torque: Vec3, work_j: f64) -> ChannelWrench {
    ChannelWrench {
        force_world_n: force,
        torque_world_nm: torque,
        work_j,
    }
}

fn sample(
    interval_start_time_s: f64,
    time_s: f64,
    disposition: RenderSampleDisposition,
    retain_forcing: bool,
) -> RenderTrajectorySampleInput {
    let state = state(time_s);
    let orientation = state.pose().orientation();
    let duration_s = time_s - interval_start_time_s;
    let channels = if retain_forcing && duration_s > 0.0 {
        ChannelOwnership {
            gravity: wrench(Vec3::new(0.0, 0.0, -19.62), Vec3::ZERO, -0.05 * duration_s),
            gas: wrench(
                Vec3::new(-0.1, 0.05, 0.0),
                Vec3::new(0.0, 0.0, -0.01),
                -0.01 * duration_s,
            ),
            ..ChannelOwnership::default()
        }
    } else {
        ChannelOwnership::default()
    };
    RenderTrajectorySampleInput {
        interval_start_time_s,
        time_s,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m: state.pose().position_world(),
        orientation_body_to_world: orientation.components(),
        linear_momentum_world_kg_m_per_s: state.linear_momentum_world(),
        angular_momentum_body_kg_m2_per_s: state.angular_momentum_body(),
        symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        contact_branch: RenderContactBranch::Open,
        contact_geometry: None,
        signed_gap_m: 1.0e-3,
        interval_contact_active: false,
        interval_normal_force_n: 0.0,
        contact_transitions: Vec::new(),
        base_mode: Some(RenderBaseModeState {
            displacement_m: 0.25 * time_s,
            velocity_m_per_s: 0.25,
        }),
        channels,
        mechanical_energy_j: 10.0 - 0.06 * time_s,
        energy_defect_j: 0.0,
        qois: DerivedEulerQois::from_state(state, mass(), -0.0).unwrap(),
        disposition,
        terminal_event: None,
    }
}

fn metadata(first: &RenderTrajectorySampleInput, timestep_s: f64) -> RenderTrajectoryMetadata {
    let orientation = UnitQuaternion::new(
        first.orientation_body_to_world[0],
        first.orientation_body_to_world[1],
        first.orientation_body_to_world[2],
        first.orientation_body_to_world[3],
    )
    .unwrap();
    let initial_state = RigidBodyState::new(
        Pose::new(first.center_of_mass_world_m, orientation).unwrap(),
        first.linear_momentum_world_kg_m_per_s,
        first.angular_momentum_body_kg_m2_per_s,
    )
    .unwrap();
    RenderTrajectoryMetadata {
        schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        specimen_profile_identity: identity("profile"),
        specimen_chart_identity: identity("chart"),
        mass_properties: RenderMassProperties {
            identity: identity("mass"),
            properties: mass(),
        },
        initial_state,
        initial_base_mode: first.base_mode.unwrap(),
        base_model_identity: identity("base"),
        base_frame: RenderBaseFrame {
            origin_world_m: Vec3::ZERO,
            orientation_base_to_world: UnitQuaternion::IDENTITY,
        },
        model_identity: identity("model"),
        channel_availability: RenderChannelAvailability::ALL_AVAILABLE,
        configuration_identity: identity("configuration"),
        configuration_fingerprint: 0x434f_4445_435f_5631,
        timestep_s,
        producer_version: "render-trajectory-codec-test-v1".into(),
        applicability: "durable reduced-model visualization transport only".into(),
        no_claims: vec![
            "does not add mechanical resolution".into(),
            "does not synthesize a physical acoustic waveform".into(),
        ],
        authority: RenderTrajectoryAuthority::SimulationEvidence,
    }
}

fn trajectory(sample_count: usize, timestep_s: f64, retain_forcing: bool) -> RenderTrajectory {
    assert!(sample_count > 0);
    let mut inputs = Vec::with_capacity(sample_count);
    for index in 0..sample_count {
        let interval_start_time_s = if index == 0 {
            0.0
        } else {
            (index - 1) as f64 * timestep_s
        };
        let time_s = index as f64 * timestep_s;
        let disposition = if index + 1 == sample_count {
            RenderSampleDisposition::HorizonCensored
        } else {
            RenderSampleDisposition::Continue
        };
        inputs.push(sample(
            interval_start_time_s,
            time_s,
            disposition,
            retain_forcing,
        ));
    }
    let trajectory_metadata = metadata(&inputs[0], timestep_s);
    RenderTrajectory::try_new(trajectory_metadata, inputs).unwrap()
}

fn small_artifact(cx: &Cx<'_>) -> EulerRenderTrajectoryArtifact {
    let timestep_s = f64::from_bits(RAW_STEP_BITS);
    EulerRenderTrajectoryArtifact::try_from_trajectory(
        identity("campaign"),
        trajectory(3, timestep_s, true),
        vec![DeclaredTimelineDiscontinuity {
            time_s: timestep_s,
            kind: DeclaredDiscontinuityKind::ContinuationSeam,
        }],
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .unwrap()
}

fn alternating_quaternion_trajectory() -> RenderTrajectory {
    // This canonical quaternion is a known binary64 two-cycle under repeated
    // normalization: new(A) = B and new(B) = A. Canonical replay must retain
    // A exactly instead of applying `new` again.
    let canonical_components = [
        f64::from_bits(0x3fe4_0d80_d1de_135d),
        f64::from_bits(0xbfcf_e426_d1e9_4e92),
        f64::from_bits(0x3fe7_8115_6776_d0ee),
        f64::from_bits(0x3fb3_6315_2552_fc97),
    ];
    let orientation = UnitQuaternion::from_canonical_components(canonical_components).unwrap();
    let renormalized = UnitQuaternion::new(
        canonical_components[0],
        canonical_components[1],
        canonical_components[2],
        canonical_components[3],
    )
    .unwrap();
    assert_ne!(
        orientation.components().map(f64::to_bits),
        renormalized.components().map(f64::to_bits)
    );

    let mut inputs = Vec::new();
    for index in 0usize..3 {
        let interval_start_time_s = index.saturating_sub(1) as f64;
        let time_s = index as f64;
        let mut input = sample(
            interval_start_time_s,
            time_s,
            if index == 2 {
                RenderSampleDisposition::HorizonCensored
            } else {
                RenderSampleDisposition::Continue
            },
            false,
        );
        let state = RigidBodyState::new(
            Pose::new(input.center_of_mass_world_m, orientation).unwrap(),
            input.linear_momentum_world_kg_m_per_s,
            input.angular_momentum_body_kg_m2_per_s,
        )
        .unwrap();
        input.orientation_body_to_world = canonical_components;
        input.symmetry_axis_world = orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        input.qois = DerivedEulerQois::from_state(state, mass(), -0.0).unwrap();
        inputs.push(input);
    }
    let trajectory_metadata = metadata(&inputs[0], 1.0);
    RenderTrajectory::try_new(trajectory_metadata, inputs).unwrap()
}

fn transition_trajectory() -> RenderTrajectory {
    let first = sample(0.0, 0.0, RenderSampleDisposition::Continue, false);
    let mut second = sample(0.0, 1.0, RenderSampleDisposition::HorizonCensored, false);
    second.contact_branch = RenderContactBranch::Closed;
    second.contact_geometry = Some(RenderContactGeometry {
        point_world_m: Vec3::new(1.0, -0.5, 0.0),
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        support_feature: RenderSupportFeature::ProfileFeature(7),
    });
    second.signed_gap_m = -1.0e-8;
    second.interval_contact_active = true;
    second.contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Reimpact,
        time_s: 0.5,
        bracket_start_s: 0.49,
        bracket_end_s: 0.51,
    }];
    let trajectory_metadata = metadata(&first, 1.0);
    RenderTrajectory::try_new(trajectory_metadata, vec![first, second]).unwrap()
}

fn metadata_text_bytes(trajectory: &RenderTrajectory) -> usize {
    let metadata = trajectory.metadata();
    metadata.producer_version.len()
        + metadata.applicability.len()
        + metadata.no_claims.iter().map(String::len).sum::<usize>()
}

#[test]
fn canonical_roundtrip_is_stable_and_preserves_receipt_and_raw_f64_bits() {
    with_cx(false, |cx| {
        let artifact = small_artifact(cx);
        let bytes = artifact
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap();
        let receipt = artifact.receipt();

        assert_eq!(receipt.source_campaign_identity(), identity("campaign"));
        assert_eq!(receipt.byte_len(), bytes.len() as u64);
        assert_eq!(receipt.sample_count(), 3);
        assert_eq!(receipt.transition_count(), 0);
        assert_eq!(receipt.chunk_count(), 1);
        assert_eq!(
            receipt.artifact_identity().to_hex(),
            "3cc387401adf46b030ace241e112061a0f642ea9018f718320e398640f12c194"
        );

        let decoded = EulerRenderTrajectoryArtifact::from_canonical_bytes_verified(
            &bytes,
            receipt.artifact_identity(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(decoded.receipt(), receipt);
        assert_eq!(decoded.trajectory(), artifact.trajectory());
        assert_eq!(
            decoded.trajectory().samples()[1].input().time_s.to_bits(),
            RAW_STEP_BITS
        );
        assert_eq!(
            decoded.trajectory().samples()[1]
                .input()
                .qois
                .precession_acceleration_rad_per_s2
                .to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(
            decoded
                .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
                .unwrap(),
            bytes
        );

        let replay = small_artifact(cx);
        assert_eq!(replay.receipt(), receipt);
        assert_eq!(
            replay
                .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
                .unwrap(),
            bytes
        );
    });
}

#[test]
fn canonical_replay_does_not_renormalize_a_two_cycle_quaternion() {
    with_cx(false, |cx| {
        let trajectory = alternating_quaternion_trajectory();
        let admitted_bits = trajectory.samples()[0]
            .input()
            .orientation_body_to_world
            .map(f64::to_bits);
        let renormalized = UnitQuaternion::new(
            f64::from_bits(admitted_bits[0]),
            f64::from_bits(admitted_bits[1]),
            f64::from_bits(admitted_bits[2]),
            f64::from_bits(admitted_bits[3]),
        )
        .unwrap();
        assert_ne!(admitted_bits, renormalized.components().map(f64::to_bits));

        let artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
            identity("alternating-quaternion-campaign"),
            trajectory,
            Vec::new(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .unwrap();
        let bytes = artifact
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap();
        let decoded = EulerRenderTrajectoryArtifact::from_canonical_bytes(
            &bytes,
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(
            decoded.trajectory().samples()[0]
                .input()
                .orientation_body_to_world
                .map(f64::to_bits),
            admitted_bits
        );
        assert_eq!(
            decoded
                .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
                .unwrap(),
            bytes
        );
    });
}

#[test]
fn every_truncated_prefix_trailing_data_and_corruption_are_rejected() {
    with_cx(false, |cx| {
        let artifact = small_artifact(cx);
        let bytes = artifact
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap();

        for prefix_len in 0..bytes.len() {
            assert!(
                EulerRenderTrajectoryArtifact::from_canonical_bytes(
                    &bytes[..prefix_len],
                    RenderTrajectoryCodecBudget::DEFAULT,
                    cx,
                )
                .is_err(),
                "accepted truncated prefix of {prefix_len}/{} bytes",
                bytes.len()
            );
        }

        let mut with_trailing_byte = bytes.clone();
        with_trailing_byte.push(0xa5);
        assert!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &with_trailing_byte,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            )
            .is_err()
        );

        let mut corrupted_chunk = bytes.clone();
        let payload_byte = corrupted_chunk.len() / 2;
        corrupted_chunk[payload_byte] ^= 0x01;
        assert!(matches!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &corrupted_chunk,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::ChunkFingerprintMismatch(0))
        ));

        let mut corrupted_trailer = bytes;
        let last = corrupted_trailer.len() - 1;
        corrupted_trailer[last] ^= 0x01;
        assert!(matches!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &corrupted_trailer,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::PayloadFingerprintMismatch)
        ));
    });
}

#[test]
fn deterministic_junk_corpus_never_crosses_the_transport_boundary() {
    with_cx(false, |cx| {
        let mut generator = 0x8c3c_010c_b475_4c91_u64;
        for case in 0..64usize {
            let length = case * 19;
            let mut bytes = Vec::with_capacity(length);
            for _ in 0..length {
                generator ^= generator << 13;
                generator ^= generator >> 7;
                generator ^= generator << 17;
                bytes.push(generator as u8);
            }
            if bytes.len() >= 8 && case % 2 == 0 {
                bytes[..8].copy_from_slice(b"FSEULTRJ");
            }
            assert!(
                EulerRenderTrajectoryArtifact::from_canonical_bytes(
                    &bytes,
                    RenderTrajectoryCodecBudget::DEFAULT,
                    cx,
                )
                .is_err(),
                "accepted deterministic junk case {case} with {} bytes",
                bytes.len()
            );
        }
    });
}

#[test]
fn versions_metadata_events_samples_and_chunk_topology_fail_closed() {
    with_cx(false, |cx| {
        let small = small_artifact(cx)
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap();

        let mut unsupported_version = small.clone();
        unsupported_version[8..10].copy_from_slice(&2u16.to_le_bytes());
        assert_eq!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &unsupported_version,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::UnsupportedCodecVersion(2))
        );

        let mut metadata_bit_flip = small.clone();
        metadata_bit_flip[WIRE_HEADER_LEN + 6] ^= 0x01;
        assert!(matches!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &metadata_bit_flip,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::PayloadFingerprintMismatch)
        ));

        let metadata_len = wire_u32(&small, 36) as usize;
        let mut declared_event_tag = small.clone();
        declared_event_tag[WIRE_HEADER_LEN + metadata_len + 8] = 3;
        assert!(matches!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &declared_event_tag,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::InvalidTag {
                field: "declared_discontinuity.kind",
                tag: 3,
            })
        ));

        let chunk = first_chunk_offset(&small);
        let mut sample_bit_flip = small;
        sample_bit_flip[chunk + 56 + 4 + 16] ^= 0x01;
        assert_eq!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &sample_bit_flip,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::ChunkFingerprintMismatch(0))
        );

        let source = trajectory(EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK + 1, 1.0, false);
        let two_chunks = EulerRenderTrajectoryArtifact::try_from_trajectory(
            identity("chunk-topology-campaign"),
            source,
            Vec::new(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .unwrap()
        .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
        .unwrap();
        let first = first_chunk_offset(&two_chunks);
        let second = next_chunk_offset(&two_chunks, first);
        let trailer = two_chunks.len() - 32;

        let mut duplicate_index = two_chunks.clone();
        duplicate_index[second..second + 4].copy_from_slice(&0u32.to_le_bytes());
        assert!(matches!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &duplicate_index,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::InvalidChunk { field: "index", .. })
        ));

        let mut reordered = Vec::with_capacity(two_chunks.len());
        reordered.extend_from_slice(&two_chunks[..first]);
        reordered.extend_from_slice(&two_chunks[second..trailer]);
        reordered.extend_from_slice(&two_chunks[first..second]);
        reordered.extend_from_slice(&two_chunks[trailer..]);
        assert!(matches!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &reordered,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::InvalidChunk { field: "index", .. })
        ));

        let mut missing_final_chunk = Vec::with_capacity(second + 32);
        missing_final_chunk.extend_from_slice(&two_chunks[..second]);
        missing_final_chunk.extend_from_slice(&two_chunks[trailer..]);
        assert!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &missing_final_chunk,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            )
            .is_err()
        );
    });
}

#[test]
fn out_of_band_identity_mismatch_is_explicit_and_fail_closed() {
    with_cx(false, |cx| {
        let artifact = small_artifact(cx);
        let bytes = artifact
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap();
        let expected = identity("different-artifact-root");

        assert_eq!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes_verified(
                &bytes,
                expected,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::ArtifactIdentityMismatch {
                expected,
                actual: artifact.receipt().artifact_identity(),
            })
        );
    });
}

#[test]
fn caller_budgets_are_enforced_at_encode_and_decode_boundaries() {
    with_cx(false, |cx| {
        let artifact = small_artifact(cx);
        let default_bytes = artifact
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap();
        let exact = RenderTrajectoryCodecBudget {
            max_artifact_bytes: artifact.receipt().byte_len(),
            max_samples: artifact.trajectory().samples().len(),
            max_total_transitions: 0,
            max_total_text_bytes: metadata_text_bytes(artifact.trajectory()),
        };
        let exact_bytes = artifact.canonical_bytes(exact, cx).unwrap();
        assert_eq!(exact_bytes, default_bytes);
        assert!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(&exact_bytes, exact, cx).is_ok()
        );

        let byte_short = RenderTrajectoryCodecBudget {
            max_artifact_bytes: exact.max_artifact_bytes - 1,
            ..exact
        };
        assert!(matches!(
            artifact.canonical_bytes(byte_short, cx),
            Err(RenderTrajectoryCodecError::ArtifactTooLarge { .. })
        ));
        assert!(matches!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(&default_bytes, byte_short, cx),
            Err(RenderTrajectoryCodecError::ArtifactTooLarge { .. })
        ));

        let sample_short = RenderTrajectoryCodecBudget {
            max_samples: exact.max_samples - 1,
            ..exact
        };
        assert!(matches!(
            artifact.canonical_bytes(sample_short, cx),
            Err(RenderTrajectoryCodecError::InvalidLength {
                field: "sample_count",
                ..
            })
        ));
        assert!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(&default_bytes, sample_short, cx,)
                .is_err()
        );

        let text_short = RenderTrajectoryCodecBudget {
            max_total_text_bytes: exact.max_total_text_bytes - 1,
            ..exact
        };
        assert!(matches!(
            artifact.canonical_bytes(text_short, cx),
            Err(RenderTrajectoryCodecError::InvalidLength {
                field: "metadata text bytes",
                ..
            })
        ));
        assert!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(&default_bytes, text_short, cx,)
                .is_err()
        );

        assert!(matches!(
            artifact.canonical_bytes(
                RenderTrajectoryCodecBudget {
                    max_artifact_bytes: 0,
                    ..RenderTrajectoryCodecBudget::DEFAULT
                },
                cx,
            ),
            Err(RenderTrajectoryCodecError::InvalidBudget(
                "max_artifact_bytes"
            ))
        ));
    });

    with_cx(false, |cx| {
        let eventful = EulerRenderTrajectoryArtifact::try_from_trajectory(
            identity("transition-budget-campaign"),
            transition_trajectory(),
            Vec::new(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(eventful.receipt().transition_count(), 1);
        let bytes = eventful
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap();
        let zero_transition_budget = RenderTrajectoryCodecBudget {
            max_total_transitions: 0,
            ..RenderTrajectoryCodecBudget::DEFAULT
        };
        assert!(matches!(
            eventful.canonical_bytes(zero_transition_budget, cx),
            Err(RenderTrajectoryCodecError::InvalidLength {
                field: "transition_count",
                ..
            })
        ));
        assert!(matches!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(&bytes, zero_transition_budget, cx,),
            Err(RenderTrajectoryCodecError::InvalidLength {
                field: "header.transition_count",
                ..
            })
        ));
    });
}

#[test]
fn chunk_boundary_roundtrip_uses_a_full_chunk_and_one_sample_tail() {
    with_cx(false, |cx| {
        let sample_count = EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK + 1;
        let source = trajectory(sample_count, 1.0, false);
        let artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
            identity("chunk-boundary-campaign"),
            source,
            Vec::new(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(artifact.receipt().sample_count(), sample_count as u32);
        assert_eq!(artifact.receipt().chunk_count(), 2);

        let bytes = artifact
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap();
        let decoded = EulerRenderTrajectoryArtifact::from_canonical_bytes(
            &bytes,
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .unwrap();
        assert_eq!(decoded.receipt(), artifact.receipt());
        assert_eq!(decoded.trajectory(), artifact.trajectory());
        assert_eq!(
            decoded.trajectory().samples()[EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK]
                .input()
                .time_s,
            EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK as f64
        );
    });
}

#[test]
fn decoded_artifact_drives_declared_seam_resampling_and_audio_controls() {
    with_cx(false, |cx| {
        let artifact = small_artifact(cx);
        let bytes = artifact
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap();
        let decoded = EulerRenderTrajectoryArtifact::from_canonical_bytes(
            &bytes,
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .unwrap();

        let resampler = TimelineResampler::with_declared_discontinuities(
            decoded.trajectory(),
            decoded.declared_discontinuities().to_vec(),
        )
        .unwrap();
        let timestep_s = f64::from_bits(RAW_STEP_BITS);
        let samples = resampler
            .resample(&[timestep_s], EventEvaluationSide::RightLimit)
            .unwrap();
        assert_eq!(
            samples[0].source,
            TimelineSampleSource::ExactSample { index: 1 }
        );
        assert!(samples[0].events_at_query.iter().any(|event| matches!(
            event,
            TimelineEvent::Declared(DeclaredTimelineDiscontinuity {
                time_s,
                kind: DeclaredDiscontinuityKind::ContinuationSeam,
            }) if time_s.to_bits() == RAW_STEP_BITS
        )));

        let bridge = EulerRenderMotionBridge::with_declared_discontinuities(
            decoded.trajectory(),
            decoded.declared_discontinuities().to_vec(),
        )
        .unwrap();
        let render_sample = bridge
            .sample_at_time(timestep_s, EventEvaluationSide::RightLimit)
            .unwrap();
        assert_eq!(
            render_sample.transform().translation_m(),
            [timestep_s, -0.5 * timestep_s, 1.0]
        );

        let controls = EulerControlStream::try_derive(decoded.trajectory(), cx).unwrap();
        assert!(controls.is_bound_to(decoded.trajectory()));
        assert_eq!(controls.visualization().len(), 3);
        assert_eq!(controls.audio().len(), 2);
        assert_eq!(
            controls.audio()[0]
                .channels
                .gas
                .available()
                .unwrap()
                .signed_work_j
                .to_bits(),
            (-0.01 * timestep_s).to_bits()
        );
        let horizon = controls.audio_visual_horizon().unwrap();
        assert_eq!(horizon.start_time_s.to_bits(), 0.0_f64.to_bits());
        assert_eq!(horizon.end_time_s.to_bits(), (2.0 * timestep_s).to_bits());
    });
}

#[test]
fn zero_campaign_invalid_seams_and_precancellation_are_refused() {
    let timestep_s = f64::from_bits(RAW_STEP_BITS);
    with_cx(false, |cx| {
        let source = trajectory(3, timestep_s, false);
        assert!(matches!(
            EulerRenderTrajectoryArtifact::try_from_trajectory(
                ContentHash([0; 32]),
                source.clone(),
                Vec::new(),
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::ZeroSourceCampaignIdentity)
        ));

        let invalid_between_samples = vec![DeclaredTimelineDiscontinuity {
            time_s: 0.5 * timestep_s,
            kind: DeclaredDiscontinuityKind::ProducerDeclared,
        }];
        assert!(matches!(
            EulerRenderTrajectoryArtifact::try_from_trajectory(
                identity("campaign"),
                source.clone(),
                invalid_between_samples,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::InvalidDeclaredDiscontinuity(0))
        ));

        let duplicate_seams = vec![
            DeclaredTimelineDiscontinuity {
                time_s: timestep_s,
                kind: DeclaredDiscontinuityKind::ContinuationSeam,
            },
            DeclaredTimelineDiscontinuity {
                time_s: timestep_s,
                kind: DeclaredDiscontinuityKind::ProducerDeclared,
            },
        ];
        assert!(matches!(
            EulerRenderTrajectoryArtifact::try_from_trajectory(
                identity("campaign"),
                source,
                duplicate_seams,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::InvalidDeclaredDiscontinuity(1))
        ));
    });

    with_cx(true, |cx| {
        assert!(matches!(
            EulerRenderTrajectoryArtifact::try_from_trajectory(
                identity("campaign"),
                trajectory(3, timestep_s, false),
                Vec::new(),
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::Cancelled)
        ));
    });

    let bytes = with_cx(false, |cx| {
        small_artifact(cx)
            .canonical_bytes(RenderTrajectoryCodecBudget::DEFAULT, cx)
            .unwrap()
    });
    with_cx(true, |cx| {
        assert_eq!(
            EulerRenderTrajectoryArtifact::from_canonical_bytes(
                &bytes,
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            ),
            Err(RenderTrajectoryCodecError::Cancelled)
        );
    });
}
