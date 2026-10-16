//! G0/G3/G4/G5 integration checks for measure-first 48 kHz audio reconstruction.

use core::{f64::consts::TAU, num::NonZeroUsize};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::{
    AUDIO_EXCITATION_ALGORITHM_VERSION, AUDIO_RECONSTRUCTION_FILTER_VERSION,
    AUDIO_RESAMPLING_ALGORITHM_VERSION, ArtisticTextureConfig, AudioEventFractionalDelay,
    AudioExcitationBudget, AudioExcitationMapper, AudioExcitationModelInput,
    AudioExcitationReduction, AudioReconstructionFilterSpec, AudioResampler,
    AudioResamplingBoundaryPolicy, AudioResamplingBudget, AudioResamplingError,
    AudioResamplingModelInput, ContactEventMeasure, ContactModeShape, ContactParticipationPolicy,
    DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerControlStream,
    EulerRenderTrajectoryArtifact, ModalComponentValues, ModalSpatialParticipation,
    ModalSynthesisBudget, ModalSynthesisError, ModalSynthesisModel, ModalSynthesisModelInput,
    ModeContactParticipationRule, RenderBaseFrame, RenderBaseModeState, RenderChannelAvailability,
    RenderContactBranch, RenderContactGeometry, RenderContactTransition, RenderMassProperties,
    RenderSampleDisposition, RenderSupportFeature, RenderTrajectory, RenderTrajectoryAuthority,
    RenderTrajectoryCodecBudget, RenderTrajectoryMetadata, RenderTrajectorySampleInput,
    RenderUnitSystem, RenderWorldFrame,
    coupled_runner::{ChannelOwnership, ChannelWrench, ContactTransitionKind},
};
use fs_evidence::{
    cinematic::{CinematicClock, CinematicClockDomain, SoundAuthority},
    cinematic_config::{CinematicComponentRef, CinematicComponentRole},
    cinematic_sound::{
        ListenerFrame, ListenerPose, SOUND_MASTER_SAMPLE_RATE_HZ, SOUND_SYNTHESIS_SCHEMA_VERSION,
        SoundAmplitudeReference, SoundChannelLayout, SoundExcitationChannel,
        SoundExcitationControl, SoundModalComponent, SoundMode, SoundModeParticipation,
        SoundModelAssumption, SoundRoomResponse, SoundSynthesisConfig, SoundSynthesisInput,
        SoundTerminalPolicy, SoundTrajectoryDisposition,
    },
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const AUDIO_START_TICK: i64 = 6_000;
const AUDIO_END_TICK: i64 = 8_000;
const VIDEO_START_TICK: i64 = 3;
const VIDEO_END_TICK: i64 = 4;
const SOURCE_INTERVAL_COUNT: usize = 40;
const FILTER_HALF_LENGTH: u16 = 1_024;

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
                seed: 0x4155_4449_4f52_534d,
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
    hash_domain("org.frankensim.test.audio-resampling.v1", label.as_bytes())
}

fn mass() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap()
}

fn state_from(input: &RenderTrajectorySampleInput) -> RigidBodyState {
    let [w, x, y, z] = input.orientation_body_to_world;
    RigidBodyState::new(
        Pose::new(
            input.center_of_mass_world_m,
            UnitQuaternion::new(w, x, y, z).unwrap(),
        )
        .unwrap(),
        input.linear_momentum_world_kg_m_per_s,
        input.angular_momentum_body_kg_m2_per_s,
    )
    .unwrap()
}

fn sample(
    start_time_s: f64,
    end_time_s: f64,
    branch: RenderContactBranch,
    disposition: RenderSampleDisposition,
) -> RenderTrajectorySampleInput {
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.2).unwrap();
    let center_of_mass_world_m = Vec3::new(0.0, 0.0, 1.0);
    let point_body_m = Vec3::new(-0.5, 0.866_025_403_784_438_6, -1.0);
    let rotated_point_m = orientation.rotate_body_to_world(point_body_m);
    let point_world_m = Vec3::new(
        center_of_mass_world_m.x + rotated_point_m.x,
        center_of_mass_world_m.y + rotated_point_m.y,
        center_of_mass_world_m.z + rotated_point_m.z,
    );
    let positive_duration = end_time_s > start_time_s;
    let contact_active = positive_duration && branch == RenderContactBranch::Closed;
    let mut input = RenderTrajectorySampleInput {
        interval_start_time_s: start_time_s,
        time_s: end_time_s,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m,
        orientation_body_to_world: orientation.components(),
        linear_momentum_world_kg_m_per_s: Vec3::ZERO,
        angular_momentum_body_kg_m2_per_s: Vec3::ZERO,
        symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        contact_branch: branch,
        contact_geometry: (branch == RenderContactBranch::Closed).then_some(
            RenderContactGeometry {
                point_world_m,
                normal_world: Vec3::new(0.0, 0.0, 1.0),
                support_feature: RenderSupportFeature::ProfileFeature(1),
            },
        ),
        signed_gap_m: if branch == RenderContactBranch::Closed {
            0.0
        } else {
            1.0e-3
        },
        interval_contact_active: contact_active,
        interval_normal_force_n: if contact_active { 2.0 } else { 0.0 },
        contact_transitions: Vec::new(),
        base_mode: Some(RenderBaseModeState {
            displacement_m: 0.0,
            velocity_m_per_s: 0.0,
        }),
        channels: ChannelOwnership::default(),
        mechanical_energy_j: 8.0,
        energy_defect_j: 0.0,
        qois: DerivedEulerQois {
            inclination_rad: 0.0,
            precession_rad_per_s: 0.0,
            spin_rad_per_s: 0.0,
            precession_acceleration_rad_per_s2: 0.0,
        },
        disposition,
        terminal_event: None,
    };
    let state = state_from(&input);
    input.qois = DerivedEulerQois::from_state(state, mass(), 0.0).unwrap();
    input
}

fn wrench(force_world_n: Vec3, work_j: f64) -> ChannelWrench {
    ChannelWrench {
        force_world_n,
        torque_world_nm: Vec3::ZERO,
        work_j,
    }
}

fn reference_interval_frame_counts() -> Vec<i64> {
    (0..SOURCE_INTERVAL_COUNT)
        .map(|index| {
            if index < SOURCE_INTERVAL_COUNT / 2 {
                40
            } else {
                60
            }
        })
        .collect()
}

fn source_inputs() -> Vec<RenderTrajectorySampleInput> {
    let frame_counts = reference_interval_frame_counts();
    source_inputs_with(AUDIO_START_TICK, &frame_counts, Some(0.5), false, |_, _| {
        3.0
    })
}

fn source_inputs_with(
    source_start_tick: i64,
    interval_frame_counts: &[i64],
    first_reimpact_offset_frames: Option<f64>,
    terminal_opening: bool,
    mut gas_work_rate_w: impl FnMut(i64, i64) -> f64,
) -> Vec<RenderTrajectorySampleInput> {
    assert!(!interval_frame_counts.is_empty());
    assert!(
        interval_frame_counts
            .iter()
            .all(|frame_count| *frame_count > 0)
    );
    let source_start_s = source_start_tick as f64 / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    let initial_branch = if first_reimpact_offset_frames.is_some() {
        RenderContactBranch::Open
    } else {
        RenderContactBranch::Closed
    };
    let mut inputs = vec![sample(
        source_start_s,
        source_start_s,
        initial_branch,
        RenderSampleDisposition::Continue,
    )];
    let mut start_tick = source_start_tick;
    for (index, &frame_count) in interval_frame_counts.iter().enumerate() {
        let end_tick = start_tick + frame_count;
        let start_time_s = start_tick as f64 / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
        let end_time_s = end_tick as f64 / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
        let last = index + 1 == interval_frame_counts.len();
        let branch = if last && terminal_opening {
            RenderContactBranch::Open
        } else {
            RenderContactBranch::Closed
        };
        let mut retained = sample(
            start_time_s,
            end_time_s,
            branch,
            if last {
                RenderSampleDisposition::HorizonCensored
            } else {
                RenderSampleDisposition::Continue
            },
        );
        // A terminal opening still represents a contact-active interval; only
        // its retained endpoint is open.
        retained.interval_contact_active = true;
        retained.interval_normal_force_n = 2.0;
        retained.channels.contact = wrench(Vec3::new(0.0, 0.0, 2.0), 0.0);
        retained.channels.gas = wrench(
            Vec3::ZERO,
            gas_work_rate_w(start_tick, end_tick) * (end_time_s - start_time_s),
        );
        if index == 0
            && let Some(offset_frames) = first_reimpact_offset_frames
        {
            assert!(offset_frames >= 0.25 && offset_frames + 0.25 <= frame_count as f64);
            retained.contact_transitions.push(RenderContactTransition {
                kind: ContactTransitionKind::Reimpact,
                time_s: (source_start_tick as f64 + offset_frames)
                    / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ),
                bracket_start_s: (source_start_tick as f64 + offset_frames - 0.25)
                    / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ),
                bracket_end_s: (source_start_tick as f64 + offset_frames + 0.25)
                    / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ),
            });
        }
        if last && terminal_opening {
            retained.contact_transitions.push(RenderContactTransition {
                kind: ContactTransitionKind::Opening,
                time_s: end_time_s,
                bracket_start_s: (end_tick as f64 - 0.25) / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ),
                bracket_end_s: end_time_s,
            });
        }
        inputs.push(retained);
        start_tick = end_tick;
    }
    inputs
}

fn trajectory(inputs: Vec<RenderTrajectorySampleInput>) -> RenderTrajectory {
    let timestep_s = inputs
        .iter()
        .map(|input| input.time_s - input.interval_start_time_s)
        .fold(f64::MIN_POSITIVE, f64::max);
    let first = &inputs[0];
    let metadata = RenderTrajectoryMetadata {
        schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        specimen_profile_identity: identity("profile"),
        specimen_chart_identity: identity("chart"),
        mass_properties: RenderMassProperties {
            identity: identity("mass"),
            properties: mass(),
        },
        initial_state: state_from(first),
        initial_base_mode: first.base_mode.unwrap(),
        base_model_identity: identity("base"),
        base_frame: RenderBaseFrame {
            origin_world_m: Vec3::ZERO,
            orientation_base_to_world: UnitQuaternion::IDENTITY,
        },
        model_identity: identity("mechanics-model"),
        channel_availability: RenderChannelAvailability::ALL_AVAILABLE,
        configuration_identity: identity("configuration"),
        configuration_fingerprint: 0x4155_4449_4f52_5301,
        timestep_s,
        producer_version: "audio-resampling-test-v1".into(),
        applicability: "mapper-to-resampler integration checks".into(),
        no_claims: vec![
            "source cadence is not an audio sample rate".into(),
            "timing-only events carry no physical impulse".into(),
        ],
        authority: RenderTrajectoryAuthority::SimulationEvidence,
    };
    RenderTrajectory::try_new(metadata, inputs).unwrap()
}

fn mode() -> SoundMode {
    mode_with_frequency(200.0)
}

fn mode_with_frequency(frequency_hz: f64) -> SoundMode {
    SoundMode {
        mode_id: 1,
        component: SoundModalComponent::Disc,
        frequency_hz,
        damping_ratio: 0.02,
        modal_mass_kg: 0.2,
        source_participation: SoundModeParticipation {
            disc: 1.0,
            glass_plate: 0.0,
            base_assembly: 0.0,
        },
        radiation_gain_fs_s_per_m: 0.1,
        material_identity: identity("material"),
        base_identity: identity("modal-base"),
    }
}

fn mappings() -> Vec<SoundExcitationControl> {
    vec![
        SoundExcitationControl {
            channel: SoundExcitationChannel::ContactNormalForce,
            target_component: SoundModalComponent::Disc,
            source_scale: 1.0,
        },
        SoundExcitationControl {
            channel: SoundExcitationChannel::ExteriorGasBodySignedWorkRate,
            target_component: SoundModalComponent::Disc,
            source_scale: 1.0,
        },
    ]
}

fn video_clock() -> CinematicClock {
    CinematicClock::try_new(
        CinematicClockDomain::Video,
        24,
        1,
        VIDEO_START_TICK,
        VIDEO_END_TICK,
    )
    .unwrap()
}

fn audio_clock() -> CinematicClock {
    CinematicClock::try_new(
        CinematicClockDomain::Audio,
        SOUND_MASTER_SAMPLE_RATE_HZ,
        1,
        AUDIO_START_TICK,
        AUDIO_END_TICK,
    )
    .unwrap()
}

#[derive(Clone, Copy)]
struct FixtureOptions {
    video_clock: CinematicClock,
    audio_clock: CinematicClock,
    declared_source_bandwidth_hz: f64,
    filter: AudioReconstructionFilterSpec,
    budget: AudioResamplingBudget,
}

impl FixtureOptions {
    fn reference() -> Self {
        Self {
            video_clock: video_clock(),
            audio_clock: audio_clock(),
            declared_source_bandwidth_hz: 60.0,
            filter: AudioReconstructionFilterSpec {
                passband_edge_hz: 80.0,
                stopband_edge_hz: 350.0,
                half_length: FILTER_HALF_LENGTH,
                maximum_passband_ripple_db: 0.1,
                minimum_stopband_attenuation_db: 80.0,
                response_grid_intervals: 8_192,
            },
            budget: AudioResamplingBudget::reference_film(),
        }
    }
}

struct Fixture {
    resampler: AudioResampler,
    modal: ModalSynthesisModel,
    sound: SoundSynthesisConfig,
    expected_localized_force_n: f64,
    expected_distributed_force_n: f64,
    expected_projected_event_impulse_n_s: f64,
}

fn component_ref(
    role: CinematicComponentRole,
    component_identity: ContentHash,
    version: u32,
) -> CinematicComponentRef {
    CinematicComponentRef::try_new(role, component_identity, version).unwrap()
}

fn sound_configuration(
    trajectory_identity: ContentHash,
    mapper: &AudioExcitationMapper<'_, '_>,
    modal: &ModalSynthesisModel,
    resampler: &AudioResampler,
    controls: Vec<SoundExcitationControl>,
    video_clock: CinematicClock,
    audio_clock: CinematicClock,
) -> SoundSynthesisConfig {
    SoundSynthesisConfig::try_admit(SoundSynthesisInput {
        schema_version: SOUND_SYNTHESIS_SCHEMA_VERSION,
        authority: SoundAuthority::PhysicallyInformed,
        trajectory: component_ref(
            CinematicComponentRole::Trajectory,
            trajectory_identity,
            u32::from(EULER_RENDER_TRAJECTORY_SCHEMA_VERSION),
        ),
        excitation: component_ref(
            CinematicComponentRole::AudioExcitation,
            mapper.identity(),
            AUDIO_EXCITATION_ALGORITHM_VERSION,
        ),
        sound_model: component_ref(
            CinematicComponentRole::SoundModel,
            modal.identity(),
            fs_euler_disc_e2e::MODAL_SYNTHESIS_ALGORITHM_VERSION,
        ),
        microphone: component_ref(
            CinematicComponentRole::Microphone,
            identity("microphone"),
            1,
        ),
        room: component_ref(CinematicComponentRole::Room, identity("room"), 1),
        timeline: component_ref(CinematicComponentRole::Timeline, identity("timeline"), 1),
        video_clock,
        audio_clock,
        channel_layout: SoundChannelLayout::Stereo,
        listener: ListenerPose {
            frame: ListenerFrame::AnimatedCamera,
            position_m: [0.0, 0.0, 0.0],
            forward: [0.0, 0.0, -1.0],
            up: [0.0, 1.0, 0.0],
        },
        excitation_controls: controls,
        modes: modal.modes().to_vec(),
        room_response: SoundRoomResponse::Dry,
        amplitude_reference: SoundAmplitudeReference::DigitalFullScale { headroom_db: 6.0 },
        trajectory_disposition: SoundTrajectoryDisposition::HorizonCensored,
        terminal_policy: SoundTerminalPolicy::FadeAtLastAccepted {
            fade_sample_frames: 240,
        },
        resampler_identity: resampler.identity(),
        resampler_version: AUDIO_RESAMPLING_ALGORITHM_VERSION,
        filter_identity: resampler.filter_identity(),
        filter_version: AUDIO_RECONSTRUCTION_FILTER_VERSION,
        assumptions: vec![
            SoundModelAssumption::LinearModalSuperposition,
            SoundModelAssumption::TimeInvariantDamping,
            SoundModelAssumption::DeclaredExcitationCompleteness,
            SoundModelAssumption::DeclaredRoomResponse,
        ],
        calibration: None,
    })
    .unwrap()
}

fn test_artistic_texture() -> ArtisticTextureConfig {
    ArtisticTextureConfig {
        seed: 7,
        rolling_force_gain_n_per_w: 0.0,
        rolling_target_component: SoundModalComponent::Disc,
        band_low_hz: 100.0,
        band_high_hz: 1_000.0,
        reimpact_impulse_n_s: 0.25,
        reimpact_target_component: SoundModalComponent::Disc,
    }
}

fn modal_model(
    mode: SoundMode,
    maximum_total_sample_frames: u64,
    cx: &Cx<'_>,
) -> ModalSynthesisModel {
    ModalSynthesisModel::try_new(
        ModalSynthesisModelInput {
            sample_rate_hz: SOUND_MASTER_SAMPLE_RATE_HZ,
            modes: vec![mode],
            budget: ModalSynthesisBudget::reference_film(maximum_total_sample_frames),
        },
        cx,
    )
    .unwrap()
}

fn excitation_model_input(
    mappings: Vec<SoundExcitationControl>,
    artistic_texture: Option<ArtisticTextureConfig>,
    interval_count: usize,
) -> AudioExcitationModelInput {
    AudioExcitationModelInput {
        mappings,
        reduction: AudioExcitationReduction::RawIntervals,
        spatial_policy: ContactParticipationPolicy::ContactCoordinates {
            rules: vec![ModeContactParticipationRule {
                mode_id: 1,
                shape: ContactModeShape::AzimuthalCosine {
                    harmonic: 1,
                    phase_rad: 0.0,
                },
            }],
        },
        artistic_texture,
        budget: AudioExcitationBudget::reference_film(interval_count),
    }
}

fn resampling_model_input(options: FixtureOptions) -> AudioResamplingModelInput {
    AudioResamplingModelInput {
        video_clock: options.video_clock,
        audio_clock: options.audio_clock,
        declared_source_bandwidth_hz: options.declared_source_bandwidth_hz,
        filter: options.filter,
        boundary_policy: AudioResamplingBoundaryPolicy::HalfSampleEvenReflectionV1,
        event_fractional_delay: AudioEventFractionalDelay::LinearTwoBoundaryV1,
        budget: options.budget,
    }
}

fn try_fixture(options: FixtureOptions, cx: &Cx<'_>) -> Result<Fixture, AudioResamplingError> {
    try_fixture_with_inputs(
        options,
        source_inputs(),
        Some(test_artistic_texture()),
        true,
        cx,
    )
}

fn try_fixture_with_inputs(
    options: FixtureOptions,
    inputs: Vec<RenderTrajectorySampleInput>,
    artistic_texture: Option<ArtisticTextureConfig>,
    expect_constant_distributed_drive: bool,
    cx: &Cx<'_>,
) -> Result<Fixture, AudioResamplingError> {
    let interval_count = inputs.len().checked_sub(1).unwrap();
    let total_audio_frames =
        u64::try_from(options.audio_clock.end_tick_exclusive() - options.audio_clock.start_tick())
            .unwrap();
    let artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
        identity("campaign"),
        trajectory(inputs),
        Vec::new(),
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .unwrap();
    let trajectory_identity = artifact.receipt().artifact_identity();
    let controls = EulerControlStream::try_derive(artifact.trajectory(), cx).unwrap();
    let modal = modal_model(mode(), total_audio_frames, cx);
    let mappings = mappings();
    let mapper = AudioExcitationMapper::try_new(
        &artifact,
        &controls,
        &modal,
        excitation_model_input(mappings.clone(), artistic_texture, interval_count),
        cx,
    )
    .unwrap();
    let mapped = mapper
        .map_next_chunk(
            &mapper.initial_checkpoint(cx).unwrap(),
            NonZeroUsize::new(interval_count).unwrap(),
            cx,
        )
        .unwrap();
    assert_eq!(mapped.intervals.len(), interval_count);
    // Source sample zero is the retained zero-duration trajectory origin; the
    // first positive control interval therefore canonically keeps index one.
    assert_eq!(mapped.intervals[0].first_source_sample_index, 1);
    let representative = &mapped.intervals[0];
    let expected_spatial_factor = representative.spatial_envelopes[0].start_factor;
    let expected_localized_force_n =
        representative.localized_mean_generalized_force_n().disc * expected_spatial_factor;
    let expected_distributed_force_n = representative.distributed_mean_generalized_force_n().disc;
    let expected_projected_event_impulse_n_s = representative
        .events
        .first()
        .and_then(|event| {
            assert_eq!(event.measure, ContactEventMeasure::TimingOnly);
            assert_eq!(event.physical_impulse_n_s, ModalComponentValues::ZERO);
            event.artistic
        })
        .map_or(0.0, |artistic| {
            artistic.impulse_n_s.disc * expected_spatial_factor
        });
    if artistic_texture.is_some() {
        assert_ne!(expected_projected_event_impulse_n_s, 0.0);
    }
    for interval in &mapped.intervals {
        assert_close(
            interval.localized_mean_generalized_force_n().disc,
            representative.localized_mean_generalized_force_n().disc,
        );
        if expect_constant_distributed_drive {
            assert_close(
                interval.distributed_mean_generalized_force_n().disc,
                expected_distributed_force_n,
            );
        }
        assert_close(
            interval.spatial_envelopes[0].start_factor,
            expected_spatial_factor,
        );
        assert_close(
            interval.spatial_envelopes[0].end_factor,
            expected_spatial_factor,
        );
    }
    assert_close(expected_localized_force_n, -1.0);
    if expect_constant_distributed_drive {
        assert_close(expected_distributed_force_n, 3.0);
    }

    let resampler = AudioResampler::try_new(
        &mapper,
        &modal,
        mapped.intervals,
        resampling_model_input(options),
        cx,
    )?;
    let sound = sound_configuration(
        trajectory_identity,
        &mapper,
        &modal,
        &resampler,
        mappings,
        options.video_clock,
        options.audio_clock,
    );
    mapper.validate_sound_configuration(&sound).unwrap();
    modal.validate_sound_configuration(&sound).unwrap();
    resampler.validate_sound_configuration(&sound)?;
    Ok(Fixture {
        resampler,
        modal,
        sound,
        expected_localized_force_n,
        expected_distributed_force_n,
        expected_projected_event_impulse_n_s,
    })
}

fn try_mismatched_modal_admission(cx: &Cx<'_>) -> Result<AudioResampler, AudioResamplingError> {
    let inputs = source_inputs();
    let interval_count = inputs.len() - 1;
    let artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
        identity("modal-mismatch-campaign"),
        trajectory(inputs),
        Vec::new(),
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .unwrap();
    let controls = EulerControlStream::try_derive(artifact.trajectory(), cx).unwrap();
    let mapper_modal = modal_model(mode_with_frequency(200.0), 2_000, cx);
    let resampler_modal = modal_model(mode_with_frequency(201.0), 2_000, cx);
    let mapper = AudioExcitationMapper::try_new(
        &artifact,
        &controls,
        &mapper_modal,
        excitation_model_input(mappings(), None, interval_count),
        cx,
    )
    .unwrap();
    let mapped = mapper
        .map_next_chunk(
            &mapper.initial_checkpoint(cx).unwrap(),
            NonZeroUsize::new(interval_count).unwrap(),
            cx,
        )
        .unwrap();
    AudioResampler::try_new(
        &mapper,
        &resampler_modal,
        mapped.intervals,
        resampling_model_input(FixtureOptions::reference()),
        cx,
    )
}

fn assert_close(actual: f64, expected: f64) {
    let scale = 1.0_f64.max(actual.abs()).max(expected.abs());
    assert!(
        (actual - expected).abs() <= 5.0e-10 * scale,
        "actual={actual:.17e}, expected={expected:.17e}, delta={:.17e}",
        actual - expected
    );
}

fn tone_amplitude(samples: &[f64], frequency_hz: f64) -> f64 {
    let sample_rate_hz = f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    let (in_phase, quadrature) =
        samples
            .iter()
            .enumerate()
            .fold((0.0, 0.0), |(in_phase, quadrature), (index, sample)| {
                let phase = TAU * frequency_hz * index as f64 / sample_rate_hz;
                (
                    sample.mul_add(phase.cos(), in_phase),
                    sample.mul_add(phase.sin(), quadrature),
                )
            });
    2.0 * in_phase.hypot(quadrature) / samples.len() as f64
}

#[test]
fn g0_g3_mapper_to_resampler_preserves_dc_projection_events_and_exact_markers() {
    with_cx(false, |cx| {
        let fixture = try_fixture(FixtureOptions::reference(), cx).unwrap();
        let diagnostics = fixture.resampler.filter_diagnostics();
        assert_eq!(
            diagnostics.tap_count,
            usize::from(FILTER_HALF_LENGTH) * 2 + 1
        );
        assert!(diagnostics.measured_passband_ripple_db <= 0.1);
        assert!(diagnostics.measured_stopband_attenuation_db >= 80.0);
        assert_eq!(
            diagnostics.intrinsic_group_delay_frames,
            u32::from(FILTER_HALF_LENGTH)
        );
        assert_eq!(
            diagnostics.group_delay_compensation_frames,
            u32::from(FILTER_HALF_LENGTH)
        );
        assert_eq!(
            diagnostics.required_lookahead_frames,
            u32::from(FILTER_HALF_LENGTH)
        );
        assert_eq!(diagnostics.published_alignment_offset_frames, 0);
        let coefficients = fixture.resampler.filter_coefficients();
        assert_close(coefficients.iter().sum::<f64>(), 1.0);
        for pair in coefficients.iter().zip(coefficients.iter().rev()) {
            assert_eq!(pair.0.to_bits(), pair.1.to_bits());
        }

        assert_eq!(fixture.resampler.total_audio_frames(), 2_000);
        let alignment = fixture.resampler.alignment();
        assert_eq!(alignment.audio_frames_per_video_frame, 2_000);
        assert_eq!(alignment.endpoint_drift_audio_frames, 0);
        assert_eq!(alignment.markers.len(), 2);
        assert_eq!(alignment.markers[0].video_tick, VIDEO_START_TICK);
        assert_eq!(alignment.markers[0].audio_tick, AUDIO_START_TICK);
        assert_eq!(alignment.markers[0].audio_frame_offset, 0);
        assert_eq!(alignment.markers[1].video_tick, VIDEO_END_TICK);
        assert_eq!(alignment.markers[1].audio_tick, AUDIO_END_TICK);
        assert_eq!(alignment.markers[1].audio_frame_offset, 2_000);

        let chunk = fixture
            .resampler
            .resample_next_chunk(
                &fixture.sound,
                &fixture.resampler.initial_checkpoint(cx).unwrap(),
                NonZeroUsize::new(2_000).unwrap(),
                cx,
            )
            .unwrap();
        assert_eq!(chunk.drive_frames.len(), 2_000);
        assert_eq!(chunk.preparticipated_localized_force_n.len(), 2_000);
        assert_eq!(chunk.preparticipated_localized_impulse_n_s.len(), 2_000);
        for (frame, localized) in chunk
            .drive_frames
            .iter()
            .zip(&chunk.preparticipated_localized_force_n)
        {
            assert_eq!(
                frame.localized_generalized_force_n,
                ModalComponentValues::ZERO
            );
            assert_eq!(
                frame.localized_boundary_impulse_n_s,
                ModalComponentValues::ZERO
            );
            assert_eq!(
                frame.distributed_boundary_impulse_n_s,
                ModalComponentValues::ZERO
            );
            assert_close(
                frame.distributed_generalized_force_n.disc,
                fixture.expected_distributed_force_n,
            );
            assert_eq!(frame.distributed_generalized_force_n.glass_plate, 0.0);
            assert_eq!(frame.distributed_generalized_force_n.base_assembly, 0.0);
            assert_close(*localized, fixture.expected_localized_force_n);
        }
        let period_s = 1.0 / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
        assert_close(
            chunk
                .drive_frames
                .iter()
                .map(|frame| frame.distributed_generalized_force_n.disc)
                .sum::<f64>()
                * period_s,
            fixture.expected_distributed_force_n / 24.0,
        );
        assert_close(
            chunk.preparticipated_localized_force_n.iter().sum::<f64>() * period_s,
            fixture.expected_localized_force_n / 24.0,
        );

        assert_eq!(chunk.events.len(), 1);
        let event = &chunk.events[0];
        assert_eq!(event.source.measure, ContactEventMeasure::TimingOnly);
        assert_eq!(
            event.source.physical_impulse_n_s,
            ModalComponentValues::ZERO
        );
        assert_close(event.requested_sample_position, 0.5);
        assert_eq!(event.left_frame_offset, Some(0));
        assert_eq!(event.right_frame_offset, Some(1));
        assert_close(event.left_weight, 0.5);
        assert_close(event.right_weight, 0.5);
        assert_close(event.centroid_error_frames, 0.0);
        assert_close(event.bracket_start_sample_position, 0.25);
        assert_close(event.bracket_end_sample_position, 0.75);
        assert_close(
            chunk.preparticipated_localized_impulse_n_s[0],
            0.5 * fixture.expected_projected_event_impulse_n_s,
        );
        assert_close(
            chunk.preparticipated_localized_impulse_n_s[1],
            0.5 * fixture.expected_projected_event_impulse_n_s,
        );
        assert_close(
            chunk
                .preparticipated_localized_impulse_n_s
                .iter()
                .sum::<f64>(),
            fixture.expected_projected_event_impulse_n_s,
        );
        assert_eq!(
            chunk.sync_markers.as_slice(),
            fixture.resampler.alignment().markers.as_slice()
        );

        // Exercise the public safe handoff: callers must not be able to forget
        // that localized drive is already projected and filtered per mode.
        let modal_initial = fixture.modal.initial_checkpoint(cx).unwrap();
        let modal_output = chunk
            .synthesize_modal(&fixture.modal, &modal_initial, cx)
            .unwrap();
        let silently_discarded_localized_drive = fixture
            .modal
            .synthesize_chunk(
                &modal_initial,
                &chunk.drive_frames,
                ModalSpatialParticipation::Declared,
                cx,
            )
            .unwrap();
        assert_ne!(
            modal_output.mixed_samples_fs,
            silently_discarded_localized_drive.mixed_samples_fs
        );
        assert_eq!(modal_output.mixed_samples_fs.len(), 2_000);
        assert_eq!(modal_output.stem_frames.len(), 2_000);
        assert_eq!(modal_output.total_modal_energy_j.len(), 2_000);
        assert_eq!(modal_output.successor.next_sample_frame(), 2_000);
        assert!(
            modal_output
                .mixed_samples_fs
                .iter()
                .all(|sample| sample.is_finite())
        );
        assert!(
            modal_output
                .mixed_samples_fs
                .iter()
                .any(|sample| *sample != 0.0)
        );
    });
}

#[test]
fn g0_g3_integral_and_terminal_censored_events_obey_measure_and_ownership_rules() {
    with_cx(false, |cx| {
        let frame_counts = reference_interval_frame_counts();
        let inputs =
            source_inputs_with(AUDIO_START_TICK, &frame_counts, Some(1.0), true, |_, _| 3.0);
        let fixture = try_fixture_with_inputs(
            FixtureOptions::reference(),
            inputs,
            Some(test_artistic_texture()),
            true,
            cx,
        )
        .unwrap();
        let chunk = fixture
            .resampler
            .resample_next_chunk(
                &fixture.sound,
                &fixture.resampler.initial_checkpoint(cx).unwrap(),
                NonZeroUsize::new(2_000).unwrap(),
                cx,
            )
            .unwrap();

        assert_eq!(chunk.events.len(), 2);
        let integral = &chunk.events[0];
        assert_eq!(integral.source.kind, ContactTransitionKind::Reimpact);
        assert_eq!(integral.source.measure, ContactEventMeasure::TimingOnly);
        assert_eq!(
            integral.source.physical_impulse_n_s,
            ModalComponentValues::ZERO
        );
        assert_close(integral.requested_sample_position, 1.0);
        assert_eq!(integral.left_frame_offset, Some(1));
        assert_eq!(integral.right_frame_offset, None);
        assert_close(integral.left_weight, 1.0);
        assert_close(integral.right_weight, 0.0);
        assert_close(integral.centroid_error_frames, 0.0);
        assert_close(
            chunk.preparticipated_localized_impulse_n_s[1],
            fixture.expected_projected_event_impulse_n_s,
        );
        assert_close(
            chunk
                .preparticipated_localized_impulse_n_s
                .iter()
                .sum::<f64>(),
            fixture.expected_projected_event_impulse_n_s,
        );

        let terminal = &chunk.events[1];
        assert_eq!(terminal.source.kind, ContactTransitionKind::Opening);
        assert_eq!(terminal.source.measure, ContactEventMeasure::TimingOnly);
        assert_eq!(
            terminal.source.physical_impulse_n_s,
            ModalComponentValues::ZERO
        );
        assert!(terminal.source.artistic.is_none());
        assert_close(terminal.requested_sample_position, 2_000.0);
        assert_eq!(terminal.left_frame_offset, None);
        assert_eq!(terminal.right_frame_offset, None);
        assert_close(terminal.left_weight, 0.0);
        assert_close(terminal.right_weight, 0.0);
        assert_close(terminal.centroid_error_frames, 0.0);
        assert_close(terminal.bracket_start_sample_position, 1_999.75);
        assert_close(terminal.bracket_end_sample_position, 2_000.0);
        assert_eq!(chunk.diagnostics.owned_event_count, 2);
    });
}

#[test]
fn g3_g5_long_clock_markers_remain_exact_and_drift_free() {
    with_cx(false, |cx| {
        const VIDEO_FRAMES: i64 = 8 * 24;
        const AUDIO_FRAMES: i64 = 8 * SOUND_MASTER_SAMPLE_RATE_HZ as i64;
        let interval_frame_counts = vec![50; (AUDIO_FRAMES / 50) as usize];
        let inputs = source_inputs_with(0, &interval_frame_counts, None, false, |_, _| 3.0);
        let mut options = FixtureOptions::reference();
        options.video_clock =
            CinematicClock::try_new(CinematicClockDomain::Video, 24, 1, 0, VIDEO_FRAMES).unwrap();
        options.audio_clock = CinematicClock::try_new(
            CinematicClockDomain::Audio,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            1,
            0,
            AUDIO_FRAMES,
        )
        .unwrap();
        let fixture = try_fixture_with_inputs(options, inputs, None, true, cx).unwrap();
        let alignment = fixture.resampler.alignment();
        assert_eq!(fixture.resampler.total_audio_frames(), AUDIO_FRAMES as u64);
        assert_eq!(alignment.audio_frames_per_video_frame, 2_000);
        assert_eq!(alignment.endpoint_drift_audio_frames, 0);
        assert_eq!(alignment.markers.len(), VIDEO_FRAMES as usize + 1);
        for (frame, marker) in alignment.markers.iter().enumerate() {
            assert_eq!(marker.video_tick, frame as i64);
            assert_eq!(marker.audio_tick, frame as i64 * 2_000);
            assert_eq!(marker.audio_frame_offset, frame as u64 * 2_000);
        }
    });
}

#[test]
fn g1_g3_mapper_reconstruction_passes_source_tone_and_suppresses_sampling_image() {
    with_cx(false, |cx| {
        const AUDIO_FRAMES: usize = 12_000;
        const SOURCE_INTERVAL_FRAMES: usize = 10;
        const SPECTRAL_FILTER_HALF_LENGTH: u16 = 1_025;
        const INTERIOR_FRAMES: usize = AUDIO_FRAMES - 2 * SPECTRAL_FILTER_HALF_LENGTH as usize;
        let source_rate_hz = f64::from(SOUND_MASTER_SAMPLE_RATE_HZ) / SOURCE_INTERVAL_FRAMES as f64;
        // Ten exact DFT cycles over the FIR-halo-free interior. The source
        // samples repeat every 199 intervals, so its first ZOH image is also
        // an exact DFT bin and leakage cannot masquerade as attenuation.
        let source_frequency_hz =
            10.0 * f64::from(SOUND_MASTER_SAMPLE_RATE_HZ) / INTERIOR_FRAMES as f64;
        let image_frequency_hz = source_rate_hz - source_frequency_hz;
        let interval_frame_counts =
            vec![SOURCE_INTERVAL_FRAMES as i64; AUDIO_FRAMES / SOURCE_INTERVAL_FRAMES];
        let inputs = source_inputs_with(
            0,
            &interval_frame_counts,
            None,
            false,
            |start_tick, end_tick| {
                let midpoint_frame = 0.5 * (start_tick + end_tick) as f64;
                (TAU * source_frequency_hz * midpoint_frame
                    / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ))
                .sin()
            },
        );
        let mut options = FixtureOptions::reference();
        options.video_clock =
            CinematicClock::try_new(CinematicClockDomain::Video, 24, 1, 0, 6).unwrap();
        options.audio_clock = CinematicClock::try_new(
            CinematicClockDomain::Audio,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            1,
            0,
            AUDIO_FRAMES as i64,
        )
        .unwrap();
        options.declared_source_bandwidth_hz = source_frequency_hz;
        options.filter.half_length = SPECTRAL_FILTER_HALF_LENGTH;
        options.filter.response_grid_intervals = 8_200;
        let fixture = try_fixture_with_inputs(options, inputs, None, false, cx).unwrap();
        let chunk = fixture
            .resampler
            .resample_next_chunk(
                &fixture.sound,
                &fixture.resampler.initial_checkpoint(cx).unwrap(),
                NonZeroUsize::new(AUDIO_FRAMES).unwrap(),
                cx,
            )
            .unwrap();

        let mut held_source = Vec::with_capacity(AUDIO_FRAMES);
        for interval in 0..AUDIO_FRAMES / SOURCE_INTERVAL_FRAMES {
            let midpoint_frame =
                (interval * SOURCE_INTERVAL_FRAMES) as f64 + 0.5 * SOURCE_INTERVAL_FRAMES as f64;
            let value = (TAU * source_frequency_hz * midpoint_frame
                / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ))
            .sin();
            held_source.extend(core::iter::repeat_n(value, SOURCE_INTERVAL_FRAMES));
        }
        let start = usize::from(SPECTRAL_FILTER_HALF_LENGTH);
        let end = AUDIO_FRAMES - start;
        let reconstructed = chunk.drive_frames[start..end]
            .iter()
            .map(|frame| frame.distributed_generalized_force_n.disc)
            .collect::<Vec<_>>();
        let held = &held_source[start..end];
        let held_source_tone = tone_amplitude(held, source_frequency_hz);
        let reconstructed_source_tone = tone_amplitude(&reconstructed, source_frequency_hz);
        let held_image = tone_amplitude(held, image_frequency_hz);
        let reconstructed_image = tone_amplitude(&reconstructed, image_frequency_hz);

        assert!(
            held_source_tone > 0.9,
            "held source tone={held_source_tone:.6e}"
        );
        assert!(held_image > 5.0e-3, "held sampling image={held_image:.6e}");
        let pass_ratio = reconstructed_source_tone / held_source_tone;
        let stop_ratio = reconstructed_image / held_image;
        assert!(
            (0.99..=1.01).contains(&pass_ratio),
            "source-tone magnitude ratio={pass_ratio:.6e}"
        );
        assert!(
            stop_ratio <= 2.0e-4,
            "sampling-image magnitude ratio={stop_ratio:.6e}"
        );
    });
}

#[test]
fn g4_g5_chunk_resume_replay_and_public_cancellation_are_exact() {
    let fixture = with_cx(false, |cx| {
        try_fixture(FixtureOptions::reference(), cx).unwrap()
    });
    with_cx(false, |cx| {
        let initial = fixture.resampler.initial_checkpoint(cx).unwrap();
        let one_shot = fixture
            .resampler
            .resample_next_chunk(
                &fixture.sound,
                &initial,
                NonZeroUsize::new(2_000).unwrap(),
                cx,
            )
            .unwrap();
        let replay = fixture
            .resampler
            .resample_next_chunk(
                &fixture.sound,
                &initial,
                NonZeroUsize::new(2_000).unwrap(),
                cx,
            )
            .unwrap();
        assert_eq!(one_shot, replay);

        let first = fixture
            .resampler
            .resample_next_chunk(&fixture.sound, &initial, NonZeroUsize::new(1).unwrap(), cx)
            .unwrap();
        let second = fixture
            .resampler
            .resample_next_chunk(
                &fixture.sound,
                &first.successor,
                NonZeroUsize::new(2_000).unwrap(),
                cx,
            )
            .unwrap();
        let mut joined_drive = first.drive_frames.clone();
        joined_drive.extend_from_slice(&second.drive_frames);
        assert_eq!(joined_drive, one_shot.drive_frames);
        let mut joined_force = first.preparticipated_localized_force_n.clone();
        joined_force.extend_from_slice(&second.preparticipated_localized_force_n);
        assert_eq!(joined_force, one_shot.preparticipated_localized_force_n);
        let mut joined_impulse = first.preparticipated_localized_impulse_n_s.clone();
        joined_impulse.extend_from_slice(&second.preparticipated_localized_impulse_n_s);
        assert_eq!(
            joined_impulse,
            one_shot.preparticipated_localized_impulse_n_s
        );
        let mut joined_events = first.events.clone();
        joined_events.extend_from_slice(&second.events);
        assert_eq!(joined_events, one_shot.events);
        // The half-sample event receipt belongs to the chunk containing its
        // left boundary, while its conserved right-half impulse is published
        // by the next chunk. Concatenation must still be exactly one-shot.
        assert_eq!(first.events, one_shot.events);
        assert!(second.events.is_empty());
        assert_close(
            first.preparticipated_localized_impulse_n_s[0],
            0.5 * fixture.expected_projected_event_impulse_n_s,
        );
        assert_close(
            second.preparticipated_localized_impulse_n_s[0],
            0.5 * fixture.expected_projected_event_impulse_n_s,
        );
        let mut joined_markers = first.sync_markers.clone();
        joined_markers.extend_from_slice(&second.sync_markers);
        assert_eq!(joined_markers, one_shot.sync_markers);
        assert_eq!(second.successor, one_shot.successor);

        let modal_initial = fixture.modal.initial_checkpoint(cx).unwrap();
        let modal_first = first
            .synthesize_modal(&fixture.modal, &modal_initial, cx)
            .unwrap();
        assert!(matches!(
            second.synthesize_modal(&fixture.modal, &modal_initial, cx),
            Err(ModalSynthesisError::InvalidCheckpoint)
        ));
        let modal_second = second
            .synthesize_modal(&fixture.modal, &modal_first.successor, cx)
            .unwrap();
        assert_eq!(modal_second.successor.next_sample_frame(), 2_000);

        assert!(matches!(
            fixture.resampler.resample_next_chunk(
                &fixture.sound,
                &one_shot.successor,
                NonZeroUsize::new(1).unwrap(),
                cx,
            ),
            Err(AudioResamplingError::Complete)
        ));

        let mut mismatched_input = fixture.sound.input().clone();
        mismatched_input.resampler_identity = identity("wrong-resampler");
        let mismatched = SoundSynthesisConfig::try_admit(mismatched_input).unwrap();
        assert!(matches!(
            fixture.resampler.resample_next_chunk(
                &mismatched,
                &initial,
                NonZeroUsize::new(1).unwrap(),
                cx,
            ),
            Err(AudioResamplingError::SoundConfigurationMismatch(
                "resampler or filter identity/version"
            ))
        ));
    });
    let cancellation_checkpoint = with_cx(false, |cx| {
        fixture.resampler.initial_checkpoint(cx).unwrap()
    });
    with_cx(true, |cx| {
        assert!(matches!(
            fixture.resampler.resample_next_chunk(
                &fixture.sound,
                &cancellation_checkpoint,
                NonZeroUsize::new(2_000).unwrap(),
                cx,
            ),
            Err(AudioResamplingError::Cancelled)
        ));
    });
}

#[test]
fn g0_invalid_bandwidth_clock_filter_and_budget_refuse_typed() {
    with_cx(false, |cx| {
        assert!(matches!(
            try_mismatched_modal_admission(cx),
            Err(AudioResamplingError::ExcitationModalMismatch)
        ));

        let mut options = FixtureOptions::reference();
        options.declared_source_bandwidth_hz = 10_000.0;
        assert!(matches!(
            try_fixture(options, cx),
            Err(AudioResamplingError::UnsupportedSourceBandwidth { .. })
        ));

        let mut options = FixtureOptions::reference();
        options.audio_clock = CinematicClock::try_new(
            CinematicClockDomain::Audio,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            1,
            AUDIO_START_TICK,
            AUDIO_END_TICK - 1,
        )
        .unwrap();
        assert!(matches!(
            try_fixture(options, cx),
            Err(AudioResamplingError::InvalidMasterClock(
                "audio/video endpoints differ"
            ))
        ));

        let mut options = FixtureOptions::reference();
        options.budget.maximum_sync_markers = 1;
        assert!(matches!(
            try_fixture(options, cx),
            Err(AudioResamplingError::BudgetExceeded {
                artifact: "A/V synchronization markers",
                requested: 2,
                limit: 1,
            })
        ));

        let mut options = FixtureOptions::reference();
        options.filter.stopband_edge_hz = options.filter.passband_edge_hz;
        assert!(matches!(
            try_fixture(options, cx),
            Err(AudioResamplingError::InvalidFilter(
                "ordered pass/stop edges within source and target Nyquist"
            ))
        ));

        let mut options = FixtureOptions::reference();
        options.budget.maximum_filter_taps = 128;
        assert!(matches!(
            try_fixture(options, cx),
            Err(AudioResamplingError::BudgetExceeded {
                artifact: "filter taps",
                requested: 2_049,
                limit: 128,
            })
        ));
    });
}
