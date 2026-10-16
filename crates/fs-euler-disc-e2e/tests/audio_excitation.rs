//! G0/G3/G4/G5 public integration checks for source-clock audio excitation.

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::{
    ArtisticTextureConfig, AudioExcitationBudget, AudioExcitationError, AudioExcitationMapper,
    AudioExcitationModelInput, AudioExcitationReconstructionStatus, AudioExcitationReduction,
    AudioExcitationStems, ChannelControl, ContactEventMeasure, ContactModeShape,
    ContactParticipationPolicy, DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
    EulerControlStream, EulerRenderTrajectoryArtifact, ExcitationSourceAvailability,
    ModalComponentValues, ModalSynthesisBudget, ModalSynthesisModel, ModalSynthesisModelInput,
    ModeContactParticipationRule, RenderBaseFrame, RenderBaseModeState, RenderChannelAvailability,
    RenderContactBranch, RenderContactGeometry, RenderContactTransition, RenderMassProperties,
    RenderNumericalRefusalReason, RenderSampleDisposition, RenderSupportFeature, RenderTrajectory,
    RenderTrajectoryAuthority, RenderTrajectoryCodecBudget, RenderTrajectoryMetadata,
    RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame, SpatialEnvelopeSource,
    coupled_runner::{ChannelOwnership, ChannelWrench, ContactTransitionKind},
};
use fs_evidence::cinematic_sound::{
    SoundExcitationChannel, SoundExcitationControl, SoundModalComponent, SoundMode,
    SoundModeParticipation,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

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
                seed: 0x4155_4449_4f45_5843,
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
    hash_domain("org.frankensim.test.audio-excitation.v1", label.as_bytes())
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

fn refresh_pose_diagnostics(input: &mut RenderTrajectorySampleInput) {
    let state = state_from(input);
    input.symmetry_axis_world = state
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    input.qois =
        DerivedEulerQois::from_state(state, mass(), input.qois.precession_acceleration_rad_per_s2)
            .unwrap();
}

fn sample(
    start_time_s: f64,
    end_time_s: f64,
    branch: RenderContactBranch,
    disposition: RenderSampleDisposition,
) -> RenderTrajectorySampleInput {
    let positive_duration = end_time_s > start_time_s;
    let contact_active = positive_duration && branch == RenderContactBranch::Closed;
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.2).unwrap();
    let mut input = RenderTrajectorySampleInput {
        interval_start_time_s: start_time_s,
        time_s: end_time_s,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m: Vec3::new(0.0, 0.0, 1.0),
        orientation_body_to_world: orientation.components(),
        linear_momentum_world_kg_m_per_s: Vec3::ZERO,
        angular_momentum_body_kg_m2_per_s: Vec3::ZERO,
        symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        contact_branch: branch,
        contact_geometry: (branch == RenderContactBranch::Closed).then_some(
            RenderContactGeometry {
                point_world_m: Vec3::new(1.0, 0.0, 0.0),
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
        interval_normal_force_n: if contact_active { 1.0 } else { 0.0 },
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
    refresh_pose_diagnostics(&mut input);
    input
}

fn metadata(
    first: &RenderTrajectorySampleInput,
    availability: RenderChannelAvailability,
    timestep_s: f64,
) -> RenderTrajectoryMetadata {
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
        initial_state: state_from(first),
        initial_base_mode: first.base_mode.unwrap(),
        base_model_identity: identity("base"),
        base_frame: RenderBaseFrame {
            origin_world_m: Vec3::ZERO,
            orientation_base_to_world: UnitQuaternion::IDENTITY,
        },
        model_identity: identity("mechanics-model"),
        channel_availability: availability,
        configuration_identity: identity("configuration"),
        configuration_fingerprint: 0x4155_4449_4f45_5801,
        timestep_s,
        producer_version: "audio-excitation-test-v1".into(),
        applicability: "source-clock excitation integration tests".into(),
        no_claims: vec![
            "source cadence is not an audio sample rate".into(),
            "timing-only events carry no physical impulse".into(),
        ],
        authority: RenderTrajectoryAuthority::SimulationEvidence,
    }
}

fn trajectory(
    inputs: Vec<RenderTrajectorySampleInput>,
    availability: RenderChannelAvailability,
) -> RenderTrajectory {
    let timestep_s = inputs
        .iter()
        .map(|input| input.time_s - input.interval_start_time_s)
        .fold(f64::MIN_POSITIVE, f64::max);
    let trajectory_metadata = metadata(&inputs[0], availability, timestep_s);
    RenderTrajectory::try_new(trajectory_metadata, inputs).unwrap()
}

fn artifact(
    inputs: Vec<RenderTrajectorySampleInput>,
    availability: RenderChannelAvailability,
    cx: &Cx<'_>,
) -> EulerRenderTrajectoryArtifact {
    EulerRenderTrajectoryArtifact::try_from_trajectory(
        identity("campaign"),
        trajectory(inputs, availability),
        Vec::new(),
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .unwrap()
}

fn wrench(force_world_n: Vec3, work_j: f64) -> ChannelWrench {
    ChannelWrench {
        force_world_n,
        torque_world_nm: Vec3::ZERO,
        work_j,
    }
}

fn values(disc: f64, glass_plate: f64, base_assembly: f64) -> ModalComponentValues {
    ModalComponentValues {
        disc,
        glass_plate,
        base_assembly,
    }
}

fn assert_close(actual: f64, expected: f64) {
    let scale = 1.0_f64.max(actual.abs()).max(expected.abs());
    assert!(
        (actual - expected).abs() <= 2.0e-12 * scale,
        "actual={actual:.17e}, expected={expected:.17e}, delta={:.17e}",
        actual - expected
    );
}

fn assert_values(actual: ModalComponentValues, expected: ModalComponentValues) {
    assert_close(actual.disc, expected.disc);
    assert_close(actual.glass_plate, expected.glass_plate);
    assert_close(actual.base_assembly, expected.base_assembly);
}

fn mode(mode_id: u32, component: SoundModalComponent) -> SoundMode {
    SoundMode {
        mode_id,
        component,
        frequency_hz: 400.0 + 100.0 * f64::from(mode_id),
        damping_ratio: 0.02,
        modal_mass_kg: 0.2,
        source_participation: match component {
            SoundModalComponent::Disc => SoundModeParticipation {
                disc: 1.0,
                glass_plate: 0.0,
                base_assembly: 0.0,
            },
            SoundModalComponent::GlassPlate => SoundModeParticipation {
                disc: 0.0,
                glass_plate: 1.0,
                base_assembly: 0.0,
            },
            SoundModalComponent::BaseAssembly => SoundModeParticipation {
                disc: 0.0,
                glass_plate: 0.0,
                base_assembly: 1.0,
            },
        },
        radiation_gain_fs_s_per_m: 0.1,
        material_identity: identity("material"),
        base_identity: identity("modal-base"),
    }
}

fn modal_model(components: &[SoundModalComponent], cx: &Cx<'_>) -> ModalSynthesisModel {
    let modes = components
        .iter()
        .copied()
        .enumerate()
        .map(|(index, component)| mode(u32::try_from(index + 1).unwrap(), component))
        .collect();
    ModalSynthesisModel::try_new(
        ModalSynthesisModelInput {
            sample_rate_hz: 48_000,
            modes,
            budget: ModalSynthesisBudget::reference_film(1_024),
        },
        cx,
    )
    .unwrap()
}

fn mapping(
    channel: SoundExcitationChannel,
    target_component: SoundModalComponent,
    source_scale: f64,
) -> SoundExcitationControl {
    SoundExcitationControl {
        channel,
        target_component,
        source_scale,
    }
}

fn all_mappings() -> Vec<SoundExcitationControl> {
    vec![
        mapping(
            SoundExcitationChannel::ContactNormalForce,
            SoundModalComponent::Disc,
            2.0,
        ),
        mapping(
            SoundExcitationChannel::ContactSignedWorkRate,
            SoundModalComponent::GlassPlate,
            -3.0,
        ),
        mapping(
            SoundExcitationChannel::RollingSignedWorkRate,
            SoundModalComponent::BaseAssembly,
            4.0,
        ),
        mapping(
            SoundExcitationChannel::BaseDampingSignedWorkRate,
            SoundModalComponent::Disc,
            -5.0,
        ),
        mapping(
            SoundExcitationChannel::ExteriorGasBodySignedWorkRate,
            SoundModalComponent::GlassPlate,
            6.0,
        ),
    ]
}

fn mapper_input(
    mappings: Vec<SoundExcitationControl>,
    raw_interval_count: usize,
) -> AudioExcitationModelInput {
    AudioExcitationModelInput {
        mappings,
        reduction: AudioExcitationReduction::RawIntervals,
        spatial_policy: ContactParticipationPolicy::DeclaredStatic,
        artistic_texture: None,
        budget: AudioExcitationBudget::reference_film(raw_interval_count),
    }
}

fn alternating_inputs() -> Vec<RenderTrajectorySampleInput> {
    let mut inputs = vec![sample(
        0.0,
        0.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    )];
    for index in 0..4 {
        let mut retained = sample(
            f64::from(index),
            f64::from(index + 1),
            RenderContactBranch::Open,
            if index == 3 {
                RenderSampleDisposition::HorizonCensored
            } else {
                RenderSampleDisposition::Continue
            },
        );
        let sign = if index % 2 == 0 { 1.0 } else { -1.0 };
        retained.channels.gas = wrench(Vec3::ZERO, sign);
        inputs.push(retained);
    }
    inputs
}

fn event_inputs() -> Vec<RenderTrajectorySampleInput> {
    let first = sample(
        0.0,
        0.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let mut reimpact = sample(
        0.0,
        1.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::Continue,
    );
    reimpact.channels.rolling = wrench(Vec3::ZERO, -0.5);
    reimpact.contact_transitions.push(RenderContactTransition {
        kind: ContactTransitionKind::Reimpact,
        time_s: 0.5,
        bracket_start_s: 0.49,
        bracket_end_s: 0.51,
    });
    let mut opening = sample(
        1.0,
        2.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    opening.interval_contact_active = true;
    opening.channels.rolling = wrench(Vec3::ZERO, -0.25);
    opening.contact_transitions.push(RenderContactTransition {
        kind: ContactTransitionKind::Opening,
        time_s: 1.5,
        bracket_start_s: 1.49,
        bracket_end_s: 1.51,
    });
    vec![first, reimpact, opening]
}

#[test]
fn g0_all_supported_mappings_preserve_signed_units_stems_and_measures() {
    let first = sample(
        0.0,
        0.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::Continue,
    );
    let mut retained = sample(
        0.0,
        1.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::HorizonCensored,
    );
    retained.interval_normal_force_n = 30.0;
    retained.channels = ChannelOwnership {
        contact: wrench(Vec3::new(4.0, -3.0, 30.0), -2.0),
        rolling: wrench(Vec3::ZERO, -0.5),
        base: wrench(Vec3::ZERO, -0.25),
        gas: wrench(Vec3::ZERO, 0.125),
        ..ChannelOwnership::default()
    };

    with_cx(false, |cx| {
        let artifact = artifact(
            vec![first, retained],
            RenderChannelAvailability::ALL_AVAILABLE,
            cx,
        );
        let controls = EulerControlStream::try_derive(artifact.trajectory(), cx).unwrap();
        let modal = modal_model(
            &[
                SoundModalComponent::Disc,
                SoundModalComponent::GlassPlate,
                SoundModalComponent::BaseAssembly,
            ],
            cx,
        );
        let mapper = AudioExcitationMapper::try_new(
            &artifact,
            &controls,
            &modal,
            mapper_input(all_mappings(), 1),
            cx,
        )
        .unwrap();
        let initial = mapper.initial_checkpoint(cx).unwrap();
        let chunk = mapper
            .map_next_chunk(&initial, NonZeroUsize::new(1).unwrap(), cx)
            .unwrap();
        let interval = &chunk.intervals[0];

        assert_eq!(
            SoundExcitationChannel::ContactNormalForce.source_unit(),
            "N"
        );
        for channel in [
            SoundExcitationChannel::ContactSignedWorkRate,
            SoundExcitationChannel::RollingSignedWorkRate,
            SoundExcitationChannel::BaseDampingSignedWorkRate,
            SoundExcitationChannel::ExteriorGasBodySignedWorkRate,
        ] {
            assert_eq!(channel.source_unit(), "W");
        }
        assert_values(interval.mean_force_stems_n.contact, values(60.0, 6.0, 0.0));
        assert_values(interval.mean_force_stems_n.rolling, values(0.0, 0.0, -2.0));
        assert_values(interval.mean_force_stems_n.base, values(1.25, 0.0, 0.0));
        assert_values(interval.mean_force_stems_n.gas, values(0.0, 0.75, 0.0));
        assert_eq!(interval.force_time_stems_n_s, interval.mean_force_stems_n);
        assert_values(interval.mean_generalized_force_n, values(61.25, 6.75, -2.0));
        assert_eq!(
            interval.generalized_force_time_n_s,
            interval.mean_generalized_force_n
        );
        assert_values(
            interval.localized_mean_generalized_force_n(),
            values(60.0, 6.0, -2.0),
        );
        assert_values(
            interval.distributed_mean_generalized_force_n(),
            values(1.25, 0.75, 0.0),
        );
        assert_eq!(
            interval.localized_force_time_measure_n_s(),
            interval.localized_mean_generalized_force_n()
        );
        assert_eq!(
            interval.distributed_force_time_measure_n_s(),
            interval.distributed_mean_generalized_force_n()
        );
        assert_eq!(
            interval.measure_residual_stems_n_s,
            AudioExcitationStems::ZERO
        );
        assert_eq!(interval.measure_residual_n_s, ModalComponentValues::ZERO);
        assert_eq!(
            interval.availability.contact,
            ExcitationSourceAvailability::Available
        );
        assert_eq!(
            interval.availability.rolling,
            ExcitationSourceAvailability::Available
        );
        assert_eq!(
            interval.availability.base,
            ExcitationSourceAvailability::Available
        );
        assert_eq!(
            interval.availability.gas,
            ExcitationSourceAvailability::Available
        );
        assert_eq!(chunk.diagnostics.maximum_abs_measure_residual_n_s, 0.0);
        assert_eq!(mapper.grid().interval_count, 1);
        assert_eq!(mapper.grid().nominal_source_nyquist_ceiling_hz, 0.5);
        assert_eq!(
            mapper.grid().reconstruction,
            AudioExcitationReconstructionStatus::RequiresBandLimitedResampling
        );
    });
}

#[test]
fn g0_normal_load_only_authority_maps_force_and_measure_without_contact_wrench() {
    let first = sample(
        0.0,
        0.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::Continue,
    );
    let mut retained = sample(
        0.0,
        1.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::HorizonCensored,
    );
    retained.interval_normal_force_n = 7.0;

    with_cx(false, |cx| {
        for sampling in [
            fs_euler_disc_e2e::RenderNormalForceSampling::IntervalMean,
            fs_euler_disc_e2e::RenderNormalForceSampling::AppliedSubstepZeroOrderHold,
        ] {
            let availability = RenderChannelAvailability {
                gravity: false,
                contact: false,
                normal_force_sampling: sampling,
                rolling: false,
                base: false,
                gas: false,
            };
            let artifact = artifact(vec![first.clone(), retained.clone()], availability, cx);
            let controls = EulerControlStream::try_derive(artifact.trajectory(), cx).unwrap();
            assert!(matches!(
                controls.audio()[0].channels.contact,
                ChannelControl::Unavailable
            ));
            assert_eq!(controls.audio()[0].normal_force_sampling, sampling);
            assert_eq!(
                controls.audio()[0].mean_base_normal_contact_force_n,
                Some(7.0)
            );
            let modal = modal_model(&[SoundModalComponent::Disc], cx);
            let mapper = AudioExcitationMapper::try_new(
                &artifact,
                &controls,
                &modal,
                mapper_input(
                    vec![mapping(
                        SoundExcitationChannel::ContactNormalForce,
                        SoundModalComponent::Disc,
                        1.0,
                    )],
                    1,
                ),
                cx,
            )
            .unwrap();
            let chunk = mapper
                .map_next_chunk(
                    &mapper.initial_checkpoint(cx).unwrap(),
                    NonZeroUsize::new(1).unwrap(),
                    cx,
                )
                .unwrap();
            let interval = &chunk.intervals[0];
            assert_eq!(
                interval.availability.contact,
                ExcitationSourceAvailability::Available
            );
            assert_eq!(interval.mean_force_stems_n.contact.disc, 7.0);
            assert_eq!(interval.force_time_stems_n_s.contact.disc, 7.0);
            assert_eq!(interval.measure_residual_stems_n_s.contact.disc, 0.0);
        }
    });
}

#[test]
fn g0_normal_only_midpoint_fails_closed_for_contact_normal_force_sound() {
    let first = sample(
        0.0,
        0.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::Continue,
    );
    let mut retained = sample(
        0.0,
        1.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::HorizonCensored,
    );
    retained.interval_normal_force_n = 7.0;
    let availability = RenderChannelAvailability {
        gravity: false,
        contact: false,
        normal_force_sampling:
            fs_euler_disc_e2e::RenderNormalForceSampling::FirstAcceptedSubintervalMidpoint,
        rolling: false,
        base: false,
        gas: false,
    };

    with_cx(false, |cx| {
        let artifact = artifact(vec![first, retained], availability, cx);
        let controls = EulerControlStream::try_derive(artifact.trajectory(), cx).unwrap();
        assert_eq!(controls.audio()[0].mean_base_normal_contact_force_n, None);
        assert_eq!(
            controls.audio()[0].normal_force_sampling,
            fs_euler_disc_e2e::RenderNormalForceSampling::FirstAcceptedSubintervalMidpoint
        );
        let modal = modal_model(&[SoundModalComponent::Disc], cx);
        let mapper = AudioExcitationMapper::try_new(
            &artifact,
            &controls,
            &modal,
            mapper_input(
                vec![mapping(
                    SoundExcitationChannel::ContactNormalForce,
                    SoundModalComponent::Disc,
                    1.0,
                )],
                1,
            ),
            cx,
        )
        .unwrap();
        let chunk = mapper
            .map_next_chunk(
                &mapper.initial_checkpoint(cx).unwrap(),
                NonZeroUsize::new(1).unwrap(),
                cx,
            )
            .unwrap();
        assert_eq!(
            chunk.intervals[0].availability.contact,
            ExcitationSourceAvailability::Unavailable
        );
        assert_eq!(chunk.intervals[0].mean_force_stems_n.contact.disc, 0.0);
        assert_eq!(chunk.intervals[0].force_time_stems_n_s.contact.disc, 0.0);
    });
}

#[test]
fn g0_available_zero_and_unavailable_sources_remain_distinct() {
    let first = sample(
        0.0,
        0.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let retained = sample(
        0.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    let availability = RenderChannelAvailability {
        gravity: true,
        contact: true,
        normal_force_sampling:
            fs_euler_disc_e2e::RenderNormalForceSampling::FirstAcceptedSubintervalMidpoint,
        rolling: true,
        base: false,
        gas: false,
    };

    with_cx(false, |cx| {
        let artifact = artifact(vec![first, retained], availability, cx);
        let controls = EulerControlStream::try_derive(artifact.trajectory(), cx).unwrap();
        let modal = modal_model(&[SoundModalComponent::Disc], cx);
        let mapper = AudioExcitationMapper::try_new(
            &artifact,
            &controls,
            &modal,
            mapper_input(all_mappings(), 1),
            cx,
        )
        .unwrap();
        let initial = mapper.initial_checkpoint(cx).unwrap();
        let chunk = mapper
            .map_next_chunk(&initial, NonZeroUsize::new(1).unwrap(), cx)
            .unwrap();
        let interval = &chunk.intervals[0];

        assert_eq!(interval.mean_force_stems_n, AudioExcitationStems::ZERO);
        assert_eq!(interval.force_time_stems_n_s, AudioExcitationStems::ZERO);
        assert_eq!(
            interval.availability.contact,
            ExcitationSourceAvailability::Available
        );
        assert_eq!(
            interval.availability.rolling,
            ExcitationSourceAvailability::Available
        );
        assert_eq!(
            interval.availability.base,
            ExcitationSourceAvailability::Unavailable
        );
        assert_eq!(
            interval.availability.gas,
            ExcitationSourceAvailability::Unavailable
        );
    });
}

#[test]
fn g3_boxcar_cancels_alternating_work_before_decimation() {
    with_cx(false, |cx| {
        let source_artifact = artifact(
            alternating_inputs(),
            RenderChannelAvailability::ALL_AVAILABLE,
            cx,
        );
        let controls = EulerControlStream::try_derive(source_artifact.trajectory(), cx).unwrap();
        let modal = modal_model(&[SoundModalComponent::GlassPlate], cx);
        let mappings = vec![mapping(
            SoundExcitationChannel::ExteriorGasBodySignedWorkRate,
            SoundModalComponent::GlassPlate,
            1.0,
        )];

        let raw = AudioExcitationMapper::try_new(
            &source_artifact,
            &controls,
            &modal,
            mapper_input(mappings.clone(), 4),
            cx,
        )
        .unwrap();
        let raw_initial = raw.initial_checkpoint(cx).unwrap();
        let raw_chunk = raw
            .map_next_chunk(&raw_initial, NonZeroUsize::new(4).unwrap(), cx)
            .unwrap();
        assert_eq!(raw.grid().interval_count, 4);
        assert_eq!(
            raw_chunk
                .intervals
                .iter()
                .map(|interval| interval.mean_generalized_force_n.glass_plate)
                .collect::<Vec<_>>(),
            vec![1.0, -1.0, 1.0, -1.0]
        );

        let mut coarsened_input = mapper_input(mappings.clone(), 4);
        coarsened_input.reduction = AudioExcitationReduction::WholeIntervalBoxcarV1 {
            intervals_per_bin: NonZeroUsize::new(2).unwrap(),
        };
        let coarsened = AudioExcitationMapper::try_new(
            &source_artifact,
            &controls,
            &modal,
            coarsened_input.clone(),
            cx,
        )
        .unwrap();
        let coarsened_initial = coarsened.initial_checkpoint(cx).unwrap();
        let coarsened_chunk = coarsened
            .map_next_chunk(&coarsened_initial, NonZeroUsize::new(2).unwrap(), cx)
            .unwrap();
        assert_ne!(raw.identity(), coarsened.identity());
        assert_eq!(coarsened.grid().interval_count, 2);
        assert_eq!(coarsened.grid().minimum_interval_duration_s, 2.0);
        assert_eq!(coarsened.grid().maximum_interval_duration_s, 2.0);
        for interval in &coarsened_chunk.intervals {
            assert_eq!(
                interval.mean_generalized_force_n,
                ModalComponentValues::ZERO
            );
            assert_eq!(
                interval.generalized_force_time_n_s,
                ModalComponentValues::ZERO
            );
            assert_eq!(
                interval.force_time_stems_n_s.gas,
                ModalComponentValues::ZERO
            );
        }
        assert_eq!(
            coarsened_chunk.successor.cumulative_force_time_stems_n_s(),
            AudioExcitationStems::ZERO
        );

        coarsened_input.spatial_policy = ContactParticipationPolicy::ContactCoordinates {
            rules: vec![ModeContactParticipationRule {
                mode_id: 1,
                shape: ContactModeShape::AzimuthalCosine {
                    harmonic: 1,
                    phase_rad: 0.0,
                },
            }],
        };
        assert!(matches!(
            AudioExcitationMapper::try_new(
                &source_artifact,
                &controls,
                &modal,
                coarsened_input,
                cx
            ),
            Err(AudioExcitationError::InvalidSpatialPolicy(
                "source coarsening cannot preserve time-varying contact participation"
            ))
        ));
    });
}

#[test]
fn g3_timing_only_events_keep_zero_physical_impulse_and_artistic_reimpact_is_replayable() {
    with_cx(false, |cx| {
        let artifact = artifact(event_inputs(), RenderChannelAvailability::ALL_AVAILABLE, cx);
        let controls = EulerControlStream::try_derive(artifact.trajectory(), cx).unwrap();
        let modal = modal_model(
            &[
                SoundModalComponent::Disc,
                SoundModalComponent::GlassPlate,
                SoundModalComponent::BaseAssembly,
            ],
            cx,
        );
        let mappings = vec![mapping(
            SoundExcitationChannel::RollingSignedWorkRate,
            SoundModalComponent::Disc,
            1.0,
        )];
        let texture = ArtisticTextureConfig {
            seed: 7,
            rolling_force_gain_n_per_w: 2.0,
            rolling_target_component: SoundModalComponent::GlassPlate,
            band_low_hz: 100.0,
            band_high_hz: 1_000.0,
            reimpact_impulse_n_s: 0.25,
            reimpact_target_component: SoundModalComponent::BaseAssembly,
        };
        let mut input = mapper_input(mappings.clone(), 2);
        input.artistic_texture = Some(texture);
        let mapper =
            AudioExcitationMapper::try_new(&artifact, &controls, &modal, input.clone(), cx)
                .unwrap();
        let initial = mapper.initial_checkpoint(cx).unwrap();
        let chunk = mapper
            .map_next_chunk(&initial, NonZeroUsize::new(2).unwrap(), cx)
            .unwrap();
        let replay = mapper
            .map_next_chunk(&initial, NonZeroUsize::new(2).unwrap(), cx)
            .unwrap();
        assert_eq!(chunk, replay);
        assert_eq!(chunk.diagnostics.event_count, 2);

        let reimpact = &chunk.intervals[0];
        assert!(reimpact.event_barrier);
        assert_eq!(reimpact.events.len(), 1);
        assert_eq!(reimpact.events[0].kind, ContactTransitionKind::Reimpact);
        assert_eq!(reimpact.events[0].measure, ContactEventMeasure::TimingOnly);
        assert_eq!(
            reimpact.events[0].physical_impulse_n_s,
            ModalComponentValues::ZERO
        );
        let artistic = reimpact.events[0].artistic.unwrap();
        assert_eq!(artistic.impulse_n_s.disc, 0.0);
        assert_eq!(artistic.impulse_n_s.glass_plate, 0.0);
        assert!(artistic.impulse_n_s.base_assembly.abs() <= 0.25);
        let envelope = reimpact.artistic_texture.unwrap();
        assert_values(envelope.peak_force_envelope_n, values(0.0, 1.0, 0.0));
        assert_eq!(
            envelope.rolling_availability,
            ExcitationSourceAvailability::Available
        );

        let opening = &chunk.intervals[1];
        assert!(opening.event_barrier);
        assert_eq!(opening.events[0].kind, ContactTransitionKind::Opening);
        assert_eq!(
            opening.events[0].physical_impulse_n_s,
            ModalComponentValues::ZERO
        );
        assert!(opening.events[0].artistic.is_none());
        assert_values(
            opening.artistic_texture.unwrap().peak_force_envelope_n,
            values(0.0, 0.5, 0.0),
        );

        let mut coarsened_input = input;
        coarsened_input.reduction = AudioExcitationReduction::WholeIntervalBoxcarV1 {
            intervals_per_bin: NonZeroUsize::new(8).unwrap(),
        };
        let coarsened =
            AudioExcitationMapper::try_new(&artifact, &controls, &modal, coarsened_input, cx)
                .unwrap();
        assert_eq!(coarsened.grid().interval_count, 2);
        let coarse_initial = coarsened.initial_checkpoint(cx).unwrap();
        let coarse_chunk = coarsened
            .map_next_chunk(&coarse_initial, NonZeroUsize::new(2).unwrap(), cx)
            .unwrap();
        assert!(
            coarse_chunk
                .intervals
                .iter()
                .all(|interval| interval.event_barrier)
        );
        assert_eq!(
            coarse_chunk
                .intervals
                .iter()
                .map(|interval| interval.events[0].kind)
                .collect::<Vec<_>>(),
            vec![
                ContactTransitionKind::Reimpact,
                ContactTransitionKind::Opening
            ]
        );
    });
}

#[test]
fn g3_spatial_envelopes_use_exact_and_event_side_contact_coordinates() {
    with_cx(false, |cx| {
        let first = sample(
            0.0,
            0.0,
            RenderContactBranch::Closed,
            RenderSampleDisposition::Continue,
        );
        let mut retained = sample(
            0.0,
            1.0,
            RenderContactBranch::Closed,
            RenderSampleDisposition::HorizonCensored,
        );
        retained.contact_geometry.as_mut().unwrap().point_world_m = Vec3::new(0.0, 1.0, 0.0);
        retained.channels.contact = wrench(Vec3::new(0.0, 0.0, 1.0), 0.0);
        let exact_artifact = artifact(
            vec![first, retained],
            RenderChannelAvailability::ALL_AVAILABLE,
            cx,
        );
        let exact_controls =
            EulerControlStream::try_derive(exact_artifact.trajectory(), cx).unwrap();
        let modal = modal_model(&[SoundModalComponent::BaseAssembly], cx);
        let mappings = vec![mapping(
            SoundExcitationChannel::ContactNormalForce,
            SoundModalComponent::BaseAssembly,
            1.0,
        )];
        let spatial_policy = ContactParticipationPolicy::ContactCoordinates {
            rules: vec![ModeContactParticipationRule {
                mode_id: 1,
                shape: ContactModeShape::AzimuthalCosine {
                    harmonic: 1,
                    phase_rad: 0.0,
                },
            }],
        };
        let mut exact_input = mapper_input(mappings.clone(), 1);
        exact_input.spatial_policy = spatial_policy.clone();
        let exact_mapper = AudioExcitationMapper::try_new(
            &exact_artifact,
            &exact_controls,
            &modal,
            exact_input,
            cx,
        )
        .unwrap();
        let exact_initial = exact_mapper.initial_checkpoint(cx).unwrap();
        let exact_chunk = exact_mapper
            .map_next_chunk(&exact_initial, NonZeroUsize::new(1).unwrap(), cx)
            .unwrap();
        let envelope = exact_chunk.intervals[0].spatial_envelopes[0];
        assert_eq!(
            envelope.source,
            SpatialEnvelopeSource::ExactEndpointInterpolation
        );
        assert_close(envelope.start_factor, 1.0);
        assert_close(envelope.end_factor, 0.0);
        assert_close(
            exact_chunk.intervals[0].spatial_factors_at(0.5).unwrap()[0],
            0.5,
        );
        assert_eq!(
            exact_chunk.intervals[0].spatial_factors_at(-0.1),
            Err(AudioExcitationError::InvalidSpatialEvaluation)
        );
        assert_eq!(
            exact_chunk.intervals[0].spatial_factors_at(f64::NAN),
            Err(AudioExcitationError::InvalidSpatialEvaluation)
        );

        let event_artifact = artifact(event_inputs(), RenderChannelAvailability::ALL_AVAILABLE, cx);
        let event_controls =
            EulerControlStream::try_derive(event_artifact.trajectory(), cx).unwrap();
        let mut event_input = mapper_input(mappings, 2);
        event_input.spatial_policy = spatial_policy;
        let event_mapper = AudioExcitationMapper::try_new(
            &event_artifact,
            &event_controls,
            &modal,
            event_input.clone(),
            cx,
        )
        .unwrap();
        let event_initial = event_mapper.initial_checkpoint(cx).unwrap();
        let event_chunk = event_mapper
            .map_next_chunk(&event_initial, NonZeroUsize::new(2).unwrap(), cx)
            .unwrap();
        assert_eq!(
            event_chunk.intervals[0].spatial_envelopes[0].source,
            SpatialEnvelopeSource::HeldEndEndpoint
        );
        assert_eq!(
            event_chunk.intervals[1].spatial_envelopes[0].source,
            SpatialEnvelopeSource::HeldStartEndpoint
        );
        for interval in &event_chunk.intervals {
            let held = interval.spatial_envelopes[0];
            assert_eq!(held.start_factor.to_bits(), held.end_factor.to_bits());
        }

        event_input.reduction = AudioExcitationReduction::WholeIntervalBoxcarV1 {
            intervals_per_bin: NonZeroUsize::new(2).unwrap(),
        };
        assert!(matches!(
            AudioExcitationMapper::try_new(
                &event_artifact,
                &event_controls,
                &modal,
                event_input,
                cx
            ),
            Err(AudioExcitationError::InvalidSpatialPolicy(
                "source coarsening cannot preserve time-varying contact participation"
            ))
        ));
    });
}

#[test]
fn g5_checkpoint_split_replay_complete_and_wrong_mapper_are_exact() {
    with_cx(false, |cx| {
        let source_artifact = artifact(
            alternating_inputs(),
            RenderChannelAvailability::ALL_AVAILABLE,
            cx,
        );
        let controls = EulerControlStream::try_derive(source_artifact.trajectory(), cx).unwrap();
        let modal = modal_model(&[SoundModalComponent::GlassPlate], cx);
        let mappings = vec![mapping(
            SoundExcitationChannel::ExteriorGasBodySignedWorkRate,
            SoundModalComponent::GlassPlate,
            1.0,
        )];
        let mapper = AudioExcitationMapper::try_new(
            &source_artifact,
            &controls,
            &modal,
            mapper_input(mappings.clone(), 4),
            cx,
        )
        .unwrap();
        let initial = mapper.initial_checkpoint(cx).unwrap();
        let one_shot = mapper
            .map_next_chunk(&initial, NonZeroUsize::new(4).unwrap(), cx)
            .unwrap();
        let replay = mapper
            .map_next_chunk(&initial, NonZeroUsize::new(4).unwrap(), cx)
            .unwrap();
        assert_eq!(one_shot, replay);

        let first = mapper
            .map_next_chunk(&initial, NonZeroUsize::new(1).unwrap(), cx)
            .unwrap();
        let second = mapper
            .map_next_chunk(&first.successor, NonZeroUsize::new(3).unwrap(), cx)
            .unwrap();
        let mut joined = first.intervals;
        joined.extend_from_slice(&second.intervals);
        assert_eq!(joined, one_shot.intervals);
        assert_eq!(second.successor, one_shot.successor);
        assert_eq!(
            second.diagnostics.cumulative_force_time_stems_n_s,
            one_shot.diagnostics.cumulative_force_time_stems_n_s
        );
        assert_eq!(
            mapper.map_next_chunk(&one_shot.successor, NonZeroUsize::new(1).unwrap(), cx),
            Err(AudioExcitationError::Complete)
        );

        let other_mapper = AudioExcitationMapper::try_new(
            &source_artifact,
            &controls,
            &modal,
            mapper_input(
                vec![mapping(
                    SoundExcitationChannel::ExteriorGasBodySignedWorkRate,
                    SoundModalComponent::GlassPlate,
                    2.0,
                )],
                4,
            ),
            cx,
        )
        .unwrap();
        assert_eq!(
            other_mapper.map_next_chunk(&initial, NonZeroUsize::new(1).unwrap(), cx),
            Err(AudioExcitationError::InvalidCheckpoint)
        );

        with_cx(true, |cancelled_cx| {
            assert_eq!(
                mapper.map_next_chunk(&initial, NonZeroUsize::new(1).unwrap(), cancelled_cx),
                Err(AudioExcitationError::Cancelled)
            );
        });
    });
}

#[test]
fn g0_unsupported_selectors_caps_and_numerical_refusal_fail_closed() {
    with_cx(false, |cx| {
        let source_artifact = artifact(
            alternating_inputs(),
            RenderChannelAvailability::ALL_AVAILABLE,
            cx,
        );
        let controls = EulerControlStream::try_derive(source_artifact.trajectory(), cx).unwrap();
        let modal = modal_model(&[SoundModalComponent::Disc], cx);

        let unsupported = mapper_input(
            vec![mapping(
                SoundExcitationChannel::ContactTangentialForce,
                SoundModalComponent::Disc,
                1.0,
            )],
            4,
        );
        assert!(matches!(
            AudioExcitationMapper::try_new(&source_artifact, &controls, &modal, unsupported, cx),
            Err(AudioExcitationError::UnsupportedMapping(
                SoundExcitationChannel::ContactTangentialForce
            ))
        ));

        let mappings = vec![mapping(
            SoundExcitationChannel::ExteriorGasBodySignedWorkRate,
            SoundModalComponent::Disc,
            1.0,
        )];
        let mut total_cap = mapper_input(mappings.clone(), 3);
        total_cap.budget.maximum_total_intervals = 3;
        assert!(matches!(
            AudioExcitationMapper::try_new(&source_artifact, &controls, &modal, total_cap, cx),
            Err(AudioExcitationError::LimitExceeded {
                field: "raw source intervals",
                magnitude: 4.0,
                limit: 3.0,
                ..
            })
        ));

        let mut chunk_cap = mapper_input(mappings.clone(), 4);
        chunk_cap.budget.maximum_chunk_intervals = 1;
        let chunk_capped =
            AudioExcitationMapper::try_new(&source_artifact, &controls, &modal, chunk_cap, cx)
                .unwrap();
        let chunk_initial = chunk_capped.initial_checkpoint(cx).unwrap();
        assert_eq!(
            chunk_capped.map_next_chunk(&chunk_initial, NonZeroUsize::new(2).unwrap(), cx),
            Err(AudioExcitationError::ChunkIntervalBudgetExceeded {
                requested: 2,
                limit: 1,
            })
        );

        let mut spatial_cap = mapper_input(mappings, 4);
        spatial_cap.budget.maximum_chunk_spatial_envelopes = 3;
        let spatial_capped =
            AudioExcitationMapper::try_new(&source_artifact, &controls, &modal, spatial_cap, cx)
                .unwrap();
        let spatial_initial = spatial_capped.initial_checkpoint(cx).unwrap();
        assert_eq!(
            spatial_capped.map_next_chunk(&spatial_initial, NonZeroUsize::new(4).unwrap(), cx),
            Err(AudioExcitationError::ChunkSpatialEnvelopeBudgetExceeded {
                requested: 4,
                limit: 3,
            })
        );

        let first = sample(
            0.0,
            0.0,
            RenderContactBranch::Open,
            RenderSampleDisposition::Continue,
        );
        let refused = sample(
            0.0,
            1.0,
            RenderContactBranch::Open,
            RenderSampleDisposition::NumericalRefusal(
                RenderNumericalRefusalReason::NonFiniteEnergyOrBaseState,
            ),
        );
        let refused_artifact = artifact(
            vec![first, refused],
            RenderChannelAvailability::ALL_AVAILABLE,
            cx,
        );
        let refused_controls =
            EulerControlStream::try_derive(refused_artifact.trajectory(), cx).unwrap();
        let refused_input = mapper_input(
            vec![mapping(
                SoundExcitationChannel::ContactNormalForce,
                SoundModalComponent::Disc,
                1.0,
            )],
            1,
        );
        assert!(matches!(
            AudioExcitationMapper::try_new(
                &refused_artifact,
                &refused_controls,
                &modal,
                refused_input,
                cx
            ),
            Err(AudioExcitationError::NumericalRefusalSource)
        ));
    });
}
