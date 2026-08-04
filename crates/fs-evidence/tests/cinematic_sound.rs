//! G0/G3 admission and identity checks for cinematic sound configuration.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::cinematic::{
    CinematicClock, CinematicClockDomain, DeclaredAcousticCalibrationReceipt, SoundAuthority,
};
use fs_evidence::cinematic_config::{CinematicComponentRef, CinematicComponentRole};
use fs_evidence::cinematic_sound::{
    ListenerFrame, ListenerPose, SOUND_MASTER_SAMPLE_RATE_HZ, SOUND_SYNTHESIS_SCHEMA_VERSION,
    SoundAmplitudeReference, SoundChannelLayout, SoundExcitationChannel, SoundExcitationControl,
    SoundMode, SoundModelAssumption, SoundOutputKind, SoundRoomResponse, SoundSynthesisConfig,
    SoundSynthesisError, SoundSynthesisInput, SoundTerminalPolicy, SoundTrajectoryDisposition,
    synthesis_class_code,
};

fn identity(label: &str) -> ContentHash {
    hash_domain("org.frankensim.test.sound-config", label.as_bytes())
}

fn component(role: CinematicComponentRole) -> CinematicComponentRef {
    CinematicComponentRef::try_new(role, identity(&format!("{role:?}")), 1).unwrap()
}

fn calibration() -> DeclaredAcousticCalibrationReceipt {
    DeclaredAcousticCalibrationReceipt::try_new(
        identity("calibration-data"),
        identity("calibration-method"),
        identity("calibration-domain"),
        1,
    )
    .unwrap()
}

fn informed_input() -> SoundSynthesisInput {
    SoundSynthesisInput {
        schema_version: SOUND_SYNTHESIS_SCHEMA_VERSION,
        authority: SoundAuthority::PhysicallyInformed,
        trajectory: component(CinematicComponentRole::Trajectory),
        excitation: component(CinematicComponentRole::AudioExcitation),
        sound_model: component(CinematicComponentRole::SoundModel),
        microphone: component(CinematicComponentRole::Microphone),
        room: component(CinematicComponentRole::Room),
        timeline: component(CinematicComponentRole::Timeline),
        video_clock: CinematicClock::try_new(CinematicClockDomain::Video, 24, 1, 0, 240).unwrap(),
        audio_clock: CinematicClock::try_new(
            CinematicClockDomain::Audio,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            1,
            0,
            480_000,
        )
        .unwrap(),
        channel_layout: SoundChannelLayout::Stereo,
        listener: ListenerPose {
            frame: ListenerFrame::AnimatedCamera,
            position_m: [0.0, -0.08, 0.0],
            forward: [0.0, 1.0, 0.0],
            up: [0.0, 0.0, 1.0],
        },
        excitation_controls: vec![
            SoundExcitationControl {
                channel: SoundExcitationChannel::ContactNormalForce,
                source_scale: 1.0e-3,
            },
            SoundExcitationControl {
                channel: SoundExcitationChannel::BaseReactionForce,
                source_scale: -2.0e-4,
            },
        ],
        modes: vec![
            SoundMode {
                mode_id: 1,
                frequency_hz: 733.0,
                damping_ratio: 0.012,
                gain: 0.18,
                material_identity: identity("steel-material"),
                base_identity: identity("glass-base"),
            },
            SoundMode {
                mode_id: 2,
                frequency_hz: 1_421.0,
                damping_ratio: 0.019,
                gain: -0.07,
                material_identity: identity("steel-material"),
                base_identity: identity("glass-base"),
            },
        ],
        room_response: SoundRoomResponse::DeclaredImpulseResponse { wet_mix: 0.08 },
        amplitude_reference: SoundAmplitudeReference::DigitalFullScale { headroom_db: 6.0 },
        trajectory_disposition: SoundTrajectoryDisposition::PhysicalTerminal,
        terminal_policy: SoundTerminalPolicy::FadeAtLastAccepted {
            fade_sample_frames: 2_400,
        },
        resampler_identity: identity("windowed-sinc-resampler"),
        resampler_version: 1,
        filter_identity: identity("anti-alias-filter"),
        filter_version: 3,
        assumptions: vec![
            SoundModelAssumption::LinearModalSuperposition,
            SoundModelAssumption::TimeInvariantDamping,
            SoundModelAssumption::DeclaredExcitationCompleteness,
            SoundModelAssumption::DeclaredRoomResponse,
        ],
        calibration: None,
    }
}

#[test]
fn admits_physically_informed_reference_and_metadata_matches() {
    let config = SoundSynthesisConfig::try_admit(informed_input()).unwrap();
    assert_eq!(config.authority(), SoundAuthority::PhysicallyInformed);
    assert_eq!(
        synthesis_class_code(config.authority()),
        "physically-informed-modal-synthesis"
    );
    let manifest = config.to_manifest_json();
    assert!(manifest.contains("\"authority\":\"physically-informed\""));
    assert!(manifest.contains("\"synthesis_class\":\"physically-informed-modal-synthesis\""));
    assert!(manifest.contains("\"sample_rate_hz\":48000"));
    assert!(manifest.contains("\"channels\":2"));
    assert!(manifest.contains(&config.identity().to_hex()));
    assert!(manifest.contains("\"calibrated_acoustic_prediction\":false"));
    let receipt = config.receipt();
    assert_eq!(receipt.configuration_identity, config.identity());
    assert_eq!(receipt.authority, config.authority());
    assert_eq!(
        receipt.trajectory_identity,
        config.input().trajectory.identity()
    );
    for kind in [
        SoundOutputKind::Waveform,
        SoundOutputKind::SpectralMetric,
        SoundOutputKind::LevelMetric,
    ] {
        let declaration = config.declare_output(kind);
        assert_eq!(declaration.authority, config.authority());
        assert_eq!(
            declaration.assumptions,
            config.input().assumptions.as_slice()
        );
        assert!(!declaration.calibrated_acoustic_prediction);
    }
}

#[test]
fn calibrated_authority_requires_both_pressure_reference_and_receipt() {
    let mut input = informed_input();
    input.authority = SoundAuthority::Calibrated;
    assert_eq!(
        SoundSynthesisConfig::try_admit(input.clone()),
        Err(SoundSynthesisError::InvalidAmplitudeReference)
    );
    input.amplitude_reference = SoundAmplitudeReference::CalibratedPressure {
        pascal_rms_at_full_scale: 2.0,
        headroom_db: 12.0,
    };
    assert_eq!(
        SoundSynthesisConfig::try_admit(input.clone()),
        Err(SoundSynthesisError::MissingCalibration)
    );
    input.calibration = Some(calibration());
    let admitted = SoundSynthesisConfig::try_admit(input).unwrap();
    assert!(
        admitted
            .to_manifest_json()
            .contains("\"synthesis_class\":\"calibrated-structural-acoustics\"")
    );
    assert!(
        admitted
            .to_manifest_json()
            .contains("\"calibrated_acoustic_prediction\":true")
    );
}

#[test]
fn lower_tiers_cannot_smuggle_calibration_or_pressure_scale() {
    let mut input = informed_input();
    input.calibration = Some(calibration());
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::UnexpectedCalibration)
    );

    let mut input = informed_input();
    input.authority = SoundAuthority::Artistic;
    input.excitation_controls.clear();
    input.modes.clear();
    input.assumptions.clear();
    input.amplitude_reference = SoundAmplitudeReference::CalibratedPressure {
        pascal_rms_at_full_scale: 1.0,
        headroom_db: 6.0,
    };
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidAmplitudeReference)
    );
}

#[test]
fn exact_rational_av_endpoints_are_required() {
    let mut input = informed_input();
    input.audio_clock = CinematicClock::try_new(
        CinematicClockDomain::Audio,
        SOUND_MASTER_SAMPLE_RATE_HZ,
        1,
        0,
        479_999,
    )
    .unwrap();
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::AudioVideoTimelineMismatch)
    );

    let mut input = informed_input();
    input.audio_clock =
        CinematicClock::try_new(CinematicClockDomain::Audio, 44_100, 1, 0, 441_000).unwrap();
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidAudioClock)
    );

    let mut input = informed_input();
    input.video_clock =
        CinematicClock::try_new(CinematicClockDomain::Video, 48, 1, 0, 480).unwrap();
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidVideoClock)
    );
}

#[test]
fn stereo_and_camera_relative_orthonormal_listener_are_mandatory() {
    let mut input = informed_input();
    input.channel_layout = SoundChannelLayout::Mono;
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidChannelLayout)
    );

    let mut input = informed_input();
    input.listener.frame = ListenerFrame::World;
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidListenerFrame)
    );

    let mut input = informed_input();
    input.listener.up = [0.0, 2.0, 0.0];
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidListenerPose)
    );
}

#[test]
fn modal_frequency_damping_gain_and_order_are_checked() {
    let mutations: [fn(&mut SoundMode); 4] = [
        |mode: &mut SoundMode| mode.frequency_hz = 24_000.0,
        |mode: &mut SoundMode| mode.damping_ratio = 0.0,
        |mode: &mut SoundMode| mode.gain = f64::NAN,
        |mode: &mut SoundMode| mode.material_identity = ContentHash([0; 32]),
    ];
    for mutate in mutations {
        let mut input = informed_input();
        mutate(&mut input.modes[0]);
        assert_eq!(
            SoundSynthesisConfig::try_admit(input),
            Err(SoundSynthesisError::InvalidMode)
        );
    }
    let mut input = informed_input();
    input.modes[1].mode_id = 1;
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidMode)
    );
}

#[test]
fn informed_source_requires_channels_modes_and_complete_assumptions() {
    let mut input = informed_input();
    input.excitation_controls.clear();
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::MissingExcitationChannels)
    );

    let mut input = informed_input();
    input.excitation_controls[0].source_scale = 0.0;
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::NonCanonicalExcitationChannels)
    );

    let mut input = informed_input();
    input.modes.clear();
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::MissingModes)
    );

    let mut input = informed_input();
    input.assumptions.pop();
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidAssumptions)
    );
}

#[test]
fn terminal_censor_and_numerical_refusal_policies_are_distinct() {
    let mut input = informed_input();
    input.trajectory_disposition = SoundTrajectoryDisposition::HorizonCensored;
    assert!(SoundSynthesisConfig::try_admit(input).is_ok());

    let mut input = informed_input();
    input.trajectory_disposition = SoundTrajectoryDisposition::NumericalRefusal;
    assert_eq!(
        SoundSynthesisConfig::try_admit(input.clone()),
        Err(SoundSynthesisError::InvalidTerminalPolicy)
    );
    input.terminal_policy = SoundTerminalPolicy::Silence;
    assert!(SoundSynthesisConfig::try_admit(input).is_ok());

    let mut input = informed_input();
    input.terminal_policy = SoundTerminalPolicy::FadeAtLastAccepted {
        fade_sample_frames: 0,
    };
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidTerminalPolicy)
    );
}

#[test]
fn wrong_roles_and_missing_algorithm_versions_refuse() {
    let mut input = informed_input();
    input.schema_version += 1;
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::UnsupportedSchemaVersion(2))
    );

    let mut input = informed_input();
    input.room = component(CinematicComponentRole::Microphone);
    assert!(matches!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::WrongComponentRole {
            expected: CinematicComponentRole::Room,
            got: CinematicComponentRole::Microphone,
        })
    ));

    let mut input = informed_input();
    input.resampler_version = 0;
    assert_eq!(
        SoundSynthesisConfig::try_admit(input),
        Err(SoundSynthesisError::InvalidAlgorithmReference)
    );
}

#[test]
fn identity_is_deterministic_and_material_changes_invalidate_it() {
    let first = SoundSynthesisConfig::try_admit(informed_input()).unwrap();
    let replay = SoundSynthesisConfig::try_admit(informed_input()).unwrap();
    assert_eq!(first.identity(), replay.identity());

    let mutations: [fn(&mut SoundSynthesisInput); 6] = [
        |input| input.modes[0].frequency_hz += 1.0,
        |input| input.excitation_controls[0].source_scale *= 2.0,
        |input| input.listener.position_m[0] += 0.001,
        |input| input.filter_version += 1,
        |input| {
            input.room_response = SoundRoomResponse::DeclaredImpulseResponse { wet_mix: 0.09 };
        },
        |input| {
            input.terminal_policy = SoundTerminalPolicy::FadeAtLastAccepted {
                fade_sample_frames: 2_401,
            };
        },
    ];
    for mutate in mutations {
        let mut input = informed_input();
        mutate(&mut input);
        let changed = SoundSynthesisConfig::try_admit(input).unwrap();
        assert_ne!(first.identity(), changed.identity());
    }
}
