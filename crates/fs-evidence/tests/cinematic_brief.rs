//! G0/G3 tests for the reference shot list, exact clocks, and alias diagnostics.

use core::f64::consts::TAU;

use fs_evidence::cinematic_brief::{
    ApparentDirection, ApparentRotationAssessment, ApparentRotationIntent, AudioPerspective,
    BackgroundPreset, BriefCameraKeyframe, BriefInterpolationDomain, BriefOptics, BriefVec3,
    CensoredTrajectoryPolicy, CinematicBrief, CinematicBriefError, CinematicBriefInput,
    DiscMaterialPreset, ExposureIntent, ExposureLegibility, FocusTarget, FrameRange,
    FrameTimeConvention, LightingPreset, ReferenceShotInput, ReferenceShotRole,
    SafeAreaInsetsPermille, ShotTimeMapping, ShotTransition, ShutterBoundaryPolicy, SpinCue,
    VisualizationHoldLabel, assess_apparent_rotation, assess_apparent_rotation_at_rate,
};

fn keyframe(frame: u32) -> BriefCameraKeyframe {
    BriefCameraKeyframe {
        frame,
        eye_m: BriefVec3 {
            x: 0.2,
            y: -0.2,
            z: 0.1,
        },
        target_m: BriefVec3 {
            x: 0.0,
            y: 0.0,
            z: 0.01,
        },
        up: BriefVec3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        },
        focus_distance_m: 0.25,
        focus_target: FocusTarget::DiscCenter,
    }
}

fn shot(start: u32, end: u32, start_tick: u64, end_tick: u64) -> ReferenceShotInput {
    ReferenceShotInput {
        role: ReferenceShotRole::EstablishingProduct,
        frames: FrameRange::try_new(start, end).expect("nonempty test range"),
        camera: vec![keyframe(start), keyframe(end - 1)],
        optics: BriefOptics {
            focal_length_um: 50_000,
            f_number_milli: 4_000,
            shutter_open_microframes: -180_000,
            shutter_close_microframes: 180_000,
            exposure: ExposureIntent::ProtectMetalHighlights,
        },
        lighting: LightingPreset::StudioSoftboxRim,
        disc_material: DiscMaterialPreset::BrushedTungsten,
        background: BackgroundPreset::NeutralBlackSweep,
        glass_visible: true,
        base_visible: true,
        audio_perspective: AudioPerspective::StudioObserver,
        transition: ShotTransition::HardCut,
        time_mapping: ShotTimeMapping::PhysicalLinear {
            start_tick,
            end_tick_exclusive: end_tick,
        },
        apparent_rotation_intent: ApparentRotationIntent::ForwardCueAndInclinationReadable,
    }
}

fn input() -> CinematicBriefInput {
    CinematicBriefInput {
        total_frames: 240,
        total_audio_sample_frames: 480_000,
        simulation_ticks_per_second: 1_000,
        trajectory_start_tick: 0,
        trajectory_end_tick_exclusive: 10_000,
        shots: vec![shot(0, 240, 0, 10_000)],
        safe_areas: SafeAreaInsetsPermille {
            action: 50,
            title: 100,
        },
        spin_cue: SpinCue {
            brushing_marking_frequency: 2,
            engraved_radial_mark: false,
        },
        censored_policy: CensoredTrajectoryPolicy::LabeledHoldAndAudioTaper {
            maximum_hold_frames: 12,
        },
        frame_time_convention: FrameTimeConvention::IntegerFrameCentersHalfOpenMaster,
        interpolation_domain: BriefInterpolationDomain::StudioFrameRigidPoseAndScalarOptics,
        shutter_boundary_policy: ShutterBoundaryPolicy::ClipAtCutsAndMaster,
        audio_lead_samples: 0,
        muted_review_required: true,
        essential_context_in_audio_only: false,
    }
}

#[test]
fn brief_identity_is_deterministic_and_binds_every_retained_semantic() {
    let first = CinematicBrief::try_new(input()).expect("admitted brief");
    let replay = CinematicBrief::try_new(input()).expect("replayed brief");
    assert_eq!(first.canonical_bytes(), replay.canonical_bytes());
    assert_eq!(first.identity(), replay.identity());
    assert!(first.identity().as_bytes().iter().any(|byte| *byte != 0));

    let mut changed_input = input();
    changed_input.spin_cue.engraved_radial_mark = true;
    let changed = CinematicBrief::try_new(changed_input).expect("changed admitted brief");
    assert_ne!(first.identity(), changed.identity());

    let mut camera_changed_input = input();
    camera_changed_input.shots[0].camera[0].focus_distance_m = 0.3;
    let camera_changed =
        CinematicBrief::try_new(camera_changed_input).expect("camera-changed admitted brief");
    assert_ne!(first.identity(), camera_changed.identity());

    let mut signed_zero_input = input();
    signed_zero_input.shots[0].camera[0].up.x = -0.0;
    let signed_zero =
        CinematicBrief::try_new(signed_zero_input).expect("signed-zero-equivalent brief");
    assert_eq!(
        first.identity(),
        signed_zero.identity(),
        "canonical identity must not distinguish IEEE signed zero"
    );
}

#[test]
fn brief_admission_bounds_identity_table_cardinality() {
    let mut too_many_shots = input();
    let template = too_many_shots.shots[0].clone();
    too_many_shots.shots.resize(241, template);
    assert!(matches!(
        CinematicBrief::try_new(too_many_shots),
        Err(CinematicBriefError::TooManyShots {
            maximum: 240,
            got: 241
        })
    ));

    let mut too_many_keyframes = input();
    let template = too_many_keyframes.shots[0].camera[0];
    too_many_keyframes.shots[0].camera.resize(241, template);
    assert!(matches!(
        CinematicBrief::try_new(too_many_keyframes),
        Err(CinematicBriefError::TooManyCameraKeyframes {
            maximum: 240,
            got: 241
        })
    ));
}

#[test]
fn reference_brief_freezes_four_contiguous_readable_shots() {
    let brief = CinematicBrief::euler_disc_v1().expect("reference brief must admit");
    assert_eq!(brief.total_frames(), 240);
    assert_eq!(brief.total_audio_sample_frames(), 480_000);
    assert_eq!(brief.simulation_ticks_per_second(), 1_000);
    assert_eq!(brief.shots().len(), 4);
    assert_eq!(
        brief
            .shots()
            .iter()
            .map(fs_evidence::cinematic_brief::ReferenceShot::role)
            .collect::<Vec<_>>(),
        vec![
            ReferenceShotRole::EstablishingProduct,
            ReferenceShotRole::InclinationAndContactOrbit,
            ReferenceShotRole::MacroPrecession,
            ReferenceShotRole::TerminalCloseUp,
        ]
    );
    assert_eq!(
        brief
            .shots()
            .iter()
            .map(|shot| (shot.frames().start(), shot.frames().end_exclusive()))
            .collect::<Vec<_>>(),
        vec![(0, 60), (60, 120), (120, 192), (192, 240)]
    );
    assert!(
        brief
            .shots()
            .iter()
            .all(fs_evidence::cinematic_brief::ReferenceShot::glass_visible)
    );
    assert!(
        brief
            .shots()
            .iter()
            .all(|shot| shot.background() == BackgroundPreset::NeutralBlackSweep)
    );
    assert_eq!(brief.audio_lead_samples(), 0);
    assert_eq!(
        brief
            .shots()
            .iter()
            .map(fs_evidence::cinematic_brief::ReferenceShot::apparent_rotation_intent)
            .collect::<Vec<_>>(),
        vec![
            ApparentRotationIntent::ForwardCueAndInclinationReadable,
            ApparentRotationIntent::ContactOrbitPrimaryNoFalseReversal,
            ApparentRotationIntent::WobblePrimaryNoFrozenCue,
            ApparentRotationIntent::TerminalSlowdownNoAliasReversal,
        ]
    );
    assert_eq!(
        brief.frame_time_convention(),
        FrameTimeConvention::IntegerFrameCentersHalfOpenMaster
    );
    assert_eq!(
        brief.shutter_boundary_policy(),
        ShutterBoundaryPolicy::ClipAtCutsAndMaster
    );
}

#[test]
fn golden_frame_sample_and_simulation_timestamps_are_exact() {
    let brief = CinematicBrief::euler_disc_v1().expect("reference brief must admit");
    for frame in 0..=240 {
        assert_eq!(
            brief.audio_sample_for_frame(frame),
            Ok(u64::from(frame) * 2_000)
        );
    }
    assert_eq!(
        brief
            .simulation_tick_for_frame(0)
            .expect("mapped")
            .numerator,
        0
    );
    let at_cut = brief.simulation_tick_for_frame(60).expect("mapped");
    assert_eq!(at_cut.numerator, 2_500 * 60);
    assert_eq!(at_cut.denominator, 60);
    let last = brief.simulation_tick_for_frame(239).expect("mapped");
    assert_eq!(last.numerator, 8_000 * 48 + 2_000 * 47);
    assert_eq!(last.denominator, 48);
    assert_eq!(
        brief.audio_sample_for_frame(241),
        Err(CinematicBriefError::FrameOutsideMaster(241))
    );

    for sample in 0..brief.total_audio_sample_frames() {
        assert!(brief.simulation_tick_for_audio_sample(sample).is_ok());
    }
    assert_eq!(
        brief.simulation_tick_for_audio_sample(480_000),
        Err(CinematicBriefError::AudioSampleOutsideMaster(480_000))
    );
    let audio_at_cut = brief
        .simulation_tick_for_audio_sample(120_000)
        .expect("cut sample maps to entering shot");
    assert_eq!(audio_at_cut.numerator, 2_500 * 120_000);
    assert_eq!(audio_at_cut.denominator, 120_000);
}

#[test]
fn shutter_support_clips_at_master_and_cut_boundaries() {
    let brief = CinematicBrief::euler_disc_v1().expect("reference brief must admit");
    let first = brief.effective_shutter_window(0).expect("first frame");
    assert_eq!(first.start_microframes, 0);
    assert!(first.clipped_at_boundary);

    let entering_cut = brief
        .effective_shutter_window(60)
        .expect("first frame after cut");
    assert_eq!(entering_cut.start_microframes, 60_000_000);
    assert!(entering_cut.clipped_at_boundary);

    let interior = brief.effective_shutter_window(61).expect("interior frame");
    assert_eq!(interior.start_microframes, 60_820_000);
    assert_eq!(interior.end_microframes, 61_180_000);
    assert!(!interior.clipped_at_boundary);
}

#[test]
fn storyboard_proxy_covers_every_frame_and_is_useful_when_muted() {
    let brief = CinematicBrief::euler_disc_v1().expect("reference brief must admit");
    let proxy = brief.storyboard_proxy().expect("whole brief must map");
    assert_eq!(proxy.len(), 240);
    assert_eq!(proxy.first().expect("first").frame, 0);
    assert_eq!(proxy.last().expect("last").frame, 239);
    assert!(proxy.iter().all(|frame| frame.audio_muted_for_review));
    assert!(
        brief
            .muted_review_manifest_json()
            .contains("\"audio_muted\":true")
    );
    assert!(
        brief
            .muted_review_manifest_json()
            .contains("\"essential_context_in_audio_only\":false")
    );
}

#[test]
fn shot_gaps_overlaps_and_uncovered_tail_refuse() {
    let mut candidate = input();
    candidate.shots = vec![shot(0, 100, 0, 4_000), shot(101, 240, 4_000, 10_000)];
    assert!(matches!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::ShotGapOrOverlap {
            expected_start: 100,
            got: 101
        })
    ));

    let mut candidate = input();
    candidate.shots = vec![shot(0, 100, 0, 4_000), shot(99, 240, 4_000, 10_000)];
    assert!(matches!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::ShotGapOrOverlap { .. })
    ));

    let mut candidate = input();
    candidate.shots = vec![shot(0, 239, 0, 10_000)];
    assert!(matches!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::ShotGapOrOverlap { .. })
    ));
}

#[test]
fn camera_singularities_missing_focus_and_keyframe_order_refuse() {
    let mut candidate = input();
    candidate.shots[0].camera[0].target_m = candidate.shots[0].camera[0].eye_m;
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::CameraSingularity)
    );

    let mut candidate = input();
    candidate.shots[0].camera[0].focus_distance_m = 0.0;
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::InvalidCameraValue)
    );

    let mut candidate = input();
    candidate.shots[0].camera[1].frame = 0;
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::InvalidCameraKeyframeOrder)
    );
}

#[test]
fn av_duration_trajectory_shutter_and_audio_only_meaning_fail_closed() {
    let mut candidate = input();
    candidate.total_audio_sample_frames -= 1;
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::InvalidMasterTimeline)
    );

    let mut candidate = input();
    candidate.shots[0].time_mapping = ShotTimeMapping::PhysicalLinear {
        start_tick: 0,
        end_tick_exclusive: 10_001,
    };
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::TimeOutsideTrajectory)
    );

    let mut candidate = input();
    candidate.shots[0].optics.shutter_close_microframes = 500_001;
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::InvalidOptics)
    );

    let mut candidate = input();
    candidate.essential_context_in_audio_only = true;
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::AudioOnlyMeaning)
    );
}

#[test]
fn explicit_censored_hold_is_mapped_and_labeled_visualization_only() {
    let mut candidate = input();
    let mut hold = shot(228, 240, 9_999, 10_000);
    hold.time_mapping = ShotTimeMapping::VisualizationHold {
        source_tick: 9_999,
        label: VisualizationHoldLabel::CensoredTrajectory,
    };
    candidate.shots = vec![shot(0, 228, 0, 9_999), hold];
    let brief = CinematicBrief::try_new(candidate).expect("declared in-range hold must admit");
    let mapped = brief.simulation_tick_for_frame(239).expect("hold maps");
    assert_eq!(mapped.numerator, 9_999);
    assert_eq!(mapped.denominator, 1);
    assert!(mapped.visualization_only);

    let mut candidate = input();
    candidate.shots[0].time_mapping = ShotTimeMapping::VisualizationHold {
        source_tick: 9_998,
        label: VisualizationHoldLabel::CensoredTrajectory,
    };
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::HoldIsNotTerminalState)
    );
}

#[test]
fn censored_hold_bounds_order_and_audio_offset_fail_closed() {
    let mut candidate = input();
    candidate.censored_policy = CensoredTrajectoryPolicy::LabeledHoldAndAudioTaper {
        maximum_hold_frames: 0,
    };
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::InvalidCensoredPolicy)
    );

    let mut candidate = input();
    candidate.shots[0].time_mapping = ShotTimeMapping::VisualizationHold {
        source_tick: 9_999,
        label: VisualizationHoldLabel::CensoredTrajectory,
    };
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::CensoredHoldTooLong {
            maximum: 12,
            got: 240,
        })
    );

    let mut candidate = input();
    let mut hold = shot(0, 12, 9_999, 10_000);
    hold.time_mapping = ShotTimeMapping::VisualizationHold {
        source_tick: 9_999,
        label: VisualizationHoldLabel::CensoredTrajectory,
    };
    candidate.shots = vec![hold, shot(12, 240, 0, 9_999)];
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::PhysicalShotAfterCensoredHold)
    );

    let mut candidate = input();
    candidate.audio_lead_samples = 1;
    assert_eq!(
        CinematicBrief::try_new(candidate),
        Err(CinematicBriefError::UnsupportedAudioLeadLag)
    );
}

#[test]
fn apparent_rotation_detects_reversal_frozen_markings_and_shutter_smear() {
    let false_reverse =
        assess_apparent_rotation(0.75 * TAU * 24.0 / 2.0, 2, 180_000).expect("finite alias case");
    assert_eq!(false_reverse.direction, ApparentDirection::FalseReversal);
    assert!((false_reverse.apparent_cycles_per_frame + 0.25).abs() < 1.0e-12);

    let frozen = assess_apparent_rotation(TAU * 24.0 / 2.0, 2, 1_000_000)
        .expect("one marking cycle per frame");
    assert_eq!(
        frozen,
        ApparentRotationAssessment {
            direction: ApparentDirection::FrozenMarking,
            exposure_legibility: ExposureLegibility::ShutterSmear,
            apparent_cycles_per_frame: 0.0,
        }
    );
    let flicker = assess_apparent_rotation(0.49 * TAU * 24.0, 1, 10_000).expect("near-Nyquist cue");
    assert_eq!(
        flicker.exposure_legibility,
        ExposureLegibility::ObjectionableFlickerRisk
    );

    let at_24 =
        assess_apparent_rotation_at_rate(TAU * 24.0, 1, 10_000, 24, 1).expect("24 fps sample");
    let at_30 =
        assess_apparent_rotation_at_rate(TAU * 24.0, 1, 10_000, 30, 1).expect("30 fps sample");
    assert_eq!(at_24.direction, ApparentDirection::FrozenMarking);
    assert_eq!(at_30.direction, ApparentDirection::FalseReversal);
    assert_eq!(
        assess_apparent_rotation(f64::NAN, 2, 100_000),
        Err(CinematicBriefError::InvalidApparentRotationInput)
    );
}
