//! G0/G3 checks for deterministic event-aware Euler timeline reconstruction.

use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::coupled_runner::{ChannelOwnership, ContactTransitionKind};
use fs_euler_disc_e2e::{
    DeclaredDiscontinuityKind, DeclaredTimelineDiscontinuity, DerivedEulerQois,
    EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EventEvaluationSide, ExposureEventPolicy,
    RenderBaseModeState, RenderChannelAvailability, RenderContactBranch, RenderContactGeometry,
    RenderContactTransition, RenderMassProperties, RenderSampleDisposition, RenderSupportFeature,
    RenderTerminalEvent, RenderTrajectory, RenderTrajectoryAuthority, RenderTrajectoryMetadata,
    RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame, TimelineEvent,
    TimelineResampler, TimelineResamplingError, TimelineSampleSource,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn identity(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.test.timeline-resampling.v1",
        label.as_bytes(),
    )
}

fn mass() -> MassProperties {
    MassProperties::new(1.0, Vec3::ZERO, Vec3::new(0.2, 0.2, 0.2)).unwrap()
}

fn assert_close(left: f64, right: f64) {
    let scale = 1.0_f64.max(left.abs()).max(right.abs());
    assert!((left - right).abs() <= 2.0e-12 * scale, "{left} != {right}");
}

fn state(time_s: f64, angle_rad: f64) -> RigidBodyState {
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), angle_rad).unwrap();
    RigidBodyState::new(
        Pose::new(Vec3::new(2.0 * time_s, -time_s, 0.5 * time_s), orientation).unwrap(),
        Vec3::new(2.0, -1.0, 0.5),
        Vec3::new(0.1, 0.0, 0.0),
    )
    .unwrap()
}

fn tilted_spin(spin_rad: f64) -> UnitQuaternion {
    let tilt = 0.35_f64;
    let (tilt_sine, tilt_cosine) = (0.5 * tilt).sin_cos();
    let (spin_sine, spin_cosine) = (0.5 * spin_rad).sin_cos();
    UnitQuaternion::new(
        tilt_cosine * spin_cosine,
        tilt_sine * spin_sine,
        tilt_sine * spin_cosine,
        tilt_cosine * spin_sine,
    )
    .unwrap()
}

fn replace_orientation(input: &mut RenderTrajectorySampleInput, orientation: UnitQuaternion) {
    input.orientation_body_to_world = orientation.components();
    input.symmetry_axis_world = orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let state = RigidBodyState::new(
        Pose::new(input.center_of_mass_world_m, orientation).unwrap(),
        input.linear_momentum_world_kg_m_per_s,
        input.angular_momentum_body_kg_m2_per_s,
    )
    .unwrap();
    input.qois = DerivedEulerQois::from_state(state, mass(), 0.0).unwrap();
}

fn input(
    time_s: f64,
    angle_rad: f64,
    branch: RenderContactBranch,
    disposition: RenderSampleDisposition,
) -> RenderTrajectorySampleInput {
    let state = state(time_s, angle_rad);
    let orientation = state.pose().orientation();
    RenderTrajectorySampleInput {
        interval_start_time_s: 0.0,
        time_s,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m: state.pose().position_world(),
        orientation_body_to_world: orientation.components(),
        linear_momentum_world_kg_m_per_s: state.linear_momentum_world(),
        angular_momentum_body_kg_m2_per_s: state.angular_momentum_body(),
        symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        contact_branch: branch,
        contact_geometry: (branch == RenderContactBranch::Closed).then_some(
            RenderContactGeometry {
                point_world_m: Vec3::new(2.0 * time_s, -time_s, 0.0),
                normal_world: Vec3::new(0.0, 0.0, 1.0),
                support_feature: RenderSupportFeature::CylinderRim,
            },
        ),
        signed_gap_m: if branch == RenderContactBranch::Closed {
            0.0
        } else {
            1.0e-3
        },
        interval_normal_force_n: if branch == RenderContactBranch::Closed {
            1.0
        } else {
            0.0
        },
        contact_transitions: Vec::new(),
        base_mode: Some(RenderBaseModeState {
            displacement_m: 2.0 * time_s,
            velocity_m_per_s: 2.0,
        }),
        channels: ChannelOwnership::default(),
        mechanical_energy_j: 1.0,
        energy_defect_j: 0.0,
        qois: DerivedEulerQois::from_state(state, mass(), 0.0).unwrap(),
        disposition,
        terminal_event: None,
    }
}

fn metadata(initial: &RenderTrajectorySampleInput) -> RenderTrajectoryMetadata {
    let orientation = UnitQuaternion::new(
        initial.orientation_body_to_world[0],
        initial.orientation_body_to_world[1],
        initial.orientation_body_to_world[2],
        initial.orientation_body_to_world[3],
    )
    .unwrap();
    let initial_state = RigidBodyState::new(
        Pose::new(initial.center_of_mass_world_m, orientation).unwrap(),
        initial.linear_momentum_world_kg_m_per_s,
        initial.angular_momentum_body_kg_m2_per_s,
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
        initial_base_mode: initial.base_mode.unwrap(),
        base_model_identity: identity("base"),
        model_identity: identity("model"),
        channel_availability: RenderChannelAvailability::ALL_AVAILABLE,
        configuration_identity: identity("configuration"),
        configuration_fingerprint: 0x7265_7361_6d70_6c65,
        timestep_s: 1.0,
        producer_version: "timeline-resampling-test-v1".into(),
        applicability: "deterministic visualization reconstruction only".into(),
        no_claims: vec!["does not add mechanical resolution".into()],
        authority: RenderTrajectoryAuthority::SimulationEvidence,
    }
}

fn analytic_trajectory() -> RenderTrajectory {
    let first = input(
        0.0,
        0.2,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let mut second = input(
        1.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    second.interval_start_time_s = first.time_s;
    RenderTrajectory::try_new(metadata(&first), vec![first, second]).unwrap()
}

fn event_trajectory() -> RenderTrajectory {
    let first = input(
        0.0,
        0.2,
        RenderContactBranch::Closed,
        RenderSampleDisposition::Continue,
    );
    let mut second = input(
        1.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    second.interval_start_time_s = first.time_s;
    second.contact_transitions.push(RenderContactTransition {
        kind: ContactTransitionKind::Opening,
        time_s: 0.5,
        bracket_start_s: 0.49,
        bracket_end_s: 0.51,
    });
    RenderTrajectory::try_new(metadata(&first), vec![first, second]).unwrap()
}

fn terminal_trajectory() -> RenderTrajectory {
    let first = input(
        0.0,
        0.2,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let mut second = input(
        1.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::TerminalInclination,
    );
    second.interval_start_time_s = first.time_s;
    second.terminal_event = Some(RenderTerminalEvent {
        time_s: 1.0,
        bracket_start_s: 0.99,
        bracket_end_s: 1.01,
    });
    RenderTrajectory::try_new(metadata(&first), vec![first, second]).unwrap()
}

#[test]
fn exact_endpoints_and_analytic_constant_motion_are_reproduced() {
    let trajectory = analytic_trajectory();
    let samples = TimelineResampler::new(&trajectory)
        .resample(&[0.0, 0.5, 1.0], EventEvaluationSide::RightLimit)
        .unwrap();

    assert_eq!(samples[0].state, trajectory.samples()[0].state());
    assert_eq!(samples[2].state, trajectory.samples()[1].state());
    assert_eq!(
        samples[0].source,
        TimelineSampleSource::ExactSample { index: 0 }
    );
    assert_eq!(
        samples[2].source,
        TimelineSampleSource::ExactSample { index: 1 }
    );
    let midpoint = &samples[1];
    assert_close(midpoint.state.pose().position_world().x, 1.0);
    assert_close(midpoint.state.pose().position_world().y, -0.5);
    assert_close(midpoint.base_mode.displacement_m, 1.0);
    assert_close(midpoint.base_mode.velocity_m_per_s, 2.0);
    let axis = midpoint
        .state
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    assert_close(axis.y, -0.6_f64.sin());
    assert_close(axis.z, 0.6_f64.cos());
}

#[test]
fn q_double_cover_tiny_and_near_pi_rotations_take_the_shortest_arc() {
    for (start, end) in [(0.2, 0.2 + 1.0e-13), (0.0, core::f64::consts::PI - 1.0e-10)] {
        let mut first = input(
            0.0,
            0.35,
            RenderContactBranch::Open,
            RenderSampleDisposition::Continue,
        );
        let mut second = input(
            1.0,
            0.35,
            RenderContactBranch::Open,
            RenderSampleDisposition::HorizonCensored,
        );
        replace_orientation(&mut first, tilted_spin(start));
        replace_orientation(&mut second, tilted_spin(end));
        for component in &mut second.orientation_body_to_world {
            *component = -*component;
        }
        let trajectory = RenderTrajectory::try_new(metadata(&first), vec![first, second]).unwrap();
        let midpoint = TimelineResampler::new(&trajectory)
            .resample(&[0.5], EventEvaluationSide::RightLimit)
            .unwrap()
            .pop()
            .unwrap();
        let expected = tilted_spin(0.5 * (start + end));
        let observed = midpoint.state.pose().orientation().components();
        let expected = expected.components();
        for index in 0..4 {
            assert_close(observed[index], expected[index]);
        }
    }
}

#[test]
fn contact_events_use_explicit_one_sided_discrete_semantics() {
    let trajectory = event_trajectory();
    let resampler = TimelineResampler::new(&trajectory);
    let left = resampler
        .resample(&[0.5], EventEvaluationSide::LeftLimit)
        .unwrap();
    let right = resampler
        .resample(&[0.5], EventEvaluationSide::RightLimit)
        .unwrap();
    assert_eq!(left[0].state, right[0].state);
    assert_eq!(left[0].contact_branch, RenderContactBranch::Closed);
    assert_eq!(right[0].contact_branch, RenderContactBranch::Open);
    assert_eq!(right[0].events_at_query.len(), 1);
    assert!(matches!(
        right[0].events_at_query[0],
        TimelineEvent::Contact(RenderContactTransition {
            bracket_start_s: 0.49,
            bracket_end_s: 0.51,
            ..
        })
    ));
}

#[test]
fn terminal_root_has_explicit_left_and_right_dispositions() {
    let trajectory = terminal_trajectory();
    let resampler = TimelineResampler::new(&trajectory);
    let left = resampler
        .resample(&[1.0], EventEvaluationSide::LeftLimit)
        .unwrap();
    let right = resampler
        .resample(&[1.0], EventEvaluationSide::RightLimit)
        .unwrap();
    assert_eq!(left[0].state, right[0].state);
    assert_eq!(left[0].disposition, RenderSampleDisposition::Continue);
    assert_eq!(
        right[0].disposition,
        RenderSampleDisposition::TerminalInclination
    );
    assert!(matches!(
        right[0].events_at_query.as_slice(),
        [TimelineEvent::TerminalInclination(RenderTerminalEvent {
            bracket_start_s: 0.99,
            bracket_end_s: 1.01,
            ..
        })]
    ));
}

#[test]
fn shutter_intervals_subdivide_or_refuse_at_events() {
    let trajectory = event_trajectory();
    let resampler = TimelineResampler::new(&trajectory);
    let partition = resampler
        .partition_exposure(0.25, 0.75, ExposureEventPolicy::Subdivide)
        .unwrap();
    assert_eq!(partition.interior_events.len(), 1);
    assert_eq!(partition.segments.len(), 2);
    assert_eq!(partition.segments[0].start_s.to_bits(), 0.25_f64.to_bits());
    assert_eq!(partition.segments[0].end_s.to_bits(), 0.5_f64.to_bits());
    assert_eq!(partition.segments[1].start_s.to_bits(), 0.5_f64.to_bits());
    assert_eq!(partition.segments[1].end_s.to_bits(), 0.75_f64.to_bits());
    assert_eq!(
        resampler.partition_exposure(0.25, 0.75, ExposureEventPolicy::Refuse),
        Err(TimelineResamplingError::ExposureSpansEvent)
    );
}

#[test]
fn query_and_declared_seam_validation_refuse_ambiguous_timelines() {
    let trajectory = analytic_trajectory();
    let resampler = TimelineResampler::new(&trajectory);
    assert_eq!(
        resampler.resample(&[], EventEvaluationSide::RightLimit),
        Err(TimelineResamplingError::EmptyQueries)
    );
    assert_eq!(
        resampler.resample(&[0.5, 0.5], EventEvaluationSide::RightLimit),
        Err(TimelineResamplingError::NonIncreasingQuery(1))
    );
    assert_eq!(
        resampler.resample(&[-0.1], EventEvaluationSide::RightLimit),
        Err(TimelineResamplingError::QueryOutOfRange {
            index: 0,
            time_s: -0.1,
        })
    );
    assert_eq!(
        resampler.resample(&[f64::NAN], EventEvaluationSide::RightLimit),
        Err(TimelineResamplingError::NonFiniteQuery(0))
    );
    let negative_zero = resampler
        .resample(&[-0.0], EventEvaluationSide::RightLimit)
        .unwrap();
    assert_eq!(
        negative_zero[0].source,
        TimelineSampleSource::ExactSample { index: 0 }
    );
    assert_eq!(
        TimelineResampler::with_declared_discontinuities(
            &trajectory,
            vec![DeclaredTimelineDiscontinuity {
                time_s: 0.7,
                kind: DeclaredDiscontinuityKind::ContinuationSeam,
            }],
        )
        .err(),
        Some(TimelineResamplingError::InvalidDeclaredDiscontinuity(0))
    );
    assert_eq!(
        TimelineResampler::with_declared_discontinuities(
            &trajectory,
            vec![
                DeclaredTimelineDiscontinuity {
                    time_s: 1.0,
                    kind: DeclaredDiscontinuityKind::ContinuationSeam,
                },
                DeclaredTimelineDiscontinuity {
                    time_s: 0.0,
                    kind: DeclaredDiscontinuityKind::ProducerDeclared,
                },
            ],
        )
        .err(),
        Some(TimelineResamplingError::InvalidDeclaredDiscontinuity(1))
    );
}

#[test]
fn time_translation_frame_rate_changes_and_rigid_translation_are_equivariant() {
    let trajectory = analytic_trajectory();
    let coarse = TimelineResampler::new(&trajectory)
        .resample(&[0.25, 0.5, 0.75], EventEvaluationSide::RightLimit)
        .unwrap();
    let fine = TimelineResampler::new(&trajectory)
        .resample(
            &[0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875],
            EventEvaluationSide::RightLimit,
        )
        .unwrap();
    assert_eq!(coarse[0], fine[1]);
    assert_eq!(coarse[1], fine[3]);
    assert_eq!(coarse[2], fine[5]);

    let offset = Vec3::new(10.0, -4.0, 3.0);
    let mut first = input(
        10.0,
        0.2,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let mut second = input(
        11.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    first.center_of_mass_world_m = first.center_of_mass_world_m.add(offset);
    second.center_of_mass_world_m = second.center_of_mass_world_m.add(offset);
    let shifted = RenderTrajectory::try_new(metadata(&first), vec![first, second]).unwrap();
    let shifted_midpoint = TimelineResampler::new(&shifted)
        .resample(&[10.5], EventEvaluationSide::RightLimit)
        .unwrap();
    let shifted_position = shifted_midpoint[0].state.pose().position_world();
    assert_close(shifted_position.x, 21.0 + offset.x);
    assert_close(shifted_position.y, -10.5 + offset.y);
    assert_close(shifted_position.z, 5.25 + offset.z);

    let mut scaled_first = input(
        0.0,
        0.2,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let mut scaled_last = input(
        1.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    let time_scale = 2.5;
    let time_offset = 7.0;
    scaled_first.time_s = time_offset;
    scaled_last.time_s = time_offset + time_scale;
    scaled_first.linear_momentum_world_kg_m_per_s = scaled_first
        .linear_momentum_world_kg_m_per_s
        .scale(time_scale.recip());
    scaled_last.linear_momentum_world_kg_m_per_s = scaled_last
        .linear_momentum_world_kg_m_per_s
        .scale(time_scale.recip());
    scaled_first.base_mode.as_mut().unwrap().velocity_m_per_s /= time_scale;
    scaled_last.base_mode.as_mut().unwrap().velocity_m_per_s /= time_scale;
    let mut scaled_metadata = metadata(&scaled_first);
    scaled_metadata.timestep_s = time_scale;
    let scaled =
        RenderTrajectory::try_new(scaled_metadata, vec![scaled_first, scaled_last]).unwrap();
    let scaled_midpoint = TimelineResampler::new(&scaled)
        .resample(
            &[time_offset + 0.5 * time_scale],
            EventEvaluationSide::RightLimit,
        )
        .unwrap();
    let canonical_midpoint = &coarse[1];
    let scaled_midpoint = &scaled_midpoint[0];
    assert_eq!(
        scaled_midpoint.state.pose().position_world(),
        canonical_midpoint.state.pose().position_world()
    );
    assert_eq!(
        scaled_midpoint.state.pose().orientation(),
        canonical_midpoint.state.pose().orientation()
    );
    assert_close(
        scaled_midpoint.base_mode.displacement_m,
        canonical_midpoint.base_mode.displacement_m,
    );
    assert_close(
        scaled_midpoint.base_mode.velocity_m_per_s * time_scale,
        canonical_midpoint.base_mode.velocity_m_per_s,
    );
}

#[test]
fn unrepresentable_finite_time_interval_refuses_instead_of_emitting_nan() {
    let mut first = input(
        0.0,
        0.2,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let mut second = input(
        1.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    first.linear_momentum_world_kg_m_per_s = Vec3::ZERO;
    second.linear_momentum_world_kg_m_per_s = Vec3::ZERO;
    first.time_s = 0.0;
    second.time_s = f64::MAX;
    let trajectory = RenderTrajectory::try_new(metadata(&first), vec![first, second]).unwrap();
    assert_eq!(
        TimelineResampler::new(&trajectory)
            .resample(&[0.5 * f64::MAX], EventEvaluationSide::RightLimit)
            .unwrap_err(),
        TimelineResamplingError::InvalidReconstruction(
            "base interpolation produced non-finite state".into()
        )
    );
}

#[test]
fn declared_continuation_seams_partition_exposure_without_changing_pose() {
    let first = input(
        0.0,
        0.2,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let middle = input(
        0.4,
        0.52,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let mut last = input(
        1.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    last.interval_start_time_s = middle.time_s;
    let trajectory =
        RenderTrajectory::try_new(metadata(&first), vec![first, middle, last]).unwrap();
    let resampler = TimelineResampler::with_declared_discontinuities(
        &trajectory,
        vec![DeclaredTimelineDiscontinuity {
            time_s: 0.4,
            kind: DeclaredDiscontinuityKind::ContinuationSeam,
        }],
    )
    .unwrap();
    let sample = resampler
        .resample(&[0.4], EventEvaluationSide::RightLimit)
        .unwrap();
    assert!(matches!(
        sample[0].events_at_query.as_slice(),
        [TimelineEvent::Declared(DeclaredTimelineDiscontinuity {
            kind: DeclaredDiscontinuityKind::ContinuationSeam,
            ..
        })]
    ));
    let partition = resampler
        .partition_exposure(0.1, 0.9, ExposureEventPolicy::Subdivide)
        .unwrap();
    assert_eq!(partition.segments.len(), 2);
}
