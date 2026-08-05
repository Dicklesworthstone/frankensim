//! G0/G3 controls: frame/unit semantics, conservative boxcar mitigation, replay,
//! cancellation, and explicit no-data behavior.

use core::{f64::consts::FRAC_PI_2, num::NonZeroUsize};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::{
    AudioControlFilter, ChannelControl, ContactEventMeasure, ControlStreamError, DerivedEulerQois,
    EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerControlStream, RenderBaseFrame,
    RenderBaseModeState, RenderChannelAvailability, RenderContactBranch, RenderContactGeometry,
    RenderContactTransition, RenderMassProperties, RenderSampleDisposition, RenderSupportFeature,
    RenderTrajectory, RenderTrajectoryAuthority, RenderTrajectoryError, RenderTrajectoryMetadata,
    RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame,
    coupled_runner::{ChannelOwnership, ChannelWrench, ContactTransitionKind},
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
                seed: 0x434f_4e54_524f_4c53,
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
    hash_domain("org.frankensim.test.euler-controls.v1", label.as_bytes())
}

fn mass() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap()
}

fn state_from(input: &RenderTrajectorySampleInput) -> RigidBodyState {
    let q = input.orientation_body_to_world;
    RigidBodyState::new(
        Pose::new(
            input.center_of_mass_world_m,
            UnitQuaternion::new(q[0], q[1], q[2], q[3]).unwrap(),
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

fn compose(left: UnitQuaternion, right: UnitQuaternion) -> UnitQuaternion {
    let [lw, lx, ly, lz] = left.components();
    let [rw, rx, ry, rz] = right.components();
    UnitQuaternion::new(
        lw.mul_add(rw, -lx.mul_add(rx, ly.mul_add(ry, lz * rz))),
        lw.mul_add(rx, lx.mul_add(rw, ly.mul_add(rz, -(lz * ry)))),
        lw.mul_add(ry, (-lx).mul_add(rz, ly.mul_add(rw, lz * rx))),
        lw.mul_add(rz, lx.mul_add(ry, (-ly).mul_add(rx, lz * rw))),
    )
    .unwrap()
}

fn sample(
    start_time_s: f64,
    end_time_s: f64,
    branch: RenderContactBranch,
    disposition: RenderSampleDisposition,
) -> RenderTrajectorySampleInput {
    let positive_duration = end_time_s > start_time_s;
    let contact_active = positive_duration && branch == RenderContactBranch::Closed;
    let mut input = RenderTrajectorySampleInput {
        interval_start_time_s: start_time_s,
        time_s: end_time_s,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m: Vec3::new(10.0, -5.0, 3.25),
        orientation_body_to_world: UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.2)
            .unwrap()
            .components(),
        linear_momentum_world_kg_m_per_s: Vec3::new(4.0, -2.0, 1.0),
        angular_momentum_body_kg_m2_per_s: Vec3::new(2.0, 6.0, 12.0),
        symmetry_axis_world: Vec3::new(0.0, 0.0, 1.0),
        contact_branch: branch,
        contact_geometry: (branch == RenderContactBranch::Closed).then_some(
            RenderContactGeometry {
                point_world_m: Vec3::new(11.0, -5.0, 2.25),
                normal_world: Vec3::new(0.0, 0.0, 1.0),
                support_feature: RenderSupportFeature::ProfileFeature(4),
            },
        ),
        signed_gap_m: if branch == RenderContactBranch::Closed {
            0.0
        } else {
            1.0e-3
        },
        interval_contact_active: contact_active,
        interval_normal_force_n: if contact_active { 31.0 } else { 0.0 },
        contact_transitions: Vec::new(),
        base_mode: Some(RenderBaseModeState {
            displacement_m: 0.25,
            velocity_m_per_s: 0.2,
        }),
        channels: ChannelOwnership::default(),
        mechanical_energy_j: 8.0,
        energy_defect_j: 1.0e-10,
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
    base_frame: RenderBaseFrame,
    channel_availability: RenderChannelAvailability,
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
        base_frame,
        model_identity: identity("model"),
        channel_availability,
        configuration_identity: identity("configuration"),
        configuration_fingerprint: 0x434f_4e54_524f_4c01,
        timestep_s: (first.time_s - first.interval_start_time_s).max(f64::MIN_POSITIVE),
        producer_version: "control-stream-test-v1".into(),
        applicability: "reduced-model visualization and audio controls only".into(),
        no_claims: vec![
            "control frequency is not a physical acoustic frequency".into(),
            "event brackets do not contain an impulse magnitude".into(),
        ],
        authority: RenderTrajectoryAuthority::SimulationEvidence,
    }
}

fn default_base_frame() -> RenderBaseFrame {
    RenderBaseFrame {
        origin_world_m: Vec3::new(10.0, -5.0, 2.0),
        orientation_base_to_world: UnitQuaternion::from_axis_angle(
            Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_2,
        )
        .unwrap(),
    }
}

fn trajectory(
    inputs: Vec<RenderTrajectorySampleInput>,
    base_frame: RenderBaseFrame,
    availability: RenderChannelAvailability,
) -> RenderTrajectory {
    let mut trajectory_metadata = metadata(&inputs[0], base_frame, availability);
    trajectory_metadata.timestep_s = inputs
        .iter()
        .map(|input| input.time_s - input.interval_start_time_s)
        .fold(f64::MIN_POSITIVE, f64::max);
    RenderTrajectory::try_new(trajectory_metadata, inputs).unwrap()
}

fn wrench(force: Vec3, torque: Vec3, work_j: f64) -> ChannelWrench {
    ChannelWrench {
        force_world_n: force,
        torque_world_nm: torque,
        work_j,
    }
}

fn assert_close(left: f64, right: f64) {
    let scale = 1.0_f64.max(left.abs()).max(right.abs());
    assert!(
        (left - right).abs() <= 2.0e-12 * scale,
        "left={left:.17e}, right={right:.17e}, delta={:.17e}",
        left - right
    );
}

fn assert_vec_close(left: Vec3, right: Vec3) {
    assert_close(left.x, right.x);
    assert_close(left.y, right.y);
    assert_close(left.z, right.z);
}

#[test]
fn resumed_interval_derives_exact_frames_velocities_signed_power_and_source_binding() {
    let mut retained = sample(
        5.0,
        5.1,
        RenderContactBranch::Closed,
        RenderSampleDisposition::HorizonCensored,
    );
    retained.channels = ChannelOwnership {
        gravity: wrench(Vec3::new(0.0, 0.0, -19.6), Vec3::ZERO, 0.1),
        contact: wrench(Vec3::new(1.0, 2.0, 30.0), Vec3::new(3.0, 4.0, 5.0), -0.2),
        rolling: wrench(Vec3::ZERO, Vec3::new(0.0, 0.0, -0.1), -0.03),
        base: wrench(Vec3::ZERO, Vec3::ZERO, -0.04),
        gas: wrench(Vec3::new(-0.2, 0.1, 0.0), Vec3::new(0.0, 0.0, -0.01), 0.01),
    };
    let source = trajectory(
        vec![retained],
        default_base_frame(),
        RenderChannelAvailability::ALL_AVAILABLE,
    );
    let source_clone = source.clone();
    with_cx(false, |cx| {
        let controls = EulerControlStream::try_derive(&source, cx).unwrap();
        assert!(controls.is_bound_to(&source));
        assert!(!controls.is_bound_to(&source_clone));
        assert_eq!(
            controls.authority(),
            RenderTrajectoryAuthority::SimulationEvidence
        );
        assert_eq!(controls.visualization().len(), 1);
        assert_eq!(controls.audio().len(), 1);

        let visual = &controls.visualization()[0];
        assert_vec_close(
            visual.center_of_mass_velocity_world_m_per_s,
            Vec3::new(2.0, -1.0, 0.5),
        );
        assert_vec_close(
            visual.angular_velocity_body_rad_per_s,
            Vec3::new(1.0, 2.0, 3.0),
        );
        let orientation = visual.disc_pose.orientation();
        let expected_angular_world = orientation.rotate_body_to_world(Vec3::new(1.0, 2.0, 3.0));
        assert_vec_close(
            visual.angular_velocity_world_rad_per_s,
            expected_angular_world,
        );
        let contact = visual.contact.expect("closed endpoint contact");
        let contact_arm_world = Vec3::new(1.0, 0.0, -1.0);
        assert_vec_close(
            contact.point_body_m,
            orientation.rotate_world_to_body(contact_arm_world),
        );
        assert_vec_close(contact.point_base_m, Vec3::new(0.0, -1.0, 0.0));
        let expected_disc_point_velocity =
            Vec3::new(2.0, -1.0, 0.5).add(expected_angular_world.cross(contact_arm_world));
        assert_vec_close(
            contact.disc_point_velocity_world_m_per_s,
            expected_disc_point_velocity,
        );
        assert_vec_close(
            contact.base_point_velocity_world_m_per_s,
            Vec3::new(0.0, 0.0, 0.2),
        );
        assert_vec_close(
            contact.relative_point_velocity_world_m_per_s,
            expected_disc_point_velocity.sub(Vec3::new(0.0, 0.0, 0.2)),
        );

        let interval = &controls.audio()[0];
        assert!(!interval.visual_coverage.is_fully_bracketed());
        assert!(controls.fully_bracketed_audio().is_empty());
        assert!(controls.audio_visual_horizon().is_none());
        assert_eq!(interval.start_time_s.to_bits(), 5.0_f64.to_bits());
        assert_close(interval.duration_s, 0.1);
        assert_close(interval.mean_base_normal_contact_force_n.unwrap(), 30.0);
        let contact_channel = interval.channels.contact.available().unwrap();
        assert_close(contact_channel.signed_work_j, -0.2);
        assert_close(contact_channel.signed_mean_work_rate_w, -2.0);
        assert_vec_close(
            contact_channel.force_time_measure_world_n_s,
            Vec3::new(0.1, 0.2, 3.0),
        );
        let gas = interval.channels.gas.available().unwrap();
        assert_close(gas.signed_mean_work_rate_w, 0.1);
        assert!(
            controls
                .work_integral_checks()
                .within_tolerance(1.0e-15)
                .unwrap()
        );
        let replay = EulerControlStream::try_derive(&source, cx).unwrap();
        assert_eq!(replay.visualization(), controls.visualization());
        assert_eq!(replay.audio(), controls.audio());
        assert_eq!(
            replay.work_integral_checks(),
            controls.work_integral_checks()
        );
        let coarse = controls
            .boxcar_coarsen(NonZeroUsize::new(1).unwrap(), cx)
            .unwrap();
        assert!(coarse.is_bound_to(&source));
        assert!(!coarse.is_bound_to(&source_clone));
        assert!(coarse.fully_bracketed_bins().is_empty());
        assert!(coarse.audio_visual_horizon().is_none());
    });
}

#[test]
fn point_only_open_state_emits_no_false_interval_contact_or_geometry() {
    let point = sample(
        5.0,
        5.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    let source = trajectory(
        vec![point],
        default_base_frame(),
        RenderChannelAvailability::NONE_AVAILABLE,
    );
    with_cx(false, |cx| {
        let controls = EulerControlStream::try_derive(&source, cx).unwrap();
        assert!(controls.audio().is_empty());
        assert!(controls.visualization()[0].contact.is_none());
        assert_eq!(
            controls.visualization()[0].preceding_audio_interval_index,
            None
        );
        let checks = controls.work_integral_checks();
        assert!(checks.gravity.is_none());
        assert!(checks.contact.is_none());
        assert!(checks.rolling.is_none());
        assert!(checks.base.is_none());
        assert!(checks.gas.is_none());
        assert_eq!(
            controls
                .boxcar_coarsen(NonZeroUsize::new(2).unwrap(), cx)
                .unwrap_err(),
            ControlStreamError::NoPositiveDurationIntervals
        );
    });
}

#[test]
fn available_zero_and_unavailable_base_gas_remain_distinct() {
    let retained = sample(
        0.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    let availability = RenderChannelAvailability {
        base: false,
        gas: false,
        ..RenderChannelAvailability::ALL_AVAILABLE
    };
    let source = trajectory(vec![retained.clone()], default_base_frame(), availability);
    with_cx(false, |cx| {
        let controls = EulerControlStream::try_derive(&source, cx).unwrap();
        let channels = controls.audio()[0].channels;
        assert!(matches!(channels.base, ChannelControl::Unavailable));
        assert!(matches!(channels.gas, ChannelControl::Unavailable));
        assert!(matches!(channels.contact, ChannelControl::Available(_)));
        let checks = controls.work_integral_checks();
        assert!(checks.base.is_none());
        assert!(checks.gas.is_none());
        assert_eq!(
            checks.contact.unwrap().retained_work_j.to_bits(),
            0.0_f64.to_bits()
        );
    });

    let mut contradictory = retained;
    contradictory.channels.gas.work_j = 1.0;
    assert_eq!(
        RenderTrajectory::try_new(
            metadata(&contradictory, default_base_frame(), availability),
            vec![contradictory]
        )
        .unwrap_err(),
        RenderTrajectoryError::UnavailableChannelHasData {
            sample: 0,
            channel: "gas",
        }
    );
}

#[test]
fn opening_and_zero_force_reimpact_retain_timing_only_events_and_barriers() {
    let first = sample(
        0.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    let mut reimpact = sample(
        1.0,
        2.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::Continue,
    );
    reimpact.interval_contact_active = true;
    reimpact.interval_normal_force_n = 0.0;
    reimpact.channels.contact = ChannelWrench::default();
    reimpact.contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Reimpact,
        time_s: 1.5,
        bracket_start_s: 1.49,
        bracket_end_s: 1.51,
    }];
    let mut opening = sample(
        2.0,
        3.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    opening.interval_contact_active = true;
    opening.contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Opening,
        time_s: 2.25,
        bracket_start_s: 2.24,
        bracket_end_s: 2.26,
    }];
    let source = trajectory(
        vec![first, reimpact, opening],
        default_base_frame(),
        RenderChannelAvailability::ALL_AVAILABLE,
    );
    let mut eventful_preroll = sample(
        0.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    eventful_preroll.interval_contact_active = true;
    eventful_preroll.contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Opening,
        time_s: 0.5,
        bracket_start_s: 0.49,
        bracket_end_s: 0.51,
    }];
    let eventful_preroll_source = trajectory(
        vec![eventful_preroll],
        default_base_frame(),
        RenderChannelAvailability::ALL_AVAILABLE,
    );
    with_cx(false, |cx| {
        let controls = EulerControlStream::try_derive(&source, cx).unwrap();
        assert!(controls.audio()[1].interval_contact_active);
        assert_eq!(
            controls.audio()[1].mean_base_normal_contact_force_n,
            Some(0.0)
        );
        assert_eq!(
            controls.audio()[1].events[0].measure,
            ContactEventMeasure::TimingOnly
        );
        assert_close(controls.audio()[1].events[0].localization_width_s, 0.02);
        assert!(controls.visualization()[2].contact.is_none());

        let coarse = controls
            .boxcar_coarsen(NonZeroUsize::new(3).unwrap(), cx)
            .unwrap();
        assert_eq!(coarse.bins().len(), 3);
        assert!(!coarse.bins()[0].event_barrier);
        assert!(coarse.bins()[1].event_barrier);
        assert!(coarse.bins()[2].event_barrier);
        assert_eq!(
            coarse.bins()[1].events[0].kind,
            ContactTransitionKind::Reimpact
        );
        assert_eq!(
            coarse.bins()[2].events[0].kind,
            ContactTransitionKind::Opening
        );

        let preroll_controls =
            EulerControlStream::try_derive(&eventful_preroll_source, cx).unwrap();
        let preroll_coarse = preroll_controls
            .boxcar_coarsen(NonZeroUsize::new(8).unwrap(), cx)
            .unwrap();
        assert!(preroll_coarse.bins()[0].event_barrier);
        assert!(
            !preroll_coarse.bins()[0]
                .visual_coverage
                .is_fully_bracketed()
        );
    });
}

#[test]
fn whole_interval_boxcar_cancels_alternation_before_decimation_and_conserves_work() {
    let mut inputs = vec![sample(
        0.0,
        0.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    )];
    for index in 0..4 {
        let mut retained = sample(
            index as f64,
            (index + 1) as f64,
            RenderContactBranch::Open,
            if index == 3 {
                RenderSampleDisposition::HorizonCensored
            } else {
                RenderSampleDisposition::Continue
            },
        );
        let sign = if index % 2 == 0 { 1.0 } else { -1.0 };
        retained.channels.gas = wrench(Vec3::new(sign, 0.0, 0.0), Vec3::ZERO, sign);
        inputs.push(retained);
    }
    let source = trajectory(
        inputs,
        default_base_frame(),
        RenderChannelAvailability::ALL_AVAILABLE,
    );
    with_cx(false, |cx| {
        let controls = EulerControlStream::try_derive(&source, cx).unwrap();
        assert_eq!(
            controls
                .audio()
                .iter()
                .step_by(2)
                .map(|interval| {
                    interval
                        .channels
                        .gas
                        .available()
                        .unwrap()
                        .signed_mean_work_rate_w
                })
                .collect::<Vec<_>>(),
            vec![1.0, 1.0],
            "raw decimation aliases the alternating source"
        );
        let coarse = controls
            .boxcar_coarsen(NonZeroUsize::new(2).unwrap(), cx)
            .unwrap();
        assert_eq!(coarse.filter(), AudioControlFilter::WholeIntervalBoxcarV1);
        assert_eq!(coarse.bins().len(), 2);
        assert_eq!(coarse.fully_bracketed_bins(), coarse.bins());
        assert_eq!(
            coarse.audio_visual_horizon(),
            Some(fs_euler_disc_e2e::AudioVisualHorizon {
                start_time_s: 0.0,
                end_time_s: 4.0,
            })
        );
        for bin in coarse.bins() {
            let gas = bin.channels.gas.available().unwrap();
            assert_eq!(gas.signed_work_j.to_bits(), 0.0_f64.to_bits());
            assert_eq!(gas.signed_mean_work_rate_w.to_bits(), 0.0_f64.to_bits());
            assert_vec_close(gas.mean_force_world_n, Vec3::ZERO);
        }
        assert!(coarse.work_integral_checks().within_tolerance(0.0).unwrap());
    });
}

#[test]
fn preroll_isolated_and_unequal_duration_controls_are_weighted_conservatively() {
    let mut preroll = sample(
        0.0,
        2.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::Continue,
    );
    preroll.channels.contact = wrench(Vec3::new(0.0, 0.0, 1.0), Vec3::ZERO, -2.0);
    let mut short = sample(
        2.0,
        2.5,
        RenderContactBranch::Closed,
        RenderSampleDisposition::Continue,
    );
    short.channels.contact = wrench(Vec3::new(0.0, 0.0, 3.0), Vec3::ZERO, -1.0);
    let mut long = sample(
        2.5,
        4.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::HorizonCensored,
    );
    long.channels.contact = wrench(Vec3::new(0.0, 0.0, 5.0), Vec3::ZERO, -3.0);
    let source = trajectory(
        vec![preroll, short, long],
        default_base_frame(),
        RenderChannelAvailability::ALL_AVAILABLE,
    );

    with_cx(false, |cx| {
        let controls = EulerControlStream::try_derive(&source, cx).unwrap();
        assert!(!controls.audio()[0].visual_coverage.is_fully_bracketed());
        assert_eq!(controls.fully_bracketed_audio().len(), 2);
        assert_eq!(
            controls.audio_visual_horizon(),
            Some(fs_euler_disc_e2e::AudioVisualHorizon {
                start_time_s: 2.0,
                end_time_s: 4.0,
            })
        );

        let coarse = controls
            .boxcar_coarsen(NonZeroUsize::new(16).unwrap(), cx)
            .unwrap();
        assert_eq!(coarse.bins().len(), 2, "preroll must be its own bin");
        assert!(!coarse.bins()[0].visual_coverage.is_fully_bracketed());
        assert!(coarse.bins()[1].visual_coverage.is_fully_bracketed());
        assert_eq!(coarse.fully_bracketed_bins(), &coarse.bins()[1..]);
        assert_eq!(
            coarse.audio_visual_horizon(),
            Some(fs_euler_disc_e2e::AudioVisualHorizon {
                start_time_s: 2.0,
                end_time_s: 4.0,
            })
        );
        let contact = coarse.bins()[1].channels.contact.available().unwrap();
        assert_close(contact.mean_force_world_n.z, 4.5);
        assert_close(
            coarse.bins()[1].mean_base_normal_contact_force_n.unwrap(),
            4.5,
        );
        assert_close(contact.force_time_measure_world_n_s.z, 9.0);
        assert_close(contact.signed_work_j, -4.0);
        assert_close(contact.signed_mean_work_rate_w, -2.0);
        assert!(
            coarse
                .work_integral_checks()
                .within_tolerance(1.0e-15)
                .unwrap()
        );
    });
}

#[test]
fn z_up_rigid_transform_rotates_world_controls_and_preserves_local_controls() {
    let mut original = sample(
        0.0,
        1.0,
        RenderContactBranch::Closed,
        RenderSampleDisposition::HorizonCensored,
    );
    original.channels.contact = wrench(Vec3::new(1.0, 2.0, 3.0), Vec3::new(4.0, 5.0, 6.0), -0.75);
    let original_base = RenderBaseFrame {
        origin_world_m: Vec3::new(10.0, -5.0, 2.0),
        orientation_base_to_world: UnitQuaternion::IDENTITY,
    };
    let original_source = trajectory(
        vec![original.clone()],
        original_base,
        RenderChannelAvailability::ALL_AVAILABLE,
    );

    let yaw = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.37).unwrap();
    let translation = Vec3::new(-7.0, 4.0, 2.0);
    let transform_point = |point: Vec3| yaw.rotate_body_to_world(point).add(translation);
    let mut transformed = original;
    transformed.center_of_mass_world_m = transform_point(transformed.center_of_mass_world_m);
    let original_orientation = UnitQuaternion::new(
        transformed.orientation_body_to_world[0],
        transformed.orientation_body_to_world[1],
        transformed.orientation_body_to_world[2],
        transformed.orientation_body_to_world[3],
    )
    .unwrap();
    transformed.orientation_body_to_world = compose(yaw, original_orientation).components();
    transformed.linear_momentum_world_kg_m_per_s =
        yaw.rotate_body_to_world(transformed.linear_momentum_world_kg_m_per_s);
    let mut geometry = transformed.contact_geometry.unwrap();
    geometry.point_world_m = transform_point(geometry.point_world_m);
    geometry.normal_world = yaw.rotate_body_to_world(geometry.normal_world);
    transformed.contact_geometry = Some(geometry);
    for channel in [
        &mut transformed.channels.gravity,
        &mut transformed.channels.contact,
        &mut transformed.channels.rolling,
        &mut transformed.channels.base,
        &mut transformed.channels.gas,
    ] {
        channel.force_world_n = yaw.rotate_body_to_world(channel.force_world_n);
        channel.torque_world_nm = yaw.rotate_body_to_world(channel.torque_world_nm);
    }
    refresh_pose_diagnostics(&mut transformed);
    let transformed_base = RenderBaseFrame {
        origin_world_m: transform_point(original_base.origin_world_m),
        orientation_base_to_world: yaw,
    };
    let transformed_source = trajectory(
        vec![transformed],
        transformed_base,
        RenderChannelAvailability::ALL_AVAILABLE,
    );

    with_cx(false, |cx| {
        let original_controls = EulerControlStream::try_derive(&original_source, cx).unwrap();
        let transformed_controls = EulerControlStream::try_derive(&transformed_source, cx).unwrap();
        let original_visual = &original_controls.visualization()[0];
        let transformed_visual = &transformed_controls.visualization()[0];
        let original_contact = original_visual.contact.unwrap();
        let transformed_contact = transformed_visual.contact.unwrap();
        assert_vec_close(
            transformed_contact.point_body_m,
            original_contact.point_body_m,
        );
        assert_vec_close(
            transformed_contact.point_base_m,
            original_contact.point_base_m,
        );
        assert_vec_close(
            transformed_contact.normal_body,
            original_contact.normal_body,
        );
        assert_vec_close(
            transformed_contact.normal_base,
            original_contact.normal_base,
        );
        assert_vec_close(
            transformed_contact.point_world_m,
            transform_point(original_contact.point_world_m),
        );
        assert_vec_close(
            transformed_visual.center_of_mass_velocity_world_m_per_s,
            yaw.rotate_body_to_world(original_visual.center_of_mass_velocity_world_m_per_s),
        );
        assert_vec_close(
            transformed_visual.angular_velocity_body_rad_per_s,
            original_visual.angular_velocity_body_rad_per_s,
        );
        assert_vec_close(
            transformed_visual.angular_velocity_world_rad_per_s,
            yaw.rotate_body_to_world(original_visual.angular_velocity_world_rad_per_s),
        );
        let original_channel = original_controls.audio()[0]
            .channels
            .contact
            .available()
            .unwrap();
        let transformed_channel = transformed_controls.audio()[0]
            .channels
            .contact
            .available()
            .unwrap();
        assert_vec_close(
            transformed_channel.mean_force_world_n,
            yaw.rotate_body_to_world(original_channel.mean_force_world_n),
        );
        assert_eq!(
            transformed_channel.signed_work_j.to_bits(),
            original_channel.signed_work_j.to_bits()
        );
        assert_eq!(
            transformed_channel.signed_mean_work_rate_w.to_bits(),
            original_channel.signed_mean_work_rate_w.to_bits()
        );
    });
}

#[test]
fn extreme_finite_controls_are_not_saturated_and_overflow_refuses_atomically() {
    let mut large = sample(
        0.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    large.channels.gas.work_j = f64::MAX / 4.0;
    let large_source = trajectory(
        vec![large],
        default_base_frame(),
        RenderChannelAvailability::ALL_AVAILABLE,
    );
    with_cx(false, |cx| {
        let controls = EulerControlStream::try_derive(&large_source, cx).unwrap();
        let raw = controls.audio()[0].channels.gas.available().unwrap();
        assert_eq!(raw.signed_work_j.to_bits(), (f64::MAX / 4.0).to_bits());
        let coarse = controls
            .boxcar_coarsen(NonZeroUsize::new(1).unwrap(), cx)
            .unwrap();
        assert_eq!(
            coarse.bins()[0]
                .channels
                .gas
                .available()
                .unwrap()
                .signed_work_j
                .to_bits(),
            raw.signed_work_j.to_bits()
        );
    });

    let mut overflow = sample(
        0.0,
        f64::MIN_POSITIVE,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    overflow.channels.gas.work_j = f64::MAX;
    let overflow_source = trajectory(
        vec![overflow],
        default_base_frame(),
        RenderChannelAvailability::ALL_AVAILABLE,
    );
    with_cx(false, |cx| {
        assert_eq!(
            EulerControlStream::try_derive(&overflow_source, cx).unwrap_err(),
            ControlStreamError::NonFiniteDerived {
                sample: 0,
                field: "signed_mean_work_rate_w",
            }
        );
    });

    let mut first = sample(
        0.0,
        1.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::Continue,
    );
    first.channels.gas.work_j = f64::MAX * 0.75;
    let mut second = sample(
        1.0,
        2.0,
        RenderContactBranch::Open,
        RenderSampleDisposition::HorizonCensored,
    );
    second.channels.gas.work_j = f64::MAX * 0.75;
    let aggregate_overflow_source = trajectory(
        vec![first, second],
        default_base_frame(),
        RenderChannelAvailability::ALL_AVAILABLE,
    );
    with_cx(false, |cx| {
        assert_eq!(
            EulerControlStream::try_derive(&aggregate_overflow_source, cx).unwrap_err(),
            ControlStreamError::NonFiniteDerived {
                sample: 1,
                field: "work_integral_accumulator",
            }
        );
    });
}

#[test]
fn derivation_and_coarsening_observe_precancel_without_partial_publication() {
    let source = trajectory(
        vec![sample(
            0.0,
            1.0,
            RenderContactBranch::Open,
            RenderSampleDisposition::HorizonCensored,
        )],
        default_base_frame(),
        RenderChannelAvailability::ALL_AVAILABLE,
    );
    with_cx(true, |cx| {
        assert_eq!(
            EulerControlStream::try_derive(&source, cx).unwrap_err(),
            ControlStreamError::Cancelled
        );
    });
    with_cx(false, |cx| {
        let controls = EulerControlStream::try_derive(&source, cx).unwrap();
        with_cx(true, |cancelled_cx| {
            assert_eq!(
                controls
                    .boxcar_coarsen(NonZeroUsize::new(2).unwrap(), cancelled_cx)
                    .unwrap_err(),
                ControlStreamError::Cancelled
            );
        });
    });
}
