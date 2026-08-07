//! G0/G3 checks for animation-grade Euler trajectory semantics.

use core::f64::consts::{PI, TAU};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::coupled_runner::{
    ChannelOwnership, ContactTransitionKind, CoupledControls, CoupledFactors, CoupledInitialState,
    CoupledRun, run_closed_reduced,
};
use fs_euler_disc_e2e::reduced_decay::{
    ReducedDecayRun, THORNE_2026_RENDER_VALIDITY_CUTOFF_THETA_RAD, Thorne2026SteelGlassBenchmark,
    run_thorne_2026_steel_glass_benchmark,
};
use fs_euler_disc_e2e::{
    DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, REDUCED_DECAY_RENDER_TAIL_HORIZON_S,
    RenderBaseFrame, RenderBaseModeState, RenderChannelAvailability, RenderContactBranch,
    RenderContactGeometry, RenderContactTransition, RenderMassProperties,
    RenderNumericalRefusalReason, RenderSampleDisposition, RenderSupportFeature,
    RenderTerminalEvent, RenderTrajectory, RenderTrajectoryAuthority, RenderTrajectoryError,
    RenderTrajectoryMetadata, RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const THETA: f64 = 0.35;

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5448_4f52_4e45_3236,
                kernel_id: 0x5245_4e44,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn thorne_trajectory(cx: &Cx<'_>) -> (ReducedDecayRun, RenderTrajectory) {
    let benchmark = Thorne2026SteelGlassBenchmark::ambient().expect("ambient benchmark");
    let profile = benchmark
        .resolve_specimen(cx)
        .expect("source-bound specimen");
    let run = run_thorne_2026_steel_glass_benchmark(&benchmark).expect("source-bound decay");
    let trajectory = RenderTrajectory::from_reduced_decay_run(&run, &profile, cx)
        .expect("grounded cinematic trajectory");
    (run, trajectory)
}

fn identity(label: &str) -> ContentHash {
    hash_domain("org.frankensim.test.render-trajectory.v1", label.as_bytes())
}

fn assert_close(left: f64, right: f64) {
    let scale = 1.0_f64.max(left.abs()).max(right.abs());
    assert!((left - right).abs() <= 1.0e-12 * scale, "{left} != {right}");
}

fn mass() -> MassProperties {
    MassProperties::new(0.4, Vec3::ZERO, Vec3::new(0.001, 0.001, 0.002)).unwrap()
}

fn quaternion() -> [f64; 4] {
    [(THETA * 0.5).cos(), 0.0, (THETA * 0.5).sin(), 0.0]
}

fn state_from(sample: &RenderTrajectorySampleInput) -> RigidBodyState {
    let q = sample.orientation_body_to_world;
    let orientation = UnitQuaternion::new(q[0], q[1], q[2], q[3]).unwrap();
    RigidBodyState::new(
        Pose::new(sample.center_of_mass_world_m, orientation).unwrap(),
        sample.linear_momentum_world_kg_m_per_s,
        sample.angular_momentum_body_kg_m2_per_s,
    )
    .unwrap()
}

fn sample(time_s: f64, disposition: RenderSampleDisposition) -> RenderTrajectorySampleInput {
    let q = quaternion();
    let orientation = UnitQuaternion::new(q[0], q[1], q[2], q[3]).unwrap();
    let mut sample = RenderTrajectorySampleInput {
        interval_start_time_s: 0.0,
        time_s,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m: Vec3::new(0.1, -0.2, 0.01),
        orientation_body_to_world: q,
        linear_momentum_world_kg_m_per_s: Vec3::new(0.02, -0.01, 0.0),
        angular_momentum_body_kg_m2_per_s: Vec3::new(0.0, 0.0, 0.02),
        symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        contact_branch: RenderContactBranch::Open,
        contact_geometry: None,
        signed_gap_m: 0.001,
        interval_contact_active: false,
        interval_normal_force_n: 0.0,
        contact_transitions: Vec::new(),
        base_mode: Some(RenderBaseModeState {
            displacement_m: 1.0e-6,
            velocity_m_per_s: -2.0e-5,
        }),
        channels: ChannelOwnership::default(),
        mechanical_energy_j: 0.2,
        energy_defect_j: 1.0e-8,
        qois: DerivedEulerQois {
            inclination_rad: 0.0,
            precession_rad_per_s: 0.0,
            spin_rad_per_s: 0.0,
            precession_acceleration_rad_per_s2: 0.0,
        },
        disposition,
        terminal_event: None,
    };
    sample.qois = DerivedEulerQois::from_state(state_from(&sample), mass(), 0.0).unwrap();
    sample
}

fn metadata() -> RenderTrajectoryMetadata {
    let first = sample(0.01, RenderSampleDisposition::HorizonCensored);
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
        initial_state: state_from(&first),
        initial_base_mode: RenderBaseModeState {
            displacement_m: 0.0,
            velocity_m_per_s: 0.0,
        },
        base_model_identity: identity("base"),
        base_frame: RenderBaseFrame {
            origin_world_m: Vec3::ZERO,
            orientation_base_to_world: UnitQuaternion::IDENTITY,
        },
        model_identity: identity("model"),
        channel_availability: RenderChannelAvailability::ALL_AVAILABLE,
        configuration_identity: identity("configuration"),
        configuration_fingerprint: 0x0123_4567_89ab_cdef,
        timestep_s: 0.01,
        producer_version: "fs-euler-disc-e2e-test-v1".into(),
        applicability: "reduced Euler-disc simulation inside its declared validity domain".into(),
        no_claims: vec![
            "not calibrated physical truth".into(),
            "no continuum finite-contact-patch claim".into(),
        ],
        authority: RenderTrajectoryAuthority::SimulationEvidence,
    }
}

fn runner_factors() -> CoupledFactors {
    let radius_m = 0.038;
    let thickness_m = 0.006;
    let mass_kg = std::f64::consts::PI * radius_m * radius_m * thickness_m * 2_680.0;
    CoupledFactors {
        mass_kg,
        radius_m,
        thickness_m,
        transverse_inertia_kg_m2: mass_kg * (3.0 * radius_m.powi(2) + thickness_m.powi(2)) / 12.0,
        axial_inertia_kg_m2: 0.5 * mass_kg * radius_m.powi(2),
        gravity_m_per_s2: 9.806_65,
        sliding_friction_coefficient: 0.0,
        rolling_resistance_m: 0.0,
        contact_stiffness_n_per_m: 8.0e4,
        contact_damping_n_s_per_m: 3.0,
        base_effective_mass_kg: 0.25,
        base_stiffness_n_per_m: 4.0e4,
        base_damping_n_s_per_m: 4.0,
        gas_rotational_damping_n_m_s: 0.0,
        gas_translation_damping_n_s_per_m: 0.0,
    }
}

fn runner_metadata(run: &CoupledRun) -> RenderTrajectoryMetadata {
    RenderTrajectoryMetadata {
        schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        specimen_profile_identity: identity("runner-profile"),
        specimen_chart_identity: identity("runner-chart"),
        mass_properties: RenderMassProperties {
            identity: identity("runner-mass"),
            properties: run.mass_properties,
        },
        initial_state: run.configuration_initial_state,
        initial_base_mode: RenderBaseModeState {
            displacement_m: run.configuration_initial_base_deflection_m,
            velocity_m_per_s: run.configuration_initial_base_velocity_m_per_s,
        },
        base_model_identity: identity("runner-base"),
        base_frame: RenderBaseFrame {
            origin_world_m: Vec3::ZERO,
            orientation_base_to_world: UnitQuaternion::IDENTITY,
        },
        model_identity: identity("runner-model"),
        channel_availability: RenderChannelAvailability::ALL_AVAILABLE,
        configuration_identity: identity("runner-configuration"),
        configuration_fingerprint: run.checkpoint.configuration_fingerprint,
        timestep_s: run.macro_timestep_s,
        producer_version: "fs-euler-disc-e2e-test-v1".into(),
        applicability: run.applicability.into(),
        no_claims: vec![
            "reduced-model simulation evidence, not calibrated physical truth".into(),
            run.model_disagreement.into(),
        ],
        authority: RenderTrajectoryAuthority::SimulationEvidence,
    }
}

fn two_samples() -> Vec<RenderTrajectorySampleInput> {
    let first = sample(0.01, RenderSampleDisposition::Continue);
    let mut second = sample(0.02, RenderSampleDisposition::HorizonCensored);
    second.interval_start_time_s = first.time_s;
    vec![first, second]
}

#[test]
fn valid_trajectory_exposes_complete_canonical_state() {
    let trajectory = RenderTrajectory::try_new(metadata(), two_samples()).unwrap();
    assert_eq!(trajectory.samples().len(), 2);
    let retained = &trajectory.samples()[1];
    assert_eq!(
        retained.state().pose().position_world(),
        Vec3::new(0.1, -0.2, 0.01)
    );
    let expected = quaternion();
    for (left, right) in retained
        .state()
        .pose()
        .orientation()
        .components()
        .into_iter()
        .zip(expected)
    {
        assert_close(left, right);
    }
    assert!(retained.input().base_mode.is_some());
    assert_eq!(
        retained.input().disposition,
        RenderSampleDisposition::HorizonCensored
    );
    assert_eq!(
        trajectory.metadata().authority,
        RenderTrajectoryAuthority::SimulationEvidence
    );
}

#[test]
fn quaternion_double_cover_canonicalizes_to_equal_samples() {
    let positive = sample(0.01, RenderSampleDisposition::HorizonCensored);
    let mut negative = positive.clone();
    for component in &mut negative.orientation_body_to_world {
        *component = -*component;
    }
    let positive = RenderTrajectory::try_new(metadata(), vec![positive]).unwrap();
    let negative = RenderTrajectory::try_new(metadata(), vec![negative]).unwrap();
    assert_eq!(positive.samples(), negative.samples());
}

#[test]
fn non_unit_and_non_finite_quaternions_refuse() {
    let mut non_unit = sample(0.01, RenderSampleDisposition::HorizonCensored);
    non_unit.orientation_body_to_world = [2.0, 0.0, 0.0, 0.0];
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![non_unit]).unwrap_err(),
        RenderTrajectoryError::QuaternionNotUnit(0)
    );

    let mut non_finite = sample(0.01, RenderSampleDisposition::HorizonCensored);
    non_finite.orientation_body_to_world[2] = f64::NAN;
    assert!(matches!(
        RenderTrajectory::try_new(metadata(), vec![non_finite]),
        Err(RenderTrajectoryError::NonFinite {
            sample: Some(0),
            field: "orientation_body_to_world"
        })
    ));
}

#[test]
fn duplicate_and_nonmonotone_times_refuse() {
    let mut duplicate = two_samples();
    duplicate[1].time_s = duplicate[0].time_s;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), duplicate).unwrap_err(),
        RenderTrajectoryError::NonMonotoneTime(1)
    );

    let mut decreasing = two_samples();
    decreasing[1].time_s = 0.005;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), decreasing).unwrap_err(),
        RenderTrajectoryError::NonMonotoneTime(1)
    );
}

#[test]
fn interval_clock_base_frame_and_contact_activity_contracts_refuse_contradictions() {
    let mut noncontiguous = two_samples();
    noncontiguous[1].interval_start_time_s = noncontiguous[0].time_s + 1.0e-6;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), noncontiguous).unwrap_err(),
        RenderTrajectoryError::InvalidIntervalStart(1)
    );

    let overlong = sample(0.02, RenderSampleDisposition::HorizonCensored);
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![overlong]).unwrap_err(),
        RenderTrajectoryError::IntervalExceedsDeclaredTimestep(0)
    );

    let large_clock_start_s = 100.0;
    let timestep_s = 2.0e-5;
    let mut runner_shaped = sample(
        large_clock_start_s + timestep_s,
        RenderSampleDisposition::HorizonCensored,
    );
    runner_shaped.interval_start_time_s = large_clock_start_s;
    let mut runner_shaped_metadata = metadata();
    runner_shaped_metadata.timestep_s = timestep_s;
    RenderTrajectory::try_new(runner_shaped_metadata, vec![runner_shaped])
        .expect("producer-shaped full step remains admitted at a large absolute clock");

    // The cinematic bridge rebases each endpoint independently.  This pair
    // comes from the exact source-clock arithmetic used by its Thorne tail;
    // cancellation leaves a few fs above `start + timestep` at a small
    // output time.  It is a binary64 clock residue, not an overlong physical
    // interval, and must remain admissible.
    let source_run = run_thorne_2026_steel_glass_benchmark(
        &Thorne2026SteelGlassBenchmark::ambient().expect("ambient benchmark"),
    )
    .expect("source-bound reduced decay");
    let source_final_time_s = source_run.final_sample().expect("source cutoff").time_s;
    let source_origin_s = source_final_time_s - REDUCED_DECAY_RENDER_TAIL_HORIZON_S;
    let first_retained = source_run
        .samples
        .partition_point(|source| source.time_s < source_origin_s);
    let source_interval_start_s = source_run.samples[first_retained + 1].time_s - source_origin_s;
    let source_interval_end_s = source_run.samples[first_retained + 2].time_s - source_origin_s;
    assert!(
        source_interval_end_s > source_interval_start_s + source_run.parameters.timestep_s,
        "regression requires the independent-rebase cancellation residue"
    );
    let mut rebased = sample(
        source_interval_end_s,
        RenderSampleDisposition::HorizonCensored,
    );
    rebased.interval_start_time_s = source_interval_start_s;
    let mut rebased_metadata = metadata();
    rebased_metadata.timestep_s = source_run.parameters.timestep_s;
    RenderTrajectory::try_new(rebased_metadata, vec![rebased])
        .expect("bounded sub-second rebase residue remains admissible");

    let infinity_bits = f64::INFINITY.to_bits();
    let mut near_maximum = sample(f64::MAX, RenderSampleDisposition::HorizonCensored);
    near_maximum.interval_start_time_s = 0.0;
    let mut near_maximum_metadata = metadata();
    near_maximum_metadata.timestep_s = f64::from_bits(infinity_bits - 16);
    RenderTrajectory::try_new(near_maximum_metadata, vec![near_maximum.clone()])
        .expect("ULP tolerance clamps at positive infinity without producing NaN");
    let mut overlong_near_maximum_metadata = metadata();
    overlong_near_maximum_metadata.timestep_s = f64::from_bits(infinity_bits - 64);
    assert_eq!(
        RenderTrajectory::try_new(overlong_near_maximum_metadata, vec![near_maximum]).unwrap_err(),
        RenderTrajectoryError::IntervalExceedsDeclaredTimestep(0)
    );

    let mut point_with_interval_data = sample(0.01, RenderSampleDisposition::HorizonCensored);
    point_with_interval_data.interval_start_time_s = point_with_interval_data.time_s;
    point_with_interval_data.channels.gravity.work_j = 1.0;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![point_with_interval_data]).unwrap_err(),
        RenderTrajectoryError::ZeroDurationIntervalData(0)
    );

    let mut inactive_contact = sample(0.01, RenderSampleDisposition::HorizonCensored);
    inactive_contact.channels.contact.work_j = -1.0e-3;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![inactive_contact]).unwrap_err(),
        RenderTrajectoryError::InactiveContactHasIntervalData(0)
    );

    let mut negative_mean_normal = sample(0.01, RenderSampleDisposition::HorizonCensored);
    negative_mean_normal.interval_contact_active = true;
    negative_mean_normal.channels.contact.force_world_n = Vec3::new(0.0, 0.0, -1.0);
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![negative_mean_normal]).unwrap_err(),
        RenderTrajectoryError::NegativeNormalForce(0)
    );

    let mut inactive_closed_without_root = sample(0.01, RenderSampleDisposition::HorizonCensored);
    inactive_closed_without_root.contact_branch = RenderContactBranch::Closed;
    inactive_closed_without_root.signed_gap_m = 0.0;
    inactive_closed_without_root.contact_geometry = Some(RenderContactGeometry {
        point_world_m: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        support_feature: RenderSupportFeature::CylinderRim,
    });
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![inactive_closed_without_root]).unwrap_err(),
        RenderTrajectoryError::InactiveContactHasIntervalData(0)
    );

    let mut inactive_interior_reimpact = sample(0.01, RenderSampleDisposition::HorizonCensored);
    inactive_interior_reimpact.contact_branch = RenderContactBranch::Closed;
    inactive_interior_reimpact.signed_gap_m = 0.0;
    inactive_interior_reimpact.contact_geometry = Some(RenderContactGeometry {
        point_world_m: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        support_feature: RenderSupportFeature::CylinderRim,
    });
    inactive_interior_reimpact.contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Reimpact,
        time_s: 0.005,
        bracket_start_s: 0.004,
        bracket_end_s: 0.006,
    }];
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![inactive_interior_reimpact]).unwrap_err(),
        RenderTrajectoryError::InactiveContactHasIntervalData(0)
    );

    let mut active_open_without_opening = sample(0.01, RenderSampleDisposition::HorizonCensored);
    active_open_without_opening.interval_contact_active = true;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![active_open_without_opening]).unwrap_err(),
        RenderTrajectoryError::ContactTransitionBranchMismatch(0)
    );

    let mut open_to_closed_without_reimpact = two_samples();
    open_to_closed_without_reimpact[1].contact_branch = RenderContactBranch::Closed;
    open_to_closed_without_reimpact[1].contact_geometry = Some(RenderContactGeometry {
        point_world_m: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        support_feature: RenderSupportFeature::CylinderRim,
    });
    open_to_closed_without_reimpact[1].signed_gap_m = 0.0;
    open_to_closed_without_reimpact[1].interval_contact_active = true;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), open_to_closed_without_reimpact).unwrap_err(),
        RenderTrajectoryError::ContactTransitionBranchMismatch(1)
    );

    let mut closed_to_open_without_opening = two_samples();
    closed_to_open_without_opening[0].contact_branch = RenderContactBranch::Closed;
    closed_to_open_without_opening[0].contact_geometry = Some(RenderContactGeometry {
        point_world_m: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        support_feature: RenderSupportFeature::CylinderRim,
    });
    closed_to_open_without_opening[0].signed_gap_m = 0.0;
    closed_to_open_without_opening[0].interval_contact_active = true;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), closed_to_open_without_opening).unwrap_err(),
        RenderTrajectoryError::ContactTransitionBranchMismatch(1)
    );

    let half_angle = 0.05_f64;
    let mut tilted_base = metadata();
    tilted_base.base_frame.orientation_base_to_world =
        UnitQuaternion::new(half_angle.cos(), half_angle.sin(), 0.0, 0.0).unwrap();
    assert_eq!(
        RenderTrajectory::try_new(
            tilted_base,
            vec![sample(0.01, RenderSampleDisposition::HorizonCensored)]
        )
        .unwrap_err(),
        RenderTrajectoryError::UnsupportedBaseFrame
    );
}

#[test]
fn frame_and_unit_cross_wiring_refuses() {
    let mut wrong_frame = sample(0.01, RenderSampleDisposition::HorizonCensored);
    wrong_frame.world_frame = RenderWorldFrame::RightHandedYUp;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![wrong_frame]).unwrap_err(),
        RenderTrajectoryError::FrameMismatch(0)
    );

    let mut wrong_units = sample(0.01, RenderSampleDisposition::HorizonCensored);
    wrong_units.units = RenderUnitSystem::SiDegrees;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![wrong_units]).unwrap_err(),
        RenderTrajectoryError::UnitMismatch(0)
    );
}

#[test]
fn contact_geometry_and_base_state_are_branch_consistent_and_complete() {
    let mut closed = sample(0.01, RenderSampleDisposition::HorizonCensored);
    closed.contact_branch = RenderContactBranch::Closed;
    closed.signed_gap_m = -1.0e-8;
    closed.interval_contact_active = true;
    closed.interval_normal_force_n = 3.0;
    closed.contact_geometry = Some(RenderContactGeometry {
        point_world_m: Vec3::new(0.0, 0.0, 0.0),
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        support_feature: RenderSupportFeature::ProfileFeature(7),
    });
    closed.contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Reimpact,
        time_s: 0.005,
        bracket_start_s: 0.004,
        bracket_end_s: 0.006,
    }];
    let closed = RenderTrajectory::try_new(metadata(), vec![closed]).unwrap();
    assert_eq!(
        closed.samples()[0]
            .input()
            .contact_geometry
            .unwrap()
            .support_feature,
        RenderSupportFeature::ProfileFeature(7)
    );

    let mut open_with_geometry = sample(0.01, RenderSampleDisposition::HorizonCensored);
    open_with_geometry.contact_geometry = Some(RenderContactGeometry {
        point_world_m: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        support_feature: RenderSupportFeature::CylinderRim,
    });
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![open_with_geometry]).unwrap_err(),
        RenderTrajectoryError::ContactGeometryMismatch(0)
    );

    let mut closed_without_geometry = sample(0.01, RenderSampleDisposition::HorizonCensored);
    closed_without_geometry.contact_branch = RenderContactBranch::Closed;
    closed_without_geometry.signed_gap_m = -1.0e-8;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![closed_without_geometry]).unwrap_err(),
        RenderTrajectoryError::ContactGeometryMismatch(0)
    );

    let mut missing_base = sample(0.01, RenderSampleDisposition::HorizonCensored);
    missing_base.base_mode = None;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![missing_base]).unwrap_err(),
        RenderTrajectoryError::MissingBaseState(0)
    );

    let mut wrong_transition_branch = sample(0.01, RenderSampleDisposition::HorizonCensored);
    wrong_transition_branch.contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Reimpact,
        time_s: 0.005,
        bracket_start_s: 0.004,
        bracket_end_s: 0.006,
    }];
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![wrong_transition_branch]).unwrap_err(),
        RenderTrajectoryError::ContactTransitionBranchMismatch(0)
    );
}

#[test]
fn localized_transition_brackets_and_terminal_placement_are_checked() {
    let mut inputs = two_samples();
    inputs[1].interval_contact_active = true;
    inputs[1].contact_transitions = vec![
        RenderContactTransition {
            kind: ContactTransitionKind::Reimpact,
            time_s: 0.012,
            bracket_start_s: 0.011,
            bracket_end_s: 0.013,
        },
        RenderContactTransition {
            kind: ContactTransitionKind::Opening,
            time_s: 0.015,
            bracket_start_s: 0.014,
            bracket_end_s: 0.016,
        },
    ];
    RenderTrajectory::try_new(metadata(), inputs).unwrap();

    let mut invalid = two_samples();
    invalid[1].contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Reimpact,
        time_s: 0.03,
        bracket_start_s: 0.014,
        bracket_end_s: 0.031,
    }];
    assert_eq!(
        RenderTrajectory::try_new(metadata(), invalid).unwrap_err(),
        RenderTrajectoryError::InvalidTransition {
            sample: 1,
            transition: 0
        }
    );

    let early_terminal = vec![
        sample(0.01, RenderSampleDisposition::HorizonCensored),
        sample(0.02, RenderSampleDisposition::HorizonCensored),
    ];
    assert_eq!(
        RenderTrajectory::try_new(metadata(), early_terminal).unwrap_err(),
        RenderTrajectoryError::TerminalBeforeFinalSample(0)
    );
    assert_eq!(
        RenderTrajectory::try_new(
            metadata(),
            vec![sample(0.01, RenderSampleDisposition::Continue)]
        )
        .unwrap_err(),
        RenderTrajectoryError::MissingFinalDisposition
    );

    let terminal_without_event = sample(0.01, RenderSampleDisposition::TerminalInclination);
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![terminal_without_event]).unwrap_err(),
        RenderTrajectoryError::TerminalEventMismatch(0)
    );
    let mut terminal = sample(0.01, RenderSampleDisposition::TerminalInclination);
    terminal.terminal_event = Some(RenderTerminalEvent {
        time_s: 0.0095,
        bracket_start_s: 0.009,
        bracket_end_s: 0.01,
    });
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![terminal]).unwrap_err(),
        RenderTrajectoryError::TerminalEventMismatch(0)
    );
    let mut terminal_at_retained_root = sample(0.01, RenderSampleDisposition::TerminalInclination);
    terminal_at_retained_root.terminal_event = Some(RenderTerminalEvent {
        time_s: 0.01,
        bracket_start_s: 0.009,
        bracket_end_s: 0.011,
    });
    RenderTrajectory::try_new(metadata(), vec![terminal_at_retained_root]).unwrap();

    let mut excessive_terminal_overhang =
        sample(0.01, RenderSampleDisposition::TerminalInclination);
    excessive_terminal_overhang.terminal_event = Some(RenderTerminalEvent {
        time_s: 0.01,
        bracket_start_s: 0.009,
        bracket_end_s: 0.020_001,
    });
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![excessive_terminal_overhang]).unwrap_err(),
        RenderTrajectoryError::TerminalEventMismatch(0)
    );

    let invalid_refusal = sample(
        0.01,
        RenderSampleDisposition::NumericalRefusal(RenderNumericalRefusalReason::BackendSpecific(0)),
    );
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![invalid_refusal]).unwrap_err(),
        RenderTrajectoryError::InvalidNumericalRefusalCode(0)
    );

    let mut excessive_reimpact_overhang = sample(
        0.01,
        RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::ReimpactLimitExceeded,
        ),
    );
    excessive_reimpact_overhang.contact_branch = RenderContactBranch::Closed;
    excessive_reimpact_overhang.signed_gap_m = -1.0e-8;
    excessive_reimpact_overhang.contact_geometry = Some(RenderContactGeometry {
        point_world_m: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        support_feature: RenderSupportFeature::CylinderRim,
    });
    excessive_reimpact_overhang.interval_contact_active = true;
    excessive_reimpact_overhang.contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Reimpact,
        time_s: 0.01,
        bracket_start_s: 0.009,
        bracket_end_s: 0.020_001,
    }];
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![excessive_reimpact_overhang]).unwrap_err(),
        RenderTrajectoryError::InvalidTransition {
            sample: 0,
            transition: 0,
        }
    );

    let mut interior_reimpact_refusal = sample(
        0.01,
        RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::ReimpactLimitExceeded,
        ),
    );
    interior_reimpact_refusal.contact_branch = RenderContactBranch::Closed;
    interior_reimpact_refusal.signed_gap_m = -1.0e-8;
    interior_reimpact_refusal.contact_geometry = Some(RenderContactGeometry {
        point_world_m: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        support_feature: RenderSupportFeature::CylinderRim,
    });
    interior_reimpact_refusal.interval_contact_active = true;
    interior_reimpact_refusal.contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Reimpact,
        time_s: 0.0095,
        bracket_start_s: 0.009,
        bracket_end_s: 0.01,
    }];
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![interior_reimpact_refusal]).unwrap_err(),
        RenderTrajectoryError::InvalidTransition {
            sample: 0,
            transition: 0,
        }
    );
}

#[test]
fn redundant_qoi_disagreement_refuses() {
    let mut wrong = sample(0.01, RenderSampleDisposition::HorizonCensored);
    wrong.qois.spin_rad_per_s += 0.1;
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![wrong]).unwrap_err(),
        RenderTrajectoryError::DerivedQoiMismatch(0)
    );
}

#[test]
fn z_axis_rigid_world_transforms_preserve_intrinsic_qois() {
    for translation in [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(3.0, -2.0, 0.5),
        Vec3::new(-1.0e3, 2.0e3, -4.0),
    ] {
        let original = sample(0.01, RenderSampleDisposition::HorizonCensored);
        let mut transformed = original.clone();
        transformed.center_of_mass_world_m = Vec3::new(
            -original.center_of_mass_world_m.x + translation.x,
            -original.center_of_mass_world_m.y + translation.y,
            original.center_of_mass_world_m.z + translation.z,
        );
        transformed.linear_momentum_world_kg_m_per_s = Vec3::new(
            -original.linear_momentum_world_kg_m_per_s.x,
            -original.linear_momentum_world_kg_m_per_s.y,
            original.linear_momentum_world_kg_m_per_s.z,
        );
        transformed.symmetry_axis_world = Vec3::new(
            -original.symmetry_axis_world.x,
            -original.symmetry_axis_world.y,
            original.symmetry_axis_world.z,
        );
        let q = original.orientation_body_to_world;
        // Left-compose a world-frame pi rotation about +z, then use its
        // canonical q/-q representative.
        transformed.orientation_body_to_world = [0.0, -q[2], q[1], q[0]];
        transformed.qois =
            DerivedEulerQois::from_state(state_from(&transformed), mass(), 0.0).unwrap();

        let original = RenderTrajectory::try_new(metadata(), vec![original]).unwrap();
        let transformed = RenderTrajectory::try_new(metadata(), vec![transformed]).unwrap();
        let original = original.samples()[0].input().qois;
        let transformed = transformed.samples()[0].input().qois;
        assert_close(original.inclination_rad, transformed.inclination_rad);
        assert_close(
            original.precession_rad_per_s,
            transformed.precession_rad_per_s,
        );
        assert_close(original.spin_rad_per_s, transformed.spin_rad_per_s);
        assert_close(
            original.precession_acceleration_rad_per_s2,
            transformed.precession_acceleration_rad_per_s2,
        );
    }
}

#[test]
fn unknown_schema_and_zero_component_identity_refuse() {
    let mut wrong_version = metadata();
    wrong_version.schema_version += 1;
    assert_eq!(
        RenderTrajectory::try_new(
            wrong_version,
            vec![sample(0.01, RenderSampleDisposition::HorizonCensored)]
        )
        .unwrap_err(),
        RenderTrajectoryError::UnsupportedSchemaVersion(2)
    );

    let mut zero_identity = metadata();
    zero_identity.model_identity = ContentHash([0; 32]);
    assert_eq!(
        RenderTrajectory::try_new(
            zero_identity,
            vec![sample(0.01, RenderSampleDisposition::HorizonCensored)]
        )
        .unwrap_err(),
        RenderTrajectoryError::ZeroIdentity("model_identity")
    );
}

#[test]
fn accepted_coupled_runner_output_converts_without_reconstructing_state() {
    let initial = CoupledInitialState {
        inclination_rad: 0.08,
        precession_rad_per_s: 16.0,
        spin_rad_per_s: 120.0,
    };
    let run = run_closed_reduced(
        runner_factors(),
        CoupledControls {
            timestep_s: 2.0e-5,
            maximum_steps: 20,
            terminal_inclination_rad: 0.002,
            reimpact_limit: 32,
        },
        initial,
        None,
    )
    .expect("accepted reduced run");
    let trajectory = RenderTrajectory::from_coupled_run(runner_metadata(&run), &run)
        .expect("runner trajectory admission");
    assert_eq!(trajectory.samples().len(), run.samples.len());
    for (retained, source) in trajectory.samples().iter().zip(&run.samples) {
        assert_eq!(retained.state(), source.state);
        assert_eq!(retained.input().time_s, source.time_s);
        assert_eq!(
            retained.input().base_mode,
            Some(RenderBaseModeState {
                displacement_m: source.base_deflection_m,
                velocity_m_per_s: source.base_velocity_m_per_s,
            })
        );
    }

    let mut wrong_metadata = runner_metadata(&run);
    wrong_metadata.timestep_s *= 2.0;
    assert_eq!(
        RenderTrajectory::from_coupled_run(wrong_metadata, &run).unwrap_err(),
        RenderTrajectoryError::RunnerConfigurationMismatch
    );

    let mut inconsistent_sample = run.clone();
    inconsistent_sample.samples[0]
        .center_of_mass_velocity_world_m_per_s
        .x += 1.0e-9;
    assert_eq!(
        RenderTrajectory::from_coupled_run(
            runner_metadata(&inconsistent_sample),
            &inconsistent_sample
        )
        .unwrap_err(),
        RenderTrajectoryError::RunnerSampleMismatch(0)
    );

    let mut inconsistent_checkpoint = run.clone();
    inconsistent_checkpoint.checkpoint.base_deflection_m += 1.0e-9;
    assert_eq!(
        RenderTrajectory::from_coupled_run(
            runner_metadata(&inconsistent_checkpoint),
            &inconsistent_checkpoint,
        )
        .unwrap_err(),
        RenderTrajectoryError::RunnerCheckpointMismatch
    );
}

#[test]
fn localized_terminal_bracket_from_runner_is_admitted_at_the_retained_root() {
    let run = run_closed_reduced(
        runner_factors(),
        CoupledControls {
            timestep_s: 1.0e-3,
            maximum_steps: 32,
            terminal_inclination_rad: 0.079,
            reimpact_limit: 8,
        },
        CoupledInitialState {
            inclination_rad: 0.08,
            precession_rad_per_s: 0.0,
            spin_rad_per_s: 0.0,
        },
        None,
    )
    .expect("terminal reduced run");
    let source = run.samples.last().expect("terminal sample");
    let event = source
        .terminal_inclination_event
        .expect("localized terminal bracket");
    assert_eq!(source.time_s, event.time_s);
    assert!(event.bracket_start_s <= event.time_s);
    assert!(event.time_s <= event.bracket_end_s);
    let trajectory = RenderTrajectory::from_coupled_run(runner_metadata(&run), &run)
        .expect("terminal runner trajectory admission");
    assert_eq!(
        trajectory.samples().last().unwrap().input().disposition,
        RenderSampleDisposition::TerminalInclination
    );
}

#[test]
fn prohibited_reimpact_root_is_a_complete_numerical_refusal_sample() {
    let controls = CoupledControls {
        timestep_s: 1.0e-4,
        maximum_steps: 100,
        terminal_inclination_rad: 0.002,
        reimpact_limit: 0,
    };
    let initial = CoupledInitialState {
        inclination_rad: 0.03,
        precession_rad_per_s: 0.0,
        spin_rad_per_s: 0.0,
    };
    let run = run_closed_reduced(runner_factors(), controls, initial, None)
        .expect("reimpact-refusal run");
    let trajectory = RenderTrajectory::from_coupled_run(runner_metadata(&run), &run)
        .expect("refusal-boundary trajectory admission");
    let final_sample = trajectory.samples().last().expect("refusal sample");
    assert_eq!(final_sample.input().time_s, run.checkpoint.time_s);
    assert_eq!(final_sample.state(), run.checkpoint.state);
    assert_eq!(
        final_sample.input().disposition,
        RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::ReimpactLimitExceeded
        )
    );
    assert_eq!(
        final_sample
            .input()
            .contact_transitions
            .last()
            .map(|event| event.kind),
        Some(ContactTransitionKind::Reimpact)
    );
    assert_eq!(
        final_sample
            .input()
            .contact_transitions
            .last()
            .map(|event| event.time_s.to_bits()),
        Some(final_sample.input().time_s.to_bits())
    );

    let resumed = run_closed_reduced(
        runner_factors(),
        controls,
        initial,
        Some(run.checkpoint.clone()),
    )
    .expect("refusal checkpoint restart");
    assert!(resumed.samples.is_empty());
    assert_eq!(resumed.checkpoint, run.checkpoint);
    assert_eq!(
        RenderTrajectory::from_coupled_run(runner_metadata(&resumed), &resumed).unwrap_err(),
        RenderTrajectoryError::EmptyTrajectory
    );
}

#[test]
fn thorne_tail_is_grounded_horizon_censored_and_contains_the_terminal_chirp() {
    with_cx(|cx| {
        let (run, trajectory) = thorne_trajectory(cx);
        let source_final = run.final_sample().expect("source cutoff");
        assert!(source_final.time_s > REDUCED_DECAY_RENDER_TAIL_HORIZON_S);
        let samples = trajectory.samples();
        let first = samples.first().expect("cinematic start").input();
        let final_sample = samples.last().expect("cinematic cutoff").input();
        assert_eq!(first.time_s.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            final_sample.time_s.to_bits(),
            REDUCED_DECAY_RENDER_TAIL_HORIZON_S.to_bits()
        );
        assert_eq!(
            final_sample.qois.inclination_rad.to_bits(),
            THORNE_2026_RENDER_VALIDITY_CUTOFF_THETA_RAD.to_bits()
        );
        assert_eq!(
            final_sample.disposition,
            RenderSampleDisposition::HorizonCensored
        );
        assert!(final_sample.terminal_event.is_none());

        let last_shutter_end_s = REDUCED_DECAY_RENDER_TAIL_HORIZON_S - 1.0 / 96.0;
        let visible_late_sample = samples
            .iter()
            .rev()
            .find(|sample| sample.input().time_s <= last_shutter_end_s)
            .expect("sample visible in the final shutter")
            .input();
        assert!(visible_late_sample.qois.inclination_rad < 0.006);
        assert!(visible_late_sample.qois.precession_rad_per_s > 400.0);

        for (index, retained) in samples.iter().enumerate() {
            let input = retained.input();
            assert_eq!(input.contact_branch, RenderContactBranch::Closed);
            let contact = input.contact_geometry.expect("kinematic ground support");
            assert!(contact.point_world_m.z.abs() < 1.0e-10);
            assert!(input.signed_gap_m <= 0.0);
            assert!(input.signed_gap_m.abs() < 1.0e-10);
            assert_eq!(input.interval_contact_active, index != 0);
            assert_eq!(input.interval_normal_force_n, 0.0);
            for channel in [input.channels.rolling, input.channels.gas] {
                assert_eq!(channel.force_world_n, Vec3::ZERO);
                assert_eq!(channel.torque_world_nm, Vec3::ZERO);
                assert!(channel.work_j <= 0.0);
            }
        }
        assert!(trajectory.metadata().no_claims.iter().any(|claim| {
            claim.contains("not a raw measured trajectory")
                && !claim.contains("validated raw trajectory")
        }));
        assert!(
            trajectory
                .metadata()
                .no_claims
                .iter()
                .any(|claim| claim.contains("full-FSI"))
        );
        assert!(
            trajectory
                .metadata()
                .no_claims
                .iter()
                .any(|claim| claim.contains("no theta-zero"))
        );
        assert!(
            trajectory
                .metadata()
                .no_claims
                .iter()
                .any(|claim| claim.contains("not a configuration") && claim.contains("ranking"))
        );
    });
}

#[test]
fn thorne_bridge_refuses_an_intermediate_state_that_does_not_follow_the_decay_step() {
    with_cx(|cx| {
        let benchmark = Thorne2026SteelGlassBenchmark::ambient().expect("ambient benchmark");
        let profile = benchmark
            .resolve_specimen(cx)
            .expect("source-bound specimen");
        let mut run =
            run_thorne_2026_steel_glass_benchmark(&benchmark).expect("source-bound decay");
        let index = run.samples.len() / 2;
        assert!(
            index >= 313,
            "the source run must exercise the accumulated-clock horizon that previously rejected sample 313"
        );
        RenderTrajectory::from_reduced_decay_run(&run, &profile, cx)
            .expect("the unmodified source run must admit across the accumulated-clock horizon");
        let previous_theta_rad = run.samples[index - 1].theta_rad;
        let original_theta_rad = run.samples[index].theta_rad;
        let corrupted_theta_rad = 0.5 * (previous_theta_rad + original_theta_rad);
        let sample = &mut run.samples[index];
        let original_omega_rad_s = sample.omega_rad_s;
        sample.theta_rad = corrupted_theta_rad;
        sample.omega_rad_s = (4.0 * run.parameters.gravity_m_per_s2
            / (run.parameters.radius_m * corrupted_theta_rad.sin()))
        .sqrt();
        sample.energy_j = 1.5
            * run.parameters.mass_kg
            * run.parameters.gravity_m_per_s2
            * run.parameters.radius_m
            * corrupted_theta_rad;
        sample.powers.published_rolling_w *= (sample.omega_rad_s * corrupted_theta_rad.cos())
            / (original_omega_rad_s * original_theta_rad.cos());
        sample.powers.bildsten_boundary_layer_w *=
            (corrupted_theta_rad / original_theta_rad).powf(-1.25);

        let error = RenderTrajectory::from_reduced_decay_run(&run, &profile, cx)
            .expect_err("a forged intermediate inclination must be refused");
        assert!(matches!(
            error,
            RenderTrajectoryError::ReducedDecayBridgeRefusal {
                field: "sample.integration_step",
                ..
            }
        ));
    });
}

#[test]
fn thorne_bridge_replays_exactly_and_phase_rates_match_rigid_state() {
    with_cx(|cx| {
        let benchmark = Thorne2026SteelGlassBenchmark::ambient().expect("ambient benchmark");
        let profile = benchmark
            .resolve_specimen(cx)
            .expect("source-bound specimen");
        let run = run_thorne_2026_steel_glass_benchmark(&benchmark).expect("source-bound decay");
        let first = RenderTrajectory::from_reduced_decay_run(&run, &profile, cx)
            .expect("first deterministic bridge");
        let replay = RenderTrajectory::from_reduced_decay_run(&run, &profile, cx)
            .expect("replayed deterministic bridge");
        assert_eq!(first, replay);

        let mass = first.metadata().mass_properties.properties;
        let mut precession_phase_rad = 0.0_f64;
        let mut spin_phase_rad = 0.0_f64;
        let mut previous: Option<&RenderTrajectorySampleInput> = None;
        for (index, retained) in first.samples().iter().enumerate() {
            let input = retained.input();
            if let Some(previous_input) = previous {
                let dt_s = input.time_s - previous_input.time_s;
                precession_phase_rad = (precession_phase_rad
                    + 0.5
                        * (previous_input.qois.precession_rad_per_s
                            + input.qois.precession_rad_per_s)
                        * dt_s)
                    .rem_euclid(TAU);
                spin_phase_rad = (spin_phase_rad
                    + 0.5
                        * (previous_input.qois.spin_rad_per_s + input.qois.spin_rad_per_s)
                        * dt_s)
                    .rem_euclid(TAU);
            }
            let expected_orientation =
                UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), precession_phase_rad)
                    .and_then(|orientation| {
                        orientation.right_exp(Vec3::new(0.0, input.qois.inclination_rad, 0.0))
                    })
                    .and_then(|orientation| {
                        orientation.right_exp(Vec3::new(0.0, 0.0, spin_phase_rad))
                    })
                    .expect("independent Euler composition");
            if index % 10_000 == 0 || index + 1 == first.samples().len() {
                for (actual, expected) in retained
                    .state()
                    .pose()
                    .orientation()
                    .components()
                    .into_iter()
                    .zip(expected_orientation.components())
                {
                    assert_close(actual, expected);
                }
            }

            let orientation = retained.state().pose().orientation();
            let angular_velocity_world = orientation.rotate_body_to_world(
                mass.angular_velocity_body(retained.state().angular_momentum_body()),
            );
            let symmetry_axis_world = input.symmetry_axis_world;
            let tilt_axis_world =
                Vec3::new(-precession_phase_rad.sin(), precession_phase_rad.cos(), 0.0);
            let residual = angular_velocity_world
                .sub(Vec3::new(0.0, 0.0, input.qois.precession_rad_per_s))
                .sub(symmetry_axis_world.scale(input.qois.spin_rad_per_s));
            let theta_rate_from_state = residual.dot(tilt_axis_world);
            let theta_rate_from_qoi = -2.0
                * input.qois.precession_acceleration_rad_per_s2
                * input.qois.inclination_rad.tan()
                / input.qois.precession_rad_per_s;
            assert_close(theta_rate_from_state, theta_rate_from_qoi);
            assert!(
                residual
                    .sub(tilt_axis_world.scale(theta_rate_from_state))
                    .norm_squared()
                    < 1.0e-18
            );

            let contact = input.contact_geometry.expect("closed support");
            let radius_body = retained
                .state()
                .pose()
                .point_body_from_world(contact.point_world_m)
                .expect("body support point");
            let body_contact_azimuth = radius_body.y.atan2(radius_body.x);
            let expected_body_contact_azimuth = -spin_phase_rad;
            let phase_error =
                (body_contact_azimuth - expected_body_contact_azimuth + PI).rem_euclid(TAU) - PI;
            assert!(phase_error.abs() < 1.0e-9);
            previous = Some(input);
        }
    });
}
