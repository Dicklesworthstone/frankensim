//! G0/G3 checks for animation-grade Euler trajectory semantics.

use fs_blake3::{ContentHash, hash_domain};
use fs_euler_disc_e2e::coupled_runner::{ChannelOwnership, ContactTransitionKind};
use fs_euler_disc_e2e::{
    DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, RenderBaseModeState,
    RenderContactBranch, RenderContactGeometry, RenderContactTransition, RenderMassProperties,
    RenderNumericalRefusalReason, RenderSampleDisposition, RenderSupportFeature,
    RenderTerminalEvent, RenderTrajectory, RenderTrajectoryAuthority, RenderTrajectoryError,
    RenderTrajectoryMetadata, RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const THETA: f64 = 0.35;

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
        model_identity: identity("model"),
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

fn two_samples() -> Vec<RenderTrajectorySampleInput> {
    vec![
        sample(0.01, RenderSampleDisposition::Continue),
        sample(0.02, RenderSampleDisposition::HorizonCensored),
    ]
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
    inputs[1].contact_transitions = vec![RenderContactTransition {
        kind: ContactTransitionKind::Opening,
        time_s: 0.015,
        bracket_start_s: 0.014,
        bracket_end_s: 0.016,
    }];
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
    RenderTrajectory::try_new(metadata(), vec![terminal]).unwrap();

    let invalid_refusal = sample(
        0.01,
        RenderSampleDisposition::NumericalRefusal(RenderNumericalRefusalReason::BackendSpecific(0)),
    );
    assert_eq!(
        RenderTrajectory::try_new(metadata(), vec![invalid_refusal]).unwrap_err(),
        RenderTrajectoryError::InvalidNumericalRefusalCode(0)
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
