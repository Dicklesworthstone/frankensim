//! G0/G3 tests for the atomic Euler mechanics composition rung.

#[path = "../src/mechanics.rs"]
mod mechanics;
#[path = "../src/ports.rs"]
mod ports;

use fs_couple::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp, StableId};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use mechanics::{
    BodyAction, EnergyAcceptanceBound, InactiveMechanicsChannel, InactiveMechanicsChannels,
    MechanicsContribution, MechanicsMacroStepError, MechanicsMacroStepInput, MechanicsReference,
    run_mechanics_macro_step,
};
use ports::{
    ChannelActivity, ContributionDomain, ContributionOwnership, DecompositionReceipt, EulerChannel,
    GeneralizedVelocityCoordinate, PatchRegion, PortDeclaration, PortInterval, SurfacePair,
};

fn stable(value: &str) -> StableId {
    StableId::new(value).expect("canonical test identity")
}

fn timestamp(tick: u64) -> PortTimestamp {
    PortTimestamp::new(stable("mechanics-clock"), tick)
}

fn domain(identity: &str, frame: &str, orientation: PortOrientation) -> ContributionDomain {
    ContributionDomain::new(
        SurfacePair::try_new(stable("disc"), stable("ground")).expect("distinct surfaces"),
        PatchRegion::try_new(stable(identity), 0, 1).expect("nonempty patch"),
        PortInterval::try_new(timestamp(0), timestamp(1)).expect("nonempty interval"),
        GeneralizedVelocityCoordinate::new(
            stable(&format!("coordinate-{identity}")),
            CoordinateBinding::new(stable("world-basis"), stable(frame), orientation),
        ),
    )
}

fn port(identity: &str, channel: EulerChannel) -> PortDeclaration {
    PortDeclaration::new(
        stable(identity),
        channel,
        if channel == EulerChannel::RollingContourSpin {
            PortKind::RotationalTorqueAngularVelocity
        } else {
            PortKind::MechanicalForceVelocity
        },
        ChannelActivity::Active,
        stable(&format!("law-{identity}")),
        stable(&format!("source-{identity}")),
        domain(identity, "world", PortOrientation::AlongFrame),
        ContributionOwnership::Exclusive,
    )
}

fn unavailable(channel: EulerChannel, name: &str) -> InactiveMechanicsChannel {
    InactiveMechanicsChannel {
        channel,
        activity: ChannelActivity::Unavailable {
            model_identity: stable(&format!("model-{name}")),
            reason_identity: stable(&format!("reason-{name}")),
        },
    }
}

fn inactive_channels() -> InactiveMechanicsChannels {
    InactiveMechanicsChannels {
        base: unavailable(EulerChannel::Base, "base"),
        external_gas: unavailable(EulerChannel::ExternalGas, "external-gas"),
        gas_film: unavailable(EulerChannel::GasFilm, "gas-film"),
    }
}

fn properties() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).expect("valid properties")
}

fn state(orientation: UnitQuaternion) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(Vec3::ZERO, orientation).expect("valid pose"),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .expect("valid state")
}

fn contribution(identity: &str, channel: EulerChannel) -> MechanicsContribution {
    MechanicsContribution {
        port: port(identity, channel),
        body_action: BodyAction::ActionOnIntegratedBody,
        reference: MechanicsReference::CenterOfMass,
        application_arm_body_m: Vec3::ZERO,
        force_world_n: Vec3::ZERO,
        free_torque_world_nm: Vec3::ZERO,
        claimed_discrete_work_j: 0.0,
        claimed_storage_j: 0.0,
        claimed_dissipation_j: 0.0,
        claimed_heat_j: 0.0,
        uncertainty_j: 0.0,
    }
}

fn input(contributions: Vec<MechanicsContribution>) -> MechanicsMacroStepInput {
    MechanicsMacroStepInput {
        state: state(UnitQuaternion::IDENTITY),
        mass_properties: properties(),
        gravity: Gravity::ZERO,
        duration_seconds: 0.1,
        expected_contribution_keys: contributions
            .iter()
            .map(|contribution| contribution.port.identity().clone())
            .collect(),
        contributions,
        inactive_channels: inactive_channels(),
        world_frame: stable("world"),
        energy_bound: EnergyAcceptanceBound {
            absolute_j: 1.0e-12,
            relative: 1.0e-12,
        },
        cancelled_before_step: false,
    }
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:?}"
    );
}

#[test]
fn force_arm_is_recentered_in_world_then_rotated_once_to_body() {
    let rotation =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), core::f64::consts::FRAC_PI_2)
            .expect("rotation");
    let mut force = contribution("normal-arm", EulerChannel::NormalContact);
    force.application_arm_body_m = Vec3::new(1.0, 0.0, 0.0);
    force.force_world_n = Vec3::new(2.0, 0.0, 0.0);
    force.free_torque_world_nm = Vec3::new(0.0, 0.0, 0.5);
    let mut request = input(vec![force]);
    request.state = state(rotation);
    request.energy_bound.absolute_j = 10.0;

    let receipt = run_mechanics_macro_step(request).expect("admitted force-arm step");
    assert_close(receipt.resultant.force_world_n.x, 2.0, 1.0e-14);
    assert_close(
        receipt.resultant.torque_center_of_mass_world_nm.z,
        -1.5,
        1.0e-12,
    );
    assert_close(
        receipt.resultant.torque_center_of_mass_body_nm.z,
        -1.5,
        1.0e-12,
    );
}

#[test]
fn action_reaction_fixture_admits_only_the_action_on_this_body() {
    let mut action = contribution("normal-action", EulerChannel::NormalContact);
    action.force_world_n = Vec3::new(4.0, 0.0, 0.0);
    action.claimed_discrete_work_j = 0.04;
    action.claimed_storage_j = -0.04;
    let receipt =
        run_mechanics_macro_step(input(vec![action.clone()])).expect("action is admitted");
    assert_close(
        receipt
            .rigid_body_step
            .state_after
            .linear_momentum_world()
            .x,
        0.4,
        1.0e-14,
    );

    action.body_action = BodyAction::ReactionOnOtherBody;
    let error =
        run_mechanics_macro_step(input(vec![action])).expect_err("reaction must not double count");
    assert!(matches!(
        error,
        MechanicsMacroStepError::ReactionCannotActOnIntegratedBody { .. }
    ));
}

#[test]
fn expected_keys_are_complete_and_exactly_once() {
    let normal = contribution("normal-key", EulerChannel::NormalContact);
    let tangent = contribution("tangent-key", EulerChannel::TangentialContact);
    let mut missing = input(vec![normal.clone()]);
    missing
        .expected_contribution_keys
        .push(tangent.port.identity().clone());
    assert!(matches!(
        run_mechanics_macro_step(missing),
        Err(MechanicsMacroStepError::MissingContribution { .. })
    ));

    let mut duplicate = input(vec![normal.clone(), normal]);
    duplicate.expected_contribution_keys = vec![stable("normal-key")];
    assert!(matches!(
        run_mechanics_macro_step(duplicate),
        Err(MechanicsMacroStepError::DuplicateContribution { .. })
    ));

    let mut extra = input(vec![tangent]);
    extra.expected_contribution_keys = vec![stable("normal-key")];
    assert!(matches!(
        run_mechanics_macro_step(extra),
        Err(MechanicsMacroStepError::UnexpectedContribution { .. })
    ));
}

#[test]
fn gravity_and_zero_wrench_preserve_mechanical_energy() {
    let mut request = input(Vec::new());
    request.gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.81)).expect("finite gravity");
    request.state = RigidBodyState::new(
        Pose::new(Vec3::new(0.0, 0.0, 2.0), UnitQuaternion::IDENTITY).expect("pose"),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::ZERO,
    )
    .expect("state");
    let receipt = run_mechanics_macro_step(request).expect("gravity-only step");
    assert_close(receipt.energy.mechanical_energy_change_j, 0.0, 1.0e-12);
    assert_close(receipt.energy.body_energy_residual_j, 0.0, 1.0e-12);
    assert_close(receipt.energy.interface_energy_residual_j, 0.0, 1.0e-12);
    assert_close(receipt.energy.combined_energy_residual_j, 0.0, 1.0e-12);

    let mut malformed_heat = contribution("tangent-heat-bound", EulerChannel::TangentialContact);
    malformed_heat.claimed_dissipation_j = 0.35;
    malformed_heat.claimed_heat_j = 0.36;
    assert!(matches!(
        run_mechanics_macro_step(input(vec![malformed_heat])),
        Err(MechanicsMacroStepError::InvalidInput(
            "mechanics contribution"
        ))
    ));
}

#[test]
fn conservative_contact_storage_release_closes_all_three_identities() {
    let mut manufactured = contribution("impact-manufactured", EulerChannel::Impact);
    manufactured.force_world_n = Vec3::new(4.0, 0.0, 0.0);
    manufactured.free_torque_world_nm = Vec3::new(0.0, 0.0, 3.0);
    // At rest, midpoint translation work is F² dt²/(2m) = 0.04 J and
    // midpoint spin work is tau² dt²/(2I) = 0.045 J.
    manufactured.claimed_discrete_work_j = 0.085;
    manufactured.claimed_storage_j = -0.085;
    let receipt = run_mechanics_macro_step(input(vec![manufactured])).expect("manufactured step");
    assert_close(receipt.energy.reconstructed_body_work_j, 0.085, 1.0e-12);
    assert_close(receipt.energy.mechanical_energy_change_j, 0.085, 1.0e-12);
    assert_close(receipt.energy.recoverable_storage_change_j, -0.085, 1.0e-12);
    assert_close(receipt.energy.irreversible_work_j, 0.0, 1.0e-12);
    assert_close(receipt.energy.heat_j, 0.0, 1.0e-12);
    assert_close(receipt.energy.body_energy_residual_j, 0.0, 1.0e-12);
    assert_close(receipt.energy.interface_energy_residual_j, 0.0, 1.0e-12);
    assert_close(receipt.energy.combined_energy_residual_j, 0.0, 1.0e-12);
}

#[test]
fn unloading_permits_negative_recoverable_storage_change() {
    let mut unloading = contribution("normal-unloading", EulerChannel::NormalContact);
    unloading.force_world_n = Vec3::new(4.0, 0.0, 0.0);
    unloading.claimed_discrete_work_j = 0.04;
    unloading.claimed_storage_j = -0.04;

    let receipt = run_mechanics_macro_step(input(vec![unloading])).expect("unloading is admitted");
    assert_close(receipt.energy.recoverable_storage_change_j, -0.04, 1.0e-12);
    assert_close(receipt.energy.interface_energy_residual_j, 0.0, 1.0e-12);
}

#[test]
fn dissipative_heat_equals_irreversible_work_without_double_counting() {
    let mut dissipative = contribution("tangent-dissipative", EulerChannel::TangentialContact);
    dissipative.force_world_n = Vec3::new(-4.0, 0.0, 0.0);
    dissipative.claimed_discrete_work_j = -0.36;
    dissipative.claimed_dissipation_j = 0.36;
    dissipative.claimed_heat_j = 0.36;
    let mut request = input(vec![dissipative]);
    request.state = RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).expect("pose"),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::ZERO,
    )
    .expect("state");

    let receipt = run_mechanics_macro_step(request).expect("dissipative step");
    assert_close(receipt.energy.mechanical_energy_change_j, -0.36, 1.0e-12);
    assert_close(receipt.energy.irreversible_work_j, 0.36, 1.0e-12);
    assert_close(receipt.energy.heat_j, 0.36, 1.0e-12);
    assert_close(receipt.energy.body_energy_residual_j, 0.0, 1.0e-12);
    assert_close(receipt.energy.interface_energy_residual_j, 0.0, 1.0e-12);
    assert_close(receipt.energy.combined_energy_residual_j, 0.0, 1.0e-12);
}

#[test]
fn force_only_refinement_is_state_and_energy_consistent() {
    let mut full_force = contribution("normal-refine", EulerChannel::NormalContact);
    full_force.force_world_n = Vec3::new(4.0, 0.0, 0.0);
    full_force.claimed_discrete_work_j = 0.04;
    full_force.claimed_storage_j = -0.04;
    let full = run_mechanics_macro_step(input(vec![full_force])).expect("full step");

    let mut half_one_force = contribution("normal-refine", EulerChannel::NormalContact);
    half_one_force.force_world_n = Vec3::new(4.0, 0.0, 0.0);
    half_one_force.claimed_discrete_work_j = 0.01;
    half_one_force.claimed_storage_j = -0.01;
    let mut half_one = input(vec![half_one_force]);
    half_one.duration_seconds = 0.05;
    let half_one = run_mechanics_macro_step(half_one).expect("first half");
    let mut half_two_force = contribution("normal-refine", EulerChannel::NormalContact);
    half_two_force.force_world_n = Vec3::new(4.0, 0.0, 0.0);
    half_two_force.claimed_discrete_work_j = 0.03;
    half_two_force.claimed_storage_j = -0.03;
    let mut half_two = input(vec![half_two_force]);
    half_two.duration_seconds = 0.05;
    half_two.state = half_one.rigid_body_step.state_after;
    let half_two = run_mechanics_macro_step(half_two).expect("second half");

    assert_close(
        full.rigid_body_step.state_after.pose().position_world().x,
        half_two
            .rigid_body_step
            .state_after
            .pose()
            .position_world()
            .x,
        1.0e-14,
    );
    assert_close(
        full.rigid_body_step.state_after.linear_momentum_world().x,
        half_two
            .rigid_body_step
            .state_after
            .linear_momentum_world()
            .x,
        1.0e-14,
    );
    assert!(
        half_two.energy.combined_energy_residual_j.abs()
            <= full.energy.combined_energy_residual_j.abs() + 1.0e-14
    );
}

#[test]
fn all_supported_channels_are_individually_ablated_and_base_gas_stay_inactive() {
    for (index, channel) in [
        EulerChannel::NormalContact,
        EulerChannel::TangentialContact,
        EulerChannel::RollingContourSpin,
        EulerChannel::Impact,
    ]
    .into_iter()
    .enumerate()
    {
        let receipt = run_mechanics_macro_step(input(vec![contribution(
            &format!("channel-{index}"),
            channel,
        )]))
        .expect("each supported channel can be independently active");
        assert_close(receipt.energy.reconstructed_body_work_j, 0.0, 1.0e-14);
    }
    let mut request = input(Vec::new());
    request.inactive_channels.base.activity = ChannelActivity::Active;
    assert!(matches!(
        run_mechanics_macro_step(request),
        Err(MechanicsMacroStepError::UnsupportedActiveChannel {
            channel: EulerChannel::Base
        })
    ));
}

#[test]
fn rigid_transform_and_cancellation_are_no_hidden_publication_paths() {
    let rotation =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.7).expect("rotation");
    let force_body = Vec3::new(3.0, -2.0, 0.0);
    let mut first = contribution("tangent-rigid", EulerChannel::TangentialContact);
    first.application_arm_body_m = Vec3::new(0.2, 0.4, 0.0);
    first.force_world_n = force_body;
    let mut plain = input(vec![first]);
    plain.energy_bound.absolute_j = 10.0;
    let plain = run_mechanics_macro_step(plain).expect("plain frame");

    let mut second = contribution("tangent-rigid", EulerChannel::TangentialContact);
    second.application_arm_body_m = Vec3::new(0.2, 0.4, 0.0);
    second.force_world_n = rotation.rotate_body_to_world(force_body);
    let mut rotated = input(vec![second]);
    rotated.state = state(rotation);
    rotated.energy_bound.absolute_j = 10.0;
    let rotated = run_mechanics_macro_step(rotated).expect("rotated frame");
    let expected_rotated_force = rotation.rotate_body_to_world(plain.resultant.force_world_n);
    assert_close(
        rotated.resultant.force_world_n.x,
        expected_rotated_force.x,
        1.0e-12,
    );
    assert_close(
        rotated.resultant.force_world_n.y,
        expected_rotated_force.y,
        1.0e-12,
    );

    let mut cancelled = input(Vec::new());
    cancelled.cancelled_before_step = true;
    assert!(matches!(
        run_mechanics_macro_step(cancelled),
        Err(MechanicsMacroStepError::CancelledBeforeStep)
    ));
}

#[test]
fn malformed_additive_ownership_refuses_before_integration() {
    let identity = stable("additive-ownership");
    let shared_domain = domain("additive-domain", "world", PortOrientation::AlongFrame);
    let decomposition = DecompositionReceipt::try_new(
        stable("additive-receipt"),
        shared_domain.clone(),
        [identity.clone(), stable("other-contributor")],
    )
    .expect("well-formed structural receipt");
    let mut malformed = contribution(identity.as_str(), EulerChannel::NormalContact);
    malformed.port = PortDeclaration::new(
        identity,
        EulerChannel::NormalContact,
        PortKind::MechanicalForceVelocity,
        ChannelActivity::Active,
        stable("law-additive-ownership"),
        stable("source-additive-ownership"),
        shared_domain,
        ContributionOwnership::AdditiveWithProof {
            decomposition_receipt: decomposition,
        },
    );

    assert!(matches!(
        run_mechanics_macro_step(input(vec![malformed])),
        Err(MechanicsMacroStepError::OverlappingOwnership { .. })
    ));
}

#[test]
fn deterministic_replay_publishes_identical_energy_receipts() {
    let mut conservative = contribution("replay-conservative", EulerChannel::NormalContact);
    conservative.force_world_n = Vec3::new(4.0, 0.0, 0.0);
    conservative.claimed_discrete_work_j = 0.04;
    conservative.claimed_storage_j = -0.04;
    let request = input(vec![conservative]);

    let first = run_mechanics_macro_step(request.clone()).expect("first replay");
    let second = run_mechanics_macro_step(request).expect("second replay");
    assert_eq!(first, second);
}
