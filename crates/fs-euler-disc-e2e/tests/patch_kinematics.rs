//! G0/G3 coverage for the unintegrated pre-constitutive patch kinematics leaf.

#[path = "../src/contact_dynamics.rs"]
mod contact_dynamics;
#[path = "../src/patch_kinematics.rs"]
mod patch_kinematics;

use fs_couple::StableId;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_rep_frep::AxisymmetricSupportAuthority;
use fs_tribo::InputAuthority;
use patch_kinematics::{
    Creepage, CurvatureMetadata, MovingOneModeBaseState, MovingOneModePatchBridgeInput,
    OrderedSurfacePair, PatchContactStatus, PatchGeometryMetadata, PatchKinematicThresholds,
    PatchKinematicsInput, ProfileSupportKinematics, SurfaceOrder, TangentGaugeInput,
    TangentGaugeSource, bridge_moving_one_mode_patch_kinematics, compute_patch_kinematics,
};

fn id(value: &str) -> StableId {
    StableId::new(value).expect("test identity")
}

fn properties() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).expect("valid inertia")
}

fn state(position: Vec3, velocity: Vec3, angular_velocity: Vec3) -> RigidBodyState {
    let props = properties();
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).expect("finite pose"),
        velocity.scale(props.mass()),
        angular_velocity,
    )
    .expect("finite state")
}

fn thresholds() -> PatchKinematicThresholds {
    PatchKinematicThresholds {
        threshold_identity: id("thresholds-v1"),
        tie_break_identity: id("ties-v1"),
        separation_gap_m: 1.0e-2,
        touching_gap_m: 1.0e-4,
        support_point_coincidence_tolerance_m: 1.0e-9,
        tangent_counterpart_coincidence_tolerance_m: 1.0e-9,
        approach_speed_m_per_s: 1.0e-2,
        impact_candidate_speed_m_per_s: 1.0,
        stationary_tangent_speed_m_per_s: 1.0e-9,
        minimum_reference_rolling_speed_m_per_s: 1.0e-6,
        gauge_degeneracy_norm: 1.0e-12,
    }
}

fn input(
    disc_velocity: Vec3,
    disc_angular_velocity: Vec3,
    base_velocity: Vec3,
) -> PatchKinematicsInput {
    PatchKinematicsInput {
        surfaces: OrderedSurfacePair::try_new(id("disc"), id("base"), SurfaceOrder::DiscThenBase)
            .expect("ordered surfaces"),
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        profile_support: ProfileSupportKinematics {
            disc_arm_world_m: Vec3::new(0.0, 0.0, -1.0),
            disc_point_world_m: Vec3::ZERO,
            gap_m: 0.0,
            source_feature: 7,
            support_authority: AxisymmetricSupportAuthority::Estimate,
        },
        patch: PatchGeometryMetadata {
            patch_identity: id("patch-rim"),
            source_feature: 7,
            gap_uncertainty_m: 0.0,
            curvature: CurvatureMetadata::Known {
                curvature_identity: id("curvature-rim"),
                authority: InputAuthority::SyntheticFixture,
                first_principal_m_inverse: 1.0,
                second_principal_m_inverse: 0.0,
                uncertainty_m_inverse: 1.0e-3,
            },
        },
        disc_state: state(
            Vec3::new(0.0, 0.0, 1.0),
            disc_velocity,
            disc_angular_velocity,
        ),
        disc_mass_properties: properties(),
        base_state: state(Vec3::ZERO, base_velocity, Vec3::ZERO),
        base_mass_properties: properties(),
        base_contact_arm_body_m: Vec3::ZERO,
        tangent_gauge: TangentGaugeInput {
            reference_world: Vec3::new(1.0, 0.0, 0.0),
            rotation_rad: 0.0,
        },
        thresholds: thresholds(),
        tangent_effort_probe_world_n: Some(Vec3::new(3.0, -2.0, 0.0)),
    }
}

fn assert_close(left: f64, right: f64) {
    assert!((left - right).abs() <= 1.0e-12, "{left} != {right}");
}

fn assert_vec_close(left: Vec3, right: Vec3) {
    assert_close(left.x, right.x);
    assert_close(left.y, right.y);
    assert_close(left.z, right.z);
}

#[test]
fn moving_one_mode_bridge_uses_disc_material_motion_without_faking_a_base_body() {
    let bridge = bridge_moving_one_mode_patch_kinematics(MovingOneModePatchBridgeInput {
        profile_support: ProfileSupportKinematics {
            disc_arm_world_m: Vec3::new(0.0, 0.0, -1.0),
            disc_point_world_m: Vec3::ZERO,
            gap_m: 0.0,
            source_feature: 7,
            support_authority: AxisymmetricSupportAuthority::Estimate,
        },
        disc_state: state(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
        ),
        disc_mass_properties: properties(),
        base_mode: MovingOneModeBaseState {
            undeformed_contact_point_world_m: Vec3::ZERO,
            vertical_displacement_m: 0.25,
            vertical_velocity_m_per_s: -0.5,
        },
        normal_world: Vec3::new(0.0, 0.0, 4.0),
        tangent_gauge: TangentGaugeInput {
            reference_world: Vec3::new(1.0, 0.0, 0.0),
            rotation_rad: 0.0,
        },
        thresholds: thresholds(),
    })
    .expect("moving-one-mode bridge");
    // omega x r cancels the disc COM velocity at the selected material point.
    assert_vec_close(bridge.disc_point.point_velocity_world, Vec3::ZERO);
    assert_vec_close(bridge.base_contact_point_world_m, Vec3::new(0.0, 0.0, 0.25));
    assert_vec_close(
        bridge.base_contact_velocity_world_m_per_s,
        Vec3::new(0.0, 0.0, -0.5),
    );
    assert_vec_close(bridge.normal_world, Vec3::new(0.0, 0.0, 1.0));
    assert_close(bridge.normal_gap_m, -0.25);
    assert_close(bridge.tangent_counterpart_residual_m, 0.0);
    assert_vec_close(
        bridge.relative_velocity_world_m_per_s,
        Vec3::new(0.0, 0.0, 0.5),
    );
    assert_close(bridge.normal_relative_velocity_m_per_s, 0.5);
    assert_vec_close(
        bridge.tangential_relative_velocity_world_m_per_s,
        Vec3::ZERO,
    );
    assert_vec_close(
        bridge
            .tangent_basis
            .first_world
            .cross(bridge.tangent_basis.second_world),
        bridge.normal_world,
    );
}

#[test]
fn moving_one_mode_bridge_refuses_a_tangentially_unrelated_base_point() {
    let result = bridge_moving_one_mode_patch_kinematics(MovingOneModePatchBridgeInput {
        profile_support: ProfileSupportKinematics {
            disc_arm_world_m: Vec3::new(0.0, 0.0, -1.0),
            disc_point_world_m: Vec3::ZERO,
            gap_m: 0.0,
            source_feature: 7,
            support_authority: AxisymmetricSupportAuthority::Estimate,
        },
        disc_state: state(Vec3::new(0.0, 0.0, 1.0), Vec3::ZERO, Vec3::ZERO),
        disc_mass_properties: properties(),
        base_mode: MovingOneModeBaseState {
            undeformed_contact_point_world_m: Vec3::new(1.0e-2, 0.0, 0.0),
            vertical_displacement_m: 0.0,
            vertical_velocity_m_per_s: 0.0,
        },
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_gauge: TangentGaugeInput {
            reference_world: Vec3::new(1.0, 0.0, 0.0),
            rotation_rad: 0.0,
        },
        thresholds: thresholds(),
    });
    assert!(matches!(
        result,
        Err(patch_kinematics::PatchKinematicsError::MovingOneModeCounterpartMismatch { .. })
    ));
}

#[test]
fn independent_point_velocity_reconstruction_and_pure_roll_hold() {
    // At r=(0,0,-1), omega=(0,-1,0) produces omega x r=(1,0,0),
    // cancelling the center-of-mass velocity. This is a kinematic fixture,
    // not a no-slip/contact-law conclusion.
    let result = compute_patch_kinematics(input(
        Vec3::new(-1.0, 0.0, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        Vec3::ZERO,
    ))
    .expect("pure roll kinematics");
    let explicit = result.disc_point.center_of_mass_velocity_world.add(
        result
            .disc_point
            .angular_velocity_world
            .cross(result.disc_point.arm_world),
    );
    assert_vec_close(result.disc_point.point_velocity_world, explicit);
    assert_vec_close(result.disc_point.point_velocity_world, Vec3::ZERO);
    assert_eq!(result.status, PatchContactStatus::Touching);
    assert!(matches!(result.creepage, Creepage::Unavailable { .. }));
}

#[test]
fn pure_spin_mixed_creepage_and_base_motion_are_distinguished_kinematically() {
    let spin = compute_patch_kinematics(input(Vec3::ZERO, Vec3::new(0.0, 0.0, 2.0), Vec3::ZERO))
        .expect("pure spin kinematics");
    assert_close(spin.normal_spin_rad_per_s, 2.0);
    assert_eq!(spin.status, PatchContactStatus::Touching);

    let mixed = compute_patch_kinematics(input(
        Vec3::new(4.0, 3.0, 0.0),
        Vec3::ZERO,
        Vec3::new(2.0, 0.0, 0.0),
    ))
    .expect("mixed/base-motion kinematics");
    assert_vec_close(
        mixed.relative_velocity_world_m_per_s,
        Vec3::new(2.0, 3.0, 0.0),
    );
    assert_eq!(mixed.status, PatchContactStatus::Grazing);
    assert!(matches!(mixed.creepage, Creepage::Available { .. }));
    assert_close(
        mixed.reference_rolling_speed_m_per_s,
        (3.0_f64 * 3.0 + 1.5 * 1.5).sqrt(),
    );
}

#[test]
fn separated_approaching_grazing_and_impact_candidates_use_boundaries() {
    let mut separated = input(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    separated.profile_support.gap_m = 0.02;
    separated.profile_support.disc_point_world_m.z = 0.02;
    separated.disc_state = state(Vec3::new(0.0, 0.0, 1.02), Vec3::ZERO, Vec3::ZERO);
    assert_eq!(
        compute_patch_kinematics(separated)
            .expect("separated")
            .status,
        PatchContactStatus::Separated
    );

    let approaching =
        compute_patch_kinematics(input(Vec3::new(0.0, 0.0, -0.2), Vec3::ZERO, Vec3::ZERO))
            .expect("approach");
    assert_eq!(approaching.status, PatchContactStatus::Approaching);

    let impact = compute_patch_kinematics(input(Vec3::new(0.0, 0.0, -1.0), Vec3::ZERO, Vec3::ZERO))
        .expect("impact threshold tie");
    assert_eq!(impact.status, PatchContactStatus::ImpactCandidate);

    let grazing = compute_patch_kinematics(input(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, Vec3::ZERO))
        .expect("grazing");
    assert_eq!(grazing.status, PatchContactStatus::Grazing);
}

#[test]
fn swapped_surface_order_reverses_relative_vector_without_changing_opening_sign() {
    let direct = compute_patch_kinematics(input(Vec3::new(2.0, 0.0, -0.2), Vec3::ZERO, Vec3::ZERO))
        .expect("direct order");
    let mut swapped = input(Vec3::new(2.0, 0.0, -0.2), Vec3::ZERO, Vec3::ZERO);
    swapped.surfaces =
        OrderedSurfacePair::try_new(id("base"), id("disc"), SurfaceOrder::BaseThenDisc)
            .expect("swapped order");
    swapped.normal_world = Vec3::new(0.0, 0.0, -1.0);
    let swapped = compute_patch_kinematics(swapped).expect("swapped order");
    assert_vec_close(
        swapped.relative_velocity_world_m_per_s,
        direct.relative_velocity_world_m_per_s.scale(-1.0),
    );
    assert_close(
        swapped.normal_relative_velocity_m_per_s,
        direct.normal_relative_velocity_m_per_s,
    );
}

#[test]
fn counterpart_point_gate_refuses_unrelated_base_arm_and_accepts_gap_uncertainty_boundary() {
    let mut unrelated = input(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    unrelated.base_contact_arm_body_m = Vec3::new(0.1, 0.0, 0.0);
    assert!(matches!(
        compute_patch_kinematics(unrelated),
        Err(patch_kinematics::PatchKinematicsError::CounterpartPointMismatch { .. })
    ));

    let mut bounded = input(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    bounded.profile_support.gap_m = 1.0e-2;
    bounded.profile_support.disc_point_world_m.z = 9.5e-3;
    bounded.disc_state = state(Vec3::new(0.0, 0.0, 1.0095), Vec3::ZERO, Vec3::ZERO);
    bounded.patch.gap_uncertainty_m = 5.0e-4;
    assert_eq!(
        compute_patch_kinematics(bounded)
            .expect("normal counterpart uncertainty boundary")
            .status,
        PatchContactStatus::Unknown
    );
}

#[test]
fn tangent_gauge_rotation_is_covariant_while_reconstruction_and_power_are_invariant() {
    let first = compute_patch_kinematics(input(Vec3::new(2.0, 3.0, 0.0), Vec3::ZERO, Vec3::ZERO))
        .expect("base gauge");
    let mut rotated_input = input(Vec3::new(2.0, 3.0, 0.0), Vec3::ZERO, Vec3::ZERO);
    rotated_input.tangent_gauge.rotation_rad = core::f64::consts::FRAC_PI_2;
    let rotated = compute_patch_kinematics(rotated_input).expect("rotated gauge");
    assert_vec_close(
        rotated
            .tangential_relative_velocity
            .reconstruct(rotated.tangent_basis),
        first.tangential_relative_velocity_world_m_per_s,
    );
    assert_close(
        rotated.tangential_relative_velocity.first,
        first.tangential_relative_velocity.second,
    );
    assert_close(
        rotated.tangential_relative_velocity.second,
        -first.tangential_relative_velocity.first,
    );
    assert_eq!(rotated.tangential_power_w, first.tangential_power_w);
}

#[test]
fn common_rigid_world_rotation_preserves_status_speed_and_power() {
    let original_input = input(Vec3::new(2.0, 3.0, 0.0), Vec3::ZERO, Vec3::ZERO);
    let original = compute_patch_kinematics(original_input.clone()).expect("original kinematics");
    let rotation =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), core::f64::consts::FRAC_PI_2)
            .expect("world rotation");
    let mut transformed = original_input;
    transformed.disc_state = RigidBodyState::new(
        Pose::new(
            rotation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
            rotation,
        )
        .expect("rotated disc pose"),
        rotation.rotate_body_to_world(Vec3::new(4.0, 6.0, 0.0)),
        Vec3::ZERO,
    )
    .expect("rotated disc state");
    transformed.base_state = RigidBodyState::new(
        Pose::new(Vec3::ZERO, rotation).expect("rotated base pose"),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .expect("rotated base state");
    transformed.normal_world = rotation.rotate_body_to_world(transformed.normal_world);
    transformed.profile_support.disc_arm_world_m =
        rotation.rotate_body_to_world(transformed.profile_support.disc_arm_world_m);
    transformed.profile_support.disc_point_world_m =
        rotation.rotate_body_to_world(transformed.profile_support.disc_point_world_m);
    transformed.tangent_gauge.reference_world =
        rotation.rotate_body_to_world(transformed.tangent_gauge.reference_world);
    transformed.tangent_effort_probe_world_n = transformed
        .tangent_effort_probe_world_n
        .map(|value| rotation.rotate_body_to_world(value));
    let transformed = compute_patch_kinematics(transformed).expect("rotated kinematics");
    assert_vec_close(
        transformed.relative_velocity_world_m_per_s,
        rotation.rotate_body_to_world(original.relative_velocity_world_m_per_s),
    );
    assert_close(
        transformed.reference_rolling_speed_m_per_s,
        original.reference_rolling_speed_m_per_s,
    );
    assert_eq!(transformed.status, original.status);
    assert_eq!(transformed.tangential_power_w, original.tangential_power_w);
}

#[test]
fn deterministic_degenerate_gauge_zero_speed_creepage_and_uncertainty_refuse_overclaim() {
    let mut degenerate = input(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    degenerate.tangent_gauge.reference_world = Vec3::new(0.0, 0.0, 5.0);
    let degenerate = compute_patch_kinematics(degenerate).expect("deterministic fallback");
    assert_eq!(
        degenerate.tangent_basis.source,
        TangentGaugeSource::DeterministicFallback
    );
    assert!(matches!(degenerate.creepage, Creepage::Unavailable { .. }));

    let mut uncertain = input(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    uncertain.profile_support.gap_m = 0.01;
    uncertain.profile_support.disc_point_world_m.z = 0.01;
    uncertain.disc_state = state(Vec3::new(0.0, 0.0, 1.01), Vec3::ZERO, Vec3::ZERO);
    uncertain.patch.gap_uncertainty_m = 1.0e-3;
    assert_eq!(
        compute_patch_kinematics(uncertain)
            .expect("uncertain boundary")
            .status,
        PatchContactStatus::Unknown
    );
}
