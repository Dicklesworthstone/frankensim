//! G0/G3 fixtures for the Euler finite-patch tangential adapter.

#[path = "../src/contact_dynamics.rs"]
mod contact_dynamics;
#[path = "../src/patch_kinematics.rs"]
mod patch_kinematics;
#[path = "../src/tangential_contact.rs"]
mod tangential_contact;

use fs_contact::tangential::smooth::{
    SmoothAuthorityPolicy, SmoothRegularization, SmoothTangentialAdapter,
};
use fs_couple::StableId;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_rep_frep::AxisymmetricSupportAuthority;
use fs_tribo::partial_slip::{
    GeneralizedWorkOwnership, NormalPatchAuthority, NormalPatchView, PARTIAL_SLIP_MODEL_ID,
    PartialSlipInterface, PartialSlipLaw, PartialSlipParameters,
};
use patch_kinematics::{
    OrderedSurfacePair, PatchContactStatus, PatchGeometryMetadata, PatchKinematicThresholds,
    PatchKinematicsInput, ProfileSupportKinematics, SurfaceOrder, TangentGaugeInput,
    compute_patch_kinematics,
};
use tangential_contact::{
    EulerTangentialContactAdapter, TangentialContactError, TangentialContactLane,
    TangentialContactRequest,
};

fn id(value: &str) -> StableId {
    StableId::new(value).expect("test identity")
}

fn properties() -> MassProperties {
    MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).expect("test inertia")
}

fn state(velocity: Vec3, angular_velocity: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).expect("test pose"),
        velocity.scale(properties().mass()),
        angular_velocity,
    )
    .expect("test state")
}

fn thresholds() -> PatchKinematicThresholds {
    PatchKinematicThresholds {
        threshold_identity: id("tangent-thresholds"),
        tie_break_identity: id("tangent-ties"),
        separation_gap_m: 1.0e-3,
        touching_gap_m: 1.0e-6,
        support_point_coincidence_tolerance_m: 1.0e-12,
        tangent_counterpart_coincidence_tolerance_m: 1.0e-12,
        approach_speed_m_per_s: 1.0e-3,
        impact_candidate_speed_m_per_s: 1.0,
        stationary_tangent_speed_m_per_s: 1.0e-12,
        minimum_reference_rolling_speed_m_per_s: 1.0e-6,
        gauge_degeneracy_norm: 1.0e-12,
    }
}

fn patch_kinematics(
    order: SurfaceOrder,
    creepage: f64,
    spin_rad_per_s: f64,
    gauge_rotation_rad: f64,
) -> patch_kinematics::PatchKinematics {
    let mean = 10.0;
    let relative = mean * creepage;
    let (disc_velocity, base_velocity) = match order {
        SurfaceOrder::DiscThenBase => (mean + relative / 2.0, mean - relative / 2.0),
        SurfaceOrder::BaseThenDisc => (mean - relative / 2.0, mean + relative / 2.0),
    };
    compute_patch_kinematics(PatchKinematicsInput {
        surfaces: OrderedSurfacePair::try_new(id("disc"), id("base"), order).expect("ordered pair"),
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        profile_support: ProfileSupportKinematics {
            disc_arm_world_m: Vec3::ZERO,
            disc_point_world_m: Vec3::ZERO,
            gap_m: 0.0,
            source_feature: 0,
            support_authority: AxisymmetricSupportAuthority::Estimate,
        },
        patch: PatchGeometryMetadata {
            patch_identity: id("patch-a"),
            source_feature: 0,
            gap_uncertainty_m: 0.0,
            curvature: patch_kinematics::CurvatureMetadata::Unavailable {
                curvature_identity: id("curvature-unavailable"),
                reason_identity: id("curvature-not-needed"),
            },
        },
        disc_state: state(
            Vec3::new(disc_velocity, 0.0, 0.0),
            Vec3::new(0.0, 0.0, spin_rad_per_s),
        ),
        disc_mass_properties: properties(),
        base_state: state(Vec3::new(base_velocity, 0.0, 0.0), Vec3::ZERO),
        base_mass_properties: properties(),
        base_contact_arm_body_m: Vec3::ZERO,
        tangent_gauge: TangentGaugeInput {
            reference_world: Vec3::new(1.0, 0.0, 0.0),
            rotation_rad: gauge_rotation_rad,
        },
        thresholds: thresholds(),
        tangent_effort_probe_world_n: None,
    })
    .expect("admitted patch kinematics")
}

fn normal_patch() -> NormalPatchView {
    NormalPatchView::new(
        "patch-a",
        "fixture/finite-normal-patch",
        "fixture/normal-source",
        NormalPatchAuthority::SyntheticFixture,
        100.0,
        0.02,
        0.01,
        1.0e-4,
    )
    .expect("finite patch")
}

fn interface() -> PartialSlipInterface {
    PartialSlipInterface::new(
        "disc->base",
        "fixture/history-a",
        "fixture/interface-source",
        NormalPatchAuthority::SyntheticFixture,
    )
    .expect("interface")
}

fn law() -> PartialSlipLaw {
    PartialSlipLaw::new(
        PARTIAL_SLIP_MODEL_ID,
        "fixture/friction-card",
        PartialSlipParameters {
            static_mu: 0.8,
            kinetic_mu: 0.4,
            tangential_stiffness_n_per_m: 10_000.0,
            torsional_stiffness_nm_per_rad: 100.0,
            torsional_capacity_factor: 0.5,
            partial_slip_onset_fraction: 0.5,
            partial_slip_hardening_fraction: 0.4,
        },
    )
    .expect("explicit synthetic law")
}

fn direct_adapter() -> EulerTangentialContactAdapter {
    EulerTangentialContactAdapter::new(
        "euler-direct-adapter",
        "fixture/euler-direct",
        TangentialContactLane::PartialSlip { law: law() },
    )
    .expect("direct adapter")
}

fn smooth_adapter() -> EulerTangentialContactAdapter {
    EulerTangentialContactAdapter::new(
        "euler-smooth-adapter",
        "fixture/euler-smooth",
        TangentialContactLane::Smooth {
            law: law(),
            adapter: SmoothTangentialAdapter::new(
                "generic-smooth-adapter",
                "fixture/generic-smooth",
                SmoothRegularization {
                    creepage_scale: 1.0e-8,
                    torsional_spin_scale_rad_per_s: 1.0e-8,
                    tangent_probe_creepage: 1.0e-5,
                    tangent_probe_spin_rad_per_s: 1.0e-5,
                },
                SmoothAuthorityPolicy::test_only(),
            )
            .expect("smooth adapter"),
        },
    )
    .expect("smooth Euler adapter")
}

fn request(
    interval: &str,
    version: u64,
    order: SurfaceOrder,
    creepage: f64,
    spin_rad_per_s: f64,
    gauge_rotation_rad: f64,
) -> TangentialContactRequest {
    TangentialContactRequest {
        request_id: format!("request-{interval}"),
        expected_state_version: version,
        patch_kinematics: patch_kinematics(order, creepage, spin_rad_per_s, gauge_rotation_rad),
        normal_patch: normal_patch(),
        interface: interface(),
        work_ownership: GeneralizedWorkOwnership::new(
            "patch-a",
            interval,
            "creep-long",
            "creep-lat",
            "spin",
        )
        .expect("work ownership"),
        dt_s: 1.0e-3,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1.0e-11 * actual.abs().max(expected.abs()).max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn assert_vec_close(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}

#[test]
fn stick_partial_slip_and_gross_slide_keep_mode_work_and_loss_explicit() {
    let adapter = direct_adapter();
    let expectations = [
        (
            "stick",
            0.2,
            fs_tribo::partial_slip::PartialSlipStateKind::Sticking,
        ),
        (
            "partial",
            0.6,
            fs_tribo::partial_slip::PartialSlipStateKind::PartialSlip,
        ),
        (
            "gross",
            1.2,
            fs_tribo::partial_slip::PartialSlipStateKind::GrossSlide,
        ),
    ];
    for (index, (interval, creepage, mode)) in expectations.into_iter().enumerate() {
        let state = adapter
            .initial_state(&normal_patch(), &interface(), 8)
            .expect("independent zero-history state");
        let receipt = adapter
            .prepare(
                &state,
                &request(
                    interval,
                    index as u64,
                    SurfaceOrder::DiscThenBase,
                    creepage,
                    0.0,
                    0.0,
                ),
            )
            .expect("admitted direct lane");
        assert_eq!(receipt.mode, mode);
        assert!(receipt.irreversible_loss_j >= 0.0);
        assert!(receipt.heat_j >= 0.0);
        close(
            receipt.exact_relative_body_work_j,
            -receipt.signed_storage_change_j - receipt.irreversible_loss_j,
        );
        assert_eq!(
            adapter
                .commit(&state, &receipt)
                .expect("exactly once commit")
                .committed_version(),
            1
        );
    }
}

#[test]
fn pure_spin_has_free_torsion_and_nonnegative_loss_without_creepage() {
    let adapter = direct_adapter();
    let state = adapter
        .initial_state(&normal_patch(), &interface(), 4)
        .expect("initial state");
    let receipt = adapter
        .prepare(
            &state,
            &request("pure-spin", 0, SurfaceOrder::DiscThenBase, 0.0, 5.0, 0.0),
        )
        .expect("pure spin is a valid finite-patch input");
    assert_eq!(
        receipt.mode,
        fs_tribo::partial_slip::PartialSlipStateKind::GrossSlide
    );
    assert_eq!(receipt.force_on_disc_world_n, Vec3::ZERO);
    assert!(receipt.free_torsional_torque_on_disc_world_nm.z.abs() > 0.0);
    assert!(receipt.irreversible_loss_j >= 0.0);
    assert!(receipt.heat_j >= 0.0);
}

#[test]
fn tangent_basis_rotation_and_surface_reversal_preserve_disc_wrench() {
    let adapter = direct_adapter();
    let state = adapter
        .initial_state(&normal_patch(), &interface(), 8)
        .expect("initial state");
    let original = adapter
        .prepare(
            &state,
            &request("basis-a", 0, SurfaceOrder::DiscThenBase, 0.6, 0.3, 0.0),
        )
        .expect("base gauge");
    let rotated = adapter
        .prepare(
            &state,
            &request(
                "basis-b",
                0,
                SurfaceOrder::DiscThenBase,
                0.6,
                0.3,
                core::f64::consts::FRAC_PI_2,
            ),
        )
        .expect("rotated gauge");
    assert_vec_close(
        rotated.force_on_disc_world_n,
        original.force_on_disc_world_n,
    );
    assert_vec_close(
        rotated.free_torsional_torque_on_disc_world_nm,
        original.free_torsional_torque_on_disc_world_nm,
    );
    let reversed = adapter
        .prepare(
            &state,
            &request("reversed", 0, SurfaceOrder::BaseThenDisc, 0.6, 0.3, 0.0),
        )
        .expect("reversed ordered pair");
    assert_vec_close(
        reversed.force_on_disc_world_n,
        original.force_on_disc_world_n,
    );
    assert_vec_close(
        reversed.free_torsional_torque_on_disc_world_nm,
        original.free_torsional_torque_on_disc_world_nm,
    );
}

#[test]
fn smooth_lane_preserves_checkpoint_recontact_and_exact_work_identity() {
    let adapter = smooth_adapter();
    let initial = adapter
        .initial_state(&normal_patch(), &interface(), 4)
        .expect("smooth initial state");
    let first = adapter
        .prepare(
            &initial,
            &request("smooth-first", 0, SurfaceOrder::DiscThenBase, 0.6, 0.2, 0.0),
        )
        .expect("smooth receipt");
    let committed = adapter.commit(&initial, &first).expect("smooth commit");
    let restored = adapter
        .restore_state(&normal_patch(), &interface(), first.checkpoint.clone())
        .expect("receipt checkpoint restores");
    assert_eq!(restored, committed);
    let next_request = request(
        "smooth-second",
        1,
        SurfaceOrder::DiscThenBase,
        0.2,
        0.0,
        0.0,
    );
    let replay_a = adapter
        .prepare(&restored, &next_request)
        .expect("recontact a");
    let replay_b = adapter
        .prepare(&committed, &next_request)
        .expect("recontact b");
    assert_eq!(replay_a, replay_b);
    assert!(replay_a.irreversible_loss_j >= 0.0);
    close(
        replay_a.exact_relative_body_work_j,
        -replay_a.signed_storage_change_j - replay_a.irreversible_loss_j,
    );
}

#[test]
fn missing_friction_unsupported_state_and_wrong_work_owner_refuse() {
    assert!(
        PartialSlipLaw::new(
            PARTIAL_SLIP_MODEL_ID,
            "fixture/missing-friction",
            PartialSlipParameters {
                static_mu: 0.0,
                ..law().parameters()
            },
        )
        .is_err()
    );
    let adapter = direct_adapter();
    let state = adapter
        .initial_state(&normal_patch(), &interface(), 4)
        .expect("initial state");
    let mut unsupported = request("separated", 0, SurfaceOrder::DiscThenBase, 0.2, 0.0, 0.0);
    unsupported.patch_kinematics.status = PatchContactStatus::Separated;
    assert!(matches!(
        adapter.prepare(&state, &unsupported),
        Err(TangentialContactError::UnsupportedPatchStatus {
            status: PatchContactStatus::Separated
        })
    ));
    let mut wrong_owner = request("wrong-owner", 0, SurfaceOrder::DiscThenBase, 0.2, 0.0, 0.0);
    wrong_owner.work_ownership =
        GeneralizedWorkOwnership::new("other-patch", "wrong-owner", "x", "y", "spin")
            .expect("syntactically valid wrong owner");
    assert!(matches!(
        adapter.prepare(&state, &wrong_owner),
        Err(TangentialContactError::PartialSlip(
            fs_tribo::partial_slip::PartialSlipError::WorkOwnershipMismatch
        ))
    ));
}
