//! G0/G3 coverage for the thin rolling/contour world-wrench adapter.

#[path = "../src/rolling_contact.rs"]
mod rolling_contact;

use fs_mbd::Vec3;
use fs_tribo::{
    ApplicabilityRange, InputAuthority, InterfaceMedium, InterfaceSystemRef,
    partial_slip::GeneralizedWorkOwnership,
    rolling_loss::{
        CoulombContourCard, HystereticRollingCard, LEINE_STYLE_CONTOUR_LAW_ID, PatchCurvature,
        RollingLossApplicability, RollingLossChannel, RollingLossLaw, RollingLossState,
        RollingPatchReceipt, RollingWorkOwnership, ViscousContourCard,
    },
};
use rolling_contact::{
    ROLLING_CONTACT_ADAPTER_ID, RollingContactError, RollingContactIdentity, RollingContactInput,
    SpinMicroslipAvailability, evaluate_rolling_contact,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1.0e-12 * actual.abs().max(expected.abs()).max(1.0),
        "{actual} != {expected}"
    );
}

fn close_vec(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}

fn applicability() -> RollingLossApplicability {
    RollingLossApplicability::new(
        ApplicabilityRange::new(250.0, 350.0).unwrap(),
        ApplicabilityRange::new(0.0, 100.0).unwrap(),
    )
    .unwrap()
}

fn patch() -> RollingPatchReceipt {
    RollingPatchReceipt::new(
        "patch/rim-v1",
        "normal/hertz-caller-v1",
        "source/patch-v1",
        InputAuthority::SyntheticFixture,
        10.0,
        2.0e-4,
        PatchCurvature::Principal {
            first_per_m: 40.0,
            second_per_m: 20.0,
        },
    )
    .unwrap()
}

fn interface(history: &str) -> InterfaceSystemRef {
    InterfaceSystemRef::new(
        "disc/rim->base/track",
        history,
        "source/interface-v1",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap()
}

fn contour_ownership() -> RollingWorkOwnership {
    RollingWorkOwnership::new(
        "patch/rim-v1",
        "interval/rolling-v1",
        "coordinate/contour-speed",
        RollingLossChannel::ContourDeformation,
    )
    .unwrap()
}

fn rolling_ownership() -> RollingWorkOwnership {
    RollingWorkOwnership::new(
        "patch/rim-v1",
        "interval/rolling-v1",
        "coordinate/rolling-rate",
        RollingLossChannel::RollingHysteresis,
    )
    .unwrap()
}

fn coulomb(coefficient: f64) -> RollingLossLaw {
    RollingLossLaw::CoulombContour(
        CoulombContourCard::new(
            LEINE_STYLE_CONTOUR_LAW_ID,
            "coupon/leine-v1",
            InputAuthority::SyntheticFixture,
            coefficient,
            applicability(),
        )
        .unwrap(),
    )
}

fn viscous(coefficient: f64) -> RollingLossLaw {
    RollingLossLaw::ViscousContour(
        ViscousContourCard::new(
            "law/viscous-contour-v1",
            "coupon/viscous-v1",
            InputAuthority::SyntheticFixture,
            coefficient,
            applicability(),
        )
        .unwrap(),
    )
}

fn hysteretic(loss_length_m: f64) -> RollingLossLaw {
    RollingLossLaw::HystereticRollingMoment(
        HystereticRollingCard::new(
            "law/hysteretic-v1",
            "coupon/hysteretic-v1",
            InputAuthority::SyntheticFixture,
            loss_length_m,
            0.5,
            applicability(),
        )
        .unwrap(),
    )
}

fn input(law: RollingLossLaw, ownership: RollingWorkOwnership) -> RollingContactInput {
    RollingContactInput {
        identity: RollingContactIdentity {
            case_id: "case/adapter-v1".into(),
            adapter_id: ROLLING_CONTACT_ADAPTER_ID.into(),
            world_frame_id: "frame/world-v1".into(),
            port_id: "port/rolling-v1".into(),
            domain_id: "domain/rim-interval-v1".into(),
        },
        patch: patch(),
        interface: interface("history/rolling-v1"),
        law,
        state: RollingLossState::zero(),
        checkpoint: None,
        ownership,
        partial_slip_ownership: None,
        contact_arm_world_m: Vec3::new(0.0, 0.0, -0.02),
        contour_tangent_axis_world: Vec3::new(1.0, 0.0, 0.0),
        rolling_axis_world: Vec3::new(0.0, 1.0, 0.0),
        contour_speed_mps: 2.0,
        rolling_rate_rad_s: 3.0,
        spin_rate_rad_s: 0.25,
        temperature_kelvin: 300.0,
        excitation_frequency_hz: 10.0,
        interval_s: 0.5,
    }
}

#[test]
fn g0_zero_reversal_and_scale_keep_scalar_law_signs() {
    let mut zero = input(coulomb(0.1), contour_ownership());
    zero.contour_speed_mps = 0.0;
    let quiescent = evaluate_rolling_contact(&zero).unwrap();
    close_vec(quiescent.body_wrench.contour_force_world_n, Vec3::ZERO);
    close(quiescent.generic.dissipation.total_heat_j, 0.0);

    let forward = evaluate_rolling_contact(&input(coulomb(0.1), contour_ownership())).unwrap();
    let mut reverse_input = input(coulomb(0.1), contour_ownership());
    reverse_input.contour_speed_mps = -2.0;
    let reverse = evaluate_rolling_contact(&reverse_input).unwrap();
    close(
        reverse.body_wrench.contour_force_world_n.x,
        -forward.body_wrench.contour_force_world_n.x,
    );
    close(
        forward.generic.dissipation.total_heat_j,
        reverse.generic.dissipation.total_heat_j,
    );

    let scaled = evaluate_rolling_contact(&input(coulomb(0.2), contour_ownership())).unwrap();
    close(
        scaled.body_wrench.contour_force_world_n.x,
        2.0 * forward.body_wrench.contour_force_world_n.x,
    );
}

#[test]
fn g0_world_force_recenter_witness_and_rolling_free_couple_hold() {
    let contour = evaluate_rolling_contact(&input(coulomb(0.1), contour_ownership())).unwrap();
    close_vec(
        contour.body_wrench.contour_force_world_n,
        Vec3::new(-1.0, 0.0, 0.0),
    );
    close_vec(
        contour.body_wrench.contour_force_moment_about_com_world_nm,
        Vec3::new(0.0, 0.02, 0.0),
    );

    let rolling = evaluate_rolling_contact(&input(hysteretic(0.02), rolling_ownership())).unwrap();
    close_vec(rolling.body_wrench.contour_force_world_n, Vec3::ZERO);
    close_vec(
        rolling.body_wrench.rolling_free_couple_world_nm,
        Vec3::new(0.0, -0.1, 0.0),
    );
    close_vec(
        rolling.body_wrench.total_moment_about_com_world_nm,
        Vec3::new(0.0, -0.1, 0.0),
    );
}

#[test]
fn g3_rigid_rotation_preserves_power_and_wrench_objectivity() {
    let original = evaluate_rolling_contact(&input(coulomb(0.1), contour_ownership())).unwrap();
    let mut rotated_input = input(coulomb(0.1), contour_ownership());
    rotated_input.contour_tangent_axis_world = Vec3::new(0.0, 1.0, 0.0);
    rotated_input.rolling_axis_world = Vec3::new(-1.0, 0.0, 0.0);
    let rotated = evaluate_rolling_contact(&rotated_input).unwrap();
    close_vec(
        rotated.body_wrench.contour_force_world_n,
        Vec3::new(0.0, -1.0, 0.0),
    );
    close_vec(
        rotated.body_wrench.contour_force_moment_about_com_world_nm,
        Vec3::new(-0.02, 0.0, 0.0),
    );
    close(
        rotated.body_wrench.body_power_w,
        original.body_wrench.body_power_w,
    );
    close(
        rotated.body_wrench.body_power_w,
        rotated.generic.generalized_work.endpoint_body_power_w,
    );
}

#[test]
fn g0_all_rival_ids_are_retained_without_blending() {
    let contour = evaluate_rolling_contact(&input(coulomb(0.1), contour_ownership())).unwrap();
    let viscous = evaluate_rolling_contact(&input(viscous(0.3), contour_ownership())).unwrap();
    let hysteretic =
        evaluate_rolling_contact(&input(hysteretic(0.02), rolling_ownership())).unwrap();
    assert_eq!(
        contour.generic.applicability,
        fs_tribo::rolling_loss::RollingLossApplicabilityKind::LeineStyleCoulombContour
    );
    assert_eq!(
        viscous.generic.applicability,
        fs_tribo::rolling_loss::RollingLossApplicabilityKind::ViscousContour
    );
    assert_eq!(
        hysteretic.generic.applicability,
        fs_tribo::rolling_loss::RollingLossApplicabilityKind::FinitePatchHystereticRollingMoment
    );
    assert_ne!(
        contour.generic.wrench.contour_force_n,
        viscous.generic.wrench.contour_force_n
    );
    assert_eq!(contour.law_model_id, LEINE_STYLE_CONTOUR_LAW_ID);
    assert_eq!(viscous.law_model_id, "law/viscous-contour-v1");
    assert_eq!(hysteretic.law_model_id, "law/hysteretic-v1");
    assert_ne!(contour.law_model_id, hysteretic.law_model_id);
}

#[test]
fn g0_checkpoint_mutation_and_work_overlap_refuse() {
    let first_input = input(coulomb(0.1), contour_ownership());
    let first = evaluate_rolling_contact(&first_input).unwrap();
    let mut replay = first_input.clone();
    replay.state = first.generic.next_state.clone();
    replay.checkpoint = Some(first.generic.checkpoint.clone());
    evaluate_rolling_contact(&replay).unwrap();
    replay.interface = interface("history/mutated-v1");
    assert!(matches!(
        evaluate_rolling_contact(&replay),
        Err(RollingContactError::GenericRefusal(_))
    ));

    let mut overlap = input(coulomb(0.1), contour_ownership());
    overlap.partial_slip_ownership = Some(
        GeneralizedWorkOwnership::new(
            "patch/rim-v1",
            "interval/rolling-v1",
            "coordinate/longitudinal",
            "coordinate/lateral",
            "coordinate/torsion",
        )
        .unwrap(),
    );
    assert!(matches!(
        evaluate_rolling_contact(&overlap),
        Err(RollingContactError::GenericRefusal(_))
    ));
}

#[test]
fn g0_spin_only_is_typed_no_mechanism_and_bad_axes_refuse() {
    let mut spin = input(coulomb(0.1), contour_ownership());
    spin.contour_speed_mps = 0.0;
    spin.rolling_rate_rad_s = 0.0;
    spin.spin_rate_rad_s = 9.0;
    let result = evaluate_rolling_contact(&spin).unwrap();
    close_vec(result.body_wrench.contour_force_world_n, Vec3::ZERO);
    close_vec(result.body_wrench.rolling_free_couple_world_nm, Vec3::ZERO);
    assert_eq!(
        result.spin_microslip,
        SpinMicroslipAvailability::Unavailable
    );

    spin.contour_tangent_axis_world = Vec3::new(2.0, 0.0, 0.0);
    assert_eq!(
        evaluate_rolling_contact(&spin),
        Err(RollingContactError::NonOrthonormalAxes)
    );
}

#[test]
fn g0_source_has_no_target_fields() {
    let source = include_str!("../src/rolling_contact.rs").to_ascii_lowercase();
    for forbidden in ["mould", "video", "stop_time", "edge_optimum"] {
        assert!(
            !source.contains(forbidden),
            "forbidden target field: {forbidden}"
        );
    }
}
