//! G0/G1/G3 public conformance for independent rolling-loss rungs.
//!
//! Every source and coefficient in this file is synthetic. These tests prove
//! algebra, refusal, replay, and work closure only.

use fs_tribo::{
    ApplicabilityRange, InputAuthority, InterfaceMedium, InterfaceSystemRef,
    partial_slip::GeneralizedWorkOwnership,
    rolling_loss::{
        CoulombContourCard, HystereticRollingCard, LEINE_STYLE_CONTOUR_LAW_ID, PatchCurvature,
        RollingKinematics, RollingLossApplicability, RollingLossChannel, RollingLossError,
        RollingLossLaw, RollingLossState, RollingLossStateKind, RollingLossUncertainty,
        RollingPatchReceipt, RollingWorkOwnership, ViscousContourCard,
    },
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1.0e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn interface() -> InterfaceSystemRef {
    InterfaceSystemRef::new(
        "fixture/roller-a->track-b",
        "fixture/rolling-history-v1",
        "fixture/interface-source-v1",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap()
}

fn patch(load_n: f64) -> RollingPatchReceipt {
    RollingPatchReceipt::new(
        "fixture/patch-a",
        "fixture/normal-patch-v1",
        "fixture/patch-source-v1",
        InputAuthority::SyntheticFixture,
        load_n,
        if load_n == 0.0 { 0.0 } else { 2.0e-4 },
        PatchCurvature::Principal {
            first_per_m: 50.0,
            second_per_m: 25.0,
        },
    )
    .unwrap()
}

fn applicability() -> RollingLossApplicability {
    RollingLossApplicability::new(
        ApplicabilityRange::new(250.0, 350.0).unwrap(),
        ApplicabilityRange::new(0.0, 100.0).unwrap(),
    )
    .unwrap()
}

fn kinematics(
    contour_speed_mps: f64,
    rolling_rate_rad_s: f64,
    spin_rate_rad_s: f64,
) -> RollingKinematics {
    RollingKinematics::new(
        contour_speed_mps,
        rolling_rate_rad_s,
        spin_rate_rad_s,
        300.0,
        10.0,
        0.5,
    )
    .unwrap()
}

fn contour_ownership() -> RollingWorkOwnership {
    RollingWorkOwnership::new(
        "fixture/patch-a",
        "fixture/rolling-interval-a",
        "fixture/contour-speed",
        RollingLossChannel::ContourDeformation,
    )
    .unwrap()
}

fn rolling_ownership() -> RollingWorkOwnership {
    RollingWorkOwnership::new(
        "fixture/patch-a",
        "fixture/rolling-interval-a",
        "fixture/rolling-rate",
        RollingLossChannel::RollingHysteresis,
    )
    .unwrap()
}

fn coulomb(coefficient: f64) -> RollingLossLaw {
    RollingLossLaw::CoulombContour(
        CoulombContourCard::new(
            LEINE_STYLE_CONTOUR_LAW_ID,
            "fixture/leine-coupon-v1",
            InputAuthority::SyntheticFixture,
            coefficient,
            applicability(),
        )
        .unwrap(),
    )
}

fn viscous(coefficient_n_s_per_m: f64) -> RollingLossLaw {
    RollingLossLaw::ViscousContour(
        ViscousContourCard::new(
            "fixture/viscous-contour-v1",
            "fixture/viscous-coupon-v1",
            InputAuthority::SyntheticFixture,
            coefficient_n_s_per_m,
            applicability(),
        )
        .unwrap(),
    )
}

fn hysteretic(length_m: f64, factor: f64) -> RollingLossLaw {
    RollingLossLaw::HystereticRollingMoment(
        HystereticRollingCard::new(
            "fixture/hysteretic-rolling-v1",
            "fixture/hysteretic-coupon-v1",
            InputAuthority::SyntheticFixture,
            length_m,
            factor,
            applicability(),
        )
        .unwrap(),
    )
}

#[test]
fn g0_zero_load_speed_and_loss_are_quiescent() {
    let zero_load = coulomb(0.2)
        .advance(
            &patch(0.0),
            &interface(),
            kinematics(4.0, 0.0, 0.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    assert_eq!(zero_load.state, RollingLossStateKind::Quiescent);
    close(zero_load.wrench.contour_force_n, 0.0);
    close(zero_load.dissipation.total_heat_j, 0.0);

    let zero_speed = viscous(7.0)
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(0.0, 0.0, 3.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    assert_eq!(zero_speed.state, RollingLossStateKind::Quiescent);
    close(zero_speed.wrench.spin_moment_nm, 0.0);

    let zero_loss = hysteretic(0.02, 0.0)
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(0.0, 4.0, 0.0),
            &rolling_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    assert_eq!(zero_loss.state, RollingLossStateKind::Quiescent);
    close(zero_loss.wrench.rolling_moment_nm, 0.0);
}

#[test]
fn g1_contour_rungs_reverse_and_preserve_declared_units() {
    let forward = coulomb(0.2)
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(3.0, 0.0, 0.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    let reverse = coulomb(0.2)
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(-3.0, 0.0, 0.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    let doubled_load = coulomb(0.2)
        .advance(
            &patch(200.0),
            &interface(),
            kinematics(3.0, 0.0, 0.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    close(forward.wrench.contour_force_n, -20.0);
    close(reverse.wrench.contour_force_n, 20.0);
    close(
        doubled_load.wrench.contour_force_n,
        2.0 * forward.wrench.contour_force_n,
    );
    close(forward.generalized_work.endpoint_body_power_w, -60.0);

    let viscous_forward = viscous(5.0)
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(3.0, 0.0, 0.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    let viscous_reverse = viscous(5.0)
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(-3.0, 0.0, 0.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    close(viscous_forward.wrench.contour_force_n, -15.0);
    close(viscous_reverse.wrench.contour_force_n, 15.0);
    close(
        viscous_forward.generalized_work.endpoint_body_power_w,
        -45.0,
    );
}

#[test]
fn g1_hysteretic_roll_spin_mixed_and_independent_work_integral_are_passive() {
    let law = hysteretic(0.01, 0.5);
    let pure_roll = law
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(0.0, 4.0, 0.0),
            &rolling_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    assert_eq!(pure_roll.state, RollingLossStateKind::RollingHysteresis);
    close(pure_roll.wrench.rolling_moment_nm, -0.5);
    close(pure_roll.generalized_work.endpoint_body_power_w, -2.0);
    close(pure_roll.generalized_work.work_into_interface_j, 1.0);
    close(
        pure_roll.dissipation.rolling_hysteresis_heat_j,
        pure_roll.generalized_work.work_into_interface_j,
    );
    let doubled_length = hysteretic(0.02, 0.5)
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(0.0, 4.0, 0.0),
            &rolling_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    close(
        doubled_length.wrench.rolling_moment_nm,
        2.0 * pure_roll.wrench.rolling_moment_nm,
    );

    let pure_spin = law
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(0.0, 0.0, 12.0),
            &rolling_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    assert_eq!(pure_spin.state, RollingLossStateKind::Quiescent);
    close(pure_spin.generalized_work.work_into_interface_j, 0.0);

    let mixed = law
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(7.0, -4.0, 12.0),
            &rolling_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    close(mixed.wrench.rolling_moment_nm, 0.5);
    close(mixed.generalized_work.endpoint_body_power_w, -2.0);
    close(mixed.generalized_work.work_into_interface_j, 1.0);
    assert!(mixed.generalized_work.endpoint_body_power_w <= 0.0);
    assert_eq!(
        mixed.uncertainty,
        RollingLossUncertainty::ModelFormNoNumericalCertificate {
            curvature_approximation_authority: None,
        }
    );
}

#[test]
fn g0_temperature_frequency_and_ownership_refuse() {
    let law = coulomb(0.2);
    let too_hot = RollingKinematics::new(1.0, 0.0, 0.0, 351.0, 10.0, 0.1).unwrap();
    assert!(matches!(
        law.advance(
            &patch(100.0),
            &interface(),
            too_hot,
            &contour_ownership(),
            &RollingLossState::zero(),
        ),
        Err(RollingLossError::OutsideApplicability {
            field: "temperature_kelvin",
            ..
        })
    ));
    let too_fast = RollingKinematics::new(1.0, 0.0, 0.0, 300.0, 101.0, 0.1).unwrap();
    assert!(matches!(
        law.advance(
            &patch(100.0),
            &interface(),
            too_fast,
            &contour_ownership(),
            &RollingLossState::zero(),
        ),
        Err(RollingLossError::OutsideApplicability {
            field: "excitation_frequency_hz",
            ..
        })
    ));
    let partial = GeneralizedWorkOwnership::new(
        "fixture/patch-a",
        "fixture/rolling-interval-a",
        "fixture/longitudinal",
        "fixture/lateral",
        "fixture/torsion",
    )
    .unwrap();
    assert_eq!(
        contour_ownership().require_disjoint_from_partial_slip(&partial),
        Err(RollingLossError::WorkOwnershipOverlap)
    );
    assert!(matches!(
        hysteretic(0.01, 0.5).advance(
            &patch(100.0),
            &interface(),
            kinematics(0.0, 1.0, 0.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        ),
        Err(RollingLossError::WorkOwnershipMismatch { field: "channel" })
    ));
}

#[test]
fn g3_checkpoint_replay_binds_source_cards_and_patch_receipts() {
    let law = coulomb(0.2);
    let step = law
        .advance(
            &patch(100.0),
            &interface(),
            kinematics(2.0, 0.0, 0.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    assert_eq!(
        law.restore_checkpoint(
            &patch(100.0),
            &interface(),
            &contour_ownership(),
            &step.checkpoint,
        )
        .unwrap(),
        step.next_state
    );

    let mutated_source = RollingLossLaw::CoulombContour(
        CoulombContourCard::new(
            LEINE_STYLE_CONTOUR_LAW_ID,
            "fixture/changed-coupon-source",
            InputAuthority::SyntheticFixture,
            0.2,
            applicability(),
        )
        .unwrap(),
    );
    assert_eq!(
        mutated_source.restore_checkpoint(
            &patch(100.0),
            &interface(),
            &contour_ownership(),
            &step.checkpoint,
        ),
        Err(RollingLossError::CheckpointMismatch { field: "law" })
    );

    let changed_patch = RollingPatchReceipt::new(
        "fixture/patch-a",
        "fixture/normal-patch-v1",
        "fixture/changed-patch-source",
        InputAuthority::SyntheticFixture,
        100.0,
        2.0e-4,
        PatchCurvature::EquivalentRadiusApproximation {
            radius_m: 0.02,
            authority: InputAuthority::Estimated,
        },
    )
    .unwrap();
    assert_eq!(
        law.restore_checkpoint(
            &changed_patch,
            &interface(),
            &contour_ownership(),
            &step.checkpoint,
        ),
        Err(RollingLossError::CheckpointMismatch { field: "patch" })
    );
    let approximation_step = law
        .advance(
            &changed_patch,
            &interface(),
            kinematics(2.0, 0.0, 0.0),
            &contour_ownership(),
            &RollingLossState::zero(),
        )
        .unwrap();
    assert_eq!(
        approximation_step.uncertainty,
        RollingLossUncertainty::ModelFormNoNumericalCertificate {
            curvature_approximation_authority: Some(InputAuthority::Estimated),
        }
    );

    let changed_interface = InterfaceSystemRef::new(
        "fixture/roller-a->track-b",
        "fixture/changed-rolling-history",
        "fixture/interface-source-v1",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap();
    assert_eq!(
        law.restore_checkpoint(
            &patch(100.0),
            &changed_interface,
            &contour_ownership(),
            &step.checkpoint,
        ),
        Err(RollingLossError::CheckpointMismatch { field: "interface" })
    );
}

#[test]
fn g0_source_has_no_target_shaping_fields() {
    let source = include_str!("../src/rolling_loss.rs");
    for forbidden in ["Euler", "Mould", "one_millimetre", "stop_time", "ring_cone"] {
        assert!(!source.contains(forbidden), "forbidden field: {forbidden}");
    }
}
