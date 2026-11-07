//! Physical path references and transactional adaptive depth-wear checks.

use super::*;
use fs_contact::normal_patch::plastic_source::depth_wear::integration::{
    DepthHardnessIntegratedWearAdvance, DepthHardnessIntegrationError, WearIntegrationControls,
};
use fs_contact::normal_patch::{FiniteGapPlasticDrive, FiniteGapPlasticLoadControls};

fn accuracy(scale: f64) -> WearIntegrationControls {
    WearIntegrationControls {
        absolute_depth_tolerance_m: 1.0e-10 * scale,
        absolute_pressure_tolerance_pa: 10.0 * scale,
        relative_tolerance: 1.0e-5 * scale,
        maximum_accepted_substeps: 16384,
        maximum_step_evaluations: 60000,
    }
}
fn load_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}
fn compliance(source: &SourceDepthHardnessWear) -> f64 {
    let config = source.model().config();
    4.0 * config.grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (core::f64::consts::PI * config.reduced_modulus_pa)
}
fn slip(distance: f64) -> Vec<f64> {
    let mut values = vec![0.0; 25];
    values[12] = distance;
    values
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0e-30),
        "{actual:e} != {expected:e}"
    );
}
fn accepted_accounts(
    before: &SourceDepthHardnessWear,
    advance: &DepthHardnessIntegratedWearAdvance,
) {
    let next = &advance.candidate;
    let receipt = &advance.receipt.path;
    assert_eq!(advance.source_identity, before.source_identity());
    assert!(receipt.maximum_accepted_normalized_error <= 1.0);
    for side in 0..2 {
        close(
            receipt.removed_volume_m3[side],
            next.volumes_m3()[side] - before.volumes_m3()[side],
        );
    }
    close(
        receipt.plastic_work_j,
        next.contact_state().cumulative_plastic_work_j()
            - before.contact_state().cumulative_plastic_work_j(),
    );
    let energy_scale = receipt.initial_elastic_energy_j.abs()
        + receipt.final_elastic_energy_j.abs()
        + receipt.plastic_work_j.abs()
        + receipt.normal_work_j.abs()
        + receipt.parameter_work_j.abs();
    assert!(
        (receipt.elastic_energy_change_j
            - (receipt.final_elastic_energy_j - receipt.initial_elastic_energy_j))
            .abs()
            <= 128.0 * f64::EPSILON * energy_scale
    );
    assert!(receipt.closure_residual_j.abs() <= 128.0 * f64::EPSILON * energy_scale);
}

#[test]
fn g1_adaptive_graded_cap_feedback_converges_to_independent_continuous_wear() {
    with_cx(|cx, _| {
        let before = Fixture::new(GRADIENT)
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
            .advance_normal(5.0e-5, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let distance = 1.0e-4;
        // The cap remains active: p=H0+g*z_A, dz_A/ds=k_A and
        // dz_B/ds=k_B*p/H_B. These are continuous-path references, not the
        // frozen-pressure depth integral evaluated by a numerical substep.
        let exact = [
            K[0] * distance,
            K[1] / HB * (H0 * distance + 0.5 * GRADIENT * K[0] * distance * distance),
        ];
        let loose = before
            .integrate_wear(&slip(distance), accuracy(10.0), cx)
            .unwrap();
        let tight = before
            .integrate_wear(&slip(distance), accuracy(0.01), cx)
            .unwrap();
        let loose_error = (loose.candidate.recession_m()[0][12] - exact[0]).abs();
        let tight_error = (tight.candidate.recession_m()[0][12] - exact[0]).abs();
        assert!(tight_error < 0.25 * loose_error);
        for side in 0..2 {
            near(tight.candidate.recession_m()[side][12], exact[side], 0.002);
            near(
                tight.candidate.volumes_m3()[side],
                AREA * exact[side],
                0.002,
            );
        }
        near(
            tight.candidate.contact_state().pressure_pa()[12],
            H0 + GRADIENT * exact[0],
            0.002,
        );
        let receipt = &tight.receipt.path;
        near(
            receipt.pressure_weighted_sliding_j,
            AREA * (H0 * distance + 0.5 * GRADIENT * K[0] * distance * distance),
            0.002,
        );
        // At held approach d_plastic = approach - z_A - z_B - C*p.
        // Integrating p*d(d_plastic) independently also checks that rejected
        // coarse-trial plastic work never entered the returned accounts.
        let c = compliance(&before);
        let a = -K[0] - K[1] * H0 / HB - c * GRADIENT * K[0];
        let b = -K[1] * GRADIENT * K[0] / HB;
        let plastic = AREA
            * (H0 * a * distance
                + 0.5 * (H0 * b + GRADIENT * K[0] * a) * distance * distance
                + (GRADIENT * K[0] * b / 3.0) * distance.powi(3));
        near(receipt.plastic_work_j, plastic, 0.003);
        assert!(receipt.rejected_intervals > 0);
        assert!(receipt.step_evaluations > receipt.accepted_substeps);
        assert!(tight.receipt.maximum_accepted_local_hardness_error_pa > 0.0);
        assert_eq!(receipt.normal_work_j, 0.0);
        accepted_accounts(&before, &tight);
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_adaptive_load_retains_exact_depth_integration_and_zero_path_boundary() {
    with_cx(|cx, _| {
        let force = 0.005;
        let before = Fixture::new(GRADIENT)
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
            .advance_normal_at_load(force, load_controls(), cx)
            .unwrap()
            .normal
            .candidate;
        let distance = 1.0e-4;
        let direct = before
            .advance_wear_integrated_at_load(&slip(distance), force, load_controls(), cx)
            .unwrap();
        let integrated = before
            .integrate_wear_at_load(&slip(distance), force, load_controls(), accuracy(0.01), cx)
            .unwrap();
        for side in 0..2 {
            close(
                integrated.candidate.recession_m()[side][12],
                direct.wear.candidate.recession_m()[side][12],
            );
        }
        let receipt = &integrated.receipt.path;
        assert_eq!(receipt.accepted_substeps, 2);
        assert_eq!(receipt.step_evaluations, 3);
        close(receipt.pressure_weighted_sliding_j, force * distance);
        close(
            receipt.normal_work_j,
            force
                * (integrated.candidate.contact_state().approach_m()
                    - before.contact_state().approach_m()),
        );
        accepted_accounts(&before, &integrated);
        let zero = before
            .integrate_wear_at_load(&[0.0; 25], force, load_controls(), accuracy(1.0), cx)
            .unwrap();
        assert_eq!(zero.candidate, before);
        assert_eq!(zero.receipt.path.accepted_substeps, 1);
        assert_eq!(zero.receipt.path.step_evaluations, 1);
        assert_eq!(zero.receipt.path.removed_volume_m3, [0.0; 2]);
        assert!(matches!(
            before.integrate_wear_at_load(&[0.0; 25], 0.004, load_controls(), accuracy(1.0), cx),
            Err(DepthHardnessIntegrationError::Step(
                DepthHardnessWearError::Load(_)
            ))
        ));
    });
}

#[test]
fn g1_adaptive_spring_recovers_supported_path_after_coarse_profile_overshoot() {
    with_cx(|cx, _| {
        let pristine =
            Fixture::new(0.0).owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx);
        let c = compliance(&pristine);
        let initial_pressure = 0.1e6;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: AREA / c,
            position_m: 2.0 * c * initial_pressure,
        };
        let loaded = pristine
            .advance_normal_with_drive(drive, load_controls(), cx)
            .unwrap();
        let initial_spring_energy = loaded.drive.spring_energy_j;
        let before = loaded.normal.candidate;
        let saved = before.clone();
        let distance = 0.004;
        assert!(matches!(
            before.advance_wear_integrated_with_drive(&slip(distance), drive, load_controls(), cx),
            Err(DepthHardnessWearError::Hardness(_))
        ));
        let integrated = before
            .integrate_wear_with_drive(&slip(distance), drive, load_controls(), accuracy(0.01), cx)
            .unwrap();
        // Original constant body profiles give dp/ds = -gamma*p/C_eff.
        let effective = c + AREA / drive.stiffness_n_per_m;
        let gamma = K[0] / H0 + K[1] / HB;
        let pressure = initial_pressure * (-gamma * distance / effective).exp();
        let removed = effective * (initial_pressure - pressure);
        for (side, h) in [H0, HB].into_iter().enumerate() {
            near(
                integrated.candidate.recession_m()[side][12],
                (K[side] / h) / gamma * removed,
                0.001,
            );
        }
        near(
            integrated.candidate.contact_state().pressure_pa()[12],
            pressure,
            0.005,
        );
        let compression = drive.position_m - integrated.candidate.contact_state().approach_m();
        close(
            integrated.receipt.path.normal_work_j,
            initial_spring_energy - 0.5 * drive.stiffness_n_per_m * compression * compression,
        );
        assert!(integrated.receipt.path.rejected_intervals > 0);
        accepted_accounts(&before, &integrated);
        assert_eq!(before, saved);
    });
}

#[test]
fn g1_body_hardness_refines_even_when_depth_pressure_and_caps_already_pass() {
    with_cx(|cx, _| {
        let pristine =
            Fixture::new(4.0e12).owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx);
        let before = pristine
            .advance_normal(compliance(&pristine) * 0.5e6, cx)
            .unwrap()
            .candidate;
        let distance = 1.0e-4;
        let coarse = before.advance_wear_integrated(&slip(distance), cx).unwrap();
        let first = before
            .advance_wear_integrated(&slip(0.5 * distance), cx)
            .unwrap();
        let fine = first
            .candidate
            .advance_wear_integrated(&slip(0.5 * distance), cx)
            .unwrap();
        let mut controls = accuracy(1.0);
        controls.absolute_depth_tolerance_m = 1.0e-3;
        controls.absolute_pressure_tolerance_pa = 5000.0;
        controls.relative_tolerance = 0.0;
        assert!(
            (coarse.candidate.contact_state().pressure_pa()[12]
                - fine.candidate.contact_state().pressure_pa()[12])
                .abs()
                < controls.absolute_pressure_tolerance_pa
        );
        assert_eq!(
            coarse.candidate.model().config().hardness_pa[12],
            fine.candidate.model().config().hardness_pa[12]
        );
        let profile = &before.hardness_profiles()[0];
        assert!(
            (profile
                .hardness_pa_at(coarse.candidate.recession_m()[0][12])
                .unwrap()
                - profile
                    .hardness_pa_at(fine.candidate.recession_m()[0][12])
                    .unwrap())
            .abs()
                > controls.absolute_pressure_tolerance_pa
        );
        let integrated = before
            .integrate_wear(&slip(distance), controls, cx)
            .unwrap();
        assert!(integrated.receipt.path.accepted_substeps > 2);
        assert!(
            integrated.receipt.maximum_accepted_local_hardness_error_pa
                <= controls.absolute_pressure_tolerance_pa
        );
        accepted_accounts(&before, &integrated);
    });
}

#[test]
fn g4_adaptive_path_budget_resolution_and_cancellation_refuse_atomically() {
    with_cx(|cx, gate| {
        let before = Fixture::new(GRADIENT)
            .owner(false, FiniteGapPlasticLinearSolver::Dense, 1.0e-4, cx)
            .advance_normal(5.0e-5, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let mut accepted_limit = accuracy(1.0);
        accepted_limit.maximum_accepted_substeps = 2;
        assert!(matches!(
            before.integrate_wear(&slip(1.0e-4), accepted_limit, cx),
            Err(DepthHardnessIntegrationError::BudgetExceeded {
                accepted_substeps: 2,
                ..
            })
        ));
        let mut evaluation_limit = accuracy(1.0);
        evaluation_limit.maximum_step_evaluations = 2;
        assert!(matches!(
            before.integrate_wear(&slip(1.0e-4), evaluation_limit, cx),
            Err(DepthHardnessIntegrationError::BudgetExceeded {
                step_evaluations: 2,
                ..
            })
        ));
        assert_eq!(before, saved);
        assert!(matches!(
            before.integrate_wear(&slip(f64::from_bits(1)), accuracy(1.0), cx),
            Err(DepthHardnessIntegrationError::ResolutionLimit { .. })
        ));
        for invalid in [vec![0.0; 24], slip(-1.0), slip(f64::NAN)] {
            assert!(matches!(
                before.integrate_wear(&invalid, accuracy(1.0), cx),
                Err(DepthHardnessIntegrationError::Invalid(_))
            ));
        }
        let mut invalid_controls = accuracy(1.0);
        invalid_controls.absolute_depth_tolerance_m = 0.0;
        assert!(matches!(
            before.integrate_wear(&[0.0; 25], invalid_controls, cx),
            Err(DepthHardnessIntegrationError::Invalid(_))
        ));
        // These fixed request errors cannot be repaired by reducing slip.
        // A deliberately oversized slip must not conceal them behind retries.
        assert!(matches!(
            before
                .integrate_wear_at_load(&slip(1.0), f64::NAN, load_controls(), accuracy(1.0), cx,),
            Err(DepthHardnessIntegrationError::Invalid(_))
        ));
        let mut outside = load_controls();
        outside.maximum_approach_m = 2.0e-4;
        assert!(matches!(
            before.integrate_wear_at_load(&slip(1.0), 0.01, outside, accuracy(1.0), cx),
            Err(DepthHardnessIntegrationError::Step(
                DepthHardnessWearError::Contact(_)
            ))
        ));
        assert!(matches!(
            before.integrate_wear_with_drive(
                &slip(1.0),
                FiniteGapPlasticDrive {
                    stiffness_n_per_m: 0.0,
                    position_m: 0.0
                },
                load_controls(),
                accuracy(1.0),
                cx,
            ),
            Err(DepthHardnessIntegrationError::Invalid(_))
        ));
        assert_eq!(before, saved);
        gate.request();
        assert!(matches!(
            before.integrate_wear(&slip(1.0e-4), accuracy(1.0), cx),
            Err(DepthHardnessIntegrationError::Step(_))
        ));
        assert_eq!(before, saved);
    });
}
