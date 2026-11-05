//! Finite-distance pressure feedback and independent scalar/matrix wear oracles.

use super::*;
use fs_contact::normal_patch::plastic_source::load::SourcePlasticLoadError;
use fs_contact::normal_patch::plastic_source::wear::integration::{
    WearIntegrationControls, WearIntegrationError, WearIntegrationReceipt,
};
use fs_contact::normal_patch::{FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError};

fn accurate() -> WearIntegrationControls {
    WearIntegrationControls {
        absolute_depth_tolerance_m: 1.0e-10,
        absolute_pressure_tolerance_pa: 10.0,
        relative_tolerance: 1.0e-5,
        maximum_accepted_substeps: 4096,
        maximum_step_evaluations: 12000,
    }
}
fn load_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-5,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-12,
        relative_force_tolerance: 1.0e-10,
        approach_tolerance_m: 1.0e-14,
        maximum_equilibrium_evaluations: 128,
    }
}
fn self_compliance(model: &SourcePlasticArchardWear) -> f64 {
    let width = model.contact().model().config().grid.cell_width_m;
    let modulus = model
        .contact()
        .normal_source()
        .elastic()
        .reduced_modulus_pa();
    4.0 * width * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * modulus)
}
fn neighbor_compliance(model: &SourcePlasticArchardWear) -> f64 {
    // Exact rectangular Boussinesq integral at the adjacent square's centre:
    // 2*[I(3w/2,w/2)-I(w/2,w/2)], I(a,b)=a asinh(b/a)+b asinh(a/b).
    let width = model.contact().model().config().grid.cell_width_m;
    let modulus = model
        .contact()
        .normal_source()
        .elastic()
        .reduced_modulus_pa();
    width * (3.0 * (1.0_f64 / 3.0).asinh() + 3.0_f64.asinh() - 2.0 * 1.0_f64.asinh())
        / (core::f64::consts::PI * modulus)
}
fn closes(receipt: &WearIntegrationReceipt) {
    let scale = receipt.elastic_energy_change_j.abs()
        + receipt.plastic_work_j.abs()
        + receipt.normal_work_j.abs()
        + receipt.parameter_work_j.abs();
    assert!(receipt.closure_residual_j.abs() <= 128.0 * f64::EPSILON * scale.max(1.0e-30));
    assert!(receipt.maximum_accepted_normalized_error <= 1.0);
}

#[test]
fn adaptive_held_approach_preserves_finite_distance_contact_and_converges() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let before = model
            .advance_normal(&model.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let original = before.clone();
        let c = self_compliance(&model);
        let beta = 0.01 / 1.0e6 + 0.03 / 2.0e6;
        let distance = 2.0 * c / beta;
        let mut slip = vec![0.0; 25];
        slip[12] = distance;
        let one_step = model.advance_wear(&before, &slip, cx).unwrap();
        assert_eq!(
            one_step
                .candidate
                .contact_state()
                .contact_state()
                .pressure_pa()[12],
            0.0
        );

        let mut loose_controls = accurate();
        loose_controls.absolute_depth_tolerance_m = 1.0e-8;
        loose_controls.absolute_pressure_tolerance_pa = 1000.0;
        loose_controls.relative_tolerance = 1.0e-3;
        let loose = model
            .integrate_wear(&before, &slip, loose_controls, cx)
            .unwrap();
        let fine = model
            .integrate_wear(&before, &slip, accurate(), cx)
            .unwrap();
        let p0 = before.contact_state().contact_state().pressure_pa()[12];
        let exact_pressure = p0 * (-2.0_f64).exp();
        let loose_pressure = loose
            .candidate
            .contact_state()
            .contact_state()
            .pressure_pa()[12];
        let fine_pressure = fine.candidate.contact_state().contact_state().pressure_pa()[12];
        assert!(fine_pressure > 0.0);
        assert!((fine_pressure - exact_pressure).abs() < 0.01 * exact_pressure);
        assert!((fine_pressure - exact_pressure).abs() < (loose_pressure - exact_pressure).abs());
        assert!(fine.receipt.accepted_substeps > loose.receipt.accepted_substeps);
        assert!(fine.receipt.rejected_intervals > 0);
        let exact_relative_recession = c * p0 * (1.0 - (-2.0_f64).exp());
        for side in 0..2 {
            let fraction = [0.01 / 1.0e6, 0.03 / 2.0e6][side] / beta;
            let exact_volume = fraction * exact_relative_recession * 1.0e-8;
            assert!((fine.candidate.volumes_m3()[side] - exact_volume).abs() < 0.01 * exact_volume);
            close(
                fine.receipt.removed_volume_m3[side],
                fine.candidate.volumes_m3()[side],
            );
            close(
                fine.receipt.removed_volume_m3[side],
                [0.01 / 1.0e6, 0.03 / 2.0e6][side] * fine.receipt.pressure_weighted_sliding_j,
            );
        }
        close(
            fine.candidate.volumes_m3()[1] / fine.candidate.volumes_m3()[0],
            1.5,
        );
        assert_eq!(
            fine.candidate
                .contact_state()
                .contact_state()
                .plastic_depth_m(),
            before.contact_state().contact_state().plastic_depth_m()
        );
        assert_eq!(fine.receipt.plastic_work_j, 0.0);
        assert_eq!(fine.receipt.normal_work_j, 0.0);
        assert!(fine.receipt.parameter_work_j < 0.0);
        assert_eq!(fine.source_identity, model.identity());
        assert_eq!(fine.receipt.controls, accurate());
        closes(&fine.receipt);
        assert_eq!(before, original);
    });
}

#[test]
fn held_load_nonuniform_sliding_matches_two_cell_pressure_redistribution() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], true, cx);
        let force = 0.006;
        let before = model
            .advance_normal_at_load(&model.initial_state(), force, load_controls(), cx)
            .unwrap()
            .candidate;
        let difference = self_compliance(&model) - neighbor_compliance(&model);
        let beta = 0.01 / 1.0e6 + 0.03 / 2.0e6;
        let distance = 4.0 * difference / beta;
        let mut slip = vec![0.0; 25];
        slip[12] = distance;
        let p0 = force / (2.0 * 1.0e-8);
        let exact_recession = 2.0 * difference * p0 * (1.0 - (-2.0_f64).exp());
        let initial_approach = before.contact_state().contact_state().approach_m();
        let mut bounded_load = load_controls();
        // The exact final travel is half the recession. A five-percent travel
        // margin admits it while rejecting the oversized beginning-pressure step.
        bounded_load.maximum_approach_m = initial_approach + 0.525 * exact_recession;
        assert!(matches!(
            model.advance_wear_at_load(&before, &slip, force, bounded_load, cx),
            Err(SourcePlasticWearError::Load(SourcePlasticLoadError::Load(
                FiniteGapPlasticLoadError::TargetOutsideApproachBounds { .. }
            )))
        ));
        let result = model
            .integrate_wear_at_load(&before, &slip, force, bounded_load, accurate(), cx)
            .unwrap();
        // Equality of the two approach equations gives
        // p_12' = -beta*p_12/[2*(C_self-C_neighbor)] at fixed total load.
        let exact_p12 = p0 * (-2.0_f64).exp();
        let native = result.candidate.contact_state().contact_state();
        assert!(native.pressure_pa()[12] > 0.0);
        assert!((native.pressure_pa()[12] - exact_p12).abs() < 0.01 * exact_p12);
        close(
            (native.pressure_pa()[12] + native.pressure_pa()[13]) * 1.0e-8,
            force,
        );
        assert!(native.pressure_pa()[13] > p0);
        assert!(result.receipt.accepted_substeps > 2);
        assert_eq!(result.candidate.recession_m()[0][13], 0.0);
        assert_eq!(result.candidate.recession_m()[1][13], 0.0);
        assert!((native.recession_m()[12] - exact_recession).abs() < 0.01 * exact_recession);
        let travel = native.approach_m() - initial_approach;
        assert!((travel - 0.5 * exact_recession).abs() < 0.005 * exact_recession);
        assert!(result.receipt.rejected_intervals > 0);
        assert_eq!(result.receipt.controls, accurate());
        let actuator_work =
            force * (native.approach_m() - before.contact_state().contact_state().approach_m());
        close(result.receipt.normal_work_j, actuator_work);
        assert!(result.receipt.normal_work_j > 0.0);
        assert_eq!(result.receipt.plastic_work_j, 0.0);
        closes(&result.receipt);
        assert!(matches!(
            model.integrate_wear_at_load(
                &before,
                &slip,
                force * 1.1,
                load_controls(),
                accurate(),
                cx,
            ),
            Err(WearIntegrationError::Step(SourcePlasticWearError::Load(
                SourcePlasticLoadError::Load(
                    FiniteGapPlasticLoadError::PredecessorLoadMismatch { .. }
                )
            )))
        ));
        for side in 0..2 {
            close(
                result.receipt.removed_volume_m3[side],
                result.candidate.volumes_m3()[side],
            );
        }
    });
}

#[test]
fn uniform_single_cell_load_has_exact_archard_volume_and_actuator_work() {
    with_cx(|cx, _| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let force = 0.005;
        let before = model
            .advance_normal_at_load(&model.initial_state(), force, load_controls(), cx)
            .unwrap()
            .candidate;
        let mut slip = vec![0.0; 25];
        slip[12] = 1.0e-4;
        let result = model
            .integrate_wear_at_load(&before, &slip, force, load_controls(), accurate(), cx)
            .unwrap();
        let native = result.candidate.contact_state().contact_state();
        close(native.pressure_pa()[12] * 1.0e-8, force);
        close(result.receipt.pressure_weighted_sliding_j, force * slip[12]);
        for side in 0..2 {
            close(
                result.receipt.removed_volume_m3[side],
                [0.01 / 1.0e6, 0.03 / 2.0e6][side] * force * slip[12],
            );
        }
        close(
            result.receipt.normal_work_j,
            force * native.recession_m()[12],
        );
        close(
            result.receipt.parameter_work_j,
            -result.receipt.normal_work_j,
        );
        assert_eq!(result.receipt.plastic_work_j, 0.0);
        closes(&result.receipt);
    });
}

#[test]
fn refinement_recovers_from_oversized_trial_and_budgets_discard_partial_paths() {
    with_cx(|cx, _| {
        let reference = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let c = self_compliance(&reference);
        let closure = c * 1.0e6;
        let model = fixture(
            [0.01, 0.03],
            [1.0e6, 2.0e6],
            [0.5 * closure, 1.0e-4],
            false,
            cx,
        );
        let before = model
            .advance_normal(&model.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let mut slip = vec![0.0; 25];
        slip[12] = 2.0 * c / (0.01 / 1.0e6 + 0.03 / 2.0e6);
        assert!(model.advance_wear(&before, &slip, cx).is_err());
        let result = model
            .integrate_wear(&before, &slip, accurate(), cx)
            .unwrap();
        assert!(result.candidate.recession_m()[0][12] < 0.5 * closure);
        assert!(result.receipt.rejected_intervals > 0);
        let mut limited = accurate();
        limited.maximum_accepted_substeps = 2;
        assert!(matches!(
            model.integrate_wear(&before, &slip, limited, cx),
            Err(WearIntegrationError::BudgetExceeded {
                accepted_substeps: 2,
                ..
            })
        ));
        let mut limited = accurate();
        limited.maximum_step_evaluations = 2;
        assert!(matches!(
            model.integrate_wear(&before, &slip, limited, cx),
            Err(WearIntegrationError::BudgetExceeded {
                step_evaluations: 2,
                ..
            })
        ));
        assert_eq!(before, saved);
    });
}

#[test]
fn zero_path_source_mismatch_invalid_controls_and_cancellation_keep_history() {
    with_cx(|cx, gate| {
        let model = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let before = model
            .advance_normal(&model.initial_state(), 5.0e-5, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let result = model
            .integrate_wear(&before, &[0.0; 25], accurate(), cx)
            .unwrap();
        assert_eq!(result.candidate, before);
        assert_eq!(result.receipt.removed_volume_m3, [0.0; 2]);
        assert_eq!(result.receipt.normal_work_j, 0.0);
        assert_eq!(result.receipt.parameter_work_j, 0.0);
        let other = fixture([0.02, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        assert!(
            other
                .integrate_wear(&before, &[0.0; 25], accurate(), cx)
                .is_err()
        );
        let mut invalid = accurate();
        invalid.absolute_depth_tolerance_m = f64::NAN;
        assert!(matches!(
            model.integrate_wear(&before, &[0.0; 25], invalid, cx),
            Err(WearIntegrationError::Invalid(_))
        ));
        assert!(
            model
                .integrate_wear(&before, &[0.0; 24], accurate(), cx)
                .is_err()
        );
        gate.request();
        assert!(
            model
                .integrate_wear(&before, &[0.0; 25], accurate(), cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
