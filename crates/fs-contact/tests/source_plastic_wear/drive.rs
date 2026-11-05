//! G1/G3: finite-stiffness wear and spring-energy release against scalar oracles.

use super::*;
use fs_contact::normal_patch::plastic_source::drive::SourcePlasticDriveError;
use fs_contact::normal_patch::plastic_source::wear::integration::{
    WearIntegrationControls, WearIntegrationError,
};
use fs_contact::normal_patch::{
    FiniteGapPlasticDrive, FiniteGapPlasticDriveError, FiniteGapPlasticLoadControls,
};

fn controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

fn accuracy() -> WearIntegrationControls {
    WearIntegrationControls {
        absolute_depth_tolerance_m: 1.0e-10,
        absolute_pressure_tolerance_pa: 10.0,
        relative_tolerance: 1.0e-5,
        maximum_accepted_substeps: 4096,
        maximum_step_evaluations: 12000,
    }
}

fn compliance(source: &SourcePlasticArchardWear) -> f64 {
    let width = source.contact().model().config().grid.cell_width_m;
    let modulus = source
        .contact()
        .normal_source()
        .elastic()
        .reduced_modulus_pa();
    4.0 * width * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * modulus)
}

#[test]
fn wear_against_a_held_spring_matches_force_depth_and_energy_references() {
    with_cx(|cx, _| {
        let source = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let c = compliance(&source);
        let p0 = 0.5e6;
        let area = 1.0e-8;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: area / c,
            position_m: 2.0 * c * p0,
        };
        let loaded = source
            .advance_normal_with_drive(&source.initial_state(), drive, controls(), cx)
            .unwrap();
        let initial_spring_energy = loaded.drive.spring_energy_j;
        let before = loaded.candidate;
        let saved = before.clone();
        let result = source
            .advance_wear_with_drive(&before, &[1.0e-4; 25], drive, controls(), cx)
            .unwrap();
        let next = result.wear.candidate.contact_state().contact_state();
        let removal = 1.25e-6;
        let pressure = p0 - removal / (2.0 * c);
        close(result.drive.achieved_force_n, area * pressure);
        close(next.pressure_pa()[12], pressure);
        close(next.approach_m(), c * p0 + 0.5 * removal);
        close(next.recession_m()[12], removal);
        for side in 0..2 {
            let depth = [0.5e-6, 0.75e-6][side];
            close(result.wear.candidate.recession_m()[side][12], depth);
            close(result.wear.wear[side].increment_m3, area * depth);
            close(result.wear.candidate.volumes_m3()[side], area * depth);
        }
        let spring_release = 0.5 * area * c * (p0 * p0 - pressure * pressure);
        close(
            initial_spring_energy - result.drive.spring_energy_j,
            spring_release,
        );
        close(result.wear.receipt.normal_work_j, spring_release);
        close(result.wear.receipt.elastic_energy_change_j, -spring_release);
        close(result.wear.receipt.parameter_work_j, -2.0 * spring_release);
        assert_eq!(result.wear.receipt.response.plastic_work_j, 0.0);
        assert_eq!(
            next.plastic_depth_m(),
            before.contact_state().contact_state().plastic_depth_m()
        );
        assert_eq!(before, saved);

        let retained = source
            .advance_normal_with_drive(&result.wear.candidate, drive, controls(), cx)
            .unwrap();
        assert_eq!(retained.candidate, result.wear.candidate);
        let moved = source
            .advance_normal_with_drive(
                &result.wear.candidate,
                FiniteGapPlasticDrive {
                    position_m: drive.position_m * 0.75,
                    ..drive
                },
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            moved.candidate.volumes_m3(),
            result.wear.candidate.volumes_m3()
        );
        assert_eq!(
            moved.candidate.recession_m(),
            result.wear.candidate.recession_m()
        );
    });
}

#[test]
fn adaptive_spring_feedback_matches_finite_distance_exponential_decay() {
    with_cx(|cx, _| {
        let source = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let c = compliance(&source);
        let area = 1.0e-8;
        let beta = 0.01 / 1.0e6 + 0.03 / 2.0e6;
        let p0 = 0.5e6;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: area / c,
            position_m: 2.0 * c * p0,
        };
        let loaded = source
            .advance_normal_with_drive(&source.initial_state(), drive, controls(), cx)
            .unwrap();
        let initial_spring_energy = loaded.drive.spring_energy_j;
        let before = loaded.candidate;
        let saved = before.clone();
        let mut slip = [0.0; 25];
        // At fixed position, dp/ds = -beta*p/(C + A/K). With K=A/C,
        // s=4C/beta gives p=p0*exp(-2), including the load-frame compliance.
        slip[12] = 4.0 * c / beta;
        let coarse = source
            .advance_wear_with_drive(&before, &slip, drive, controls(), cx)
            .unwrap();
        assert_eq!(coarse.drive.achieved_force_n, 0.0);
        let expected_pressure = p0 * (-2.0_f64).exp();
        let expected_removal = 2.0 * c * p0 * (1.0 - (-2.0_f64).exp());
        let bounded = FiniteGapPlasticLoadControls {
            maximum_approach_m: c * p0 + 0.525 * expected_removal,
            ..controls()
        };
        assert!(matches!(
            source.advance_wear_with_drive(&before, &slip, drive, bounded, cx),
            Err(SourcePlasticWearError::Drive(
                SourcePlasticDriveError::Drive(
                    FiniteGapPlasticDriveError::TargetOutsideApproachBounds { .. }
                )
            ))
        ));
        let fine = source
            .integrate_wear_with_drive(&before, &slip, drive, bounded, accuracy(), cx)
            .unwrap();
        let next = fine.candidate.contact_state().contact_state();
        assert!((next.pressure_pa()[12] - expected_pressure).abs() < 0.01 * expected_pressure);
        assert!((next.recession_m()[12] - expected_removal).abs() < 0.01 * expected_removal);
        close(next.approach_m(), c * p0 + 0.5 * next.recession_m()[12]);
        assert!(fine.receipt.rejected_intervals > 0 && fine.receipt.accepted_substeps > 2);
        assert!(fine.receipt.maximum_accepted_normalized_error <= 1.0);
        for side in 0..2 {
            close(
                fine.receipt.removed_volume_m3[side],
                fine.candidate.volumes_m3()[side],
            );
            close(
                fine.candidate.volumes_m3()[side],
                [0.01 / 1.0e6, 0.03 / 2.0e6][side] * fine.receipt.pressure_weighted_sliding_j,
            );
        }
        close(
            fine.candidate.volumes_m3()[1] / fine.candidate.volumes_m3()[0],
            1.5,
        );
        let final_extension = drive.position_m - next.approach_m();
        let final_spring_energy = 0.5 * drive.stiffness_n_per_m * final_extension * final_extension;
        close(
            fine.receipt.normal_work_j,
            initial_spring_energy - final_spring_energy,
        );
        close(
            fine.receipt.elastic_energy_change_j + final_spring_energy - initial_spring_energy,
            fine.receipt.parameter_work_j,
        );
        assert_eq!(fine.receipt.plastic_work_j, 0.0);
        assert_eq!(fine.source_identity, source.identity());
        assert_eq!(before, saved);
    });
}

#[test]
fn spring_wear_refusals_and_exhausted_path_budget_are_transactional() {
    with_cx(|cx, gate| {
        let source = fixture([0.01, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        let c = compliance(&source);
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: 1.0e-8 / c,
            position_m: c * 1.0e6,
        };
        let before = source
            .advance_normal_with_drive(&source.initial_state(), drive, controls(), cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let changed = FiniteGapPlasticDrive {
            position_m: drive.position_m * 1.1,
            ..drive
        };
        assert!(matches!(
            source.advance_wear_with_drive(&before, &[1.0e-4; 25], changed, controls(), cx),
            Err(SourcePlasticWearError::Drive(
                SourcePlasticDriveError::Drive(
                    FiniteGapPlasticDriveError::PredecessorDriveMismatch { .. }
                )
            ))
        ));
        assert!(matches!(
            source.integrate_wear_with_drive(
                &before,
                &[1.0e-4; 25],
                drive,
                controls(),
                WearIntegrationControls {
                    maximum_step_evaluations: 1,
                    ..accuracy()
                },
                cx
            ),
            Err(WearIntegrationError::BudgetExceeded { .. })
        ));
        let zero = source
            .integrate_wear_with_drive(&before, &[0.0; 25], drive, controls(), accuracy(), cx)
            .unwrap();
        assert_eq!(zero.candidate, before);
        assert_eq!(zero.receipt.normal_work_j, 0.0);
        let foreign = fixture([0.02, 0.03], [1.0e6, 2.0e6], [1.0e-4; 2], false, cx);
        assert!(
            foreign
                .advance_normal_with_drive(&before, drive, controls(), cx)
                .is_err()
        );
        assert!(
            source
                .advance_wear_with_drive(&before, &[1.0e-4; 24], drive, controls(), cx)
                .is_err()
        );
        gate.request();
        assert!(
            source
                .integrate_wear_with_drive(
                    &before,
                    &[1.0e-4; 25],
                    drive,
                    controls(),
                    accuracy(),
                    cx
                )
                .is_err()
        );
        assert_eq!(before, saved);
    });
}
