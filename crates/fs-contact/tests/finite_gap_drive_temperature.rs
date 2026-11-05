//! G1/G3 material continuation under a finite, stationary loading device.
//! Original permanent indentation and removed surface depth remain independent.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapContactError, FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticContinuationKind,
    FiniteGapPlasticDrive, FiniteGapPlasticDriveDisposition, FiniteGapPlasticDriveError,
    FiniteGapPlasticError, FiniteGapPlasticLoadControls, FiniteGapPlasticModel,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 103,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn config() -> FiniteGapPlasticConfig {
    let mut gap = vec![0.017; 25];
    gap[12] = 0.007;
    FiniteGapPlasticConfig {
        grid: FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: 1e-3,
            cell_depth_m: 1e-3,
        },
        undeformed_gap_m: gap,
        reduced_modulus_pa: 1e9,
        hardness_pa: vec![1e6; 25],
        maximum_active_set_iterations: 256,
        complementarity_tolerance_m: 1e-12,
        boundary_clearance_cells: 1,
    }
}

fn controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1e-5,
        maximum_approach_m: 1e-5,
        absolute_force_tolerance_n: 1e-13,
        relative_force_tolerance: 1e-12,
        approach_tolerance_m: 1e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

fn compliance(config: &FiniteGapPlasticConfig) -> f64 {
    4.0 * config.grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (core::f64::consts::PI * config.reduced_modulus_pa)
}

fn drive(position_m: f64) -> FiniteGapPlasticDrive {
    FiniteGapPlasticDrive {
        stiffness_n_per_m: 2e5,
        position_m,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

#[test]
fn compliance_change_uses_both_springs_and_preserves_dent_recession_through_reversal() {
    with_cx(|cx, _| {
        let original = config();
        let c = compliance(&original);
        let a = 1e-6;
        let previous = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let dented = previous
            .advance(&previous.initial_state(), c * 1e6 + 2e-6, cx)
            .unwrap();
        let mut removal = vec![0.0; 25];
        removal[12] = 3e-7;
        let worn = previous.recede(&dented.candidate, &removal, cx).unwrap();
        let old_gap = worn.candidate.plastic_depth_m()[12] + worn.candidate.recession_m()[12];
        let boundary = drive(old_gap + 0.4 * (c / a + 1.0 / 2e5));
        let loaded = previous
            .advance_with_drive(&worn.candidate, boundary, controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let saved = before.clone();
        let mut softer = original;
        softer.reduced_modulus_pa *= 0.5;
        let successor = FiniteGapPlasticModel::new(softer, cx).unwrap();
        let next = successor
            .refresh_from_with_drive(&previous, &before, boundary, controls(), cx)
            .unwrap();
        let expected_force =
            (boundary.position_m - old_gap) / (2.0 * c / a + 1.0 / boundary.stiffness_n_per_m);
        let expected_approach = old_gap + 2.0 * c * expected_force / a;
        let expected_work =
            0.5 * (0.4 + expected_force) * (expected_approach - before.approach_m());
        let old_contact_energy = 0.5 * c * 0.4 * 0.4 / a;
        let new_contact_energy = c * expected_force * expected_force / a;
        let receipt = &next.continuation.receipt;
        close(next.drive.achieved_force_n, expected_force, 4e-12);
        close(
            next.continuation.candidate.approach_m(),
            expected_approach,
            4e-12,
        );
        close(receipt.normal_work_j, expected_work, 4e-11);
        close(
            receipt.elastic_energy_change_j,
            new_contact_energy - old_contact_energy,
            4e-11,
        );
        close(
            receipt.parameter_work_j,
            new_contact_energy - old_contact_energy - expected_work,
            4e-11,
        );
        close(
            loaded.drive.spring_energy_j - next.drive.spring_energy_j,
            expected_work,
            4e-11,
        );
        assert_eq!(
            next.continuation.candidate.plastic_depth_m(),
            before.plastic_depth_m()
        );
        assert_eq!(
            next.continuation.candidate.recession_m(),
            before.recession_m()
        );
        assert_eq!(
            next.continuation.candidate.cumulative_plastic_work_j(),
            before.cumulative_plastic_work_j()
        );
        assert_eq!(receipt.response.plastic_work_j, 0.0);
        assert_eq!(receipt.recession_increment_m, vec![0.0; 25]);
        assert_eq!(
            receipt.kind,
            FiniteGapPlasticContinuationKind::MaterialRefreshWithDrive
        );
        assert_eq!(receipt.predecessor_model_identity, previous.identity());
        assert_eq!(receipt.successor_model_identity, successor.identity());
        assert_eq!(receipt.predecessor_state_identity, before.identity());
        assert_eq!(receipt.response.predecessor_identity, before.identity());
        assert_eq!(next.drive.predecessor_state_identity, before.identity());
        assert_eq!(next.drive.model_identity, successor.identity());
        assert_eq!(next.drive.response_identity, receipt.response.identity);
        assert_eq!(receipt.response.normal.undeformed_gap_datum_m, 0.007);
        assert_eq!(before, saved);
        let restored = previous
            .refresh_from_with_drive(
                &successor,
                &next.continuation.candidate,
                boundary,
                controls(),
                cx,
            )
            .unwrap();
        close(restored.drive.achieved_force_n, 0.4, 4e-12);
        close(
            restored.continuation.candidate.approach_m(),
            before.approach_m(),
            4e-12,
        );
        assert_eq!(
            restored.continuation.candidate.plastic_depth_m(),
            before.plastic_depth_m()
        );
        assert_eq!(
            restored.continuation.candidate.recession_m(),
            before.recession_m()
        );
        assert_eq!(
            restored.continuation.candidate.cumulative_plastic_work_j(),
            before.cumulative_plastic_work_j()
        );
    });
}

#[test]
fn lower_hardness_yields_under_the_loading_spring_and_cooling_does_not_heal() {
    with_cx(|cx, _| {
        let original = config();
        let c = compliance(&original);
        let previous = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let boundary = drive(0.6 * (c / 1e-6 + 1.0 / 2e5));
        let loaded = previous
            .advance_with_drive(&previous.initial_state(), boundary, controls(), cx)
            .unwrap();
        let mut softened = original;
        softened.hardness_pa[12] = 3e5;
        let successor = FiniteGapPlasticModel::new(softened, cx).unwrap();
        let next = successor
            .refresh_from_with_drive(
                &previous,
                &loaded.contact.candidate,
                boundary,
                controls(),
                cx,
            )
            .unwrap();
        let expected_force = 0.3;
        let expected_approach = boundary.position_m - expected_force / boundary.stiffness_n_per_m;
        let expected_dent = expected_approach - c * 3e5;
        assert!(expected_dent > 0.0);
        close(next.drive.achieved_force_n, expected_force, 4e-12);
        close(
            next.continuation.candidate.approach_m(),
            expected_approach,
            4e-12,
        );
        close(
            next.continuation.candidate.plastic_depth_m()[12],
            expected_dent,
            4e-12,
        );
        close(
            next.continuation.receipt.response.plastic_work_j,
            0.3 * expected_dent,
            4e-12,
        );
        let cooled = previous
            .refresh_from_with_drive(
                &successor,
                &next.continuation.candidate,
                boundary,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            cooled.continuation.candidate.plastic_depth_m(),
            next.continuation.candidate.plastic_depth_m()
        );
        assert_eq!(
            cooled.continuation.candidate.cumulative_plastic_work_j(),
            next.continuation.candidate.cumulative_plastic_work_j()
        );
        assert_eq!(cooled.continuation.receipt.response.plastic_work_j, 0.0);
        close(cooled.drive.achieved_force_n, expected_force, 4e-12);
        close(
            cooled.continuation.candidate.approach_m(),
            expected_approach,
            4e-12,
        );
    });
}

#[test]
fn direct_material_refresh_does_not_preserve_a_phantom_fixed_approach_dent() {
    with_cx(|cx, _| {
        let original = config();
        let c = compliance(&original);
        let previous = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let boundary = drive(0.4 * (c / 1e-6 + 1.0 / 2e5));
        let loaded = previous
            .advance_with_drive(&previous.initial_state(), boundary, controls(), cx)
            .unwrap();
        let mut changed = original;
        changed.reduced_modulus_pa *= 4.0;
        changed.hardness_pa[12] = 5e5;
        let successor = FiniteGapPlasticModel::new(changed, cx).unwrap();
        let fixed_approach = successor
            .refresh_from(&previous, &loaded.contact.candidate, cx)
            .unwrap();
        assert!(fixed_approach.candidate.plastic_depth_m()[12] > 0.0);
        let next = successor
            .refresh_from_with_drive(
                &previous,
                &loaded.contact.candidate,
                boundary,
                controls(),
                cx,
            )
            .unwrap();
        let expected_force =
            boundary.position_m / (0.25 * c / 1e-6 + 1.0 / boundary.stiffness_n_per_m);
        assert!(expected_force < 0.5);
        close(next.drive.achieved_force_n, expected_force, 4e-12);
        close(
            next.continuation.candidate.approach_m(),
            0.25 * c * expected_force / 1e-6,
            4e-12,
        );
        assert_eq!(next.continuation.candidate.plastic_depth_m(), &[0.0; 25]);
        assert_eq!(next.continuation.receipt.response.plastic_work_j, 0.0);
        assert_eq!(next.continuation.candidate.cumulative_plastic_work_j(), 0.0);
        assert!(next.continuation.receipt.normal_work_j < 0.0);
        assert!(next.drive.spring_energy_j > loaded.drive.spring_energy_j);
    });
}

#[test]
fn maintained_drive_material_refresh_retains_noops_and_refuses_foreign_or_changed_boundaries() {
    with_cx(|cx, gate| {
        let original = config();
        let previous = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let boundary = drive(3e-6);
        let loaded = previous
            .advance_with_drive(&previous.initial_state(), boundary, controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let saved = before.clone();
        let unchanged = previous
            .refresh_from_with_drive(&previous, &before, boundary, controls(), cx)
            .unwrap();
        assert_eq!(unchanged.continuation.candidate, before);
        assert_eq!(unchanged.drive.equilibrium_evaluations, 0);
        assert_eq!(
            unchanged.drive.disposition,
            FiniteGapPlasticDriveDisposition::RetainedExistingEquilibrium
        );
        assert_eq!(unchanged.continuation.receipt.normal_work_j, 0.0);
        assert_eq!(unchanged.continuation.receipt.elastic_energy_change_j, 0.0);
        assert_eq!(unchanged.continuation.receipt.parameter_work_j, 0.0);
        // Force underflow must not turn a finite spring into a fictitious slack
        // interval or retain the virgin approach instead of its drive position.
        let vanishing_drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: f64::MIN_POSITIVE,
            position_m: 1e-100,
        };
        assert!(matches!(
            previous.advance_with_drive(&previous.initial_state(), vanishing_drive, controls(), cx),
            Err(FiniteGapPlasticDriveError::InvalidInput { .. })
        ));
        let mut changed = original.clone();
        changed.reduced_modulus_pa *= 2.0;
        let successor = FiniteGapPlasticModel::new(changed, cx).unwrap();
        assert!(matches!(
            successor.refresh_from_with_drive(&previous, &before, drive(4e-6), controls(), cx),
            Err(FiniteGapPlasticDriveError::PredecessorDriveMismatch { .. })
        ));
        assert_eq!(
            successor.refresh_from_with_drive(&successor, &before, boundary, controls(), cx),
            Err(FiniteGapPlasticDriveError::Contact(
                FiniteGapPlasticError::ForeignState
            ))
        );
        let mut moved = original;
        moved.undeformed_gap_m[0] += 1e-4;
        let wrong_geometry = FiniteGapPlasticModel::new(moved, cx).unwrap();
        assert!(matches!(
            wrong_geometry.refresh_from_with_drive(&previous, &before, boundary, controls(), cx),
            Err(FiniteGapPlasticDriveError::Contact(
                FiniteGapPlasticError::Contact(FiniteGapContactError::InvalidInput { .. })
            ))
        ));
        let mut exhausted = controls();
        exhausted.maximum_equilibrium_evaluations = 1;
        assert!(matches!(
            successor.refresh_from_with_drive(&previous, &before, boundary, exhausted, cx),
            Err(FiniteGapPlasticDriveError::DidNotConverge { .. })
        ));
        assert_eq!(before, saved);
        gate.request();
        assert!(matches!(
            successor.refresh_from_with_drive(&previous, &before, boundary, controls(), cx),
            Err(FiniteGapPlasticDriveError::Contact(
                FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled)
            ))
        ));
        assert_eq!(before, saved);
    });
}
