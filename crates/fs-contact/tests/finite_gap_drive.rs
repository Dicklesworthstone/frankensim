//! G1 series-spring contact references and G3 irreversible-history checks.
//! The cell compliance is evaluated independently from the rectangular integral.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapContactError, FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticContinuationKind,
    FiniteGapPlasticDrive, FiniteGapPlasticDriveDisposition, FiniteGapPlasticDriveError,
    FiniteGapPlasticError, FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError,
    FiniteGapPlasticModel,
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
                seed: 37,
                kernel_id: 91,
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
        maximum_approach_m: 2e-5,
        absolute_force_tolerance_n: 1e-13,
        relative_force_tolerance: 1e-12,
        approach_tolerance_m: 1e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

fn drive(position_m: f64) -> FiniteGapPlasticDrive {
    FiniteGapPlasticDrive {
        stiffness_n_per_m: 2e5,
        position_m,
    }
}

fn self_compliance(config: &FiniteGapPlasticConfig) -> f64 {
    4.0 * config.grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (core::f64::consts::PI * config.reduced_modulus_pa)
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

#[test]
fn elastic_series_spring_discards_the_yielding_upper_trial() {
    with_cx(|cx, _| {
        let config = config();
        let c = self_compliance(&config);
        let a = config.grid.cell_width_m * config.grid.cell_depth_m;
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let before = model.initial_state();
        let saved = before.clone();
        let loading = drive(3e-6);
        assert!(
            model
                .advance(&before, loading.position_m, cx)
                .unwrap()
                .receipt
                .plastic_work_j
                > 0.0
        );
        let next = model
            .advance_with_drive(&before, loading, controls(), cx)
            .unwrap();
        let p =
            loading.stiffness_n_per_m * loading.position_m / (a + loading.stiffness_n_per_m * c);
        close(next.contact.candidate.pressure_pa()[12], p, 3e-12);
        close(next.contact.candidate.approach_m(), c * p, 3e-12);
        close(next.drive.achieved_force_n, a * p, 3e-12);
        close(next.drive.drive_force_n, a * p, 3e-12);
        close(
            next.drive.spring_energy_j,
            0.5 * a * p * (loading.position_m - c * p),
            3e-12,
        );
        assert_eq!(next.contact.candidate.plastic_depth_m(), &[0.0; 25]);
        assert_eq!(next.contact.receipt.plastic_work_j, 0.0);
        assert_eq!(next.contact.candidate.cumulative_plastic_work_j(), 0.0);
        assert_eq!(next.contact.receipt.normal.undeformed_gap_datum_m, 0.007);
        assert_eq!(next.drive.predecessor_state_identity, before.identity());
        assert_eq!(
            next.drive.candidate_state_identity,
            next.contact.candidate.identity()
        );
        assert_eq!(next.drive.response_identity, next.contact.receipt.identity);
        assert!(next.drive.force_residual_n.abs() <= next.drive.force_tolerance_n);
        assert_eq!(before, saved);
    });
}

#[test]
fn finite_stiffness_selects_motion_through_the_perfectly_plastic_plateau() {
    with_cx(|cx, _| {
        let config = config();
        let c = self_compliance(&config);
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let before = model.initial_state();
        assert!(matches!(
            model.advance_at_load(&before, 1.0, controls(), cx),
            Err(FiniteGapPlasticLoadError::NonUniqueLoad { .. })
        ));
        let loading = drive(1e-5);
        let next = model
            .advance_with_drive(&before, loading, controls(), cx)
            .unwrap();
        let approach = loading.position_m - 1.0 / loading.stiffness_n_per_m;
        close(next.contact.candidate.approach_m(), approach, 3e-12);
        assert_eq!(next.contact.candidate.pressure_pa()[12], 1e6);
        close(
            next.contact.candidate.plastic_depth_m()[12],
            approach - c * 1e6,
            3e-12,
        );
        let further = model
            .advance_with_drive(&next.contact.candidate, drive(1.2e-5), controls(), cx)
            .unwrap();
        close(
            further.contact.candidate.approach_m() - approach,
            2e-6,
            3e-12,
        );
        close(further.contact.receipt.plastic_increment_m[12], 2e-6, 3e-12);
        close(further.contact.receipt.plastic_work_j, 2e-6, 3e-12);
        assert_eq!(further.drive.achieved_force_n, 1.0);
        assert_eq!(
            further.contact.receipt.predecessor_identity,
            next.contact.candidate.identity()
        );
    });
}

#[test]
fn unloading_retains_the_dent_and_detachment_selects_the_unloaded_drive_position() {
    with_cx(|cx, _| {
        let config = config();
        let c = self_compliance(&config);
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let loaded = model
            .advance_with_drive(&model.initial_state(), drive(1e-5), controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let approach = before.plastic_depth_m()[12] + c * 0.4 / 1e-6;
        let unloading = drive(approach + 0.4 / 2e5);
        let next = model
            .advance_with_drive(&before, unloading, controls(), cx)
            .unwrap();
        close(next.drive.achieved_force_n, 0.4, 3e-12);
        close(next.contact.candidate.approach_m(), approach, 3e-12);
        assert_eq!(
            next.contact.candidate.plastic_depth_m(),
            before.plastic_depth_m()
        );
        assert_eq!(
            next.contact.candidate.cumulative_plastic_work_j(),
            before.cumulative_plastic_work_j()
        );
        let open = model
            .advance_with_drive(&next.contact.candidate, drive(-2e-6), controls(), cx)
            .unwrap();
        assert_eq!(open.contact.candidate.approach_m(), -2e-6);
        assert_eq!(open.drive.achieved_force_n, 0.0);
        assert_eq!(open.drive.drive_force_n, 0.0);
        assert_eq!(open.drive.spring_energy_j, 0.0);
        assert_eq!(
            open.contact.candidate.plastic_depth_m(),
            before.plastic_depth_m()
        );
    });
}

#[test]
fn repeated_drive_and_zero_recession_preserve_the_complete_plastic_state() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let boundary = drive(1e-5);
        let first = model
            .advance_with_drive(&model.initial_state(), boundary, controls(), cx)
            .unwrap();
        let before = first.contact.candidate;
        let repeated = model
            .advance_with_drive(&before, boundary, controls(), cx)
            .unwrap();
        assert_eq!(repeated.contact.candidate, before);
        assert_eq!(repeated.contact.receipt.plastic_work_j, 0.0);
        assert_eq!(repeated.drive.equilibrium_evaluations, 0);
        assert_eq!(
            repeated.drive.disposition,
            FiniteGapPlasticDriveDisposition::RetainedExistingEquilibrium
        );
        let zero = model
            .recede_with_drive(&before, &[0.0; 25], boundary, controls(), cx)
            .unwrap();
        assert_eq!(zero.continuation.candidate, before);
        assert_eq!(zero.drive.equilibrium_evaluations, 0);
        assert_eq!(zero.continuation.receipt.normal_work_j, 0.0);
        assert_eq!(zero.continuation.receipt.elastic_energy_change_j, 0.0);
        assert_eq!(zero.continuation.receipt.parameter_work_j, 0.0);
    });
}

#[test]
fn fixed_drive_wear_releases_spring_storage_without_actuator_work_or_healing() {
    with_cx(|cx, _| {
        let config = config();
        let c = self_compliance(&config);
        let a = 1e-6;
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let dented = model
            .advance_with_drive(&model.initial_state(), drive(1e-5), controls(), cx)
            .unwrap();
        let dent = dented.contact.candidate.plastic_depth_m()[12];
        let boundary = drive(dent + c * 0.5 / a + 0.5 / 2e5);
        let loaded = model
            .advance_with_drive(&dented.contact.candidate, boundary, controls(), cx)
            .unwrap();
        let before = &loaded.contact.candidate;
        let mut removal = vec![0.0; 25];
        removal[12] = 2e-7;
        let next = model
            .recede_with_drive(before, &removal, boundary, controls(), cx)
            .unwrap();
        let expected_force =
            (boundary.position_m - dent - removal[12]) / (c / a + 1.0 / boundary.stiffness_n_per_m);
        let expected_approach = dent + removal[12] + c * expected_force / a;
        let expected_work =
            0.5 * (0.5 + expected_force) * (expected_approach - before.approach_m());
        close(next.drive.achieved_force_n, expected_force, 4e-12);
        close(
            next.continuation.candidate.approach_m(),
            expected_approach,
            4e-12,
        );
        close(
            next.continuation.receipt.normal_work_j,
            expected_work,
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
            next.continuation.candidate.cumulative_plastic_work_j(),
            before.cumulative_plastic_work_j()
        );
        assert_eq!(next.continuation.candidate.recession_m(), removal);
        assert_eq!(
            next.continuation.receipt.kind,
            FiniteGapPlasticContinuationKind::SurfaceRecessionWithDrive
        );
        let receipt = &next.continuation.receipt;
        close(
            receipt.parameter_work_j,
            receipt.elastic_energy_change_j - expected_work,
            4e-11,
        );
        close(
            receipt.parameter_work_j,
            receipt.elastic_energy_change_j + next.drive.spring_energy_j
                - loaded.drive.spring_energy_j,
            4e-11,
        );
        assert_eq!(receipt.response.plastic_work_j, 0.0);
        assert_eq!(receipt.predecessor_state_identity, before.identity());
    });
}

#[test]
fn complete_recession_can_unload_both_contact_and_the_finite_loading_spring() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let boundary = drive(3e-6);
        let loaded = model
            .advance_with_drive(&model.initial_state(), boundary, controls(), cx)
            .unwrap();
        let mut removal = vec![0.0; 25];
        removal[12] = 4e-6;
        let next = model
            .recede_with_drive(
                &loaded.contact.candidate,
                &removal,
                boundary,
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            next.continuation.candidate.approach_m(),
            boundary.position_m
        );
        assert_eq!(next.drive.achieved_force_n, 0.0);
        assert_eq!(next.drive.drive_force_n, 0.0);
        assert_eq!(next.drive.spring_energy_j, 0.0);
        close(
            next.continuation.receipt.normal_work_j,
            loaded.drive.spring_energy_j,
            4e-12,
        );
        close(
            next.continuation.receipt.parameter_work_j,
            -loaded.contact.receipt.normal.reversible_energy_j - loaded.drive.spring_energy_j,
            4e-12,
        );
    });
}

#[test]
fn direct_pressure_solve_resolves_drive_changes_smaller_than_closure_tolerance() {
    with_cx(|cx, _| {
        let mut config = config();
        config.complementarity_tolerance_m = 1e-8;
        let c = self_compliance(&config);
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let initial = model
            .advance_with_drive(&model.initial_state(), drive(3e-6), controls(), cx)
            .unwrap();
        let boundary = drive(3e-6 + 1e-9);
        let next = model
            .advance_with_drive(&initial.contact.candidate, boundary, controls(), cx)
            .unwrap();
        let expected = 2e5 * boundary.position_m / (1e-6 + 2e5 * c);
        close(next.contact.candidate.pressure_pa()[12], expected, 4e-12);
        assert!(next.drive.achieved_force_n > initial.drive.achieved_force_n);
        assert!(
            (next.contact.candidate.approach_m() - initial.contact.candidate.approach_m()).abs()
                < 1e-8
        );
        assert_eq!(next.contact.receipt.plastic_work_j, 0.0);
        assert_eq!(next.drive.equilibrium_evaluations, 1);
    });
}

#[test]
fn drive_failures_are_transactional_and_tensile_equilibria_are_refused() {
    with_cx(|cx, gate| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let boundary = drive(3e-6);
        let loaded = model
            .advance_with_drive(&model.initial_state(), boundary, controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let saved = before.clone();
        let mut bad_drive = boundary;
        bad_drive.stiffness_n_per_m = 0.0;
        assert!(matches!(
            model.advance_with_drive(&before, bad_drive, controls(), cx),
            Err(FiniteGapPlasticDriveError::InvalidInput { .. })
        ));
        assert!(matches!(
            model.recede_with_drive(&before, &[0.0; 25], drive(4e-6), controls(), cx),
            Err(FiniteGapPlasticDriveError::PredecessorDriveMismatch { .. })
        ));
        assert!(matches!(
            model.recede_with_drive(&before, &[0.0; 24], boundary, controls(), cx),
            Err(FiniteGapPlasticDriveError::InvalidInput { .. })
        ));
        let mut exhausted = controls();
        exhausted.maximum_equilibrium_evaluations = 0;
        assert!(matches!(
            model.advance_with_drive(&before, drive(4e-6), exhausted, cx),
            Err(FiniteGapPlasticDriveError::InvalidInput { .. })
        ));
        let mut outside = controls();
        outside.minimum_approach_m = 1e-5;
        assert!(matches!(
            model.advance_with_drive(&before, boundary, outside, cx),
            Err(FiniteGapPlasticDriveError::TargetOutsideApproachBounds { .. })
        ));
        let mut foreign_config = config();
        foreign_config.hardness_pa[12] *= 0.5;
        let foreign = FiniteGapPlasticModel::new(foreign_config, cx).unwrap();
        assert_eq!(
            foreign.advance_with_drive(&before, boundary, controls(), cx),
            Err(FiniteGapPlasticDriveError::Contact(
                FiniteGapPlasticError::ForeignState
            ))
        );
        assert_eq!(before, saved);
        gate.request();
        assert!(matches!(
            model.advance_with_drive(&before, boundary, controls(), cx),
            Err(FiniteGapPlasticDriveError::Contact(
                FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled)
            ))
        ));
        assert_eq!(before, saved);
    });
}
