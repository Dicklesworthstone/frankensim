//! G1/G3 native material continuation at maintained normal force.
//! The independent one/two-cell formulas distinguish true loaded equilibria
//! from a held-approach refresh followed by unloading an artificial dent.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapContactError, FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticContinuationKind,
    FiniteGapPlasticError, FiniteGapPlasticLoadControls, FiniteGapPlasticLoadDisposition,
    FiniteGapPlasticLoadError, FiniteGapPlasticModel,
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
                kernel_id: 97,
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

fn cross_compliance(config: &FiniteGapPlasticConfig) -> f64 {
    let s = config.grid.cell_width_m;
    let b = 0.5 * s;
    let primitive = |x: f64| x * (b / x).asinh() + b * (x / b).asinh();
    2.0 * (primitive(2.5 * s) - primitive(1.5 * s))
        / (core::f64::consts::PI * config.reduced_modulus_pa)
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

#[test]
fn compliance_change_at_load_preserves_dent_recession_and_separates_work() {
    with_cx(|cx, _| {
        let original = config();
        let c = compliance(&original);
        let previous = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let indented = previous
            .advance(&previous.initial_state(), c * 1e6 + 2e-6, cx)
            .unwrap();
        let mut recession = vec![0.0; 25];
        recession[12] = 3e-7;
        let worn = previous
            .recede(&indented.candidate, &recession, cx)
            .unwrap();
        let loaded = previous
            .advance_at_load(&worn.candidate, 0.4, controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let saved = before.clone();
        let mut softer = original;
        softer.reduced_modulus_pa *= 0.5;
        let model = FiniteGapPlasticModel::new(softer, cx).unwrap();
        let next = model
            .refresh_from_at_load(&previous, &before, 0.4, controls(), cx)
            .unwrap();
        let state = &next.continuation.candidate;
        let receipt = &next.continuation.receipt;
        let pressure = 0.4 / 1e-6;
        let expected_travel = c * pressure;
        let expected_storage_change = 0.5 * 1e-6 * c * pressure * pressure;
        close(
            state.approach_m() - before.approach_m(),
            expected_travel,
            3e-11,
        );
        close(next.load.achieved_force_n, 0.4, 3e-12);
        assert_eq!(state.plastic_depth_m(), before.plastic_depth_m());
        assert_eq!(state.recession_m(), before.recession_m());
        assert_eq!(
            state.cumulative_plastic_work_j(),
            before.cumulative_plastic_work_j()
        );
        assert_eq!(receipt.response.plastic_work_j, 0.0);
        close(
            receipt.elastic_energy_change_j,
            expected_storage_change,
            3e-11,
        );
        close(receipt.normal_work_j, 0.4 * expected_travel, 3e-11);
        close(receipt.parameter_work_j, -expected_storage_change, 3e-11);
        assert_eq!(
            receipt.kind,
            FiniteGapPlasticContinuationKind::MaterialRefreshAtLoad
        );
        assert_eq!(receipt.predecessor_model_identity, previous.identity());
        assert_eq!(receipt.successor_model_identity, model.identity());
        assert_eq!(receipt.predecessor_state_identity, before.identity());
        assert_eq!(receipt.response.predecessor_identity, before.identity());
        assert_eq!(next.load.predecessor_state_identity, before.identity());
        assert_eq!(next.load.response_identity, receipt.response.identity);
        assert_eq!(receipt.response.normal.undeformed_gap_datum_m, 0.007);
        assert_eq!(receipt.recession_increment_m, vec![0.0; 25]);
        assert_eq!(before, saved);

        // Reversing an elastic material change returns the original geometry,
        // without healing either independently accumulated irreversible field.
        let restored = previous
            .refresh_from_at_load(&model, state, 0.4, controls(), cx)
            .unwrap();
        close(
            restored.continuation.candidate.approach_m(),
            before.approach_m(),
            3e-11,
        );
        assert_eq!(
            restored.continuation.candidate.plastic_depth_m(),
            before.plastic_depth_m()
        );
        assert_eq!(
            restored.continuation.candidate.recession_m(),
            before.recession_m()
        );
        close(
            restored.continuation.receipt.normal_work_j,
            -receipt.normal_work_j,
            3e-11,
        );
        close(
            restored.continuation.receipt.parameter_work_j,
            expected_storage_change,
            3e-11,
        );
    });
}

#[test]
fn direct_loaded_refresh_never_retains_a_fictitious_held_approach_dent() {
    with_cx(|cx, _| {
        let original = config();
        let previous = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let loaded = previous
            .advance_at_load(&previous.initial_state(), 0.4, controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let mut changed = original;
        changed.reduced_modulus_pa *= 4.0;
        changed.hardness_pa.fill(0.5e6);
        let c = compliance(&changed);
        let model = FiniteGapPlasticModel::new(changed, cx).unwrap();
        let fixed_position = model.refresh_from(&previous, &before, cx).unwrap();
        assert!(fixed_position.candidate.plastic_depth_m()[12] > 0.0);
        assert!(fixed_position.receipt.response.plastic_work_j > 0.0);
        let loaded_refresh = model
            .refresh_from_at_load(&previous, &before, 0.4, controls(), cx)
            .unwrap();
        close(
            loaded_refresh.continuation.candidate.approach_m(),
            c * 0.4e6,
            3e-12,
        );
        close(loaded_refresh.load.achieved_force_n, 0.4, 3e-12);
        assert_eq!(
            loaded_refresh.continuation.candidate.plastic_depth_m(),
            &[0.0; 25]
        );
        assert_eq!(
            loaded_refresh.continuation.receipt.response.plastic_work_j,
            0.0
        );
        assert_eq!(
            loaded_refresh
                .continuation
                .candidate
                .cumulative_plastic_work_j(),
            0.0
        );
        assert_eq!(before.plastic_depth_m(), &[0.0; 25]);
        assert!(matches!(
            model.advance_at_load(&before, 0.4, controls(), cx),
            Err(FiniteGapPlasticLoadError::Contact(
                FiniteGapPlasticError::ForeignState
            ))
        ));
    });
}

#[test]
fn softening_one_region_redistributes_maintained_force_and_adds_only_required_flow() {
    with_cx(|cx, _| {
        let mut original = config();
        original.undeformed_gap_m[12] = 0.017;
        original.undeformed_gap_m[11] = 0.007;
        original.undeformed_gap_m[13] = 0.007;
        let old_c = compliance(&original);
        let old_d = cross_compliance(&original);
        let previous = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let loaded = previous
            .advance_at_load(&previous.initial_state(), 1.0, controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let saved = before.clone();
        let mut changed = original;
        changed.hardness_pa[11] = 0.25e6;
        changed.reduced_modulus_pa *= 0.75;
        let c = compliance(&changed);
        let d = cross_compliance(&changed);
        let model = FiniteGapPlasticModel::new(changed, cx).unwrap();
        let result = model
            .refresh_from_at_load(&previous, &before, 1.0, controls(), cx)
            .unwrap();
        let state = &result.continuation.candidate;
        let receipt = &result.continuation.receipt;
        let p1 = 0.25e6;
        let p2 = 0.75e6;
        let expected_approach = d * p1 + c * p2;
        let expected_dent = (c - d) * (p2 - p1);
        let old_energy = 1e-6 * 0.5e6_f64.powi(2) * (old_c + old_d);
        let new_energy = 0.5 * 1e-6 * (c * (p1 * p1 + p2 * p2) + 2.0 * d * p1 * p2);
        close(state.pressure_pa()[11], p1, 3e-12);
        close(state.pressure_pa()[13], p2, 3e-12);
        close(state.approach_m(), expected_approach, 3e-12);
        close(state.plastic_depth_m()[11], expected_dent, 3e-11);
        assert_eq!(state.plastic_depth_m()[13], 0.0);
        close(
            receipt.response.plastic_work_j,
            1e-6 * p1 * expected_dent,
            3e-11,
        );
        close(receipt.previous_elastic_energy_j, old_energy, 3e-11);
        close(receipt.elastic_energy_j, new_energy, 3e-11);
        close(
            receipt.elastic_energy_change_j,
            new_energy - old_energy,
            3e-11,
        );
        close(
            receipt.normal_work_j,
            expected_approach - 0.5e6 * (old_c + old_d),
            3e-11,
        );
        close(
            receipt.parameter_work_j,
            new_energy - old_energy + 1e-6 * p1 * expected_dent - receipt.normal_work_j,
            3e-10,
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn identical_loaded_refresh_retains_the_complete_accepted_plateau_state() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let loaded = model.advance(&model.initial_state(), 5e-6, cx).unwrap();
        let before = loaded.candidate;
        let result = model
            .refresh_from_at_load(&model, &before, 1.0, controls(), cx)
            .unwrap();
        assert_eq!(result.continuation.candidate, before);
        assert_eq!(result.load.equilibrium_evaluations, 0);
        assert_eq!(
            result.load.disposition,
            FiniteGapPlasticLoadDisposition::RetainedExistingEquilibrium
        );
        assert_eq!(result.continuation.receipt.normal_work_j, 0.0);
        assert_eq!(result.continuation.receipt.elastic_energy_change_j, 0.0);
        assert_eq!(result.continuation.receipt.parameter_work_j, 0.0);
        assert_eq!(result.continuation.receipt.response.plastic_work_j, 0.0);
    });
}

#[test]
fn material_load_refresh_refuses_mismatched_history_geometry_load_and_lost_capacity() {
    with_cx(|cx, gate| {
        let original = config();
        let previous = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let loaded = previous
            .advance_at_load(&previous.initial_state(), 0.4, controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let saved = before.clone();
        let mut changed = original.clone();
        changed.reduced_modulus_pa *= 0.5;
        let model = FiniteGapPlasticModel::new(changed.clone(), cx).unwrap();
        assert!(matches!(
            model.refresh_from_at_load(&model, &before, 0.4, controls(), cx),
            Err(FiniteGapPlasticLoadError::Contact(
                FiniteGapPlasticError::ForeignState
            ))
        ));
        assert!(matches!(
            model.refresh_from_at_load(&previous, &before, 0.5, controls(), cx),
            Err(FiniteGapPlasticLoadError::PredecessorLoadMismatch { .. })
        ));
        let mut short = controls();
        short.maximum_approach_m = before.approach_m();
        assert!(matches!(
            model.refresh_from_at_load(&previous, &before, 0.4, short, cx),
            Err(FiniteGapPlasticLoadError::TargetOutsideApproachBounds { .. })
        ));
        let mut limited = controls();
        limited.maximum_equilibrium_evaluations = 1;
        assert!(matches!(
            model.refresh_from_at_load(&previous, &before, 0.4, limited, cx),
            Err(FiniteGapPlasticLoadError::DidNotConverge { .. })
        ));
        let mut wrong_gap = changed.clone();
        wrong_gap.undeformed_gap_m[11] += 1e-7;
        let wrong = FiniteGapPlasticModel::new(wrong_gap, cx).unwrap();
        assert!(matches!(
            wrong.refresh_from_at_load(&previous, &before, 0.4, controls(), cx),
            Err(FiniteGapPlasticLoadError::Contact(
                FiniteGapPlasticError::Contact(FiniteGapContactError::InvalidInput { .. })
            ))
        ));
        changed.hardness_pa.fill(1e3);
        let insufficient = FiniteGapPlasticModel::new(changed.clone(), cx).unwrap();
        assert!(matches!(
            insufficient.refresh_from_at_load(&previous, &before, 0.4, controls(), cx),
            Err(FiniteGapPlasticLoadError::TargetExceedsCapacity { .. })
        ));
        changed.hardness_pa.fill(0.4e6);
        let plateau = FiniteGapPlasticModel::new(changed, cx).unwrap();
        assert!(matches!(
            plateau.refresh_from_at_load(&previous, &before, 0.4e6 * (1e-3 * 1e-3), controls(), cx),
            Err(FiniteGapPlasticLoadError::NonUniqueLoad { .. })
        ));
        assert_eq!(before, saved);
        gate.request();
        assert!(matches!(
            model.refresh_from_at_load(&previous, &before, 0.4, controls(), cx),
            Err(FiniteGapPlasticLoadError::Contact(
                FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled)
            ))
        ));
        assert_eq!(before, saved);
    });
}
