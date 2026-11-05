//! Analytical finite-gap force control and maintained-load recession checks.
//! Every trial must retain the original accepted plastic/recession/work state.

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
                seed: 31,
                kernel_id: 89,
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

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn self_compliance(config: &FiniteGapPlasticConfig) -> f64 {
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

#[test]
fn one_cell_force_control_discards_yielding_bracket_trials() {
    with_cx(|cx, _| {
        let config = config();
        let c = self_compliance(&config);
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let before = model.initial_state();
        let saved = before.clone();
        assert!(
            model
                .advance(&before, controls().maximum_approach_m, cx)
                .unwrap()
                .receipt
                .plastic_work_j
                > 0.0
        );
        let result = model.advance_at_load(&before, 0.5, controls(), cx).unwrap();
        close(result.contact.candidate.approach_m(), c * 0.5 / 1e-6, 3e-12);
        close(result.contact.receipt.normal.normal_force_n, 0.5, 3e-12);
        assert_eq!(result.contact.receipt.plastic_work_j, 0.0);
        assert_eq!(result.contact.candidate.plastic_depth_m(), &[0.0; 25]);
        assert_eq!(result.contact.candidate.recession_m(), &[0.0; 25]);
        assert_eq!(result.contact.candidate.cumulative_plastic_work_j(), 0.0);
        assert_eq!(
            result.contact.receipt.predecessor_identity,
            before.identity()
        );
        assert_eq!(result.load.predecessor_state_identity, before.identity());
        assert_eq!(
            result.load.response_identity,
            result.contact.receipt.identity
        );
        assert!(result.load.equilibrium_evaluations > 2);
        assert!(result.load.force_residual_n.abs() <= 1e-13 + 1e-12 * 0.5);
        assert_eq!(
            result.load.disposition,
            FiniteGapPlasticLoadDisposition::ForceBalanced
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn force_control_retains_dent_recession_and_original_datum() {
    with_cx(|cx, _| {
        let config = config();
        let c = self_compliance(&config);
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let loaded = model
            .advance(&model.initial_state(), c * 1e6 + 2e-6, cx)
            .unwrap();
        let mut removal = vec![0.0; 25];
        removal[12] = 3e-7;
        let worn = model.recede(&loaded.candidate, &removal, cx).unwrap();
        let before = worn.candidate;
        let result = model.advance_at_load(&before, 0.4, controls(), cx).unwrap();
        let expected = before.plastic_depth_m()[12] + before.recession_m()[12] + c * 0.4 / 1e-6;
        close(result.contact.candidate.approach_m(), expected, 3e-12);
        assert_eq!(
            result.contact.candidate.plastic_depth_m(),
            before.plastic_depth_m()
        );
        assert_eq!(result.contact.candidate.recession_m(), before.recession_m());
        assert_eq!(
            result.contact.candidate.cumulative_plastic_work_j(),
            before.cumulative_plastic_work_j()
        );
        assert_eq!(result.contact.receipt.plastic_work_j, 0.0);
        assert_eq!(result.contact.receipt.normal.undeformed_gap_datum_m, 0.007);
    });
}

#[test]
fn maintained_load_recession_separates_actuator_work_from_parameter_work() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let loaded = model
            .advance_at_load(&model.initial_state(), 0.4, controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let saved = before.clone();
        let mut recession = vec![0.0; 25];
        recession[12] = 5e-7;
        let result = model
            .recede_at_load(&before, &recession, 0.4, controls(), cx)
            .unwrap();
        let receipt = &result.continuation.receipt;
        let candidate = &result.continuation.candidate;
        close(
            candidate.approach_m() - before.approach_m(),
            recession[12],
            3e-11,
        );
        close(result.load.achieved_force_n, 0.4, 3e-12);
        assert_eq!(candidate.plastic_depth_m(), before.plastic_depth_m());
        assert_eq!(
            candidate.cumulative_plastic_work_j(),
            before.cumulative_plastic_work_j()
        );
        assert_eq!(candidate.recession_m(), recession.as_slice());
        assert_eq!(
            receipt.kind,
            FiniteGapPlasticContinuationKind::SurfaceRecessionAtLoad
        );
        assert_eq!(receipt.previous_approach_m, before.approach_m());
        assert_eq!(receipt.approach_m, candidate.approach_m());
        close(receipt.normal_work_j, 0.4 * recession[12], 3e-11);
        assert!(receipt.elastic_energy_change_j.abs() <= 1e-18);
        close(receipt.parameter_work_j, -receipt.normal_work_j, 3e-11);
        assert_eq!(receipt.response.plastic_work_j, 0.0);
        assert_eq!(result.load.response_identity, receipt.response.identity);
        assert_eq!(before, saved);
    });
}

#[test]
fn two_cell_load_control_passes_an_intermediate_plateau() {
    with_cx(|cx, _| {
        let mut config = config();
        config.undeformed_gap_m[12] = 0.017;
        config.undeformed_gap_m[11] = 0.007;
        config.undeformed_gap_m[13] = 0.007 + 5e-6;
        let c = self_compliance(&config);
        let d = cross_compliance(&config);
        let relative_second_gap = config.undeformed_gap_m[13] - 0.007;
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let before = model.initial_state();
        let result = model.advance_at_load(&before, 1.5, controls(), cx).unwrap();
        let expected_approach = relative_second_gap + d * 1e6 + c * 0.5e6;
        let expected_dent = expected_approach - c * 1e6 - d * 0.5e6;
        close(result.contact.candidate.pressure_pa()[11], 1e6, 3e-12);
        close(result.contact.candidate.pressure_pa()[13], 0.5e6, 3e-12);
        close(
            result.contact.candidate.approach_m(),
            expected_approach,
            3e-12,
        );
        close(
            result.contact.candidate.plastic_depth_m()[11],
            expected_dent,
            3e-12,
        );
        assert_eq!(result.contact.candidate.plastic_depth_m()[13], 0.0);
        close(
            result.contact.receipt.plastic_work_j,
            1e6 * expected_dent * 1e-6,
            3e-12,
        );
    });
}

#[test]
fn maintained_load_two_cell_wear_redistributes_force_and_stores_elastic_energy() {
    with_cx(|cx, _| {
        let mut config = config();
        config.undeformed_gap_m[12] = 0.017;
        config.undeformed_gap_m[11] = 0.007;
        config.undeformed_gap_m[13] = 0.007;
        let c = self_compliance(&config);
        let d = cross_compliance(&config);
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let loaded = model
            .advance_at_load(&model.initial_state(), 0.8, controls(), cx)
            .unwrap();
        let before = loaded.contact.candidate;
        let w = 2e-7;
        let mut recession = vec![0.0; 25];
        recession[11] = w;
        let changed = model
            .recede_at_load(&before, &recession, 0.8, controls(), cx)
            .unwrap();
        let pressure_change = w / (2.0 * (c - d));
        close(
            changed.continuation.candidate.pressure_pa()[11],
            0.4e6 - pressure_change,
            2e-10,
        );
        close(
            changed.continuation.candidate.pressure_pa()[13],
            0.4e6 + pressure_change,
            2e-10,
        );
        close(
            changed.continuation.candidate.approach_m() - before.approach_m(),
            0.5 * w,
            2e-10,
        );
        let expected_energy_change = 1e-6 * w * w / (4.0 * (c - d));
        close(
            changed.continuation.receipt.elastic_energy_change_j,
            expected_energy_change,
            2e-10,
        );
        close(
            changed.continuation.receipt.normal_work_j,
            0.8 * 0.5 * w,
            2e-10,
        );
        close(
            changed.continuation.receipt.parameter_work_j,
            expected_energy_change - 0.8 * 0.5 * w,
            2e-10,
        );
        assert_eq!(changed.continuation.receipt.response.plastic_work_j, 0.0);
        assert_eq!(changed.continuation.candidate.plastic_depth_m(), &[0.0; 25]);
    });
}

#[test]
fn positive_plastic_plateau_capacity_and_bounds_are_distinct_refusals() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let before = model.initial_state();
        let saved = before.clone();
        assert!(matches!(
            model.advance_at_load(&before, 1.0, controls(), cx),
            Err(FiniteGapPlasticLoadError::NonUniqueLoad { .. })
        ));
        assert!(matches!(
            model.advance_at_load(&before, 1.1, controls(), cx),
            Err(FiniteGapPlasticLoadError::TargetOutsideApproachBounds { .. })
        ));
        assert!(matches!(
            model.advance_at_load(&before, 10.0, controls(), cx),
            Err(FiniteGapPlasticLoadError::TargetExceedsCapacity { .. })
        ));
        assert_eq!(before, saved);
        let saturated = model.advance(&before, 5e-6, cx).unwrap();
        let held = model
            .advance_at_load(&saturated.candidate, 1.0, controls(), cx)
            .unwrap();
        assert_eq!(held.contact.candidate, saturated.candidate);
        assert_eq!(
            held.load.disposition,
            FiniteGapPlasticLoadDisposition::RetainedExistingEquilibrium
        );
        assert_eq!(held.load.equilibrium_evaluations, 0);
        let mut removal = vec![0.0; 25];
        removal[12] = 1e-7;
        assert!(matches!(
            model.recede_at_load(&saturated.candidate, &removal, 1.0, controls(), cx),
            Err(FiniteGapPlasticLoadError::NonUniqueLoad { .. })
        ));
        assert_eq!(saturated.candidate.recession_m(), &[0.0; 25]);
    });
}

#[test]
fn zero_load_uses_the_explicit_lower_position_and_then_reuses_open_history() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let loaded = model.advance(&model.initial_state(), 5e-6, cx).unwrap();
        let open = model
            .advance_at_load(&loaded.candidate, 0.0, controls(), cx)
            .unwrap();
        assert_eq!(
            open.contact.candidate.approach_m(),
            controls().minimum_approach_m
        );
        assert_eq!(open.contact.receipt.normal.normal_force_n, 0.0);
        assert_eq!(open.contact.receipt.normal.reversible_energy_j, 0.0);
        assert_eq!(
            open.contact.candidate.plastic_depth_m(),
            loaded.candidate.plastic_depth_m()
        );
        assert_eq!(
            open.contact.candidate.cumulative_plastic_work_j(),
            loaded.candidate.cumulative_plastic_work_j()
        );
        assert_eq!(
            open.load.disposition,
            FiniteGapPlasticLoadDisposition::UnloadedAtLowerApproach
        );
        let repeat = model
            .advance_at_load(&open.contact.candidate, 0.0, controls(), cx)
            .unwrap();
        assert_eq!(repeat.contact.candidate, open.contact.candidate);
        assert_eq!(repeat.load.equilibrium_evaluations, 0);
    });
}

#[test]
fn load_refusals_and_cancellation_never_mutate_accepted_histories() {
    with_cx(|cx, gate| {
        let config = config();
        let model = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        let before = model.initial_state();
        let saved = before.clone();
        let mut bounded = controls();
        bounded.maximum_equilibrium_evaluations = 1;
        assert!(matches!(
            model.advance_at_load(&before, 0.5, bounded, cx),
            Err(FiniteGapPlasticLoadError::DidNotConverge { .. })
        ));
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                model.advance_at_load(&before, invalid, controls(), cx),
                Err(FiniteGapPlasticLoadError::InvalidInput { .. })
            ));
        }
        assert!(matches!(
            model.recede_at_load(&before, &[0.0; 25], 0.5, controls(), cx),
            Err(FiniteGapPlasticLoadError::PredecessorLoadMismatch { .. })
        ));
        let mut foreign_config = config;
        foreign_config.hardness_pa.fill(2e6);
        let foreign = FiniteGapPlasticModel::new(foreign_config, cx).unwrap();
        assert!(matches!(
            foreign.advance_at_load(&before, 0.5, controls(), cx),
            Err(FiniteGapPlasticLoadError::Contact(
                FiniteGapPlasticError::ForeignState
            ))
        ));
        assert_eq!(before, saved);
        gate.request();
        assert!(matches!(
            model.advance_at_load(&before, 0.5, controls(), cx),
            Err(FiniteGapPlasticLoadError::Contact(
                FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled)
            ))
        ));
        assert!(
            model
                .recede_at_load(&before, &[0.0; 25], 0.0, controls(), cx)
                .is_err()
        );
        assert_eq!(before, saved);
    });
}

#[test]
fn inadmissible_bracket_endpoint_keeps_the_original_domain_refusal() {
    with_cx(|cx, _| {
        let mut config = config();
        config.undeformed_gap_m[0] = 0.007 + 2e-6;
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let before = model.initial_state();
        let saved = before.clone();
        assert!(matches!(
            model.advance_at_load(&before, 0.5, controls(), cx),
            Err(FiniteGapPlasticLoadError::Contact(
                FiniteGapPlasticError::Contact(
                    FiniteGapContactError::ContactReachedDomainBoundary { .. }
                )
            ))
        ));
        assert_eq!(before, saved);
    });
}

#[test]
fn loose_force_tolerance_does_not_accept_the_far_end_of_a_plastic_plateau() {
    with_cx(|cx, _| {
        let mut config = config();
        config.undeformed_gap_m[12] = 0.017;
        config.undeformed_gap_m[11] = 0.007;
        config.undeformed_gap_m[13] = 0.007 + 5e-6;
        let c = self_compliance(&config);
        let d = cross_compliance(&config);
        let second_gap = config.undeformed_gap_m[13] - 0.007;
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let before = model.initial_state();
        let mut search = controls();
        search.maximum_approach_m = second_gap + d * 1e6 + c * 10.0;
        search.absolute_force_tolerance_n = 2e-5;
        search.relative_force_tolerance = 0.0;
        let upper = model
            .advance(&before, search.maximum_approach_m, cx)
            .unwrap();
        assert!(upper.candidate.pressure_pa()[13] > 0.0);
        assert!(
            (upper.receipt.normal.normal_force_n - 1.0).abs() < search.absolute_force_tolerance_n
        );
        assert!(matches!(
            model.advance_at_load(&before, 1.0, search, cx),
            Err(FiniteGapPlasticLoadError::NonUniqueLoad { .. })
                | Err(FiniteGapPlasticLoadError::DidNotConverge { .. })
        ));
        assert_eq!(before.plastic_depth_m(), &[0.0; 25]);
    });
}
