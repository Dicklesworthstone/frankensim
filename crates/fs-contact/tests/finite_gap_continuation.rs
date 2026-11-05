//! Independent one/two-cell checks of held-approach geometry/material changes.
//! These validate the pressure-saturation approximation, not full plasticity.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapContactError, FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticContinuationKind,
    FiniteGapPlasticError, FiniteGapPlasticModel,
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
                seed: 29,
                kernel_id: 83,
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

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn compliance(config: &FiniteGapPlasticConfig) -> f64 {
    4.0 * config.grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (core::f64::consts::PI * config.reduced_modulus_pa)
}

fn two_cells() -> FiniteGapPlasticConfig {
    let mut config = config();
    config.undeformed_gap_m[12] = 0.017;
    config.undeformed_gap_m[11] = 0.007;
    config.undeformed_gap_m[13] = 0.007;
    config.hardness_pa.fill(1e9);
    config
}

fn cross_compliance(config: &FiniteGapPlasticConfig) -> f64 {
    // Independent closed rectangle integral for cells two widths apart.
    let s = config.grid.cell_width_m;
    let b = 0.5 * s;
    let primitive = |x: f64| x * (b / x).asinh() + b * (x / b).asinh();
    2.0 * (primitive(2.5 * s) - primitive(1.5 * s))
        / (core::f64::consts::PI * config.reduced_modulus_pa)
}

#[test]
fn recession_unloads_elastic_storage_without_erasing_dent_or_gap_datum() {
    with_cx(|cx, _| {
        let config = config();
        let c = compliance(&config);
        let h = config.hardness_pa[12];
        let area = 1e-6;
        let approach = c * h + 3e-6;
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let loaded = model.advance(&model.initial_state(), approach, cx).unwrap();
        let saved = loaded.candidate.clone();
        let mut recession = vec![0.0; 25];
        recession[12] = 0.25 * c * h;
        let worn = model.recede(&loaded.candidate, &recession, cx).unwrap();
        assert_eq!(loaded.candidate, saved);
        assert_eq!(worn.candidate.plastic_depth_m(), saved.plastic_depth_m());
        assert_eq!(worn.candidate.recession_m(), recession.as_slice());
        assert_eq!(worn.candidate.approach_m(), approach);
        assert_eq!(
            worn.receipt.kind,
            FiniteGapPlasticContinuationKind::SurfaceRecession
        );
        assert_eq!(worn.receipt.response.predecessor_identity, saved.identity());
        assert_eq!(worn.receipt.predecessor_state_identity, saved.identity());
        assert_eq!(
            worn.receipt.successor_state_identity,
            worn.candidate.identity()
        );
        assert_eq!(worn.receipt.response.normal.undeformed_gap_datum_m, 0.007);
        close(worn.candidate.pressure_pa()[12], 0.75 * h, 2e-13);
        close(
            worn.receipt.previous_elastic_energy_j,
            0.5 * c * h * h * area,
            2e-13,
        );
        close(
            worn.receipt.elastic_energy_j,
            0.5 * c * (0.75 * h).powi(2) * area,
            2e-13,
        );
        close(
            worn.receipt.parameter_work_j,
            (0.75_f64.powi(2) - 1.0) * 0.5 * c * h * h * area,
            2e-13,
        );
        assert_eq!(worn.receipt.response.plastic_work_j, 0.0);
        assert_eq!(
            worn.candidate.cumulative_plastic_work_j(),
            saved.cumulative_plastic_work_j()
        );
        recession[12] = c * h;
        let open = model.recede(&worn.candidate, &recession, cx).unwrap();
        assert_eq!(open.receipt.response.normal.normal_force_n, 0.0);
        assert_eq!(open.receipt.elastic_energy_j, 0.0);
        let threshold = open.candidate.plastic_depth_m()[12] + open.candidate.recession_m()[12];
        let reloaded = model
            .advance(&open.candidate, threshold + 0.5 * c * h, cx)
            .unwrap();
        close(
            reloaded.receipt.normal.normal_force_n,
            0.5 * h * area,
            2e-13,
        );
        assert_eq!(
            reloaded.candidate.recession_m(),
            open.candidate.recession_m()
        );
        assert_eq!(
            reloaded.candidate.plastic_depth_m(),
            saved.plastic_depth_m()
        );
    });
}

#[test]
fn material_refresh_changes_elastic_stiffness_at_fixed_approach_and_recession() {
    with_cx(|cx, _| {
        let mut original = config();
        original.hardness_pa.fill(1e9);
        let c = compliance(&original);
        let model = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let loaded = model.advance(&model.initial_state(), 1e-6, cx).unwrap();
        let mut recession = vec![0.0; 25];
        recession[12] = 2e-7;
        let worn = model.recede(&loaded.candidate, &recession, cx).unwrap();
        let before = worn.candidate.clone();
        let mut changed = original.clone();
        changed.reduced_modulus_pa *= 2.0;
        let stiffer = FiniteGapPlasticModel::new(changed, cx).unwrap();
        let refresh = stiffer.refresh_from(&model, &before, cx).unwrap();
        close(
            refresh.candidate.pressure_pa()[12],
            2.0 * (1e-6 - 2e-7) / c,
            2e-13,
        );
        close(
            refresh.receipt.elastic_energy_j,
            2.0 * worn.receipt.elastic_energy_j,
            2e-13,
        );
        close(
            refresh.receipt.parameter_work_j,
            worn.receipt.elastic_energy_j,
            2e-13,
        );
        assert_eq!(refresh.receipt.response.plastic_work_j, 0.0);
        assert_eq!(refresh.candidate.recession_m(), before.recession_m());
        assert_eq!(
            refresh.candidate.plastic_depth_m(),
            before.plastic_depth_m()
        );
        assert_eq!(refresh.candidate.approach_m(), before.approach_m());
        assert_eq!(refresh.candidate.owner_identity(), stiffer.identity());
        assert_eq!(refresh.receipt.predecessor_model_identity, model.identity());
        assert_eq!(refresh.receipt.successor_model_identity, stiffer.identity());
        assert_eq!(
            refresh.receipt.response.predecessor_identity,
            before.identity()
        );
        assert!(matches!(
            stiffer.advance(&before, before.approach_m(), cx),
            Err(FiniteGapPlasticError::ForeignState)
        ));
        assert!(matches!(
            model.advance(&refresh.candidate, before.approach_m(), cx),
            Err(FiniteGapPlasticError::ForeignState)
        ));
        let softened = model
            .refresh_from(&stiffer, &refresh.candidate, cx)
            .unwrap();
        close(
            softened.candidate.pressure_pa()[12],
            before.pressure_pa()[12],
            2e-13,
        );
        close(
            softened.receipt.parameter_work_j,
            -worn.receipt.elastic_energy_j,
            2e-13,
        );
        assert_eq!(before, worn.candidate);
    });
}

#[test]
fn softer_hardness_adds_plastic_depth_without_repricing_prior_work() {
    with_cx(|cx, _| {
        let mut original = config();
        original.hardness_pa.fill(4e6);
        let c = compliance(&original);
        let model = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let loaded = model
            .advance(&model.initial_state(), 4e6 * c + 7e-7, cx)
            .unwrap();
        let mut recession = vec![0.0; 25];
        recession[12] = 1e6 * c;
        let worn = model.recede(&loaded.candidate, &recession, cx).unwrap();
        close(worn.candidate.pressure_pa()[12], 3e6, 2e-13);
        let mut updated = original;
        updated.hardness_pa.fill(1e6);
        let softer = FiniteGapPlasticModel::new(updated.clone(), cx).unwrap();
        let refresh = softer.refresh_from(&model, &worn.candidate, cx).unwrap();
        let added_depth = 2e6 * c;
        close(refresh.candidate.pressure_pa()[12], 1e6, 2e-13);
        close(
            refresh.receipt.response.plastic_increment_m[12],
            added_depth,
            2e-13,
        );
        close(
            refresh.candidate.plastic_depth_m()[12],
            worn.candidate.plastic_depth_m()[12] + added_depth,
            2e-13,
        );
        close(
            refresh.receipt.response.plastic_work_j,
            1e6 * added_depth * 1e-6,
            2e-13,
        );
        close(
            refresh.candidate.cumulative_plastic_work_j(),
            loaded.candidate.cumulative_plastic_work_j() + 1e6 * added_depth * 1e-6,
            2e-13,
        );
        close(
            refresh.receipt.parameter_work_j,
            -0.5 * c * (2e6_f64).powi(2) * 1e-6,
            2e-13,
        );
        assert_eq!(
            refresh.candidate.recession_m(),
            worn.candidate.recession_m()
        );
        updated.hardness_pa.fill(1e8);
        let harder = FiniteGapPlasticModel::new(updated, cx).unwrap();
        let hardened = harder
            .refresh_from(&softer, &refresh.candidate, cx)
            .unwrap();
        assert_eq!(
            hardened.candidate.plastic_depth_m(),
            refresh.candidate.plastic_depth_m()
        );
        assert_eq!(
            hardened.candidate.recession_m(),
            refresh.candidate.recession_m()
        );
        assert_eq!(hardened.receipt.response.plastic_work_j, 0.0);
    });
}

#[test]
fn two_cell_recession_can_yield_an_untouched_neighbor() {
    with_cx(|cx, _| {
        let mut config = two_cells();
        config.hardness_pa[13] = 1.6e6;
        let c = compliance(&config);
        let d = cross_compliance(&config);
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let loaded = model.advance(&model.initial_state(), 2e-6, cx).unwrap();
        close(loaded.candidate.pressure_pa()[11], 2e-6 / (c + d), 2e-11);
        assert_eq!(loaded.receipt.plastic_work_j, 0.0);
        let mut recession = vec![0.0; 25];
        recession[11] = 8e-7;
        let changed = model.recede(&loaded.candidate, &recession, cx).unwrap();
        let expected_p1 = (1.2e-6 - d * 1.6e6) / c;
        let expected_u2 = d * expected_p1 + c * 1.6e6;
        let expected_plastic = 2e-6 - expected_u2;
        assert!(expected_plastic > 0.0);
        close(changed.candidate.pressure_pa()[11], expected_p1, 2e-11);
        close(changed.candidate.pressure_pa()[13], 1.6e6, 2e-13);
        close(
            changed.candidate.plastic_depth_m()[13],
            expected_plastic,
            2e-10,
        );
        assert_eq!(changed.candidate.plastic_depth_m()[11], 0.0);
        assert_eq!(changed.candidate.recession_m()[13], 0.0);
        close(
            changed.receipt.response.plastic_work_j,
            1.6e6 * expected_plastic * 1e-6,
            2e-10,
        );
        let expected_energy = 0.5 * 1e-6 * (expected_p1 * 1.2e-6 + 1.6e6 * expected_u2);
        close(changed.receipt.elastic_energy_j, expected_energy, 2e-11);
        close(
            changed.receipt.parameter_work_j,
            expected_energy - loaded.receipt.normal.reversible_energy_j
                + 1.6e6 * expected_plastic * 1e-6,
            2e-10,
        );
        assert!(changed.receipt.parameter_work_j < 0.0);
    });
}

#[test]
fn two_cell_material_refresh_matches_coupled_capped_solution() {
    with_cx(|cx, _| {
        let original = two_cells();
        let c = compliance(&original);
        let d = cross_compliance(&original);
        let previous = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let loaded = previous
            .advance(&previous.initial_state(), 2e-6, cx)
            .unwrap();
        let mut updated = original;
        updated.hardness_pa[11] = 1e6;
        updated.hardness_pa[13] = 3e6;
        let model = FiniteGapPlasticModel::new(updated, cx).unwrap();
        let result = model
            .refresh_from(&previous, &loaded.candidate, cx)
            .unwrap();
        let expected_p2 = (2e-6 - d * 1e6) / c;
        let expected_plastic = 2e-6 - c * 1e6 - d * expected_p2;
        close(result.candidate.pressure_pa()[11], 1e6, 2e-13);
        close(result.candidate.pressure_pa()[13], expected_p2, 2e-11);
        close(
            result.receipt.response.plastic_increment_m[11],
            expected_plastic,
            2e-11,
        );
        assert_eq!(result.receipt.response.plastic_increment_m[13], 0.0);
        assert_eq!(
            result.receipt.predecessor_state_identity,
            loaded.candidate.identity()
        );
        assert_eq!(result.candidate.recession_m(), &[0.0; 25]);
    });
}

#[test]
fn unchanged_inputs_preserve_the_complete_physical_state_exactly() {
    with_cx(|cx, _| {
        let config = config();
        let model = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        let loaded = model.advance(&model.initial_state(), 5e-6, cx).unwrap();
        let mut recession = vec![0.0; 25];
        recession[12] = 1e-7;
        let before = model
            .recede(&loaded.candidate, &recession, cx)
            .unwrap()
            .candidate;
        let identical = FiniteGapPlasticModel::new(config, cx).unwrap();
        let refreshed = identical.refresh_from(&model, &before, cx).unwrap();
        assert_eq!(refreshed.candidate, before);
        assert_eq!(refreshed.receipt.parameter_work_j, 0.0);
        assert_eq!(refreshed.receipt.elastic_energy_change_j, 0.0);
        assert_eq!(refreshed.receipt.response.plastic_work_j, 0.0);
        let unworn = model.recede(&before, &[0.0; 25], cx).unwrap();
        assert_eq!(unworn.candidate, before);
        assert_eq!(unworn.receipt.parameter_work_j, 0.0);
        assert_ne!(unworn.receipt.identity, refreshed.receipt.identity);
    });
}

#[test]
fn small_recession_work_retains_the_one_cell_energy_difference() {
    with_cx(|cx, _| {
        let mut config = config();
        config.hardness_pa.fill(1e9);
        // Resolve this parameter change instead of accepting the previous
        // pressure within a looser displacement tolerance.
        config.complementarity_tolerance_m = 1e-20;
        let c = compliance(&config);
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let approach = 1e-5;
        let increment = 1e-15;
        let loaded = model.advance(&model.initial_state(), approach, cx).unwrap();
        let mut recession = vec![0.0; 25];
        recession[12] = increment;
        let result = model.recede(&loaded.candidate, &recession, cx).unwrap();
        let expected = (-approach * increment + 0.5 * increment * increment) * 1e-6 / c;
        close(result.receipt.parameter_work_j, expected, 2e-5);
        assert!(result.receipt.parameter_work_j < 0.0);
        assert_eq!(result.receipt.response.plastic_work_j, 0.0);
    });
}

#[test]
fn continuation_refusals_preserve_accepted_history() {
    with_cx(|cx, gate| {
        let original = config();
        let model = FiniteGapPlasticModel::new(original.clone(), cx).unwrap();
        let before = model
            .advance(&model.initial_state(), 4e-6, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            let mut increments = vec![0.0; 25];
            increments[12] = invalid;
            assert!(model.recede(&before, &increments, cx).is_err());
        }
        assert!(model.recede(&before, &[], cx).is_err());
        let mut changed = original.clone();
        changed.hardness_pa.fill(2e6);
        let other = FiniteGapPlasticModel::new(changed, cx).unwrap();
        assert!(matches!(
            other.recede(&before, &[0.0; 25], cx),
            Err(FiniteGapPlasticError::ForeignState)
        ));
        assert!(matches!(
            model.refresh_from(&other, &before, cx),
            Err(FiniteGapPlasticError::ForeignState)
        ));
        for which in 0..3 {
            let mut changed = original.clone();
            match which {
                0 => {
                    for gap in &mut changed.undeformed_gap_m {
                        *gap += 0.001;
                    }
                }
                1 => changed.grid.cell_width_m *= 2.0,
                _ => changed.maximum_active_set_iterations += 1,
            }
            let target = FiniteGapPlasticModel::new(changed, cx).unwrap();
            assert!(target.refresh_from(&model, &before, cx).is_err());
        }
        assert_eq!(before, saved);
        gate.request();
        assert!(matches!(
            model.recede(&before, &[0.0; 25], cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::Cancelled
            ))
        ));
        assert!(matches!(
            other.refresh_from(&model, &before, cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::Cancelled
            ))
        ));
        assert_eq!(before, saved);
    });
}

#[test]
fn recession_that_moves_contact_into_the_boundary_is_refused_atomically() {
    with_cx(|cx, _| {
        let mut config = config();
        config.hardness_pa.fill(1e9);
        let reference = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        let approach = 1e-6;
        let reference_load = reference
            .advance(&reference.initial_state(), approach, cx)
            .unwrap();
        let edge_opening = reference_load.receipt.normal.elastic_opening_m[0];
        config.undeformed_gap_m[0] = 0.007 + approach - 0.5 * edge_opening;
        let model = FiniteGapPlasticModel::new(config, cx).unwrap();
        let loaded = model.advance(&model.initial_state(), approach, cx).unwrap();
        assert_eq!(loaded.candidate.pressure_pa()[0], 0.0);
        let saved = loaded.candidate.clone();
        let mut recession = vec![0.0; 25];
        recession[12] = 2.0 * approach;
        assert!(matches!(
            model.recede(&loaded.candidate, &recession, cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::ContactReachedDomainBoundary { .. }
            ))
        ));
        assert_eq!(loaded.candidate, saved);
    });
}
