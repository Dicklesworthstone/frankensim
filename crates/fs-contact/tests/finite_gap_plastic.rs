//! Numerical/model-limit checks for hardness-capped half-space contact.
//! These are not experimental validation of pressure saturation or J2 plasticity.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapContactError, FiniteGapContactRequest, FiniteGapGrid, FiniteGapPlasticAdvance,
    FiniteGapPlasticConfig, FiniteGapPlasticError, FiniteGapPlasticModel, FiniteGapPlasticState,
    solve_finite_gap_half_space,
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
                seed: 19,
                kernel_id: 81,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn one_cell_config() -> FiniteGapPlasticConfig {
    let mut gap = vec![0.01; 25];
    gap[12] = 0.0;
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
        maximum_active_set_iterations: 512,
        complementarity_tolerance_m: 1e-12,
        boundary_clearance_cells: 1,
    }
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn square_self_compliance(side: f64, reduced_modulus: f64) -> f64 {
    // Independent integral over [-s/2,s/2]^2 of 1/r, divided by pi E*.
    4.0 * side * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * reduced_modulus)
}

fn assert_kkt(
    config: &FiniteGapPlasticConfig,
    previous: &FiniteGapPlasticState,
    result: &FiniteGapPlasticAdvance,
    approach: f64,
) {
    let datum = config
        .undeformed_gap_m
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let tolerance = config.complementarity_tolerance_m;
    for i in 0..config.undeformed_gap_m.len() {
        let pressure = result.receipt.normal.pressure_pa[i];
        let cap = config.hardness_pa[i];
        assert!(pressure >= 0.0 && pressure <= cap);
        let gradient = result.receipt.normal.elastic_opening_m[i]
            + previous.plastic_depth_m()[i]
            + config.undeformed_gap_m[i]
            - datum
            - approach;
        if pressure == 0.0 {
            assert!(
                gradient >= -tolerance,
                "inactive cell {i}: gradient {gradient}"
            );
        } else if pressure == cap {
            assert!(gradient <= tolerance, "capped cell {i}: gradient {gradient}");
        } else {
            assert!(
                gradient.abs() <= tolerance,
                "elastic cell {i}: gradient {gradient}"
            );
        }
        let increment = result.receipt.plastic_increment_m[i];
        assert!(increment >= 0.0);
        if increment > 0.0 {
            assert_eq!(pressure, cap);
        }
        let separation = gradient + increment;
        if pressure > 0.0 {
            assert!(separation.abs() <= tolerance);
        } else {
            assert!(separation >= -tolerance);
        }
    }
    assert!(result.receipt.kkt_residual_m <= tolerance);
    assert!(result.receipt.normal.complementarity_residual_m <= tolerance);
}

#[test]
fn analytical_single_cell_cap_force_storage_and_plastic_work() {
    with_cx(|cx, _| {
        let config = one_cell_config();
        let compliance = square_self_compliance(config.grid.cell_width_m, config.reduced_modulus_pa);
        let cap = config.hardness_pa[12];
        let area = config.grid.cell_width_m * config.grid.cell_depth_m;
        let yield_approach = compliance * cap;
        let plastic_increment = 3e-6;
        let approach = yield_approach + plastic_increment;
        let model = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        let initial = model.initial_state();
        let before = initial.clone();
        let step = model.advance(&initial, approach, cx).unwrap();
        assert_eq!(initial, before);
        assert_eq!(step.receipt.normal.active_cells, 1);
        assert_eq!(step.receipt.yielded_cells, 1);
        assert_eq!(step.receipt.normal.pressure_pa[12], cap);
        close(step.receipt.normal.normal_force_n, cap * area, 1e-14, 0.0);
        close(
            step.receipt.normal.elastic_opening_m[12],
            yield_approach,
            2e-14,
            0.0,
        );
        close(
            step.candidate.plastic_depth_m()[12],
            plastic_increment,
            2e-14,
            0.0,
        );
        close(
            step.receipt.plastic_work_j,
            cap * plastic_increment * area,
            2e-14,
            0.0,
        );
        close(
            step.receipt.normal.reversible_energy_j,
            0.5 * cap * yield_approach * area,
            2e-14,
            0.0,
        );
        close(step.receipt.contact_area_m2, area, 1e-14, 0.0);
        assert_eq!(step.receipt.predecessor_identity, initial.identity());
        assert_eq!(step.candidate.owner_identity(), model.identity());
        assert_eq!(step.receipt.normal.identity, step.receipt.identity);
        assert_kkt(&config, &initial, &step, approach);
        assert_eq!(step, model.advance(&initial, approach, cx).unwrap());
    });
}

#[test]
fn unloading_and_reloading_preserve_residual_dent_and_original_gap_datum() {
    with_cx(|cx, _| {
        let mut config = one_cell_config();
        for gap in &mut config.undeformed_gap_m {
            *gap += 0.007;
        }
        let compliance = square_self_compliance(config.grid.cell_width_m, config.reduced_modulus_pa);
        let cap = config.hardness_pa[12];
        let yield_approach = compliance * cap;
        let dent = 3e-6;
        let model = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        assert_eq!(model.undeformed_gap_datum_m(), 0.007);
        let initial = model.initial_state();
        let loaded = model
            .advance(&initial, yield_approach + dent, cx)
            .unwrap();
        let half = model
            .advance(&loaded.candidate, dent + 0.5 * yield_approach, cx)
            .unwrap();
        close(half.receipt.normal.pressure_pa[12], 0.5 * cap, 2e-13, 0.0);
        assert_eq!(half.receipt.plastic_work_j, 0.0);
        assert_eq!(
            half.candidate.plastic_depth_m(),
            loaded.candidate.plastic_depth_m()
        );
        let open = model.advance(&half.candidate, -1e-6, cx).unwrap();
        assert_eq!(open.receipt.normal.normal_force_n, 0.0);
        assert_eq!(open.receipt.normal.reversible_energy_j, 0.0);
        assert_eq!(open.receipt.plastic_work_j, 0.0);
        assert_eq!(
            open.candidate.plastic_depth_m(),
            loaded.candidate.plastic_depth_m()
        );
        assert_eq!(
            open.candidate.cumulative_plastic_work_j(),
            loaded.receipt.plastic_work_j
        );
        let still_open = model.advance(&open.candidate, dent * 0.5, cx).unwrap();
        assert_eq!(
            still_open.receipt.normal.normal_force_n, 0.0,
            "reload must not renormalize the residual indentation away"
        );
        let reloaded = model
            .advance(&still_open.candidate, dent + 0.75 * yield_approach, cx)
            .unwrap();
        close(
            reloaded.receipt.normal.pressure_pa[12],
            0.75 * cap,
            2e-13,
            0.0,
        );
        assert_eq!(reloaded.receipt.plastic_work_j, 0.0);
        let deeper = model
            .advance(&reloaded.candidate, dent + yield_approach + 2e-6, cx)
            .unwrap();
        close(deeper.candidate.plastic_depth_m()[12], 5e-6, 2e-13, 0.0);
        close(
            deeper.candidate.cumulative_plastic_work_j(),
            cap * 5e-6 * 1e-6,
            2e-13,
            0.0,
        );
        assert_eq!(deeper.receipt.normal.undeformed_gap_datum_m, 0.007);
        assert_kkt(
            &config,
            &reloaded.candidate,
            &deeper,
            dent + yield_approach + 2e-6,
        );
    });
}

#[test]
fn high_hardness_recovers_existing_coupled_elastic_half_space_solution() {
    with_cx(|cx, _| {
        let grid = FiniteGapGrid {
            cells_x: 15,
            cells_y: 15,
            cell_width_m: 2e-4,
            cell_depth_m: 2e-4,
        };
        let gap: Vec<f64> = (0..225)
            .map(|i| {
                let x = (i % 15) as f64 - 7.0;
                let y = (i / 15) as f64 - 7.0;
                ((x * grid.cell_width_m).powi(2) + (y * grid.cell_depth_m).powi(2)) / (2.0 * 0.02)
            })
            .collect();
        let request = FiniteGapContactRequest {
            grid,
            undeformed_gap_m: gap.clone(),
            approach_m: 1e-5,
            reduced_modulus_pa: 2e9,
            maximum_active_set_iterations: 512,
            complementarity_tolerance_m: 1e-11,
            boundary_clearance_cells: 2,
        };
        let reference = solve_finite_gap_half_space(&request, cx).unwrap();
        let config = FiniteGapPlasticConfig {
            grid,
            undeformed_gap_m: gap,
            reduced_modulus_pa: request.reduced_modulus_pa,
            hardness_pa: vec![1e12; 225],
            maximum_active_set_iterations: 512,
            complementarity_tolerance_m: request.complementarity_tolerance_m,
            boundary_clearance_cells: request.boundary_clearance_cells,
        };
        let model = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        let initial = model.initial_state();
        let step = model.advance(&initial, request.approach_m, cx).unwrap();
        assert_eq!(step.receipt.yielded_cells, 0);
        assert_eq!(step.receipt.plastic_work_j, 0.0);
        assert_eq!(step.receipt.normal.active_cells, reference.active_cells);
        for (actual, expected) in step
            .receipt
            .normal
            .pressure_pa
            .iter()
            .zip(&reference.pressure_pa)
        {
            close(*actual, *expected, 1e-9, 1e-5);
        }
        close(
            step.receipt.normal.normal_force_n,
            reference.normal_force_n,
            1e-10,
            1e-12,
        );
        close(
            step.receipt.normal.reversible_energy_j,
            reference.reversible_energy_j,
            1e-10,
            1e-14,
        );
        assert_kkt(&config, &initial, &step, request.approach_m);
    });
}

#[test]
fn heterogeneous_caps_yield_locally_and_preserve_coupled_kkt_conditions() {
    with_cx(|cx, _| {
        let mut config = one_cell_config();
        config.undeformed_gap_m[12] = 0.01;
        config.undeformed_gap_m[11] = 0.0;
        config.undeformed_gap_m[13] = 0.0;
        config.hardness_pa[13] = 3e6;
        let model = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        let initial = model.initial_state();
        let mixed = model.advance(&initial, 2e-6, cx).unwrap();
        assert_eq!(mixed.receipt.normal.active_cells, 2);
        assert_eq!(mixed.receipt.yielded_cells, 1);
        assert_eq!(
            mixed.receipt.normal.pressure_pa[11],
            config.hardness_pa[11]
        );
        assert!(mixed.receipt.normal.pressure_pa[13] > 0.0);
        assert!(mixed.receipt.normal.pressure_pa[13] < config.hardness_pa[13]);
        assert!(mixed.candidate.plastic_depth_m()[11] > 0.0);
        assert_eq!(mixed.candidate.plastic_depth_m()[13], 0.0);
        assert_kkt(&config, &initial, &mixed, 2e-6);
        let capped = model.advance(&mixed.candidate, 1e-5, cx).unwrap();
        assert_eq!(capped.receipt.yielded_cells, 2);
        close(
            capped.receipt.normal.normal_force_n,
            (1e6 + 3e6) * 1e-6,
            1e-14,
            0.0,
        );
        let independently_summed_work: f64 = config
            .hardness_pa
            .iter()
            .zip(&capped.receipt.plastic_increment_m)
            .map(|(h, dz)| h * dz * 1e-6)
            .sum();
        close(
            capped.receipt.plastic_work_j,
            independently_summed_work,
            1e-14,
            0.0,
        );
        assert_kkt(&config, &mixed.candidate, &capped, 1e-5);
    });
}

#[test]
fn foreign_history_nonconvergence_and_boundary_failures_do_not_change_state() {
    with_cx(|cx, _| {
        let config = one_cell_config();
        let model = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        let initial = model.initial_state();
        let before = initial.clone();
        let mut other_config = config.clone();
        other_config.hardness_pa[12] *= 2.0;
        let other = FiniteGapPlasticModel::new(other_config, cx).unwrap();
        assert_eq!(
            other.advance(&initial, 1e-6, cx).unwrap_err(),
            FiniteGapPlasticError::ForeignState
        );
        assert_eq!(initial, before);
        let mut bounded = config.clone();
        bounded.maximum_active_set_iterations = 1;
        let bounded = FiniteGapPlasticModel::new(bounded, cx).unwrap();
        let state = bounded.initial_state();
        let unchanged = state.clone();
        assert!(matches!(
            bounded.advance(&state, 1e-5, cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::ActiveSetDidNotConverge { .. }
            ))
        ));
        assert_eq!(state, unchanged);
        let mut edge_config = config;
        edge_config.undeformed_gap_m[0] = 0.0;
        let edge_model = FiniteGapPlasticModel::new(edge_config, cx).unwrap();
        let state = edge_model.initial_state();
        let unchanged = state.clone();
        assert!(matches!(
            edge_model.advance(&state, 1e-5, cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::ContactReachedDomainBoundary { .. }
            ))
        ));
        assert_eq!(state, unchanged);
    });
}

#[test]
fn cancellation_is_typed_and_does_not_publish_plastic_history() {
    with_cx(|cx, gate| {
        let model = FiniteGapPlasticModel::new(one_cell_config(), cx).unwrap();
        let state = model.initial_state();
        let unchanged = state.clone();
        gate.request();
        assert_eq!(
            model.advance(&state, 1e-5, cx).unwrap_err(),
            FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled)
        );
        assert_eq!(state, unchanged);
    });
}

#[test]
fn invalid_hardness_and_nonfinite_trial_inputs_are_refused() {
    with_cx(|cx, _| {
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut config = one_cell_config();
            config.hardness_pa[12] = invalid;
            assert!(matches!(
                FiniteGapPlasticModel::new(config, cx),
                Err(FiniteGapPlasticError::Contact(
                    FiniteGapContactError::InvalidInput { .. }
                ))
            ));
        }
        let model = FiniteGapPlasticModel::new(one_cell_config(), cx).unwrap();
        let state = model.initial_state();
        let unchanged = state.clone();
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                model.advance(&state, invalid, cx),
                Err(FiniteGapPlasticError::Contact(
                    FiniteGapContactError::InvalidInput { .. }
                ))
            ));
        }
        assert_eq!(state, unchanged);
    });
}
