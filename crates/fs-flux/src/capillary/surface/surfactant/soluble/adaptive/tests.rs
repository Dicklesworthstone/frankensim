use super::*;
use fs_material::osmosis::nonideal::VirialSolutionModel;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::sync::Arc;

const TEMPERATURE: f64 = 300.0;
const RT: f64 = R * TEMPERATURE;
const VOLUME: f64 = 1e-6;

fn law(saturation: f64) -> SurfactantLaw {
    SurfactantLaw::Langmuir {
        clean_tension_n_m: 0.2,
        thermal_energy_j_mol: RT,
        saturation_concentration_mol_m2: saturation,
    }
}

fn bulk(reference: f64, amounts: Vec<f64>, interactions: Vec<f64>) -> AqueousSolution {
    let model = VirialSolutionModel::new(
        TEMPERATURE,
        18e-6,
        reference,
        vec![20.0; amounts.len()],
        interactions,
    )
    .unwrap();
    AqueousSolution::nonideal(Arc::new(model), VOLUME, amounts).unwrap()
}

fn physical(step: f64) -> SolubleSurfaceControls {
    SolubleSurfaceControls {
        maximum_step_s: step,
        maximum_substeps: 100_000,
        maximum_root_evaluations: 2_000_000,
        maximum_surface_pair_exchanges: 2_000_000,
        maximum_faces: 3,
        maximum_edges: 3,
        root_tolerance_mol: 1e-20,
        amount_tolerance_mol: 1e-15,
        energy_tolerance_j: 1e-11,
    }
}

fn accuracy(step: f64, absolute: f64) -> AdaptiveSolubleSurfaceControls {
    AdaptiveSolubleSurfaceControls {
        initial_step_s: step,
        minimum_step_s: step * 1e-9,
        absolute_surface_amount_tolerance_mol: absolute,
        absolute_bulk_amount_tolerance_mol: absolute,
        absolute_bulk_activity_tolerance: absolute / VOLUME,
        relative_tolerance: 0.0,
        maximum_trials: 10_000,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance {tolerance:e}"
    );
}

fn edge(i: usize, j: usize, conductance: f64) -> SurfacePartitionDiffusionEdge {
    SurfacePartitionDiffusionEdge {
        faces: [i, j],
        conductance_mol_s: conductance,
    }
}

fn riccati(mobile: f64, occupied: f64, capacity: f64, kd: f64, duration: f64) -> f64 {
    // n'=a*(M-n)*(C-n)-kd*n. Derive the exact scalar trajectory from
    // its two roots; no original numerical transport map enters this oracle.
    let total = mobile + occupied;
    let a = kd / VOLUME;
    let b = kd + a * (total + capacity);
    let c = a * total * capacity;
    let gap = (b * b - 4.0 * a * c).sqrt();
    let low = 2.0 * c / (b + gap);
    let high = (b + gap) / (2.0 * a);
    let ratio = (occupied - low) / (occupied - high) * (-gap * duration).exp();
    (low - ratio * high) / (1.0 - ratio)
}

#[test]
fn finite_depletion_and_desorption_converge_to_the_independent_riccati_flow() {
    let (area, saturation, kd, duration) = (0.2, 1e-5, 2.0, 0.7);
    for (mobile, occupied) in [(3e-6, 0.0), (0.0, 1.5e-6)] {
        let exact = riccati(mobile, occupied, area * saturation, kd, duration);
        let mut coarse_error = None;
        for tolerance in [2e-10, 2e-12] {
            let mut bath = bulk(1.0, vec![mobile, 1e-6], vec![0.0; 4]);
            let mut amounts = [occupied];
            let report = advance_soluble_surface_adaptive(
                &mut bath,
                &[law(saturation)],
                &[area],
                &mut amounts,
                SolubleSurfaceKinetics {
                    bulk_species: 0,
                    standard_potential_offset_j_mol: 0.0,
                    desorption_rates_s_inverse: &[kd],
                },
                &[],
                duration,
                physical(duration),
                accuracy(duration, tolerance),
                || false,
            )
            .unwrap();
            let error = (amounts[0] - exact).abs();
            if let Some(previous) = coarse_error {
                assert!(error * 5.0 < previous, "{error:e} vs {previous:e}");
                assert!(error < 3e-9, "{error:e}");
            }
            coarse_error = Some(error);
            assert!(report.rejected_intervals > 0);
            assert!(report.maximum_accepted_normalized_error <= 1.0);
            assert_eq!(
                report.trials,
                report.accepted_intervals + report.rejected_intervals
            );
            assert_eq!(report.evaluated_substeps, 3 * report.trials);
            assert_eq!(report.transport.substeps, 2 * report.accepted_intervals);
            assert!(report.evaluated_root_evaluations > report.transport.root_evaluations);
            assert_eq!(report.evaluated_surface_pair_exchanges, 0);
            assert_eq!(bath.solute_moles()[1], 1e-6);
            near(
                bath.solute_moles()[0] + amounts[0],
                mobile + occupied,
                1e-15,
            );
            near(report.transport.amount_defect_mol, 0.0, 1e-15);
            near(report.transport.energy_account_residual_j, 0.0, 1e-11);
            assert!(amounts[0] > 0.0 && amounts[0] < area * saturation);
            assert!(bath.solute_moles()[0] > 0.0);
            assert!(report.transport.free_energy_change_j < 0.0);
            assert!(report.transport.exchange_free_energy_dissipation_j > 0.0);
            near(
                report.transport.gross_adsorption_mol - report.transport.gross_desorption_mol,
                amounts[0] - occupied,
                1e-15,
            );
        }
    }
}

fn nonlinear_reference(initial: [f64; 4], steps: usize) -> [f64; 4] {
    let capacity = [0.8e-6, 1.2e-6, 2e-6];
    let kd = [0.5, 3.0, 1.5];
    let edges = [(0, 1, 8e-6), (1, 2, 0.4e-6), (0, 2, 1e-6)];
    let rate = |state: [f64; 4]| {
        // Authored virial row [0.2,0.05], co-solute concentration 1,
        // activity reference 1, and offset RT*ln(2).
        let concentration = state[0] / VOLUME;
        let activity = concentration * (0.2 * concentration + 0.05).exp();
        let mut derivative = [0.0; 4];
        for i in 0..3 {
            let exchange = kd[i] * (0.5 * activity * (capacity[i] - state[i + 1]) - state[i + 1]);
            derivative[i + 1] += exchange;
            derivative[0] -= exchange;
        }
        for (i, j, conductance) in edges {
            let flux = conductance * (state[j + 1] / capacity[j] - state[i + 1] / capacity[i]);
            derivative[i + 1] += flux;
            derivative[j + 1] -= flux;
        }
        derivative
    };
    let h = 0.4 / steps as f64;
    let mut state = initial;
    for _ in 0..steps {
        let a = rate(state);
        let b = rate(core::array::from_fn(|i| state[i] + 0.5 * h * a[i]));
        let c = rate(core::array::from_fn(|i| state[i] + 0.5 * h * b[i]));
        let d = rate(core::array::from_fn(|i| state[i] + h * c[i]));
        state =
            core::array::from_fn(|i| state[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    state
}

#[test]
fn nonideal_finite_bath_and_noncommuting_heterogeneous_graph_match_independent_ode() {
    let areas = [0.1, 0.2, 0.15];
    let capacities = [0.8e-6, 1.2e-6, 2e-6];
    let laws = core::array::from_fn::<_, 3, _>(|i| law(capacities[i] / areas[i]));
    let edges = [edge(0, 1, 8e-6), edge(1, 2, 0.4e-6), edge(0, 2, 1e-6)];
    let initial = [3e-6, 0.6e-6, 0.08e-6, 0.2e-6];
    let oracle = nonlinear_reference(initial, 4096);
    let finer_oracle = nonlinear_reference(initial, 8192);
    for i in 0..4 {
        near(oracle[i], finer_oracle[i], 2e-17);
    }
    let mut previous = None;
    for tolerance in [2e-9, 2e-12] {
        let mut bath = bulk(1.0, vec![initial[0], 1e-6], vec![0.2, 0.05, 0.05, 0.1]);
        let original = bath.clone();
        let mut amounts = [initial[1], initial[2], initial[3]];
        let report = advance_soluble_surface_adaptive(
            &mut bath,
            &laws,
            &areas,
            &mut amounts,
            SolubleSurfaceKinetics {
                bulk_species: 0,
                standard_potential_offset_j_mol: RT * 2.0_f64.ln(),
                desorption_rates_s_inverse: &[0.5, 3.0, 1.5],
            },
            &edges,
            0.4,
            physical(0.4),
            accuracy(0.4, tolerance),
            || false,
        )
        .unwrap();
        let mut error = (bath.solute_moles()[0] - oracle[0]).abs();
        for i in 0..3 {
            error = error.max((amounts[i] - oracle[i + 1]).abs());
            assert!(amounts[i] > 0.0 && amounts[i] < capacities[i]);
        }
        if let Some(coarse) = previous {
            assert!(error * 10.0 < coarse, "{error:e} vs {coarse:e}");
            assert!(error < 3e-9, "{error:e}");
        }
        previous = Some(error);
        assert_eq!(bath.nonideal_model(), original.nonideal_model());
        assert_eq!(bath.solute_moles()[1], 1e-6);
        assert_eq!(report.evaluated_surface_pair_exchanges, 30 * report.trials);
        assert_eq!(
            report.transport.surface_pair_exchanges,
            20 * report.accepted_intervals
        );
        assert!(report.transport.exchange_free_energy_dissipation_j > 0.0);
        assert!(report.transport.surface_diffusion_free_energy_dissipation_j > 0.0);
        near(report.transport.amount_defect_mol, 0.0, 1e-15);
        near(report.transport.energy_account_residual_j, 0.0, 1e-11);
    }
}

#[test]
fn first_order_error_factor_and_unextrapolated_fine_account_are_observable() {
    let duration = 0.3;
    let original = bulk(1.0, vec![3e-6], vec![0.0]);
    let initial = [0.1e-6];
    let laws = [law(1e-5)];
    let areas = [0.2];
    let kinetics = SolubleSurfaceKinetics {
        bulk_species: 0,
        standard_potential_offset_j_mol: 0.0,
        desorption_rates_s_inverse: &[2.0],
    };
    let mut coarse_bulk = original.clone();
    let mut coarse = initial;
    advance_soluble_surface(
        &mut coarse_bulk,
        &laws,
        &areas,
        &mut coarse,
        kinetics,
        &[],
        duration,
        physical(duration),
        || false,
    )
    .unwrap();
    let mut fine_bulk = original.clone();
    let mut fine = initial;
    let first = advance_soluble_surface(
        &mut fine_bulk,
        &laws,
        &areas,
        &mut fine,
        kinetics,
        &[],
        0.5 * duration,
        physical(duration),
        || false,
    )
    .unwrap();
    let second = advance_soluble_surface(
        &mut fine_bulk,
        &laws,
        &areas,
        &mut fine,
        kinetics,
        &[],
        0.5 * duration,
        physical(duration),
        || false,
    )
    .unwrap();
    let discrepancy = (fine[0] - coarse[0]).abs();
    assert!(discrepancy > 0.0);
    let mut limits = accuracy(duration, 0.5 * discrepancy);
    limits.absolute_bulk_amount_tolerance_mol = 1.0;
    limits.absolute_bulk_activity_tolerance = 1.0;
    limits.minimum_step_s = duration;
    let mut bath = original.clone();
    let mut amounts = initial;
    match advance_soluble_surface_adaptive(
        &mut bath,
        &laws,
        &areas,
        &mut amounts,
        kinetics,
        &[],
        duration,
        physical(duration),
        limits,
        || false,
    ) {
        Err(AdaptiveSolubleSurfaceError::TemporalAccuracy {
            normalized_error, ..
        }) => near(normalized_error, 2.0, 2e-12),
        other => panic!("first-order fine error must refuse this floor: {other:?}"),
    }
    assert_eq!(bath, original);
    assert_eq!(amounts, initial);
    limits.absolute_surface_amount_tolerance_mol = 2.0 * discrepancy;
    let report = advance_soluble_surface_adaptive(
        &mut bath,
        &laws,
        &areas,
        &mut amounts,
        kinetics,
        &[],
        duration,
        physical(duration),
        limits,
        || false,
    )
    .unwrap();
    assert_eq!(report.trials, 1);
    assert_eq!(bath, fine_bulk);
    assert_eq!(amounts, fine);
    near(
        report.transport.gross_adsorption_mol,
        first.gross_adsorption_mol + second.gross_adsorption_mol,
        1e-20,
    );
    near(
        report.transport.gross_desorption_mol,
        first.gross_desorption_mol + second.gross_desorption_mol,
        1e-20,
    );
    near(
        report.transport.free_energy_change_j,
        first.free_energy_change_j + second.free_energy_change_j,
        1e-14,
    );
    assert_eq!(
        report.transport.root_evaluations,
        first.root_evaluations + second.root_evaluations
    );
    assert!(report.evaluated_root_evaluations > report.transport.root_evaluations);
}

#[test]
fn partition_equilibrium_retains_exact_turnover_and_zero_duration_admission() {
    let capacities = [1e-6, 2e-6, 3e-6];
    let areas = [0.1, 0.2, 0.3];
    let laws = [law(1e-5); 3];
    let initial = capacities.map(|capacity| 0.5 * capacity);
    let original = bulk(1.0, vec![1e-6], vec![0.0]);
    let edges = [edge(0, 1, 1e-3), edge(1, 2, 2e-3)];
    let rates = [1.0, 0.0, 3.0];
    let kinetics = SolubleSurfaceKinetics {
        bulk_species: 0,
        standard_potential_offset_j_mol: 0.0,
        desorption_rates_s_inverse: &rates,
    };
    for duration in [0.0, 1.0] {
        let mut bath = original.clone();
        let mut amounts = initial;
        let mut work = physical(1.0);
        work.maximum_root_evaluations = 3;
        let report = advance_soluble_surface_adaptive(
            &mut bath,
            &laws,
            &areas,
            &mut amounts,
            kinetics,
            &edges,
            duration,
            work,
            accuracy(1.0, 1e-14),
            || false,
        )
        .unwrap();
        near(report.transport.free_energy_change_j, 0.0, 1e-17);
        near(report.transport.amount_defect_mol, 0.0, 1e-20);
        near(bath.solute_moles()[0], original.solute_moles()[0], 1e-20);
        for i in 0..3 {
            near(amounts[i], initial[i], 1e-20);
        }
        let turnover = duration * (rates[0] * initial[0] + rates[2] * initial[2]);
        near(report.transport.gross_adsorption_mol, turnover, 1e-20);
        near(report.transport.gross_desorption_mol, turnover, 1e-20);
        assert_eq!(
            report.evaluated_root_evaluations,
            if duration == 0.0 { 0 } else { 3 }
        );
        assert_eq!(
            report.transport.root_evaluations,
            if duration == 0.0 { 0 } else { 2 }
        );
    }
    let mut bath = original.clone();
    let mut amounts = initial;
    assert!(
        advance_soluble_surface_adaptive(
            &mut bath,
            &laws,
            &areas,
            &mut amounts,
            kinetics,
            &[edge(0, 3, 1e-6)],
            0.0,
            physical(1.0),
            accuracy(1.0, 1e-14),
            || false
        )
        .is_err()
    );
    assert_eq!(bath, original);
    assert_eq!(amounts, initial);
}

#[test]
fn accepted_prefix_rejected_probe_and_last_cancellation_never_publish_partial_state() {
    let original = bulk(1.0, vec![1e-6], vec![0.0]);
    let initial = [0.5e-6, 1e-6];
    let laws = [law(1e-5); 2];
    let areas = [0.1, 0.2];
    let edges = [edge(0, 1, 1e-5)];
    let kinetics = SolubleSurfaceKinetics {
        bulk_species: 0,
        standard_potential_offset_j_mol: 0.0,
        desorption_rates_s_inverse: &[1.0, 2.0],
    };
    for budget in 0..4 {
        let mut work = physical(0.25);
        let mut time = accuracy(0.25, 1e-14);
        match budget {
            0 => time.maximum_trials = 1,
            1 => work.maximum_substeps = 4,
            2 => work.maximum_surface_pair_exchanges = 8,
            _ => work.maximum_root_evaluations = 5,
        }
        let mut bath = original.clone();
        let mut amounts = initial;
        let error = advance_soluble_surface_adaptive(
            &mut bath,
            &laws,
            &areas,
            &mut amounts,
            kinetics,
            &edges,
            1.0,
            work,
            time,
            || false,
        )
        .unwrap_err();
        match (budget, error) {
            (
                0,
                AdaptiveSolubleSurfaceError::TrialBudget {
                    completed_trials: 1,
                },
            )
            | (
                1,
                AdaptiveSolubleSurfaceError::SubstepBudget {
                    completed: 3,
                    required: 3,
                },
            )
            | (
                2,
                AdaptiveSolubleSurfaceError::PairExchangeBudget {
                    completed: 6,
                    required: 6,
                },
            )
            | (
                3,
                AdaptiveSolubleSurfaceError::Transport(SolubleSurfaceError::RootBudget {
                    completed_evaluations: 5,
                }),
            ) => {}
            (_, other) => panic!("unexpected global budget refusal: {other:?}"),
        }
        assert_eq!(bath, original);
        assert_eq!(amounts, initial);
    }
    let mut calls = 0;
    let mut bath = original.clone();
    let mut amounts = initial;
    let report = advance_soluble_surface_adaptive(
        &mut bath,
        &laws,
        &areas,
        &mut amounts,
        kinetics,
        &edges,
        1.0,
        physical(0.25),
        accuracy(0.25, 1e-14),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(report.accepted_intervals, 4);
    for stop in [calls / 2, calls] {
        let mut bath = original.clone();
        let mut amounts = initial;
        let mut visits = 0;
        let error = advance_soluble_surface_adaptive(
            &mut bath,
            &laws,
            &areas,
            &mut amounts,
            kinetics,
            &edges,
            1.0,
            physical(0.25),
            accuracy(0.25, 1e-14),
            || {
                visits += 1;
                visits == stop
            },
        )
        .unwrap_err();
        assert!(matches!(
            error,
            AdaptiveSolubleSurfaceError::Transport(
                SolubleSurfaceError::Cancelled
                    | SolubleSurfaceError::Diffusion(
                        super::super::SurfaceDiffusionError::Cancelled
                    )
            )
        ));
        assert_eq!(bath, original);
        assert_eq!(amounts, initial);
    }
    // A rejected full comparison also consumes its three original substeps.
    let original = bulk(1.0, vec![3e-6], vec![0.0]);
    let initial = [0.0];
    let mut bath = original.clone();
    let mut amounts = initial;
    let mut work = physical(1.0);
    work.maximum_substeps = 3;
    let error = advance_soluble_surface_adaptive(
        &mut bath,
        &[law(1e-5)],
        &[0.2],
        &mut amounts,
        SolubleSurfaceKinetics {
            bulk_species: 0,
            standard_potential_offset_j_mol: 0.0,
            desorption_rates_s_inverse: &[2.0],
        },
        &[],
        1.0,
        work,
        accuracy(1.0, 1e-14),
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        AdaptiveSolubleSurfaceError::SubstepBudget {
            completed: 3,
            required: 3
        }
    ));
    assert_eq!(bath, original);
    assert_eq!(amounts, initial);
}

#[test]
fn original_activity_gate_resolves_nonideal_changes_and_extreme_reference_datums() {
    let mut reference_amount = None;
    for reference in [1.0, 1e-310] {
        let mut bath = bulk(reference, vec![0.0, 1e-6], vec![0.2, 0.05, 0.05, 0.1]);
        let mut amounts = [1.5e-6];
        let mut limits = accuracy(0.5, 1.0);
        // Amount allowances are deliberately loose: only the actual source
        // activity gate resolves desorption from the exactly vacant bath.
        limits.absolute_bulk_activity_tolerance = 1e-4 / reference;
        let report = advance_soluble_surface_adaptive(
            &mut bath,
            &[law(1e-5)],
            &[0.2],
            &mut amounts,
            SolubleSurfaceKinetics {
                bulk_species: 0,
                standard_potential_offset_j_mol: -RT * reference.ln(),
                desorption_rates_s_inverse: &[2.0],
            },
            &[],
            0.5,
            physical(0.5),
            limits,
            || false,
        )
        .unwrap();
        assert!(report.rejected_intervals > 0);
        assert!(report.accepted_intervals > 2);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert!(bath.solute_moles()[0] > 0.0);
        if let Some(expected) = reference_amount {
            near(amounts[0], expected, 3e-15);
        }
        reference_amount = Some(amounts[0]);
        assert_eq!(bath.solute_moles()[1], 1e-6);
        near(bath.solute_moles()[0] + amounts[0], 1.5e-6, 1e-15);
    }
}
