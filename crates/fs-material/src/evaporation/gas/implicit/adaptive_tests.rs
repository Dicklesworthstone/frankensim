//! Time-error admission and whole-request accounting on both original owners.
use super::*;
use crate::evaporation::gas::implicit::adaptive::{
    AdaptiveImplicitGasControls, AdaptiveImplicitGasError, advance_deposition_with_gas_adaptive,
};

fn policy(temperature_tolerance: f64, mass_tolerance: f64) -> AdaptiveImplicitGasControls {
    let mut solve = controls();
    solve.thermal.maximum_iterations = 10_000;
    solve.thermal.maximum_evaluations = 100_000;
    solve.thermal.newton.absolute_tolerance = 1e-10;
    let mut controls = AdaptiveImplicitGasControls::new(solve);
    controls.initial_dt_s = 2.0;
    controls.minimum_dt_s = 1e-5;
    controls.maximum_dt_s = 2.0;
    controls.maximum_comparisons = 2048;
    controls.maximum_reductions = 16;
    controls.maximum_gas_iterations = 100_000;
    controls.enthalpy_absolute_tolerance_j_kg = R * 4.3 * temperature_tolerance;
    controls.temperature_absolute_tolerance_k = temperature_tolerance;
    controls.mass_absolute_tolerance_kg = mass_tolerance;
    controls.relative_tolerance = 0.0;
    controls
}

fn thermal_with_budget(
    source_w: f64,
    max_substeps: usize,
    later: Vec<TimedEvent>,
) -> DepositionSimulation<'static> {
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(cell(1, 300.0, source_w, Continuous::AREA)),
    }];
    events.extend(later);
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &ORIGINAL_LIQUID)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-6,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn first_comparison_work(
    thermal: &DepositionSimulation<'_>,
    surfaces: &[ImplicitGasSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    controls: ImplicitGasControls,
) -> (usize, usize) {
    let coarse = solve_deposition_with_gas(
        thermal,
        surfaces,
        gas,
        catalog,
        2.0,
        controls,
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let first = solve_deposition_with_gas(
        thermal,
        surfaces,
        gas,
        catalog,
        1.0,
        controls,
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    let second = solve_deposition_with_gas(
        first.thermal(),
        surfaces,
        first.gas(),
        catalog,
        2.0,
        controls,
        empty_heat,
        || Ok(()),
    )
    .unwrap();
    (
        coarse.report().residual_evaluations
            + first.report().residual_evaluations
            + second.report().residual_evaluations,
        coarse.report().gas_iterations
            + first.report().gas_iterations
            + second.report().gas_iterations,
    )
}

#[test]
fn tighter_adaptive_accuracy_improves_all_physical_coordinates_and_retains_only_fine_exchange() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = thermal_with_budget(2.0, 10_000, vec![]);
    let initial_gas = gas(&source, 340.0, 5000.0, 0.001);
    let snapshot = original.clone();
    let gas_snapshot = initial_gas.clone();
    let surfaces = [surface(
        1,
        &binding,
        Continuous::MASS_COEFFICIENT,
        Continuous::HEAT_COEFFICIENT,
    )];
    let oracle = Continuous {
        carrier_mass: initial_gas.species_mass_kg(&carrier()).unwrap(),
        initial_vapor_mass: initial_gas.species_mass_kg(&vapor()).unwrap(),
        volume: initial_gas.volume_m3(),
        source_w: 2.0,
    };
    let reference = oracle.integrate([300.0, oracle.initial_vapor_mass, 340.0], 2048);
    let mut outcomes = Vec::new();
    for config in [policy(0.05, 1e-8), policy(0.0005, 1e-10)] {
        let candidate = advance_deposition_with_gas_adaptive(
            &original,
            &surfaces,
            &initial_gas,
            &source,
            Continuous::DURATION,
            config,
            |thermal, liquids, gas| {
                Ok(coefficients((
                    thermal.time_s(),
                    gas.state().state().temperature,
                    gas.species_mass_kg(&vapor())?,
                    liquids[&CellId(1)].mass_kg,
                )))
            },
            || Ok(()),
        )
        .unwrap();
        let report = candidate.report();
        let thermal = candidate.thermal();
        let gas = candidate.gas();
        let state = [
            thermal.temperature_k(CellId(1)).unwrap(),
            gas.species_mass_kg(&vapor()).unwrap(),
            gas.state().state().temperature,
        ];
        let error: [f64; 3] = std::array::from_fn(|i| (state[i] - reference[i]).abs());
        let totals = &report.surface_totals[&CellId(1)];
        let mass = original.cells()[&CellId(1)].mass_kg - thermal.cells()[&CellId(1)].mass_kg;
        assert!(mass < 0.0);
        near(totals.transferred_mass_kg, mass, 1e-13);
        near(
            totals.transferred_mass_kg,
            state[1] - oracle.initial_vapor_mass,
            1e-13,
        );
        assert_eq!(totals.evaporated_mass_kg, 0.0);
        near(totals.condensed_mass_kg, -mass, 1e-13);
        near(totals.latent_energy_j, mass * LATENT, 1e-8);
        assert!(totals.vapor_enthalpy_correction_j < 0.0);
        assert!(totals.sensible_heat_into_liquid_j > 0.0);
        near(
            totals.recession_m,
            mass / (Continuous::AREA * 1000.0),
            1e-13,
        );
        let gas_energy = gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
            - initial_gas.mass_kg() * initial_gas.state().specific_internal_energy_j_kg();
        near(
            totals.vapor_energy_j - totals.sensible_heat_into_liquid_j,
            gas_energy,
            1e-7,
        );
        near(
            thermal.energy_j() - original.energy_j() + totals.vapor_energy_j,
            report.external_heat_j + totals.sensible_heat_into_liquid_j,
            1e-7,
        );
        near(
            thermal.energy_j() - original.energy_j() + gas_energy,
            20.0,
            1e-7,
        );
        near(report.external_heat_j, 20.0, 1e-9);
        near(report.mass_residual_kg, 0.0, 1e-13);
        near(report.energy_residual_j, 0.0, 1e-7);
        assert!(report.absolute_energy_residual_j <= config.solve.thermal.total_energy_tolerance_j);
        assert!(report.absolute_mass_residual_kg <= config.solve.phase.gas.mass_tolerance_kg);
        near(
            gas.species_mass_kg(&carrier()).unwrap(),
            oracle.carrier_mass,
            1e-14,
        );
        assert_eq!(totals.intervals, report.accepted_intervals);
        assert_eq!(report.accepted_intervals, 2 * report.accepted_comparisons);
        assert_eq!(
            report.comparisons,
            report.accepted_comparisons + report.accuracy_rejections
        );
        assert!(report.accepted_intervals <= original.work_budget().max_substeps);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert!(
            report
                .last_error
                .as_ref()
                .is_some_and(|error| error.normalized <= 1.0)
        );
        assert!(report.last_step.is_some());
        assert!(report.newton_iterations <= config.solve.thermal.maximum_iterations);
        assert!(report.residual_evaluations <= config.solve.thermal.maximum_evaluations);
        assert!(report.gas_iterations <= config.maximum_gas_iterations);
        assert_eq!(
            candidate.payload(),
            Some(&(
                Continuous::DURATION,
                state[2],
                state[1],
                thermal.cells()[&CellId(1)].mass_kg
            ))
        );
        assert!(std::ptr::eq(candidate.catalog(), &source));
        outcomes.push((
            error,
            report.accepted_comparisons,
            report.accuracy_rejections,
        ));
    }
    for coordinate in 0..3 {
        assert!(
            outcomes[1].0[coordinate] < 0.35 * outcomes[0].0[coordinate],
            "adaptive errors: {outcomes:?}"
        );
    }
    assert!(outcomes[1].1 > outcomes[0].1);
    assert!(outcomes[1].2 > 0);

    // Make both thermal tolerances ineffective, then refuse the unchanged
    // full span solely because the actual liquid/species inventories disagree.
    let mut mass_only = policy(1000.0, 1e-12);
    mass_only.minimum_dt_s = 2.0;
    let error = advance_deposition_with_gas_adaptive(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        2.0,
        mass_only,
        empty_heat,
        || Ok(()),
    )
    .err()
    .expect("mass discrepancy must refuse the indivisible comparison span");
    match error {
        AdaptiveImplicitGasError::TimeAccuracy {
            time_s,
            span_s,
            error,
        } => {
            assert_eq!(time_s, 0.0);
            assert_eq!(span_s, 2.0);
            assert!(
                error.maximum_temperature_difference_k < mass_only.temperature_absolute_tolerance_k
            );
            assert!(
                error.gas_temperature_difference_k < mass_only.temperature_absolute_tolerance_k
            );
            assert!(
                error.maximum_enthalpy_difference_j_kg < mass_only.enthalpy_absolute_tolerance_j_kg
            );
            assert!(error.maximum_liquid_mass_difference_kg > mass_only.mass_absolute_tolerance_kg);
            assert!(
                error.maximum_species_mass_difference_kg > mass_only.mass_absolute_tolerance_kg
            );
            assert!(error.normalized > 1.0);
        }
        other => panic!("unexpected mass-only refusal: {other:?}"),
    }
    unchanged(&original, &snapshot, &initial_gas, &gas_snapshot);
}

#[test]
fn rejected_comparisons_spend_shared_work_and_every_refusal_keeps_both_original_owners() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let original = thermal_with_budget(2.0, 10_000, vec![]);
    let initial_gas = gas(&source, 340.0, 5000.0, 0.001);
    let snapshot = original.clone();
    let gas_snapshot = initial_gas.clone();
    let surfaces = [surface(
        1,
        &binding,
        Continuous::MASS_COEFFICIENT,
        Continuous::HEAT_COEFFICIENT,
    )];
    let config = policy(0.0005, 1e-10);
    let (comparison_evaluations, comparison_gas_iterations) =
        first_comparison_work(&original, &surfaces, &initial_gas, &source, config.solve);

    let calls = Cell::new(0usize);
    let mut polls = 0usize;
    let accepted = advance_deposition_with_gas_adaptive(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        2.0,
        config,
        |_, _, _| {
            calls.set(calls.get() + 1);
            Ok(coefficients(calls.get()))
        },
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(accepted.payload(), Some(&calls.get()));
    assert_eq!(calls.get(), accepted.report().residual_evaluations);
    assert!(accepted.report().accuracy_rejections > 0);
    assert!(accepted.report().residual_evaluations > comparison_evaluations);
    assert!(accepted.report().gas_iterations > comparison_gas_iterations);

    // Independently measured work of coarse + half + half is spent even when
    // that entire physical comparison is rejected and the controller retries.
    for budget in 0..3 {
        let mut limited = config;
        match budget {
            0 => limited.maximum_comparisons = 1,
            1 => limited.solve.thermal.maximum_evaluations = comparison_evaluations,
            _ => limited.maximum_gas_iterations = comparison_gas_iterations,
        }
        calls.set(0);
        let result = advance_deposition_with_gas_adaptive(
            &original,
            &surfaces,
            &initial_gas,
            &source,
            2.0,
            limited,
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok(coefficients(()))
            },
            || Ok(()),
        );
        match budget {
            0 => assert!(matches!(
                result,
                Err(AdaptiveImplicitGasError::ComparisonBudget)
            )),
            1 => assert!(matches!(
                result,
                Err(AdaptiveImplicitGasError::Implicit(
                    ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
                        DepositionError::Budget("adaptive finite-gas whole-window solver work")
                    ))
                ))
            )),
            _ => assert!(matches!(
                result,
                Err(AdaptiveImplicitGasError::Implicit(
                    ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
                        DepositionError::Budget("adaptive finite-gas whole-window gas iterations")
                    ))
                ))
            )),
        }
        assert_eq!(calls.get(), comparison_evaluations);
        unchanged(&original, &snapshot, &initial_gas, &gas_snapshot);
    }

    calls.set(0);
    let result = advance_deposition_with_gas_adaptive(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        2.0,
        config,
        |thermal, _, _| {
            calls.set(calls.get() + 1);
            if calls.get() == comparison_evaluations + 1 {
                assert!(thermal.time_s() > 0.0 && thermal.time_s() < 2.0);
                Err(DepositionError::InvalidInput("adaptive original heat source refused").into())
            } else {
                Ok(coefficients(()))
            }
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(AdaptiveImplicitGasError::Implicit(
            ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
                DepositionError::InvalidInput("adaptive original heat source refused")
            ))
        ))
    ));
    assert_eq!(calls.get(), comparison_evaluations + 1);
    unchanged(&original, &snapshot, &initial_gas, &gas_snapshot);

    // Reserve both fine intervals against the ORIGINAL thermal work ceiling.
    // A coarse solve may precede the refusal, but an unreservable fine pair
    // must not start, and accepted pairs cannot renew the window allowance.
    for (thermal_budget, target, forbidden_half_endpoint) in [(1, 2.0, 1.0), (2, 4.0, 3.0)] {
        let limited_thermal = thermal_with_budget(2.0, thermal_budget, vec![]);
        let limited_snapshot = limited_thermal.clone();
        let mut queried_times = Vec::new();
        let result = advance_deposition_with_gas_adaptive(
            &limited_thermal,
            &surfaces,
            &initial_gas,
            &source,
            target,
            policy(1000.0, 1e-3),
            |thermal, _, _| {
                queried_times.push(thermal.time_s());
                Ok(coefficients(()))
            },
            || Ok(()),
        );
        assert!(matches!(
            result,
            Err(AdaptiveImplicitGasError::Implicit(
                ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
                    DepositionError::Budget(_)
                ))
            ))
        ));
        assert!(!queried_times.contains(&forbidden_half_endpoint));
        unchanged(
            &limited_thermal,
            &limited_snapshot,
            &initial_gas,
            &gas_snapshot,
        );
    }

    for cancel_at in [1, polls] {
        let mut count = 0usize;
        let result = advance_deposition_with_gas_adaptive(
            &original,
            &surfaces,
            &initial_gas,
            &source,
            2.0,
            config,
            empty_heat,
            || {
                count += 1;
                if count == cancel_at {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert!(matches!(
            result,
            Err(AdaptiveImplicitGasError::Implicit(
                ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
                    DepositionError::Cancelled
                ))
            ))
        ));
        assert_eq!(count, cancel_at);
        unchanged(&original, &snapshot, &initial_gas, &gas_snapshot);
    }
}

#[test]
fn adaptive_window_rejects_interior_or_unprocessed_events_before_sources_and_preserves_target_events()
 {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial_gas = gas(&source, 300.0, 1000.0, 0.001);
    let mut changed = cell(1, 300.0, 1.0, Continuous::AREA).boundary;
    changed.source_w = 3.0;
    let original = thermal_with_budget(
        1.0,
        10_000,
        vec![TimedEvent {
            time_s: 1.0,
            event: Event::Boundary(CellId(1), changed),
        }],
    );
    let snapshot = original.clone();
    let gas_snapshot = initial_gas.clone();
    let surfaces = [surface(
        1,
        &binding,
        Continuous::MASS_COEFFICIENT,
        Continuous::HEAT_COEFFICIENT,
    )];
    let mut config = policy(0.05, 1e-8);
    config.initial_dt_s = 0.25;
    config.maximum_dt_s = 0.25;
    let calls = Cell::new(0usize);
    let result = advance_deposition_with_gas_adaptive(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        2.0,
        config,
        |_, _, _| {
            calls.set(calls.get() + 1);
            Ok(coefficients(()))
        },
        || Ok(()),
    );
    assert!(result.is_err());
    assert_eq!(
        calls.get(),
        0,
        "the entire event-free window must be checked before advancing its first subinterval"
    );

    let candidate = advance_deposition_with_gas_adaptive(
        &original,
        &surfaces,
        &initial_gas,
        &source,
        1.0,
        config,
        |thermal, _, _| Ok(coefficients(thermal.time_s())),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(candidate.thermal().time_s(), 1.0);
    assert_eq!(candidate.payload(), Some(&1.0));
    assert_eq!(
        candidate.thermal().pending_events(),
        original.pending_events()
    );
    assert_eq!(
        candidate.thermal().cells()[&CellId(1)].boundary.source_w,
        1.0
    );
    near(candidate.report().external_heat_j, 1.0, 1e-9);
    assert!(candidate.report().accepted_intervals > 2);
    assert_eq!(
        candidate.report().surface_totals[&CellId(1)].intervals,
        candidate.report().accepted_intervals
    );
    calls.set(0);
    let result = advance_deposition_with_gas_adaptive(
        candidate.thermal(),
        &surfaces,
        candidate.gas(),
        &source,
        2.0,
        config,
        |_, _, _| {
            calls.set(calls.get() + 1);
            Ok(coefficients(()))
        },
        || Ok(()),
    );
    assert!(result.is_err());
    assert_eq!(calls.get(), 0);
    assert_eq!(original.pending_events(), snapshot.pending_events());
    unchanged(&original, &snapshot, &initial_gas, &gas_snapshot);
}
