//! G1 independent heat/fusion balances; G3 split refinement and owner equivalence.
use super::*;
use crate::evaporation::gas::{
    equilibrium::tests::{CP, FUSION, MELTING, TRANSFER, VAPORIZATION, bulk, catalog, curve, near},
    heat::advance_surface_with_gas_and_heat,
    tests::{Phase, binding, controls, gas},
};
use crate::phase::EquilibriumEnthalpyPhaseCurve;

fn heat_controls(chart: &EquilibriumEnthalpyPhaseCurve) -> EquilibriumGasHeatSolve {
    EquilibriumGasHeatSolve {
        bulk_enthalpy_range_j_kg: [
            chart.knots()[0].specific_enthalpy_j_kg,
            chart.knots().last().unwrap().specific_enthalpy_j_kg,
        ],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}

#[test]
fn finite_gas_heat_melts_and_freezes_an_isothermal_body_with_one_energy_account() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    for gas_t in [340.0, 260.0] {
        let gas = gas(&source, gas_t, 1000.0, 0.01);
        let result = advance_equilibrium_surface_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &source,
            0.0,
            60.0,
            2.0,
            controls(),
            heat_controls(&chart),
            || false,
        )
        .unwrap();
        // P*V/T = nR. The fixture's fixed mole fractions are 0.99 carrier
        // (Cv/R=2.5), 0.01 volatile (Cv/R=3), independent of the heat solver.
        let gas_capacity = 100_000.0 * 0.01 / gas_t * (0.99 * 2.5 + 0.01 * 3.0);
        let exposure = 60.0 * initial.area_m2 * 2.0;
        let expected_heat = exposure * (gas_t - MELTING) / (1.0 + exposure / gas_capacity);
        let heat = result.heat();
        let step = result.step();
        near(heat.heat_into_liquid_j, expected_heat, 5e-9);
        near(
            step.gas_temperature_k[1],
            gas_t - expected_heat / gas_capacity,
            1e-8,
        );
        assert_eq!(step.surface_temperature_k, [MELTING, MELTING]);
        assert_eq!(step.transferred_mass_kg, 0.0);
        assert_eq!(step.energy_into_gas_j, 0.0);
        assert_eq!(result.bulk().mass_kg, initial.mass_kg);
        assert_eq!(result.gas().mass_kg(), gas.mass_kg());
        assert_eq!(result.gas().model(), gas.model());
        assert_eq!(result.gas().point().axes().get("local-lot"), Some(&11.0));
        near(
            result.phase().freezing_mass_kg,
            -heat.heat_into_liquid_j / FUSION,
            2e-17,
        );
        near(
            result.phase().fusion_energy_j,
            -heat.heat_into_liquid_j,
            2e-12,
        );
        assert_eq!(result.transfer().phase().freezing_mass_kg, 0.0);
        near(
            step.liquid_enthalpy_change_j + step.gas_internal_energy_change_j,
            0.0,
            1e-8,
        );
        assert!(result.sum_absolute_energy_residual_j() <= controls().gas.energy_tolerance_j);
        near(
            result.phase().liquid_partition_residual_kg,
            0.0,
            controls().gas.mass_tolerance_kg,
        );
    }
}

#[test]
fn simultaneous_heat_and_vapor_streams_share_the_original_bulk_fusion_storage() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    for (gas_t, vapor_pressure) in [(330.0, 0.0), (280.0, 5000.0)] {
        let gas = gas(&source, gas_t, vapor_pressure, 0.01);
        let result = advance_equilibrium_surface_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &source,
            TRANSFER,
            30.0,
            1.0,
            controls(),
            heat_controls(&chart),
            || false,
        )
        .unwrap();
        let step = result.step();
        let dm = step.transferred_mass_kg;
        let q = result.heat().heat_into_liquid_j;
        assert!(q != 0.0 && dm != 0.0);
        assert_eq!(step.surface_temperature_k, [MELTING, MELTING]);
        let phase_gas_t = result.transfer().exchange().step().gas_temperature_k[0];
        let carried_latent_gap = VAPORIZATION
            + if dm < 0.0 {
                CP * (phase_gas_t - MELTING)
            } else {
                0.0
            };
        near(
            result.phase().fusion_energy_j,
            dm * carried_latent_gap - q,
            3e-12,
        );
        near(
            result.transfer().phase().fusion_energy_j,
            dm * carried_latent_gap,
            3e-12,
        );
        near(
            step.gas_internal_energy_change_j,
            step.energy_into_gas_j - q,
            1e-8,
        );
        near(
            step.liquid_enthalpy_change_j + step.gas_internal_energy_change_j,
            0.0,
            1e-8,
        );
        assert!(result.sum_absolute_energy_residual_j() <= controls().gas.energy_tolerance_j);
        let original_total = initial.mass_kg * initial.specific_enthalpy_j_kg
            + gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
        let final_total = result.bulk().mass_kg * result.bulk().specific_enthalpy_j_kg
            + result.gas().mass_kg() * result.gas().state().specific_internal_energy_j_kg();
        near(final_total, original_total, 1e-7);
        near(step.mass_residual_kg, 0.0, controls().gas.mass_tolerance_kg);
    }
}

#[test]
fn zero_heat_coefficient_exactly_delegates_even_with_unused_invalid_heat_controls() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let gas = gas(&source, 310.0, 1000.0, 0.01);
    let phase = advance_equilibrium_surface(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &gas,
        &source,
        TRANSFER,
        0.5,
        controls(),
        || false,
    )
    .unwrap();
    let result = advance_equilibrium_surface_with_heat(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &gas,
        &source,
        TRANSFER,
        0.0,
        0.5,
        controls(),
        EquilibriumGasHeatSolve {
            bulk_enthalpy_range_j_kg: [f64::NAN; 2],
            temperature_tolerance_k: f64::NAN,
            maximum_iterations: 0,
        },
        || false,
    )
    .unwrap();
    assert_eq!(result.transfer(), &phase);
    assert_eq!(result.phase(), phase.phase());
    assert_eq!(result.step(), phase.exchange().step());
    assert_eq!(result.heat(), &GasFilmHeatStep::default());
    assert_eq!(result.phase_initial_bulk(), initial);
}

#[test]
fn fully_liquid_body_reuses_the_same_reciprocal_heat_and_vapor_owner() {
    let source = catalog(0.0);
    let chart = curve();
    let mut initial = bulk(&chart, 1.0);
    initial.specific_enthalpy_j_kg = CP * 340.0 - VAPORIZATION;
    let gas = gas(&source, 300.0, 1000.0, 0.01);
    let mut heat_budget = heat_controls(&chart);
    heat_budget.bulk_enthalpy_range_j_kg[0] = chart.knots()[3].specific_enthalpy_j_kg;
    let result = advance_equilibrium_surface_with_heat(
        EquilibriumEvaporation::new(&chart).unwrap(),
        &Phase,
        &binding(0.0),
        initial,
        &gas,
        &source,
        TRANSFER,
        60.0,
        2.0,
        controls(),
        heat_budget,
        || false,
    )
    .unwrap();
    let legacy = advance_surface_with_gas_and_heat(
        &Phase,
        &binding(0.0),
        initial,
        &gas,
        &source,
        |h| {
            chart
                .temperature_at_specific_enthalpy(h)
                .map_err(|_| EvaporationError::Invalid("original liquid source"))
        },
        TRANSFER,
        60.0,
        2.0,
        controls(),
        GasFilmHeatSolve {
            pure_liquid_enthalpy_range_j_kg: heat_budget.bulk_enthalpy_range_j_kg,
            temperature_tolerance_k: heat_budget.temperature_tolerance_k,
            maximum_iterations: heat_budget.maximum_iterations,
        },
        || false,
    )
    .unwrap();
    assert_eq!(result.bulk(), legacy.liquid());
    assert_eq!(result.gas(), legacy.gas());
    assert_eq!(result.heat(), legacy.heat());
    assert_eq!(result.step(), legacy.step());
    assert_eq!(result.phase_initial_bulk(), legacy.phase_initial_liquid());
    assert_eq!(result.phase().freezing_mass_kg, 0.0);
}

#[test]
fn combined_temperature_departure_is_checked_after_both_individually_admitted_substeps() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let mut initial = bulk(&chart, 1.0);
    initial.specific_enthalpy_j_kg = CP * 330.0 - VAPORIZATION;
    let gas = gas(&source, 300.0, 0.0, 1.0);
    let run = |budget| {
        advance_equilibrium_surface_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &source,
            TRANSFER,
            60.0,
            1.0,
            budget,
            heat_controls(&chart),
            || false,
        )
    };
    let baseline = run(controls()).unwrap();
    let temperatures = baseline.step().surface_temperature_k;
    let phase_start = baseline
        .transfer()
        .phase()
        .initial
        .equilibrium
        .temperature_k();
    let heat_departure = (phase_start - temperatures[0]).abs();
    let phase_departure = (temperatures[1] - phase_start).abs();
    let complete_departure = (temperatures[1] - temperatures[0]).abs();
    assert!(complete_departure > heat_departure && heat_departure > phase_departure);
    let mut budget = controls();
    budget.maximum_temperature_increment_k = 0.5 * (heat_departure + complete_departure);
    assert!(matches!(
        run(budget),
        Err(GasEvaporationError::StepSizeLimit {
            reason: "combined heat and phase temperature departure",
            ..
        })
    ));
}

#[test]
fn complete_heat_and_phase_steps_converge_under_physical_time_refinement() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let original_gas = gas(&source, 330.0, 1000.0, 0.01);
    let run = |count: usize| {
        let mut bulk = initial;
        let mut gas = original_gas.clone();
        for _ in 0..count {
            let result = advance_equilibrium_surface_with_heat(
                inventory,
                &Phase,
                &binding(0.0),
                bulk,
                &gas,
                &source,
                0.2,
                80.0,
                8.0 / count as f64,
                controls(),
                heat_controls(&chart),
                || false,
            )
            .unwrap();
            assert_eq!(result.step().surface_temperature_k, [MELTING, MELTING]);
            assert!(result.step().transferred_mass_kg > 0.0);
            assert!(result.heat().heat_into_liquid_j > 0.0);
            let (next_bulk, next_gas, _, _, _, _) = result.into_parts();
            bulk = next_bulk;
            gas = next_gas;
        }
        (
            gas.state().state().temperature,
            inventory.state(bulk).unwrap().liquid_mass_kg,
        )
    };
    let reference = run(128);
    let errors = [4, 8, 16].map(|count| {
        let endpoint = run(count);
        (endpoint.0 - reference.0).abs() + 1e5 * (endpoint.1 - reference.1).abs()
    });
    assert!(
        errors[0] > 1e-5,
        "inactive coupled time refinement: {errors:?}"
    );
    assert!(
        errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn original_source_brackets_liquid_seed_and_late_cancellation_remain_atomic() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let gas = gas(&source, 260.0, 1000.0, 0.01);
    let saved = gas.clone();
    let mut wrong = heat_controls(&chart);
    wrong.bulk_enthalpy_range_j_kg[0] -= 1.0;
    assert!(matches!(
        advance_equilibrium_surface_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &source,
            0.0,
            60.0,
            1.0,
            controls(),
            wrong,
            || false,
        ),
        Err(GasEvaporationError::Invalid(_))
    ));
    assert!(
        advance_equilibrium_surface_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            bulk(&chart, 1e-6),
            &gas,
            &source,
            0.0,
            60.0,
            1.0,
            controls(),
            heat_controls(&chart),
            || false,
        )
        .is_err()
    );

    let mut calls = 0;
    advance_equilibrium_surface_with_heat(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &gas,
        &source,
        TRANSFER,
        30.0,
        1.0,
        controls(),
        heat_controls(&chart),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    assert!(calls > 1);
    let mut attempt = 0;
    assert!(matches!(
        advance_equilibrium_surface_with_heat(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &source,
            TRANSFER,
            30.0,
            1.0,
            controls(),
            heat_controls(&chart),
            || {
                attempt += 1;
                attempt == calls
            },
        ),
        Err(GasEvaporationError::Cancelled)
    ));
    assert_eq!(gas, saved);
    assert_eq!(initial, bulk(&chart, 0.5));
}
