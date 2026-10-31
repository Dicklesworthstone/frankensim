//! Independent polynomial energy accounts for the manufactured NASA-9 sources.
use super::*;
use crate::evaporation::gas::tests::{Phase, binding, catalog, controls, gas, liquid, temperature};

fn heat_controls() -> GasFilmHeatSolve {
    GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: [
            liquid(200.0).specific_enthalpy_j_kg,
            liquid(900.0).specific_enthalpy_j_kg,
        ],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn gas_capacity_coefficients() -> (f64, f64) {
    // Fixture mole fractions are [0.99,0.01]; Cp/R is 3.5 and 4+0.001*T.
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / (0.99 * 0.02897 + 0.01 * 0.018015);
    (r * (0.99 * 2.5 + 0.01 * 3.0), r * 0.01 * 0.001)
}
fn independent_gas_du(a: f64, b: f64) -> f64 {
    let (cv0, slope) = gas_capacity_coefficients();
    (b - a) * (cv0 + 0.5 * slope * (a + b))
}
fn independent_liquid_dh(a: f64, b: f64) -> f64 {
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    r * (b - a) * (4.0 + 0.5 * 0.001 * (a + b))
}

#[test]
fn nonlinear_original_calorics_exchange_reciprocal_heat_at_zero_mass_flux() {
    let source = catalog(0.0);
    for (liquid_t, gas_t) in [(300.0, 340.0), (340.0, 300.0)] {
        let original_liquid = liquid(liquid_t);
        let original_gas = gas(&source, gas_t, 1000.0, 0.01);
        let candidate = advance_surface_with_gas_and_heat(
            &Phase,
            &binding(0.0),
            original_liquid,
            &original_gas,
            &source,
            temperature,
            0.0,
            60.0,
            2.0,
            controls(),
            heat_controls(),
            || false,
        )
        .unwrap();
        let heat = candidate.heat();
        let step = candidate.step();
        let final_liquid_t = step.surface_temperature_k[1];
        let final_gas_t = step.gas_temperature_k[1];
        assert_eq!(step.transferred_mass_kg, 0.0);
        assert_eq!(step.energy_into_gas_j, 0.0);
        assert_eq!(candidate.liquid().mass_kg, original_liquid.mass_kg);
        assert_eq!(candidate.gas().mass_kg(), original_gas.mass_kg());
        assert_eq!(candidate.gas().volume_m3(), original_gas.volume_m3());
        assert_eq!(candidate.gas().model(), original_gas.model());
        assert_eq!(candidate.gas().point().axes().get("local-lot"), Some(&11.0));
        let liquid_energy =
            original_liquid.mass_kg * independent_liquid_dh(liquid_t, final_liquid_t);
        let gas_energy = original_gas.mass_kg() * independent_gas_du(gas_t, final_gas_t);
        near(heat.heat_into_liquid_j, liquid_energy, 2e-9);
        near(heat.gas_internal_energy_change_j, gas_energy, 2e-9);
        near(liquid_energy + gas_energy, 0.0, 4e-9);
        near(
            step.liquid_enthalpy_change_j + step.gas_internal_energy_change_j,
            0.0,
            1e-8,
        );
        near(
            (final_gas_t - final_liquid_t) - heat.heat_into_liquid_j / (60.0 * 0.001 * 2.0),
            0.0,
            heat_controls().temperature_tolerance_k,
        );
        assert!(heat.heat_into_liquid_j * (gas_t - liquid_t) > 0.0);
        assert!((final_gas_t - final_liquid_t) * (gas_t - liquid_t) > 0.0);
        assert!((final_gas_t - final_liquid_t).abs() < (gas_t - liquid_t).abs());
        // Both original capacities vary with temperature; the gas EOS also
        // updates pressure at fixed mass/volume instead of leaving it stale.
        near(
            candidate.gas().state().state().pressure / original_gas.state().state().pressure,
            final_gas_t / gas_t,
            3e-15,
        );
    }
}

#[test]
fn zero_heat_coefficient_exactly_delegates_to_the_original_phase_path() {
    let source = catalog(0.0);
    let original_gas = gas(&source, 310.0, 1000.0, 0.01);
    let original_liquid = liquid(300.0);
    let phase = advance_surface_with_gas(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        0.2,
        controls(),
        || false,
    )
    .unwrap();
    let coupled = advance_surface_with_gas_and_heat(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        0.0,
        0.2,
        controls(),
        heat_controls(),
        || false,
    )
    .unwrap();
    assert_eq!(coupled.phase(), &phase);
    assert_eq!(coupled.step(), phase.step());
    assert_eq!(coupled.heat(), &GasFilmHeatStep::default());
    assert_eq!(coupled.gas(), phase.gas());
    assert_eq!(coupled.liquid(), phase.liquid());
}

#[test]
fn stiff_implicit_heat_approaches_the_nonlinear_common_temperature_without_crossing() {
    let source = catalog(0.0);
    let original_gas = gas(&source, 340.0, 1000.0, 0.01);
    let original_liquid = liquid(300.0);
    let mut phase = controls();
    phase.maximum_temperature_increment_k = 100.0;
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let (a_g, b_g) = gas_capacity_coefficients();
    let a = original_liquid.mass_kg * 4.0 * r + original_gas.mass_kg() * a_g;
    let b = original_liquid.mass_kg * 0.001 * r + original_gas.mass_kg() * b_g;
    let capacity = a + 300.0 * b;
    let excess = original_gas.mass_kg() * independent_gas_du(300.0, 340.0);
    let equilibrium =
        300.0 + 2.0 * excess / (capacity + (capacity * capacity + 2.0 * b * excess).sqrt());
    let mut previous_gap = 40.0;
    for duration in [100.0, 1000.0, 10000.0] {
        let candidate = advance_surface_with_gas_and_heat(
            &Phase,
            &binding(0.0),
            original_liquid,
            &original_gas,
            &source,
            temperature,
            0.0,
            1000.0,
            duration,
            phase,
            heat_controls(),
            || false,
        )
        .unwrap();
        let tl = candidate.step().surface_temperature_k[1];
        let tg = candidate.step().gas_temperature_k[1];
        assert!(tl > 300.0 && tl <= equilibrium && equilibrium <= tg && tg < 340.0);
        assert!(tg - tl < previous_gap);
        near(
            candidate.step().energy_residual_j,
            0.0,
            phase.gas.energy_tolerance_j,
        );
        previous_gap = tg - tl;
    }
    assert!(previous_gap < 0.03);
}

#[test]
fn complete_heat_and_evaporation_split_converges_under_physical_time_refinement() {
    let source = catalog(0.0);
    let original_gas = gas(&source, 330.0, 1000.0, 0.005);
    let run = |count: usize| {
        let mut gas = original_gas.clone();
        let mut liquid = liquid(300.0);
        let mut transferred = 0.0;
        for _ in 0..count {
            let result = advance_surface_with_gas_and_heat(
                &Phase,
                &binding(0.0),
                liquid,
                &gas,
                &source,
                temperature,
                0.2,
                80.0,
                8.0 / count as f64,
                controls(),
                heat_controls(),
                || false,
            )
            .unwrap();
            let (next_liquid, next_gas, _, step, heat) = result.into_parts();
            assert!(heat.heat_into_liquid_j > 0.0);
            transferred += step.transferred_mass_kg;
            liquid = next_liquid;
            gas = next_gas;
        }
        assert!(transferred > 0.0);
        (
            temperature(liquid.specific_enthalpy_j_kg).unwrap(),
            gas.state().state().temperature,
            liquid.mass_kg,
        )
    };
    let reference = run(256);
    let errors: Vec<_> = [8, 16, 32]
        .into_iter()
        .map(|count| {
            let endpoint = run(count);
            (endpoint.0 - reference.0).abs()
                + (endpoint.1 - reference.1).abs()
                + 1e6 * (endpoint.2 - reference.2).abs()
        })
        .collect();
    assert!(errors[0] > 1e-5, "inactive coupled refinement: {errors:?}");
    assert!(
        errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn heat_source_limits_and_late_cancellation_publish_neither_original_owner() {
    let source = catalog(0.0);
    let original_gas = gas(&source, 340.0, 1000.0, 0.01);
    let original_liquid = liquid(300.0);
    let before = original_gas.clone();
    let mut small = controls();
    small.maximum_temperature_increment_k = 0.01;
    assert!(matches!(
        advance_surface_with_gas_and_heat(
            &Phase,
            &binding(0.0),
            original_liquid,
            &original_gas,
            &source,
            temperature,
            0.0,
            1000.0,
            1.0,
            small,
            heat_controls(),
            || false,
        ),
        Err(GasEvaporationError::StepSizeLimit { .. })
    ));
    let narrow = GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: [
            original_liquid.specific_enthalpy_j_kg - 1.0,
            original_liquid.specific_enthalpy_j_kg + 1.0,
        ],
        ..heat_controls()
    };
    assert!(matches!(
        advance_surface_with_gas_and_heat(
            &Phase,
            &binding(0.0),
            original_liquid,
            &original_gas,
            &source,
            temperature,
            0.0,
            1000.0,
            1.0,
            controls(),
            narrow,
            || false,
        ),
        Err(GasEvaporationError::HeatOutsideBracket { .. })
    ));
    // A small liquid inventory makes its admitted h endpoint bound the search.
    // This datum and mass reproduce a one-ulp excursion if m*Delta-h is divided
    // back to reconstruct that endpoint; the real heat root remains interior.
    let shift = liquid(280.0).specific_enthalpy_j_kg - 29000.123;
    let small_liquid = LiquidSurfaceState {
        mass_kg: 0.00002345678991,
        specific_enthalpy_j_kg: liquid(280.0).specific_enthalpy_j_kg - shift,
        ..liquid(280.0)
    };
    let branch = GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: [
            small_liquid.specific_enthalpy_j_kg,
            liquid(330.0).specific_enthalpy_j_kg - shift,
        ],
        ..heat_controls()
    };
    let upper = branch.pure_liquid_enthalpy_range_j_kg[1];
    let gap = upper - small_liquid.specific_enthalpy_j_kg;
    assert!(
        small_liquid.specific_enthalpy_j_kg + (small_liquid.mass_kg * gap) / small_liquid.mass_kg
            > upper
    );
    let warm = gas(&source, 320.0, 1000.0, 0.01);
    let mut wide = controls();
    wide.maximum_temperature_increment_k = 100.0;
    let admitted = advance_surface_with_gas_and_heat(
        &Phase,
        &GasLiquidBinding {
            liquid_enthalpy_offset_j_kg: shift,
            ..binding(0.0)
        },
        small_liquid,
        &warm,
        &source,
        |h| temperature(h + shift),
        0.0,
        100.0,
        1.0,
        wide,
        branch,
        || false,
    )
    .unwrap();
    assert!(admitted.liquid().specific_enthalpy_j_kg < upper);
    assert!(admitted.heat().heat_into_liquid_j > 0.0);
    let wrong = GasLiquidBinding {
        liquid_enthalpy_offset_j_kg: 100.0,
        ..binding(0.0)
    };
    assert!(matches!(
        advance_surface_with_gas_and_heat(
            &Phase,
            &wrong,
            original_liquid,
            &original_gas,
            &source,
            temperature,
            0.0,
            1000.0,
            1.0,
            controls(),
            heat_controls(),
            || false,
        ),
        Err(GasEvaporationError::CaloricMismatch { .. })
    ));
    let mut calls = 0;
    advance_surface_with_gas_and_heat(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        60.0,
        0.5,
        controls(),
        heat_controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        let rejected = advance_surface_with_gas_and_heat(
            &Phase,
            &binding(0.0),
            original_liquid,
            &original_gas,
            &source,
            temperature,
            0.02,
            60.0,
            0.5,
            controls(),
            heat_controls(),
            || {
                count += 1;
                count == stop
            },
        );
        assert!(matches!(
            rejected,
            Err(GasEvaporationError::Cancelled)
                | Err(GasEvaporationError::Gas(
                    super::super::GasEnergyError::Cancelled
                ))
        ));
        assert_eq!(original_gas, before);
        assert_eq!(original_liquid, liquid(300.0));
    }
}

#[test]
fn native_mass_selection_runs_once_on_the_phase_flux_after_sensible_heat() {
    let source = catalog(0.0);
    let original_gas = gas(&source, 340.0, 1000.0, 0.01);
    let original_liquid = liquid(300.0);
    let initial_phase = advance_surface_with_gas(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        0.5,
        controls(),
        || false,
    )
    .unwrap();
    let baseline = advance_surface_with_gas_and_heat(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        60.0,
        0.5,
        controls(),
        heat_controls(),
        || false,
    )
    .unwrap();
    let requested = baseline.step().requested_transferred_mass_kg;
    assert!(requested > 0.0 && baseline.heat().heat_into_liquid_j > 0.0);
    assert_ne!(
        requested,
        initial_phase.step().requested_transferred_mass_kg
    );
    assert!(
        baseline.phase().step().surface_temperature_k[0] > baseline.step().surface_temperature_k[0]
    );
    let selected = requested.next_down();
    let mut calls = 0;
    let candidate = advance_surface_with_gas_and_heat_mass_selection(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        60.0,
        0.5,
        controls(),
        heat_controls(),
        |mass| {
            calls += 1;
            assert_eq!(mass, requested);
            Ok(selected)
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls, 1);
    assert_eq!(candidate.heat(), baseline.heat());
    assert_eq!(
        temperature(candidate.phase_initial_liquid().specific_enthalpy_j_kg).unwrap(),
        candidate.phase().step().surface_temperature_k[0]
    );
    assert_eq!(candidate.step().requested_transferred_mass_kg, requested);
    assert_eq!(candidate.step().transferred_mass_kg, selected);
    assert_eq!(candidate.phase().step().transferred_mass_kg, selected);
    assert_eq!(
        candidate.liquid().mass_kg,
        original_liquid.mass_kg - selected
    );
    assert_eq!(
        candidate.step().gas_temperature_k[0],
        original_gas.state().state().temperature
    );
    near(
        candidate.step().energy_residual_j,
        0.0,
        controls().gas.energy_tolerance_j,
    );
    near(
        candidate.step().liquid_enthalpy_change_j + candidate.step().gas_internal_energy_change_j,
        0.0,
        controls().gas.energy_tolerance_j,
    );
}

#[test]
fn zero_heat_exactly_delegates_to_selected_phase_without_using_heat_controls() {
    let source = catalog(0.0);
    let original_gas = gas(&source, 310.0, 1000.0, 0.01);
    let original_liquid = liquid(300.0);
    let phase = advance_surface_with_gas_mass_selection(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        0.2,
        controls(),
        |requested| Ok(requested.next_up()),
        || false,
    )
    .unwrap();
    let mut calls = 0;
    let candidate = advance_surface_with_gas_and_heat_mass_selection(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        0.0,
        0.2,
        controls(),
        GasFilmHeatSolve {
            pure_liquid_enthalpy_range_j_kg: [f64::NAN; 2],
            temperature_tolerance_k: f64::NAN,
            maximum_iterations: 0,
        },
        |requested| {
            calls += 1;
            Ok(requested.next_up())
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls, 1);
    assert_eq!(candidate.phase(), &phase);
    assert_eq!(candidate.step(), phase.step());
    assert_eq!(candidate.heat(), &GasFilmHeatStep::default());
    assert_eq!(candidate.gas(), phase.gas());
    assert_eq!(candidate.liquid(), phase.liquid());
}

#[test]
fn selector_failure_invalid_mass_and_early_cancellation_publish_no_heat_or_phase() {
    let source = catalog(0.0);
    let original_gas = gas(&source, 340.0, 1000.0, 0.01);
    let original_liquid = liquid(300.0);
    let before = original_gas.clone();
    let refusal = GasEvaporationError::Invalid("native inventory refuses transfer");
    let mut calls = 0;
    let rejected = advance_surface_with_gas_and_heat_mass_selection(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        60.0,
        0.5,
        controls(),
        heat_controls(),
        |_| {
            calls += 1;
            Err(refusal.clone())
        },
        || false,
    );
    assert_eq!(calls, 1);
    assert_eq!(rejected, Err(refusal));
    assert!(matches!(
        advance_surface_with_gas_and_heat_mass_selection(
            &Phase,
            &binding(0.0),
            original_liquid,
            &original_gas,
            &source,
            temperature,
            0.02,
            60.0,
            0.5,
            controls(),
            heat_controls(),
            |_| Ok(f64::NAN),
            || false,
        ),
        Err(GasEvaporationError::Invalid(_))
    ));
    let mut calls = 0;
    let cancelled = advance_surface_with_gas_and_heat_mass_selection(
        &Phase,
        &binding(0.0),
        original_liquid,
        &original_gas,
        &source,
        temperature,
        0.02,
        60.0,
        0.5,
        controls(),
        heat_controls(),
        |requested| {
            calls += 1;
            Ok(requested)
        },
        || true,
    );
    assert_eq!(cancelled, Err(GasEvaporationError::Cancelled));
    assert_eq!(calls, 0);
    assert_eq!(original_gas, before);
    assert_eq!(original_liquid, liquid(300.0));
}
