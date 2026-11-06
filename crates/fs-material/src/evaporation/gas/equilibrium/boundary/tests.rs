//! G1 original mass/enthalpy event balances and G3 legacy/scale comparisons.
use super::*;
use crate::evaporation::gas::{
    equilibrium::{
        advance_equilibrium_surface,
        tests::{CP, FUSION, MELTING, TRANSFER, VAPORIZATION, bulk, catalog, curve, near},
    },
    tests::{Phase, binding, controls, gas},
};
use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_blake3::ContentHash;
use fs_thermal::evaporation::{EvaporationError, EvaporationProperties};

fn event_controls() -> GasEvaporationSolve {
    let mut result = controls();
    // This fixture freezes after evaporating one third of its original liquid.
    // The event must fit this explicit donor budget; the owner cannot waive it.
    result.maximum_inventory_fraction = 0.5;
    result
}

#[test]
fn complete_freezing_conserves_the_original_system_and_returns_unconsumed_time() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let gas = gas(&source, MELTING, 0.0, 1.0);
    let maximum_duration = 20_000.0;
    // Spending this entire requested interval would exceed the donor budget.
    assert!(matches!(
        advance_equilibrium_surface(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &source,
            TRANSFER,
            maximum_duration,
            event_controls(),
            || false,
        ),
        Err(GasEvaporationError::StepSizeLimit { .. })
    ));
    let result = advance_to_freezing(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &gas,
        &source,
        TRANSFER,
        maximum_duration,
        event_controls(),
        || false,
    )
    .unwrap();
    let boundary = result.boundary();
    let exchange = result.exchange();
    let step = exchange.step();
    // Independent isothermal constituent balance and integrated Stefan film.
    let dm = 0.005 * FUSION / (FUSION + VAPORIZATION);
    let mass_rate = initial.area_m2 * TRANSFER * 100_000.0 / ((CP / 4.0) * MELTING)
        * (100_000.0_f64 / (100_000.0 - 1500.0)).ln();
    near(step.transferred_mass_kg, dm, 1e-15);
    near(boundary.elapsed_duration_s, dm / mass_rate, 1e-8);
    assert!(boundary.freezing_completed);
    assert!(boundary.remaining_duration_s > 0.0);
    near(
        boundary.elapsed_duration_s + boundary.remaining_duration_s,
        maximum_duration,
        1e-11,
    );
    assert_eq!(step.duration_s, boundary.elapsed_duration_s);
    assert_eq!(
        exchange.liquid().specific_enthalpy_j_kg,
        chart.knots()[1].specific_enthalpy_j_kg
    );
    assert_eq!(result.phase().final_state.liquid_mass_kg, 0.0);
    assert_eq!(
        result.phase().final_state.solid_mass_kg,
        exchange.liquid().mass_kg
    );
    assert_eq!(step.surface_temperature_k, [MELTING, MELTING]);
    near(result.phase().freezing_mass_kg, 0.005 - dm, 1e-15);
    near(result.phase().fusion_energy_j, dm * VAPORIZATION, 1e-9);
    near(
        result.phase().liquid_partition_residual_kg,
        0.0,
        controls().gas.mass_tolerance_kg,
    );
    near(
        step.transferred_mass_kg,
        step.requested_transferred_mass_kg,
        controls().gas.mass_tolerance_kg,
    );
    let total_before = initial.mass_kg * initial.specific_enthalpy_j_kg
        + gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let total_after = exchange.liquid().mass_kg * exchange.liquid().specific_enthalpy_j_kg
        + exchange.gas().mass_kg() * exchange.gas().state().specific_internal_energy_j_kg();
    near(total_after, total_before, 1e-7);
    assert!(
        boundary.condensed_energy_residual_j.abs()
            + exchange.gas_exchange().step().energy_residual_j.abs()
            <= controls().gas.energy_tolerance_j
    );
    assert!(inventory.state(exchange.liquid()).is_err());
    assert_eq!(
        result
            .phase()
            .final_state
            .equilibrium
            .phase_curve_identity(),
        chart.identity()
    );
    assert!(std::ptr::eq(exchange.catalog(), &source));
}

#[test]
fn no_event_and_inactive_intervals_match_the_existing_scalar_owner_exactly() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    for (temperature, pressure, coefficient, duration) in [
        (300.0, 0.0, TRANSFER, 1.0),
        (320.0, 5000.0, TRANSFER, 0.5),
        (300.0, 0.0, 0.0, 20_000.0),
        (300.0, 0.0, TRANSFER, 0.0),
    ] {
        let gas = gas(&source, temperature, pressure, 0.01);
        let old = advance_equilibrium_surface(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &source,
            coefficient,
            duration,
            controls(),
            || false,
        )
        .unwrap();
        let new = advance_to_freezing(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &source,
            coefficient,
            duration,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(new.exchange(), old.exchange());
        assert_eq!(new.phase(), old.phase());
        assert!(!new.boundary().freezing_completed);
        assert_eq!(new.boundary().elapsed_duration_s, duration);
        assert_eq!(new.boundary().remaining_duration_s, 0.0);
    }
}

#[test]
fn finite_early_event_survives_overflow_of_the_unused_full_exposure() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let gas = gas(&source, MELTING, 0.0, 1.0);
    let normal = advance_to_freezing(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &gas,
        &source,
        TRANSFER,
        20_000.0,
        event_controls(),
        || false,
    )
    .unwrap();
    let scaled = LiquidSurfaceState {
        area_m2: 1e200,
        ..initial
    };
    let horizon = 1e200;
    assert!(!(scaled.area_m2 * horizon).is_finite());
    let result = advance_to_freezing(
        inventory,
        &Phase,
        &binding(0.0),
        scaled,
        &gas,
        &source,
        TRANSFER,
        horizon,
        event_controls(),
        || false,
    )
    .unwrap();
    assert!(result.boundary().freezing_completed);
    assert!(result.boundary().elapsed_duration_s > 0.0);
    assert!(result.boundary().elapsed_duration_s < 1e-190);
    assert_eq!(result.exchange().gas(), normal.exchange().gas());
    assert_eq!(result.phase(), normal.phase());
    near(
        result.boundary().elapsed_duration_s * scaled.area_m2,
        normal.boundary().elapsed_duration_s * initial.area_m2,
        1e-12,
    );
    near(
        result.exchange().step().transferred_mass_kg,
        result.exchange().step().requested_transferred_mass_kg,
        controls().gas.mass_tolerance_kg,
    );
}

#[test]
fn sufficiently_cold_admitted_condensation_can_complete_freezing_with_signed_transfer() {
    struct LowLatent;
    impl EvaporationLaw for LowLatent {
        fn temperature_range_k(&self) -> [f64; 2] {
            Phase.temperature_range_k()
        }
        fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
            let mut value = Phase.properties(t)?;
            value.latent_heat_j_kg = 40_000.0;
            value.saturation_pressure_pa =
                1500.0 * fs_math::det::exp(40_000.0 / (CP / 4.0) * (1.0 / MELTING - 1.0 / t));
            value.saturation_pressure_derivative_pa_k =
                value.saturation_pressure_pa * 40_000.0 / ((CP / 4.0) * t * t);
            Ok(value)
        }
    }
    let mut knots = curve().knots().to_vec();
    for knot in &mut knots {
        knot.specific_enthalpy_j_kg += VAPORIZATION - 40_000.0;
    }
    let chart = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x78; 32]), knots).unwrap();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let source = catalog(0.0);
    let gas = gas(&source, 220.0, 5000.0, 10.0);
    let result = advance_to_freezing(
        inventory,
        &LowLatent,
        &binding(0.0),
        initial,
        &gas,
        &source,
        TRANSFER,
        1e6,
        event_controls(),
        || false,
    )
    .unwrap();
    let solid_h = chart.knots()[1].specific_enthalpy_j_kg;
    let expected_dm =
        -initial.mass_kg * (initial.specific_enthalpy_j_kg - solid_h) / (solid_h - CP * 220.0);
    let step = result.exchange().step();
    assert!(result.boundary().freezing_completed && step.transferred_mass_kg < 0.0);
    near(step.transferred_mass_kg, expected_dm, 1e-12);
    assert_eq!(result.phase().final_state.liquid_mass_kg, 0.0);
    near(result.phase().freezing_mass_kg, 0.005 - expected_dm, 1e-12);
    assert_eq!(result.exchange().liquid().specific_enthalpy_j_kg, solid_h);
    near(step.energy_into_gas_j, expected_dm * CP * 220.0, 1e-6);
    near(
        step.liquid_enthalpy_change_j + step.gas_internal_energy_change_j,
        0.0,
        1e-7,
    );
    near(step.mass_residual_kg, 0.0, controls().gas.mass_tolerance_kg);
    assert!(step.vapor_partial_pressure_pa[1] > 1500.0);
}

#[test]
fn terminal_knot_is_the_original_value_not_a_reconstructed_latent_subtraction() {
    let solid_h = -0.1;
    let liquid_h = 100_000.0;
    assert_ne!(liquid_h - (liquid_h - solid_h), solid_h);
    let knot = |h, t, fraction: f64| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: fraction,
        bulk_density_kg_m3: 1.0 / ((1.0 - fraction) / 900.0 + fraction / 1000.0),
    };
    let chart = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x79; 32]),
        vec![
            knot(solid_h - CP * 100.0, 200.0, 0.0),
            knot(solid_h, MELTING, 0.0),
            knot(0.5 * (solid_h + liquid_h), MELTING, 0.5),
            knot(liquid_h, MELTING, 1.0),
            knot(liquid_h + CP * 600.0, 900.0, 1.0),
        ],
    )
    .unwrap();
    let source = catalog(0.0);
    let initial = LiquidSurfaceState {
        mass_kg: 0.01,
        area_m2: 1e-3,
        specific_enthalpy_j_kg: chart.knots()[2].specific_enthalpy_j_kg,
    };
    let mut binding = binding(0.0);
    binding.liquid_enthalpy_offset_j_kg = CP * MELTING - VAPORIZATION - liquid_h;
    let gas = gas(&source, MELTING, 0.0, 1.0);
    let result = advance_to_freezing(
        EquilibriumEvaporation::new(&chart).unwrap(),
        &Phase,
        &binding,
        initial,
        &gas,
        &source,
        TRANSFER,
        20_000.0,
        event_controls(),
        || false,
    )
    .unwrap();
    assert!(result.boundary().freezing_completed);
    assert_eq!(result.exchange().liquid().specific_enthalpy_j_kg, solid_h);
    assert_eq!(result.phase().final_state.liquid_mass_kg, 0.0);
    assert!(
        result.boundary().condensed_energy_residual_j.abs()
            + result
                .exchange()
                .gas_exchange()
                .step()
                .energy_residual_j
                .abs()
            <= controls().gas.energy_tolerance_j
    );
}

#[test]
fn terminal_events_keep_donor_temperature_saturation_and_cancellation_admission() {
    let source = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let original_gas = gas(&source, MELTING, 0.0, 1.0);
    let saved = original_gas.clone();
    let run = |gas: &ClosedGasParcel, budget, cancelled: &mut dyn FnMut() -> bool| {
        advance_to_freezing(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            gas,
            &source,
            TRANSFER,
            20_000.0,
            budget,
            cancelled,
        )
    };
    assert!(matches!(
        run(&original_gas, controls(), &mut || false),
        Err(GasEvaporationError::StepSizeLimit {
            reason: "initial liquid or vapor inventory withdrawal limit",
            ..
        })
    ));
    let mut temperature_budget = event_controls();
    temperature_budget.maximum_temperature_increment_k = 1e-6;
    assert!(matches!(
        run(&original_gas, temperature_budget, &mut || false),
        Err(GasEvaporationError::StepSizeLimit {
            reason: "gas temperature departure",
            ..
        })
    ));
    let small_gas = gas(&source, MELTING, 0.0, 0.01);
    let mut saturation_budget = event_controls();
    saturation_budget.maximum_temperature_increment_k = 100.0;
    assert!(matches!(
        run(&small_gas, saturation_budget, &mut || false),
        Err(GasEvaporationError::StepSizeLimit {
            reason: "finite-atmosphere equilibrium overshoot",
            ..
        })
    ));
    let mut calls = 0;
    run(&original_gas, event_controls(), &mut || {
        calls += 1;
        false
    })
    .unwrap();
    let mut attempt = 0;
    assert!(matches!(
        run(&original_gas, event_controls(), &mut || {
            attempt += 1;
            attempt == calls
        }),
        Err(GasEvaporationError::Cancelled)
    ));
    assert_eq!(original_gas, saved);
    assert_eq!(initial, bulk(&chart, 0.5));
}
