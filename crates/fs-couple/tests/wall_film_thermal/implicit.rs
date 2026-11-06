//! Independent physical endpoints and atomic admission on the joint solve.
use super::*;
use fs_couple::finite_gas_thermal::wall_film::implicit::{
    WallFilmImplicitControls, WallFilmImplicitError, WallFilmImplicitReport,
};
use fs_material::evaporation::gas::GasEvaporationError;
use fs_material::phase::{UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

fn balances(
    input: UniformEnthalpyStepInput<'_>,
    gas: &ClosedGasParcel,
    film: GasWallFilmInventory,
    binding: &fs_material::evaporation::gas::GasLiquidBinding,
    next: &UniformEnthalpyStep<WallFilmImplicitReport>,
) {
    let report = &next.report;
    let wall_change = input.mass_kg
        * (next.state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
    let gas_change = report.next_gas.mass_kg()
        * report.next_gas.state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let film_change = film_energy(report.next_film, binding) - film_energy(film, binding);
    near(
        wall_change + gas_change + film_change,
        input.internal_heat_j + report.external_gas_heat_j + report.radiation_heat_j,
        1e-8,
    );
    near(
        wall_change,
        input.internal_heat_j + next.external_heat_j,
        1e-8,
    );
    near(
        report.next_gas.mass_kg() + report.next_film.mass_kg(),
        gas.mass_kg() + film.mass_kg(),
        1e-14,
    );
    near(
        report.next_gas.species_mass_kg(&water()).unwrap() + report.next_film.mass_kg(),
        gas.species_mass_kg(&water()).unwrap() + film.mass_kg(),
        1e-14,
    );
    near(
        report.next_gas.species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-18,
    );
    assert_eq!(report.next_gas.volume_m3(), gas.volume_m3());
    assert_eq!(
        report.next_gas.point().axes().get("wall-film-lot"),
        Some(&13.0)
    );
    assert!(report.sum_absolute_energy_residual_j <= controls().wall.energy_tolerance_j);
    assert!(report.maximum_mass_residual_kg <= controls().phase.gas.mass_tolerance_kg);
    assert!(report.sum_absolute_mass_residual_kg <= controls().phase.gas.mass_tolerance_kg);
}

#[test]
fn dry_condensation_uses_actual_final_gas_for_flux_and_incoming_vapor_capacity() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let old_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let input = input(&curve, 300.0, 0.1);
    let film = GasWallFilmInventory::Dry;
    let environment = environment(&gas, &source, film, &binding, &body);
    let next = scope(|cx, _| {
        environment
            .advance_implicit(cx, input, WallFilmImplicitControls::default())
            .unwrap()
    });
    let report = &next.report;
    let wall_t = next.state.temperature_k();
    let gas_t = report.next_gas.state().state().temperature;
    let condensed = report.next_film.mass_kg();
    assert!(condensed > 0.0 && report.phase.film_formed);
    assert!(wall_t > 300.0 && gas_t < 320.0);
    assert_eq!(report.solid_to_gas_heat_j, 0.0);
    assert_eq!(report.phase.transferred_mass_kg, -condensed);
    // The final gas, not its old 320 K value, supplies incoming vapor Cp.
    let expected_wall = 300.0
        + condensed * (LATENT + 3.5 * R * (gas_t - 300.0))
            / (input.mass_kg * CP_WALL + condensed * CP_LIQUID);
    near(wall_t, expected_wall, 2e-8);
    let pressure = report.next_gas.state().state().pressure;
    let vapor = report.next_gas.species_mass_kg(&water()).unwrap() * R * gas_t / GAS_VOLUME;
    let saturation = PHASE.properties(wall_t).unwrap().saturation_pressure_pa;
    let endpoint_flux = environment.boundary.mass_transfer_coefficient_m_s * pressure / (R * gas_t)
        * ((pressure - vapor) / (pressure - saturation)).ln();
    near(
        report.phase.transferred_mass_kg,
        input.duration_s * AREA * endpoint_flux,
        1e-14,
    );
    near(
        liquid_t(report.next_film.specific_enthalpy_j_kg().unwrap()).unwrap(),
        wall_t,
        1e-7,
    );
    balances(input, &gas, film, &binding, &next);
    assert_eq!(gas, old_gas);
}

#[test]
fn complete_dryout_keeps_vapor_capacity_and_reduces_dry_limit_to_original_heat_owner() {
    let source = source();
    let gas = gas(&source, 290.0, 0.0);
    let binding = binding(&source);
    let curve = curve();
    let body = body(20.0);
    let input = input(&curve, 300.0, 0.1);
    let dm = 1e-7;
    let film = GasWallFilmInventory::Wet {
        mass_kg: dm,
        specific_enthalpy_j_kg: liquid_h(300.0).unwrap(),
    };
    let environment = environment(&gas, &source, film, &binding, &body);
    let next = scope(|cx, _| {
        environment
            .advance_implicit(cx, input, WallFilmImplicitControls::default())
            .unwrap()
    });
    let report = &next.report;
    assert_eq!(report.next_film, GasWallFilmInventory::Dry);
    assert_eq!(report.phase.transferred_mass_kg, dm);
    assert_eq!(report.phase.duration_s, input.duration_s);
    let exhaustion = report.phase.dry_out_time_s.unwrap();
    assert!(exhaustion > 0.0 && exhaustion < input.duration_s);
    let cw = input.mass_kg * CP_WALL;
    let cg = 2.5 * 100000.0 * GAS_VOLUME / 290.0;
    let exposure = body.convection_w_per_m2_k * AREA * input.duration_s;
    let a = cw + dm * 3.5 * R + exposure;
    let b = -exposure;
    let c = -(exposure + dm * 3.5 * R);
    let d = cg + dm * 2.5 * R + exposure;
    let rw = (cw + dm * 3.5 * R) * 300.0 - dm * LATENT;
    let rg = cg * 290.0;
    let tw = (rw * d - b * rg) / (a * d - b * c);
    let tg = (a * rg - c * rw) / (a * d - b * c);
    near(next.state.temperature_k(), tw, 1e-7);
    near(report.next_gas.state().state().temperature, tg, 1e-7);
    near(report.solid_to_gas_heat_j, exposure * (tw - tg), 2e-9);
    balances(input, &gas, film, &binding, &next);

    let mut dry_environment = environment.clone();
    dry_environment.film = GasWallFilmInventory::Dry;
    let dry = scope(|cx, _| {
        dry_environment
            .advance_implicit(cx, input, WallFilmImplicitControls::default())
            .unwrap()
    });
    let old = scope(|cx, _| dry_environment.advance(cx, input).unwrap());
    assert_eq!(dry.report.next_film, GasWallFilmInventory::Dry);
    assert_eq!(dry.report.phase.transferred_mass_kg, 0.0);
    near(dry.state.temperature_k(), old.state.temperature_k(), 1e-7);
    near(
        dry.report.next_gas.state().state().temperature,
        old.report.next_gas.state().state().temperature,
        1e-7,
    );
    balances(input, &gas, GasWallFilmInventory::Dry, &binding, &dry);
}

#[test]
fn work_final_cancellation_original_source_and_biot_refusals_are_atomic() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let old_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let input = input(&curve, 300.0, 0.1);
    let initial = input.initial;
    let base = environment(&gas, &source, GasWallFilmInventory::Dry, &binding, &body);
    let queries = Cell::new(0usize);
    let counted_h = |t| {
        queries.set(queries.get() + 1);
        liquid_h(t)
    };
    let mut environment = base.clone();
    environment.boundary.liquid_enthalpy_at_temperature = &counted_h;
    let accepted = scope(|cx, _| {
        environment
            .advance_implicit(cx, input, WallFilmImplicitControls::default())
            .unwrap()
    });
    let complete_queries = queries.get();
    assert!(complete_queries > 4);
    for limit in [1, accepted.report.residual_evaluations - 1] {
        let controls = WallFilmImplicitControls {
            maximum_evaluations: limit,
            ..WallFilmImplicitControls::default()
        };
        assert!(matches!(
            scope(|cx, _| base.advance_implicit(cx, input, controls)),
            Err(WallFilmImplicitError::EvaluationBudget)
        ));
    }
    scope(|cx, gate| {
        queries.set(0);
        let cancel_h = |t| {
            queries.set(queries.get() + 1);
            if queries.get() == complete_queries {
                gate.request();
            }
            liquid_h(t)
        };
        let mut environment = base.clone();
        environment.boundary.liquid_enthalpy_at_temperature = &cancel_h;
        assert!(matches!(
            environment.advance_implicit(cx, input, WallFilmImplicitControls::default()),
            Err(WallFilmImplicitError::Physical(
                WallFilmThermalError::Cancelled
            ))
        ));
        assert_eq!(queries.get(), complete_queries);
    });
    let refuse_h = |_: f64| {
        Err(EvaporationError::Invalid(
            "original liquid admission refused",
        ))
    };
    let mut refused = base.clone();
    refused.boundary.liquid_enthalpy_at_temperature = &refuse_h;
    assert!(matches!(
        scope(|cx, _| refused.advance_implicit(cx, input, WallFilmImplicitControls::default())),
        Err(WallFilmImplicitError::Physical(
            WallFilmThermalError::Phase(GasEvaporationError::Surface(EvaporationError::Invalid(
                "original liquid admission refused"
            )))
        ))
    ));
    let mut invalid_body = body.clone();
    invalid_body.convection_w_per_m2_k = 1e9;
    refused = base.clone();
    refused.body = &invalid_body;
    assert!(matches!(
        scope(|cx, _| refused.advance_implicit(cx, input, WallFilmImplicitControls::default())),
        Err(WallFilmImplicitError::Physical(
            WallFilmThermalError::Thermal(_)
        ))
    ));
    assert_eq!(gas, old_gas);
    assert_eq!(base.film, GasWallFilmInventory::Dry);
    assert_eq!(input.initial, initial);
    let mut tiny_body = body.clone();
    tiny_body.convection_w_per_m2_k = f64::from_bits(1);
    let mut tiny = base.clone();
    tiny.body = &tiny_body;
    assert!(matches!(
        scope(|cx, _| tiny.advance_implicit(cx, input, WallFilmImplicitControls::default())),
        Err(WallFilmImplicitError::Physical(
            WallFilmThermalError::Invalid("nonzero wall convection exposure underflows")
        ))
    ));
}
