use super::*;
use crate::evaporation::gas::tests::{Phase, binding, catalog, controls, gas, liquid, temperature};
use fs_thermochem::SpeciesId;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn solve() -> GasEvaporationSolve {
    let mut solve = controls();
    solve.maximum_temperature_increment_k = 200.0;
    solve
}
fn surface<'a>(
    binding: &'a GasLiquidBinding,
    film: GasWallFilmInventory,
    temperature_k: f64,
) -> GasWallFilmSurface<'a> {
    GasWallFilmSurface {
        law: &Phase,
        binding,
        film,
        area_m2: 0.001,
        interface_liquid_enthalpy_j_kg: liquid(temperature_k).specific_enthalpy_j_kg,
        temperature_from_enthalpy: &temperature,
        mass_transfer_coefficient_m_s: 0.02,
    }
}
fn wet(mass_kg: f64, temperature_k: f64) -> GasWallFilmInventory {
    GasWallFilmInventory::Wet {
        mass_kg,
        specific_enthalpy_j_kg: liquid(temperature_k).specific_enthalpy_j_kg,
    }
}
fn vapor() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
fn carrier() -> SpeciesId {
    SpeciesId::new("A").unwrap()
}

#[test]
fn dry_supersaturated_wall_grows_real_film_and_receives_upstream_condensation_heat() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 320.0, 5000.0, 0.01);
    let input = surface(&binding, GasWallFilmInventory::Dry, 300.0);
    let candidate =
        advance_wall_film_with_gas(input, &initial, &source, 0.5, solve(), || false).unwrap();
    let step = candidate.step();
    assert!(step.transferred_mass_kg < 0.0);
    assert!(step.film_formed);
    assert_eq!(step.dry_out_time_s, None);
    assert_eq!(step.phase_active_duration_s, 0.5);
    assert_eq!(candidate.film().mass_kg(), -step.transferred_mass_kg);
    assert_eq!(
        candidate.film().specific_enthalpy_j_kg(),
        Some(input.interface_liquid_enthalpy_j_kg)
    );
    assert_eq!(step.film_sensible_energy_j, 0.0);
    assert!(step.wall_heat_into_system_j < 0.0);
    // Manufactured vapor Cp/R=4+0.001*T and L=200 kJ/kg. A cooling wall
    // receives both latent heat and the upstream vapor's sensible energy.
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let upstream_sensible = r * (4.0 * 20.0 + 0.5 * 0.001 * (320.0 * 320.0 - 300.0 * 300.0));
    near(
        step.wall_heat_into_system_j,
        step.transferred_mass_kg * (200000.0 + upstream_sensible),
        1e-10,
    );
    near(
        step.film_enthalpy_change_j + step.gas_internal_energy_change_j
            - step.wall_heat_into_system_j,
        0.0,
        solve().gas.energy_tolerance_j,
    );
    near(
        initial.species_mass_kg(&vapor()).unwrap(),
        candidate.gas().species_mass_kg(&vapor()).unwrap() + candidate.film().mass_kg(),
        solve().gas.mass_tolerance_kg,
    );
    near(
        initial.species_mass_kg(&carrier()).unwrap(),
        candidate.gas().species_mass_kg(&carrier()).unwrap(),
        solve().gas.mass_tolerance_kg,
    );
    assert!(step.sum_absolute_energy_residual_j <= solve().gas.energy_tolerance_j);
    assert!(std::ptr::eq(candidate.catalog(), &source));
}

#[test]
fn complete_evaporation_has_exact_frozen_rate_event_and_stays_dry_afterward() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 320.0, 500.0, 0.01);
    let initial_film = wet(1e-7, 300.0);
    let input = surface(&binding, initial_film, 340.0);
    let wall_t = temperature(input.interface_liquid_enthalpy_j_kg).unwrap();
    let rate = input.area_m2
        * surface_flux(
            &Phase,
            wall_t,
            GasFilm {
                total_pressure_pa: initial.state().state().pressure,
                vapor_partial_pressure_pa: initial.species_partial_pressure_pa(&vapor()).unwrap(),
                temperature_k: initial.state().state().temperature,
                mass_transfer_coefficient_m_s: input.mass_transfer_coefficient_m_s,
            },
        )
        .unwrap()
        .mass_flux_kg_m2_s;
    assert!(rate > 0.0);
    let event = initial_film.mass_kg() / rate;
    let duration = 2.0 * event;
    let candidate =
        advance_wall_film_with_gas(input, &initial, &source, duration, solve(), || false).unwrap();
    let step = candidate.step();
    assert_eq!(step.duration_s, duration);
    assert_eq!(step.phase_active_duration_s, event);
    assert_eq!(step.dry_out_time_s, Some(event));
    assert_eq!(step.transferred_mass_kg, initial_film.mass_kg());
    assert!(step.requested_transferred_mass_kg > step.transferred_mass_kg);
    assert_eq!(candidate.film(), GasWallFilmInventory::Dry);
    assert_eq!(candidate.film().specific_enthalpy_j_kg(), None);
    assert_eq!(step.liquid_mass_kg, [1e-7, 0.0]);
    assert_eq!(
        step.film_enthalpy_change_j,
        -1e-7 * initial_film.specific_enthalpy_j_kg().unwrap()
    );
    let sensible = 1e-7
        * (input.interface_liquid_enthalpy_j_kg - initial_film.specific_enthalpy_j_kg().unwrap());
    near(
        step.wall_heat_into_system_j,
        sensible + 1e-7 * 200000.0,
        1e-10,
    );
    near(
        step.film_enthalpy_change_j + step.gas_internal_energy_change_j
            - step.wall_heat_into_system_j,
        0.0,
        solve().gas.energy_tolerance_j,
    );
    let partial =
        advance_wall_film_with_gas(input, &initial, &source, 0.5 * event, solve(), || false)
            .unwrap();
    assert!(matches!(partial.film(), GasWallFilmInventory::Wet { .. }));
    assert_eq!(partial.step().dry_out_time_s, None);
    near(partial.film().mass_kg(), 0.5e-7, 1e-18);
    let idle = advance_wall_film_with_gas(
        GasWallFilmSurface {
            film: candidate.film(),
            ..input
        },
        candidate.gas(),
        &source,
        duration,
        solve(),
        || false,
    )
    .unwrap();
    assert_eq!(idle.film(), GasWallFilmInventory::Dry);
    assert_eq!(idle.gas(), candidate.gas());
    assert!(idle.step().requested_transferred_mass_kg > 0.0);
    assert_eq!(idle.step().transferred_mass_kg, 0.0);
    assert_eq!(idle.step().phase_active_duration_s, 0.0);
    assert_eq!(idle.step().wall_heat_into_system_j, 0.0);
}

#[test]
fn existing_film_equilibrates_through_explicit_wall_port_and_zero_time_is_identity() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 320.0, 5000.0, 0.01);
    let initial_film = wet(0.01, 300.0);
    let mut input = surface(&binding, initial_film, 310.0);
    input.mass_transfer_coefficient_m_s = 0.0;
    let expected_heat = initial_film.mass_kg()
        * (input.interface_liquid_enthalpy_j_kg - initial_film.specific_enthalpy_j_kg().unwrap());
    let candidate =
        advance_wall_film_with_gas(input, &initial, &source, 1.0, solve(), || false).unwrap();
    assert_eq!(candidate.gas(), &initial);
    assert_eq!(candidate.film().mass_kg(), initial_film.mass_kg());
    assert_eq!(
        candidate.film().specific_enthalpy_j_kg(),
        Some(input.interface_liquid_enthalpy_j_kg)
    );
    assert_eq!(candidate.step().film_sensible_energy_j, expected_heat);
    assert_eq!(candidate.step().wall_heat_into_system_j, expected_heat);
    assert_eq!(candidate.step().film_enthalpy_change_j, expected_heat);
    assert_eq!(candidate.step().transferred_mass_kg, 0.0);
    let zero =
        advance_wall_film_with_gas(input, &initial, &source, 0.0, solve(), || false).unwrap();
    assert_eq!(zero.film(), initial_film);
    assert_eq!(zero.gas(), &initial);
    assert_eq!(zero.step().wall_heat_into_system_j, 0.0);
    assert_eq!(zero.step().film_enthalpy_change_j, 0.0);
    let dry = advance_wall_film_with_gas(
        GasWallFilmSurface {
            film: GasWallFilmInventory::Dry,
            ..input
        },
        &initial,
        &source,
        1.0,
        solve(),
        || false,
    )
    .unwrap();
    assert_eq!(dry.film(), GasWallFilmInventory::Dry);
    assert_eq!(dry.step().wall_heat_into_system_j, 0.0);
}

#[test]
fn fixed_caloric_datum_changes_reports_without_changing_dryout_or_wall_heat() {
    let mut physical = Vec::new();
    for reference in [0.0, -1e8, 1e12] {
        let source = catalog(reference);
        let binding = binding(reference);
        let initial = gas(&source, 320.0, 500.0, 0.01);
        let mut controls = solve();
        controls.latent_compatibility_tolerance_j_kg = 0.01;
        let result = advance_wall_film_with_gas(
            surface(&binding, wet(1e-7, 300.0), 340.0),
            &initial,
            &source,
            2.0,
            controls,
            || false,
        )
        .unwrap();
        assert_eq!(result.film(), GasWallFilmInventory::Dry);
        assert!(result.step().sum_absolute_energy_residual_j <= controls.gas.energy_tolerance_j);
        physical.push((
            result.step().dry_out_time_s.unwrap(),
            result.step().wall_heat_into_system_j,
            result.gas().state().state().temperature,
        ));
    }
    for value in &physical[1..] {
        near(value.0, physical[0].0, 0.0);
        near(value.1, physical[0].1, 1e-8);
        near(value.2, physical[0].2, 1e-9);
    }
}

#[test]
fn invalid_wet_inventory_vapor_overdraw_and_late_cancellation_publish_nothing() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 320.0, 5000.0, 0.01);
    let before = initial.clone();
    for mass in [0.0, -1.0, f64::INFINITY] {
        assert!(
            advance_wall_film_with_gas(
                surface(&binding, wet(mass, 300.0), 300.0),
                &initial,
                &source,
                0.2,
                solve(),
                || false,
            )
            .is_err()
        );
    }
    let tiny = gas(&source, 320.0, 5000.0, 1e-6);
    let error = advance_wall_film_with_gas(
        surface(&binding, GasWallFilmInventory::Dry, 300.0),
        &tiny,
        &source,
        1.0,
        solve(),
        || false,
    )
    .unwrap_err();
    assert!(
        matches!(
            error,
            GasEvaporationError::StepSizeLimit {
                reason: "initial wall-film vapor inventory withdrawal limit",
                ..
            }
        ),
        "{error:?}"
    );
    let input = surface(&binding, GasWallFilmInventory::Dry, 300.0);
    let mut calls = 0;
    advance_wall_film_with_gas(input, &initial, &source, 0.2, solve(), || {
        calls += 1;
        false
    })
    .unwrap();
    let final_gate = calls;
    let mut calls = 0;
    assert_eq!(
        advance_wall_film_with_gas(input, &initial, &source, 0.2, solve(), || {
            calls += 1;
            calls == final_gate
        })
        .unwrap_err(),
        GasEvaporationError::Cancelled
    );
    assert_eq!(initial, before);
    assert_eq!(input.film, GasWallFilmInventory::Dry);
}
