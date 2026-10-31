use super::*;
use crate::evaporation::gas::tests::{Phase, binding, catalog, controls, gas, liquid, temperature};
use fs_thermal::evaporation::{GasFilm, surface_flux};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}; tolerance {tolerance}"
    );
}
fn surface<'a>(
    id: u64,
    binding: &'a GasLiquidBinding,
    liquid: LiquidSurfaceState,
) -> GasLiquidChamberSurface<'a> {
    GasLiquidChamberSurface {
        id,
        law: &Phase,
        binding,
        liquid,
        temperature_from_enthalpy: &temperature,
        mass_transfer_coefficient_m_s: 0.02,
    }
}
fn solve() -> GasLiquidChamberSolve {
    GasLiquidChamberSolve {
        surface: controls(),
        maximum_surfaces: 8,
    }
}
fn vapor() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
fn carrier() -> SpeciesId {
    SpeciesId::new("A").unwrap()
}

#[test]
fn one_surface_preserves_original_gas_liquid_and_physical_accounts_and_empty_admits() {
    let source = catalog(-1e7);
    let binding = binding(-1e7);
    let initial = gas(&source, 320.0, 5000.0, 0.01);
    let liquid = liquid(300.0);
    let old = super::super::advance_surface_with_gas(
        &Phase,
        &binding,
        liquid,
        &initial,
        &source,
        temperature,
        0.02,
        0.5,
        controls(),
        || false,
    )
    .unwrap();
    let joint = advance_chamber_with_gas(
        &[surface(9, &binding, liquid)],
        &initial,
        &source,
        0.5,
        solve(),
        || false,
    )
    .unwrap();
    let local = &joint.surfaces()[0];
    assert_eq!(local.id, 9);
    assert_eq!(local.liquid, old.liquid());
    assert_eq!(joint.gas_exchange(), old.gas_exchange());
    assert_eq!(local.transferred_mass_kg, old.step().transferred_mass_kg);
    assert_eq!(local.latent_energy_j, old.step().latent_energy_j);
    assert_eq!(local.energy_into_gas_j, old.step().energy_into_gas_j);
    assert_eq!(
        local.liquid_enthalpy_change_j,
        old.step().liquid_enthalpy_change_j
    );
    assert_eq!(
        local.surface_temperature_k,
        old.step().surface_temperature_k
    );
    assert_eq!(
        joint.step().gas_internal_energy_change_j,
        old.step().gas_internal_energy_change_j
    );
    near(
        joint.step().energy_residual_j,
        old.step().energy_residual_j,
        0.0,
    );
    let empty = advance_chamber_with_gas(&[], &initial, &source, 0.5, solve(), || false).unwrap();
    assert_eq!(empty.gas(), &initial);
    assert!(empty.surfaces().is_empty());
    assert!(empty.step().species.is_empty());
    assert_eq!(empty.step().transferred_mass_kg, 0.0);
    assert_eq!(empty.step().sum_absolute_energy_residual_j, 0.0);
    assert!(std::ptr::eq(empty.catalog(), &source));
}

#[test]
fn opposing_equal_species_streams_keep_heat_and_are_order_independent_with_native_selection() {
    let mut physical = Vec::new();
    for reference in [0.0, 1e12] {
        let source = catalog(reference);
        let binding = binding(reference);
        let initial = gas(&source, 320.0, 1600.0, 0.01);
        let hot = liquid(360.0);
        let mut cold = liquid(280.0);
        let duration = 0.2;
        let film = GasFilm {
            total_pressure_pa: initial.state().state().pressure,
            vapor_partial_pressure_pa: initial.species_partial_pressure_pa(&vapor()).unwrap(),
            temperature_k: 320.0,
            mass_transfer_coefficient_m_s: 0.02,
        };
        let target =
            (hot.area_m2 * duration) * surface_flux(&Phase, 360.0, film).unwrap().mass_flux_kg_m2_s;
        cold.area_m2 =
            target / (-surface_flux(&Phase, 280.0, film).unwrap().mass_flux_kg_m2_s * duration);
        let cold_surface = surface(3, &binding, cold);
        let hot_surface = surface(8, &binding, hot);
        let mut config = solve();
        config.surface.latent_compatibility_tolerance_j_kg = 0.01;
        let mut selected_ids = Vec::new();
        let selected = advance_chamber_with_gas_mass_selection(
            &[hot_surface, cold_surface],
            &initial,
            &source,
            duration,
            config,
            |id, requested| {
                selected_ids.push(id);
                let represented = if id == 3 { -target } else { target };
                near(represented, requested, config.surface.gas.mass_tolerance_kg);
                Ok(represented)
            },
            || false,
        )
        .unwrap();
        assert_eq!(selected_ids, [3, 8]);
        let reordered = advance_chamber_with_gas_mass_selection(
            &[cold_surface, hot_surface],
            &initial,
            &source,
            duration,
            config,
            |id, _| Ok(if id == 3 { -target } else { target }),
            || false,
        )
        .unwrap();
        assert_eq!(selected, reordered);
        assert_eq!(selected.step().transferred_mass_kg, 0.0);
        assert_eq!(selected.step().species.len(), 1);
        assert_eq!(selected.step().species[0].transferred_mass_kg, 0.0);
        assert!(selected.step().evaporated_mass_kg > 0.0);
        assert_eq!(
            selected.step().evaporated_mass_kg,
            selected.step().condensed_mass_kg
        );
        assert!(selected.step().energy_into_gas_j > 0.0);
        assert!(selected.step().gas_temperature_k[1] > 320.0);
        assert!(selected.surfaces()[0].liquid.mass_kg > cold.mass_kg);
        assert!(selected.surfaces()[1].liquid.mass_kg < hot.mass_kg);
        near(
            selected.gas().species_mass_kg(&carrier()).unwrap(),
            initial.species_mass_kg(&carrier()).unwrap(),
            config.surface.gas.mass_tolerance_kg,
        );
        near(
            selected.step().mass_residual_kg,
            0.0,
            config.surface.gas.mass_tolerance_kg,
        );
        assert!(
            selected.step().sum_absolute_energy_residual_j <= config.surface.gas.energy_tolerance_j
        );
        physical.push((
            selected.gas().state().state().temperature,
            selected.surfaces()[0].liquid,
            selected.surfaces()[1].liquid,
        ));
    }
    near(physical[0].0, physical[1].0, 1e-9);
    for (a, b) in [
        (physical[0].1, physical[1].1),
        (physical[0].2, physical[1].2),
    ] {
        near(a.mass_kg, b.mass_kg, 0.0);
        near(a.specific_enthalpy_j_kg, b.specific_enthalpy_j_kg, 1e-6);
    }
}

#[test]
fn gross_shared_vapor_is_checked_before_one_joint_gas_reconstruction() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 300.0, 5000.0, 1e-6);
    let liquid = liquid(300.0);
    let old_vapor = initial.species_mass_kg(&vapor()).unwrap();
    let flux = surface_flux(
        &Phase,
        300.0,
        GasFilm {
            total_pressure_pa: 100000.0,
            vapor_partial_pressure_pa: 5000.0,
            temperature_k: 300.0,
            mass_transfer_coefficient_m_s: 0.02,
        },
    )
    .unwrap()
    .mass_flux_kg_m2_s;
    let duration = 0.2 * old_vapor / (-flux * liquid.area_m2);
    // Each film alone is inside the 25% donor allowance and remains condensing.
    advance_chamber_with_gas(
        &[surface(1, &binding, liquid)],
        &initial,
        &source,
        duration,
        solve(),
        || false,
    )
    .unwrap();
    let error = advance_chamber_with_gas(
        &[surface(2, &binding, liquid), surface(1, &binding, liquid)],
        &initial,
        &source,
        duration,
        solve(),
        || false,
    )
    .unwrap_err();
    assert!(
        matches!(
            error,
            GasEvaporationError::StepSizeLimit {
                reason: "shared initial vapor inventory withdrawal limit",
                ..
            }
        ),
        "{error:?}"
    );
    assert_eq!(initial.species_mass_kg(&vapor()).unwrap(), old_vapor);
}

#[test]
fn all_films_check_equilibrium_against_the_common_endpoint() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 300.0, 1700.0, 1e-4);
    let liquid = liquid(300.0);
    let old_vapor = initial.species_mass_kg(&vapor()).unwrap();
    let flux = surface_flux(
        &Phase,
        300.0,
        GasFilm {
            total_pressure_pa: 100000.0,
            vapor_partial_pressure_pa: 1700.0,
            temperature_k: 300.0,
            mass_transfer_coefficient_m_s: 0.02,
        },
    )
    .unwrap()
    .mass_flux_kg_m2_s;
    let duration = 0.08 * old_vapor / (-flux * liquid.area_m2);
    advance_chamber_with_gas(
        &[surface(1, &binding, liquid)],
        &initial,
        &source,
        duration,
        solve(),
        || false,
    )
    .unwrap();
    let error = advance_chamber_with_gas(
        &[surface(1, &binding, liquid), surface(2, &binding, liquid)],
        &initial,
        &source,
        duration,
        solve(),
        || false,
    )
    .unwrap_err();
    assert!(
        matches!(
            error,
            GasEvaporationError::StepSizeLimit {
                reason: "finite-atmosphere equilibrium overshoot",
                ..
            }
        ),
        "{error:?}"
    );
}

#[test]
fn duplicate_work_selection_and_late_cancellation_refuse_without_partial_publication() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 320.0, 5000.0, 0.01);
    let before = initial.clone();
    let first = surface(1, &binding, liquid(300.0));
    let second = surface(2, &binding, liquid(305.0));
    let mut selections = 0;
    assert!(
        advance_chamber_with_gas_mass_selection(
            &[first, first],
            &initial,
            &source,
            0.2,
            solve(),
            |_, dm| {
                selections += 1;
                Ok(dm)
            },
            || false,
        )
        .is_err()
    );
    assert_eq!(selections, 0);
    assert!(
        advance_chamber_with_gas(
            &[first, second],
            &initial,
            &source,
            0.2,
            GasLiquidChamberSolve {
                maximum_surfaces: 1,
                ..solve()
            },
            || false,
        )
        .is_err()
    );
    let requested_error = GasEvaporationError::Invalid("native inventory refused");
    assert_eq!(
        advance_chamber_with_gas_mass_selection(
            &[first, second],
            &initial,
            &source,
            0.2,
            solve(),
            |_, _| Err(requested_error.clone()),
            || false,
        )
        .unwrap_err(),
        requested_error,
    );
    let mut calls = 0;
    advance_chamber_with_gas(&[first, second], &initial, &source, 0.2, solve(), || {
        calls += 1;
        false
    })
    .unwrap();
    let final_gate = calls;
    let mut calls = 0;
    assert_eq!(
        advance_chamber_with_gas(&[first, second], &initial, &source, 0.2, solve(), || {
            calls += 1;
            calls == final_gate
        },)
        .unwrap_err(),
        GasEvaporationError::Cancelled,
    );
    assert_eq!(initial, before);
}
