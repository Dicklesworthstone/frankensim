//! Manufactured compatible charts exercise the original molecular/gas owners.
use super::*;
use crate::sublimation::gas::{UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, tests as fixture};
use fs_thermochem::SpeciesId;

const SUBSTRATE_CP: f64 = 500.0;
const SUBSTRATE_MASS: f64 = 0.02;
const AREA: f64 = 1e-4;
const LATENT: f64 = 200000.0;

fn controls() -> GasSupportedChamberSolve {
    GasSupportedChamberSolve {
        chamber: GasChamberSolve {
            surface: fixture::controls(),
            maximum_surfaces: 16,
        },
        storage: SolverConfig {
            relative_tolerance: 1e-12,
            absolute_energy_tolerance_j: 1e-10,
            max_iterations: 128,
        },
    }
}
fn state(t: f64, deposit: f64, datum: f64) -> SupportedSurfaceState {
    SupportedSurfaceState {
        substrate_mass_kg: SUBSTRATE_MASS,
        substrate_specific_enthalpy_j_kg: datum + SUBSTRATE_CP * (t - 400.0),
        deposit_mass_kg: deposit,
        area_m2: AREA,
    }
}
fn substrate(h: f64) -> Result<f64, SublimationError> {
    Ok(400.0 + h / SUBSTRATE_CP)
}
fn surface<'a>(
    id: u64,
    law: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    state: SupportedSurfaceState,
    substrate_temperature: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    deposit_enthalpy: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasSupportedSurface<'a> {
    GasSupportedSurface {
        id,
        law,
        binding,
        state,
        substrate_temperature,
        deposit_enthalpy,
    }
}
fn r() -> f64 {
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015
}
fn incoming_rate(pressure: f64, temperature: f64, area: f64) -> f64 {
    0.1 * pressure / (2.0 * core::f64::consts::PI * r() * temperature).sqrt() * area
}

#[test]
fn bare_substrate_captures_finite_vapor_and_shares_real_sensible_storage() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let initial = state(300.0, 0.0, 0.0);
    let duration = 0.1;
    let accepted = advance_supported_chamber(
        &[surface(4, &law, &binding, initial, &substrate, &deposit_h)],
        None,
        &gas,
        &catalog,
        duration,
        controls(),
        || false,
    )
    .unwrap();
    let row = accepted.surfaces()[0];
    let mass = incoming_rate(3000.0, 400.0, AREA) * duration;
    let heat = mass * (LATENT + 3.5 * r() * 400.0 - 4.0 * r() * 300.0);
    let temperature = 300.0 + heat / (SUBSTRATE_MASS * SUBSTRATE_CP + mass * chart.1);
    assert_eq!(row.exchange.outgoing_mass_kg, 0.0);
    fixture::near(row.exchange.incoming_mass_kg, mass, 1e-18);
    fixture::near(row.state.deposit_mass_kg, mass, 1e-18);
    assert_eq!(row.state.substrate_mass_kg, initial.substrate_mass_kg);
    fixture::near(row.exchange.surface_temperature_k[1], temperature, 1e-9);
    assert!(temperature > 300.0);
    assert!(accepted.gas().state().state().temperature < 400.0);
    fixture::near(row.exchange.recession_m, -mass / 917.0 / AREA, 1e-18);
    fixture::near(
        row.substrate_energy_change_j
            + row.deposit_energy_change_j
            + accepted.step().gas_energy_change_j,
        0.0,
        1e-7,
    );
    fixture::near(accepted.step().mass_residual_kg, 0.0, 1e-14);
    fixture::near(accepted.step().energy_residual_j, 0.0, 1e-7);
    assert_eq!(accepted.catalog(), &catalog);
    assert_eq!(accepted.gas().point().axes().get("local-lot"), Some(&9.0));
}

#[test]
fn bare_surfaces_share_species_budget_and_one_order_independent_gas_endpoint() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let mut second = state(350.0, 0.0, 0.0);
    second.area_m2 *= 2.0;
    let rows = [
        surface(
            20,
            &law,
            &binding,
            state(300.0, 0.0, 0.0),
            &substrate,
            &deposit_h,
        ),
        surface(10, &law, &binding, second, &substrate, &deposit_h),
    ];
    let accepted =
        advance_supported_chamber(&rows, None, &gas, &catalog, 0.1, controls(), || false).unwrap();
    let permuted = advance_supported_chamber(
        &[rows[1], rows[0]],
        None,
        &gas,
        &catalog,
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(accepted, permuted);
    assert_eq!(
        accepted.surfaces().iter().map(|s| s.id).collect::<Vec<_>>(),
        [10, 20]
    );
    let mass = 3.0 * incoming_rate(3000.0, 400.0, AREA) * 0.1;
    fixture::near(accepted.step().transferred_mass_kg, -mass, 1e-18);
    let carrier = SpeciesId::new("A").unwrap();
    let carrier_mass = gas.species_mass_kg(&carrier).unwrap();
    fixture::near(
        accepted.gas().species_mass_kg(&carrier).unwrap(),
        carrier_mass,
        1e-14,
    );
    let vapor_mass = gas.species_mass_kg(&fixture::species()).unwrap();
    fixture::near(
        accepted.gas().species_mass_kg(&fixture::species()).unwrap(),
        vapor_mass - mass,
        1e-14,
    );
    let capacity = carrier_mass * 2.5 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897
        + (vapor_mass - mass) * 3.0 * r();
    let expected_t = 400.0 - mass * 0.5 * r() * 400.0 / capacity;
    fixture::near(accepted.gas().state().state().temperature, expected_t, 2e-8);
    fixture::near(accepted.step().energy_residual_j, 0.0, 1e-7);
    fixture::near(accepted.step().species[0].mass_residual_kg, 0.0, 1e-14);
}

#[test]
fn warm_coating_supplies_an_originally_bare_cold_substrate_through_carrier_gas() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let warm = state(600.0, 0.01, 0.0);
    let cold = state(220.0, 0.0, 0.0);
    let rows = [
        surface(1, &law, &binding, warm, &substrate, &deposit_h),
        surface(2, &law, &binding, cold, &substrate, &deposit_h),
    ];
    let first =
        advance_supported_chamber(&rows, None, &gas, &catalog, 0.1, controls(), || false).unwrap();
    assert!(first.surfaces()[0].state.deposit_mass_kg < warm.deposit_mass_kg);
    assert!(first.surfaces()[0].exchange.surface_temperature_k[1] < 600.0);
    assert_eq!(first.surfaces()[1].state, cold);
    assert_eq!(first.surfaces()[1].exchange.incoming_mass_kg, 0.0);
    assert!(first.gas().species_mass_kg(&fixture::species()).unwrap() > 0.0);
    let next = [
        GasSupportedSurface {
            state: first.surfaces()[0].state,
            ..rows[0]
        },
        GasSupportedSurface {
            state: first.surfaces()[1].state,
            ..rows[1]
        },
    ];
    let second =
        advance_supported_chamber(&next, None, first.gas(), &catalog, 0.1, controls(), || {
            false
        })
        .unwrap();
    assert!(second.surfaces()[1].state.deposit_mass_kg > 0.0);
    assert!(second.surfaces()[1].exchange.surface_temperature_k[1] > 220.0);
    for row in second.surfaces() {
        assert_eq!(row.state.substrate_mass_kg, SUBSTRATE_MASS);
    }
    let deposit_mass = second
        .surfaces()
        .iter()
        .map(|row| row.state.deposit_mass_kg)
        .sum::<f64>();
    fixture::near(
        deposit_mass + second.gas().species_mass_kg(&fixture::species()).unwrap(),
        warm.deposit_mass_kg,
        1e-14,
    );
    let final_solid = second
        .surfaces()
        .iter()
        .map(|row| {
            row.state.substrate_mass_kg * row.state.substrate_specific_enthalpy_j_kg
                + row.state.deposit_mass_kg * row.exchange.solid.specific_enthalpy_j_kg
        })
        .sum::<f64>();
    let initial_solid = SUBSTRATE_MASS
        * (warm.substrate_specific_enthalpy_j_kg + cold.substrate_specific_enthalpy_j_kg)
        + warm.deposit_mass_kg * deposit_h(600.0).unwrap();
    fixture::near(
        final_solid - initial_solid
            + first.step().gas_energy_change_j
            + second.step().gas_energy_change_j,
        0.0,
        2e-7,
    );
    assert_ne!(second.gas().state().state().temperature, 400.0);
}

#[test]
fn bare_saturation_is_not_stationary_and_competing_incoming_streams_share_old_inventory() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let rows = [
        surface(
            1,
            &law,
            &binding,
            state(400.0, 0.0, 0.0),
            &substrate,
            &deposit_h,
        ),
        surface(
            2,
            &law,
            &binding,
            state(400.0, 0.0, 0.0),
            &substrate,
            &deposit_h,
        ),
    ];
    let duration = 0.15 * gas.species_mass_kg(&fixture::species()).unwrap()
        / incoming_rate(1000.0, 400.0, AREA);
    advance_supported_chamber(
        &rows[..1],
        None,
        &gas,
        &catalog,
        duration,
        controls(),
        || false,
    )
    .unwrap();
    let error =
        advance_supported_chamber(&rows, None, &gas, &catalog, duration, controls(), || false)
            .unwrap_err();
    let GasSurfaceError::InventoryLimit { maximum_duration_s } = error else {
        panic!("expected aggregate original vapor allowance, got {error:?}");
    };
    assert!(maximum_duration_s > 0.0 && maximum_duration_s < duration);
    advance_supported_chamber(
        &rows,
        None,
        &gas,
        &catalog,
        maximum_duration_s,
        controls(),
        || false,
    )
    .unwrap();
    assert!(matches!(
        advance_supported_chamber(
            &rows[..1],
            None,
            &gas,
            &catalog,
            2.0 * duration,
            controls(),
            || false
        ),
        Err(GasSurfaceError::InventoryLimit { .. })
    ));
}

#[test]
fn carrier_reflection_heats_bare_substrates_in_the_same_conservative_endpoint() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let rows = [
        surface(
            1,
            &law,
            &binding,
            state(300.0, 0.0, 0.0),
            &substrate,
            &deposit_h,
        ),
        surface(
            2,
            &law,
            &binding,
            state(500.0, 0.0, 0.0),
            &substrate,
            &deposit_h,
        ),
    ];
    let reflect = ThermalAccommodation::new(
        &catalog,
        &[
            (SpeciesId::new("A").unwrap(), 0.5),
            (fixture::species(), 0.0),
        ],
        [200.0, 900.0],
    )
    .unwrap();
    let maps = BTreeMap::from([(1, reflect.clone()), (2, reflect)]);
    let accepted = advance_supported_chamber(
        &rows,
        Some(&maps),
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(accepted.gas(), &gas);
    let ra = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let heat = 0.5 * 100000.0 / (2.0 * core::f64::consts::PI * ra * 400.0).sqrt()
        * 3.0
        * ra
        * 100.0
        * AREA
        * 0.001;
    for (row, sign) in accepted.surfaces().iter().zip([-1.0, 1.0]) {
        assert_eq!(row.state.deposit_mass_kg, 0.0);
        assert_eq!(row.exchange.transferred_mass_kg, 0.0);
        fixture::near(row.reflected_heat_into_gas_j, sign * heat, 1e-12);
        assert_eq!(
            row.exchange.energy_into_gas_j,
            row.reflected_heat_into_gas_j
        );
        fixture::near(row.substrate_energy_change_j, -sign * heat, 1e-9);
        fixture::near(
            row.exchange.surface_temperature_k[1],
            row.exchange.surface_temperature_k[0] - sign * heat / (SUBSTRATE_MASS * SUBSTRATE_CP),
            1e-9,
        );
    }
    assert_eq!(accepted.step().mass_residual_kg, 0.0);
    fixture::near(accepted.step().energy_residual_j, 0.0, 1e-7);
}

#[test]
fn independent_substrate_and_large_species_datums_preserve_supported_physics() {
    let base_catalog = fixture::catalog(0.0);
    let chart = fixture::chart(&base_catalog, &fixture::binding(0.0));
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let mut baseline = None;
    for (gas_datum, substrate_datum) in [(0.0, 0.0), (1e14, 1e8), (-1e14, -1e8)] {
        let catalog = fixture::catalog(gas_datum);
        let binding = fixture::binding(gas_datum);
        let gas = fixture::gas(&catalog, 400.0, 3000.0);
        let caloric = |h| Ok(400.0 + (h - substrate_datum) / SUBSTRATE_CP);
        let initial = state(300.0, 0.0, substrate_datum);
        let mut controls = controls();
        controls.chamber.surface.latent_compatibility_tolerance_j_kg = 0.125;
        controls.storage.absolute_energy_tolerance_j = 1e-9;
        let accepted = advance_supported_chamber(
            &[surface(1, &law, &binding, initial, &caloric, &deposit_h)],
            None,
            &gas,
            &catalog,
            0.1,
            controls,
            || false,
        )
        .unwrap();
        let row = accepted.surfaces()[0];
        let alignment = (catalog
            .species_specific_enthalpy_j_kg(&fixture::species(), 300.0)
            .unwrap()
            - (deposit_h(300.0).unwrap() + gas_datum))
            - LATENT;
        let correction = row.exchange.incoming_mass_kg * alignment
            / (SUBSTRATE_MASS * SUBSTRATE_CP + row.state.deposit_mass_kg * chart.1);
        let state = (
            row.state.deposit_mass_kg,
            row.exchange.surface_temperature_k[1] - correction,
            accepted.gas().state().state().temperature,
        );
        if let Some(first) = baseline {
            let (mass, solid_t, gas_t): (f64, f64, f64) = first;
            assert_eq!(state.0, mass);
            fixture::near(state.1, solid_t, 1e-9);
            fixture::near(state.2, gas_t, 2e-8);
        } else {
            baseline = Some(state);
        }
        fixture::near(accepted.step().energy_residual_j, 0.0, 1e-7);
        fixture::near(accepted.step().mass_residual_kg, 0.0, 1e-14);
    }
}

#[test]
fn empty_zero_exposure_and_refusals_preserve_original_supported_and_gas_inputs() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let saved = gas.clone();
    let initial = state(300.0, 0.0, 0.0);
    let rows = [surface(1, &law, &binding, initial, &substrate, &deposit_h)];
    let empty =
        advance_supported_chamber(&[], None, &gas, &catalog, 1.0, controls(), || false).unwrap();
    assert!(empty.surfaces().is_empty());
    assert_eq!(empty.gas(), &gas);
    for (duration, area) in [(0.0, AREA), (100.0, 0.0)] {
        let mut row = rows[0];
        row.state.area_m2 = area;
        let accepted =
            advance_supported_chamber(&[row], None, &gas, &catalog, duration, controls(), || false)
                .unwrap();
        assert_eq!(accepted.surfaces()[0].state, row.state);
        assert_eq!(accepted.gas(), &gas);
    }
    let broken = |h| {
        if h == initial.substrate_specific_enthalpy_j_kg {
            substrate(h)
        } else {
            Err(SublimationError::Invalid(
                "outside original substrate chart",
            ))
        }
    };
    let late_failure = [GasSupportedSurface {
        substrate_temperature: &broken,
        ..rows[0]
    }];
    assert!(
        advance_supported_chamber(&late_failure, None, &gas, &catalog, 0.1, controls(), || {
            false
        })
        .is_err()
    );
    assert!(
        advance_supported_chamber(
            &[rows[0], rows[0]],
            None,
            &gas,
            &catalog,
            0.1,
            controls(),
            || false
        )
        .is_err()
    );
    let mut limited = controls();
    limited.chamber.maximum_surfaces = 0;
    assert!(
        advance_supported_chamber(&rows, None, &gas, &catalog, 0.1, limited, || false).is_err()
    );
    let mut storage = controls();
    storage.storage.max_iterations = 0;
    assert!(advance_supported_chamber(&[], None, &gas, &catalog, 0.1, storage, || false).is_err());
    let wrong = fixture::binding(1000.0);
    assert!(matches!(
        advance_supported_chamber(
            &[GasSupportedSurface {
                binding: &wrong,
                ..rows[0]
            }],
            None,
            &gas,
            &catalog,
            0.1,
            controls(),
            || false
        ),
        Err(GasSurfaceError::CaloricMismatch { .. })
    ));
    let mut polls = 0;
    advance_supported_chamber(&rows, None, &gas, &catalog, 0.1, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            advance_supported_chamber(&rows, None, &gas, &catalog, 0.1, controls(), || {
                count += 1;
                count == stop
            }),
            Err(GasSurfaceError::Surface(SublimationError::Cancelled))
        ));
        assert_eq!(gas, saved);
        assert_eq!(rows[0].state, initial);
    }
}
