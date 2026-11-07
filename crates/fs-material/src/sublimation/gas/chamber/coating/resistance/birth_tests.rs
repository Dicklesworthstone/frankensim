//! G0/G1/G3: finite-G film birth, independent molecular/storage balances.
use super::*;
use crate::sublimation::gas::{chamber::GasChamberSolve, tests as fixture};
use fs_qty::chemistry::SpeciesId;
use fs_thermal::sublimation::SolverConfig;

const SUPPORT_MASS: f64 = 0.02;
const SUPPORT_CP: f64 = 500.0;
const AREA: f64 = 1e-4;
const SUPPORT_T: f64 = 350.0;
const GAS_T: f64 = 400.0;
const VAPOR_P: f64 = 3000.0;

fn controls() -> GasSupportedChamberSolve {
    GasSupportedChamberSolve {
        chamber: GasChamberSolve {
            surface: fixture::controls(),
            maximum_surfaces: 8,
        },
        storage: SolverConfig {
            relative_tolerance: 1e-12,
            absolute_energy_tolerance_j: 1e-10,
            max_iterations: 128,
        },
    }
}

fn state(chart: (f64, f64), support_datum: f64) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: SUPPORT_MASS,
        substrate_specific_enthalpy_j_kg: support_datum + SUPPORT_CP * SUPPORT_T,
        deposit: SolidSurfaceState {
            mass_kg: 0.0,
            specific_enthalpy_j_kg: chart.0 + chart.1 * (SUPPORT_T - 400.0),
            area_m2: AREA,
        },
    }
}

fn balance(result: &GasResistiveChamberCandidate<'_>) {
    fixture::near(result.step().mass_residual_kg, 0.0, 1e-14);
    fixture::near(result.step().energy_residual_j, 0.0, 1e-7);
    fixture::near(
        result.step().solid_energy_change_j + result.step().gas_energy_change_j,
        0.0,
        1e-7,
    );
    for report in result.surfaces() {
        fixture::near(report.substrate_energy_residual_j, 0.0, 1e-7);
        fixture::near(report.deposit_energy_residual_j, 0.0, 1e-7);
        assert_eq!(report.state.substrate_mass_kg, SUPPORT_MASS);
    }
}

fn birth<'catalog>(
    rows: &[GasResistiveSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    dt: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveChamberCandidate<'catalog>, GasSurfaceError> {
    advance_resistive_chamber_with_births(rows, None, None, gas, catalog, dt, controls(), cancelled)
}

#[test]
fn incoming_molecules_create_actual_film_storage_without_seed_or_datum_dependence() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, GAS_T, VAPOR_P);
    let law = fixture::phase();
    let dt = 0.1;
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let incoming = dt * AREA * 0.1 * VAPOR_P / (2.0 * core::f64::consts::PI * r * GAS_T).sqrt();
    let incoming_energy = 3.5 * r * GAS_T;
    for conductance in [0.0, 2.0, 2000.0] {
        let mut reference = None;
        for (support_datum, deposit_offset) in [(0.0, 0.0), (1e6, 5e5)] {
            let binding = fixture::binding(deposit_offset);
            let chart = fixture::chart(&catalog, &binding);
            let support = |h| Ok((h - support_datum) / SUPPORT_CP);
            let deposit = |h| fixture::temperature(h, chart);
            let initial = state(chart, support_datum);
            let rows = [GasResistiveSurface {
                id: 1,
                law: &law,
                binding: &binding,
                state: initial,
                conductance_w_k: conductance,
                substrate_temperature: &support,
                deposit_temperature: &deposit,
            }];
            let born = birth(&rows, &gas, &catalog, dt, || false).unwrap();
            let report = &born.surfaces()[0];
            // Independent backward-Euler two-body balance using the incoming
            // mass as the ENTIRE final deposit heat capacity.
            let heat = incoming
                * (incoming_energy - initial.deposit.specific_enthalpy_j_kg - deposit_offset);
            let cs = SUPPORT_MASS * SUPPORT_CP;
            let cd = incoming * chart.1;
            let q =
                -dt * conductance * heat / cd / (1.0 + dt * conductance * (1.0 / cs + 1.0 / cd));
            fixture::near(report.state.deposit.mass_kg, incoming, 1e-16);
            assert_eq!(report.exchange.outgoing_mass_kg, 0.0);
            fixture::near(report.exchange.incoming_mass_kg, incoming, 1e-16);
            fixture::near(report.exchange.transferred_mass_kg, -incoming, 1e-16);
            fixture::near(report.interface_heat_j, q, 1e-8);
            fixture::near(report.substrate_temperature_k[1], SUPPORT_T - q / cs, 1e-8);
            fixture::near(
                report.exchange.surface_temperature_k[1],
                SUPPORT_T + (heat + q) / cd,
                1e-7,
            );
            fixture::near(
                report.exchange.energy_into_gas_j,
                -incoming * incoming_energy,
                1e-10,
            );
            let carrier = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap();
            let vapor = gas.species_mass_kg(&fixture::species()).unwrap();
            let carrier_cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
            let capacity = carrier * carrier_cv + vapor * 3.0 * r;
            let gas_t =
                (capacity * GAS_T - incoming * incoming_energy) / (capacity - incoming * 3.0 * r);
            fixture::near(born.gas().state().state().temperature, gas_t, 1e-8);
            fixture::near(
                born.gas().species_mass_kg(&fixture::species()).unwrap(),
                vapor - incoming,
                1e-16,
            );
            balance(&born);
            let temperatures = (
                report.substrate_temperature_k[1],
                report.exchange.surface_temperature_k[1],
            );
            if let Some((ts, td)) = reference {
                fixture::near(temperatures.0, ts, 1e-8);
                fixture::near(temperatures.1, td, 1e-7);
            }
            reference = Some(temperatures);
            // Re-entry uses the ordinary positive-film owner and the SAME gas
            // endpoint. No special born-film state or mass floor persists.
            let next = [GasResistiveSurface {
                state: report.state,
                ..rows[0]
            }];
            let continued = advance_resistive_chamber(
                &next,
                None,
                born.gas(),
                &catalog,
                1e-5,
                controls(),
                || false,
            )
            .unwrap();
            assert!(continued.surfaces()[0].state.deposit.mass_kg > 0.0);
            balance(&continued);
            assert_eq!(rows[0].state.deposit.mass_kg, 0.0);
        }
    }
}

#[test]
fn film_resistance_and_reflected_carrier_heat_share_one_newborn_energy_account() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, GAS_T, VAPOR_P);
    let binding = fixture::binding(0.0);
    let law = fixture::phase();
    let chart = fixture::chart(&catalog, &binding);
    let support = |h| Ok(h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let rows = [GasResistiveSurface {
        id: 1,
        law: &law,
        binding: &binding,
        state: state(chart, 0.0),
        conductance_w_k: 20.0,
        substrate_temperature: &support,
        deposit_temperature: &deposit,
    }];
    let reflection = ThermalAccommodation::new(
        &catalog,
        &[
            (SpeciesId::new("A").unwrap(), 0.5),
            (fixture::species(), 0.0),
        ],
        [200.0, 900.0],
    )
    .unwrap();
    let dt = 0.001;
    let film = 0.2;
    let result = advance_resistive_chamber_with_births(
        &rows,
        Some(&BTreeMap::from([(1, reflection)])),
        Some(&BTreeMap::from([(1, film)])),
        &gas,
        &catalog,
        dt,
        controls(),
        || false,
    )
    .unwrap();
    let report = &result.surfaces()[0];
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let ra = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let kinetic = 0.1 / (2.0 * core::f64::consts::PI * r * GAS_T).sqrt();
    let gas_film = film / (r * GAS_T);
    let incoming = AREA * dt * kinetic * VAPOR_P * gas_film / (gas_film + kinetic);
    let reflected = AREA * dt * 0.5 * (100000.0 - VAPOR_P)
        / (2.0 * core::f64::consts::PI * ra * GAS_T).sqrt()
        * 3.0
        * ra
        * (SUPPORT_T - GAS_T);
    fixture::near(report.state.deposit.mass_kg, incoming, 1e-16);
    fixture::near(report.reflected_heat_into_gas_j, reflected, 1e-12);
    fixture::near(
        report.exchange.energy_into_gas_j,
        -incoming * 3.5 * r * GAS_T + reflected,
        1e-10,
    );
    assert_eq!(report.exchange.outgoing_mass_kg, 0.0);
    assert!(report.substrate_temperature_k[1] > SUPPORT_T);
    assert!(report.exchange.surface_temperature_k[1] > report.substrate_temperature_k[1]);
    balance(&result);
}

#[test]
fn birth_admission_refuses_unmatched_calorics_missing_incoming_inventory_and_cancellation() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, GAS_T, VAPOR_P);
    let saved = gas.clone();
    let binding = fixture::binding(0.0);
    let law = fixture::phase();
    let chart = fixture::chart(&catalog, &binding);
    let support = |h| Ok(h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let rows = [GasResistiveSurface {
        id: 1,
        law: &law,
        binding: &binding,
        state: state(chart, 0.0),
        conductance_w_k: 2.0,
        substrate_temperature: &support,
        deposit_temperature: &deposit,
    }];
    assert!(
        advance_resistive_chamber(&rows, None, &gas, &catalog, 0.1, controls(), || false).is_err()
    );
    let mut mismatch = rows;
    mismatch[0].state.deposit.specific_enthalpy_j_kg += chart.1;
    assert!(birth(&mismatch, &gas, &catalog, 0.1, || false).is_err());
    assert!(birth(&rows, &gas, &catalog, 0.0, || false).is_err());
    let dry = fixture::gas(&catalog, GAS_T, 0.0);
    assert!(birth(&rows, &dry, &catalog, 0.1, || false).is_err());
    assert!(matches!(
        birth(&rows, &gas, &catalog, 100.0, || false),
        Err(GasSurfaceError::InventoryLimit { .. })
    ));
    let mut calls = 0;
    birth(&rows, &gas, &catalog, 0.1, || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(
            birth(&rows, &gas, &catalog, 0.1, || {
                count += 1;
                count == stop
            })
            .is_err()
        );
        assert_eq!(gas, saved);
        assert_eq!(rows[0].state.deposit.mass_kg, 0.0);
    }
}
