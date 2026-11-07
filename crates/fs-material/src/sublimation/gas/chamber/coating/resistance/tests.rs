//! Independent two-body and molecular/ideal-gas balances on the original owners.
use super::*;
use crate::sublimation::gas::{chamber::GasChamberSolve, tests as fixture};
use fs_qty::chemistry::SpeciesId;
use fs_thermal::sublimation::SolverConfig;

const SUPPORT_CP: f64 = 500.0;
const SUPPORT_MASS: f64 = 0.02;
const FILM_MASS: f64 = 0.002;
const AREA: f64 = 1e-4;

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
fn state(support_t: f64, deposit_t: f64, datum: f64, chart: (f64, f64)) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: SUPPORT_MASS,
        substrate_specific_enthalpy_j_kg: datum + SUPPORT_CP * (support_t - 400.0),
        deposit: SolidSurfaceState {
            mass_kg: FILM_MASS,
            area_m2: AREA,
            specific_enthalpy_j_kg: chart.0 + chart.1 * (deposit_t - 400.0),
        },
    }
}
fn support(h: f64) -> Result<f64, SublimationError> {
    Ok(400.0 + h / SUPPORT_CP)
}
fn row<'a>(
    id: u64,
    state: ResistiveSurfaceState,
    conductance: f64,
    law: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    substrate_temperature: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    deposit_temperature: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasResistiveSurface<'a> {
    GasResistiveSurface {
        id,
        state,
        conductance_w_k: conductance,
        law,
        binding,
        substrate_temperature,
        deposit_temperature,
    }
}
fn balanced(candidate: &GasResistiveChamberCandidate<'_>) {
    fixture::near(candidate.step().mass_residual_kg, 0.0, 1e-14);
    fixture::near(candidate.step().energy_residual_j, 0.0, 1e-7);
    fixture::near(
        candidate.step().solid_energy_change_j + candidate.step().gas_energy_change_j,
        0.0,
        1e-7,
    );
    for surface in candidate.surfaces() {
        assert_eq!(surface.state.substrate_mass_kg, SUPPORT_MASS);
        fixture::near(surface.substrate_energy_residual_j, 0.0, 1e-7);
        fixture::near(surface.deposit_energy_residual_j, 0.0, 1e-7);
    }
}

#[test]
fn sealed_mass_film_retains_two_temperatures_and_original_caloric_datums() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let law = fixture::phase();
    let mut endpoints = Vec::new();
    for (datum, offset) in [(0.0, 0.0), (1e6, 5e5)] {
        let binding = fixture::binding(offset);
        let chart = fixture::chart(&catalog, &binding);
        let substrate = |h| support(h - datum);
        let deposit = |h| fixture::temperature(h, chart);
        let initial = state(300.0, 350.0, datum, chart);
        let accepted = advance_resistive_chamber_with_films(
            &[row(8, initial, 2.0, &law, &binding, &substrate, &deposit)],
            None,
            &BTreeMap::from([(8, 0.0)]),
            &gas,
            &catalog,
            0.2,
            controls(),
            || false,
        )
        .unwrap();
        let report = &accepted.surfaces()[0];
        // Backward Euler of two independent constant-capacity bodies. No
        // production heat function is used in this reference calculation.
        let cs = SUPPORT_MASS * SUPPORT_CP;
        let cd = FILM_MASS * chart.1;
        let q = 0.2 * 2.0 * (300.0 - 350.0) / (1.0 + 0.2 * 2.0 * (1.0 / cs + 1.0 / cd));
        fixture::near(report.interface_heat_j, q, 1e-8);
        fixture::near(report.substrate_temperature_k[1], 300.0 - q / cs, 1e-8);
        fixture::near(
            report.exchange.surface_temperature_k[1],
            350.0 + q / cd,
            1e-8,
        );
        assert!(report.substrate_temperature_k[1] < report.exchange.surface_temperature_k[1]);
        assert_eq!(report.state.deposit.mass_kg, FILM_MASS);
        assert_eq!(report.exchange.transferred_mass_kg, 0.0);
        assert_eq!(accepted.gas(), &gas);
        balanced(&accepted);
        endpoints.push((
            report.substrate_temperature_k[1],
            report.exchange.surface_temperature_k[1],
        ));
    }
    fixture::near(endpoints[0].0, endpoints[1].0, 1e-8);
    fixture::near(endpoints[0].1, endpoints[1].1, 1e-8);
}

#[test]
fn finite_gas_molecular_energy_heats_the_deposit_before_resistive_sharing() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let initial = state(300.0, 350.0, 0.0, chart);
    let dt = 0.1;
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let saturation = 1000.0 * (200000.0 / r * (1.0 / 400.0 - 1.0 / 350.0)).exp();
    let outgoing = dt * AREA * 0.1 * saturation / (2.0 * core::f64::consts::PI * r * 350.0).sqrt();
    let incoming = dt * AREA * 0.1 * 3000.0 / (2.0 * core::f64::consts::PI * r * 400.0).sqrt();
    let dm = outgoing - incoming;
    let stream = outgoing * 3.5 * r * 350.0 - incoming * 3.5 * r * 400.0;
    let next_mass = FILM_MASS - dm;
    let sensible_heat = dm * initial.deposit.specific_enthalpy_j_kg - stream;
    let cs = SUPPORT_MASS * SUPPORT_CP;
    let cd = next_mass * chart.1;
    let carrier = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap();
    let vapor = gas.species_mass_kg(&fixture::species()).unwrap();
    let carrier_cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let expected_gas_t = ((carrier * carrier_cv + vapor * 3.0 * r) * 400.0 + stream)
        / (carrier * carrier_cv + (vapor + dm) * 3.0 * r);
    let mut difference = Vec::new();
    for conductance in [0.0, 2.0, 2e3] {
        let accepted = advance_resistive_chamber(
            &[row(
                3,
                initial,
                conductance,
                &law,
                &binding,
                &support,
                &deposit,
            )],
            None,
            &gas,
            &catalog,
            dt,
            controls(),
            || false,
        )
        .unwrap();
        let report = &accepted.surfaces()[0];
        let q = dt * conductance * (300.0 - 350.0 - sensible_heat / cd)
            / (1.0 + dt * conductance * (1.0 / cs + 1.0 / cd));
        fixture::near(report.state.deposit.mass_kg, next_mass, 1e-16);
        fixture::near(report.exchange.outgoing_mass_kg, outgoing, 1e-16);
        fixture::near(report.exchange.incoming_mass_kg, incoming, 1e-16);
        fixture::near(report.interface_heat_j, q, 1e-8);
        fixture::near(report.substrate_temperature_k[1], 300.0 - q / cs, 1e-8);
        fixture::near(
            report.exchange.surface_temperature_k[1],
            350.0 + (sensible_heat + q) / cd,
            1e-8,
        );
        fixture::near(
            accepted.gas().state().state().temperature,
            expected_gas_t,
            1e-8,
        );
        fixture::near(
            accepted.gas().species_mass_kg(&fixture::species()).unwrap(),
            vapor + dm,
            1e-16,
        );
        assert!(report.state.deposit.mass_kg > initial.deposit.mass_kg);
        assert_eq!(accepted.gas().point().axes().get("local-lot"), Some(&9.0));
        balanced(&accepted);
        difference.push(
            (report.exchange.surface_temperature_k[1] - report.substrate_temperature_k[1]).abs(),
        );
    }
    assert!(difference[2] < 0.03 * difference[1] && difference[1] < difference[0]);
}

#[test]
fn one_gas_is_shared_canonically_and_reflection_reaches_the_actual_deposit() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let rows = [
        row(
            42,
            state(300.0, 350.0, 0.0, chart),
            0.5,
            &law,
            &binding,
            &support,
            &deposit,
        ),
        row(
            7,
            state(330.0, 375.0, 0.0, chart),
            3.0,
            &law,
            &binding,
            &support,
            &deposit,
        ),
    ];
    let accommodation = ThermalAccommodation::new(
        &catalog,
        &[
            (SpeciesId::new("A").unwrap(), 0.25),
            (fixture::species(), 0.5),
        ],
        [200.0, 900.0],
    )
    .unwrap();
    let reflections = BTreeMap::from([(7, accommodation.clone()), (42, accommodation.clone())]);
    let accepted = advance_resistive_chamber(
        &rows,
        Some(&reflections),
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    let replay = advance_resistive_chamber(
        &[rows[1], rows[0]],
        Some(&reflections),
        &gas,
        &catalog,
        0.001,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(accepted, replay);
    assert_eq!(
        accepted
            .surfaces()
            .iter()
            .map(|row| row.id)
            .collect::<Vec<_>>(),
        vec![7, 42]
    );
    for (initial, report) in [rows[1], rows[0]].iter().zip(accepted.surfaces()) {
        let reference = accommodation
            .evaluate(
                &binding,
                &gas,
                deposit(initial.state.deposit.specific_enthalpy_j_kg).unwrap(),
                AREA,
                0.001,
                || false,
            )
            .unwrap();
        assert_eq!(report.reflected_heat_into_gas_j, reference.heat_into_gas_j);
        assert!(report.reflected_heat_into_gas_j < 0.0);
    }
    balanced(&accepted);
}

#[test]
fn internal_heat_disables_stationary_turnover_exception_and_refusals_publish_nothing() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let original = state(300.0, 400.0, 0.0, chart);
    let insulated = row(8, original, 0.0, &law, &binding, &support, &deposit);
    let retained =
        advance_resistive_chamber(&[insulated], None, &gas, &catalog, 1e6, controls(), || {
            false
        })
        .unwrap();
    assert_eq!(retained.surfaces()[0].state, original);
    assert_eq!(retained.gas(), &gas);
    let conducting = GasResistiveSurface {
        conductance_w_k: 1.0,
        ..insulated
    };
    assert!(
        matches!(advance_resistive_chamber(&[conducting], None, &gas, &catalog,
        1e6, controls(), || false), Err(GasSurfaceError::InventoryLimit { maximum_duration_s })
        if maximum_duration_s > 0.0 && maximum_duration_s < 1e6)
    );
    let bare = GasResistiveSurface {
        state: ResistiveSurfaceState {
            deposit: SolidSurfaceState {
                mass_kg: 0.0,
                ..original.deposit
            },
            ..original
        },
        ..conducting
    };
    assert!(
        advance_resistive_chamber(&[bare], None, &gas, &catalog, 0.001, controls(), || false)
            .is_err()
    );
    assert!(
        advance_resistive_chamber_with_films(
            &[conducting],
            None,
            &BTreeMap::from([(99, 0.1)]),
            &gas,
            &catalog,
            0.001,
            controls(),
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    advance_resistive_chamber(
        &[conducting],
        None,
        &gas,
        &catalog,
        0.001,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, polls] {
        let mut seen = 0;
        assert!(matches!(
            advance_resistive_chamber(
                &[conducting],
                None,
                &gas,
                &catalog,
                0.001,
                controls(),
                || {
                    seen += 1;
                    seen == stop
                }
            ),
            Err(GasSurfaceError::Surface(SublimationError::Cancelled))
        ));
    }
    assert_eq!(conducting.state, original);
    assert_eq!(gas, fixture::gas(&catalog, 400.0, 1000.0));
}
