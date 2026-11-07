//! Independent driven heat limits and accepted finite-gas source histories.
use super::super::{
    GasSupportedChamberSolve, advance_resistive_chamber, advance_resistive_chamber_driven,
};
use super::*;
use crate::sublimation::gas::{GasSurfaceBinding, chamber::GasChamberSolve, tests as fixture};
use fs_thermal::sublimation::{
    SolidSurfaceState, SolverConfig, SublimationError, SublimationLaw,
    coating::resistance::ResistiveSurfaceState,
};

const SUPPORT_MASS: f64 = 0.02;
const SUPPORT_CP: f64 = 500.0;
const FILM_MASS: f64 = 0.002;

fn controls() -> GasResistiveAdvanceSolve {
    GasResistiveAdvanceSolve {
        step: GasSupportedChamberSolve {
            chamber: GasChamberSolve {
                surface: fixture::controls(),
                maximum_surfaces: 8,
            },
            storage: SolverConfig {
                relative_tolerance: 1e-12,
                absolute_energy_tolerance_j: 1e-10,
                max_iterations: 128,
            },
        },
        absolute_mass_tolerance_kg: 1e-9,
        absolute_temperature_tolerance_k: 0.01,
        absolute_specific_enthalpy_tolerance_j_kg: 2.0,
        absolute_heat_tolerance_j: 0.001,
        relative_tolerance: 1e-6,
        minimum_dt_s: 1e-10,
        maximum_dt_s: 2.0,
        maximum_trials: 20000,
    }
}
fn state(ts: f64, td: f64, chart: (f64, f64), datum: f64) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: SUPPORT_MASS,
        substrate_specific_enthalpy_j_kg: datum + SUPPORT_CP * (ts - 400.0),
        deposit: SolidSurfaceState {
            mass_kg: FILM_MASS,
            area_m2: 1e-4,
            specific_enthalpy_j_kg: chart.0 + chart.1 * (td - 400.0),
        },
    }
}
fn row<'a>(
    state: ResistiveSurfaceState,
    conductance: f64,
    law: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    support: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    deposit: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasResistiveSurface<'a> {
    GasResistiveSurface {
        id: 7,
        state,
        conductance_w_k: conductance,
        law,
        binding,
        substrate_temperature: support,
        deposit_temperature: deposit,
    }
}

#[test]
fn independent_heating_cooling_and_datums_match_the_two_body_linear_balance() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let law = fixture::phase();
    let films = BTreeMap::from([(7, 0.0)]);
    let qs = -12.0;
    let qd = 5.0;
    let g = 2.0;
    for (support_datum, deposit_offset) in [(0.0, 0.0), (1e7, -3e6)] {
        let binding = fixture::binding(deposit_offset);
        let chart = fixture::chart(&catalog, &binding);
        let support = |h| Ok(400.0 + (h - support_datum) / SUPPORT_CP);
        let deposit = |h| fixture::temperature(h, chart);
        let initial = state(350.0, 400.0, chart, support_datum);
        let rows = [row(initial, g, &law, &binding, &support, &deposit)];
        let drive = BTreeMap::from([(
            7,
            GasResistiveDrive {
                substrate_heat_j: qs,
                deposit_heat_j: qd,
                interface: None,
            },
        )]);
        // A zero-duration prescribed-energy call is a real independent impulse.
        for dt in [0.0, 0.2] {
            let result = advance_resistive_chamber_driven(
                &rows,
                None,
                Some(&films),
                &drive,
                &gas,
                &catalog,
                dt,
                controls().step,
                || false,
            )
            .unwrap();
            let accepted = result.surfaces()[0];
            let cs = SUPPORT_MASS * SUPPORT_CP;
            let cd = FILM_MASS * chart.1;
            let q = dt * g * (350.0 + qs / cs - 400.0 - qd / cd)
                / (1.0 + dt * g * (1.0 / cs + 1.0 / cd));
            fixture::near(accepted.interface_heat_j, q, 1e-8);
            fixture::near(
                accepted.substrate_temperature_k[1],
                350.0 + (qs - q) / cs,
                2e-9,
            );
            fixture::near(
                accepted.exchange.surface_temperature_k[1],
                400.0 + (qd + q) / cd,
                2e-9,
            );
            assert_eq!(accepted.supplied_substrate_heat_j, qs);
            assert_eq!(accepted.supplied_deposit_heat_j, qd);
            assert_eq!(accepted.exchange.supplied_heat_j, qs + qd);
            assert_eq!(result.step().supplied_heat_j, qs + qd);
            assert_eq!(result.gas(), &gas);
            assert_eq!(accepted.state.deposit.mass_kg, initial.deposit.mass_kg);
            fixture::near(result.step().solid_energy_change_j - (qs + qd), 0.0, 1e-7);
            fixture::near(accepted.substrate_energy_residual_j, 0.0, 1e-7);
            fixture::near(accepted.deposit_energy_residual_j, 0.0, 1e-7);
            fixture::near(result.step().energy_residual_j, 0.0, 1e-7);
        }
    }
}

#[test]
fn endpoint_interface_replaces_scalar_contact_inside_the_same_driven_transaction() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let law = fixture::phase();
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = state(300.0, 350.0, chart, 0.0);
    let films = BTreeMap::from([(7, 0.0)]);
    // This passive nonlinear law has a quadratic primitive, independent of the
    // authored scalar conductance. Check its constitutive endpoint, not a
    // frozen-conductance approximation.
    let interface = |ts: f64, td: f64| Ok(0.002 * (ts - td) * (ts + td));
    let drives = BTreeMap::from([(
        7,
        GasResistiveDrive {
            substrate_heat_j: -2.0,
            deposit_heat_j: 0.5,
            interface: Some(&interface),
        },
    )]);
    let dt = 0.2;
    let accepted = advance_resistive_chamber_driven(
        &[row(initial, f64::NAN, &law, &binding, &support, &deposit)],
        None,
        Some(&films),
        &drives,
        &gas,
        &catalog,
        dt,
        controls().step,
        || false,
    )
    .unwrap();
    let endpoint = accepted.surfaces()[0];
    let ts = endpoint.substrate_temperature_k[1];
    let td = endpoint.exchange.surface_temperature_k[1];
    fixture::near(
        endpoint.interface_heat_j,
        dt * 0.002 * (ts - td) * (ts + td),
        1e-8,
    );
    fixture::near(
        endpoint.substrate_energy_change_j + endpoint.interface_heat_j,
        -2.0,
        1e-8,
    );
    fixture::near(
        endpoint.deposit_energy_change_j - endpoint.interface_heat_j,
        0.5,
        1e-8,
    );
    assert_eq!(accepted.gas(), &gas);
    assert_eq!(accepted.step().supplied_heat_j, -1.5);
}

#[test]
fn rejected_trials_cannot_supply_heat_to_a_finite_gas_deposition_history() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let original_gas = gas.clone();
    let law = fixture::phase();
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = state(300.0, 350.0, chart, 0.0);
    let rows = [row(initial, 2.0, &law, &binding, &support, &deposit)];
    let powers = BTreeMap::from([(
        7,
        GasResistivePower {
            substrate_power_w: -20.0,
            deposit_power_w: 3.0,
            interface: None,
        },
    )]);
    let duration = 2.0;
    let result = advance_resistive_chamber_adaptive_driven(
        &rows,
        None,
        None,
        &powers,
        &gas,
        &catalog,
        duration,
        controls(),
        || false,
    )
    .unwrap();
    let audit = result.audit();
    let thermal = audit.thermal[0];
    assert!(audit.chamber.rejected_trials > 0);
    assert!(audit.chamber.accepted_intervals > 1);
    assert_eq!(audit.chamber.elapsed_s, duration);
    fixture::near(thermal.supplied_substrate_heat_j, -20.0 * duration, 1e-10);
    fixture::near(thermal.supplied_deposit_heat_j, 3.0 * duration, 1e-10);
    fixture::near(audit.chamber.supplied_heat_j, -17.0 * duration, 1e-10);
    fixture::near(
        audit.chamber.surfaces[0].supplied_heat_j,
        -17.0 * duration,
        1e-10,
    );
    assert!(result.endpoint().surfaces()[0].state.deposit.mass_kg > FILM_MASS);
    assert!(
        result
            .endpoint()
            .gas()
            .species_mass_kg(&fixture::species())
            .unwrap()
            < gas.species_mass_kg(&fixture::species()).unwrap()
    );
    fixture::near(
        audit.chamber.solid_energy_change_j + audit.chamber.gas_energy_change_j
            - audit.chamber.supplied_heat_j,
        0.0,
        1e-7,
    );
    fixture::near(audit.chamber.energy_residual_j, 0.0, 1e-7);
    fixture::near(thermal.substrate_energy_residual_j, 0.0, 1e-7);
    fixture::near(thermal.deposit_energy_residual_j, 0.0, 1e-7);
    fixture::near(audit.chamber.mass_residual_kg, 0.0, 1e-14);
    assert_eq!(gas, original_gas);
    assert_eq!(rows[0].state, initial);
}

#[test]
fn opposite_nonzero_sources_disable_the_stationary_inventory_shortcut() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let law = fixture::phase();
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let rows = [row(
        state(400.0, 400.0, chart, 0.0),
        2.0,
        &law,
        &binding,
        &support,
        &deposit,
    )];
    let duration = 1e6;
    let closed = advance_resistive_chamber(
        &rows,
        None,
        &gas,
        &catalog,
        duration,
        controls().step,
        || false,
    )
    .unwrap();
    assert_eq!(closed.gas(), &gas);
    assert_eq!(closed.surfaces()[0].state, rows[0].state);
    assert_eq!(closed.step().supplied_heat_j, 0.0);
    let drives = BTreeMap::from([(
        7,
        GasResistiveDrive {
            substrate_heat_j: -1.0,
            deposit_heat_j: 1.0,
            interface: None,
        },
    )]);
    assert!(matches!(
        advance_resistive_chamber_driven(
            &rows,
            None,
            None,
            &drives,
            &gas,
            &catalog,
            duration,
            controls().step,
            || false,
        ),
        Err(GasSurfaceError::InventoryLimit { .. })
    ));
}

#[test]
fn invalid_sources_and_late_cancellation_leave_original_states_untouched() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let original_gas = gas.clone();
    let law = fixture::phase();
    let support = |h| Ok(400.0 + h / SUPPORT_CP);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = state(300.0, 350.0, chart, 0.0);
    let rows = [row(initial, 2.0, &law, &binding, &support, &deposit)];
    let valid = GasResistivePower {
        substrate_power_w: -20.0,
        deposit_power_w: 3.0,
        interface: None,
    };
    for (id, power) in [
        (8, valid),
        (
            7,
            GasResistivePower {
                deposit_power_w: f64::NAN,
                ..valid
            },
        ),
    ] {
        assert!(
            advance_resistive_chamber_adaptive_driven(
                &rows,
                None,
                None,
                &BTreeMap::from([(id, power)]),
                &gas,
                &catalog,
                0.1,
                controls(),
                || false,
            )
            .is_err()
        );
    }
    assert!(
        GasResistivePower {
            substrate_power_w: f64::MIN_POSITIVE,
            deposit_power_w: 0.0,
            interface: None
        }
        .for_interval(f64::MIN_POSITIVE)
        .is_err()
    );
    assert_eq!(valid.for_interval(0.0).unwrap().substrate_heat_j, 0.0);
    let powers = BTreeMap::from([(7, valid)]);
    let mut polls = 0;
    advance_resistive_chamber_adaptive_driven(
        &rows,
        None,
        None,
        &powers,
        &gas,
        &catalog,
        0.1,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let cancel_at = polls - 1;
    let mut replay_polls = 0;
    let result = advance_resistive_chamber_adaptive_driven(
        &rows,
        None,
        None,
        &powers,
        &gas,
        &catalog,
        0.1,
        controls(),
        || {
            replay_polls += 1;
            replay_polls == cancel_at
        },
    );
    assert!(matches!(
        result,
        Err(GasResistiveAdvanceError::Surface(GasSurfaceError::Surface(
            SublimationError::Cancelled
        )))
    ));
    assert_eq!(gas, original_gas);
    assert_eq!(rows[0].state, initial);
}
