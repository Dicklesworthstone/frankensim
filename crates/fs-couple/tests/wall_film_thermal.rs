//! Dry onset and exact film exhaustion on the original finite thermal wall.
use std::cell::Cell;

use fs_couple::finite_gas_thermal::{FiniteGasThermalEnvironment, wall_film::WallFilmThermalError};
use fs_material::evaporation::gas::wall_film::GasWallFilmInventory;
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, EvaporationProperties};

#[path = "wall_film_thermal/fixture.rs"]
mod fixture;
use fixture::*;

#[test]
fn dry_condensation_warms_finite_wall_and_retains_actual_film_capacity() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let input = input(&curve, 300.0, 0.1);
    let film = GasWallFilmInventory::Dry;
    let next = scope(|cx, _| {
        environment(&gas, &source, film, &binding, &body)
            .advance(cx, input)
            .unwrap()
    });
    let report = &next.report;
    assert!(report.phase.film_formed);
    assert!(report.phase.transferred_mass_kg < 0.0);
    assert!(report.next_film.mass_kg() > 0.0);
    assert!(next.state.temperature_k() > input.initial.temperature_k());
    assert!(report.wall_heat_into_phase_system_j < 0.0);
    assert!(report.phase.vapor_partial_pressure_pa[1] < report.phase.vapor_partial_pressure_pa[0]);
    assert_eq!(report.solid_to_gas_heat_j, 0.0);
    assert_eq!(report.phase.dry_out_time_s, None);

    // Independent finite-capacity balance. New liquid has its own actual
    // heat capacity; no film seed or capacity is added to the wall material.
    let condensed = report.next_film.mass_kg();
    let capacity = input.mass_kg * CP_WALL + condensed * CP_LIQUID;
    let source_gap = LATENT + 3.5 * R * (320.0 - 300.0);
    near(
        next.state.temperature_k(),
        300.0 + condensed * source_gap / capacity,
        1e-7,
    );
    near(
        liquid_t(report.next_film.specific_enthalpy_j_kg().unwrap()).unwrap(),
        next.state.temperature_k(),
        controls().temperature_tolerance_k,
    );
    complete_inventory(input, &gas, film, &binding, 0.0, &next);
    assert_eq!(gas, original_gas);
}

#[test]
fn complete_evaporation_keeps_outgoing_vapor_capacity_and_full_sensible_interval() {
    let source = source();
    let gas = gas(&source, 290.0, 0.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(20.0);
    let input = input(&curve, 300.0, 0.1);
    let evaporated = 1e-7;
    let film = GasWallFilmInventory::Wet {
        mass_kg: evaporated,
        specific_enthalpy_j_kg: liquid_h(300.0).unwrap(),
    };
    let next = scope(|cx, _| {
        environment(&gas, &source, film, &binding, &body)
            .advance(cx, input)
            .unwrap()
    });
    let report = &next.report;
    assert_eq!(report.next_film, GasWallFilmInventory::Dry);
    assert_eq!(report.next_film.specific_enthalpy_j_kg(), None);
    assert_eq!(report.phase.transferred_mass_kg, evaporated);
    assert_eq!(report.phase.duration_s, input.duration_s);
    let dry_out_time = report.phase.dry_out_time_s.unwrap();
    assert!(dry_out_time > 0.0 && dry_out_time < input.duration_s);
    assert_eq!(report.phase.phase_active_duration_s, dry_out_time);
    assert!(report.solid_to_gas_heat_j > 0.0);
    assert!(next.state.temperature_k() < input.initial.temperature_k());
    assert!(report.next_gas.state().state().temperature > 290.0);

    // Independent constant-capacity 2x2 endpoint system. Evaporated vapor
    // leaves at the accepted wall temperature: its Cp remains in the wall
    // equation even though final liquid mass and liquid capacity are zero.
    let cp_vapor = 3.5 * R;
    let cv_vapor = 2.5 * R;
    let wall_capacity = input.mass_kg * CP_WALL;
    let gas_capacity = 2.5 * 100000.0 * GAS_VOLUME / 290.0;
    let conductance_time = body.convection_w_per_m2_k * AREA * input.duration_s;
    let a = wall_capacity + evaporated * cp_vapor + conductance_time;
    let b = -conductance_time;
    let c = -(conductance_time + evaporated * cp_vapor);
    let d = gas_capacity + evaporated * cv_vapor + conductance_time;
    let wall_rhs = (wall_capacity + evaporated * cp_vapor) * 300.0 - evaporated * LATENT;
    let gas_rhs = gas_capacity * 290.0;
    let determinant = a * d - b * c;
    let expected_wall = (wall_rhs * d - b * gas_rhs) / determinant;
    let expected_gas = (a * gas_rhs - c * wall_rhs) / determinant;
    near(next.state.temperature_k(), expected_wall, 1e-7);
    near(
        report.next_gas.state().state().temperature,
        expected_gas,
        1e-7,
    );
    near(
        report.solid_to_gas_heat_j,
        conductance_time * (expected_wall - expected_gas),
        2e-9,
    );
    complete_inventory(input, &gas, film, &binding, 0.0, &next);
    assert_eq!(gas, original_gas);
}

#[test]
fn dry_undersaturated_wall_reduces_to_original_reciprocal_sensible_owner() {
    let source = source();
    let gas = gas(&source, 290.0, 0.0);
    let binding = binding(&source);
    let curve = curve();
    let body = body(20.0);
    let input = input(&curve, 300.0, 0.1);
    let film = GasWallFilmInventory::Dry;
    let original_isochoric = gas.clone().into_isochoric();
    let (next, baseline) = scope(|cx, _| {
        let next = environment(&gas, &source, film, &binding, &body)
            .advance(cx, input)
            .unwrap();
        let baseline = FiniteGasThermalEnvironment {
            gas: &original_isochoric,
            body: &body,
            gas_heat_j: 0.0,
            gas_solve: controls().phase.gas.gas,
            exchange: controls().sensible,
            maximum_gas_temperature_increment_k: 10.0,
        }
        .advance(cx, input)
        .unwrap();
        (next, baseline)
    });
    assert_eq!(next.report.next_film, GasWallFilmInventory::Dry);
    assert!(!next.report.phase.film_formed);
    assert_eq!(next.report.phase.transferred_mass_kg, 0.0);
    assert_eq!(next.report.phase.film_enthalpy_change_j, 0.0);
    assert_eq!(next.report.wall_heat_into_phase_system_j, 0.0);
    near(
        next.state.temperature_k(),
        baseline.state.temperature_k(),
        1e-7,
    );
    near(
        next.report.next_gas.state().state().temperature,
        baseline.report.next_gas.state().state().temperature,
        1e-7,
    );
    near(
        next.report.solid_to_gas_heat_j,
        baseline.report.solid_to_gas_heat_j,
        2e-9,
    );
    complete_inventory(input, &gas, film, &binding, 0.0, &next);
}

struct CountedPhase(Cell<usize>);
impl EvaporationLaw for CountedPhase {
    fn temperature_range_k(&self) -> [f64; 2] {
        PHASE.temperature_range_k()
    }
    fn properties(&self, temperature: f64) -> Result<EvaporationProperties, EvaporationError> {
        self.0.set(self.0.get() + 1);
        PHASE.properties(temperature)
    }
}

#[test]
fn late_original_wall_refusal_and_cancellation_do_not_publish_phase_candidates() {
    let source = source();
    let gas = gas(&source, 320.0, 6000.0);
    let original_gas = gas.clone();
    let binding = binding(&source);
    let curve = curve();
    let body = body(0.0);
    let film = GasWallFilmInventory::Dry;
    let valid_input = input(&curve, 300.0, 0.1);
    let original_wall = valid_input.initial;
    let law = CountedPhase(Cell::new(0));
    let mut environment = environment(&gas, &source, film, &binding, &body);
    environment.boundary.law = &law;
    let mut refused_input = valid_input;
    // Volume is validated by the original wall thermal owner after the
    // otherwise valid gas/film trial has been assembled.
    refused_input.volume_m3 = 0.0;
    scope(|cx, gate| {
        assert!(matches!(
            environment.advance(cx, refused_input),
            Err(WallFilmThermalError::Thermal(_))
        ));
        assert!(law.0.get() > 0);
        let calls_after_refusal = law.0.get();
        gate.request();
        assert!(matches!(
            environment.advance(cx, valid_input),
            Err(WallFilmThermalError::Cancelled)
        ));
        assert_eq!(law.0.get(), calls_after_refusal);
    });
    assert_eq!(gas, original_gas);
    assert_eq!(environment.film, GasWallFilmInventory::Dry);
    assert_eq!(valid_input.initial, original_wall);
}
