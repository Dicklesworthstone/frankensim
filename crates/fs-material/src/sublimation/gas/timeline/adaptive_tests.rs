//! Coupled temporal accuracy on the real heat, phase and gas owners.
use super::*;
use crate::sublimation::gas::{GasSurfaceBinding, tests as fixture};
use crate::sublimation::gas::timeline::{advance_deposition_with_accommodation, advance_deposition_with_gas};
use fs_thermal::enthalpy::{Boundary, CellGeometry, DepositedCell, DepositionBudget,
    MaterialId, Reservoir, ThermalLaw, ThermalLink, TimedEvent};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

#[derive(Clone, Copy)]
struct Caloric(f64, f64);
impl Caloric {
    fn h(self, t: f64) -> f64 { self.0 + self.1 * (t - 400.0) }
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 400.0 + (h - self.0) / self.1;
        if (200.0..=900.0).contains(&t) { Ok(t) } else { Err(DepositionError::MaterialDomain) }
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 / self.1 }
}
fn chart(catalog: &ThermochemicalGasModel, binding: &GasSurfaceBinding) -> Caloric {
    let (h, cp) = fixture::chart(catalog, binding); Caloric(h, cp)
}
fn phase_controls() -> GasChamberSolve {
    GasChamberSolve { surface: fixture::controls(), maximum_surfaces: 2 }
}
fn heat(dt: f64) -> StepConfig { StepConfig { max_dt_s: dt, energy_tolerance_j: 1e-7 } }
fn controls(tolerance: f64, cp: f64) -> GasTimelineAdaptiveSolve {
    GasTimelineAdaptiveSolve {
        absolute_mass_tolerance_kg: 1e-11, absolute_temperature_tolerance_k: tolerance,
        absolute_specific_enthalpy_tolerance_j_kg: cp * tolerance,
        relative_tolerance: 0.0, minimum_dt_s: 1e-9, maximum_trials: 10000,
    }
}
fn boundary(power: f64, cover: f64) -> Boundary {
    Boundary { ambient_k: 350.0, htc_w_m2_k: 0.0, source_w: power,
        reservoir: Some(Reservoir { area_m2: cover, conductance_w_k: 0.0, temperature_k: 350.0 }) }
}
fn cell(id: u64, t: f64, chart: Caloric, power: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id), segment_id: id, body_id: id, material_id: MaterialId(1),
        mass_kg: 0.01, specific_h_j_kg: chart.h(t),
        geometry: CellGeometry { start_m: [id as f64, 0.0, 0.0], end_m: [id as f64, 0.01, 0.0],
            volume_m3: 0.01 / 917.0, surface_area_m2: 2e-4 },
        boundary: boundary(power, 1e-4),
    }
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent { TimedEvent { time_s, event: Event::Birth(cell) } }
fn schedule(chart: &Caloric, events: Vec<TimedEvent>, steps: usize) -> DepositionSimulation<'_> {
    DepositionSimulation::new(events, &[(MaterialId(1), chart)], DepositionBudget {
        max_cells: 2, max_links: 1, max_events: 16, max_substeps: steps,
    }).unwrap()
}
fn selection<'a>(phase: &'a dyn fs_thermal::sublimation::SublimationLaw,
    binding: &'a GasSurfaceBinding, chart: Caloric) -> [GasTimelineSurface<'a>; 2] {
    std::array::from_fn(|i| GasTimelineSurface { cell_id: CellId(i as u64 + 1), law: phase, binding,
        pure_solid_enthalpy_range_j_kg: [chart.h(200.0), chart.h(900.0)] })
}
fn laws(catalog: &ThermochemicalGasModel, alpha: f64) -> BTreeMap<CellId, ThermalAccommodation> {
    let law = ThermalAccommodation::new(catalog, &[
        (SpeciesId::new("A").unwrap(), alpha), (fixture::species(), alpha),
    ], [200.0, 900.0]).unwrap();
    BTreeMap::from([(CellId(1), law.clone()), (CellId(2), law)])
}
fn temperatures(candidate: &GasTimelineCandidate<'_, '_>) -> [f64; 3] {
    let t = |id| {
        let cell = &candidate.thermal().cells()[&CellId(id)];
        candidate.thermal().temperature_at_specific_enthalpy_k(CellId(id), cell.specific_h_j_kg).unwrap()
    };
    [t(1), t(2), candidate.gas().state().state().temperature]
}
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

// Independent ODE reference: two finite solid stores, original specified
// conductance, authored heater, and ideal-diatomic reflected heat into one gas.
// This integrates simultaneous rates with RK4, NOT the production split kernel.
fn reference(cp: f64, gas_mass: f64, duration: f64) -> [f64; 3] {
    let r = R / 0.02897;
    let cs = 0.01 * cp;
    let cg = gas_mass * 2.5 * r;
    let factor = 0.5 * 9e-5 * gas_mass / 0.01 * (r / (2.0 * core::f64::consts::PI)).sqrt() * 3.0 * r;
    let rates = |y: [f64; 3]| {
        let conductance = 10.0 * (y[0] - y[1]);
        let q1 = factor * y[2].sqrt() * (y[0] - y[2]);
        let q2 = factor * y[2].sqrt() * (y[1] - y[2]);
        [(200.0 - conductance - q1) / cs, (conductance - q2) / cs, (q1 + q2) / cg]
    };
    let mut y = [400.0, 300.0, 350.0];
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let k1 = rates(y);
        let k2 = rates(std::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]));
        let k3 = rates(std::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]));
        let k4 = rates(std::array::from_fn(|i| y[i] + dt * k3[i]));
        y = std::array::from_fn(|i| y[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}

#[test]
fn complete_heat_contact_reflection_split_refines_against_independent_simultaneous_rates() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0); binding.accommodation_coefficient = 0.0;
    let chart = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, chart);
    let thermal = schedule(&chart, vec![birth(0.0, cell(1, 400.0, chart, 200.0)),
        birth(0.0, cell(2, 300.0, chart, 0.0)),
        TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: 9, a: CellId(1), b: CellId(2), area_m2: 1e-5, conductance_w_k: 10.0,
        }) }], 100000);
    let gas = fixture::gas(&catalog, 350.0, 0.0);
    let coefficients = laws(&catalog, 0.5);
    let expected = reference(chart.1, gas.mass_kg(), 0.25);
    let run = |tolerance| advance_deposition_adaptive_with_accommodation(&thermal, &selected,
        &coefficients, &gas, &catalog, 0.25, heat(0.125), phase_controls(),
        controls(tolerance, chart.1), || Ok(())).unwrap();
    let loose = run(0.02);
    let tight = run(1e-5);
    let errors = [distance(temperatures(loose.endpoint()), expected), distance(temperatures(tight.endpoint()), expected)];
    assert!(errors[1] < 0.1 * errors[0] && errors[1] < 0.005, "{errors:?}");
    assert!(tight.audit().accepted_intervals > loose.audit().accepted_intervals);
    assert!(tight.audit().rejected_trials > 0);
    assert!(tight.audit().maximum_accepted_error_ratio <= 1.0);
    let report = tight.endpoint().audit();
    fixture::near(report.thermal.external_heat_j, 50.0, 1e-10);
    fixture::near(report.thermal.incoming_mass_kg, 0.02, 1e-14);
    assert_eq!(report.transferred_mass_kg, 0.0);
    assert_eq!(tight.endpoint().gas().mass_kg(), gas.mass_kg());
    fixture::near(report.energy_into_gas_j, report.reflected_heat_into_gas_j, 1e-9);
    fixture::near(report.energy_residual_j, 0.0, 1e-7);
    fixture::near(tight.audit().gas_stream_energy_residual_j, 0.0, 1e-7);
    fixture::near(tight.audit().solid_stream_energy_residual_j, 0.0, 1e-7);
    assert!(tight.audit().thermal_substeps_evaluated > report.thermal.substeps);
    assert!(thermal.cells().is_empty());
}

#[test]
fn phase_inventory_retries_and_trace_species_share_the_full_coupled_controller() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, chart);
    let thermal = schedule(&chart, vec![birth(0.0, cell(1, 400.0, chart, 1.0)),
        birth(0.0, cell(2, 300.0, chart, 1.0))], 100000);
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let coefficients = laws(&catalog, 0.3);
    let accepted = advance_deposition_adaptive_with_accommodation(&thermal, &selected,
        &coefficients, &gas, &catalog, 1.0, heat(1.0), phase_controls(), controls(0.001, chart.1), || Ok(())).unwrap();
    assert!(accepted.audit().rejected_trials > 0);
    assert!(accepted.audit().maximum_accepted_error_ratio <= 1.0);
    assert!(accepted.endpoint().audit().cells[&CellId(1)].transferred_mass_kg > 0.0);
    assert!(accepted.endpoint().audit().cells[&CellId(2)].transferred_mass_kg < 0.0);
    fixture::near(accepted.audit().maximum_species_mass_residual_kg, 0.0, 1e-14);
    fixture::near(accepted.endpoint().audit().mass_residual_kg, 0.0, 1e-14);
    fixture::near(accepted.endpoint().audit().thermal.external_heat_j, 2.0, 1e-10);
    fixture::near(accepted.endpoint().gas().species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(),
        gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(), 1e-14);
    let reference = advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 1.0, heat(1.0 / 2048.0), phase_controls(), || Ok(())).unwrap();
    assert!(distance(temperatures(accepted.endpoint()), temperatures(&reference)) < 0.03);
}

#[test]
fn births_and_endpoint_coverage_are_not_replayed_into_the_accepted_heat_history() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0); binding.accommodation_coefficient = 0.0;
    let chart = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, chart);
    let thermal = schedule(&chart, vec![birth(0.03125, cell(1, 400.0, chart, 0.0)),
        birth(0.03125, cell(2, 340.0, chart, 0.0)),
        TimedEvent { time_s: 0.0625, event: Event::Boundary(CellId(1), boundary(0.0, 2e-4)) },
        TimedEvent { time_s: 0.0625, event: Event::Boundary(CellId(2), boundary(0.0, 2e-4)) }], 100);
    let gas = fixture::gas(&catalog, 300.0, 0.0);
    let coefficients = laws(&catalog, 0.5);
    let accepted = advance_deposition_adaptive_with_accommodation(&thermal, &selected,
        &coefficients, &gas, &catalog, 0.0625, heat(1.0), phase_controls(), controls(100.0, chart.1), || Ok(())).unwrap();
    let direct = advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 0.0625, heat(0.015625), phase_controls(), || Ok(())).unwrap();
    assert_eq!(accepted.endpoint().gas(), direct.gas());
    assert_eq!(accepted.endpoint().thermal().cells(), direct.thermal().cells());
    assert_eq!(accepted.endpoint().audit().thermal.events, 4);
    assert_eq!(accepted.endpoint().audit().phase_steps, 2);
    assert_eq!(accepted.audit().dormant_advances, 1);
    assert_eq!(accepted.audit().accepted_intervals, 1);
    assert_eq!(accepted.audit().thermal_substeps_evaluated, 3);
    for id in [CellId(1), CellId(2)] {
        assert_eq!(accepted.endpoint().thermal().exposed_area_m2(id), Ok(0.0));
        fixture::near(accepted.endpoint().audit().cells[&id].reflected_heat_into_gas_j,
            direct.audit().cells[&id].reflected_heat_into_gas_j, 1e-10);
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert!(thermal.cells().is_empty());
}

#[test]
fn specular_defaults_and_surface_permutations_share_identical_adaptive_decisions() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, chart);
    let thermal = schedule(&chart, vec![birth(0.0, cell(1, 400.0, chart, 1.0)),
        birth(0.0, cell(2, 300.0, chart, 1.0))], 100000);
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let original = advance_deposition_adaptive(&thermal, &selected, &gas, &catalog,
        0.125, heat(0.0625), phase_controls(), controls(0.001, chart.1), || Ok(())).unwrap();
    let zeros = advance_deposition_adaptive_with_accommodation(&thermal, &selected, &laws(&catalog, 0.0),
        &gas, &catalog, 0.125, heat(0.0625), phase_controls(), controls(0.001, chart.1), || Ok(())).unwrap();
    let reversed = advance_deposition_adaptive(&thermal, &[selected[1], selected[0]], &gas, &catalog,
        0.125, heat(0.0625), phase_controls(), controls(0.001, chart.1), || Ok(())).unwrap();
    for other in [&zeros, &reversed] {
        assert_eq!(original.endpoint().gas(), other.endpoint().gas());
        assert_eq!(original.endpoint().thermal().cells(), other.endpoint().thermal().cells());
        assert_eq!(original.endpoint().audit(), other.endpoint().audit());
        assert_eq!(original.audit(), other.audit());
    }
}

#[test]
fn all_trial_heat_work_spends_one_budget_and_late_cancellation_publishes_nothing() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, chart);
    let events = vec![birth(0.0, cell(1, 400.0, chart, 1.0)), birth(0.0, cell(2, 300.0, chart, 1.0))];
    let thermal = schedule(&chart, events.clone(), 2);
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let saved = gas.clone();
    assert!(matches!(advance_deposition_adaptive(&thermal, &selected, &gas, &catalog,
        0.03125, heat(0.03125), phase_controls(), controls(100.0, chart.1), || Ok(())),
        Err(GasTimelineAdaptiveError::Timeline(GasTimelineError::Thermal(DepositionError::Budget(_))))));
    let thermal = schedule(&chart, events, 10000);
    let mut calls = 0;
    advance_deposition_adaptive(&thermal, &selected, &gas, &catalog,
        0.03125, heat(0.03125), phase_controls(), controls(0.001, chart.1), || { calls += 1; Ok(()) }).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(advance_deposition_adaptive(&thermal, &selected, &gas, &catalog,
            0.03125, heat(0.03125), phase_controls(), controls(0.001, chart.1), || {
                count += 1;
                if count == stop { Err(DepositionError::Cancelled) } else { Ok(()) }
            }).is_err());
        assert!(thermal.cells().is_empty());
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(gas, saved);
    }
}

#[test]
fn unselected_isothermal_plateau_inventory_is_part_of_the_error_norm() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if (0.0..=1000.0).contains(&h) { Ok(300.0) } else { Err(DepositionError::MaterialDomain) }
        }
        fn maximum_temperature_slope(&self) -> f64 { 0.0 }
    }
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, chart);
    let gas = fixture::gas(&catalog, 350.0, 0.0);
    let plateau = Plateau;
    let make = |h| {
        let mut second = cell(2, 300.0, chart, 0.0);
        second.material_id = MaterialId(2); second.specific_h_j_kg = h;
        DepositionSimulation::new(vec![birth(0.0, cell(1, 400.0, chart, 0.0)), birth(0.0, second)],
            &[(MaterialId(1), &chart as &dyn ThermalLaw), (MaterialId(2), &plateau)],
            DepositionBudget { max_cells: 2, max_links: 1, max_events: 2, max_substeps: 100 }).unwrap()
    };
    let base = make(100.0);
    let changed = make(102.0);
    let a = advance_deposition_with_gas(&base, &selected[..1], &gas, &catalog, 0.0,
        heat(1.0), phase_controls(), || Ok(())).unwrap();
    let b = advance_deposition_with_gas(&changed, &selected[..1], &gas, &catalog, 0.0,
        heat(1.0), phase_controls(), || Ok(())).unwrap();
    let mut scales = controls(100.0, chart.1);
    scales.absolute_specific_enthalpy_tolerance_j_kg = 1.0;
    assert_eq!(compare(a.thermal(), &gas, &a, &b, scales, &mut || Ok(())).unwrap(), 2.0);
}

#[test]
fn dormant_time_zero_duration_admission_trial_and_minimum_limits_are_bounded() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, chart);
    let thermal = schedule(&chart, vec![birth(2e9, cell(1, 400.0, chart, 0.0)),
        birth(2e9, cell(2, 300.0, chart, 0.0))], 1);
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let mut bounded = controls(0.001, chart.1); bounded.maximum_trials = 1;
    let dormant = advance_deposition_adaptive(&thermal, &selected, &gas, &catalog,
        1e9, heat(1e-6), phase_controls(), bounded, || Ok(())).unwrap();
    assert_eq!(dormant.endpoint().thermal().time_s(), 1e9);
    assert_eq!(dormant.endpoint().gas(), &gas);
    assert_eq!(dormant.audit().trials, 0);
    assert_eq!(dormant.audit().interval_evaluations, 2);
    assert_eq!(dormant.audit().thermal_substeps_evaluated, 0);
    assert!(advance_deposition_adaptive_with_accommodation(&thermal, &selected, &BTreeMap::new(),
        &gas, &catalog, 0.0, heat(1.0), phase_controls(), bounded, || Ok(())).is_err());
    let active = schedule(&chart, vec![birth(0.0, cell(1, 400.0, chart, 500.0)),
        birth(0.0, cell(2, 300.0, chart, 0.0))], 10000);
    assert!(matches!(advance_deposition_adaptive(&active, &selected, &gas, &catalog,
        1.0, heat(1.0), phase_controls(), bounded, || Ok(())),
        Err(GasTimelineAdaptiveError::TrialBudget { .. })));
    bounded.minimum_dt_s = 1.0;
    assert!(matches!(advance_deposition_adaptive(&active, &selected, &gas, &catalog,
        1.0, heat(1.0), phase_controls(), bounded, || Ok(())),
        Err(GasTimelineAdaptiveError::MinimumStep { .. })));
}
