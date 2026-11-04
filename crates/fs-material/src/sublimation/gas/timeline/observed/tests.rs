//! Actual molecular, gas, thermal and conduction owners; manufactured calorics.
use super::*;
use crate::sublimation::gas::{GasSurfaceBinding, tests as fixture};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, Reservoir,
    ThermalLaw, TimedEvent,
    conduction::{ConductionObserver, ConductiveReservoir, ConductivityLaw, PlanarPhaseLayer},
};

#[derive(Clone, Copy)]
struct Caloric(f64, f64);
impl Caloric {
    fn h(self, t: f64) -> f64 {
        self.0 + self.1 * (t - 400.0)
    }
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 400.0 + (h - self.0) / self.1;
        if (200.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.1
    }
}
struct Conductivity;
impl ConductivityLaw for Conductivity {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(1.0)
    }
}
const AREA: f64 = 1e-4;
const MASS: f64 = 1e-4;
const DENSITY: f64 = 917.0;
fn thermal(caloric: &Caloric) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(1),
                body_id: 1,
                segment_id: 1,
                material_id: MaterialId(1),
                mass_kg: MASS,
                specific_h_j_kg: caloric.h(300.0),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.01, 0.0, 0.0],
                    volume_m3: MASS / DENSITY,
                    surface_area_m2: 2.0 * AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: Some(Reservoir {
                        area_m2: AREA,
                        temperature_k: 280.0,
                        conductance_w_k: 0.0,
                    }),
                },
            }),
        }],
        &[(MaterialId(1), caloric)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, heat(0.001), || Ok(()))
    .unwrap()
    .0
}
fn heat(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-7,
    }
}
fn controls() -> GasChamberSolve {
    GasChamberSolve {
        surface: fixture::controls(),
        maximum_surfaces: 1,
    }
}
fn selected<'a>(
    law: &'a dyn fs_thermal::sublimation::SublimationLaw,
    binding: &'a GasSurfaceBinding,
    caloric: Caloric,
) -> [GasTimelineSurface<'a>; 1] {
    [GasTimelineSurface {
        cell_id: CellId(1),
        law,
        binding,
        pure_solid_enthalpy_range_j_kg: [caloric.h(200.0), caloric.h(900.0)],
    }]
}

#[test]
fn empty_observer_preserves_original_trajectory_and_audits_bitwise() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (h, cp) = fixture::chart(&catalog, &binding);
    let caloric = Caloric(h, cp);
    let thermal = thermal(&caloric);
    let gas = fixture::gas(&catalog, 350.0, 3000.0);
    let phase = fixture::phase();
    let surfaces = selected(&phase, &binding, caloric);
    let plain = super::super::advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas,
        &catalog,
        0.01,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let (observed, ()) = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        None,
        &(),
        &gas,
        &catalog,
        0.01,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(observed.gas(), plain.gas());
    assert_eq!(observed.thermal().cells(), plain.thermal().cells());
    assert_eq!(observed.audit(), plain.audit());
}

#[test]
fn deposited_mass_changes_actual_layer_resistance_on_the_next_heat_interval() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (h, cp) = fixture::chart(&catalog, &binding);
    let caloric = Caloric(h, cp);
    let thermal = thermal(&caloric);
    let original = thermal.cells().clone();
    let gas = fixture::gas(&catalog, 350.0, 3000.0);
    let phase = fixture::phase();
    let surfaces = selected(&phase, &binding, caloric);
    let initial_distance = MASS / (2.0 * DENSITY * AREA);
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(1), &Conductivity)],
        &[],
        &[ConductiveReservoir {
            cell_id: CellId(1),
            half_distance_m: initial_distance,
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_planar_phase_layers(
        &thermal,
        &[PlanarPhaseLayer {
            cell_id: CellId(1),
            footprint_area_m2: AREA,
            density_kg_m3: DENSITY,
        }],
        || Ok(()),
    )
    .unwrap();
    let (first, observed) = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        None,
        &conduction,
        &gas,
        &catalog,
        0.001,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let mass = first.thermal().cells()[&CellId(1)].mass_kg;
    assert!(mass > MASS);
    fixture::near(
        observed.last_planar_half_distance_m(CellId(1)).unwrap(),
        initial_distance,
        1e-16,
    );
    let (second, final_observer) = advance_deposition_with_gas(
        first.thermal(),
        &surfaces,
        None,
        &observed,
        first.gas(),
        &catalog,
        0.002,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    fixture::near(
        final_observer
            .last_planar_half_distance_m(CellId(1))
            .unwrap(),
        mass / (2.0 * DENSITY * AREA),
        1e-16,
    );
    assert!(
        final_observer
            .last_planar_half_distance_m(CellId(1))
            .unwrap()
            > initial_distance
    );
    let before_t = first.thermal().temperature_k(CellId(1)).unwrap();
    let expected_heat = AREA / (mass / (2.0 * DENSITY * AREA)) * (280.0 - before_t) * 0.001;
    fixture::near(second.audit().thermal.external_heat_j, expected_heat, 1e-13);
    for candidate in [&first, &second] {
        let report = candidate.audit();
        fixture::near(report.mass_residual_kg, 0.0, 1e-14);
        fixture::near(
            report.solid_energy_change_j + report.gas_energy_change_j
                - report.incoming_energy_j
                - report.thermal.external_heat_j,
            0.0,
            1e-7,
        );
    }
    assert_eq!(thermal.cells(), &original);
    fixture::near(
        conduction.last_planar_half_distance_m(CellId(1)).unwrap(),
        initial_distance,
        1e-16,
    );
}

#[derive(Clone, Debug, PartialEq)]
enum SourceError {
    Thermal(DepositionError),
    Rejected,
}
impl From<DepositionError> for SourceError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
#[derive(Clone, Debug, Default, PartialEq)]
struct Pulse {
    begin_calls: usize,
    accepted_duration_s: f64,
    accepted_steps: usize,
    fail_after: Option<usize>,
    interval_heat_j: f64,
}
impl DepositionObserver for Pulse {
    type Error = SourceError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.begin_calls += 1;
        self.accepted_duration_s = 0.0;
        self.accepted_steps = 0;
        self.interval_heat_j = 0.0;
        Ok(())
    }
    fn next_transition_time_s(
        &mut self,
        state: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        Ok((state.time_s() < 0.002).then_some(0.002))
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        _: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let watts = if start < 0.002 { 10.0 } else { -5.0 };
        self.interval_heat_j += watts * dt;
        Ok(BTreeMap::from([(CellId(1), watts)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.accepted_duration_s += dt;
        self.accepted_steps += 1;
        if self
            .fail_after
            .is_some_and(|limit| self.accepted_steps >= limit)
        {
            return Err(SourceError::Rejected);
        }
        Ok(())
    }
}

#[test]
fn interval_heat_and_transition_hooks_share_the_physical_clock_and_are_counted_once() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (h, cp) = fixture::chart(&catalog, &binding);
    let caloric = Caloric(h, cp);
    let thermal = thermal(&caloric);
    let gas = fixture::gas(&catalog, 350.0, 3000.0);
    let phase = fixture::phase();
    let surfaces = selected(&phase, &binding, caloric);
    let pulse = Pulse::default();
    let (candidate, accepted) = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        None,
        &pulse,
        &gas,
        &catalog,
        0.004,
        heat(0.004),
        controls(),
        || Ok(()),
    )
    .unwrap();
    fixture::near(accepted.accepted_duration_s, 0.004, 1e-15);
    fixture::near(accepted.interval_heat_j, 0.01, 1e-14);
    fixture::near(candidate.audit().thermal.external_heat_j, 0.01, 1e-14);
    assert_eq!(candidate.audit().phase_steps, 2);
    assert_eq!(accepted.accepted_steps, 2);
    assert_eq!(accepted.begin_calls, 1);
    fixture::near(candidate.audit().energy_residual_j, 0.0, 1e-7);
    assert_eq!(pulse, Pulse::default());
}

#[test]
fn a_late_source_refusal_preserves_heat_gas_and_observer_state() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let (h, cp) = fixture::chart(&catalog, &binding);
    let caloric = Caloric(h, cp);
    let thermal = thermal(&caloric);
    let original_cells = thermal.cells().clone();
    let gas = fixture::gas(&catalog, 350.0, 3000.0);
    let original_gas = gas.clone();
    let phase = fixture::phase();
    let surfaces = selected(&phase, &binding, caloric);
    let pulse = Pulse {
        fail_after: Some(3),
        ..Pulse::default()
    };
    let saved_pulse = pulse.clone();
    let result = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        None,
        &pulse,
        &gas,
        &catalog,
        0.004,
        heat(0.001),
        controls(),
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(GasTimelineObservedError::Observer(SourceError::Rejected))
    ));
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(gas, original_gas);
    assert_eq!(pulse, saved_pulse);
}

#[path = "adaptive_tests.rs"]
mod adaptive;
