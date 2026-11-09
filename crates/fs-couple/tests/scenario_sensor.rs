//! Real finite-cell heat drives the same scenario instruments at accepted steps.

use fs_couple::scenario_sensor::{
    TemperatureSensorObserver, ThermalSensorError, ThermalSensorSupport,
};
use fs_scenario::{
    CompiledSensorOperator, EntityDeclaration, EntityId, EntityRef, KindExpectation,
    ObservationSupport, ObservationTerm, PlacementUncertainty, ScenarioSensor, SensorCalibration,
    SensorDynamics, SensorKind, SensorLocation, SensorMount,
};
use fs_thermal::enthalpy::*;
use fs_thermal::observers::CoupledObservers;

struct Law;
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > 0.0 {
            Ok(h / 10.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.1
    }
}
static LAW: Law = Law;

fn entity() -> EntityId {
    let assembly = EntityDeclaration::assembly("instrumented-heating").identity();
    let part = EntityDeclaration::part(assembly, "finite-body").identity();
    EntityDeclaration::region(part, "thermal-cells").identity()
}

fn sensor(
    kind: SensorKind,
    support: ObservationSupport,
    tau: Option<f64>,
    gain: f64,
    offset: f64,
) -> CompiledSensorOperator {
    ScenarioSensor::new(
        "mounted-temperature",
        kind,
        SensorLocation::new(
            EntityRef::new(entity(), KindExpectation::Domain),
            [0.0; 3],
            PlacementUncertainty::declared_exact("explicit fixture mapping").unwrap(),
        )
        .unwrap(),
        support,
        SensorMount::affine(gain, offset, "fixture mount").unwrap(),
        match tau {
            Some(tau) => SensorDynamics::first_order(tau, "fixture time constant").unwrap(),
            None => SensorDynamics::declared_instantaneous("ideal reference").unwrap(),
        },
        SensorCalibration::physical("fixture-calibration", "2026-07-01", "synthetic", 0.04)
            .unwrap(),
        false,
    )
    .unwrap()
    .compile()
    .unwrap()
}

fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        source_w,
        reservoir: None,
    }
}
fn cell(id: u64, temperature: f64, mass: f64, watts: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: mass,
        specific_h_j_kg: temperature * 10.0,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 1.0, 0.0],
            volume_m3: 0.001,
            surface_area_m2: 1.0,
        },
        boundary: boundary(watts),
    }
}
fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.125,
        energy_tolerance_j: 1e-9,
    }
}
fn initial(change_source: bool) -> DepositionSimulation<'static> {
    let mut schedule = vec![
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(10, 300.0, 2.0, 40.0)),
        },
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(20, 350.0, 1.0, 30.0)),
        },
    ];
    if change_source {
        schedule.push(TimedEvent {
            time_s: 0.3,
            event: Event::Boundary(CellId(10), boundary(-20.0)),
        });
    }
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), &LAW)],
        DepositionBudget {
            max_cells: 3,
            max_links: 0,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-11, "{a:.15e} != {b:.15e}");
}
fn response(old: f64, target: f64, slope: f64, dt: f64, tau: f64) -> f64 {
    let weight = 1.0 - (-dt / tau).exp();
    old + (target - old) * weight + slope * (dt - tau * weight)
}
fn support(cells: &[CellId]) -> ThermalSensorSupport<'_> {
    ThermalSensorSupport {
        entity: entity(),
        cells,
    }
}

#[test]
fn native_ramp_drives_mount_lag_and_patch_probe_without_altering_heat() {
    let heat = initial(false);
    // Slot ordering deliberately differs from the BTreeMap ordering of cells.
    let cells = [CellId(20), CellId(10)];
    let instrument = sensor(
        SensorKind::Thermocouple,
        ObservationSupport::point(2, 1).unwrap(),
        Some(0.25),
        2.0,
        5.0,
    );
    let patch = sensor(
        SensorKind::Rtd,
        ObservationSupport::patch_average(
            2,
            vec![ObservationTerm::new(0, 0.25), ObservationTerm::new(1, 0.75)],
        )
        .unwrap(),
        None,
        1.0,
        0.0,
    );
    let lagged =
        TemperatureSensorObserver::new(&instrument, support(&cells), &heat, 605.0, || Ok(()))
            .unwrap();
    let average =
        TemperatureSensorObserver::new(&patch, support(&cells), &heat, 312.5, || Ok(())).unwrap();
    let observers = CoupledObservers::new(lagged, average);
    let (observed, next, audit) = heat
        .advance_observed(&observers, 1.0, config(), || Ok(()))
        .unwrap();
    let (plain, plain_audit) = heat.advance_to(1.0, config(), || Ok(())).unwrap();
    assert_eq!(observed.cells(), plain.cells());
    assert_eq!(audit, plain_audit);
    near(
        next.first.state().reading(),
        response(605.0, 605.0, 4.0, 1.0, 0.25),
    );
    near(next.first.state().equilibrium_reading(), 609.0);
    near(next.second.state().reading(), 0.25 * 353.0 + 0.75 * 302.0);
    assert_eq!(
        next.first.state().operator().sensor_identity(),
        instrument.sensor_identity()
    );
    assert_eq!(next.first.cells(), cells);
    assert_eq!(observers.first.state().time_s(), 0.0);
    assert_eq!(heat.temperature_k(CellId(10)).unwrap(), 300.0);

    // An implicit/multiphysics owner may already have solved the endpoint.
    // For this affine ramp, reading that candidate agrees with substep observation.
    let endpoint_sensor = observers
        .first
        .observe_endpoint_candidate(&heat, &plain, 1.0, || Ok(()))
        .unwrap();
    near(
        endpoint_sensor.state().reading(),
        next.first.state().reading(),
    );
    assert_eq!(
        endpoint_sensor.state().equilibrium_reading(),
        next.first.state().equilibrium_reading()
    );
    assert_eq!(plain.cells(), observed.cells());

    let (middle, middle_sensors, _) = heat
        .advance_observed(&observers, 0.5, config(), || Ok(()))
        .unwrap();
    let (restarted, final_sensors, _) = middle
        .advance_observed(&middle_sensors, 1.0, config(), || Ok(()))
        .unwrap();
    assert_eq!(restarted.cells(), observed.cells());
    assert_eq!(final_sensors, next);
}

#[test]
fn native_boundary_event_preserves_sensor_memory_and_resolves_both_ramps() {
    let heat = initial(true);
    let instrument = sensor(
        SensorKind::Thermocouple,
        ObservationSupport::point(1, 0).unwrap(),
        Some(0.2),
        1.0,
        0.0,
    );
    let observer = TemperatureSensorObserver::new(
        &instrument,
        support(&[CellId(10)]),
        &heat,
        300.0,
        || Ok(()),
    )
    .unwrap();
    // Decimal dt exercises native clocks for which (start+dt)-start != dt.
    let config = StepConfig {
        max_dt_s: 0.1,
        ..config()
    };
    let (next, sensor, audit) = observer
        .advance_candidate(&heat, 0.8, config, || Ok(()))
        .unwrap();
    let at_event = response(300.0, 300.0, 2.0, 0.3, 0.2);
    near(
        sensor.state().reading(),
        response(at_event, 300.6, -1.0, 0.5, 0.2),
    );
    near(next.temperature_k(CellId(10)).unwrap(), 300.1);
    assert_eq!(sensor.state().time_s(), next.time_s());
    assert_eq!(audit.events, 1);
    let (_, plain_audit) = heat.advance_to(0.8, config, || Ok(())).unwrap();
    assert_eq!(audit, plain_audit);
}

#[test]
fn cell_mapping_quantity_and_stale_history_refuse_before_publication() {
    let heat = initial(false);
    let instrument = sensor(
        SensorKind::Thermocouple,
        ObservationSupport::point(2, 0).unwrap(),
        Some(0.2),
        1.0,
        0.0,
    );
    for ids in [
        vec![CellId(10)],
        vec![CellId(10), CellId(10)],
        vec![CellId(10), CellId(99)],
    ] {
        assert!(
            TemperatureSensorObserver::new(&instrument, support(&ids), &heat, 300.0, || Ok(()))
                .is_err()
        );
    }
    let other = EntityDeclaration::assembly("unrelated-entity").identity();
    assert!(
        TemperatureSensorObserver::new(
            &instrument,
            ThermalSensorSupport {
                entity: other,
                cells: &[CellId(10), CellId(20)]
            },
            &heat,
            300.0,
            || Ok(())
        )
        .is_err()
    );
    let pressure = sensor(
        SensorKind::PressureTap,
        ObservationSupport::point(1, 0).unwrap(),
        Some(0.2),
        1.0,
        0.0,
    );
    assert!(
        TemperatureSensorObserver::new(&pressure, support(&[CellId(10)]), &heat, 300.0, || Ok(()))
            .is_err()
    );

    let observer = TemperatureSensorObserver::new(
        &instrument,
        support(&[CellId(10), CellId(20)]),
        &heat,
        300.0,
        || Ok(()),
    )
    .unwrap();
    let (later_heat, _, _) = observer
        .advance_candidate(&heat, 0.5, config(), || Ok(()))
        .unwrap();
    assert!(matches!(
        observer.advance_candidate(&later_heat, 1.0, config(), || Ok(())),
        Err(ThermalSensorError::Invalid(_))
    ));
    let mut stale_temperatures = std::collections::BTreeMap::new();
    stale_temperatures.insert(CellId(10), 299.0);
    stale_temperatures.insert(CellId(20), 350.0);
    let mut candidate = observer.clone();
    assert!(
        candidate
            .step(0.0, 0.5, &stale_temperatures, &later_heat, &mut || Ok(()))
            .is_err()
    );
    assert_eq!(candidate, observer);
    assert!(
        observer
            .observe_endpoint_candidate(&heat, &later_heat, 0.25, || Ok(()))
            .is_err()
    );
    assert!(
        observer
            .observe_endpoint_candidate(&later_heat, &later_heat, 0.5, || Ok(()))
            .is_err()
    );
}

#[test]
fn late_cancellation_rolls_back_sensor_and_heat_and_retry_replays() {
    let heat = initial(true);
    let instrument = sensor(
        SensorKind::Thermocouple,
        ObservationSupport::point(1, 0).unwrap(),
        Some(0.2),
        1.0,
        0.0,
    );
    let observer = TemperatureSensorObserver::new(
        &instrument,
        support(&[CellId(10)]),
        &heat,
        299.0,
        || Ok(()),
    )
    .unwrap();
    let mut calls = 0;
    let expected = observer
        .advance_candidate(&heat, 1.0, config(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    let rejected = observer.advance_candidate(&heat, 1.0, config(), || {
        seen += 1;
        if seen == calls {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        rejected,
        Err(ThermalSensorError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(observer.state().time_s(), 0.0);
    assert_eq!(observer.state().reading(), 299.0);
    assert_eq!(heat.time_s(), 0.0);
    let retried = observer
        .advance_candidate(&heat, 1.0, config(), || Ok(()))
        .unwrap();
    assert_eq!(retried.0.cells(), expected.0.cells());
    assert_eq!(retried.1, expected.1);
    assert_eq!(retried.2, expected.2);
}

struct LatentLaw;
impl ThermalLaw for LatentLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || h <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 3000.0 {
            h / 10.0
        } else if h <= 4000.0 {
            300.0
        } else {
            (h - 1000.0) / 10.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.1
    }
}

#[test]
fn latent_heat_uses_original_caloric_temperature_instead_of_energy_as_temperature() {
    let law = LatentLaw;
    let heat = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(10, 300.0, 2.0, 200.0)),
        }],
        &[(MaterialId(1), &law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0;
    let instrument = sensor(
        SensorKind::Thermocouple,
        ObservationSupport::point(1, 0).unwrap(),
        Some(0.25),
        1.0,
        0.0,
    );
    let observer = TemperatureSensorObserver::new(
        &instrument,
        support(&[CellId(10)]),
        &heat,
        290.0,
        || Ok(()),
    )
    .unwrap();
    let (next, sensor, audit) = observer
        .advance_candidate(&heat, 1.0, config(), || Ok(()))
        .unwrap();
    assert_eq!(next.temperature_k(CellId(10)).unwrap(), 300.0);
    assert_eq!(next.cells()[&CellId(10)].specific_h_j_kg, 3100.0);
    near(audit.external_heat_j, 200.0);
    near(
        sensor.state().reading(),
        response(290.0, 300.0, 0.0, 1.0, 0.25),
    );
}
