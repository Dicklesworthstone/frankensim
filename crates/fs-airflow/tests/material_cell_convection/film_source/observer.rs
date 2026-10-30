//! Material convection shares accepted native substeps with conduction/events.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::observer::SourceConvectionObserver;
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_thermal::observers::CoupledObservers;

struct ConstantConductivity(f64);
impl ConductivityLaw for ConstantConductivity {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.0)
    }
}

fn controls(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}

// Independent explicit finite-body balance. The convection formula is the
// closed-form material/Churchill-Chu reference in the parent fixture.
fn reference_step(t: &mut [f64; 2], dt: f64, area: f64, conductance: f64) -> f64 {
    let q =
        (*t).map(|wall| area * expected_h(0.5, 0.01, wall, 0.5 * wall + 150.0) * (300.0 - wall));
    let exchange = conductance * (t[1] - t[0]);
    t[0] += dt * (q[0] + exchange) / (MASS * CP);
    t[1] += dt * (q[1] - exchange) / (MASS * CP);
    dt * (q[0] + q[1])
}

#[test]
fn g1_same_substeps_drive_material_convection_and_native_spatial_conduction() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let saved = thermal.cells().clone();
    let card = liquid_card(0.5, 0.01, 320.0);
    let law = liquid_law(&card)
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
        .unwrap();
    let bindings = IDS.map(film);
    let convection =
        SourceConvectionObserver::new(&bindings, 2, |_: &SourceConvectiveCell, t: Temperature| {
            law.at_temperature(t)
        });
    let conductivity = ConstantConductivity(0.5);
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), &conductivity as &dyn ConductivityLaw)],
        &[ConductiveLink {
            link_id: 5,
            half_distances_m: [0.01; 2],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    let observer = CoupledObservers::new(convection, conduction);
    let (next, accepted, audit) = thermal
        .advance_observed(&observer, 0.04, controls(0.02), || Ok(()))
        .unwrap();
    let mut expected = [310.0, 290.0];
    let g = 0.002 / (0.01 / 0.5 + 0.01 / 0.5);
    let q1 = reference_step(&mut expected, 0.02, 0.008, g);
    let q2 = reference_step(&mut expected, 0.02, 0.008, g);
    for (id, temperature) in IDS.into_iter().zip(expected) {
        near(next.temperature_k(id).unwrap(), temperature, 1e-9);
    }
    near(audit.external_heat_j, q1 + q2, 1e-9);
    near(
        accepted.first.external_heat_j(),
        audit.external_heat_j,
        1e-9,
    );
    near(
        next.energy_j() - thermal.energy_j(),
        audit.external_heat_j,
        1e-9,
    );
    let source = accepted.first.sample().unwrap();
    for sample in source.materials() {
        near(
            sample.property_temperature_k,
            0.5 * next.temperature_k(sample.cell).unwrap() + 150.0,
            1e-12,
        );
        assert_eq!(
            sample.material.material().card_identity(),
            card.content_hash()
        );
    }
    for sample in source.convection().samples() {
        near(
            sample.wall_temperature_k,
            next.temperature_k(sample.cell).unwrap(),
            1e-12,
        );
        near(sample.exposed_area_m2, 0.008, 1e-15);
    }
    // A second advance owns only its own interval heat, with the same sources.
    let (again, continued, audit2) = next
        .advance_observed(&accepted, 0.06, controls(0.02), || Ok(()))
        .unwrap();
    let q3 = reference_step(&mut expected, 0.02, 0.008, g);
    near(continued.first.external_heat_j(), q3, 1e-9);
    near(audit2.external_heat_j, q3, 1e-9);
    near(again.temperature_k(IDS[0]).unwrap(), expected[0], 1e-9);
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.first.sample().is_none());
    assert_eq!(observer.first.external_heat_j(), 0.0);
    assert_eq!(next.links()[&5].conductance_w_k, 0.0);
}

fn scheduled<'a>(
    chart: &'a EquilibriumEnthalpyPhaseCurve,
    heater_w: f64,
    separation_s: Option<f64>,
) -> DepositionSimulation<'a> {
    let base = initial(chart, [310.0, 290.0], 0.0, 0.002);
    let mut events = base
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            if cell.id == IDS[0] {
                cell.boundary.source_w = heater_w;
            }
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .chain(base.links().values().copied().map(|link| TimedEvent {
            time_s: 0.0,
            event: Event::Contact(link),
        }))
        .collect::<Vec<_>>();
    if let Some(time_s) = separation_s {
        events.push(TimedEvent {
            time_s,
            event: Event::Separation(5),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), chart as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 4,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(0.02), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn g1_contact_events_change_the_exposed_area_at_the_original_event_time() {
    let chart = curve();
    let thermal = scheduled(&chart, 0.0, Some(0.015));
    let card = liquid_card(0.5, 0.01, 320.0);
    let law = liquid_law(&card)
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
        .unwrap();
    let bindings = IDS.map(film);
    let observer =
        SourceConvectionObserver::new(&bindings, 2, |_: &SourceConvectiveCell, t: Temperature| {
            law.at_temperature(t)
        });
    let (next, accepted, audit) = thermal
        .advance_observed(&observer, 0.04, controls(0.02), || Ok(()))
        .unwrap();
    let mut expected = [310.0, 290.0];
    let mut heat = 0.0;
    for (dt, area) in [(0.015, 0.008), (0.02, AREA), (0.005, AREA)] {
        heat += reference_step(&mut expected, dt, area, 0.0);
    }
    for (id, temperature) in IDS.into_iter().zip(expected) {
        near(next.temperature_k(id).unwrap(), temperature, 1e-9);
    }
    near(accepted.external_heat_j(), heat, 1e-9);
    near(audit.external_heat_j, heat, 1e-9);
    assert!(!next.links().contains_key(&5));
    for sample in accepted.sample().unwrap().convection().samples() {
        near(sample.exposed_area_m2, AREA, 1e-15);
    }
    let (_, at_event, _) = thermal
        .advance_observed(&observer, 0.015, controls(0.02), || Ok(()))
        .unwrap();
    assert!(at_event.sample().is_none());
    near(
        at_event.external_heat_j(),
        reference_step(&mut [310.0, 290.0], 0.015, 0.008, 0.0),
        1e-9,
    );
    assert!(thermal.links().contains_key(&5));
}

#[test]
fn g0_endpoint_source_refusal_and_late_cancellation_publish_neither_history() {
    let chart = curve();
    let thermal = scheduled(&chart, 500.0, None);
    let saved = thermal.cells().clone();
    // Initial film is 305 K; the heater takes the candidate above the 306 K
    // material support while its native caloric chart still admits the body.
    let limited_card = liquid_card(0.5, 0.01, 306.0);
    let limited = liquid_law(&limited_card);
    let forced = [SourceConvectiveCell {
        law: SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::FlatPlateLaminarAverage,
            speed_m_s: 1.0,
            settings: ForcedConvectionSettings::default(),
        },
        ..film(IDS[0])
    }];
    let observer =
        SourceConvectionObserver::new(&forced, 1, |_: &SourceConvectiveCell, t: Temperature| {
            limited.at_temperature(t)
        });
    assert!(matches!(
        thermal.advance_observed(&observer, 0.01, controls(0.01), || Ok(())),
        Err(SourceCellConvectionError::Source {
            cell: CellId(30),
            ..
        })
    ));
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none());
    assert_eq!(observer.external_heat_j(), 0.0);

    let card = liquid_card(0.5, 0.01, 320.0);
    let law = liquid_law(&card);
    let observer =
        SourceConvectionObserver::new(&forced, 1, |_: &SourceConvectiveCell, t: Temperature| {
            law.at_temperature(t)
        });
    let mut polls = 0;
    thermal
        .advance_observed(&observer, 0.01, controls(0.01), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut at = 0;
    assert!(matches!(
        thermal.advance_observed(&observer, 0.01, controls(0.01), || {
            at += 1;
            if at == polls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Cancelled
        ))
    ));
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none());
    assert_eq!(observer.external_heat_j(), 0.0);

    // Preserve rich native budget refusals even inside the Boolean polling seam.
    let mut trial = observer.clone();
    let mut calls = 0;
    assert!(matches!(
        trial.sources(&thermal, &mut || {
            calls += 1;
            if calls == 3 {
                Err(DepositionError::Budget("convection checkpoint"))
            } else {
                Ok(())
            }
        }),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Thermal(DepositionError::Budget("convection checkpoint"))
        ))
    ));
}
