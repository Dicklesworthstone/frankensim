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

// Independent subcycling: original conduction rate plus the separately bounded
// convection rate, with event/end clipping AFTER the composed rate restriction.
fn reference_advance(
    t: &mut [f64; 2], duration: f64, maximum_dt: f64, area: f64, conductance: f64,
) -> f64 {
    let capacity = MASS * CP;
    let native_dt = if conductance == 0.0 { maximum_dt }
        else { maximum_dt.min(0.8 * capacity / conductance) };
    let mut time = 0.0;
    let mut heat = 0.0;
    while time < duration {
        let extra_rate = t.iter().map(|&wall| {
            if wall == 300.0 { 0.0 }
            else { area * expected_h(0.5, 0.01, wall, 0.5 * wall + 150.0) / capacity }
        }).fold(0.0, f64::max);
        let dt = (1.0 / (1.0 / native_dt + extra_rate / 0.8)).min(duration - time);
        heat += reference_step(t, dt, area, conductance);
        time += dt;
    }
    heat
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
    let expected_heat = reference_advance(&mut expected, 0.04, 0.02, 0.008, g);
    for (id, temperature) in IDS.into_iter().zip(expected) {
        near(next.temperature_k(id).unwrap(), temperature, 1e-9);
    }
    near(audit.external_heat_j, expected_heat, 1e-9);
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
    let q3 = reference_advance(&mut expected, 0.02, 0.02, 0.008, g);
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
    for (duration, area) in [(0.015, 0.008), (0.025, AREA)] {
        heat += reference_advance(&mut expected, duration, 0.02, area, 0.0);
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
        reference_advance(&mut [310.0, 290.0], 0.015, 0.02, 0.008, 0.0),
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

#[test]
fn forced_fluid_substitution_changes_native_subcycling_without_thermal_overshoot() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let saved = thermal.cells().clone();
    let bindings = IDS.map(|cell| SourceConvectiveCell {
        law: SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::FlatPlateLaminarAverage,
            speed_m_s: 1.0,
            settings: ForcedConvectionSettings::default(),
        },
        ..film(cell)
    });
    let conductivity = ConstantConductivity(0.5);
    let conduction = ConductionObserver::new(
        &thermal,
        &[(MaterialId(0), &conductivity as &dyn ConductivityLaw)],
        &[ConductiveLink {
            link_id: 5, half_distances_m: [0.01; 2], contact_resistance_m2_k_w: 0.0,
        }], &[], || Ok(()),
    ).unwrap();
    let mut counts = Vec::new();
    for density in [250.0, 1000.0] {
        let fluid = FluidTransportState::try_new(density, 4000.0, 0.5, 0.001).unwrap();
        let convection = SourceConvectionObserver::new(
            &bindings, 2, move |_: &SourceConvectiveCell, t: Temperature| {
                ThermalFluidTransportState::try_new(fluid, t.value(), None)
            },
        );
        let observer = CoupledObservers::new(convection, conduction.clone());
        let (next, accepted, audit) = thermal.advance_observed(
            &observer, 0.2, controls(1.0), || Ok(()),
        ).unwrap();
        // Independent laminar-plate correlation and symmetric two-cell balance.
        let h = 0.664 * (density * 0.1 / 0.001).sqrt()
            * (0.001_f64 * 4000.0 / 0.5).cbrt() * 0.5 / 0.1;
        let g = 0.008 * h;
        let link_g = 0.002 / (0.01 / 0.5 + 0.01 / 0.5);
        let dt_bound = 1.0 / (1.0 + g / (0.8 * MASS * CP));
        let mut contrast = 10.0;
        let mut time = 0.0;
        let mut steps = 0;
        while time < 0.2 {
            let dt = dt_bound.min(0.2 - time);
            contrast *= 1.0 - dt * (g + 2.0 * link_g) / (MASS * CP);
            time += dt;
            steps += 1;
        }
        near(next.temperature_k(IDS[0]).unwrap(), 300.0 + contrast, 1e-10);
        near(next.temperature_k(IDS[1]).unwrap(), 300.0 - contrast, 1e-10);
        assert_eq!(audit.substeps, steps);
        assert!(audit.substeps > 1);
        assert!((300.0..310.0).contains(&next.temperature_k(IDS[0]).unwrap()));
        assert!((290.0..300.0).contains(&next.temperature_k(IDS[1]).unwrap()));
        near(accepted.first.external_heat_j(), audit.external_heat_j, 1e-10);
        near(next.energy_j() - thermal.energy_j(), audit.external_heat_j, 1e-9);
        assert!(matches!(
            thermal.advance_observed_limited(&observer, 0.2, controls(1.0), 1, || Ok(())),
            Err(fs_thermal::observers::ObserverCouplingError::Thermal(
                DepositionError::Budget("thermal substeps")
            ))
        ));
        counts.push(audit.substeps);
    }
    assert!(counts[1] > counts[0]);
    assert_eq!(thermal.cells(), &saved);
    assert_eq!(thermal.links()[&5].conductance_w_k, 0.0);
}

#[test]
fn convection_step_bound_requires_a_current_source_and_is_invalidated_by_events() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let card = liquid_card(0.5, 0.01, 320.0);
    let law = liquid_law(&card)
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0)).unwrap();
    let bindings = IDS.map(film);
    let mut observer = SourceConvectionObserver::new(
        &bindings, 2, |_: &SourceConvectiveCell, t: Temperature| law.at_temperature(t),
    );
    assert!(observer.limit_timestep_s(&thermal, 1.0, &mut || Ok(())).is_err());
    observer.sources(&thermal, &mut || Ok(())).unwrap();
    assert!(observer.limit_timestep_s(&thermal, 1.0, &mut || Ok(())).unwrap() < 1.0);
    observer.event(&Event::Separation(5), &thermal).unwrap();
    assert!(observer.limit_timestep_s(&thermal, 1.0, &mut || Ok(())).is_err());
    observer.sources(&thermal, &mut || Ok(())).unwrap();
    observer.begin().unwrap();
    assert!(observer.limit_timestep_s(&thermal, 1.0, &mut || Ok(())).is_err());
    assert_eq!(thermal.time_s(), 0.0);
}
