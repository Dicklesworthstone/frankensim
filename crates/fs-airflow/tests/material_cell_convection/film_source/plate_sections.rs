//! One physical leading edge across the original finite-cell heat owners.

use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::observer::SourceConvectionObserver;

const PLATE_LENGTH: f64 = 0.2;
const WIDTH: f64 = 0.05;

fn section(
    cell: CellId,
    start: f64,
    end: f64,
    speed_m_s: f64,
    correlation: CorrelationId,
    reservoir_k: f64,
) -> SourceConvectiveCell {
    SourceConvectiveCell {
        cell,
        reference_temperature: Temperature::new(reservoir_k),
        characteristic_length: Length::new(end),
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        law: SourceCellConvectionLaw::ForcedPlateSection {
            correlation,
            speed_m_s,
            start: Length::new(start),
        },
    }
}

fn plate<'a>(
    chart: &'a EquilibriumEnthalpyPhaseCurve,
    edges: &[f64],
    temperature_k: f64,
) -> DepositionSimulation<'a> {
    let events = edges.windows(2).enumerate().map(|(i, edge)| {
        let width = edge[1] - edge[0];
        let fraction = width / PLATE_LENGTH;
        TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(500 + i as u64),
                body_id: 0,
                segment_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: MASS * fraction,
                specific_h_j_kg: CP * (temperature_k - 300.0),
                geometry: CellGeometry {
                    start_m: [edge[0], 0.0, 0.0],
                    end_m: [edge[1], 0.0, 0.0],
                    volume_m3: 1e-6 * fraction,
                    surface_area_m2: WIDTH * width,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: None,
                },
            }),
        }
    }).collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), chart as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: edges.len() - 1,
            max_links: 1,
            max_events: edges.len() - 1,
            max_substeps: 100,
        },
    ).unwrap().advance_to(
        0.0,
        StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-9 },
        || Ok(()),
    ).unwrap().0
}

// Independent integral of h(x), using the actual source's physical properties.
// Positive-start turbulent strips cancel the shared 871 leading-edge term.
pub(super) fn independent_h(
    start: f64,
    end: f64,
    speed: f64,
    transport: FluidTransportState,
    correlation: CorrelationId,
) -> f64 {
    let reynolds_per_m = transport.density_kg_m3() * speed
        / transport.dynamic_viscosity_pa_s();
    let k = transport.thermal_conductivity_w_m_k();
    let pr = transport.prandtl();
    match correlation {
        CorrelationId::FlatPlateLaminarAverage => {
            0.664 * k * reynolds_per_m.sqrt() * pr.cbrt() / (end.sqrt() + start.sqrt())
        }
        CorrelationId::FlatPlateTurbulentAverage => {
            let entrance = if start == 0.0 { 871.0 } else { 0.0 };
            k * (0.037 * reynolds_per_m.powf(0.8) * (end.powf(0.8) - start.powf(0.8))
                - entrance) * pr.cbrt() / (end - start)
        }
        _ => panic!("reference requires a flat-plate card"),
    }
}

#[test]
fn unequal_subdivision_preserves_whole_plate_heat_and_both_original_mean_receipts() {
    let chart = curve();
    let card = liquid_card(0.5, 0.01, 320.0);
    let law = liquid_law(&card);
    let reference = law.at_temperature(Temperature::new(300.0)).unwrap();
    let transport = reference.transport_state().unwrap();
    for (correlation, speed, edges) in [
        (CorrelationId::FlatPlateLaminarAverage, 0.1, [0.0, 0.005, 0.03, 0.11, PLATE_LENGTH]),
        (CorrelationId::FlatPlateTurbulentAverage, 20.0, [0.0, 0.04, 0.08, 0.14, PLATE_LENGTH]),
    ] {
        for (wall_k, bulk_k) in [(310.0, 290.0), (290.0, 310.0)] {
            let thermal = plate(&chart, &edges, wall_k);
            let saved = thermal.cells().clone();
            let mut bindings = edges.windows(2).enumerate().map(|(i, edge)| {
                section(CellId(500 + i as u64), edge[0], edge[1], speed, correlation, bulk_k)
            }).collect::<Vec<_>>();
            bindings.reverse();
            let mut queries = Vec::new();
            let result = evaluate_source_cell_convection(
                &thermal, &bindings, bindings.len(),
                |binding, t| {
                    queries.push((binding.cell, t.value()));
                    law.at_temperature(t)
                },
                || false,
            ).unwrap();
            assert_eq!(queries, (0..4).map(|i| (CellId(500 + i), 300.0)).collect::<Vec<_>>());
            let whole_h = independent_h(0.0, PLATE_LENGTH, speed, transport, correlation);
            let expected_power = WIDTH * PLATE_LENGTH * whole_h * (bulk_k - wall_k);
            near(result.convection().total_heat_into_cells_w(), expected_power,
                2e-11 * expected_power.abs().max(1.0));
            let mut restarted_power = 0.0;
            for sample in result.convection().samples() {
                let i = (sample.cell.0 - 500) as usize;
                let (start, end) = (edges[i], edges[i + 1]);
                let evaluation = sample.evaluation.as_ref().unwrap();
                let details = evaluation.section.as_ref().unwrap();
                let h = independent_h(start, end, speed, transport, correlation);
                near(evaluation.coefficient.value.value(), h, 2e-11 * h);
                near(evaluation.effective_nusselt(), h * end
                    / transport.thermal_conductivity_w_m_k(), 1e-9);
                assert_eq!(details.start.value(), start);
                assert_eq!(details.end.value(), end);
                assert_eq!(evaluation.nusselt.card().id, correlation);
                near(evaluation.nusselt.groups()["Re"],
                    transport.density_kg_m3() * speed * end
                        / transport.dynamic_viscosity_pa_s(), 1e-8);
                if start == 0.0 {
                    assert!(details.upstream.is_none());
                    assert_eq!(evaluation.effective_nusselt(), evaluation.nusselt.evidence().value);
                } else {
                    let upstream = details.upstream.as_ref().unwrap();
                    assert_eq!(upstream.card().id, correlation);
                    near(upstream.groups()["Re"], transport.density_kg_m3() * speed * start
                        / transport.dynamic_viscosity_pa_s(), 1e-8);
                    assert!(evaluation.coefficient.model.discrepancy_rel.is_infinite());
                    assert_ne!(evaluation.effective_nusselt(), evaluation.nusselt.evidence().value);
                }
                restarted_power += WIDTH * (end - start)
                    * independent_h(0.0, end - start, speed, transport, correlation)
                    * (bulk_k - wall_k);
            }
            // A new entrance at every cell produces a materially different heat rate.
            assert!((restarted_power - expected_power).abs() > 0.01 * expected_power.abs());
            for material in result.materials() {
                assert_eq!(material.material.material().card_identity(), card.content_hash());
                assert_eq!(material.property_temperature_k, 300.0);
            }
            assert_eq!(thermal.cells(), &saved);
        }
    }
}

#[test]
fn plate_sections_refresh_actual_films_and_use_native_area_energy_and_timestep_bounds() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let saved = thermal.cells().clone();
    let card = liquid_card(0.5, 0.01, 320.0);
    let law = liquid_law(&card);
    let bindings = [
        section(IDS[0], 0.0, 0.1, 0.1, CorrelationId::FlatPlateLaminarAverage, 300.0),
        section(IDS[1], 0.1, 0.2, 0.1, CorrelationId::FlatPlateLaminarAverage, 300.0),
    ];
    let observer = SourceConvectionObserver::new(
        &bindings, 2, |_: &SourceConvectiveCell, t: Temperature| law.at_temperature(t),
    );
    let (next, accepted, audit) = thermal.advance_observed(
        &observer, 1.0,
        StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-8 },
        || Ok(()),
    ).unwrap();
    assert!(audit.substeps > 1, "the actual strip hA must tighten the native timestep");
    assert!(next.temperature_k(IDS[0]).unwrap() < 310.0);
    assert!(next.temperature_k(IDS[1]).unwrap() > 290.0);
    near(next.energy_j() - thermal.energy_j(), audit.external_heat_j, 1e-8);
    near(accepted.external_heat_j(), audit.external_heat_j, 1e-8);
    let result = accepted.sample().unwrap();
    for (sample, material) in result.convection().samples().iter().zip(result.materials()) {
        let binding = bindings.iter().find(|binding| binding.cell == sample.cell).unwrap();
        let wall = next.temperature_k(sample.cell).unwrap();
        near(material.property_temperature_k, 0.5 * wall + 150.0, 1e-12);
        assert_eq!(material.material.material().query().axes()["T"], material.property_temperature_k);
        assert_eq!(material.material.material().card_identity(), card.content_hash());
        near(sample.exposed_area_m2, 0.008, 1e-15);
        let evaluation = sample.evaluation.as_ref().unwrap();
        let details = evaluation.section.as_ref().unwrap();
        let h = independent_h(details.start.value(), binding.characteristic_length.value(),
            0.1, material.thermal.transport(), CorrelationId::FlatPlateLaminarAverage);
        near(evaluation.coefficient.value.value(), h, 2e-11 * h);
        near(sample.heat_into_cell_w, 0.008 * h * (300.0 - wall), 1e-9);
    }
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none());
    assert_eq!(next.links()[&5].conductance_w_k, 0.0);
}

#[test]
fn section_geometry_film_convention_and_both_endpoint_domains_remain_admitted() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let saved = thermal.cells().clone();
    let card = liquid_card(0.5, 0.01, 320.0);
    let law = liquid_law(&card);
    let good = section(IDS[0], 0.1, 0.2, 0.1,
        CorrelationId::FlatPlateLaminarAverage, 300.0);
    for binding in [
        section(IDS[0], -0.1, 0.2, 0.1, CorrelationId::FlatPlateLaminarAverage, 300.0),
        section(IDS[0], 0.2, 0.2, 0.1, CorrelationId::FlatPlateLaminarAverage, 300.0),
        section(IDS[0], f64::NAN, 0.2, 0.1, CorrelationId::FlatPlateLaminarAverage, 300.0),
        section(IDS[0], 1e-300, 1e300, 0.1, CorrelationId::FlatPlateLaminarAverage, 300.0),
        section(IDS[0], 0.1, 0.2, -0.1, CorrelationId::FlatPlateLaminarAverage, 300.0),
        SourceConvectiveCell { property_temperature: ConvectionPropertyTemperature::Reservoir, ..good },
    ] {
        let mut queries = 0;
        assert!(evaluate_source_cell_convection(&thermal, &[binding], 1,
            |_, t| { queries += 1; law.at_temperature(t) }, || false).is_err());
        assert_eq!(queries, 0);
    }
    // Admitted downstream Re cannot authorize an unsupported positive start.
    for binding in [
        section(IDS[0], 1e-10, 0.2, 0.1, CorrelationId::FlatPlateLaminarAverage, 300.0),
        section(IDS[0], 0.001, 0.2, 20.0, CorrelationId::FlatPlateTurbulentAverage, 300.0),
    ] {
        let mut queries = 0;
        assert!(evaluate_source_cell_convection(&thermal, &[binding], 1,
            |_, t| { queries += 1; law.at_temperature(t) }, || false).is_err());
        assert_eq!(queries, 1);
    }
    let queries = std::cell::Cell::new(0);
    let bindings = [good, SourceConvectiveCell { cell: IDS[1], ..good }];
    assert!(matches!(
        evaluate_source_cell_convection(&thermal, &bindings, 2,
            |_, t| { queries.set(queries.get() + 1); law.at_temperature(t) },
            || queries.get() == 1),
        Err(SourceCellConvectionError::Convection(ConvectionSourceError::Cancelled))
    ));
    assert_eq!(thermal.cells(), &saved);
}
