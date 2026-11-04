//! Vertical-cavity heat on the original finite walls, material sources and CFL.

use super::*;

// Wright (1996), equations 9 and 7, evaluated independently of fs-convection.
pub(super) fn expected_nusselt(rayleigh: f64, height_over_gap: f64) -> f64 {
    let tall = if rayleigh <= 1e4 {
        1.0 + 1.75967e-10 * rayleigh.powf(2.2984755)
    } else if rayleigh <= 5e4 {
        0.028154 * rayleigh.powf(0.4134)
    } else {
        0.0673838 * rayleigh.cbrt()
    };
    tall.max(0.242 * (rayleigh / height_over_gap).powf(0.272))
}

fn vertical() -> SourceEnclosedGasLayer {
    SourceEnclosedGasLayer {
        inclination_degrees: 90.0,
        ..binding()
    }
}

#[test]
fn g1_vertical_gap_uses_rayleigh_magnitude_and_signed_reciprocal_heat() {
    let chart = gas_source::body_curve();
    // Manufactured source packs exercise the air-policy API and source ownership;
    // they do not constitute empirical validation of this gas or correlation.
    let model = gas_source::model(1.0);
    for pressure in [100000.0, 200000.0, 500000.0] {
        let reference = model.evaluate(&gas_source::point(400.0, pressure)).unwrap();
        let s = reference.state();
        for reversed in [false, true] {
            for temperatures in [[410.0, 390.0], [390.0, 410.0], [400.0, 400.0]] {
                let thermal = pair(&chart, temperatures, reversed, 0.0);
                let saved = thermal.cells().clone();
                for reference_wall in [0, 1] {
                    let layer = SourceEnclosedGasLayer {
                        lower_cell: IDS[reference_wall],
                        ..vertical()
                    };
                    let mut queries = Vec::new();
                    let result = evaluate_source_enclosed_gas_layers(
                        &thermal,
                        &[layer],
                        1,
                        |_, _, temperature| {
                            queries.push(temperature.value());
                            model.evaluate(&gas_source::point(temperature.value(), pressure))
                        },
                        || false,
                    )
                    .unwrap();
                    let sample = &result.samples()[0];
                    let difference =
                        temperatures[reference_wall] - temperatures[1 - reference_wall];
                    let ra = 9.81 * difference.abs() / s.temperature
                        * GAP.powi(3)
                        * s.density.powi(2)
                        * s.specific_heat_cp
                        / (s.dynamic_viscosity * s.thermal_conductivity);
                    let nu = expected_nusselt(ra, SPANS[1] / GAP);
                    let conductance = GAP_AREA * s.thermal_conductivity / GAP * nu;
                    assert_eq!(queries, vec![400.0]);
                    assert_eq!(sample.lower_cell, IDS[reference_wall]);
                    assert_eq!(sample.upper_cell, IDS[1 - reference_wall]);
                    assert_eq!(sample.inclination_degrees, 90.0);
                    assert_eq!(
                        sample.evaluation.nusselt.card().id,
                        CorrelationId::WrightVerticalAirLayer
                    );
                    near(sample.evaluation.nusselt.groups()["Ra"], ra, 1e-7);
                    near(
                        sample.evaluation.nusselt.groups()["aspect_ratio"],
                        SPANS[1] / GAP,
                        1e-12,
                    );
                    near(sample.evaluation.nusselt.evidence().value, nu, 1e-12);
                    near(sample.conductance_w_k, conductance, 1e-12);
                    near(sample.heat_into_upper_w, conductance * difference, 1e-11);
                    assert_eq!(sample.material.source_packs(), model.source_packs());
                    assert_eq!(sample.material.state().pressure, pressure);
                    assert_eq!(sample.material.state().temperature, 400.0);
                    if difference == 0.0 {
                        assert_eq!(sample.evaluation.nusselt.evidence().value, 1.0);
                        assert_eq!(sample.heat_into_upper_w, 0.0);
                    }
                    assert_eq!(thermal.cells(), &saved);
                    assert_eq!(thermal.links()[&5].conductance_w_k, 0.0);
                }
            }
        }
    }
}

#[test]
fn g1_vertical_height_is_the_second_whole_span_and_short_cavities_are_admitted() {
    let chart = gas_source::body_curve();
    let model = gas_source::model(1.0);
    let thermal = pair(&chart, [410.0, 390.0], false, 0.0);
    let baseline_ra = expected(&model, 410.0, 390.0).0;
    let mut nusselt = Vec::new();
    for (gap_m, spans) in [
        (GAP, SPANS),
        (GAP, [SPANS[1], SPANS[0]]),
        (GAP, [GAP_AREA / (5.0 * GAP), 5.0 * GAP]),
        (0.005, [GAP_AREA / 0.55, 0.55]),
    ] {
        let layer = SourceEnclosedGasLayer {
            gap: Length::new(gap_m),
            plate_spans: spans.map(Length::new),
            ..vertical()
        };
        let result = evaluate_source_enclosed_gas_layers(
            &thermal,
            &[layer],
            1,
            |_, _, t| model.evaluate(&gas_source::point(t.value(), PRESSURE)),
            || false,
        )
        .unwrap();
        let sample = &result.samples()[0];
        let ra = baseline_ra * (gap_m / GAP).powi(3);
        let height_over_gap = spans[1] / gap_m;
        let nu = expected_nusselt(ra, height_over_gap);
        near(sample.evaluation.nusselt.groups()["Ra"], ra, 1e-9);
        near(
            sample.evaluation.nusselt.groups()["aspect_ratio"],
            height_over_gap,
            1e-12,
        );
        near(sample.evaluation.nusselt.evidence().value, nu, 1e-12);
        near(sample.conductance_w_k, GAP_AREA * 0.026 / gap_m * nu, 1e-12);
        nusselt.push(nu);
    }
    // Swapping width and height keeps the area and minimum span unchanged.
    // Only the true vertical height changes the finite-height correction.
    assert!(nusselt[1] > nusselt[0]);
    assert!(nusselt[2] > nusselt[0]);
}

#[test]
fn g3_vertical_observer_refreshes_film_and_conserves_both_wall_energies() {
    let chart = gas_source::body_curve();
    let model = gas_source::model_transport(1.0, true);
    let thermal = pair(&chart, [390.0, 410.0], false, 0.0);
    let saved = thermal.cells().clone();
    let layers = [vertical()];
    let observer = SourceEnclosedGasLayerObserver::new(
        &layers,
        1,
        |_: &SourceEnclosedGasLayer, _: &DepositionSimulation<'_>, t: Temperature| {
            model.evaluate(&gas_source::point(t.value(), PRESSURE))
        },
    );
    let (next, accepted, audit) = thermal
        .advance_observed(&observer, 0.5, controls(0.01), || Ok(()))
        .unwrap();
    let lower = next.temperature_k(IDS[0]).unwrap();
    let upper = next.temperature_k(IDS[1]).unwrap();
    assert!(lower > 390.0 && upper < 410.0);
    near(next.energy_j(), thermal.energy_j(), 1e-8);
    near(audit.external_heat_j, 0.0, 1e-12);
    let sample = &accepted.sample().unwrap().samples()[0];
    let applied = &accepted.applied_sample().unwrap().samples()[0];
    near(sample.property_temperature_k, 0.5 * lower + 0.5 * upper, 1e-12);
    near(
        sample.conductance_w_k,
        expected_at_inclination(&model, lower, upper, 90.0).2,
        1e-12,
    );
    assert!(sample.heat_into_upper_w < 0.0);
    assert_ne!(sample.property_temperature_k, applied.property_temperature_k);
    assert_eq!(sample.material.state().temperature, sample.property_temperature_k);
    assert_eq!(sample.material.source_packs(), model.source_packs());
    assert_eq!(next.links()[&5].conductance_w_k, 0.0);
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none() && observer.applied_sample().is_none());
}

#[test]
fn g0_vertical_geometry_and_material_limits_refuse_before_mutating_walls() {
    let chart = gas_source::body_curve();
    let model = gas_source::model(1.0);
    let thermal = pair(&chart, [410.0, 390.0], false, 0.0);
    let saved = thermal.cells().clone();
    for (gap_m, spans) in [
        (GAP, [10.0 * GAP, GAP_AREA / (10.0 * GAP)]),
        (GAP, [GAP_AREA / (4.9 * GAP), 4.9 * GAP]),
        (0.005, [GAP_AREA / 0.5505, 0.5505]),
    ] {
        let mut queries = 0;
        let result = evaluate_source_enclosed_gas_layers(
            &thermal,
            &[SourceEnclosedGasLayer {
                gap: Length::new(gap_m),
                plate_spans: spans.map(Length::new),
                ..vertical()
            }],
            1,
            |_, _, t| {
                queries += 1;
                model.evaluate(&gas_source::point(t.value(), PRESSURE))
            },
            || false,
        );
        assert!(result.is_err());
        assert_eq!(queries, 0);
        assert_eq!(thermal.cells(), &saved);
    }

    let foreign_pr = gas_source::model(0.75);
    for (source, pressure, layer) in [
        (&model, 2000000.0, vertical()),
        (&foreign_pr, PRESSURE, vertical()),
        (
            &model,
            PRESSURE,
            SourceEnclosedGasLayer {
                maximum_relative_density_departure: 0.01,
                ..vertical()
            },
        ),
    ] {
        assert!(
            evaluate_source_enclosed_gas_layers(
                &thermal,
                &[layer],
                1,
                |_, _, t| source.evaluate(&gas_source::point(t.value(), pressure)),
                || false,
            )
            .is_err()
        );
        assert_eq!(thermal.cells(), &saved);
    }
}
