//! Inclined gap physics on the same finite-wall and source-state owners.

use super::*;

#[test]
fn g1_g3_inclination_changes_gap_heat_and_retains_link_orientation_and_source() {
    let chart = gas_source::body_curve();
    let model = gas_source::model(1.0);
    for reversed in [false, true] {
        let thermal = pair(&chart, [410.0, 390.0], reversed, 0.0);
        for angle in [0.0, 15.0, 30.0, 45.0, 60.0] {
            let binding = SourceEnclosedGasLayer {
                inclination_degrees: angle,
                ..binding()
            };
            let result = evaluate_source_enclosed_gas_layers(
                &thermal,
                &[binding],
                1,
                |_, _, temperature| {
                    model.evaluate(&gas_source::point(temperature.value(), PRESSURE))
                },
                || false,
            )
            .unwrap();
            let sample = &result.samples()[0];
            let (ra, nu, conductance) = expected_at_inclination(&model, 410.0, 390.0, angle);
            near(sample.evaluation.nusselt.groups()["Ra"], ra, 1e-9);
            near(sample.evaluation.nusselt.evidence().value, nu, 1e-12);
            near(sample.conductance_w_k, conductance, 1e-12);
            near(sample.heat_into_upper_w, 20.0 * conductance, 1e-11);
            assert_eq!(sample.inclination_degrees, angle);
            assert_eq!((sample.lower_cell, sample.upper_cell), (IDS[0], IDS[1]));
            assert_eq!(sample.material.source_packs(), model.source_packs());
            assert_eq!(sample.material.state().temperature, 400.0);
            assert_eq!(sample.material.state().pressure, PRESSURE);
            let expected_card = if angle == 0.0 {
                CorrelationId::HollandsHorizontalAirLayer
            } else {
                CorrelationId::HollandsInclinedAirLayer
            };
            assert_eq!(sample.evaluation.nusselt.card().id, expected_card);
            if angle > 0.0 {
                assert_eq!(
                    sample.evaluation.nusselt.groups()["inclination_degrees"],
                    angle
                );
            }
        }
    }
}

#[test]
fn g3_inclined_observer_transfers_heat_once_and_refreshes_actual_film_state() {
    let chart = gas_source::body_curve();
    let model = gas_source::model_transport(1.0, true);
    let thermal = pair(&chart, [410.0, 390.0], false, 0.0);
    let saved = thermal.cells().clone();
    let layers = [SourceEnclosedGasLayer {
        inclination_degrees: 45.0,
        ..binding()
    }];
    let observer = SourceEnclosedGasLayerObserver::new(
        &layers,
        1,
        |_: &SourceEnclosedGasLayer, _: &DepositionSimulation<'_>, temperature: Temperature| {
            model.evaluate(&gas_source::point(temperature.value(), PRESSURE))
        },
    );
    let (next, accepted, audit) = thermal
        .advance_observed(&observer, 0.5, controls(0.01), || Ok(()))
        .unwrap();
    assert!(next.temperature_k(IDS[0]).unwrap() < 410.0);
    assert!(next.temperature_k(IDS[1]).unwrap() > 390.0);
    near(next.energy_j(), thermal.energy_j(), 1e-8);
    near(audit.external_heat_j, 0.0, 1e-12);
    let sample = &accepted.sample().unwrap().samples()[0];
    let applied = &accepted.applied_sample().unwrap().samples()[0];
    let lower = next.temperature_k(IDS[0]).unwrap();
    let upper = next.temperature_k(IDS[1]).unwrap();
    near(
        sample.property_temperature_k,
        0.5 * lower + 0.5 * upper,
        1e-12,
    );
    near(
        sample.conductance_w_k,
        expected_at_inclination(&model, lower, upper, 45.0).2,
        1e-12,
    );
    assert_ne!(
        sample.property_temperature_k,
        applied.property_temperature_k
    );
    assert_eq!(
        sample.material.state().temperature,
        sample.property_temperature_k
    );
    assert_eq!(sample.inclination_degrees, 45.0);
    assert_eq!(next.links()[&5].conductance_w_k, 0.0);
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none());
}

#[test]
fn g0_inclined_geometry_and_unsupported_heated_above_states_refuse() {
    let chart = gas_source::body_curve();
    let model = gas_source::model(1.0);
    let thermal = pair(&chart, [410.0, 390.0], false, 0.0);
    for angle in [-1.0, 60.0_f64.next_up(), 90.0, f64::NAN, f64::INFINITY] {
        let mut queries = 0;
        let result = evaluate_source_enclosed_gas_layers(
            &thermal,
            &[SourceEnclosedGasLayer {
                inclination_degrees: angle,
                ..binding()
            }],
            1,
            |_, _, temperature| {
                queries += 1;
                model.evaluate(&gas_source::point(temperature.value(), PRESSURE))
            },
            || false,
        );
        assert!(result.is_err());
        assert_eq!(
            queries, 0,
            "invalid orientation must refuse before source queries"
        );
    }
    let tilted = [SourceEnclosedGasLayer {
        inclination_degrees: 45.0,
        ..binding()
    }];
    let reversed_heat = pair(&chart, [390.0, 410.0], false, 0.0);
    let saved = reversed_heat.cells().clone();
    assert!(
        evaluate_source_enclosed_gas_layers(
            &reversed_heat,
            &tilted,
            1,
            |_, _, temperature| model.evaluate(&gas_source::point(temperature.value(), PRESSURE)),
            || false,
        )
        .is_err()
    );
    assert_eq!(reversed_heat.cells(), &saved);

    // Exact equilibrium still retains the real finite Fourier coupling.
    let resting = pair(&chart, [400.0, 400.0], false, 0.0);
    let result = evaluate_source_enclosed_gas_layers(
        &resting,
        &tilted,
        1,
        |_, _, temperature| model.evaluate(&gas_source::point(temperature.value(), PRESSURE)),
        || false,
    )
    .unwrap();
    let sample = &result.samples()[0];
    assert_eq!(sample.evaluation.nusselt.evidence().value, 1.0);
    near(sample.conductance_w_k, GAP_AREA * 0.026 / GAP, 1e-14);
    assert_eq!(sample.heat_into_upper_w, 0.0);
}
