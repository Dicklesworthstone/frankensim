//! Original bulk/wall liquid queries and reciprocal turbulent tube heat.
use super::*;

fn tube_link() -> SourceConvectiveLink {
    SourceConvectiveLink {
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::SiederTateTurbulent,
            speed_m_s: 1.0,
            settings: ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(100.0),
                ..Default::default()
            },
        },
        ..link()
    }
}

// Independent authored-property evaluation, not a call back into the card.
fn coefficient(viscosity_scale: f64, wall: f64, bulk: f64) -> f64 {
    let mu = viscosity_scale * 0.001 * (1.0 - 0.01 * (bulk - 300.0));
    let wall_mu = viscosity_scale * 0.001 * (1.0 - 0.01 * (wall - 300.0));
    let k = 0.5 + 0.002 * (bulk - 300.0);
    let re = 1000.0 * 0.1 / mu;
    let pr = mu * FLUID_CP / k;
    0.027 * re.powf(0.8) * pr.cbrt() * (mu / wall_mu).powf(0.14) * k / 0.1
}

fn increments(step: &BackwardEulerStep<'_, '_>, temperatures: [f64; 2]) -> Vec<f64> {
    step.cell_ids()
        .iter()
        .map(|&id| {
            if id == IDS[0] {
                MASS * CP * (temperatures[0] - 310.0)
            } else {
                MASS * FLUID_CP * (temperatures[1] - 290.0)
            }
        })
        .collect()
}

#[test]
fn g1_tube_heating_and_cooling_use_bulk_properties_and_actual_wall_viscosity() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let law = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    let original = thermal.cells().clone();
    for [wall, bulk] in [[310.0, 290.0], [290.0, 310.0]] {
        let step = BackwardEulerStep::new(&thermal, 0.1, || Ok(())).unwrap();
        let mut queried = Vec::new();
        let trial = step
            .evaluate(
                &increments(&step, [wall, bulk]),
                |state| {
                    let result = evaluate_source_link_convection(
                        state,
                        &[tube_link()],
                        1,
                        |binding, _, temperature| {
                            assert_eq!(binding.fluid_cell, IDS[1]);
                            queried.push(temperature.value());
                            law.at_temperature(temperature)
                        },
                        || false,
                    )?;
                    Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: result.conductances_w_k().clone(),
                        payload: result,
                    })
                },
                || Ok(()),
            )
            .unwrap();
        assert_eq!(queried, [bulk, wall]);
        let sample = &trial.payload().samples()[0];
        near(sample.property_temperature_k, bulk, 1e-12);
        let surface = sample.surface.as_ref().unwrap();
        near(surface.thermal.temperature_k(), wall, 1e-12);
        assert_eq!(
            sample.material.material().card_identity(),
            card.content_hash()
        );
        assert_eq!(
            surface.material.material().card_identity(),
            card.content_hash()
        );
        let groups = sample.evaluation.as_ref().unwrap().nusselt.groups();
        let ratio = (1.0 - 0.01 * (bulk - 300.0)) / (1.0 - 0.01 * (wall - 300.0));
        near(groups["mu_bulk_over_surface"], ratio, 1e-13);
        assert_eq!(ratio > 1.0, wall > bulk);
        let expected = LINK_AREA * coefficient(1.0, wall, bulk);
        near(sample.conductance_w_k, expected, 1e-9);
        near(sample.heat_into_fluid_w, expected * (wall - bulk), 1e-8);
        // A film tuple or a unit wall correction must visibly change h.
        assert!((expected - LINK_AREA * coefficient(1.0, 300.0, 300.0)).abs() > 1e-3);
        assert!((expected - LINK_AREA * coefficient(1.0, bulk, bulk)).abs() > 1e-3);
    }
    assert_eq!(thermal.cells(), &original);
}

fn endpoint(viscosity_scale: f64, dt: f64) -> [f64; 2] {
    let (mut lo, mut hi) = (294.0, 310.0);
    for _ in 0..80 {
        let wall = 0.5 * lo + 0.5 * hi;
        let bulk = 290.0 + CP / FLUID_CP * (310.0 - wall);
        let power = LINK_AREA * coefficient(viscosity_scale, wall, bulk) * (wall - bulk);
        if MASS * CP * (wall - 310.0) + dt * power > 0.0 {
            hi = wall;
        } else {
            lo = wall;
        }
    }
    let wall = 0.5 * lo + 0.5 * hi;
    [wall, 290.0 + CP / FLUID_CP * (310.0 - wall)]
}

#[test]
fn g1_g3_tube_link_closes_both_original_stores_and_responds_to_source_viscosity() {
    let chart = curve();
    let mut wall_endpoints = Vec::new();
    for scale in [1.0, 2.0] {
        let (card, caloric) = fluid_source(scale);
        let law = liquid_law(&card);
        let dt = 0.2;
        let expected = endpoint(scale, dt);
        for reversed in [false, true] {
            let thermal = pair(&chart, &caloric, reversed, 0.0);
            let original = thermal.cells().clone();
            let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
            let trial = step
                .evaluate(
                    &increments(&step, expected),
                    |state| {
                        let source = evaluate_source_link_convection(
                            state,
                            &[tube_link()],
                            1,
                            |_, _, temperature| law.at_temperature(temperature),
                            || false,
                        )?;
                        Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(
                            StepCoefficients {
                                sources_w: BTreeMap::new(),
                                conductances_w_k: source.conductances_w_k().clone(),
                                payload: source,
                            },
                        )
                    },
                    || Ok(()),
                )
                .unwrap();
            let (next, source, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
            let sample = &source.samples()[0];
            for (id, temperature) in IDS.into_iter().zip(expected) {
                near(next.temperature_k(id).unwrap(), temperature, 1e-9);
            }
            near(sample.property_temperature_k, expected[1], 1e-12);
            near(
                sample.surface.as_ref().unwrap().thermal.temperature_k(),
                expected[0],
                1e-12,
            );
            near(
                sample.heat_into_fluid_w,
                MASS * FLUID_CP * (expected[1] - 290.0) / dt,
                1e-8,
            );
            near(
                MASS * CP * (expected[0] - 310.0) + dt * sample.heat_into_fluid_w,
                0.0,
                1e-8,
            );
            near(next.energy_j(), thermal.energy_j(), 1e-8);
            near(audit.external_heat_j, 0.0, 1e-12);
            assert_eq!(thermal.cells(), &original);
            assert_eq!(next.links()[&5].conductance_w_k, 0.0);
        }
        wall_endpoints.push(expected[0]);
    }
    assert!(wall_endpoints[1] > wall_endpoints[0] + 0.1);
}

#[test]
fn g0_g4_tube_rejects_film_substitution_stale_wall_data_and_late_cancellation() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let law = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    let original = thermal.cells().clone();
    let calls = std::cell::Cell::new(0);
    let film = SourceConvectiveLink {
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        ..tube_link()
    };
    assert!(
        evaluate_source_link_convection(
            &thermal,
            &[film],
            1,
            |_, _, temperature| {
                calls.set(calls.get() + 1);
                law.at_temperature(temperature)
            },
            || false
        )
        .is_err()
    );
    assert_eq!(calls.get(), 0);
    assert!(
        evaluate_source_link_convection(
            &thermal,
            &[tube_link()],
            1,
            |_, _, _| {
                calls.set(calls.get() + 1);
                law.at_temperature(Temperature::new(290.0))
            },
            || false
        )
        .is_err()
    );
    assert_eq!(
        calls.get(),
        2,
        "the wall query cannot reuse a stale bulk snapshot"
    );
    calls.set(0);
    assert!(matches!(
        evaluate_source_link_convection(
            &thermal,
            &[tube_link()],
            1,
            |_, _, temperature| {
                calls.set(calls.get() + 1);
                law.at_temperature(temperature)
            },
            || calls.get() == 2
        ),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Cancelled
        ))
    ));
    assert_eq!(calls.get(), 2);
    assert_eq!(thermal.cells(), &original);
    assert_eq!(thermal.time_s(), 0.0);
}
