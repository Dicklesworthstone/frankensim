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

fn developing_link(length_over_diameter: f64) -> SourceConvectiveLink {
    SourceConvectiveLink {
        law: SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::GnielinskiDeveloping,
            speed_m_s: 1.0,
            settings: ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(length_over_diameter),
                ..Default::default()
            },
        },
        ..tube_link()
    }
}

fn section_link(start_m: f64, end_m: f64) -> SourceConvectiveLink {
    SourceConvectiveLink {
        law: SourceCellConvectionLaw::ForcedTubeSection {
            speed_m_s: 1.0,
            start: Length::new(start_m),
            end: Length::new(end_m),
        },
        ..tube_link()
    }
}

// Independent difference of cumulative entrance transfer, moderate intervals.
fn section_coefficient(bulk: f64, start_m: f64, end_m: f64) -> f64 {
    let developed = developing_coefficient(bulk, 1e6) / (1.0 + 1e6_f64.powf(-2.0 / 3.0));
    let a = start_m / 0.1;
    let b = end_m / 0.1;
    developed * ((b + b.cbrt()) - (a + a.cbrt())) / (b - a)
}

#[test]
fn downstream_tube_sections_close_finite_body_energy_with_the_actual_bulk_source() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let law = liquid_law(&card);
    let dt = 0.2;
    let mut final_walls = Vec::new();
    for start_m in [0.0, 1.0] {
        let (mut lo, mut hi) = (294.0, 310.0);
        for _ in 0..80 {
            let wall = 0.5 * lo + 0.5 * hi;
            let bulk = 290.0 + CP / FLUID_CP * (310.0 - wall);
            let heat = dt * LINK_AREA * section_coefficient(bulk, start_m, 2.0) * (wall - bulk);
            if MASS * CP * (wall - 310.0) + heat > 0.0 { hi = wall; }
            else { lo = wall; }
        }
        let wall = 0.5 * lo + 0.5 * hi;
        let expected = [wall, 290.0 + CP / FLUID_CP * (310.0 - wall)];
        for reversed in [false, true] {
            let thermal = pair(&chart, &caloric, reversed, 0.0);
            let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
            let mut queried = Vec::new();
            let trial = step.evaluate(&increments(&step, expected), |state| {
                let sources = evaluate_source_link_convection(
                    state, &[section_link(start_m, 2.0)], 1,
                    |_, _, temperature| {
                        queried.push(temperature.value());
                        law.at_temperature(temperature)
                    }, || false,
                )?;
                Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: sources.conductances_w_k().clone(), payload: sources,
                })
            }, || Ok(())).unwrap();
            let (next, source, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
            assert_eq!(queried.len(), 1);
            near(queried[0], expected[1], 1e-12);
            let sample = &source.samples()[0];
            assert!(sample.surface.is_none());
            let evaluation = sample.evaluation.as_ref().unwrap();
            assert!(evaluation.section.is_none());
            let details = evaluation.tube_section.as_ref().unwrap();
            assert_eq!(details.diameter.value(), 0.1);
            assert_eq!(details.start.value(), start_m);
            assert_eq!(details.end.value(), 2.0);
            assert_eq!(details.upstream.is_some(), start_m > 0.0);
            assert_eq!(evaluation.effective_nusselt(), details.nusselt);
            assert_eq!(evaluation.coefficient.model.discrepancy_rel.is_infinite(), start_m > 0.0);
            near(sample.conductance_w_k,
                LINK_AREA * section_coefficient(expected[1], start_m, 2.0), 1e-9);
            near(next.temperature_k(IDS[0]).unwrap(), expected[0], 1e-9);
            near(next.temperature_k(IDS[1]).unwrap(), expected[1], 1e-9);
            near(next.energy_j(), thermal.energy_j(), 1e-8);
            near(audit.external_heat_j, 0.0, 1e-12);
            near(audit.residual_j, 0.0, 1e-8);
            near(MASS * CP * (expected[0] - 310.0) + dt * sample.heat_into_fluid_w,
                0.0, 1e-8);
        }
        final_walls.push(wall);
    }
    assert!(final_walls[1] > final_walls[0] + 0.1,
        "a downstream section must not repeat the physical inlet enhancement");
}

#[test]
fn tube_section_geometry_and_property_convention_refuse_before_material_queries() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let law = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    let invalid = [
        section_link(-1.0, 2.0), section_link(2.0, 2.0), section_link(3.0, 2.0),
        section_link(0.0, f64::NAN),
        SourceConvectiveLink {
            characteristic_length: Length::new(0.0), ..section_link(1.0, 2.0)
        },
        SourceConvectiveLink {
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            ..section_link(1.0, 2.0)
        },
        SourceConvectiveLink {
            law: SourceCellConvectionLaw::ForcedTubeSection {
                speed_m_s: 0.0, start: Length::new(1.0), end: Length::new(2.0),
            },
            ..section_link(1.0, 2.0)
        },
    ];
    for binding in invalid {
        let mut calls = 0;
        assert!(evaluate_source_link_convection(&thermal, &[binding], 1, |_, _, t| {
            calls += 1;
            law.at_temperature(t)
        }, || false).is_err());
        assert_eq!(calls, 0);
    }
}

// Independent Gnielinski (2013) Eq. (2) with corrected Eq. (6), K=1.
fn developing_coefficient(bulk: f64, length_over_diameter: f64) -> f64 {
    let mu = 0.001 * (1.0 - 0.01 * (bulk - 300.0));
    let k = 0.5 + 0.002 * (bulk - 300.0);
    let re = 100.0 / mu;
    let pr = mu * FLUID_CP / k;
    let f = (1.8 * re.log10() - 1.5).powi(-2);
    let nu = (f / 8.0) * (re - 1000.0) * pr
        / (1.0 + 12.7 * (f / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0));
    nu * (1.0 + length_over_diameter.powf(-2.0 / 3.0)) * k / 0.1
}

#[test]
fn g1_g3_developing_tube_entrance_changes_closed_pair_cooling_without_creating_heat() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let law = liquid_law(&card);
    let dt = 0.2;
    let mut endpoints = Vec::new();
    for length_ratio in [10.0, 1e6] {
        let (mut lo, mut hi) = (294.0, 310.0);
        for _ in 0..80 {
            let wall = 0.5 * lo + 0.5 * hi;
            let bulk = 290.0 + CP / FLUID_CP * (310.0 - wall);
            let power = LINK_AREA * developing_coefficient(bulk, length_ratio) * (wall - bulk);
            if MASS * CP * (wall - 310.0) + dt * power > 0.0 { hi = wall; }
            else { lo = wall; }
        }
        let wall = 0.5 * lo + 0.5 * hi;
        let expected = [wall, 290.0 + CP / FLUID_CP * (310.0 - wall)];
        for reversed in [false, true] {
            let thermal = pair(&chart, &caloric, reversed, 0.0);
            let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
            let mut temperatures = Vec::new();
            let trial = step.evaluate(&increments(&step, expected), |state| {
                let source = evaluate_source_link_convection(
                    state, &[developing_link(length_ratio)], 1,
                    |_, _, temperature| {
                        temperatures.push(temperature.value());
                        law.at_temperature(temperature)
                    }, || false,
                )?;
                Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: source.conductances_w_k().clone(), payload: source,
                })
            }, || Ok(())).unwrap();
            let (next, source, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
            assert_eq!(temperatures.len(), 1, "only actual bulk properties enter K=1");
            near(temperatures[0], expected[1], 1e-12);
            let sample = &source.samples()[0];
            assert!(sample.surface.is_none());
            near(sample.conductance_w_k,
                LINK_AREA * developing_coefficient(expected[1], length_ratio), 1e-9);
            near(next.temperature_k(IDS[0]).unwrap(), expected[0], 1e-9);
            near(next.temperature_k(IDS[1]).unwrap(), expected[1], 1e-9);
            near(next.energy_j(), thermal.energy_j(), 1e-8);
            near(audit.external_heat_j, 0.0, 1e-12);
            near(MASS * CP * (expected[0] - 310.0) + dt * sample.heat_into_fluid_w,
                0.0, 1e-8);
            assert_eq!(thermal.time_s(), 0.0);
        }
        endpoints.push(wall);
    }
    assert!(endpoints[0] < endpoints[1] - 0.1,
        "the actual short-tube entrance enhancement must cool the finite wall more");
}

#[test]
fn g0_developing_tube_rejects_film_substitution_and_non_circular_geometry_before_queries() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let law = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    let mut no_length = developing_link(10.0);
    let mut rectangular = no_length;
    if let SourceCellConvectionLaw::Forced { settings, .. } = &mut no_length.law {
        settings.length_over_hydraulic_diameter = None;
    }
    if let SourceCellConvectionLaw::Forced { settings, .. } = &mut rectangular.law {
        settings.aspect_ratio = Some(0.5);
    }
    for binding in [SourceConvectiveLink {
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        ..developing_link(10.0)
    }, no_length, rectangular] {
        let mut calls = 0;
        let result = evaluate_source_link_convection(&thermal, &[binding], 1, |_, _, t| {
            calls += 1;
            law.at_temperature(t)
        }, || false);
        assert!(result.is_err());
        assert_eq!(calls, 0);
    }
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
