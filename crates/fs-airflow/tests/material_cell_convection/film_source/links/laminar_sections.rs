//! Physical laminar entrances reach the original finite solid/coolant stores.
use super::*;

fn circular(start: f64, end: f64) -> SourceConvectiveLink {
    SourceConvectiveLink {
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: SourceCellConvectionLaw::ForcedLaminarTubeSection {
            speed_m_s: 0.01,
            start: Length::new(start),
            end: Length::new(end),
        },
        ..link()
    }
}

// Independent Hausen cumulative transfer at moderate endpoint separation.
fn circular_h(bulk: f64, start: f64, end: f64) -> f64 {
    let k = 0.5 + 0.002 * (bulk - 300.0);
    let peclet = 0.01 * 0.1 * 1000.0 * FLUID_CP / k;
    let cumulative = |x: f64| {
        if x == 0.0 {
            return 0.0;
        }
        let gz = peclet / (x / 0.1);
        x * (3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0)))
    };
    (cumulative(end) - cumulative(start)) / (end - start) * k / 0.1
}

fn energy_increments(
    step: &BackwardEulerStep<'_, '_>,
    temperatures: [f64; 2],
    fluid_cp: f64,
) -> Vec<f64> {
    step.cell_ids()
        .iter()
        .map(|&id| {
            MASS * if id == IDS[0] {
                CP * (temperatures[0] - 310.0)
            } else {
                fluid_cp * (temperatures[1] - 290.0)
            }
        })
        .collect()
}

#[test]
fn laminar_entrance_position_changes_closed_pair_cooling_without_restarting_flow() {
    let wall_curve = curve();
    let (card, caloric) = fluid_source(1.0);
    let liquid = liquid_law(&card);
    let dt = 2.0;
    let mut endpoints = Vec::new();
    for start in [0.0, 1.0] {
        let (mut lo, mut hi) = (294.0, 310.0);
        for _ in 0..80 {
            let wall = 0.5 * (lo + hi);
            let bulk = 290.0 + CP / FLUID_CP * (310.0 - wall);
            let heat = dt * LINK_AREA * circular_h(bulk, start, 2.0) * (wall - bulk);
            if MASS * CP * (wall - 310.0) + heat > 0.0 {
                hi = wall;
            } else {
                lo = wall;
            }
        }
        let wall = 0.5 * (lo + hi);
        let expected = [wall, 290.0 + CP / FLUID_CP * (310.0 - wall)];
        for reversed in [false, true] {
            let thermal = pair(&wall_curve, &caloric, reversed, 0.0);
            let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
            let mut queries = Vec::new();
            let trial = step
                .evaluate(
                    &energy_increments(&step, expected, FLUID_CP),
                    |state| {
                        let source = evaluate_source_link_convection(
                            state,
                            &[circular(start, 2.0)],
                            1,
                            |_, _, temperature| {
                                queries.push(temperature.value());
                                liquid.at_temperature(temperature)
                            },
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
            assert_eq!(queries.len(), 1);
            near(queries[0], expected[1], 1e-12);
            let sample = &source.samples()[0];
            assert!(sample.surface.is_none());
            let evaluation = sample.evaluation.as_ref().unwrap();
            let section = evaluation.tube_section.as_ref().unwrap();
            assert_eq!(section.start.value(), start);
            assert_eq!(section.end.value(), 2.0);
            assert_eq!(section.diameter.value(), 0.1);
            assert_eq!(section.upstream.is_some(), start > 0.0);
            assert!(evaluation.rectangular_section.is_none());
            assert_eq!(evaluation.effective_nusselt(), section.nusselt);
            assert_eq!(
                evaluation.coefficient.model.discrepancy_rel.is_infinite(),
                start > 0.0
            );
            near(
                sample.conductance_w_k,
                LINK_AREA * circular_h(expected[1], start, 2.0),
                1e-10,
            );
            near(next.temperature_k(IDS[0]).unwrap(), expected[0], 1e-9);
            near(next.temperature_k(IDS[1]).unwrap(), expected[1], 1e-9);
            near(next.energy_j(), thermal.energy_j(), 1e-8);
            near(audit.external_heat_j, 0.0, 1e-12);
            near(audit.residual_j, 0.0, 1e-8);
            near(
                MASS * CP * (expected[0] - 310.0) + dt * sample.heat_into_fluid_w,
                0.0,
                1e-8,
            );
            near(thermal.temperature_k(IDS[0]).unwrap(), 310.0, 1e-12);
        }
        endpoints.push(wall);
    }
    assert!(
        endpoints[1] > endpoints[0] + 0.8,
        "a downstream laminar section must not repeat its inlet enhancement"
    );
}

const RECT_CP: f64 = 720.0;

// Authored, constant transport with EXACT source Pr=0.72, not a corrected query.
// The same Cp and density claims also construct the native coolant caloric law.
fn rectangular_source() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let original = liquid_card(0.0, 0.0, 320.0);
    let mut claims = ClaimSet::new();
    for (_, original) in original.claims().claims_ordered() {
        let mut claim = original.clone();
        let value = match claim.key.name() {
            SPECIFIC_HEAT_CAPACITY_PROPERTY => Some(RECT_CP),
            THERMAL_CONDUCTIVITY_PROPERTY => Some(1.0),
            DYNAMIC_VISCOSITY_PROPERTY => Some(0.001),
            _ => None,
        };
        if let Some(value) = value {
            claim.value = PropertyValue::Scalar {
                value,
                dims: claim.key.quantity().dims(),
            };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &context(280.0),
            upper: &context(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        -20.0 * RECT_CP,
        1e-8,
    )
    .unwrap();
    (card, caloric)
}

fn rectangular(area: f64, aspect_ratio: f64) -> SourceConvectiveLink {
    let diameter = 0.02;
    let start = 0.08;
    let perimeter = diameter * (1.0 + aspect_ratio) * ((1.0 + aspect_ratio) / aspect_ratio);
    SourceConvectiveLink {
        characteristic_length: Length::new(diameter),
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: SourceCellConvectionLaw::ForcedRectangularDuctSection {
            speed_m_s: 0.005,
            aspect_ratio,
            start: Length::new(start),
            end: Length::new(start + area / perimeter),
        },
        ..link()
    }
}

#[test]
fn four_wall_rectangular_section_uses_original_table_transfer_and_both_finite_stores() {
    let wall_curve = curve();
    let (card, caloric) = rectangular_source();
    let liquid = liquid_law(&card);
    // Both endpoints lie on Table 52's 10<=Gz<=20 chord. The derivative of
    // cumulative transfer is its independent intercept 4.20-(4.79-4.20)=3.61.
    let h = 3.61 / 0.02;
    let dt = 2.0;
    let conductance = LINK_AREA * h;
    let difference = 20.0 / (1.0 + dt * conductance * (1.0 / (MASS * CP) + 1.0 / (MASS * RECT_CP)));
    let heat = dt * conductance * difference;
    let expected = [310.0 - heat / (MASS * CP), 290.0 + heat / (MASS * RECT_CP)];
    for reversed in [false, true] {
        let thermal = pair(&wall_curve, &caloric, reversed, 0.0);
        let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
        let mut queries = Vec::new();
        let trial = step
            .evaluate(
                &energy_increments(&step, expected, RECT_CP),
                |state| {
                    let source = evaluate_source_link_convection(
                        state,
                        &[rectangular(LINK_AREA, 0.5)],
                        1,
                        |_, _, temperature| {
                            queries.push(temperature.value());
                            liquid.at_temperature(temperature)
                        },
                        || false,
                    )?;
                    Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: source.conductances_w_k().clone(),
                        payload: source,
                    })
                },
                || Ok(()),
            )
            .unwrap();
        let (next, source, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
        assert_eq!(queries.len(), 1);
        near(queries[0], expected[1], 1e-12);
        let sample = &source.samples()[0];
        assert!(sample.surface.is_none());
        let evaluation = sample.evaluation.as_ref().unwrap();
        assert!(evaluation.tube_section.is_none());
        let section = evaluation.rectangular_section.as_ref().unwrap();
        assert_eq!(section.hydraulic_diameter.value(), 0.02);
        assert_eq!(section.aspect_ratio, 0.5);
        assert_eq!(section.start.value(), 0.08);
        assert!(section.upstream.is_some());
        assert_eq!(evaluation.nusselt.groups()["Pr"], 0.72);
        near(evaluation.effective_nusselt(), 3.61, 1e-12);
        assert!(evaluation.coefficient.model.discrepancy_rel.is_infinite());
        near(sample.conductance_w_k, conductance, 1e-12);
        near(next.temperature_k(IDS[0]).unwrap(), expected[0], 1e-9);
        near(next.temperature_k(IDS[1]).unwrap(), expected[1], 1e-9);
        near(next.energy_j(), thermal.energy_j(), 1e-8);
        near(audit.external_heat_j, 0.0, 1e-12);
        near(audit.residual_j, 0.0, 1e-8);
        near(dt * sample.heat_into_fluid_w, heat, 1e-9);
    }
}

#[test]
fn laminar_sections_refuse_wrong_source_domains_partial_walls_and_film_substitution() {
    let wall_curve = curve();
    let (card, caloric) = rectangular_source();
    let liquid = liquid_law(&card);
    let thermal = pair(&wall_curve, &caloric, false, 0.0);
    for binding in [
        rectangular(2.0 * LINK_AREA, 0.5), // half of the four-wall perimeter
        rectangular(0.5 * LINK_AREA, 0.5), // native area exceeds actual channel
        rectangular(LINK_AREA, 0.4),       // original source supports exactly 0.5
        SourceConvectiveLink {
            law: SourceCellConvectionLaw::ForcedLaminarTubeSection {
                speed_m_s: 1.0,
                start: Length::new(1.0),
                end: Length::new(2.0),
            },
            ..circular(1.0, 2.0)
        }, // source is turbulent: never select another card silently
    ] {
        assert!(
            evaluate_source_link_convection(
                &thermal,
                &[binding],
                1,
                |_, _, temperature| liquid.at_temperature(temperature),
                || false
            )
            .is_err()
        );
    }
    for binding in [circular(1.0, 2.0), rectangular(LINK_AREA, 0.5)] {
        let binding = SourceConvectiveLink {
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            ..binding
        };
        let mut calls = 0;
        assert!(
            evaluate_source_link_convection(
                &thermal,
                &[binding],
                1,
                |_, _, temperature| {
                    calls += 1;
                    liquid.at_temperature(temperature)
                },
                || false
            )
            .is_err()
        );
        assert_eq!(calls, 0);
    }
    let (other_card, _) = fluid_source(1.0);
    let other_liquid = liquid_law(&other_card);
    assert!(
        evaluate_source_link_convection(
            &thermal,
            &[rectangular(LINK_AREA, 0.5)],
            1,
            |_, _, temperature| other_liquid.at_temperature(temperature),
            || false
        )
        .is_err()
    );
}
