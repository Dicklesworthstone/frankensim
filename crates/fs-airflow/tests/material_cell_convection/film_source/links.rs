//! One liquid source drives transport and caloric storage in a closed pair.
use super::*;

#[path = "links/observer.rs"]
mod observer;
#[path = "links/natural_geometry.rs"]
mod natural_geometry;
#[path = "links/mixed_convection.rs"]
mod mixed_convection;
use fs_airflow::graph::thermal::medium::convection::deposition::source::links::*;
use fs_matdb::SelectionPolicy;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};

const FLUID_CP: f64 = 4000.0;
const LINK_AREA: f64 = 0.002;

fn context(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
        .with("p", 100000.0)
        .unwrap()
        .with("solute_mass_fraction", 0.2)
        .unwrap()
}
fn fluid_source(viscosity_scale: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    fluid_source_with_expansion(viscosity_scale, 0.0)
}
fn fluid_source_with_expansion(
    viscosity_scale: f64,
    density_slope: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let original = liquid_card(density_slope, 0.01, 320.0);
    let mut claims = ClaimSet::new();
    for (_, value) in original.claims().claims_ordered() {
        let mut claim = value.clone();
        if claim.key.name() == SPECIFIC_HEAT_CAPACITY_PROPERTY {
            claim.value = PropertyValue::Scalar {
                value: FLUID_CP,
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        if claim.key.name() == DYNAMIC_VISCOSITY_PROPERTY {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                for (_, value) in knots {
                    *value *= viscosity_scale;
                }
            }
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
        -20.0 * FLUID_CP,
        1e-8,
    )
    .unwrap();
    (card, caloric)
}
fn pair<'a>(
    wall: &'a EquilibriumEnthalpyPhaseCurve,
    fluid: &'a ResolvedHeatCapacityCurve,
    reversed: bool,
    authored_g: f64,
) -> DepositionSimulation<'a> {
    // This fixture has constant Cp. Read its enthalpy from the actual source
    // chart so changing Cp data preserves the same initial 290 K fluid.
    let knots = fluid.curve().knots();
    let (lower, upper) = (&knots[0], &knots[knots.len() - 1]);
    let fluid_h = lower.specific_enthalpy_j_kg
        + (290.0 - lower.temperature_k) / (upper.temperature_k - lower.temperature_k)
            * (upper.specific_enthalpy_j_kg - lower.specific_enthalpy_j_kg);
    let mut events = IDS
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                body_id: i as u64,
                segment_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: MASS,
                specific_h_j_kg: if i == 0 { 10.0 * CP } else { fluid_h },
                geometry: CellGeometry {
                    start_m: [0.0, i as f64, 0.0],
                    end_m: [0.1, i as f64, 0.0],
                    volume_m3: 1e-6,
                    surface_area_m2: AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: None,
                },
            }),
        })
        .collect::<Vec<_>>();
    let [a, b] = if reversed { [IDS[1], IDS[0]] } else { IDS };
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 5,
            a,
            b,
            area_m2: LINK_AREA,
            conductance_w_k: authored_g,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), wall as &dyn ThermalLaw),
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 200,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn link() -> SourceConvectiveLink {
    SourceConvectiveLink {
        link: 5,
        fluid_cell: IDS[1],
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        law: SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::FlatPlateLaminarAverage,
            speed_m_s: 1.0,
            settings: ForcedConvectionSettings::default(),
        },
    }
}
// Independent flat-plate formula, using the authored liquid curves.
fn link_h(viscosity_scale: f64, temperature: f64) -> f64 {
    let k = 0.5 + 0.002 * (temperature - 300.0);
    let mu = viscosity_scale * 0.001 * (1.0 - 0.01 * (temperature - 300.0));
    let re = 1000.0 * 0.1 / mu;
    let pr = mu * FLUID_CP / k;
    0.664 * re.sqrt() * pr.powf(1.0 / 3.0) * k / 0.1
}
fn closed_endpoint(viscosity_scale: f64, dt: f64) -> [f64; 2] {
    let (mut lower, mut upper) = (294.0, 310.0);
    for _ in 0..80 {
        let wall = 0.5 * lower + 0.5 * upper;
        let fluid = 290.0 + CP / FLUID_CP * (310.0 - wall);
        let power = LINK_AREA * link_h(viscosity_scale, 0.5 * wall + 0.5 * fluid) * (wall - fluid);
        if MASS * CP * (wall - 310.0) + dt * power > 0.0 {
            upper = wall;
        } else {
            lower = wall;
        }
    }
    let wall = 0.5 * lower + 0.5 * upper;
    [wall, 290.0 + CP / FLUID_CP * (310.0 - wall)]
}

#[test]
fn g1_g3_liquid_source_changes_both_finite_stores_without_external_heat() {
    let wall_chart = curve();
    let run = |viscosity_scale, reversed| {
        let (card, caloric) = fluid_source(viscosity_scale);
        let fluid = liquid_law(&card);
        let thermal = pair(&wall_chart, &caloric, reversed, 0.0);
        let saved = thermal.cells().clone();
        let dt = 0.2;
        let expected = closed_endpoint(viscosity_scale, dt);
        let residual = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
        let increments = residual
            .cell_ids()
            .iter()
            .map(|&id| {
                let i = if id == IDS[0] { 0 } else { 1 };
                MASS * if i == 0 {
                    CP * (expected[0] - 310.0)
                } else {
                    FLUID_CP * (expected[1] - 290.0)
                }
            })
            .collect::<Vec<_>>();
        let trial = residual
            .evaluate(
                &increments,
                |state| {
                    let source = evaluate_source_link_convection(
                        state,
                        &[link()],
                        1,
                        |binding, state, t| {
                            assert_eq!(binding.fluid_cell, IDS[1]);
                            assert!(state.temperature_k(binding.fluid_cell).unwrap() > 290.0);
                            fluid.at_temperature(t)
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
        for (id, t) in IDS.into_iter().zip(expected) {
            near(next.temperature_k(id).unwrap(), t, 1e-9);
        }
        let sample = &source.samples()[0];
        assert_eq!(sample.fluid_cell, IDS[1]);
        assert_eq!(sample.wall_cell, IDS[0]);
        assert_eq!(
            sample.material.material().card_identity(),
            card.content_hash()
        );
        near(
            sample.property_temperature_k,
            0.5 * expected[0] + 0.5 * expected[1],
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
        assert_eq!(thermal.cells(), &saved);
        let mut mapped = BTreeMap::from([(99, 3.0)]);
        source.apply_to(&mut mapped, || false).unwrap();
        assert_eq!(mapped[&99], 3.0);
        assert_eq!(mapped[&5], sample.conductance_w_k);
        let saved_map = mapped.clone();
        assert!(source.apply_to(&mut mapped, || false).is_err());
        assert_eq!(mapped, saved_map);
        expected
    };
    let reference = run(1.0, false);
    assert_eq!(run(1.0, true), reference);
    let viscous = run(4.0, false);
    assert!(viscous[0] > reference[0]);
    assert!(viscous[1] < reference[1]);
}

#[test]
fn g0_original_link_ownership_and_source_refusals_preserve_both_cells() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let law = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    let saved = thermal.cells().clone();
    let mut calls = 0;
    assert!(
        evaluate_source_link_convection(
            &thermal,
            &[link(), link()],
            2,
            |_, _, t| {
                calls += 1;
                law.at_temperature(t)
            },
            || false
        )
        .is_err()
    );
    assert_eq!(calls, 0);
    let authored = pair(&chart, &caloric, false, 0.5);
    assert!(matches!(
        evaluate_source_link_convection(
            &authored,
            &[link()],
            1,
            |_, _, t| law.at_temperature(t),
            || false
        ),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Invalid("convective links require zero authored conductance")
        ))
    ));
    let wrong_endpoint = SourceConvectiveLink {
        fluid_cell: CellId(999),
        ..link()
    };
    assert!(
        evaluate_source_link_convection(
            &thermal,
            &[wrong_endpoint],
            1,
            |_, _, t| law.at_temperature(t),
            || false
        )
        .is_err()
    );
    let natural = SourceConvectiveLink {
        law: SourceCellConvectionLaw::VerticalPlate { gravity_m_s2: 9.81 },
        ..link()
    };
    assert!(matches!(
        evaluate_source_link_convection(
            &thermal,
            &[natural],
            1,
            |_, _, t| law.at_temperature(t),
            || false
        ),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Invalid(
                "natural convection requires sourced isobaric volumetric expansion"
            )
        ))
    ));
    assert!(matches!(
        evaluate_source_link_convection(
            &thermal,
            &[link()],
            1,
            |_, _, _| law.at_temperature(Temperature::new(300.1)),
            || false
        ),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Invalid(_)
        ))
    ));
    assert!(matches!(
        evaluate_source_link_convection(
            &thermal,
            &[link()],
            1,
            |_, _, _| Err::<fs_material::liquid::ResolvedLiquidState, _>("source failure"),
            || false
        ),
        Err(SourceCellConvectionError::Source {
            cell: CellId(10),
            source: "source failure",
            ..
        })
    ));
    assert!(matches!(
        evaluate_source_link_convection(
            &thermal,
            &[link()],
            1,
            |_, _, t| law.at_temperature(t),
            || true
        ),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Cancelled
        ))
    ));
    assert_eq!(thermal.cells(), &saved);
}
