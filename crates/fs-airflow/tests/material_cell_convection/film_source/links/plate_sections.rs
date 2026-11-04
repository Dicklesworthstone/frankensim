//! Section-dependent film transfer in the original closed fluid/wall residual.

use super::*;
use super::super::plate_sections::independent_h;

fn section_link() -> SourceConvectiveLink {
    SourceConvectiveLink {
        characteristic_length: Length::new(0.2),
        law: SourceCellConvectionLaw::ForcedPlateSection {
            correlation: CorrelationId::FlatPlateLaminarAverage,
            speed_m_s: 0.1,
            start: Length::new(0.1),
        },
        ..link()
    }
}

#[test]
fn original_finite_fluid_and_wall_close_energy_with_source_refreshed_plate_sections() {
    let wall_chart = curve();
    let mut transfers = Vec::new();
    for viscosity_scale in [1.0, 2.0] {
        let (card, caloric) = fluid_source(viscosity_scale);
        let fluid = liquid_law(&card);
        for reversed in [false, true] {
            let thermal = pair(&wall_chart, &caloric, reversed, 0.0);
            let saved = thermal.cells().clone();
            let dt = 0.2;
            let capacities = [MASS * CP, MASS * FLUID_CP];
            let mut lower = 0.0;
            let mut upper = 20.0 / (1.0 / capacities[0] + 1.0 / capacities[1]);
            for _ in 0..80 {
                let heat = 0.5 * lower + 0.5 * upper;
                let wall = 310.0 - heat / capacities[0];
                let bulk = 290.0 + heat / capacities[1];
                let material = fluid.at_temperature(Temperature::new(0.5 * wall + 0.5 * bulk)).unwrap();
                let h = independent_h(0.1, 0.2, 0.1, material.transport_state().unwrap(),
                    CorrelationId::FlatPlateLaminarAverage);
                if heat - dt * LINK_AREA * h * (wall - bulk) > 0.0 {
                    upper = heat;
                } else {
                    lower = heat;
                }
            }
            let heat = 0.5 * lower + 0.5 * upper;
            let expected = [310.0 - heat / capacities[0], 290.0 + heat / capacities[1]];
            let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
            let increments = step.cell_ids().iter().map(|&id| {
                if id == IDS[0] { -heat } else { heat }
            }).collect::<Vec<_>>();
            let trial = step.evaluate(&increments, |state| {
                let source = evaluate_source_link_convection(
                    state, &[section_link()], 1,
                    |_, _, t| fluid.at_temperature(t), || false,
                )?;
                let mut conductances_w_k = BTreeMap::new();
                source.apply_to(&mut conductances_w_k, || false)?;
                Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k,
                    payload: source,
                })
            }, || Ok(())).unwrap();
            let (next, source, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
            for (cell, temperature) in IDS.into_iter().zip(expected) {
                near(next.temperature_k(cell).unwrap(), temperature, 1e-9);
            }
            let sample = &source.samples()[0];
            let evaluation = sample.evaluation.as_ref().unwrap();
            let details = evaluation.section.as_ref().unwrap();
            assert_eq!((sample.wall_cell, sample.fluid_cell), (IDS[0], IDS[1]));
            assert_eq!(details.start.value(), 0.1);
            assert_eq!(details.end.value(), 0.2);
            assert!(details.upstream.is_some());
            assert!(evaluation.coefficient.model.discrepancy_rel.is_infinite());
            assert_eq!(sample.material.material().card_identity(), card.content_hash());
            near(sample.property_temperature_k, 0.5 * expected[0] + 0.5 * expected[1], 1e-12);
            assert!((sample.property_temperature_k - 300.0).abs() > 0.01);
            near(sample.heat_into_fluid_w * dt, heat, 1e-8);
            near(next.energy_j(), thermal.energy_j(), 1e-8);
            near(audit.external_heat_j, 0.0, 1e-12);
            near(audit.residual_j, 0.0, 1e-8);
            assert_eq!(next.links()[&5].conductance_w_k, 0.0);
            assert_eq!(thermal.cells(), &saved);
            transfers.push(heat);
        }
    }
    near(transfers[0], transfers[1], 1e-12);
    near(transfers[2], transfers[3], 1e-12);
    assert!(transfers[2] < transfers[0]);
}

#[test]
fn invalid_native_link_strip_refuses_before_a_source_query() {
    let wall_chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let fluid = liquid_law(&card);
    let thermal = pair(&wall_chart, &caloric, false, 0.0);
    let saved = thermal.cells().clone();
    let binding = SourceConvectiveLink {
        characteristic_length: Length::new(0.1),
        ..section_link()
    };
    let mut queries = 0;
    assert!(evaluate_source_link_convection(
        &thermal, &[binding], 1,
        |_, _, t| { queries += 1; fluid.at_temperature(t) }, || false,
    ).is_err());
    assert_eq!(queries, 0);
    assert_eq!(thermal.cells(), &saved);
}
