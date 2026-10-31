//! Oriented natural convection reaches the original finite thermal stores.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    HorizontalPlateFacing, NaturalConvectionGeometry,
};

#[test]
fn horizontal_orientation_and_negative_expansion_select_actual_cell_heat() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let fluid = FluidTransportState::try_new(1.0, 1000.0, 0.02, 2e-5).unwrap();
    for facing in [HorizontalPlateFacing::Up, HorizontalPlateFacing::Down] {
        for beta in [0.002, -0.002] {
            let cells = IDS.map(|cell| MaterialConvectiveCell {
                cell,
                fluid: &fluid,
                reference_temperature: Temperature::new(300.0),
                characteristic_length: Length::new(0.1),
                law: CellConvectionLaw::Natural {
                    geometry: NaturalConvectionGeometry::HorizontalPlate(facing),
                    gravity_m_s2: 9.81,
                    expansion_per_k: beta,
                },
            });
            let result = evaluate_material_cell_convection(&thermal, &cells, 2, || false).unwrap();
            for sample in result.samples() {
                let difference = sample.wall_temperature_k - 300.0;
                let unstable = (beta * difference > 0.0) == (facing == HorizontalPlateFacing::Up);
                let ra = (9.81 * beta * difference * 0.1_f64.powi(3) / (2e-5 * 2e-5)).abs();
                let nu = if unstable { 0.54 } else { 0.27 } * ra.powf(0.25);
                let evaluation = sample.evaluation.as_ref().unwrap();
                assert_eq!(
                    evaluation.nusselt.card().id,
                    if unstable {
                        CorrelationId::HorizontalPlateUnstableLaminar
                    } else {
                        CorrelationId::HorizontalPlateStable
                    }
                );
                near(
                    sample.heat_into_cell_w,
                    -0.008 * nu * 0.02 / 0.1 * difference,
                    1e-11,
                );
                assert!(sample.heat_into_cell_w * difference < 0.0);
            }
        }
    }
    let equilibrium = initial(&chart, [300.0; 2], 0.0, 0.0);
    for geometry in [
        NaturalConvectionGeometry::HorizontalCylinder,
        NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Up),
    ] {
        let cell = MaterialConvectiveCell {
            law: CellConvectionLaw::Natural {
                geometry,
                gravity_m_s2: 9.81,
                expansion_per_k: 0.002,
            },
            ..natural(&fluid, IDS[0])
        };
        let zero = evaluate_material_cell_convection(&equilibrium, &[cell], 1, || false).unwrap();
        assert_eq!(zero.total_heat_into_cells_w(), 0.0);
        assert!(zero.samples()[0].evaluation.is_none());
    }
}

#[test]
fn horizontal_source_uses_original_film_state_full_contrast_and_authored_planform_length() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    // A 0.1 m square has L=A/P=0.01/0.4=0.025 m. Contact masking changes
    // exposed heat-transfer area to 0.008 m², never the full plate perimeter.
    let length = 0.01 / 0.4;
    for density_slope in [0.5, -0.5] {
        let card = liquid_card(density_slope, 0.01, 320.0);
        let law = liquid_law(&card)
            .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        let cells = IDS.map(|cell| SourceConvectiveCell {
            characteristic_length: Length::new(length),
            law: SourceCellConvectionLaw::Natural {
                geometry: NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Up),
                gravity_m_s2: 9.81,
            },
            ..film(cell)
        });
        let result = evaluate_source_cell_convection(
            &thermal,
            &cells,
            2,
            |_, t| law.at_temperature(t),
            || false,
        )
        .unwrap();
        for (source, heat) in result.materials().iter().zip(result.convection().samples()) {
            let t = 0.5 * heat.wall_temperature_k + 150.0;
            near(source.property_temperature_k, t, 1e-12);
            assert_eq!(source.material.material().query().axes()["p"], 100000.0);
            assert_eq!(
                source.material.material().query().axes()["solute_mass_fraction"],
                0.2
            );
            assert_eq!(
                source.material.material().card_identity(),
                card.content_hash()
            );
            let rho = 1000.0 - density_slope * (t - 300.0);
            let cp = 4000.0 + 2.0 * (t - 300.0);
            let k = 0.5 + 0.002 * (t - 300.0);
            let mu = 0.001 * (1.0 - 0.01 * (t - 300.0));
            let difference = heat.wall_temperature_k - 300.0;
            let ra = (9.81 * density_slope / rho * difference * length.powi(3)
                / ((mu / rho) * (k / (rho * cp))))
                .abs();
            let nu = if density_slope * difference > 0.0 {
                0.54
            } else {
                0.27
            } * ra.powf(0.25);
            let evaluation = heat.evaluation.as_ref().unwrap();
            near(evaluation.nusselt.groups()["Ra"], ra, 1e-7);
            near(heat.exposed_area_m2, 0.008, 1e-15);
            near(
                heat.heat_into_cell_w,
                -0.008 * nu * k / length * difference,
                1e-8,
            );
        }
    }
}

// Independent ASHRAE horizontal-plate algebra, with this fixture's source
// curves and constant Cp shared with its original finite fluid caloric law.
fn horizontal_link_h(
    density_slope: f64,
    wall: f64,
    fluid: f64,
    facing: HorizontalPlateFacing,
) -> f64 {
    let film = 0.5 * wall + 0.5 * fluid;
    let rho = 1000.0 - density_slope * (film - 300.0);
    let k = 0.5 + 0.002 * (film - 300.0);
    let mu = 0.001 * (1.0 - 0.01 * (film - 300.0));
    let ra = (9.81 * density_slope / rho * (wall - fluid) * 0.025_f64.powi(3)
        / ((mu / rho) * (k / (rho * FLUID_CP))))
        .abs();
    let unstable = (density_slope * (wall - fluid) > 0.0) == (facing == HorizontalPlateFacing::Up);
    let nu = if unstable {
        if ra <= 8.0e6 {
            0.54 * ra.powf(0.25)
        } else {
            0.15 * ra.cbrt()
        }
    } else {
        0.27 * ra.powf(0.25)
    };
    nu * k / 0.025
}

#[test]
fn horizontal_finite_link_exchanges_equal_heat_at_both_actual_endpoints() {
    let chart = curve();
    let dt = 0.02;
    for (density_slope, facing) in [
        (0.5, HorizontalPlateFacing::Up),
        (0.5, HorizontalPlateFacing::Down),
        (-0.5, HorizontalPlateFacing::Up),
    ] {
        let (card, caloric) = fluid_source_with_expansion(1.0, density_slope);
        let law = liquid_law(&card)
            .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        let thermal = pair(&chart, &caloric, false, 0.0);
        let binding = SourceConvectiveLink {
            characteristic_length: Length::new(0.025),
            law: SourceCellConvectionLaw::Natural {
                geometry: NaturalConvectionGeometry::HorizontalPlate(facing),
                gravity_m_s2: 9.81,
            },
            ..link()
        };
        let (mut lower, mut upper) = (294.0, 310.0);
        for _ in 0..80 {
            let wall = 0.5 * lower + 0.5 * upper;
            let fluid = 290.0 + CP / FLUID_CP * (310.0 - wall);
            let heat =
                LINK_AREA * horizontal_link_h(density_slope, wall, fluid, facing) * (wall - fluid);
            if MASS * CP * (wall - 310.0) + dt * heat > 0.0 {
                upper = wall;
            } else {
                lower = wall;
            }
        }
        let wall = 0.5 * lower + 0.5 * upper;
        let fluid = 290.0 + CP / FLUID_CP * (310.0 - wall);
        let residual = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
        let increments: Vec<_> = residual
            .cell_ids()
            .iter()
            .map(|&cell| {
                if cell == IDS[0] {
                    MASS * CP * (wall - 310.0)
                } else {
                    MASS * FLUID_CP * (fluid - 290.0)
                }
            })
            .collect();
        let trial = residual
            .evaluate(
                &increments,
                |state| {
                    let sources = evaluate_source_link_convection(
                        state,
                        &[binding],
                        1,
                        |_, _, t| law.at_temperature(t),
                        || false,
                    )?;
                    Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: sources.conductances_w_k().clone(),
                        payload: sources,
                    })
                },
                || Ok(()),
            )
            .unwrap();
        let (next, sources, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
        let sample = &sources.samples()[0];
        near(
            sample.heat_into_fluid_w * dt,
            MASS * CP * (310.0 - wall),
            1e-8,
        );
        near(
            sample.heat_into_fluid_w * dt,
            MASS * FLUID_CP * (fluid - 290.0),
            1e-8,
        );
        near(
            sample.property_temperature_k,
            0.5 * wall + 0.5 * fluid,
            1e-12,
        );
        near(next.energy_j(), thermal.energy_j(), 1e-8);
        near(audit.external_heat_j, 0.0, 1e-12);
        assert!(sample.heat_into_fluid_w > 0.0);
    }
}

#[test]
fn horizontal_source_laws_refuse_bulk_properties_before_resolving_any_material() {
    let chart = curve();
    let (card, caloric) = fluid_source_with_expansion(1.0, 0.5);
    let law = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    for geometry in [
        NaturalConvectionGeometry::HorizontalCylinder,
        NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Down),
    ] {
        let natural = SourceCellConvectionLaw::Natural {
            geometry,
            gravity_m_s2: 9.81,
        };
        let mut calls = 0;
        let cell = SourceConvectiveCell {
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            law: natural,
            ..film(IDS[0])
        };
        assert!(
            evaluate_source_cell_convection(
                &thermal,
                &[cell],
                1,
                |_, t| {
                    calls += 1;
                    law.at_temperature(t)
                },
                || false
            )
            .is_err()
        );
        let binding = SourceConvectiveLink {
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            law: natural,
            ..link()
        };
        assert!(
            evaluate_source_link_convection(
                &thermal,
                &[binding],
                1,
                |_, _, t| {
                    calls += 1;
                    law.at_temperature(t)
                },
                || false
            )
            .is_err()
        );
        assert_eq!(calls, 0);
    }
}
