//! Aiding mixed convection reaches actual source states and finite heat owners.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    AidingVerticalPlateSettings, evaluate_aiding_vertical_plate_from_material,
};

fn settings(speed: f64) -> AidingVerticalPlateSettings {
    AidingVerticalPlateSettings {
        speed_up_m_s: speed,
        gravity_m_s2: 9.81,
        maximum_relative_density_departure: 0.1,
    }
}

#[test]
fn aiding_direction_uses_signed_expansion_and_fan_stop_keeps_natural_heat() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let fluid = FluidTransportState::try_new(1.0, 1000.0, 0.02, 2e-5).unwrap();
    for beta in [0.002_f64, -0.002] {
        for (cell, wall) in [(IDS[0], 310.0), (IDS[1], 290.0)] {
            let speed = (beta * (wall - 300.0)).signum() * 0.5;
            let binding = |speed| MaterialConvectiveCell {
                cell,
                fluid: &fluid,
                reference_temperature: Temperature::new(300.0),
                characteristic_length: Length::new(0.1),
                law: CellConvectionLaw::AidingVerticalPlate {
                    settings: settings(speed),
                    expansion_per_k: beta,
                },
            };
            let run =
                |speed| evaluate_material_cell_convection(&thermal, &[binding(speed)], 1, || false);
            let mixed = run(speed).unwrap();
            let natural = run(0.0).unwrap();
            let almost_stopped = run(speed * 1e-10).unwrap();
            let sample = &mixed.samples()[0];
            assert_eq!(
                sample.evaluation.as_ref().unwrap().nusselt.card().id,
                CorrelationId::ManningQureshiAidingVerticalPlate
            );
            assert!(sample.heat_into_cell_w * (wall - 300.0) < 0.0);
            assert!(sample.heat_into_cell_w.abs() > natural.total_heat_into_cells_w().abs());
            assert!(natural.total_heat_into_cells_w() != 0.0);
            near(
                almost_stopped.total_heat_into_cells_w(),
                natural.total_heat_into_cells_w(),
                1e-10,
            );
            assert!(run(-speed).is_err());
        }
    }
    let equilibrium = initial(&chart, [300.0; 2], 0.0, 0.0);
    let binding = MaterialConvectiveCell {
        cell: IDS[0],
        fluid: &fluid,
        reference_temperature: Temperature::new(300.0),
        characteristic_length: Length::new(0.1),
        law: CellConvectionLaw::AidingVerticalPlate {
            settings: settings(0.0),
            expansion_per_k: 0.002,
        },
    };
    let resting = evaluate_material_cell_convection(&equilibrium, &[binding], 1, || false).unwrap();
    assert_eq!(resting.total_heat_into_cells_w(), 0.0);
    assert!(resting.samples()[0].evaluation.is_none());
}

#[test]
fn mixed_material_reads_one_snapshot_and_enforces_its_approximation_bound() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.0.set(self.0.get() + 1);
            FluidTransportState::try_new(1.0, 1000.0, 0.02, 2e-5)
        }
    }
    let fluid = Counted(std::cell::Cell::new(0));
    let evaluated = evaluate_aiding_vertical_plate_from_material(
        &fluid,
        settings(0.5),
        0.002,
        Temperature::new(310.0),
        Temperature::new(300.0),
        Length::new(0.1),
    )
    .unwrap();
    assert_eq!(fluid.0.get(), 1);
    near(evaluated.nusselt().groups()["Re"], 2500.0, 1e-10);
    near(evaluated.nusselt().groups()["Ra"], 490500.0, 1e-7);
    near(
        evaluated.coefficient().value.value(),
        evaluated.nusselt().evidence().value * 0.2,
        1e-12,
    );
    let strict = AidingVerticalPlateSettings {
        maximum_relative_density_departure: 0.01,
        ..settings(0.5)
    };
    assert!(
        evaluate_aiding_vertical_plate_from_material(
            &fluid,
            strict,
            0.002,
            Temperature::new(310.0),
            Temperature::new(300.0),
            Length::new(0.1),
        )
        .is_err()
    );
}

#[test]
fn aiding_source_requeries_actual_film_with_original_pressure_composition_and_heat_sign() {
    let chart = curve();
    let thermal = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    for slope in [0.5_f64, -0.5] {
        let card = liquid_card(slope, 0.01, 320.0);
        let law = liquid_law(&card)
            .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        let cells = IDS.map(|cell| SourceConvectiveCell {
            characteristic_length: Length::new(0.025),
            law: SourceCellConvectionLaw::AidingVerticalPlate {
                settings: settings(
                    (slope * (thermal.temperature_k(cell).unwrap() - 300.0)).signum() * 0.01,
                ),
            },
            ..film(cell)
        });
        let mut calls = 0;
        let result = evaluate_source_cell_convection(
            &thermal,
            &cells,
            2,
            |_, t| {
                calls += 1;
                law.at_temperature(t)
            },
            || false,
        )
        .unwrap();
        assert_eq!(calls, 2);
        for (source, heat) in result.materials().iter().zip(result.convection().samples()) {
            let film = 0.5 * heat.wall_temperature_k + 150.0;
            near(source.property_temperature_k, film, 1e-12);
            assert_eq!(source.material.material().query().axes()["p"], 100000.0);
            assert_eq!(
                source.material.material().query().axes()["solute_mass_fraction"],
                0.2
            );
            assert_eq!(
                source.material.material().card_identity(),
                card.content_hash()
            );
            let rho = 1000.0 - slope * (film - 300.0);
            let cp = 4000.0 + 2.0 * (film - 300.0);
            let k = 0.5 + 0.002 * (film - 300.0);
            let mu = 0.001 * (1.0 - 0.01 * (film - 300.0));
            let expected_ra =
                (9.81 * slope / rho * (heat.wall_temperature_k - 300.0) * 0.025_f64.powi(3)
                    / ((mu / rho) * (k / (rho * cp))))
                    .abs();
            let evaluated = heat.evaluation.as_ref().unwrap();
            near(evaluated.nusselt.groups()["Ra"], expected_ra, 1e-6);
            near(
                evaluated.nusselt.groups()["Re"],
                0.01 * 0.025 * rho / mu,
                1e-9,
            );
            near(
                heat.heat_into_cell_w,
                evaluated.nusselt.evidence().value * k / 0.025
                    * 0.008
                    * (300.0 - heat.wall_temperature_k),
                1e-9,
            );
        }
    }
}

// Independent average Eq. (10), with the fixture's actual source curves.
fn link_h(wall: f64, fluid: f64, speed: f64) -> f64 {
    let film = 0.5 * wall + 0.5 * fluid;
    let rho = 1000.0 - 0.5 * (film - 300.0);
    let k = 0.5 + 0.002 * (film - 300.0);
    let mu = 0.001 * (1.0 - 0.01 * (film - 300.0));
    let re = speed * 0.025 * rho / mu;
    let pr = mu * FLUID_CP / k;
    let gr = 9.81 * 0.5 / rho * (wall - fluid) * 0.025_f64.powi(3) / (mu / rho).powi(2);
    let forced = 2.0 * 0.339 * pr.powf(1.0 / 3.0) / (1.0 + 0.100 * pr.powf(-0.75)).powf(2.0 / 9.0)
        * re.sqrt();
    let natural = 4.0 / 3.0 * 0.503 * pr.powf(0.25)
        / (1.0 + 0.670 * pr.powf(-9.0 / 16.0)).powf(4.0 / 9.0)
        * gr.powf(0.25);
    let n = 3.5 * pr.powf(0.075);
    (forced.powf(n) + natural.powf(n)).powf(1.0 / n) * k / 0.025
}

#[test]
fn aiding_finite_link_closes_native_heat_balance_with_running_or_stopped_fan() {
    let chart = curve();
    let (card, caloric) = fluid_source_with_expansion(1.0, 0.5);
    let law = liquid_law(&card)
        .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
        .unwrap();
    let thermal = pair(&chart, &caloric, false, 0.0);
    let dt = 0.02;
    for speed in [0.0, 0.01, 0.1] {
        let binding = SourceConvectiveLink {
            characteristic_length: Length::new(0.025),
            law: SourceCellConvectionLaw::AidingVerticalPlate {
                settings: settings(speed),
            },
            ..link()
        };
        let (mut lower, mut upper) = (294.0, 310.0);
        for _ in 0..80 {
            let wall = 0.5 * lower + 0.5 * upper;
            let fluid = 290.0 + CP / FLUID_CP * (310.0 - wall);
            let flux = LINK_AREA * link_h(wall, fluid, speed) * (wall - fluid);
            if MASS * CP * (wall - 310.0) + dt * flux > 0.0 {
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
        near(next.energy_j(), thermal.energy_j(), 1e-8);
        near(audit.external_heat_j, 0.0, 1e-12);
        assert!(sample.heat_into_fluid_w > 0.0);
    }
}

#[test]
fn mixed_source_bulk_property_shortcut_refuses_before_source_queries() {
    let chart = curve();
    let (card, caloric) = fluid_source_with_expansion(1.0, 0.5);
    let law = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    let mixed = SourceCellConvectionLaw::AidingVerticalPlate {
        settings: settings(0.01),
    };
    let mut calls = 0;
    let cell = SourceConvectiveCell {
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: mixed,
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
        law: mixed,
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
