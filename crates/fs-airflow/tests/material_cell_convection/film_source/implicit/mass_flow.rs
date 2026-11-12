//! Original bulk-density velocity -> film/tube convection -> native surface heat.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    evaluate_flat_plate_section_from_material,
    implicit::{ImplicitGasConvection, PrescribedMassFlowConvection},
};
use fs_material::gas_species::input::{read_gas_record, temperature::IsobaricGasLaw};

fn gas(text: &str) -> IsobaricGasLaw {
    IsobaricGasLaw::try_from_record(
        read_gas_record(&mut text.split_ascii_whitespace(), None, None).unwrap(),
    )
    .unwrap()
}

fn plate() -> SourceCellConvectionLaw {
    SourceCellConvectionLaw::ForcedPlateSection {
        correlation: CorrelationId::FlatPlateLaminarAverage,
        speed_m_s: 9.0, // Replaced by the actual bulk-density mass-flow velocity.
        start: Length::new(0.025),
    }
}

fn surface_derivative(flow: &PrescribedMassFlowConvection<'_>, surface: f64, bath: f64) {
    let result = flow.evaluate(surface, bath).unwrap().convection;
    let delta = 1e-3;
    let h = |temperature| {
        flow.evaluate(temperature, bath)
            .unwrap()
            .convection
            .coefficient
            .value
            .value()
    };
    let numerical = (h(surface + delta) - h(surface - delta)) / (2.0 * delta);
    near(
        result.coefficient_surface_derivative_w_m2_k2,
        numerical,
        2e-7 * numerical.abs().max(1.0),
    );
    let flux = flow.flux_state_at_surface(0.0, surface, bath).unwrap();
    near(
        flux.convection_temperature_tangent_w_m2_k,
        result.coefficient.value.value()
            + (surface - bath) * result.coefficient_surface_derivative_w_m2_k2,
        1e-9,
    );
}

#[test]
fn gas_mass_flow_uses_actual_bulk_density_and_retains_distinct_film_properties() {
    let body = curve();
    let mut pressure_coefficients = Vec::new();
    for text in ["300 100000 0", "300 50000 0"] {
        let gas = gas(text);
        let before = gas.reference().clone();
        let source =
            ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &gas, settings(plate()))
                .unwrap();
        let mass = 0.01;
        let area = 0.01;
        let flow = PrescribedMassFlowConvection::for_gas(&source, mass, area).unwrap();
        assert!(core::ptr::eq(flow.caloric_law(), source.caloric_law()));
        for bath in [290.0, 300.0, 310.0] {
            let surface = 315.0;
            let bulk = gas.at_temperature(Temperature::new(bath)).unwrap();
            let film = gas
                .at_temperature(Temperature::new(0.5 * (surface + bath)))
                .unwrap();
            let velocity = mass / area / bulk.state().density;
            let expected = evaluate_flat_plate_section_from_material(
                CorrelationId::FlatPlateLaminarAverage,
                &film,
                velocity,
                Length::new(0.025),
                Length::new(0.1),
            )
            .unwrap();
            let actual = flow.evaluate(surface, bath).unwrap();
            near(actual.speed_m_s, velocity, 1e-15);
            assert_eq!(actual.bulk_transport, bulk.transport_state().unwrap());
            assert_eq!(actual.convection.transport, film.transport_state().unwrap());
            assert_ne!(
                actual.bulk_transport.density_kg_m3(),
                actual.convection.transport.density_kg_m3()
            );
            near(
                actual.convection.coefficient.value.value(),
                expected.coefficient().value.value(),
                1e-10,
            );
            // The actual laminar strip retains exponent 1/2 in mass flow.
            near(
                actual.coefficient_mass_flow_derivative_w_m2_k_s_kg.unwrap(),
                0.5 * actual.convection.coefficient.value.value() / mass,
                1e-8,
            );
            surface_derivative(&flow, surface, bath);
        }
        pressure_coefficients.push(
            flow.evaluate(315.0, 300.0)
                .unwrap()
                .convection
                .coefficient
                .value
                .value(),
        );
        assert_eq!(gas.reference(), &before);
    }
    // For this ideal gas source, rho_film/rho_bulk cancels pressure at fixed
    // mass flux. A held speed or film-density velocity would give other physics.
    near(pressure_coefficients[0], pressure_coefficients[1], 1e-12);
}

#[test]
fn original_liquid_bulk_density_and_all_forced_section_rate_partials_are_retained() {
    let body = curve();
    let card = liquid_card(0.3, 0.01, 320.0);
    let liquid = liquid_law(&card);
    let area = 0.01;
    for (law, convention, mass) in [
        (plate(), ConvectionPropertyTemperature::FilmMean, 1.0),
        (
            SourceCellConvectionLaw::Forced {
                correlation: CorrelationId::DittusBoelter,
                speed_m_s: 0.1,
                settings: ForcedConvectionSettings {
                    length_over_hydraulic_diameter: Some(50.0),
                    ..Default::default()
                },
            },
            ConvectionPropertyTemperature::Reservoir,
            2.0,
        ),
        (
            SourceCellConvectionLaw::ForcedSphere { speed_m_s: 0.1 },
            ConvectionPropertyTemperature::Reservoir,
            1.0,
        ),
        (
            SourceCellConvectionLaw::ForcedTubeSection {
                speed_m_s: 0.1,
                start: Length::new(1.0),
                end: Length::new(2.0),
            },
            ConvectionPropertyTemperature::Reservoir,
            // GnielinskiDeveloping admits Re 1e4 to 1e6: (m/A) D / mu is 18182 to 22222 over
            // the three baths here; 1.0 kg/s gave 9091 at 290 K.
            2.0,
        ),
        (
            SourceCellConvectionLaw::ForcedLaminarTubeSection {
                speed_m_s: 0.01,
                start: Length::new(0.01),
                end: Length::new(1.0),
            },
            ConvectionPropertyTemperature::Reservoir,
            0.1,
        ),
    ] {
        let mut config = settings(law);
        config.property_temperature = convention;
        let source =
            ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, config).unwrap();
        let flow = PrescribedMassFlowConvection::for_liquid(&source, mass, area).unwrap();
        for bath in [290.0, 300.0, 310.0] {
            let surface = 315.0;
            let actual = flow.evaluate(surface, bath).unwrap();
            let bulk = liquid
                .at_temperature(Temperature::new(bath))
                .unwrap()
                .transport_state()
                .unwrap();
            for (prepared, original) in [
                (actual.bulk_transport.density_kg_m3(), bulk.density_kg_m3()),
                (
                    actual.bulk_transport.specific_heat_j_kg_k(),
                    bulk.specific_heat_j_kg_k(),
                ),
                (
                    actual.bulk_transport.thermal_conductivity_w_m_k(),
                    bulk.thermal_conductivity_w_m_k(),
                ),
                (
                    actual.bulk_transport.dynamic_viscosity_pa_s(),
                    bulk.dynamic_viscosity_pa_s(),
                ),
            ] {
                near(prepared, original, 4e-14 * original.abs());
            }
            near(actual.speed_m_s, mass / area / bulk.density_kg_m3(), 1e-15);
            let delta = mass * 1e-5;
            let value = |rate| {
                PrescribedMassFlowConvection::for_liquid(&source, rate, area)
                    .unwrap()
                    .evaluate(surface, bath)
                    .unwrap()
                    .convection
                    .coefficient
                    .value
                    .value()
            };
            let numerical = (value(mass + delta) - value(mass - delta)) / (2.0 * delta);
            near(
                actual.coefficient_mass_flow_derivative_w_m2_k_s_kg.unwrap(),
                numerical,
                3e-8 * numerical.abs().max(1.0),
            );
            surface_derivative(&flow, surface, bath);
            if convention == ConvectionPropertyTemperature::Reservoir
                && !matches!(law, SourceCellConvectionLaw::ForcedSphere { .. })
            {
                assert_eq!(
                    actual.convection.coefficient_surface_derivative_w_m2_k2,
                    0.0
                );
            }
        }
    }
}

#[test]
fn prescribed_flow_refuses_unowned_geometry_sources_and_stagnant_regimes_without_a_floor() {
    let body = curve();
    let gas = gas("300 100000 0");
    let source =
        ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &gas, settings(plate()))
            .unwrap();
    for (mass, area) in [
        (-0.1, 0.01),
        (f64::NAN, 0.01),
        (1.0, 0.0),
        (1.0, f64::INFINITY),
    ] {
        assert!(PrescribedMassFlowConvection::for_gas(&source, mass, area).is_err());
    }
    let stopped = PrescribedMassFlowConvection::for_gas(&source, 0.0, 0.01).unwrap();
    assert!(stopped.evaluate(315.0, 300.0).is_err());
    let flowing = PrescribedMassFlowConvection::for_gas(&source, 0.01, 0.01).unwrap();
    // Film 307.5 K would be supported, but the ACTUAL bulk 330 K is not.
    assert!(flowing.evaluate(285.0, 330.0).is_err());
    for law in [
        SourceCellConvectionLaw::Natural {
            geometry: NaturalConvectionGeometry::Sphere,
            gravity_m_s2: 9.81,
        },
        SourceCellConvectionLaw::AidingVerticalPlate {
            settings: mixed_settings(0.1, 9.81),
        },
    ] {
        let natural =
            ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &gas, settings(law))
                .unwrap();
        assert!(PrescribedMassFlowConvection::for_gas(&natural, 0.01, 0.01).is_err());
    }
    // The original Re upper boundary is admitted at fixed rho and mu. Its
    // missing TWO-SIDED flow derivative must not become zero or block the
    // independently valid native Pr/k surface tangent.
    let card = liquid_card(0.0, 0.0, 320.0);
    let liquid = liquid_law(&card);
    let source =
        ImplicitLiquidConvection::try_new(&body, [-20000.0, 20000.0], &liquid, settings(plate()))
            .unwrap();
    let bounded = PrescribedMassFlowConvection::for_liquid(&source, 50.0, 0.01).unwrap();
    let response = bounded.evaluate(310.0, 300.0).unwrap();
    assert_eq!(response.convection.nusselt.groups()["Re"], 500000.0);
    assert!(
        response
            .coefficient_mass_flow_derivative_w_m2_k_s_kg
            .is_none()
    );
    assert!(bounded.flux_state_at_surface(0.0, 310.0, 300.0).is_ok());
}

#[test]
fn prescribed_mass_flow_drives_the_original_conducting_surface_and_heat_tangent() {
    let card = isobaric_liquid_card(0.0, 0.01, 320.0);
    let material = isobaric_liquid_law(&card);
    let point = |temperature| {
        material
            .reference_point()
            .clone()
            .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100000.0)
            .unwrap()
            .with("T", temperature)
            .unwrap()
    };
    let cp_key = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let k_key = PropertyKey::with_quantity(
        THERMAL_CONDUCTIVITY_PROPERTY,
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    );
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let caloric = IsobaricHeatCapacityLaw::try_new(
        &card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &cp_key,
            temperature_axis: "T",
            pressure_axis: "p",
            lower: &point(280.0),
            upper: &point(320.0),
            selection,
        },
        300.0,
    )
    .unwrap();
    let conduction = IsobaricConductionLaw::try_new(&card, &caloric, &k_key, selection).unwrap();
    let layer = CaloricConductionLayer::from_isobaric(&conduction).unwrap();
    let gas = gas("300 100000 0");
    let source = ImplicitGasConvection::try_new(
        layer.caloric_law(),
        layer.energy_bounds_j_kg(),
        &gas,
        settings(plate()),
    )
    .unwrap();
    let mass = 0.01;
    let area = 0.01;
    let flow = PrescribedMassFlowConvection::for_gas(&source, mass, area).unwrap();
    let boundary = IntegratedCaloricGrayBoundary::new(
        &layer,
        ImplicitBoundaryGeometry {
            area_m2: AREA,
            half_distance_m: 0.005,
            surface_resistance_m2_k_w: 0.0,
        },
        0.0,
        0.0,
        305.0,
        LayerConductionConfig {
            max_iterations: 64,
            absolute_heat_tolerance_w: 1e-10,
            relative_heat_tolerance: 1e-12,
        },
    )
    .unwrap()
    .with_exterior_law(&flow)
    .unwrap();
    let energy = |temperature: f64| 4000.0 * (temperature - 300.0) + (temperature - 300.0).powi(2);
    for bath in [290.0, 300.0, 310.0] {
        let body = 315.0;
        let response = boundary
            .evaluate(energy(body), bath, &mut || Ok(()))
            .unwrap();
        let surface = response.surface_temperature_k;
        let bulk = gas.at_temperature(Temperature::new(bath)).unwrap();
        let film = gas
            .at_temperature(Temperature::new(0.5 * (surface + bath)))
            .unwrap();
        let film = film.state();
        let speed = mass / area / bulk.state().density;
        let reynolds = film.density * speed * 0.1 / film.dynamic_viscosity;
        let prandtl = film.dynamic_viscosity * film.specific_heat_cp / film.thermal_conductivity;
        // Independent original laminar plate equation: x1/x2=1/4 strip.
        let h = 0.664 * reynolds.sqrt() * prandtl.cbrt() * film.thermal_conductivity / 0.1 / 1.5;
        let convection = AREA * h * (surface - bath);
        let radiation = AREA * 0.7 * 5.670374419e-8 * (surface.powi(4) - 305.0_f64.powi(4));
        let conducted = AREA / 0.005
            * (0.5 * (body - surface)
                + 0.001 * ((body - 300.0).powi(2) - (surface - 300.0).powi(2)));
        near(response.convection_power_w, convection, 1e-8);
        near(response.radiation_power_w, radiation, 1e-8);
        near(conducted, convection + radiation, 1e-8);
        near(response.surface_residual_w, 0.0, 1e-8);
        assert!(surface < body && surface > bath);
        let delta = 1e-3;
        let total = |temperature| {
            let state = boundary
                .evaluate(energy(temperature), bath, &mut || Ok(()))
                .unwrap();
            state.convection_power_w + state.radiation_power_w
        };
        near(
            response.boundary.power_enthalpy_derivative_w_kg_j,
            (total(body + delta) - total(body - delta))
                / (energy(body + delta) - energy(body - delta)),
            1e-9,
        );
    }
    assert!(
        boundary
            .evaluate(energy(315.0), 300.0, &mut || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
}
