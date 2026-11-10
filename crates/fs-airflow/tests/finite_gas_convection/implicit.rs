//! Original thermochemical calorics and transport on implicit gas surfaces.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    deposition::source::{ConvectionPropertyTemperature, SourceCellConvectionLaw},
    evaluate_flat_plate_section_from_material, evaluate_laminar_tube_section_from_material,
    evaluate_sphere_from_material,
    implicit::{ExteriorConvectionSettings, ImplicitGasConvection},
};
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas_species::input::temperature::IsobaricGasLaw;
use fs_qty::Temperature;
use fs_thermal::enthalpy::conduction::implicit::boundary::ImplicitGrayExteriorLaw;

#[test]
fn bulk_held_thermochemical_gas_keeps_original_sources_and_zero_wall_property_tangents() {
    let gas = IsobaricGasLaw::try_from_model(
        source::model_transport(1.0, true),
        source::point(400.0, 100000.0),
    )
    .unwrap();
    let saved = gas.reference().clone();
    let body = source::body_curve();
    let settings = ExteriorConvectionSettings {
        fluid_temperature_bounds_k: [300.0, 600.0],
        surface_temperature_bounds_k: [300.0, 600.0],
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: SourceCellConvectionLaw::ForcedLaminarTubeSection {
            speed_m_s: 0.1,
            start: Length::new(1.0),
            end: Length::new(2.0),
        },
        emissivity: 0.0,
    };
    let convection =
        ImplicitGasConvection::try_new(&body, [0.0, 360000.0], &gas, settings).unwrap();
    let mut coefficients = Vec::new();
    for bath in [350.0, 450.0] {
        let original = gas.at_temperature(Temperature::new(bath)).unwrap();
        let expected = evaluate_laminar_tube_section_from_material(
            &original,
            0.1,
            Length::new(0.1),
            Length::new(1.0),
            Length::new(2.0),
        )
        .unwrap();
        for surface in [350.0, 450.0, 550.0] {
            let result = convection.evaluate(surface, bath).unwrap();
            near(
                result.coefficient.value.value(),
                expected.coefficient().value.value(),
                1e-12,
            );
            assert_eq!(result.property_temperature_k, bath);
            assert_eq!(result.coefficient_surface_derivative_w_m2_k2, 0.0);
        }
        coefficients.push(expected.coefficient().value.value());
    }
    assert!(
        (coefficients[0] - coefficients[1]).abs() > 0.1,
        "bulk-held means independent of wall temperature, not frozen at construction"
    );
    assert_eq!(convection.gas().reference(), &saved);
    assert_eq!(gas.reference(), &saved);
}

fn film_settings() -> ExteriorConvectionSettings {
    ExteriorConvectionSettings {
        fluid_temperature_bounds_k: [300.0, 600.0],
        surface_temperature_bounds_k: [300.0, 600.0],
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        law: SourceCellConvectionLaw::ForcedPlateSection {
            correlation: CorrelationId::FlatPlateLaminarAverage,
            speed_m_s: 1.0,
            start: Length::new(0.025),
        },
        emissivity: 0.0,
    }
}

#[test]
fn active_thermochemical_film_retains_nasa_heat_capacity_and_original_mixture_tangents() {
    let body = source::body_curve();
    let cp_r_slope = 0.002;
    for fraction in [0.8, 1.0] {
        let gas = IsobaricGasLaw::try_from_model(
            source::model_transport_with_heat_capacity_slope(fraction, cp_r_slope),
            source::point(400.0, 100000.0),
        ).unwrap();
        let saved = gas.reference().clone();
        let convection = ImplicitGasConvection::try_new(
            &body, [0.0, 360000.0], &gas, film_settings(),
        ).unwrap();
        for surface in [350.0, 500.0, 550.0] {
            let film = 0.5 * (surface + 400.0);
            let original = gas.at_temperature(Temperature::new(film)).unwrap();
            let expected = evaluate_flat_plate_section_from_material(
                CorrelationId::FlatPlateLaminarAverage, &original, 1.0,
                Length::new(0.025), Length::new(0.1),
            ).unwrap();
            let actual = convection.evaluate(surface, 400.0).unwrap();
            near(actual.coefficient.value.value(), expected.coefficient().value.value(), 1e-11);
            assert_eq!(actual.transport, original.transport_state().unwrap());
            assert_eq!(actual.property_temperature_k, film);
            // In these original source packs mu and k are proportional to T.
            // Their Wilke ratios are constant, rho is proportional to 1/T,
            // and therefore flat-plate h is proportional to Cp^(1/3).
            // This independent homogeneity identity detects dropping NASA Cp'.
            let molar_mass = fraction * 0.028 + (1.0 - fraction) * 0.004;
            let cp_slope = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                * cp_r_slope / molar_mass;
            let expected_slope = actual.coefficient.value.value() * cp_slope
                / (6.0 * original.state().specific_heat_cp);
            near(actual.coefficient_surface_derivative_w_m2_k2, expected_slope, 2e-12);
            assert!(expected_slope > 0.0);
            let flux = convection.flux_state_at_surface(0.0, surface, 400.0).unwrap();
            let delta = 0.01;
            let power = |t| convection.flux_state_at_surface(0.0, t, 400.0)
                .unwrap().coefficients.convection_w_m2_k * (t - 400.0);
            near(flux.convection_temperature_tangent_w_m2_k,
                (power(surface + delta) - power(surface - delta)) / (2.0 * delta), 1e-7);
        }
        assert!(convection.evaluate(650.0, 400.0).is_err());
        assert_eq!(gas.reference(), &saved);
    }
}

#[test]
fn active_thermochemical_wall_viscosity_uses_its_original_bulk_and_surface_sources() {
    let body = source::body_curve();
    let gas = IsobaricGasLaw::try_from_model(
        source::model_transport_with_heat_capacity_slope(1.0, 0.002),
        source::point(400.0, 100000.0),
    ).unwrap();
    let settings = ExteriorConvectionSettings {
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: SourceCellConvectionLaw::ForcedSphere { speed_m_s: 1.0 },
        ..film_settings()
    };
    let convection = ImplicitGasConvection::try_new(
        &body, [0.0, 360000.0], &gas, settings,
    ).unwrap();
    let bulk = gas.at_temperature(Temperature::new(400.0)).unwrap();
    for surface in [320.0, 360.0, 390.0] {
        let wall = gas.at_temperature(Temperature::new(surface)).unwrap();
        let expected = evaluate_sphere_from_material(
            &bulk, &wall, 1.0, Length::new(0.1),
        ).unwrap();
        let actual = convection.evaluate(surface, 400.0).unwrap();
        near(actual.coefficient.value.value(), expected.bulk().coefficient().value.value(), 1e-11);
        assert_eq!(actual.property_temperature_k, 400.0);
        let delta = 0.01;
        let value = |t| convection.evaluate(t, 400.0).unwrap().coefficient.value.value();
        near(actual.coefficient_surface_derivative_w_m2_k2,
            (value(surface + delta) - value(surface - delta)) / (2.0 * delta), 1e-9);
        assert!(actual.coefficient_surface_derivative_w_m2_k2 < 0.0);
    }
}
