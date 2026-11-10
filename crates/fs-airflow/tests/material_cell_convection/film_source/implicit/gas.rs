//! Original gas transport and analytic film tangents on native solid surfaces.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    evaluate_flat_plate_section_from_material, evaluate_sphere_from_material,
    evaluate_tube_from_material,
    implicit::{ExteriorConvectionSettings, ImplicitGasConvection},
};
use fs_material::gas_species::input::{read_gas_record, temperature::IsobaricGasLaw};

#[path = "../../../finite_gas_convection/source.rs"]
mod thermochemical_source;

fn gas(text: &str) -> IsobaricGasLaw {
    IsobaricGasLaw::try_from_record(
        read_gas_record(&mut text.split_ascii_whitespace(), None, None).unwrap(),
    )
    .unwrap()
}

fn plate() -> SourceCellConvectionLaw {
    SourceCellConvectionLaw::ForcedPlateSection {
        correlation: CorrelationId::FlatPlateLaminarAverage,
        speed_m_s: 1.0,
        start: Length::new(0.025),
    }
}

fn derivative(source: &ImplicitGasConvection<'_, '_>, surface: f64, bath: f64) {
    let result = source.evaluate(surface, bath).unwrap();
    let delta = 1e-3;
    let h = |t| source.evaluate(t, bath).unwrap().coefficient.value.value();
    let numerical = (h(surface + delta) - h(surface - delta)) / (2.0 * delta);
    near(
        result.coefficient_surface_derivative_w_m2_k2,
        numerical,
        2e-7 * numerical.abs().max(1.0),
    );
    let flux = source.flux_state_at_surface(0.0, surface, bath).unwrap();
    near(
        flux.convection_temperature_tangent_w_m2_k,
        result.coefficient.value.value()
            + (surface - bath) * result.coefficient_surface_derivative_w_m2_k2,
        1e-10,
    );
}

#[test]
fn film_sections_keep_original_gas_pressure_composition_and_analytic_transport_response() {
    let body = curve();
    let mut coefficients = Vec::new();
    for text in [
        "300 100000 0",
        "300 100000 0.2",
        "mixture 300 100000 mole 2 He 0.1 N2 0.9",
    ] {
        let gas = gas(text);
        let reference = gas.reference().clone();
        let source =
            ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &gas, settings(plate()))
                .unwrap();
        for surface in [285.0, 315.0] {
            let film = 0.5 * (surface + 300.0);
            let original = gas.at_temperature(Temperature::new(film)).unwrap();
            let expected = evaluate_flat_plate_section_from_material(
                CorrelationId::FlatPlateLaminarAverage,
                &original,
                1.0,
                Length::new(0.025),
                Length::new(0.1),
            )
            .unwrap();
            let result = source.evaluate(surface, 300.0).unwrap();
            near(
                result.coefficient.value.value(),
                expected.coefficient().value.value(),
                1e-10,
            );
            assert_eq!(result.property_temperature_k, film);
            assert_eq!(result.transport, original.transport_state().unwrap());
            assert_eq!(
                result.nusselt.groups(),
                expected.section().downstream().groups()
            );
            assert_eq!(result.plate_section.as_ref().unwrap().start.value(), 0.025);
            assert!(result.coefficient.model.discrepancy_rel.is_infinite());
            assert_eq!(original.state().pressure, reference.state().pressure);
            near(
                original.state().water_mole_fraction,
                reference.state().water_mole_fraction,
                1e-16,
            );
            derivative(&source, surface, 300.0);
        }
        coefficients.push(
            source
                .evaluate(315.0, 300.0)
                .unwrap()
                .coefficient
                .value
                .value(),
        );
        assert_eq!(source.gas().reference(), &reference);
    }
    assert!(
        (coefficients[2] - coefficients[0]).abs() > 0.5,
        "an original helium mixture must change actual surface cooling"
    );
    let low = gas("300 50000 0");
    let high = gas("300 100000 0");
    let h = |gas| {
        ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], gas, settings(plate()))
            .unwrap()
            .evaluate(315.0, 300.0)
            .unwrap()
            .coefficient
            .value
            .value()
    };
    near(h(&high) / h(&low), 2.0_f64.sqrt(), 1e-12);
}

#[test]
fn natural_and_aiding_gas_tangents_include_the_original_inverse_temperature_expansion() {
    let body = curve();
    let gas = gas("mixture 300 100000 mole 1 CO2 1");
    for surface in [290.0, 310.0] {
        for law in [
            SourceCellConvectionLaw::Natural {
                geometry: NaturalConvectionGeometry::Sphere,
                gravity_m_s2: 9.81,
            },
            SourceCellConvectionLaw::AidingVerticalPlate {
                settings: mixed_settings(if surface > 300.0 { 0.2 } else { -0.2 }, 9.81),
            },
        ] {
            let mut config = settings(law);
            if let SourceCellConvectionLaw::AidingVerticalPlate { ref mut settings } = config.law {
                settings.maximum_relative_density_departure = 0.1;
            }
            let source =
                ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &gas, config).unwrap();
            let film = 0.5 * (surface + 300.0);
            let original = gas.at_temperature(Temperature::new(film)).unwrap();
            let beta = original
                .thermal_transport_state()
                .unwrap()
                .volumetric_expansion_per_k()
                .unwrap();
            near(beta, 1.0 / film, 1e-16);
            let expected = match config.law {
                SourceCellConvectionLaw::Natural {
                    geometry,
                    gravity_m_s2,
                } => evaluate_natural_from_material(
                    geometry,
                    &original,
                    gravity_m_s2,
                    beta,
                    Temperature::new(surface),
                    Temperature::new(300.0),
                    Length::new(0.1),
                )
                .unwrap(),
                SourceCellConvectionLaw::AidingVerticalPlate { settings } => {
                    evaluate_aiding_vertical_plate_from_material(
                        &original,
                        settings,
                        beta,
                        Temperature::new(surface),
                        Temperature::new(300.0),
                        Length::new(0.1),
                    )
                    .unwrap()
                }
                _ => unreachable!(),
            };
            near(
                source
                    .evaluate(surface, 300.0)
                    .unwrap()
                    .coefficient
                    .value
                    .value(),
                expected.coefficient().value.value(),
                1e-10,
            );
            derivative(&source, surface, 300.0);
        }
    }
    let sphere = ImplicitGasConvection::try_new(
        &body,
        [-20000.0, 20000.0],
        &gas,
        settings(SourceCellConvectionLaw::Natural {
            geometry: NaturalConvectionGeometry::Sphere,
            gravity_m_s2: 9.81,
        }),
    )
    .unwrap();
    let flux = sphere.flux_state_at_surface(0.0, 300.0, 300.0).unwrap();
    near(
        flux.convection_temperature_tangent_w_m2_k,
        2.0 * gas.reference().state().thermal_conductivity / 0.1,
        1e-12,
    );
    assert!(
        sphere.evaluate(300.0, 300.0).is_err(),
        "finite equilibrium heat tangent must not invent a regular coefficient derivative"
    );
}

#[test]
fn bulk_wall_viscosity_corrections_query_the_same_gas_without_a_film_substitution() {
    let body = curve();
    let gas = gas("mixture 300 100000 mole 1 CO2 1");
    for law in [
        SourceCellConvectionLaw::ForcedSphere { speed_m_s: 1.0 },
        SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::SiederTateTurbulent,
            speed_m_s: 5.0,
            settings: ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(100.0),
                ..Default::default()
            },
        },
    ] {
        let mut config = settings(law);
        config.property_temperature = ConvectionPropertyTemperature::Reservoir;
        let source =
            ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &gas, config).unwrap();
        // Colder CO2 wall gives Whitaker's supported mu_bulk/mu_wall >= 1.
        let wall = gas.at_temperature(Temperature::new(290.0)).unwrap();
        let expected = match law {
            SourceCellConvectionLaw::ForcedSphere { speed_m_s } => {
                evaluate_sphere_from_material(gas.reference(), &wall, speed_m_s, Length::new(0.1))
                    .unwrap()
                    .bulk()
                    .coefficient()
                    .value
                    .value()
            }
            SourceCellConvectionLaw::Forced { speed_m_s, .. } => evaluate_tube_from_material(
                gas.reference(),
                &wall,
                speed_m_s,
                Length::new(0.1),
                100.0,
            )
            .unwrap()
            .bulk()
            .coefficient()
            .value
            .value(),
            _ => unreachable!(),
        };
        let actual = source.evaluate(290.0, 300.0).unwrap();
        assert_eq!(actual.property_temperature_k, 300.0);
        assert_eq!(
            actual.surface_dynamic_viscosity_pa_s,
            Some(wall.state().dynamic_viscosity)
        );
        near(actual.coefficient.value.value(), expected, 1e-10);
        derivative(&source, 290.0, 300.0);
    }
}

#[test]
fn gas_phase_property_interval_and_pointwise_geometry_fail_without_mutation() {
    let body = curve();
    let humid = gas("300 100000 0.9");
    let saved = humid.reference().clone();
    assert!(
        ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &humid, settings(plate()),)
            .is_err(),
        "the declared cold film interval crosses actual saturation"
    );
    assert_eq!(humid.reference(), &saved);
    let dry = gas("300 100000 0");
    let source =
        ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &dry, settings(plate()))
            .unwrap();
    assert!(source.evaluate(315.0, 340.0).is_err()); // film exceeds declared 320 K
    let mut invalid = settings(SourceCellConvectionLaw::ForcedRectangularDuctSection {
        speed_m_s: 1.0,
        aspect_ratio: 0.5,
        start: Length::new(1.0),
        end: Length::new(2.0),
    });
    invalid.property_temperature = ConvectionPropertyTemperature::Reservoir;
    assert!(ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &dry, invalid).is_err());
}

#[path = "gas/native.rs"]
mod native;
#[path = "gas/surface.rs"]
mod surface;
#[path = "gas/turbulent_air.rs"]
mod turbulent_air;
#[path = "gas/laminar_tube.rs"]
mod laminar_tube;
