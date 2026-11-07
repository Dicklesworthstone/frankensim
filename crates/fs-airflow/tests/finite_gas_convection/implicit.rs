//! Preserve original thermochemical sources when wall-property derivatives are zero.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    deposition::source::{ConvectionPropertyTemperature, SourceCellConvectionLaw},
    evaluate_laminar_tube_section_from_material,
    implicit::{ExteriorConvectionError, ExteriorConvectionSettings, ImplicitGasConvection},
};
use fs_material::gas_species::input::temperature::IsobaricGasLaw;
use fs_qty::Temperature;

#[test]
fn bulk_held_thermochemical_gas_uses_original_sources_and_refuses_missing_active_tangents() {
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
    for configuration in [
        ExteriorConvectionSettings {
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: SourceCellConvectionLaw::ForcedPlateSection {
                correlation: CorrelationId::FlatPlateLaminarAverage,
                speed_m_s: 1.0,
                start: Length::new(0.0),
            },
            ..settings
        },
        ExteriorConvectionSettings {
            law: SourceCellConvectionLaw::ForcedSphere { speed_m_s: 1.0 },
            ..settings
        },
    ] {
        assert!(matches!(
            ImplicitGasConvection::try_new(&body, [0.0, 360000.0], &gas, configuration,),
            Err(ExteriorConvectionError::Material(_))
        ));
    }
    assert_eq!(convection.gas().reference(), &saved);
    assert_eq!(gas.reference(), &saved);
}
