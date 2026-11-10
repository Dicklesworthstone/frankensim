//! Original gas wall viscosity in the whole heated laminar-tube mean.
use super::*;
use fs_airflow::graph::thermal::medium::convection::evaluate_laminar_tube_from_material;

#[test]
fn g1_laminar_gas_mean_retains_original_bulk_wall_and_analytic_tangent() {
    let body = curve();
    for text in ["300 100000 0", "mixture 300 100000 mole 1 CO2 1"] {
        let gas = gas(text);
        let saved = gas.reference().clone();
        let bulk = gas.reference().transport_state().unwrap();
        let diameter = Length::new(0.1);
        let speed = 1000.0 * bulk.dynamic_viscosity_pa_s()
            / (bulk.density_kg_m3() * diameter.value());
        let mut config = settings(SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::SiederTateLaminar,
            speed_m_s: speed,
            settings: ForcedConvectionSettings {
                length_over_hydraulic_diameter: Some(20.0),
                ..Default::default()
            },
        });
        config.property_temperature = ConvectionPropertyTemperature::Reservoir;
        let source = ImplicitGasConvection::try_new(
            &body, [-20000.0, 20000.0], &gas, config,
        ).unwrap();
        let prepared = gas.prepare_transport_temperature_tangent().unwrap();
        for surface in [290.0, 300.0, 310.0] {
            let wall = gas.at_temperature(Temperature::new(surface)).unwrap();
            let expected = evaluate_laminar_tube_from_material(
                gas.reference(), &wall, speed, diameter, 20.0,
            ).unwrap();
            let actual = source.evaluate(surface, 300.0).unwrap();
            let h = actual.coefficient.value.value();
            near(h, expected.bulk().coefficient().value.value(), 2e-12 * h);
            assert_eq!(actual.property_temperature_k, 300.0);
            assert_eq!(actual.transport, bulk);
            assert_eq!(actual.nusselt.card().id, CorrelationId::SiederTateLaminar);
            assert_eq!(actual.nusselt.groups()["L_over_Dh"], 20.0);
            assert!(actual.tube_section.is_none());
            assert!(actual.turbulent_property_correction.is_none());
            assert_eq!(actual.surface_dynamic_viscosity_pa_s, Some(wall.state().dynamic_viscosity));
            let tangent = prepared.evaluate(Temperature::new(surface)).unwrap();
            near(actual.coefficient_surface_derivative_w_m2_k2,
                -0.14 * h * tangent.dynamic_viscosity_derivative_pa_s_k
                    / wall.state().dynamic_viscosity, 1e-12 * h);
            derivative(&source, surface, 300.0);
        }
        assert_eq!(gas.reference(), &saved);
        config.property_temperature = ConvectionPropertyTemperature::FilmMean;
        assert!(ImplicitGasConvection::try_new(
            &body, [-20000.0, 20000.0], &gas, config,
        ).is_err());
    }
}
