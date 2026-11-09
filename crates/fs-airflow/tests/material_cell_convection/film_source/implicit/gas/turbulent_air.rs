//! G1/G3: explicit air wall temperature correction on the original native heat path.
use super::*;
use fs_airflow::graph::thermal::medium::convection::implicit::PrescribedMassFlowConvection;
use fs_convection::developing::TurbulentTubePropertyCorrection;

fn configuration(start: f64) -> ExteriorConvectionSettings {
    let mut config = settings(SourceCellConvectionLaw::ForcedTubeSection {
        speed_m_s: 5.0,
        start: Length::new(start),
        end: Length::new(2.0),
    });
    config.property_temperature = ConvectionPropertyTemperature::Reservoir;
    config
}

// Independent Gnielinski algebra; only the original material supplies bulk
// properties. No production correlation, corrected section or tangent is used.
fn reference_coefficient(
    gas: &IsobaricGasLaw,
    bulk: f64,
    wall: f64,
    start: f64,
    speed: f64,
) -> f64 {
    let original = gas.at_temperature(Temperature::new(bulk)).unwrap();
    let state = original.state();
    let re = state.density * speed * 0.1 / state.dynamic_viscosity;
    let pr = state.dynamic_viscosity * state.specific_heat_cp / state.thermal_conductivity;
    let friction = (1.8 * re.log10() - 1.5).powi(-2);
    let developed = friction / 8.0 * (re - 1000.0) * pr
        / (1.0 + 12.7 * (friction / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0));
    let a = (start / 0.1).cbrt();
    let b = 20.0_f64.cbrt();
    developed * (1.0 + 1.0 / (a * a + a * b + b * b)) * state.thermal_conductivity / 0.1
        * (bulk / wall).powf(0.45)
}

#[test]
fn g1_air_correction_uses_actual_kelvin_temperatures_and_retains_partial_axes() {
    let body = curve();
    let air = gas("300 100000 0");
    let original = air.reference().clone();
    for start in [0.0, 1.0] {
        let baseline =
            ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &air, configuration(start))
                .unwrap();
        let uncorrected = baseline.evaluate(310.0, 300.0).unwrap();
        assert!(uncorrected.turbulent_property_correction.is_none());
        assert_eq!(uncorrected.coefficient_surface_derivative_w_m2_k2, 0.0);
        let corrected = baseline.with_gnielinski_air_wall_correction().unwrap();
        assert!(core::ptr::eq(
            corrected.caloric_law(),
            &body as &dyn ThermalLaw
        ));
        for bulk in [290.0, 300.0, 310.0] {
            for wall in [290.0, 300.0, 310.0] {
                let actual = corrected.evaluate(wall, bulk).unwrap();
                let expected = reference_coefficient(&air, bulk, wall, start, 5.0);
                near(actual.coefficient.value.value(), expected, 2e-12 * expected);
                near(
                    actual.coefficient_surface_derivative_w_m2_k2,
                    -0.45 * expected / wall,
                    2e-13 * expected,
                );
                assert_eq!(actual.property_temperature_k, bulk);
                assert_eq!(
                    actual.transport,
                    air.at_temperature(Temperature::new(bulk))
                        .unwrap()
                        .transport_state()
                        .unwrap()
                );
                assert!(actual.surface_transport.is_none());
                assert_eq!(
                    actual.turbulent_property_correction,
                    Some(TurbulentTubePropertyCorrection::Air {
                        bulk_temperature: Temperature::new(bulk),
                        wall_temperature: Temperature::new(wall),
                    })
                );
                let partials = &actual.coefficient.sensitivity.d_qoi;
                near(partials["T_bulk"], 0.45 * expected / bulk, 2e-13 * expected);
                near(
                    partials["T_wall"],
                    -0.45 * expected / wall,
                    2e-13 * expected,
                );
                assert!(!partials.contains_key("Pr_wall"));
                let section = actual.tube_section.as_ref().unwrap();
                near(
                    section.nusselt * actual.transport.thermal_conductivity_w_m_k() / 0.1,
                    expected,
                    2e-12 * expected,
                );
                assert_eq!(section.start.value(), start);
                assert_eq!(section.upstream.is_some(), start > 0.0);
                derivative(&corrected, wall, bulk);
            }
        }
    }
    assert_eq!(air.reference(), &original);
}

#[test]
fn g3_air_mass_flow_keeps_actual_bulk_density_correction_and_independent_tangents() {
    let body = curve();
    let air = gas("300 100000 0");
    let source =
        ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &air, configuration(1.0))
            .unwrap()
            .with_gnielinski_air_wall_correction()
            .unwrap();
    let rate = 0.06;
    let area = 0.01;
    let flow = PrescribedMassFlowConvection::for_gas(&source, rate, area).unwrap();
    for bulk in [290.0, 310.0] {
        for wall in [290.0, 310.0] {
            let result = flow.evaluate(wall, bulk).unwrap();
            let rho = air
                .at_temperature(Temperature::new(bulk))
                .unwrap()
                .state()
                .density;
            let speed = rate / (rho * area);
            let expected = reference_coefficient(&air, bulk, wall, 1.0, speed);
            near(result.speed_m_s, speed, 1e-14);
            near(
                result.convection.coefficient.value.value(),
                expected,
                2e-12 * expected,
            );
            near(
                result.convection.coefficient_surface_derivative_w_m2_k2,
                -0.45 * expected / wall,
                2e-13 * expected,
            );
            let delta = 1e-6;
            let at = |rate| {
                PrescribedMassFlowConvection::for_gas(&source, rate, area)
                    .unwrap()
                    .evaluate(wall, bulk)
                    .unwrap()
                    .convection
                    .coefficient
                    .value
                    .value()
            };
            let numerical = (at(rate + delta) - at(rate - delta)) / (2.0 * delta);
            near(
                result.coefficient_mass_flow_derivative_w_m2_k_s_kg.unwrap(),
                numerical,
                2e-7 * numerical.abs(),
            );
            let flux = flow.flux_state_at_surface(0.0, wall, bulk).unwrap();
            near(
                flux.convection_temperature_tangent_w_m2_k,
                expected * (1.0 - 0.45 * (wall - bulk) / wall),
                2e-11 * expected,
            );
        }
    }
}

#[test]
fn g0_air_exponent_requires_explicit_dry_air_geometry_and_original_source_support() {
    let body = curve();
    for text in [
        "300 100000 0.2",
        "mixture 300 100000 mole 1 CO2 1",
        "mixture 300 100000 mole 2 N2 0.79 O2 0.21",
    ] {
        let medium = gas(text);
        let original = medium.reference().clone();
        let source =
            ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &medium, configuration(1.0))
                .unwrap();
        assert!(source.evaluate(310.0, 300.0).is_ok());
        assert!(source.with_gnielinski_air_wall_correction().is_err());
        assert_eq!(medium.reference(), &original);
    }
    let air = gas("300 100000 0");
    let wrong_law =
        ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &air, settings(plate()))
            .unwrap();
    assert!(wrong_law.with_gnielinski_air_wall_correction().is_err());
    let mut unsupported_wall = configuration(1.0);
    unsupported_wall.surface_temperature_bounds_k = [281.0, 330.0];
    assert!(
        ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &air, unsupported_wall,)
            .unwrap()
            .with_gnielinski_air_wall_correction()
            .is_err()
    );
    let source =
        ImplicitGasConvection::try_new(&body, [-20000.0, 20000.0], &air, configuration(1.0))
            .unwrap()
            .with_gnielinski_air_wall_correction()
            .unwrap();
    for (wall, bulk) in [
        (330.0, 300.0),
        (300.0, 330.0),
        (f64::NAN, 300.0),
        (300.0, 0.0),
    ] {
        assert!(source.evaluate(wall, bulk).is_err());
    }
    assert!(source.state_at_surface(30000.0, 310.0, 300.0).is_err());
}

fn surface_reference(air: &IsobaricGasLaw, body: f64) -> (f64, [f64; 2]) {
    let powers = |wall: f64| {
        [
            AREA * reference_coefficient(air, 300.0, wall, 1.0, 5.0) * (wall - 300.0),
            AREA * 0.7 * 5.670_374_419e-8 * (wall.powi(4) - 305.0_f64.powi(4)),
        ]
    };
    let (mut lo, mut hi) = (300.0, body);
    for _ in 0..90 {
        let wall = 0.5 * (lo + hi);
        let conducted = AREA / 0.005
            * (0.5 * (body - wall) + 0.001 * ((body - 300.0).powi(2) - (wall - 300.0).powi(2)));
        let heat = powers(wall);
        if conducted > heat[0] + heat[1] {
            lo = wall;
        } else {
            hi = wall;
        }
    }
    let wall = 0.5 * (lo + hi);
    (wall, powers(wall))
}

#[test]
fn g3_air_correction_reaches_native_conductive_surface_storage_and_bath_entropy() {
    let air = gas("300 100000 0");
    let original = air.reference().clone();
    let final_temperature =
        super::native::advance_with(&air, configuration(1.0), true, surface_reference);
    assert!(final_temperature > 300.0 && final_temperature < 315.0);
    assert_eq!(air.reference(), &original);
}
