use super::*;

fn liquid(t: f64) -> WaterConductivityState {
    liquid_at_pressure_temperature(Pressure::new(300_000.0), Temperature::new(t)).unwrap()
}

#[test]
fn original_iapws_verification_points_and_unretuned_ice_regression() {
    // SR6-08(2011), Table 8. The printed 260 K value is outside this API's
    // deliberately non-extrapolating liquid interval and is tested as a refusal.
    for (t, published) in [(298.15, 0.606_502_308), (375.0, 0.677_913_788)] {
        let state = liquid(t);
        assert!((state.conductivity().value() - published).abs() < 5.1e-10);
        assert_eq!(
            state.source(),
            WaterConductivitySource::LiquidIapwsSr6082011
        );
        assert_eq!(state.temperature().value(), t);
        assert_eq!(state, liquid(t));
    }
    // At T*=1 the original four coefficients sum without fractional powers.
    assert!((liquid(300.0).conductivity().value() - 0.609_485).abs() < 8e-16);

    // Independent standard-library evaluation of Yen equation (33); the source
    // regression is not adjusted to the report's rounded 2.2 W/(m K) discussion.
    for t in [243.15, 250.0, 273.15, 273.16] {
        let value = ice_at_temperature(Temperature::new(t)).unwrap();
        let original = 9.828 * (-0.0057_f64 * t).exp();
        assert!((value.conductivity().value() - original).abs() < 2e-14);
        assert_eq!(value.source(), WaterConductivitySource::IceYen1981Empirical);
        assert_eq!(value, ice_at_temperature(Temperature::new(t)).unwrap());
    }
    let near_melting = ice_at_temperature(Temperature::new(273.15))
        .unwrap()
        .conductivity()
        .value();
    assert!((2.06..2.08).contains(&near_melting));
}

#[test]
fn analytic_transport_derivatives_match_temperature_changes() {
    for t in [274.0, 298.15, 330.0, 375.0, 382.0] {
        let lo = t - 0.001;
        let hi = t + 0.001;
        let difference =
            (liquid(hi).conductivity().value() - liquid(lo).conductivity().value()) / (hi - lo);
        assert!((difference - liquid(t).temperature_derivative().value()).abs() < 5e-12);
    }
    for t in [244.0, 253.15, 270.0] {
        let lo = t - 0.001;
        let hi = t + 0.001;
        let at = |t| ice_at_temperature(Temperature::new(t)).unwrap();
        let difference =
            (at(hi).conductivity().value() - at(lo).conductivity().value()) / (hi - lo);
        assert!((difference - at(t).temperature_derivative().value()).abs() < 5e-12);
        assert!(at(t).temperature_derivative().value() < 0.0);
        assert_eq!(
            at(t).temperature_derivative().value(),
            -0.0057 * at(t).conductivity().value()
        );
    }
}

#[test]
fn exact_domain_edges_pressure_applicability_and_nonfinite_refusals() {
    for t in [LIQUID_MINIMUM_TEMPERATURE_K, LIQUID_MAXIMUM_TEMPERATURE_K] {
        let temperature = Temperature::new(t);
        let saturation = if97::saturation_pressure(temperature).unwrap();
        assert_eq!(
            liquid_at_pressure_temperature(saturation, temperature).unwrap(),
            liquid_at_pressure_temperature(Pressure::new(300_000.0), temperature).unwrap()
        );
        assert!(matches!(
            liquid_at_pressure_temperature(
                Pressure::new(saturation.value().next_down()),
                temperature
            ),
            Err(WaterConductivityError::LiquidPressureOutOfRange { .. })
        ));
    }
    for t in [
        260.0,
        LIQUID_MINIMUM_TEMPERATURE_K.next_down(),
        LIQUID_MAXIMUM_TEMPERATURE_K.next_up(),
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        assert!(matches!(
            liquid_at_pressure_temperature(Pressure::new(101_325.0), Temperature::new(t)),
            Err(WaterConductivityError::TemperatureOutOfRange { .. })
        ));
    }
    for t in [ICE_MINIMUM_TEMPERATURE_K, ICE_MAXIMUM_TEMPERATURE_K] {
        assert!(ice_at_temperature(Temperature::new(t)).is_ok());
    }
    for t in [
        ICE_MINIMUM_TEMPERATURE_K.next_down(),
        ICE_MAXIMUM_TEMPERATURE_K.next_up(),
        f64::NAN,
        f64::INFINITY,
    ] {
        assert!(matches!(
            ice_at_temperature(Temperature::new(t)),
            Err(WaterConductivityError::TemperatureOutOfRange { .. })
        ));
    }
    for p in [0.0, -1.0, 300_000.0_f64.next_up(), f64::NAN, f64::INFINITY] {
        assert!(matches!(
            liquid_at_pressure_temperature(Pressure::new(p), Temperature::new(298.15)),
            Err(WaterConductivityError::LiquidPressureOutOfRange { .. })
        ));
    }
    let normal = liquid_at_pressure_temperature(Pressure::new(101_325.0), Temperature::new(298.15));
    assert_eq!(normal.unwrap(), liquid(298.15));
}
