use super::*;

fn support() -> WetPoreLiquidSupport {
    WetPoreLiquidSupport {
        liquid_saturation: (0.0, 1.0),
        temperature_k: (280.0, 350.0),
        liquid_pressure_pa: (10000.0, 200000.0),
        maximum_superficial_speed_m_s: 0.25,
    }
}
fn law(slr: f64, sgr: f64, ell: f64, scale: f64) -> WetPoreLiquidLaw {
    WetPoreLiquidLaw::try_new(
        slr,
        sgr,
        0.5,
        ell,
        scale,
        support(),
        "synthetic calibrated wetting branch",
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 8e-12 * a.abs().max(b.abs()).max(1e-300),
        "{a} != {b}"
    );
}
fn pressure(state: &WetPoreLiquidState) -> f64 {
    match state.capillary_pressure() {
        LiquidCapillaryPressure::Finite(value) => value,
        LiquidCapillaryPressure::DryLimit => panic!("expected a finite capillary pressure"),
    }
}

#[test]
fn normalized_m_half_matches_independent_wetting_and_capillary_oracle() {
    let calibration = law(0.1, 0.1, 0.5, 800.0);
    let state = calibration.sample(0.58).unwrap(); // Se=3/5; sqrt(1-Se^2)=4/5.
    near(state.effective_liquid_saturation(), 0.6);
    near(state.liquid_relative_permeability(), 0.6_f64.sqrt() * 0.04);
    near(pressure(&state), 800.0 * 4.0 / 3.0);
    assert_eq!(state.law(), &calibration);
    // Connectivity is an explicit independent calibration, not a hidden 1/2.
    let ell0 = law(0.1, 0.1, 0.0, 800.0).sample(0.58).unwrap();
    let ell2 = law(0.1, 0.1, 2.0, 800.0).sample(0.58).unwrap();
    near(
        ell2.liquid_relative_permeability() / ell0.liquid_relative_permeability(),
        0.36,
    );
    assert_eq!(pressure(&ell2), pressure(&ell0));
}

#[test]
fn forward_inverse_retention_uses_pressure_not_pressure_head() {
    let first = law(0.05, 0.1, 0.5, 1500.0);
    let twice = law(0.05, 0.1, 0.5, 3000.0);
    let mut previous = 1.0;
    for pc in [0.0, 10.0, 100.0, 1000.0, 10000.0, 1e6] {
        let state = first.at_capillary_pressure(pc).unwrap();
        let effective = 1.0 / (1.0 + (pc / 1500.0).powi(2)).sqrt();
        near(state.liquid_saturation(), 0.05 + 0.85 * effective);
        if pc > 0.0 {
            near(pressure(&state), pc);
        }
        let scaled = twice.sample(state.liquid_saturation()).unwrap();
        near(pressure(&scaled), 2.0 * pressure(&state));
        assert_eq!(
            scaled.liquid_relative_permeability(),
            state.liquid_relative_permeability()
        );
        assert!(state.liquid_relative_permeability() <= previous);
        previous = state.liquid_relative_permeability();
    }
}

#[test]
fn dry_cancellation_and_near_wet_pressure_remain_positive() {
    let calibration = law(0.0, 0.0, 0.5, 1000.0);
    for saturation in [
        1e-30,
        1e-8,
        1e-3,
        0.5,
        0.999999,
        f64::from_bits(1.0_f64.to_bits() - 1),
    ] {
        let state = calibration.sample(saturation).unwrap();
        let gap = ((1.0 - saturation) * (1.0 + saturation)).sqrt();
        let integral = saturation * saturation / (1.0 + gap);
        near(
            state.liquid_relative_permeability(),
            saturation.sqrt() * integral * integral,
        );
        near(pressure(&state), 1000.0 * gap / saturation);
        assert!(state.liquid_relative_permeability() > 0.0);
        assert!(pressure(&state) > 0.0);
    }
}

#[test]
fn exact_endpoint_extensions_and_calibration_refusals() {
    let calibration = law(0.1, 0.2, 0.5, 1000.0);
    for sl in [0.0, 0.1] {
        let state = calibration.sample(sl).unwrap();
        assert_eq!(state.liquid_relative_permeability(), 0.0);
        assert_eq!(
            state.capillary_pressure(),
            LiquidCapillaryPressure::DryLimit
        );
    }
    for sl in [0.8, 1.0] {
        let state = calibration.sample(sl).unwrap();
        assert_eq!(state.liquid_relative_permeability(), 1.0);
        assert_eq!(
            state.capillary_pressure(),
            LiquidCapillaryPressure::Finite(0.0)
        );
    }
    assert!(calibration.sample(-0.01).is_err());
    assert!(calibration.sample(f64::NAN).is_err());
    for pc in [-1.0, f64::NAN, f64::INFINITY, 1e-300, 1e300] {
        assert!(calibration.at_capillary_pressure(pc).is_err());
    }
    let limited = WetPoreLiquidLaw::try_new(
        0.1,
        0.2,
        0.5,
        0.5,
        1000.0,
        WetPoreLiquidSupport {
            liquid_saturation: (0.2, 0.7),
            ..support()
        },
        "bounded synthetic calibration",
    )
    .unwrap();
    assert!(limited.sample(0.1).is_err());
    assert!(limited.sample(0.8).is_err());
    assert!(limited.at_capillary_pressure(0.0).is_err());
    for (m, ell, scale) in [
        (0.0, 0.5, 1000.0),
        (1.0, 0.5, 1000.0),
        (0.5, -0.1, 1000.0),
        (0.5, 0.5, 0.0),
    ] {
        assert!(
            WetPoreLiquidLaw::try_new(
                0.1,
                0.2,
                m,
                ell,
                scale,
                support(),
                "invalid synthetic calibration"
            )
            .is_err()
        );
    }
}
