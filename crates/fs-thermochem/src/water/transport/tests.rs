use super::*;
use crate::water::if97::region3::If97DenseState;
use crate::water::if97::{If97Region, If97State, saturation_pressure};
use fs_qty::Density;

fn controls() -> DensePressureSolve {
    DensePressureSolve {
        pressure_tolerance_pa: 1e-4,
        density_tolerance_kg_m3: 1e-8,
        maximum_iterations: 128,
    }
}
fn state(p: f64, t: f64) -> If97TransportState {
    If97TransportState::at_pressure_temperature(
        Pressure::new(p),
        Temperature::new(t),
        controls(),
        || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.15e} != {expected:.15e}, allowance {tolerance:e}"
    );
}

#[test]
fn original_viscosity_and_background_conductivity_verification_tables() {
    // Literal R12-08 Table 4, including liquid, compressed liquid and vapor.
    // Densities are supplied by that table; these test the correlation itself,
    // not phase admission at arbitrary rho/T or general critical viscosity.
    for (t, rho, micro_pa_s) in [
        (298.15, 998.0, 889.735100),
        (298.15, 1200.0, 1437.649467),
        (373.15, 1000.0, 307.883622),
        (433.15, 1.0, 14.538324),
        (433.15, 1000.0, 217.685358),
        (873.15, 1.0, 32.619287),
        (873.15, 100.0, 35.802262),
        (873.15, 600.0, 77.430195),
        (1173.15, 1.0, 44.217245),
        (1173.15, 100.0, 47.640433),
        (1173.15, 400.0, 64.154608),
    ] {
        near(viscosity(t, rho) * 1e6, micro_pa_s, 5.1e-7);
    }
    // Literal R15-11 Table 4, whose critical contribution is zero.
    for (t, rho, mw_m_k) in [
        (298.15, 0.0, 18.4341883),
        (298.15, 998.0, 607.712868),
        (298.15, 1200.0, 799.038144),
        (873.15, 0.0, 79.1034659),
    ] {
        near(background_conductivity(t, rho) * 1e3, mw_m_k, 5.1e-7);
    }
}

#[test]
fn original_if97_state_drives_industrial_conductivity_and_critical_enhancement() {
    // R15-11 Tables 7 and 8: actual IF97 Cp, Cv and density derivative feed
    // the INDUSTRIAL law. General-IAPWS95 Table 5 is a different calculation.
    for (p, t, conductivity, critical, viscosity_micro_pa_s) in [
        (20e6, 620.0, 0.481485195, 0.0126391714, 70.9051068),
        (50e6, 620.0, 0.545038940, 0.00575816285, 84.1527945),
        (0.3e6, 650.0, 0.0522311024, 0.000000129246457, 23.4877453),
        (50e6, 800.0, 0.177709914, 0.00664341394, 39.3727534),
    ] {
        let actual = state(p, t);
        near(actual.thermal_conductivity().value(), conductivity, 6e-10);
        near(
            actual.critical_thermal_conductivity().value(),
            critical,
            6e-11,
        );
        near(
            actual.dynamic_viscosity().value() * 1e6,
            viscosity_micro_pa_s,
            6e-8,
        );
        assert_eq!(
            actual.thermal_conductivity().value(),
            actual.background_thermal_conductivity().value()
                + actual.critical_thermal_conductivity().value()
        );
        let original = If97FluidState::at_pressure_temperature(
            Pressure::new(p),
            Temperature::new(t),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(actual.thermodynamics(), &original);
    }
    // R15-11 industrial Table 9 uses the original Region 3 EOS. Invert the
    // pressure computed by that same native EOS, retaining the solver receipt.
    for (rho, conductivity) in [(222.0, 0.366879411), (322.0, 1.24182415)] {
        let original =
            If97DenseState::at_density_temperature(Density::new(rho), Temperature::new(647.35))
                .unwrap();
        let actual = state(original.pressure().value(), 647.35);
        near(actual.thermodynamics().density().value(), rho, 2e-6);
        near(actual.thermal_conductivity().value(), conductivity, 1e-6);
        assert!(actual.critical_thermal_conductivity().value() > 0.1);
    }
}

#[test]
fn pressure_effects_preserve_original_calorics_and_explicit_phase_domain() {
    let low = state(100_000.0, 300.0);
    let high = state(80e6, 300.0);
    assert!(high.thermodynamics().density().value() > low.thermodynamics().density().value());
    assert!(high.thermal_conductivity().value() > low.thermal_conductivity().value());
    assert_ne!(
        high.dynamic_viscosity().value(),
        low.dynamic_viscosity().value()
    );
    assert_ne!(
        high.thermodynamics().sound_speed().value(),
        low.thermodynamics().sound_speed().value()
    );
    let t = Temperature::new(400.0);
    let ps = saturation_pressure(t).unwrap();
    assert!(If97TransportState::at_pressure_temperature(ps, t, controls(), || false).is_err());
    let liquid = If97State::at_pressure_temperature(ps, t, If97Region::Region1).unwrap();
    let vapor = If97State::at_pressure_temperature(ps, t, If97Region::Region2).unwrap();
    let liquid = If97TransportState::from_state(If97FluidState::Gibbs(liquid)).unwrap();
    let vapor = If97TransportState::from_state(If97FluidState::Gibbs(vapor)).unwrap();
    assert!(liquid.dynamic_viscosity().value() > vapor.dynamic_viscosity().value());
    assert!(liquid.thermal_conductivity().value() > vapor.thermal_conductivity().value());
    for (p, t) in [(101e6, 300.0), (1e5, 273.0), (1e5, 1100.0), (1e5, f64::NAN)] {
        assert!(
            If97TransportState::at_pressure_temperature(
                Pressure::new(p),
                Temperature::new(t),
                controls(),
                || false
            )
            .is_err()
        );
    }
    assert_eq!(
        If97TransportState::at_pressure_temperature(
            Pressure::new(1e5),
            Temperature::new(300.0),
            controls(),
            || true
        ),
        Err(WaterTransportError::Cancelled)
    );
    assert_eq!(low, state(100_000.0, 300.0));
}
