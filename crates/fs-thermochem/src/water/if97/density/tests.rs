use super::*;
use fs_qty::{Pressure, Temperature};

fn state(p: f64, t: f64, region: If97Region) -> If97State {
    If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:.3e}"
    );
}

#[test]
fn g1_density_curvature_matches_original_scalar_eos_on_both_branches() {
    use If97Region::{Region1 as L, Region2 as V};
    let cases = [
        (0.1e6, 280.0, L),
        (3e6, 300.0, L),
        (80e6, 300.0, L),
        (5e6, 500.0, L),
        (20e6, 620.0, L),
        (1000.0, 300.0, V),
        (0.1e6, 400.0, V),
        (3e6, 600.0, V),
        (25e6, 700.0, V),
        (1e6, 1000.0, V),
    ];
    for (p, t, region) in cases {
        let original = state(p, t, region);
        let jet = original.isobaric_density_jet().unwrap();
        assert_eq!(jet.density, original.density());
        assert_eq!(
            jet.temperature_derivative.value(),
            -original.density().value() * original.thermal_expansion().value()
        );
        // Independent fourth-order differences of the existing scalar EOS,
        // plus a difference of its independently retained first Gibbs slope.
        let rho = |dt| state(p, t + dt, region).density().value();
        let h = 0.2;
        let second = ((-rho(2.0 * h) + 16.0 * rho(h)) + (16.0 * rho(-h) - rho(-2.0 * h))
            - 30.0 * rho(0.0))
            / (12.0 * h * h);
        let exact = jet.temperature_second_derivative.value();
        near(exact, second, 3e-7 * exact.abs().max(1e-4));
        let slope = |dt| {
            let s = state(p, t + dt, region);
            -s.density().value() * s.thermal_expansion().value()
        };
        let h = 0.1;
        let derivative =
            (-slope(2.0 * h) + 8.0 * slope(h) - 8.0 * slope(-h) + slope(-2.0 * h)) / (12.0 * h);
        near(exact, derivative, 2e-8 * exact.abs().max(1e-4));
    }
}

#[test]
fn g3_dilute_vapor_recovers_ideal_gas_density_derivatives() {
    let p = 1e-4;
    let t = 600.0;
    let jet = state(p, t, If97Region::Region2)
        .isobaric_density_jet()
        .unwrap();
    let rho = p / (R_SPECIFIC_J_KG_K * t);
    near(jet.density.value(), rho, 1e-10 * rho);
    near(
        jet.temperature_derivative.value(),
        -rho / t,
        1e-10 * rho / t,
    );
    near(
        jet.temperature_second_derivative.value(),
        2.0 * rho / (t * t),
        1e-10 * rho / (t * t),
    );
}

#[test]
fn g1_density_maximum_and_saturation_keep_original_phase_derivatives() {
    let cold = state(1e5, 274.0, If97Region::Region1)
        .isobaric_density_jet()
        .unwrap();
    let warm = state(1e5, 280.0, If97Region::Region1)
        .isobaric_density_jet()
        .unwrap();
    assert!(cold.temperature_derivative.value() > 0.0);
    assert!(warm.temperature_derivative.value() < 0.0);
    assert!(cold.temperature_second_derivative.value() < 0.0);
    assert!(warm.temperature_second_derivative.value() < 0.0);
    let pair =
        crate::water::if97::If97SaturationPair::at_temperature(Temperature::new(373.15)).unwrap();
    let liquid = pair.liquid().isobaric_density_jet().unwrap();
    let vapor = pair.vapor().isobaric_density_jet().unwrap();
    assert_eq!(liquid.density, pair.liquid().density());
    assert_eq!(vapor.density, pair.vapor().density());
    assert!(liquid.density.value() > 1000.0 * vapor.density.value());
    assert_ne!(
        liquid.temperature_second_derivative,
        vapor.temperature_second_derivative
    );
}
