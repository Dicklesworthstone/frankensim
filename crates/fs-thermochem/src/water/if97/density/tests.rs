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

fn derivative(f: impl Fn(f64) -> f64, h: f64) -> f64 {
    (-f(2.0 * h) + 8.0 * f(h) - 8.0 * f(-h) + f(-2.0 * h)) / (12.0 * h)
}

#[test]
fn g1_joint_density_jet_matches_original_pressure_and_temperature_responses() {
    use If97Region::{Region1 as L, Region2 as V};
    for (p, t, region) in [
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
    ] {
        let original = state(p, t, region);
        let full = original.density_jet().unwrap();
        let thermal = original.isobaric_density_jet().unwrap();
        assert_eq!(full.density, thermal.density);
        assert_eq!(full.temperature_derivative, thermal.temperature_derivative);
        assert_eq!(
            full.temperature_second_derivative,
            thermal.temperature_second_derivative
        );
        assert_eq!(
            full.pressure_derivative.value(),
            original.density().value() * original.isothermal_compressibility().value()
        );
        let hp = 0.001 * p;
        let scalar_pressure = derivative(|dp| state(p + dp, t, region).density().value(), hp);
        near(
            full.pressure_derivative.value(),
            scalar_pressure,
            3e-8 * scalar_pressure.abs(),
        );
        let original_pressure_slope = |dp, dt| {
            let s = state(p + dp, t + dt, region);
            s.density().value() * s.isothermal_compressibility().value()
        };
        let pressure_second = derivative(|dp| original_pressure_slope(dp, 0.0), hp);
        near(
            full.pressure_second_derivative.value(),
            pressure_second,
            3e-7 * pressure_second.abs().max(1e-20),
        );
        // Check the mixed derivative in both directions using only the native
        // state's independently retained first responses (Maxwell symmetry).
        let mixed = full.pressure_temperature_derivative.value();
        near(
            mixed,
            derivative(|dt| original_pressure_slope(0.0, dt), 0.05),
            3e-7 * mixed.abs().max(1e-14),
        );
        near(
            mixed,
            derivative(
                |dp| {
                    let s = state(p + dp, t, region);
                    -s.density().value() * s.thermal_expansion().value()
                },
                hp,
            ),
            3e-7 * mixed.abs().max(1e-14),
        );
    }
}

#[test]
fn g3_dilute_pressure_curvature_retains_the_original_virial_limit() {
    let t = 600.0;
    // At 600 K, independently evaluate the five i=1 Table 11 terms at
    // tau-0.5 = 0.4 with 60-decimal arithmetic. For q=1+a*pi+O(pi²),
    // rho_pp -> -2*a/(R*T*p_star), rather than cancellation to zero.
    let virial_curvature = 1.455_881_451_149_169_6e-13;
    for p in [1e-4, 1e-40, 1e-140] {
        let jet = state(p, t, If97Region::Region2).density_jet().unwrap();
        near(
            jet.pressure_derivative.value(),
            1.0 / (R_SPECIFIC_J_KG_K * t),
            1e-15,
        );
        near(
            jet.pressure_temperature_derivative.value(),
            -1.0 / (R_SPECIFIC_J_KG_K * t * t),
            1e-18,
        );
        near(
            jet.pressure_second_derivative.value(),
            virial_curvature,
            1e-10 * virial_curvature,
        );
    }
}

#[test]
fn g1_joint_density_jet_retains_the_selected_saturation_tangent() {
    let pair =
        crate::water::if97::If97SaturationPair::at_temperature(Temperature::new(373.15)).unwrap();
    let liquid = pair.liquid().density_jet().unwrap();
    let vapor = pair.vapor().density_jet().unwrap();
    assert_eq!(liquid.density, pair.liquid().density());
    assert_eq!(vapor.density, pair.vapor().density());
    assert!(liquid.pressure_derivative.value() > 0.0);
    assert!(vapor.pressure_derivative.value() > 0.0);
    assert_ne!(
        liquid.pressure_second_derivative,
        vapor.pressure_second_derivative
    );
    assert_ne!(
        liquid.pressure_temperature_derivative,
        vapor.pressure_temperature_derivative
    );
}
