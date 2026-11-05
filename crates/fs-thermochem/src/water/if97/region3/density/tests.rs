use super::*;
use crate::water::if97::region3::pressure::{DensePressureSolve, If97DenseSaturationPair};
use crate::water::if97::regions::{If97FluidRegion, If97FluidState};
use fs_qty::{Density, Pressure, Temperature};

fn source(rho: f64, t: f64) -> If97DenseState {
    If97DenseState::at_density_temperature(Density::new(rho), Temperature::new(t)).unwrap()
}
fn controls() -> DensePressureSolve {
    DensePressureSolve {
        pressure_tolerance_pa: 1e-3,
        density_tolerance_kg_m3: 1e-8,
        maximum_iterations: 128,
    }
}
fn pressure_source(p: f64, t: f64) -> If97DenseState {
    If97DenseState::at_pressure_temperature(
        Pressure::new(p),
        Temperature::new(t),
        controls(),
        || false,
    )
    .unwrap()
    .state
}
fn derivative(f: impl Fn(f64) -> f64, h: f64) -> f64 {
    (-f(2.0 * h) + 8.0 * f(h) - 8.0 * f(-h) + f(-2.0 * h)) / (12.0 * h)
}
fn close(actual: f64, expected: f64, relative: f64, floor: f64) {
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= relative * expected.abs().max(floor),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

#[test]
fn g1_dense_pressure_hessian_matches_original_stiffness_and_thermal_response() {
    // Three original Table 33 states plus stable dense liquid and near-critical
    // water. Independent fourth-order differences use the unchanged source's
    // already retained first derivatives, not the new Hessian implementation.
    for (rho, t) in [
        (500.0, 650.0),
        (200.0, 650.0),
        (500.0, 750.0),
        (600.0, 640.0),
        (322.0, 647.2),
    ] {
        let [rr, rt, tt] = pressure_hessian(rho, t);
        close(
            rr,
            derivative(
                |dr| source(rho + dr, t).dp_drho_at_temperature().value(),
                0.05,
            ),
            2e-6,
            1.0,
        );
        close(
            rt,
            derivative(
                |dt| source(rho, t + dt).dp_drho_at_temperature().value(),
                0.002,
            ),
            2e-6,
            1.0,
        );
        close(
            rt,
            derivative(|dr| source(rho + dr, t).dp_dt_at_density().value(), 0.05),
            2e-6,
            1.0,
        );
        close(
            tt,
            derivative(|dt| source(rho, t + dt).dp_dt_at_density().value(), 0.002),
            2e-6,
            1.0,
        );
    }
}

#[test]
fn g1_dense_density_hessian_matches_independent_pressure_solutions() {
    for (rho, t) in [
        (500.0, 650.0),
        (200.0, 650.0),
        (500.0, 750.0),
        (600.0, 640.0),
    ] {
        let original = source(rho, t);
        let p = original.pressure().value();
        let jet = original.density_jet().unwrap();
        let rho_p = |dp, dt| {
            pressure_source(p + dp, t + dt)
                .dp_drho_at_temperature()
                .value()
                .recip()
        };
        let rho_t = |dp, dt| {
            let s = pressure_source(p + dp, t + dt);
            -s.dp_dt_at_density().value() / s.dp_drho_at_temperature().value()
        };
        close(
            jet.pressure_derivative.value(),
            derivative(|dp| pressure_source(p + dp, t).density().value(), 200.0),
            2e-6,
            1e-10,
        );
        close(
            jet.temperature_derivative.value(),
            derivative(|dt| pressure_source(p, t + dt).density().value(), 0.005),
            2e-6,
            1e-5,
        );
        close(
            jet.pressure_second_derivative.value(),
            derivative(|dp| rho_p(dp, 0.0), 200.0),
            3e-6,
            1e-16,
        );
        close(
            jet.temperature_second_derivative.value(),
            derivative(|dt| rho_t(0.0, dt), 0.005),
            3e-6,
            1e-6,
        );
        close(
            jet.pressure_temperature_derivative.value(),
            derivative(|dt| rho_p(0.0, dt), 0.005),
            3e-6,
            1e-12,
        );
        close(
            jet.pressure_temperature_derivative.value(),
            derivative(|dp| rho_t(dp, 0.0), 200.0),
            3e-6,
            1e-12,
        );
    }
}

#[test]
fn g0_regional_density_jets_retain_each_original_source_and_solve_diagnostics() {
    for (p, t, region) in [
        (3e6, 300.0, If97FluidRegion::Region1),
        (3_500.0, 700.0, If97FluidRegion::Region2),
        (25e6, 650.0, If97FluidRegion::Region3),
    ] {
        let state = If97FluidState::at_pressure_temperature(
            Pressure::new(p),
            Temperature::new(t),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(state.region(), region);
        let before = state;
        let full = state.density_jet().unwrap();
        let thermal = state.isobaric_density_jet().unwrap();
        let expected = match state {
            If97FluidState::Gibbs(s) => s.density_jet().unwrap(),
            If97FluidState::Dense(s) => s.state.density_jet().unwrap(),
        };
        assert_eq!(full, expected);
        assert_eq!(full.density, state.density());
        assert_eq!(thermal.density, full.density);
        assert_eq!(thermal.temperature_derivative, full.temperature_derivative);
        assert_eq!(
            thermal.temperature_second_derivative,
            full.temperature_second_derivative
        );
        assert_eq!(state, before);
    }
}

#[test]
fn g3_critical_amplification_and_coexistence_tangents_are_not_smoothed() {
    let regular = source(500.0, 650.0).density_jet().unwrap();
    let critical = source(322.0, 647.2).density_jet().unwrap();
    assert!(critical.pressure_derivative.value() > 100.0 * regular.pressure_derivative.value());
    assert!(
        critical.temperature_derivative.value().abs()
            > 100.0 * regular.temperature_derivative.value().abs()
    );
    let pair =
        If97DenseSaturationPair::at_temperature(Temperature::new(640.0), controls(), || false)
            .unwrap();
    let liquid = *pair.liquid();
    let vapor = *pair.vapor();
    for s in [liquid, vapor] {
        let jet = s.density_jet().unwrap();
        assert_eq!(jet.density, s.density());
        assert_eq!(
            jet.pressure_derivative.value(),
            s.dp_drho_at_temperature().value().recip()
        );
        assert_eq!(
            jet.temperature_derivative.value(),
            -s.dp_dt_at_density().value() / s.dp_drho_at_temperature().value()
        );
    }
    assert!(liquid.density().value() > vapor.density().value());
    assert_ne!(liquid.density_jet().unwrap(), vapor.density_jet().unwrap());
}
