//! G0/G1/G3: constitutive limits, independent tangents and retained sources.

use super::*;
use crate::visco::temperature::TemperatureShiftModel;
use std::cell::Cell;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    let scale = actual.abs().max(expected.abs());
    if scale == 0.0 {
        return;
    }
    assert!(
        (actual / scale - expected / scale).abs() <= relative,
        "{actual} != {expected}",
    );
}

fn model(n: f64, a: f64, mode: ViscosityShiftMode) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(12.0, 0.2, 0.6, n, a, [0.0, 1e10]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        mode,
    )
}

#[test]
fn g0_reference_limits_and_carreau_special_case_match_the_declared_equation() {
    let curve = CarreauYasuda::new(12.0, 0.2, 0.6, 0.45, 2.0, [0.0, 1e10]).unwrap();
    assert_eq!(curve.viscosity_pa_s(0.0).unwrap(), 12.0);
    for rate in [1e-8, 0.1, 1.0, 20.0, 1e8] {
        let expected = 0.2 + 11.8 * (1.0 + (0.6_f64 * rate).powi(2)).powf((0.45 - 1.0) / 2.0);
        close(curve.viscosity_pa_s(rate).unwrap(), expected, 2e-13);
    }
    for (time, n, infinity) in [(0.0, 0.45, 0.2), (0.6, 1.0, 0.2), (0.6, 0.45, 12.0)] {
        let constant = CarreauYasuda::new(12.0, infinity, time, n, 0.4, [0.0, 1e10]).unwrap();
        for rate in [0.0, 1.0, 1e10] {
            assert_eq!(constant.viscosity_pa_s(rate).unwrap(), 12.0);
        }
    }
}

#[test]
fn g1_analytic_shear_stress_and_temperature_tangents_match_finite_differences() {
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        for a in [0.6, 1.0, 2.0, 4.0] {
            let law = model(0.45, a, mode);
            for rate in [0.1_f64, 2.0, 20.0] {
                let temperature = 320.0;
                let state = law.evaluate(temperature, rate).unwrap();
                let step = 1e-4 * rate.max(1.0);
                let lower = law
                    .evaluate(temperature, rate - step)
                    .unwrap()
                    .dynamic_viscosity_pa_s();
                let upper = law
                    .evaluate(temperature, rate + step)
                    .unwrap()
                    .dynamic_viscosity_pa_s();
                close(
                    state.shear_derivative_pa_s2().unwrap(),
                    (upper - lower) / (2.0 * step),
                    2e-6,
                );
                close(
                    state.stress_tangent_pa_s(),
                    (upper * (rate + step) - lower * (rate - step)) / (2.0 * step),
                    2e-6,
                );
                close(
                    state.scaled_shear_derivative_pa_s(),
                    rate * state.shear_derivative_pa_s2().unwrap(),
                    2e-12,
                );
                let lower = law
                    .evaluate(temperature - 1e-3, rate)
                    .unwrap()
                    .dynamic_viscosity_pa_s();
                let upper = law
                    .evaluate(temperature + 1e-3, rate)
                    .unwrap()
                    .dynamic_viscosity_pa_s();
                close(
                    state.temperature_derivative_pa_s_per_k(),
                    (upper - lower) / 2e-3,
                    2e-6,
                );
                assert!(state.stress_tangent_pa_s() > 0.0);
                assert!(state.temperature_derivative_pa_s_per_k() < 0.0);
            }
        }
    }
}

#[test]
fn g0_zero_rate_singular_scalar_slope_does_not_invent_a_stress_singularity() {
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        for a in [0.5, 1.0, 2.0] {
            let law = model(0.45, a, mode);
            let state = law.evaluate(320.0, 0.0).unwrap();
            assert_eq!(state.scaled_shear_derivative_pa_s(), 0.0);
            assert_eq!(state.stress_tangent_pa_s(), state.dynamic_viscosity_pa_s());
            if a < 1.0 {
                assert_eq!(state.shear_derivative_pa_s2(), None);
            } else if a == 1.0 {
                let shift = state.temperature().factor();
                let horizontal = if mode == ViscosityShiftMode::VerticalAndRate {
                    shift
                } else {
                    1.0
                };
                close(
                    state.shear_derivative_pa_s2().unwrap(),
                    shift * horizontal * (0.45 - 1.0) * 11.8 * 0.6,
                    2e-13,
                );
            } else {
                assert_eq!(state.shear_derivative_pa_s2(), Some(0.0));
            }
        }
    }
}

#[test]
fn g3_temperature_mode_moves_reference_support_without_clipping_the_curve() {
    let curve = CarreauYasuda::new(10.0, 0.1, 1.0, 0.5, 2.0, [0.0, 10.0]).unwrap();
    let shift = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 1.0,
            c2_k: 100.0,
        },
        300.0,
        [300.0, 400.0],
    )
    .unwrap();
    let vertical = ThermalCarreauYasuda::new(curve, shift, ViscosityShiftMode::Vertical);
    let combined = ThermalCarreauYasuda::new(curve, shift, ViscosityShiftMode::VerticalAndRate);
    assert!(vertical.evaluate(400.0, 20.0).is_err());
    let state = combined.evaluate(400.0, 20.0).unwrap();
    let factor = 10.0_f64.powf(-0.5);
    close(state.reduced_rate_per_s(), factor * 20.0, 2e-14);
    close(
        state.dynamic_viscosity_pa_s(),
        factor * curve.viscosity_pa_s(factor * 20.0).unwrap(),
        2e-13,
    );
    assert!(combined.evaluate(401.0, 1.0).is_err());
    assert!(combined.evaluate(300.0, 10.01).is_err());
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(combined.evaluate(300.0, invalid).is_err());
    }
}

#[test]
fn g0_extreme_products_preserve_representable_viscosity_and_shifted_slopes() {
    let curve = CarreauYasuda::new(1e300, 0.0, 1e300, 0.3, 2.0, [0.0, 1e300]).unwrap();
    let law = ThermalCarreauYasuda::new(
        curve,
        *model(0.5, 2.0, ViscosityShiftMode::Vertical).temperature_shift(),
        ViscosityShiftMode::Vertical,
    );
    let state = law.evaluate(300.0, 1e300).unwrap();
    // lambda*g overflows, and the shape underflows, but mu_0*shape is finite.
    let expected = (1e300_f64.ln() + (0.3 - 1.0) * (1e300_f64.ln() + 1e300_f64.ln())).exp();
    close(state.dynamic_viscosity_pa_s(), expected, 2e-12);
    close(state.stress_tangent_pa_s(), 0.3 * expected, 2e-12);
    assert_eq!(state.shear_derivative_pa_s2(), None); // Nonzero scalar derivative underflows.

    let curve = CarreauYasuda::new(1.0, 0.1, 1.0, 0.5, 2.0, [0.0, 1.0]).unwrap();
    let shift = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 300.0,
            c2_k: 100.0,
        },
        300.0,
        [250.0, 300.0],
    )
    .unwrap();
    let state = ThermalCarreauYasuda::new(curve, shift, ViscosityShiftMode::Vertical)
        .evaluate(250.0, 1e-300)
        .unwrap();
    // The reference g*mu' is below f64, but its thermal shift is representable.
    close(state.shear_derivative_pa_s2().unwrap(), -0.45, 2e-12);
    close(state.scaled_shear_derivative_pa_s(), -0.45e-300, 2e-12);
}

#[derive(Debug)]
struct CountedFluid {
    state: FluidTransportState,
    calls: Cell<usize>,
}
impl FluidTransport for CountedFluid {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.calls.set(self.calls.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn g0_binding_preserves_current_thermal_source_and_replaces_viscosity_once() {
    let source = CountedFluid {
        state: FluidTransportState::try_new(7.0, 11.0, 13.0, 101.0).unwrap(),
        calls: Cell::new(0),
    };
    let state = model(0.45, 2.0, ViscosityShiftMode::VerticalAndRate)
        .evaluate(320.0, 2.0)
        .unwrap();
    let bound = state.bind(&source).unwrap();
    assert!(std::ptr::eq(*bound.source(), &source));
    assert_eq!(bound.original_transport(), &source.state);
    assert_eq!(bound.rheology(), &state);
    for _ in 0..2 {
        let actual = bound.transport_state().unwrap();
        assert_eq!(actual.density_kg_m3(), 7.0);
        assert_eq!(actual.specific_heat_j_kg_k(), 11.0);
        assert_eq!(actual.thermal_conductivity_w_m_k(), 13.0);
        assert_eq!(
            actual.dynamic_viscosity_pa_s(),
            state.dynamic_viscosity_pa_s()
        );
        assert_eq!(
            actual.kinematic_viscosity_m2_s(),
            state.dynamic_viscosity_pa_s() / 7.0
        );
    }
    assert_eq!(source.calls.get(), 1);

    #[derive(Debug)]
    struct Refusing;
    impl FluidTransport for Refusing {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            Err(MaterialError::State {
                what: "source pressure outside support".into(),
            })
        }
    }
    assert_eq!(
        state.bind(Refusing).unwrap_err(),
        MaterialError::State {
            what: "source pressure outside support".into()
        },
    );
}

#[test]
fn g0_invalid_calibration_parameters_refuse() {
    for (zero, infinity, time, n, a, band) in [
        (0.0, 0.0, 1.0, 0.5, 2.0, [0.0, 10.0]),
        (1.0, 2.0, 1.0, 0.5, 2.0, [0.0, 10.0]),
        (1.0, -0.1, 1.0, 0.5, 2.0, [0.0, 10.0]),
        (1.0, 0.0, -1.0, 0.5, 2.0, [0.0, 10.0]),
        (1.0, 0.0, 1.0, 0.0, 2.0, [0.0, 10.0]),
        (1.0, 0.0, 1.0, 1.1, 2.0, [0.0, 10.0]),
        (1.0, 0.0, 1.0, 0.5, 0.0, [0.0, 10.0]),
        (1.0, 0.0, 1.0, 0.5, 2.0, [-1.0, 10.0]),
        (1.0, 0.0, 1.0, 0.5, 2.0, [10.0, 1.0]),
        (1.0, 0.0, 1.0, 0.5, 2.0, [0.0, f64::INFINITY]),
    ] {
        assert!(CarreauYasuda::new(zero, infinity, time, n, a, band).is_err());
    }
}
