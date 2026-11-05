//! G0/G1/G3: analytic flow stresses, material derivatives and objectivity.
use super::*;
use crate::fluid_transport::FluidTransportState;
use crate::fluid_transport::rheology::{CarreauYasuda, ViscosityShiftMode};
use crate::tensor::{matmul3, rotate, rotation};
use crate::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}");
}
fn model(a: f64, mode: ViscosityShiftMode) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(12.0, 0.2, 0.6, 0.45, a, [0.0, 1e5]).unwrap(),
        TemperatureShift::new(TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 12_000.0,
        }, 300.0, [280.0, 360.0]).unwrap(), mode,
    )
}

#[test]
fn g1_couette_stress_traction_and_heat_use_dynamic_viscosity() {
    let source = FluidTransportState::try_new(930.0, 1900.0, 0.15, 0.7).unwrap();
    let law = NewtonianViscosity::from_transport(&source, 3.0).unwrap();
    for tangent in 0..3 {
        for normal in 0..3 {
            if tangent == normal { continue; }
            let rate: f64 = -4.0;
            let mut gradient = [[0.0; 3]; 3];
            gradient[tangent][normal] = rate;
            let response = law.evaluate(DeformationRate::from_velocity_gradient(gradient).unwrap()).unwrap();
            close(response.rate().shear_rate_per_s(), rate.abs(), 1e-14);
            close(response.dissipation_w_m3(), 0.7 * rate * rate, 1e-14);
            let mut n = [0.0; 3];
            n[normal] = 1.0;
            let traction = response.traction(n).unwrap();
            close(traction.vector_pa[tangent], 0.7 * rate, 1e-14);
            assert_eq!(traction.shear_pa, traction.vector_pa);
            assert_eq!(traction.normal_pa, 0.0);
            let reversed = response.traction(n.map(|x| -x)).unwrap();
            assert_eq!(reversed.vector_pa, traction.vector_pa.map(|x| -x));
            assert_eq!(reversed.shear_pa, traction.shear_pa.map(|x| -x));
            assert_eq!(reversed.normal_pa, traction.normal_pa);
        }
    }
}

#[test]
fn g1_rigid_spin_dilation_and_extension_are_distinct() {
    let law = NewtonianViscosity::new(2.0, 5.0).unwrap();
    let spin = [[0.0, -3.0, 2.0], [3.0, 0.0, -7.0], [-2.0, 7.0, 0.0]];
    let response = law.evaluate(DeformationRate::from_velocity_gradient(spin).unwrap()).unwrap();
    assert_eq!(*response.stress_pa(), [0.0; 6]);
    assert_eq!(response.dissipation_w_m3(), 0.0);
    let expansion = law.evaluate(DeformationRate::from_tensor([4.0, 4.0, 4.0, 0.0, 0.0, 0.0]).unwrap()).unwrap();
    assert_eq!(expansion.rate().shear_rate_per_s(), 0.0);
    assert_eq!(*expansion.stress_pa(), [60.0, 60.0, 60.0, 0.0, 0.0, 0.0]);
    assert_eq!(expansion.dissipation_w_m3(), 720.0);
    let stokes = NewtonianViscosity::new(2.0, 0.0).unwrap().evaluate(*expansion.rate()).unwrap();
    assert_eq!(*stokes.stress_pa(), [0.0; 6]);
    assert_eq!(stokes.dissipation_w_m3(), 0.0);
    // Incompressible extension has a nonzero invariant despite zero xy shear.
    let extension = law.evaluate(DeformationRate::from_tensor([3.0, -1.5, -1.5, 0.0, 0.0, 0.0]).unwrap()).unwrap();
    close(extension.rate().shear_rate_per_s(), 3.0_f64.sqrt() * 3.0, 1e-14);
    assert_eq!(*extension.stress_pa(), [12.0, -6.0, -6.0, 0.0, 0.0, 0.0]);
    assert_eq!(extension.dissipation_w_m3(), 54.0);
}

#[test]
fn g3_stress_rotates_and_heating_is_objective() {
    let gradient = [[0.7, 2.0, -1.2], [0.5, -0.2, 0.3], [1.1, -0.4, 0.6]];
    let q = rotation([0.0, 0.0, 1.0], 0.63);
    let qt = core::array::from_fn(|i| core::array::from_fn(|j| q[j][i]));
    let rotated_gradient = matmul3(&matmul3(&q, &gradient), &qt);
    let law = model(1.7, ViscosityShiftMode::VerticalAndRate);
    let original = law.viscous_response(315.0, DeformationRate::from_velocity_gradient(gradient).unwrap(), 0.8).unwrap();
    let rotated = law.viscous_response(315.0, DeformationRate::from_velocity_gradient(rotated_gradient).unwrap(), 0.8).unwrap();
    for (actual, expected) in rotated.viscous.stress_pa().iter().zip(rotate(original.viscous.stress_pa(), &q)) {
        close(*actual, expected, 5e-13);
    }
    close(rotated.viscous.dissipation_w_m3(), original.viscous.dissipation_w_m3(), 5e-13);
    let mut with_spin = gradient;
    with_spin[0][1] -= 17.0;
    with_spin[1][0] += 17.0;
    let spun = law.viscous_response(315.0, DeformationRate::from_velocity_gradient(with_spin).unwrap(), 0.8).unwrap();
    for (actual, expected) in spun.viscous.stress_pa().iter().zip(original.viscous.stress_pa()) {
        close(*actual, *expected, 5e-13);
    }
}

#[test]
fn g1_full_coupled_tensor_and_temperature_tangents_match_finite_differences() {
    let tensor = [1.2, -0.3, 0.6, -1.7, 0.8, 2.1];
    for mode in [ViscosityShiftMode::Vertical, ViscosityShiftMode::VerticalAndRate] {
        for transition in [0.5, 2.0] {
            let law = model(transition, mode);
            let response = law.viscous_response(320.0, DeformationRate::from_tensor(tensor).unwrap(), 0.4).unwrap();
            let tangent = response.viscous.tangent_pa_s().unwrap();
            let h = 1e-5;
            for column in 0..6 {
                let mut upper = tensor;
                let mut lower = tensor;
                upper[column] += h;
                lower[column] -= h;
                let upper = law.viscous_response(320.0, DeformationRate::from_tensor(upper).unwrap(), 0.4).unwrap();
                let lower = law.viscous_response(320.0, DeformationRate::from_tensor(lower).unwrap(), 0.4).unwrap();
                for row in 0..6 {
                    close(tangent[row][column], (upper.viscous.stress_pa()[row] - lower.viscous.stress_pa()[row]) / (2.0 * h), 2e-9);
                    let weight_row = if row < 3 { 1.0 } else { 2.0 };
                    let weight_column = if column < 3 { 1.0 } else { 2.0 };
                    close(weight_row * tangent[row][column], weight_column * tangent[column][row], 2e-13);
                }
            }
            let upper = law.viscous_response(320.0 + h, DeformationRate::from_tensor(tensor).unwrap(), 0.4).unwrap();
            let lower = law.viscous_response(320.0 - h, DeformationRate::from_tensor(tensor).unwrap(), 0.4).unwrap();
            for row in 0..6 {
                close(response.temperature_tangent_pa_per_k[row], (upper.viscous.stress_pa()[row] - lower.viscous.stress_pa()[row]) / (2.0 * h), 3e-9);
            }
        }
    }
}

#[test]
fn g0_zero_shear_has_finite_tangent_even_with_singular_viscosity_slope() {
    let law = model(0.5, ViscosityShiftMode::VerticalAndRate);
    let response = law.viscous_response(310.0, DeformationRate::from_tensor([2.0, 2.0, 2.0, 0.0, 0.0, 0.0]).unwrap(), 0.0).unwrap();
    assert_eq!(response.rheology.shear_derivative_pa_s2(), None);
    assert_eq!(*response.viscous.stress_pa(), [0.0; 6]);
    assert_eq!(response.temperature_tangent_pa_per_k, [0.0; 6]);
    let reference = NewtonianViscosity::new(response.rheology.dynamic_viscosity_pa_s(), 0.0).unwrap()
        .evaluate(*response.viscous.rate()).unwrap();
    assert_eq!(response.viscous.tangent_pa_s().unwrap(), reference.tangent_pa_s().unwrap());
}

#[test]
fn g0_bulk_stress_does_not_pollute_resolved_surface_shear() {
    let rate = DeformationRate::from_tensor([1.0, 1.0, 1.0, 2.0, -3.0, 4.0]).unwrap();
    let shear = NewtonianViscosity::new(0.25, 0.0).unwrap().evaluate(rate).unwrap();
    let bulk = NewtonianViscosity::new(0.25, 1e16).unwrap().evaluate(rate).unwrap();
    let n = [0.6, 0.8, 0.0];
    let before = shear.traction(n).unwrap();
    let after = bulk.traction(n).unwrap();
    assert_eq!(before.shear_pa, after.shear_pa);
    close((0..3).map(|i| n[i] * after.shear_pa[i]).sum(), 0.0, 1e-14);
    close(contract(shear.stress_pa(), rate.tensor_per_s()), shear.dissipation_w_m3(), 1e-14);
}

#[test]
fn g0_source_support_and_invalid_states_refuse_without_clipping() {
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(NewtonianViscosity::new(1.0, bad).is_err());
        assert!(NewtonianViscosity::new(bad, 0.0).is_err());
    }
    assert!(NewtonianViscosity::new(0.0, 0.0).is_err());
    assert!(DeformationRate::from_tensor([f64::NAN; 6]).is_err());
    let law = model(2.0, ViscosityShiftMode::Vertical);
    let zero = DeformationRate::from_tensor([0.0; 6]).unwrap();
    assert!(law.viscous_response(400.0, zero, 0.0).is_err());
    let high = DeformationRate::from_tensor([0.0, 0.0, 0.0, 1e6, 0.0, 0.0]).unwrap();
    assert!(law.viscous_response(300.0, high, 0.0).is_err());
    let response = NewtonianViscosity::new(1.0, 0.0).unwrap().evaluate(zero).unwrap();
    for normal in [[0.0; 3], [2.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]] {
        assert!(response.traction(normal).is_err());
    }
}
