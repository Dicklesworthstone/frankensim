//! Source-consistent real-gas mechanical forces, tangents and work.

use super::{entropy_controls, parcel_at};
use crate::fluid_acoustics::FluidAcoustics;

#[test]
fn g1_dense_gas_bulk_stiffness_matches_adiabatic_pressure_differences() {
    let mut gas = parcel_at(5.0e6);
    gas.set_volume_isentropic(0.0075, entropy_controls(), || false)
        .unwrap();
    let response = gas.pressure_response().unwrap();
    let acoustic = gas.state().acoustic_state().unwrap();
    assert_eq!(response.pressure_pa, gas.state().pressure_pa());
    assert_eq!(
        response.bulk_modulus_pa,
        acoustic.density_kg_m3() * acoustic.sound_speed_m_s() * acoustic.sound_speed_m_s()
    );
    // A dense real-fluid Cp/Cv is not the ideal-gas compressibility exponent.
    assert!(
        (response.bulk_modulus_pa - gas.state().gamma() * gas.state().pressure_pa()).abs()
            > 100_000.0
    );
    let h = gas.volume_m3() * 1.0e-5;
    let mut plus = gas.clone();
    let mut minus = gas.clone();
    plus.set_volume_isentropic(gas.volume_m3() + h, entropy_controls(), || false)
        .unwrap();
    minus
        .set_volume_isentropic(gas.volume_m3() - h, entropy_controls(), || false)
        .unwrap();
    let derivative = (plus.state().pressure_pa() - minus.state().pressure_pa()) / (2.0 * h);
    assert!((derivative / response.pressure_volume_derivative_pa_m3 - 1.0).abs() < 5.0e-6);
}

#[test]
fn g1_dense_gas_pressure_pullback_includes_the_geometric_tangent() {
    let gas = parcel_at(5.0e6);
    // V(q)=0.25 m*q^2: at q=0.2 m, the chamber volume is 0.01 m^3.
    let q = 0.2;
    let gradient = [0.5 * q];
    let direction = [0.003];
    let exterior = 200_000.0;
    let load = gas
        .mechanical_load(gas.volume_m3(), &gradient, exterior)
        .unwrap();
    assert_eq!(load.state(), gas.state());
    let force = load.generalized_forces().unwrap();
    assert_eq!(
        force[0],
        (gas.state().pressure_pa() - exterior) * gradient[0]
    );
    let power = load.mechanical_power_w(&direction).unwrap();
    assert!((power - force[0] * direction[0]).abs() < 1.0e-10);
    let action = load
        .tangent_action(&direction, &[0.5 * direction[0]])
        .unwrap();
    let h = 1.0e-3;
    let force_at = |sign: f64| {
        let position = q + sign * h * direction[0];
        let mut moved = gas.clone();
        moved
            .set_volume_isentropic(0.25 * position * position, entropy_controls(), || false)
            .unwrap();
        (moved.state().pressure_pa() - exterior) * (0.5 * position)
    };
    let derivative = (force_at(1.0) - force_at(-1.0)) / (2.0 * h);
    assert!((derivative / action[0] - 1.0).abs() < 5.0e-6);
    assert!(gas.mechanical_load(0.011, &gradient, exterior).is_err());
}

#[test]
fn g3_integrated_generalized_force_matches_compression_and_exterior_work() {
    let initial = parcel_at(5.0e6);
    let final_volume = 0.0075;
    let exterior = 200_000.0;
    let mut compressed = initial.clone();
    let report = compressed
        .set_volume_isentropic(final_volume, entropy_controls(), || false)
        .unwrap();
    let expected = report.work_by_gas_j - exterior * (final_volume - initial.volume_m3());
    let q0 = 0.2;
    let q1 = fs_math::det::sqrt(4.0 * final_volume);
    let force_at = |q: f64| {
        let mut gas = initial.clone();
        if q != q0 {
            gas.set_volume_isentropic(0.25 * q * q, entropy_controls(), || false)
                .unwrap();
        }
        let gradient = [0.5 * q];
        gas.mechanical_load(gas.volume_m3(), &gradient, exterior)
            .unwrap()
            .generalized_forces()
            .unwrap()[0]
    };
    // Independent geometric-coordinate quadrature converges to caloric work;
    // it does not reuse the algebraic Delta U + (-Delta U) receipt.
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let dq = (q1 - q0) / n as f64;
        let mut sum = force_at(q0) + force_at(q1);
        for i in 1..n {
            sum += if i % 2 == 0 { 2.0 } else { 4.0 } * force_at(q0 + i as f64 * dq);
        }
        errors.push((sum * dq / 3.0 - expected).abs());
    }
    assert!(
        errors[0] > 12.0 * errors[1] && errors[1] > 12.0 * errors[2],
        "{errors:?}"
    );
    assert!(errors[2] < 0.001, "{errors:?}");
}
