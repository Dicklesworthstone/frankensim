//! Source domains belong to each physical layer, including its actual face.
use super::*;

fn chart(lo: f64, hi: f64, k: f64) -> ConductivityCurve {
    ConductivityCurve::new(vec![[lo, k], [hi, k]]).unwrap()
}

#[test]
fn perfect_contact_needs_only_the_actual_shared_face_in_both_sources() {
    let hot = chart(350.0, 700.0, 2.0);
    let cold = chart(250.0, 500.0, 3.0);
    let model = TwoLayerConduction::new([&hot, &cold], [0.02, 0.03], 0.004, 0.0).unwrap();
    let got = model.evaluate([600.0, 300.0], config(), || Ok(())).unwrap();
    close(got.power_a_to_b_w, 60.0, 1e-9);
    for face in got.interface_temperature_k {
        close(face, 450.0, 1e-10);
    }

    // The only common supported temperature is exactly the physical interface.
    let hot = chart(400.0, 700.0, 2.0);
    let cold = chart(250.0, 400.0, 3.0);
    let model = TwoLayerConduction::new([&hot, &cold], [0.02, 0.015], 0.004, 0.0).unwrap();
    let got = model.evaluate([600.0, 300.0], config(), || Ok(())).unwrap();
    assert_eq!(got.interface_temperature_k, [400.0; 2]);
    close(got.power_a_to_b_w, 80.0, 1e-10);
}

#[test]
fn finite_contact_bridges_disjoint_material_temperature_supports() {
    let hot = chart(500.0, 800.0, 8.0);
    let cold = chart(250.0, 450.0, 2.0);
    let resistance = 0.02 / 8.0 + 0.03 / 2.0 + 0.04;
    let flux = 300.0 / resistance;
    let model = TwoLayerConduction::new([&hot, &cold], [0.02, 0.03], 0.004, 0.04).unwrap();
    let got = model.evaluate([600.0, 300.0], config(), || Ok(())).unwrap();
    close(got.power_a_to_b_w, 0.004 * flux, 1e-9);
    close(
        got.interface_temperature_k[0],
        600.0 - flux * 0.02 / 8.0,
        1e-10,
    );
    close(
        got.interface_temperature_k[1],
        300.0 + flux * 0.03 / 2.0,
        1e-10,
    );
    close(got.endpoint_derivative_w_k[0], 0.004 / resistance, 1e-12);
    close(got.endpoint_derivative_w_k[1], -0.004 / resistance, 1e-12);
    assert!(got.interface_temperature_k[0] > hot.temperature_range_k()[0]);
    assert!(got.interface_temperature_k[1] < cold.temperature_range_k()[1]);
}

#[test]
fn nonlinear_local_sources_keep_integrated_heat_and_physical_derivatives() {
    // k_hot = 2 + 0.01 T; k_cold = 12 - 0.01 T. Choose geometry from
    // independent affine integrals so j=10000 and the faces are [540, 410] K.
    let hot = ConductivityCurve::new(vec![[500.0, 7.0], [700.0, 9.0]]).unwrap();
    let cold = ConductivityCurve::new(vec![[250.0, 9.5], [450.0, 7.5]]).unwrap();
    let model = TwoLayerConduction::new([&hot, &cold], [0.08745, 0.09295], 0.004, 0.013).unwrap();
    let temperatures = [650.0, 300.0];
    let got = model.evaluate(temperatures, config(), || Ok(())).unwrap();
    close(got.power_a_to_b_w, 40.0, 1e-9);
    close(got.interface_temperature_k[0], 540.0, 1e-9);
    close(got.interface_temperature_k[1], 410.0, 1e-9);
    let differential_resistance = 0.08745 / 7.4 + 0.09295 / 7.9 + 0.013;
    close(
        got.endpoint_derivative_w_k[0],
        0.004 / differential_resistance * 8.5 / 7.4,
        1e-12,
    );
    close(
        got.endpoint_derivative_w_k[1],
        -0.004 / differential_resistance * 9.0 / 7.9,
        1e-12,
    );
    for i in 0..2 {
        let mut plus = temperatures;
        let mut minus = temperatures;
        plus[i] += 0.001;
        minus[i] -= 0.001;
        let finite_difference = (model
            .evaluate(plus, config(), || Ok(()))
            .unwrap()
            .power_a_to_b_w
            - model
                .evaluate(minus, config(), || Ok(()))
                .unwrap()
                .power_a_to_b_w)
            / 0.002;
        close(got.endpoint_derivative_w_k[i], finite_difference, 1e-7);
    }
    let reversed = TwoLayerConduction::new([&cold, &hot], [0.09295, 0.08745], 0.004, 0.013)
        .unwrap()
        .evaluate([300.0, 650.0], config(), || Ok(()))
        .unwrap();
    assert_eq!(reversed.power_a_to_b_w, -got.power_a_to_b_w);
    assert_eq!(
        reversed.interface_temperature_k,
        [
            got.interface_temperature_k[1],
            got.interface_temperature_k[0]
        ]
    );
}

#[test]
fn unsupported_actual_interfaces_and_own_bulk_states_still_refuse() {
    let hot = chart(500.0, 800.0, 8.0);
    let cold = chart(250.0, 450.0, 2.0);
    // The weaker contact would require the cold material's face above 450 K.
    // Perfect contact is also impossible when the original supports are disjoint.
    for resistance in [0.0, 0.01] {
        let model =
            TwoLayerConduction::new([&hot, &cold], [0.02, 0.03], 0.004, resistance).unwrap();
        assert_eq!(
            model.evaluate([600.0, 300.0], config(), || Ok(())),
            Err(DepositionError::MaterialDomain)
        );
    }
    let model = TwoLayerConduction::new([&hot, &cold], [0.02, 0.03], 0.004, 0.04).unwrap();
    for temperatures in [[499.0, 300.0], [600.0, 451.0], [f64::NAN, 300.0]] {
        assert_eq!(
            model.evaluate(temperatures, config(), || Ok(())),
            Err(DepositionError::MaterialDomain)
        );
    }
    // No positive flux can leave a hot layer already at its lower source limit.
    assert_eq!(
        model.evaluate([500.0, 300.0], config(), || Ok(())),
        Err(DepositionError::MaterialDomain)
    );
}
