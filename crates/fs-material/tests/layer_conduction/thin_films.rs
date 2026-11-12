//! Thin-layer heat continuity must retain the solved relative material drop.
use super::*;
use fs_material::thermal::conduction::{
    ConvectiveLayerError, LayerStackConduction, ReservoirLayerGeometry,
};

#[test]
fn thin_conducting_film_and_insulation_recover_series_heat_in_both_directions() {
    let insulation = constant(0.02);
    let film = constant(10.0);
    let area = 0.004;
    let distances = [0.01, 1e-6];
    for contact in [0.0, 0.05] {
        let resistance = distances[0] / 0.02 + distances[1] / 10.0 + contact;
        let conductance = area / resistance;
        let model =
            TwoLayerConduction::new([&insulation, &film], distances, area, contact).unwrap();
        for temperatures in [[600.0, 300.0], [300.0, 600.0]] {
            let got = model.evaluate(temperatures, config(), || Ok(())).unwrap();
            let expected = conductance * (temperatures[0] - temperatures[1]);
            close(got.power_a_to_b_w, expected, 1e-11);
            close(got.layer_heat_residual_w, 0.0, 1e-12);
            close(got.endpoint_derivative_w_k[0], conductance, 1e-13);
            close(got.endpoint_derivative_w_k[1], -conductance, 1e-13);
            close(
                got.interface_temperature_k[1],
                temperatures[1] + expected / area * distances[1] / 10.0,
                1e-12,
            );

            let reversed = TwoLayerConduction::new(
                [&film, &insulation],
                [distances[1], distances[0]],
                area,
                contact,
            )
            .unwrap()
            .evaluate([temperatures[1], temperatures[0]], config(), || Ok(()))
            .unwrap();
            assert_eq!(got.power_a_to_b_w, -reversed.power_a_to_b_w);
            assert_eq!(
                got.interface_temperature_k,
                [
                    reversed.interface_temperature_k[1],
                    reversed.interface_temperature_k[0]
                ]
            );

            let doubled =
                TwoLayerConduction::new([&insulation, &film], distances, 2.0 * area, contact)
                    .unwrap()
                    .evaluate(temperatures, config(), || Ok(()))
                    .unwrap();
            close(doubled.power_a_to_b_w, 2.0 * expected, 2e-11);
            close(doubled.layer_heat_residual_w, 0.0, 2e-12);
        }
    }
}

#[test]
fn thin_nonlinear_film_retains_heat_across_the_original_insulation_knots() {
    let insulation = ConductivityCurve::new(vec![
        [200.0, 0.01],
        [400.0, 0.05],
        [600.0, 0.03],
        [900.0, 0.09],
    ])
    .unwrap();
    let film = ConductivityCurve::new(vec![[200.0, 10.0], [400.0, 20.0], [900.0, 30.0]]).unwrap();
    let area = 0.004;
    let model = TwoLayerConduction::new([&insulation, &film], [0.01, 1e-6], area, 0.0).unwrap();

    // The three original insulation segments integrate to 25.75 W/m from
    // 300 to 850 K. With face = 300 + delta, continuity is exactly
    // (250.0001)*delta^2 + (150000.03)*delta - 25.75 = 0.
    // This independent relative-coordinate quadratic avoids subtracting
    // two large absolute temperatures in the thin nonlinear film oracle.
    let linear = 150_000.03_f64;
    let quadratic = 250.0001;
    let drop = 2.0 * 25.75 / (linear + (linear * linear + 4.0 * quadratic * 25.75).sqrt());
    let expected = area / 0.01 * (25.75 - (0.03 + 0.0001 * drop) * drop);
    let got = model.evaluate([850.0, 300.0], config(), || Ok(())).unwrap();
    close(got.power_a_to_b_w, expected, 1e-10);
    close(got.layer_heat_residual_w, 0.0, 1e-11);
    for face in got.interface_temperature_k {
        close(face, 300.0 + drop, 1e-10);
    }
    let hot_face_k = 0.03 + 0.0002 * drop;
    let cold_face_k = 15.0 + 0.05 * drop;
    let differential_resistance = 0.01 / hot_face_k + 1e-6 / cold_face_k;
    close(
        got.endpoint_derivative_w_k[0],
        area / differential_resistance * 0.08 / hot_face_k,
        1e-12,
    );
    close(
        got.endpoint_derivative_w_k[1],
        -area / differential_resistance * 15.0 / cold_face_k,
        1e-12,
    );

    let reversed = TwoLayerConduction::new([&film, &insulation], [1e-6, 0.01], area, 0.0)
        .unwrap()
        .evaluate([300.0, 850.0], config(), || Ok(()))
        .unwrap();
    assert_eq!(got.power_a_to_b_w, -reversed.power_a_to_b_w);
}

#[test]
fn retained_inverse_keeps_original_knots_endpoints_and_source_refusals() {
    let curve =
        ConductivityCurve::new(vec![[200.0, 1.0], [400.0, 5.0], [600.0, 3.0], [900.0, 9.0]])
            .unwrap();
    close(
        curve.temperature_after_integral_w_m(250.0, 2700.0).unwrap(),
        850.0,
        1e-12,
    );
    close(
        curve
            .temperature_after_integral_w_m(850.0, -2700.0)
            .unwrap(),
        250.0,
        1e-12,
    );
    assert_eq!(
        curve.temperature_after_integral_w_m(200.0, 600.0),
        Ok(400.0)
    );
    assert_eq!(
        curve.temperature_after_integral_w_m(900.0, -1800.0),
        Ok(600.0)
    );
    assert_eq!(
        curve.temperature_after_integral_w_m(200.0, 3200.0),
        Ok(900.0)
    );
    assert_eq!(
        curve.temperature_after_integral_w_m(900.0, -3200.0),
        Ok(200.0)
    );
    assert_eq!(
        curve.temperature_after_integral_w_m(200.0, 3200.0_f64.next_up()),
        Err(DepositionError::MaterialDomain)
    );
    assert_eq!(
        curve.temperature_after_integral_w_m(900.0, -3200.0_f64.next_up()),
        Err(DepositionError::MaterialDomain)
    );
    // A nonzero drop that cannot change the stored absolute temperature
    // remains refused; this change does not enlarge the representable domain.
    assert_eq!(
        curve.temperature_after_integral_w_m(400.0, 1e-16),
        Err(DepositionError::Arithmetic)
    );
}

#[test]
fn three_layer_stack_conserves_heat_through_an_internal_metal_film() {
    let curves = [constant(0.02), constant(100.0), constant(0.04)];
    let lengths = [0.01, 1e-6, 0.02];
    let area = 0.004;
    for contacts in [[0.0, 0.0], [0.02, 0.03]] {
        let resistance = 0.01 / 0.02 + 1e-6 / 100.0 + 0.02 / 0.04 + contacts[0] + contacts[1];
        let stack = LayerStackConduction::new(&curves, &lengths, &contacts, area).unwrap();
        for temperatures in [[600.0, 300.0], [300.0, 600.0]] {
            let got = stack.evaluate(temperatures, config(), || Ok(())).unwrap();
            close(
                got.power_a_to_b_w,
                area * (temperatures[0] - temperatures[1]) / resistance,
                1e-10,
            );
            close(got.max_layer_heat_residual_w, 0.0, 1e-10);
            close(got.endpoint_derivative_w_k[0], area / resistance, 1e-13);
            close(got.endpoint_derivative_w_k[1], -area / resistance, 1e-13);
            let film = got.layer_temperatures_k()[1];
            assert!(
                film[0] != film[1],
                "the thin material drop remains representable"
            );

            let reversed_curves = [curves[2].clone(), curves[1].clone(), curves[0].clone()];
            let reversed = LayerStackConduction::new(
                &reversed_curves,
                &[lengths[2], lengths[1], lengths[0]],
                &[contacts[1], contacts[0]],
                area,
            )
            .unwrap()
            .evaluate([temperatures[1], temperatures[0]], config(), || Ok(()))
            .unwrap();
            assert_eq!(got.power_a_to_b_w, -reversed.power_a_to_b_w);
        }
    }
}

#[test]
fn sub_resolution_convective_face_drop_returns_the_cell_face_only_when_its_power_closes() {
    // One ulp from a 300 K bath, a 0.05 W/K film carries about 2.8e-15 W.
    // Through this 2e5 W/K layer (k 200, 0.1 m2, 0.1 mm) the face would move
    // by about 1.4e-20 K, below half an ulp of 300 K, so the inverse refuses
    // that drop. The cell face, the power-zero end of the one-bath solve's
    // bracket, closes within the declared 1e-12 W and is returned with its
    // actual residual and no interior trial.
    let curve = constant(200.0);
    let controls = LayerConductionConfig {
        absolute_heat_tolerance_w: 1e-12,
        ..config()
    };
    let layer = ReservoirLayerGeometry::planar(0.1, 0.0001).unwrap();
    let bulk: f64 = 300.0;
    for cell in [bulk.next_up(), bulk.next_down()] {
        let film_power = 0.05 * (cell - bulk);
        assert_eq!(
            curve.temperature_after_integral_w_m(cell, -film_power / 1000.0),
            Err(DepositionError::Arithmetic)
        );
        let response = layer
            .evaluate_convective_material_to_bath(
                &curve,
                [cell, bulk],
                controls,
                |surface| Ok::<_, ()>((0.05 * (surface - bulk), 0.05, surface)),
                || Ok(()),
            )
            .unwrap();
        assert_eq!(response.surface_temperature_k, cell);
        assert_eq!(response.film, cell);
        assert_eq!(response.solid_power_to_bath_w, 0.0);
        assert_eq!(response.power_to_bath_w, film_power);
        assert_eq!(response.heat_residual_w, -film_power);
        assert_eq!(response.conductance_w_k, film_power / (cell - bulk));
        assert_eq!(response.iterations, 0);
    }
    // A 1e12 m multiplier cannot represent a 3 mK cell's face drop either,
    // but omitting it would leave 1.5e-4 W unclosed: that stays refused.
    let unclosed = ReservoirLayerGeometry::planar(1.0, 1e-12)
        .unwrap()
        .evaluate_convective_material_to_bath(
            &curve,
            [300.003, bulk],
            controls,
            |surface| Ok::<_, ()>((0.05 * (surface - bulk), 0.05, ())),
            || Ok(()),
        );
    assert!(matches!(
        unclosed,
        Err(ConvectiveLayerError::Thermal(DepositionError::Arithmetic))
    ));
}
