use super::*;
use fs_qty::Pressure;

fn surface(dims: Dims, base: f64, thermal: f64, pressure: f64, reversed: bool) -> PropertyValue {
    let mut axes = ["temperature".into(), "pressure".into()];
    let mut axis_dims = [Temperature::DIMS, Pressure::DIMS];
    let samples = [
        [256.0, 65536.0],
        [512.0, 65536.0],
        [256.0, 131072.0],
        [512.0, 131072.0],
    ]
    .into_iter()
    .map(|[t, p]| {
        let value = base + thermal * (t - 256.0) + pressure * (p - 65536.0);
        (if reversed { [p, t] } else { [t, p] }, value)
    })
    .collect();
    if reversed {
        axes.swap(0, 1);
        axis_dims.swap(0, 1);
    }
    PropertyValue::Surface2D {
        axes,
        axis_dims,
        samples,
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims,
    }
}

fn transport(reversed: bool) -> (PropertyValue, PropertyValue) {
    (
        surface(
            DynViscosity::DIMS,
            2.0_f64.powi(-16),
            2.0_f64.powi(-25),
            2.0_f64.powi(-36),
            reversed,
        ),
        surface(
            THERMAL_CONDUCTIVITY_DIMS,
            2.0_f64.powi(-6),
            2.0_f64.powi(-15),
            2.0_f64.powi(-24),
            reversed,
        ),
    )
}

#[test]
fn prepared_surface_tangents_retain_original_axis_order_pressure_edges_and_inward_endpoints() {
    for reversed in [false, true] {
        let (mu, k) = transport(reversed);
        let model = gas(
            vec![source("A", mu, k, InterpolationPolicy::LinearInside, None)],
            vec![1.0],
        );
        for pressure in [65536.0, 98304.0, 131072.0] {
            let law =
                IsobaricGasLaw::try_from_model(model.clone(), point(320.0, pressure)).unwrap();
            let original = law.reference().clone();
            let prepared = law.prepare_transport_temperature_tangent().unwrap();
            for t in [256.0, 320.0, 384.0, 512.0] {
                let tangent = prepared.evaluate(Temperature::new(t)).unwrap();
                assert_eq!(
                    tangent,
                    law.transport_temperature_tangent(Temperature::new(t))
                        .unwrap()
                );
                assert_eq!(
                    tangent.gas,
                    law.at_temperature(Temperature::new(t)).unwrap()
                );
                assert_eq!(tangent.gas.state().pressure, pressure);
                near(
                    tangent.dynamic_viscosity_derivative_pa_s_k,
                    2.0_f64.powi(-25),
                    1e-14,
                );
                near(
                    tangent.conductivity_derivative_w_m_k2,
                    2.0_f64.powi(-15),
                    1e-14,
                );
                near(
                    tangent.specific_heat_derivative_j_kg_k2,
                    R * 0.001 / 0.028,
                    1e-14,
                );
                near(
                    tangent.gas.state().dynamic_viscosity,
                    2.0_f64.powi(-16)
                        + 2.0_f64.powi(-25) * (t - 256.0)
                        + 2.0_f64.powi(-36) * (pressure - 65536.0),
                    1e-14,
                );
            }
            assert!(prepared.evaluate(Temperature::new(255.0)).is_err());
            assert_eq!(law.reference(), &original);
        }
    }
}

#[test]
fn mixed_cached_surfaces_and_uncached_curves_use_original_wilke_wms_tangents() {
    let a = source(
        "A",
        curve(&[(256.0, 2.0e-5), (512.0, 3.0e-5)], DynViscosity::DIMS),
        curve(&[(256.0, 0.025), (512.0, 0.04)], THERMAL_CONDUCTIVITY_DIMS),
        InterpolationPolicy::LinearInside,
        Some([256.0, 512.0]),
    );
    let (mu, k) = transport(true);
    let b = source(
        "B",
        mu,
        k,
        InterpolationPolicy::LinearInside,
        Some([256.0, 512.0]),
    );
    let law =
        IsobaricGasLaw::try_from_model(gas(vec![a, b], vec![0.75, 0.25]), point(320.0, 98304.0))
            .unwrap();
    let prepared = law.prepare_transport_temperature_tangent().unwrap();
    for t in [300.0, 384.0, 450.0] {
        let tangent = prepared.evaluate(Temperature::new(t)).unwrap();
        let h = 0.01;
        let samples = [-2.0, -1.0, 1.0, 2.0]
            .map(|scale| values(&law.at_temperature(Temperature::new(t + scale * h)).unwrap()));
        for (i, analytic) in derivatives(&tangent).into_iter().enumerate() {
            let difference = (samples[0][i] - 8.0 * samples[1][i] + 8.0 * samples[2][i]
                - samples[3][i])
                / (12.0 * h);
            near(analytic, difference, 3e-8);
        }
        assert_eq!(
            tangent,
            law.transport_temperature_tangent(Temperature::new(t))
                .unwrap()
        );
    }
}

#[test]
fn surface_corners_refuse_except_on_the_original_claims_admitted_side() {
    let bind = |bounds| {
        let (mut mu, k) = transport(false);
        let PropertyValue::Surface2D { samples, .. } = &mut mu else {
            unreachable!()
        };
        samples[3].1 += 2.0_f64.powi(-18);
        IsobaricGasLaw::try_from_model(
            gas(
                vec![source(
                    "A",
                    mu,
                    k,
                    InterpolationPolicy::LinearInside,
                    Some(bounds),
                )],
                vec![1.0],
            ),
            point(384.0, 98304.0),
        )
        .unwrap()
    };
    let complete = bind([256.0, 512.0]);
    let original = complete.reference().clone();
    assert!(complete.at_temperature(Temperature::new(384.0)).is_ok());
    let prepared = complete.prepare_transport_temperature_tangent().unwrap();
    assert!(prepared.evaluate(Temperature::new(384.0)).is_err());
    assert!(
        complete
            .transport_temperature_tangent(Temperature::new(384.0))
            .is_err()
    );
    assert_eq!(complete.reference(), &original);
    let left = bind([256.0, 384.0])
        .transport_temperature_tangent(Temperature::new(384.0))
        .unwrap();
    let right = bind([384.0, 512.0])
        .transport_temperature_tangent(Temperature::new(384.0))
        .unwrap();
    near(
        left.dynamic_viscosity_derivative_pa_s_k,
        2.0_f64.powi(-25),
        1e-14,
    );
    near(
        right.dynamic_viscosity_derivative_pa_s_k,
        2.0_f64.powi(-25) + 2.0_f64.powi(-26),
        1e-14,
    );
    assert_eq!(left.gas.state(), right.gas.state());
}

#[test]
fn tabulated_and_isolated_surface_samples_retain_values_without_gaining_tangents() {
    for isolated in [false, true] {
        let (mut mu, mut k) = transport(false);
        if isolated {
            for value in [&mut mu, &mut k] {
                let PropertyValue::Surface2D { triangles, .. } = value else {
                    unreachable!()
                };
                triangles.pop();
            }
        }
        let interpolation = if isolated {
            InterpolationPolicy::LinearInside
        } else {
            InterpolationPolicy::TabulatedOnly
        };
        let law = IsobaricGasLaw::try_from_model(
            gas(
                vec![source("A", mu, k, interpolation, Some([256.0, 512.0]))],
                vec![1.0],
            ),
            point(512.0, 131072.0),
        )
        .unwrap();
        law.at_temperature(Temperature::new(512.0)).unwrap();
        assert!(
            law.transport_temperature_tangent(Temperature::new(512.0))
                .is_err()
        );
        if isolated {
            assert!(
                law.prepare_transport_temperature_tangent()
                    .unwrap()
                    .evaluate(Temperature::new(512.0))
                    .is_err()
            );
        } else {
            assert!(law.prepare_transport_temperature_tangent().is_err());
        }
    }
}
