//! G1/G3/G4: source data follows the actual solidus/liquidus caloric path.
use super::*;

fn source(scale: f64, reverse: bool, triangles: Vec<[u32; 3]>) -> MaterialCard {
    let mut value = PropertyValue::Surface2D {
        axes: ["T".into(), "f".into()],
        axis_dims: [Temperature::DIMS, Dims::NONE],
        samples: vec![
            ([300.0, 0.0], 2.0 * scale),
            ([800.0, 0.0], 7.0 * scale),
            ([800.0, 1.0], 12.0 * scale),
            ([300.0, 1.0], 5.0 * scale),
        ],
        triangles,
        dims: KD,
    };
    if reverse {
        if let PropertyValue::Surface2D {
            axes,
            axis_dims,
            samples,
            ..
        } = &mut value
        {
            axes.swap(0, 1);
            axis_dims.swap(0, 1);
            for (point, _) in samples {
                point.swap(0, 1);
            }
        }
    }
    card(vec![claim(value, [0.0, 1.0])])
}
fn chart(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    curve(
        card.content_hash(),
        &[
            (-50.0, 300.0, 0.0),
            (0.0, 350.0, 0.0),
            (100.0, 750.0, 1.0),
            (150.0, 800.0, 1.0),
        ],
    )
}

#[test]
fn solidus_liquidus_data_keeps_both_tangents_and_every_crossed_cell() {
    for scale in [0.5, 2.0] {
        for reverse in [false, true] {
            let card = source(scale, reverse, vec![[0, 1, 2], [0, 2, 3]]);
            let curve = chart(&card);
            let path = prepare(&card, &curve, [-50.0, 150.0], || Ok(())).unwrap();
            let knots: Vec<_> = path.knots().iter().map(|k| k[0]).collect();
            assert_eq!(knots.len(), 5);
            for (actual, expected) in knots.into_iter().zip([-50.0, 0.0, 50.0, 100.0, 150.0]) {
                close(actual, expected);
            }
            // Below the diagonal k=2+.01(T-300)+5f; above it
            // k=2+.014(T-300)+3f. Both source coordinates evolve.
            for (h, expected, dk, dt) in [
                (-25.0, 2.25, 0.01, 1.0),
                (0.0, 2.5, 0.09, 4.0),
                (25.0, 4.75, 0.09, 4.0),
                (75.0, 9.15, 0.086, 4.0),
                (125.0, 11.65, 0.014, 1.0),
                (150.0, 12.0, 0.014, 1.0),
            ] {
                let sample = path.state_at_specific_enthalpy(h).unwrap();
                close(sample.conductivity_w_m_k(), scale * expected);
                close(
                    sample.conductivity_enthalpy_derivative_kg_w_m_k_j(),
                    scale * dk,
                );
                close(sample.temperature_enthalpy_derivative_kg_k_j(), dt);
                let native = ImplicitConductionLaw::state_at_specific_enthalpy(&path, h).unwrap();
                close(native.conductivity_slope_kg_w_m_k_j, scale * dk);
                let original = curve.values_at_specific_enthalpy(h).unwrap();
                assert_eq!(sample.temperature_k(), original.temperature_k);
                assert_eq!(sample.liquid_mass_fraction(), original.liquid_mass_fraction);
            }
            for source in path.source_paths() {
                for (_, answer) in &source.knots {
                    assert_eq!(answer.receipt.selected, path.selected_claim());
                    card.claims().verify_receipt(&answer.receipt).unwrap();
                }
            }
        }
    }
}

#[test]
fn mixed_sensible_curve_and_constant_data_need_no_phase_mixing_rule() {
    for value in [
        PropertyValue::Scalar {
            value: 3.0,
            dims: KD,
        },
        phase_values(&[(0.0, 2.0), (0.25, 4.0), (1.0, 7.0)]),
    ] {
        let is_scalar = matches!(&value, PropertyValue::Scalar { .. });
        let mut response = claim(value, [0.0, 1.0]);
        if is_scalar {
            response.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        let card = card(vec![response]);
        let curve = curve(
            card.content_hash(),
            &[
                (0.0, 300.0, 0.0),
                (50.0, 400.0, 0.5),
                (100.0, 500.0, 0.5),
                (150.0, 800.0, 1.0),
            ],
        );
        let path = prepare(&card, &curve, [0.0, 150.0], || Ok(())).unwrap();
        let sample = path.state_at_specific_enthalpy(75.0).unwrap();
        close(sample.temperature_k(), 450.0);
        close(sample.liquid_mass_fraction(), 0.5);
        close(
            sample.conductivity_w_m_k(),
            if is_scalar { 3.0 } else { 5.0 },
        );
        assert_eq!(sample.conductivity_enthalpy_derivative_kg_w_m_k_j(), 0.0);
        close(sample.temperature_enthalpy_derivative_kg_k_j(), 2.0);
        if !is_scalar {
            assert!(path.knots().iter().any(|k| k[0] == 25.0));
        }
    }
}

#[test]
fn missing_mushy_support_and_cancelled_preparation_do_not_change_the_original_chart() {
    // Initial solid state belongs to cell zero; the actual alloy path exits it.
    let hole = source(1.0, false, vec![[0, 1, 2]]);
    let curve = chart(&hole);
    assert!(matches!(
        prepare(&hole, &curve, [0.0, 100.0], || Ok(())),
        Err(PhaseConductivityError::Source(_))
    ));
    let card = source(1.0, false, vec![[0, 1, 2], [0, 2, 3]]);
    let curve = chart(&card);
    let original = curve.identity();
    let mut calls = 0;
    let path = prepare(&card, &curve, [-50.0, 150.0], || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    let mut retry_calls = 0;
    assert!(matches!(
        prepare(&card, &curve, [-50.0, 150.0], || {
            retry_calls += 1;
            if retry_calls == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(PhaseConductivityError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(curve.identity(), original);
    let retry = prepare(&card, &curve, [-50.0, 150.0], || Ok(())).unwrap();
    assert_eq!(path.identity(), retry.identity());
    assert_eq!(path.knots(), retry.knots());
    assert!(path.state_at_specific_enthalpy(151.0).is_err());
}
