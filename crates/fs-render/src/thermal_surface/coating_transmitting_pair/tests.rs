use super::*;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, ConductorDataStatus,
    ConductorIorSample, ConductorOptics, ConductorSource,
};
use crate::dielectric::{
    CauchyIor, DielectricCoatingStack, DielectricIor, DielectricIorSample, TabulatedIor,
};
use crate::thermal_surface::{CoatingExchangeResponse, CoatingGapResponse};

fn source(tag: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative).unwrap()
}
fn film(bounds: (f64, f64), eta: f64, k: f64, depth: f64, tag: u8) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(
                bounds.0 + (bounds.1 - bounds.0) * i as f64 / 8.0,
                eta,
                k,
            )
            .unwrap()
        }),
        depth,
        (bounds.1 - bounds.0) / 8.0,
        source(tag),
    )
    .unwrap()
}
fn medium(bounds: (f64, f64), index: impl Fn(usize) -> f64) -> DielectricIor {
    TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(
                bounds.0 + (bounds.1 - bounds.0) * i as f64 / 8.0,
                index(i),
            )
            .unwrap()
        }),
        (bounds.1 - bounds.0) / 8.0,
    )
    .unwrap()
    .into()
}
fn emitter(
    films: &[CoherentAbsorbingCoating],
    bounds: (f64, f64),
    temperature: f64,
    reverse: bool,
    tag: u8,
) -> CoatingThermalEmission {
    CoatingThermalEmission::try_transmitting(
        &films
            .iter()
            .map(CoatingLayer::Absorbing)
            .collect::<Vec<_>>(),
        &vec![temperature; films.len()],
        medium(bounds, |_| 1.0),
        reverse,
        ContentHash([tag; 32]),
    )
    .unwrap()
}
fn budget(q: CoatingExchangeQuadrature) -> CoatingTransmittingPairBudget {
    CoatingTransmittingPairBudget {
        max_channel_spectral_nodes: CHANNELS * q.wavelength_bins,
        max_optical_evaluations: 4 * q.wavelength_bins * q.cosine_bins,
    }
}
fn response(
    first: &CoatingThermalEmission,
    second: &CoatingThermalEmission,
    area: f64,
    q: CoatingExchangeQuadrature,
) -> CoatingTransmittingPairResponse {
    CoatingTransmittingPairResponse::from_emitters(first, second, area, q, budget(q), || Ok(()))
        .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= relative * expected.abs() + 4.0 * f64::from_bits(1),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn planck(t: f64, lambda_nm: f64) -> f64 {
    let lambda = lambda_nm * 1e-9;
    let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
    2.0 * h * c * c / (lambda.powi(5) * (h * c / (lambda * k * t)).exp_m1()) * 1e-9
}
fn scalar_port(a: f64, t: f64) -> ReciprocalPortPower {
    ReciprocalPortPower {
        reflected: [1.0 - (a + t); 2],
        transmitted: [t; 2],
        absorbed: [a; 2],
    }
}

// Independent two-unknown intensity solve and direct boundary powers, without
// the reduced six-channel coefficients. Intended for ordinary, well-resolved
// R/T/A; tiny escapes have a separate analytic scaled reference below.
fn direct_intensity(ports: [ReciprocalPortPower; 4], b: [f64; 4], pol: usize) -> [f64; 4] {
    let [g1, o1, g2, o2] = ports;
    let [b1, b2, bl, br] = b;
    let input_plus = g1.transmitted[pol] * bl + g1.absorbed[pol] * b1;
    let input_minus = g2.transmitted[pol] * br + g2.absorbed[pol] * b2;
    let plus = (input_plus + g1.reflected[pol] * input_minus)
        / (1.0 - g1.reflected[pol] * g2.reflected[pol]);
    let minus = input_minus + g2.reflected[pol] * plus;
    let left_out = o1.reflected[pol] * bl + g1.transmitted[pol] * minus + o1.absorbed[pol] * b1;
    let right_out = o2.reflected[pol] * br + g2.transmitted[pol] * plus + o2.absorbed[pol] * b2;
    [
        o1.absorbed[pol] * (b1 - bl) + g1.absorbed[pol] * (b1 - minus),
        o2.absorbed[pol] * (b2 - br) + g2.absorbed[pol] * (b2 - plus),
        bl - left_out,
        br - right_out,
    ]
}
fn original_port(
    stack: &DielectricCoatingStack<'_>,
    lambda: f64,
    mu: f64,
    reverse: bool,
) -> ReciprocalPortPower {
    let work = stack
        .layer_power_balance(lambda, reverse, mu, 1.0, 1.0)
        .unwrap();
    let power = |matrix: crate::polarization::Mueller| {
        [
            matrix.0[0][0] + matrix.0[0][1],
            matrix.0[0][0] - matrix.0[0][1],
        ]
    };
    ReciprocalPortPower {
        reflected: power(work.fresnel.reflection),
        transmitted: power(work.fresnel.transmission),
        absorbed: work.absorbed,
    }
}

#[test]
fn transmitting_pair_positive_network_matches_intensity_balance_and_tiny_escape() {
    let ports = [
        scalar_port(0.2, 0.3),
        scalar_port(0.4, 0.3),
        scalar_port(0.3, 0.15),
        scalar_port(0.6, 0.15),
    ];
    let c = exchange_fractions(ports).unwrap().map(|value| value[0]);
    assert!(c.into_iter().all(|value| value > 0.0));
    for b in [[2.0, 0.5, 3.0, 0.2], [0.0, 0.0, 1.0, 0.0], [1.0; 4]] {
        let [b1, b2, bl, br] = b;
        let edges = [
            c[0] * (b1 - b2),
            c[1] * (b1 - bl),
            c[2] * (b1 - br),
            c[3] * (b2 - bl),
            c[4] * (b2 - br),
            c[5] * (bl - br),
        ];
        let reduced = [
            edges[0] + edges[1] + edges[2],
            -edges[0] + edges[3] + edges[4],
            -edges[1] - edges[3] + edges[5],
            -edges[2] - edges[4] - edges[5],
        ];
        let direct = direct_intensity(ports, b, 0);
        for i in 0..4 {
            assert!((reduced[i] - direct[i]).abs() < 2e-15);
        }
        assert!(reduced.into_iter().sum::<f64>().abs() < 5e-16);
    }
    let tiny = 1e-200;
    assert_eq!(tiny * tiny, 0.0);
    let weak = exchange_fractions([
        scalar_port(tiny, 2.0 * tiny),
        scalar_port(5.0 * tiny, 2.0 * tiny),
        scalar_port(3.0 * tiny, 4.0 * tiny),
        scalar_port(6.0 * tiny, 4.0 * tiny),
    ])
    .unwrap();
    for (actual, expected) in weak.into_iter().zip([0.3, 5.2, 0.4, 0.6, 7.2, 0.8]) {
        close(actual[0], expected * tiny, 4e-15);
    }
    let trapped = exchange_fractions([
        scalar_port(0.0, 0.0),
        scalar_port(0.2, 0.0),
        scalar_port(0.0, 0.0),
        scalar_port(0.4, 0.0),
    ])
    .unwrap();
    assert_eq!(
        trapped.map(|value| value[0]),
        [0.0, 0.2, 0.0, 0.0, 0.4, 0.0]
    );
    let mut invalid = ports;
    invalid[0].absorbed[0] = -1.0;
    assert_eq!(
        exchange_fractions(invalid),
        Err(ThermalSurfaceError::Unrepresentable)
    );
    invalid = ports;
    invalid[1] = scalar_port(0.4, 0.29);
    assert_eq!(
        exchange_fractions(invalid),
        Err(ThermalSurfaceError::Unrepresentable)
    );
}

#[test]
fn transmitting_pair_coherent_asymmetric_infrared_fields_close_all_four_terminals() {
    let bounds = (4000.0, 20_000.0);
    let first_films = [
        film(bounds, 2.7, 0.4, 0.42e-6, 1),
        film(bounds, 1.3, 0.02, 1.6e-6, 2),
    ];
    let second_films = [
        film(bounds, 1.8, 0.09, 1.0e-6, 3),
        film(bounds, 3.5, 0.6, 0.25e-6, 4),
    ];
    let first = emitter(&first_films, bounds, 390.0, false, 5);
    let second = emitter(&second_films, bounds, 280.0, true, 6);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 12,
        cosine_bins: 14,
    };
    let area = 0.025;
    let pair = response(&first, &second, area, q);
    let temperatures = [390.0, 280.0, 510.0, 220.0];
    let rates = pair
        .evaluate(390.0, 280.0, Some(510.0), Some(220.0), || Ok(()))
        .unwrap();
    let first_layers = first_films
        .iter()
        .map(CoatingLayer::Absorbing)
        .collect::<Vec<_>>();
    let second_layers = second_films
        .iter()
        .map(CoatingLayer::Absorbing)
        .collect::<Vec<_>>();
    let stacks = [
        DielectricCoatingStack::try_mixed_in_band(&first_layers, bounds).unwrap(),
        DielectricCoatingStack::try_mixed_in_band(&second_layers, bounds).unwrap(),
    ];
    let mut expected = [0.0; 4];
    let mut asymmetry = 0.0_f64;
    let width = (bounds.1 - bounds.0) / q.wavelength_bins as f64;
    for w in 0..q.wavelength_bins {
        let lambda = bounds.0 + (w as f64 + 0.5) * width;
        let b = temperatures.map(|t| planck(t, lambda));
        for angle in 0..q.cosine_bins {
            let mu = (angle as f64 + 0.5) / q.cosine_bins as f64;
            let ports = [
                original_port(&stacks[0], lambda, mu, false),
                original_port(&stacks[0], lambda, mu, true),
                original_port(&stacks[1], lambda, mu, true),
                original_port(&stacks[1], lambda, mu, false),
            ];
            for pol in 0..2 {
                asymmetry = asymmetry.max((ports[0].absorbed[pol] - ports[1].absorbed[pol]).abs());
                let direct = direct_intensity(ports, b, pol);
                for i in 0..4 {
                    expected[i] += area * PI * mu * width / q.cosine_bins as f64 * direct[i];
                }
            }
        }
    }
    assert!(asymmetry > 0.05);
    let body = rates.body_outward_power_w();
    let baths = rates.bath_outward_power_w();
    for (actual, reference) in [body[0], body[1], baths[0], baths[1]]
        .into_iter()
        .zip(expected)
    {
        close(actual, reference, 3e-12);
    }
    assert!((body.into_iter().sum::<f64>() + baths.into_iter().sum::<f64>()).abs() < 2e-14);
    assert_eq!(pair.optical_evaluations(), 4 * 12 * 14);
    assert_eq!(pair.channel_spectral_nodes(), 6 * 12);
    assert_eq!(
        pair.optical_identities(),
        [first.optical_identity(), second.optical_identity()]
    );
    assert_eq!(pair.wavelength_bounds_nm(), bounds);
    assert_eq!(pair.area_m2(), area);
    assert_eq!(pair.quadrature(), q);

    let equilibrium = pair
        .evaluate(300.0, 300.0, Some(300.0), Some(300.0), || Ok(()))
        .unwrap();
    for rate in [equilibrium.pair, equilibrium.bath_transfer]
        .into_iter()
        .chain(equilibrium.first_baths)
        .chain(equilibrium.second_baths)
    {
        assert_eq!(rate.outward_power_w, 0.0);
        assert!(rate.secant_w_k > 0.0);
    }
    let mirrored = response(&second, &first, area, q)
        .evaluate(280.0, 390.0, Some(220.0), Some(510.0), || Ok(()))
        .unwrap();
    assert_eq!(mirrored.pair.outward_power_w, -rates.pair.outward_power_w);
    assert_eq!(mirrored.pair.secant_w_k, rates.pair.secant_w_k);
    assert_eq!(
        mirrored.bath_transfer.outward_power_w,
        -rates.bath_transfer.outward_power_w
    );
    for i in 0..2 {
        assert_eq!(mirrored.first_baths[i], rates.second_baths[1 - i]);
        assert_eq!(mirrored.second_baths[i], rates.first_baths[1 - i]);
    }
    let warmer_preparation = emitter(&first_films, bounds, 900.0, false, 5);
    assert_eq!(
        pair.identity(),
        response(&warmer_preparation, &second, area, q).identity()
    );
}

#[test]
fn transmitting_pair_transparent_limit_passes_baths_without_heating_bodies() {
    let bounds = (4000.0, 80_000.0);
    let first = emitter(&[], bounds, 100.0, false, 10);
    let second = emitter(&[], (2000.0, 90_000.0), 900.0, true, 11);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 1024,
        cosine_bins: 1,
    };
    let pair = response(&first, &second, 0.01, q);
    assert_eq!(pair.wavelength_bounds_nm(), bounds);
    let rates = pair
        .evaluate(100.0, 900.0, Some(310.0), Some(280.0), || Ok(()))
        .unwrap();
    assert_eq!(rates.body_outward_power_w(), [0.0; 2]);
    assert_eq!(rates.pair.outward_power_w, 0.0);
    let n = 16_384;
    let step = (bounds.1 - bounds.0) / n as f64;
    let value = |lambda| planck(310.0, lambda) - planck(280.0, lambda);
    let mut reference = value(bounds.0) + value(bounds.1);
    for i in 1..n {
        reference += (if i % 2 == 0 { 2.0 } else { 4.0 }) * value(bounds.0 + i as f64 * step);
    }
    reference *= 0.01 * PI * step / 3.0;
    close(rates.bath_transfer.outward_power_w, reference, 2e-5);
    close(rates.bath_transfer.secant_w_k, reference / 30.0, 2e-5);
    let left = pair
        .evaluate(100.0, 900.0, Some(310.0), None, || Ok(()))
        .unwrap();
    let right = pair
        .evaluate(100.0, 900.0, None, Some(310.0), || Ok(()))
        .unwrap();
    assert!(left.bath_transfer.outward_power_w > 0.0);
    assert_eq!(
        left.bath_transfer.outward_power_w,
        -right.bath_transfer.outward_power_w
    );
    assert_eq!(
        left.bath_transfer.secant_w_k,
        right.bath_transfer.secant_w_k
    );
    let dark = pair.evaluate(100.0, 900.0, None, None, || Ok(())).unwrap();
    assert_eq!(
        dark.bath_transfer,
        CoatingExchangeRates {
            outward_power_w: 0.0,
            secant_w_k: 0.0
        }
    );
}

#[test]
fn transmitting_pair_opaque_limit_recovers_original_gap_and_exterior_owners() {
    let bounds = (4000.0, 20_000.0);
    // Stay within the original 100-micrometre film admission. These extinction
    // coefficients make transmitted power round to zero throughout this band.
    let films = [
        film(bounds, 0.8, 16.0, 1e-4, 20),
        film(bounds, 1.2, 18.0, 1e-4, 21),
    ];
    let first = emitter(&films[..1], bounds, 400.0, false, 22);
    let second = emitter(&films[1..], bounds, 300.0, true, 23);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 16,
        cosine_bins: 12,
    };
    let pair = response(&first, &second, 0.02, q);
    let rates = pair
        .evaluate(400.0, 300.0, Some(280.0), Some(310.0), || Ok(()))
        .unwrap();
    let opaque = |eta, k, t, tag| {
        let metal = ConductorOptics::try_new_in_band(
            core::array::from_fn(|i| {
                ConductorIorSample::try_new_in_band(
                    bounds.0 + (bounds.1 - bounds.0) * i as f64 / 8.0,
                    eta,
                    k,
                )
                .unwrap()
            }),
            (bounds.1 - bounds.0) / 8.0,
            source(tag),
        )
        .unwrap();
        CoatingThermalEmission::try_opaque(&[], &[], metal, t, ContentHash([tag; 32])).unwrap()
    };
    let conventional_budget = CoatingExchangeBudget {
        max_spectral_nodes: q.wavelength_bins,
        max_optical_evaluations: 2 * q.wavelength_bins * q.cosine_bins,
    };
    let gap = CoatingGapResponse::from_emitters(
        &opaque(0.8, 16.0, 400.0, 24),
        &opaque(1.2, 18.0, 300.0, 25),
        0.02,
        q,
        conventional_budget,
        || Ok(()),
    )
    .unwrap()
    .evaluate(400.0, 300.0, || Ok(()))
    .unwrap();
    close(rates.pair.outward_power_w, gap.outward_power_w, 3e-12);
    let left =
        CoatingExchangeResponse::from_emitter(&first, 0.02, q, conventional_budget, || Ok(()))
            .unwrap()
            .evaluate(400.0, Some(280.0), || Ok(()))
            .unwrap();
    let right =
        CoatingExchangeResponse::from_emitter(&second, 0.02, q, conventional_budget, || Ok(()))
            .unwrap()
            .evaluate(300.0, Some(310.0), || Ok(()))
            .unwrap();
    close(
        rates.first_baths[0].outward_power_w,
        left.outward_power_w,
        3e-12,
    );
    close(
        rates.second_baths[1].outward_power_w,
        right.outward_power_w,
        3e-12,
    );
    assert_eq!(rates.first_baths[1].outward_power_w, 0.0);
    assert_eq!(rates.second_baths[0].outward_power_w, 0.0);
    assert_eq!(rates.bath_transfer.outward_power_w, 0.0);
}

#[test]
fn transmitting_pair_retains_tiny_physical_work_and_planck_secants() {
    let bounds = (380.0, 780.0);
    let films = [film(bounds, 1e-120, 100.0, 5e-6, 30)];
    let body = emitter(&films, bounds, 300.0, false, 31);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 1,
        cosine_bins: 1,
    };
    let pair = response(&body, &body, 1e-202, q);
    let rate = pair.evaluate(1e300, 300.0, None, None, || Ok(())).unwrap();
    let a = body.reciprocal_absorbed_fraction(580.0, 0.5).unwrap();
    assert!(
        a.into_iter()
            .all(|fraction| fraction > 0.0 && 1.0 - fraction == 1.0)
    );
    assert_eq!(1e-202 * (a[0] + a[1]), 0.0);
    let slope = 2.0 * 299_792_458.0 * 1.380_649e-23 * 1e27 / 580.0_f64.powi(4);
    let expected = (slope * 1e300) * (a[0] + a[1]) * (0.5 * 1e-202) * (PI * 0.5 * 400.0);
    close(rate.pair.outward_power_w, expected, 3e-12);
    assert!(rate.pair.outward_power_w > 0.0 && rate.pair.secant_w_k > 0.0);
    close(rate.pair.secant_w_k, rate.pair.outward_power_w / 1e300, 0.0);
}

#[test]
fn transmitting_pair_admission_cancellation_and_dark_temperature_validation() {
    let bounds = (4000.0, 20_000.0);
    let body = emitter(&[], bounds, 300.0, false, 40);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 2,
        cosine_bins: 3,
    };
    for b in [
        CoatingTransmittingPairBudget {
            max_channel_spectral_nodes: 11,
            ..budget(q)
        },
        CoatingTransmittingPairBudget {
            max_optical_evaluations: 23,
            ..budget(q)
        },
    ] {
        let mut calls = 0;
        let result =
            CoatingTransmittingPairResponse::from_emitters(&body, &body, 0.01, q, b, || {
                calls += 1;
                Ok(())
            });
        assert!(matches!(result, Err(CoatingExchangeError::Budget)));
        assert_eq!(calls, 0);
    }
    let other_medium = |ior| {
        CoatingThermalEmission::try_transmitting(&[], &[], ior, false, ContentHash([41; 32]))
            .unwrap()
    };
    for other in [
        other_medium(medium(bounds, |_| 1.5)),
        // A wavelength sample at n=1 must not certify a dispersive full band.
        other_medium(medium(bounds, |i| if i % 2 == 0 { 1.0 } else { 1.01 })),
    ] {
        assert!(matches!(
            CoatingTransmittingPairResponse::from_emitters(
                &body,
                &other,
                0.01,
                q,
                budget(q),
                || Ok(())
            ),
            Err(CoatingExchangeError::VacuumTransmittingRequired)
        ));
    }
    let disjoint = emitter(&[], (30_000.0, 40_000.0), 300.0, false, 42);
    assert!(matches!(
        CoatingTransmittingPairResponse::from_emitters(
            &body,
            &disjoint,
            0.01,
            q,
            budget(q),
            || Ok(())
        ),
        Err(CoatingExchangeError::WavelengthDomain)
    ));
    let films = [
        film(bounds, 1.5, 0.1, 1e-6, 43),
        film(bounds, 2.0, 0.1, 1e-6, 44),
    ];
    let nonisothermal = CoatingThermalEmission::try_transmitting(
        &[
            CoatingLayer::Absorbing(&films[0]),
            CoatingLayer::Absorbing(&films[1]),
        ],
        &[300.0, 301.0],
        medium(bounds, |_| 1.0),
        false,
        ContentHash([45; 32]),
    )
    .unwrap();
    assert!(matches!(
        CoatingTransmittingPairResponse::from_emitters(
            &body,
            &nonisothermal,
            0.01,
            q,
            budget(q),
            || Ok(())
        ),
        Err(CoatingExchangeError::Nonisothermal)
    ));
    let mut calls = 0;
    let result =
        CoatingTransmittingPairResponse::from_emitters(&body, &body, 0.01, q, budget(q), || {
            calls += 1;
            if calls == 7 {
                Err(CoatingExchangeError::Cancelled)
            } else {
                Ok(())
            }
        });
    assert!(matches!(result, Err(CoatingExchangeError::Cancelled)));
    assert_eq!(calls, 7);
    let pair = response(&body, &body, 0.01, q);
    for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert_eq!(
            pair.evaluate(invalid, 300.0, None, None, || Ok(())),
            Err(CoatingExchangeError::Thermal(
                ThermalSurfaceError::Temperature
            ))
        );
        assert_eq!(
            pair.evaluate(300.0, 300.0, Some(invalid), None, || Ok(())),
            Err(CoatingExchangeError::Thermal(
                ThermalSurfaceError::Temperature
            ))
        );
    }
    assert_eq!(
        pair.evaluate(300.0, 300.0, None, None, || Err(
            CoatingExchangeError::Cancelled
        )),
        Err(CoatingExchangeError::Cancelled)
    );
    // The original constant Cauchy vacuum is also admitted, within its visible
    // support. No source-band extrapolation is granted by the analytic law.
    let visible = other_medium(CauchyIor::try_new(1.0, 0.0, 0.0).unwrap().into());
    assert!(visible.is_vacuum_transmitting());
}
