use super::*;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, ConductorDataStatus, ConductorSource,
};
use crate::dielectric::{DielectricIor, DielectricIorSample, TabulatedIor};
use crate::thermal_surface::{
    CoatingExchangeResponse, CoatingTransmittingPairBudget, CoatingTransmittingPairResponse,
};

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
        ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative)
            .unwrap(),
    )
    .unwrap()
}
fn medium(bounds: (f64, f64), eta: f64) -> DielectricIor {
    TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(
                bounds.0 + (bounds.1 - bounds.0) * i as f64 / 8.0,
                eta,
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
    reverse: bool,
    tag: u8,
) -> CoatingThermalEmission {
    CoatingThermalEmission::try_transmitting(
        &films
            .iter()
            .map(CoatingLayer::Absorbing)
            .collect::<Vec<_>>(),
        &vec![300.0; films.len()],
        medium(bounds, 1.0),
        reverse,
        ContentHash([tag; 32]),
    )
    .unwrap()
}
fn budget(n: usize, q: CoatingExchangeQuadrature) -> CoatingRadiativeStackBudget {
    CoatingRadiativeStackBudget {
        max_bodies: n,
        max_channel_spectral_nodes: channel_count(n) * q.wavelength_bins,
        max_optical_evaluations: 2 * n * q.wavelength_bins * q.cosine_bins,
    }
}
fn prepare(
    bodies: &[CoatingThermalEmission],
    area: f64,
    q: CoatingExchangeQuadrature,
) -> CoatingRadiativeStackResponse {
    CoatingRadiativeStackResponse::from_emitters(
        bodies,
        area,
        q,
        budget(bodies.len(), q),
        || Ok(()),
    )
    .unwrap()
}
fn scalar(a: f64, t: f64) -> ReciprocalPortPower {
    ReciprocalPortPower {
        reflected: [1.0 - (a + t); 2],
        transmitted: [t; 2],
        absorbed: [a; 2],
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs() + 32.0 * f64::from_bits(1),
        "actual={actual:.17e} expected={expected:.17e}"
    );
}
fn planck(t: f64, wavelength_nm: f64) -> f64 {
    let lambda = wavelength_nm * 1e-9;
    let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
    2.0 * h * c * c / (lambda.powi(5) * (h * c / (lambda * k * t)).exp_m1()) * 1e-9
}

// Independent dense solve of 2N incoming directional intensities, including
// the two authored bath boundary equations. No prefix or reduced-channel
// calculation appears here. This reference is for ordinary resolved R/T/A;
// the singular weak-loss limit has its own rescaled analytic reference below.
fn dense_balance(ports: &[[ReciprocalPortPower; 2]], b: &[f64], pol: usize) -> Vec<f64> {
    let n = ports.len();
    let count = 2 * n;
    let mut matrix = vec![vec![0.0; count + 1]; count];
    for i in 0..n {
        matrix[2 * i][2 * i] = 1.0;
        if i == 0 {
            matrix[0][count] = b[n];
        } else {
            let [left, right] = ports[i - 1];
            matrix[2 * i][2 * (i - 1)] = -left.transmitted[pol];
            matrix[2 * i][2 * (i - 1) + 1] = -right.reflected[pol];
            matrix[2 * i][count] = right.absorbed[pol] * b[i - 1];
        }
        matrix[2 * i + 1][2 * i + 1] = 1.0;
        if i == n - 1 {
            matrix[2 * i + 1][count] = b[n + 1];
        } else {
            let [left, _] = ports[i + 1];
            matrix[2 * i + 1][2 * (i + 1)] = -left.reflected[pol];
            matrix[2 * i + 1][2 * (i + 1) + 1] = -left.transmitted[pol];
            matrix[2 * i + 1][count] = left.absorbed[pol] * b[i + 1];
        }
    }
    for col in 0..count {
        let pivot = (col..count)
            .max_by(|&a, &b| matrix[a][col].abs().total_cmp(&matrix[b][col].abs()))
            .unwrap();
        matrix.swap(col, pivot);
        let pivot = matrix[col][col];
        assert!(pivot.abs() > 1e-12);
        for k in col..=count {
            matrix[col][k] /= pivot;
        }
        for row in 0..count {
            if row == col {
                continue;
            }
            let factor = matrix[row][col];
            for k in col..=count {
                matrix[row][k] -= factor * matrix[col][k];
            }
        }
    }
    let field = matrix.iter().map(|row| row[count]).collect::<Vec<_>>();
    let mut outward = vec![0.0; n + 2];
    for i in 0..n {
        outward[i] = ports[i][0].absorbed[pol] * (b[i] - field[2 * i])
            + ports[i][1].absorbed[pol] * (b[i] - field[2 * i + 1]);
    }
    outward[n] = b[n]
        - (ports[0][0].reflected[pol] * field[0]
            + ports[0][0].transmitted[pol] * field[1]
            + ports[0][0].absorbed[pol] * b[0]);
    outward[n + 1] = b[n + 1]
        - (ports[n - 1][0].transmitted[pol] * field[2 * n - 2]
            + ports[n - 1][1].reflected[pol] * field[2 * n - 1]
            + ports[n - 1][1].absorbed[pol] * b[n - 1]);
    outward
}
fn network_balance(c: &[[f64; 2]], b: &[f64], pol: usize) -> Vec<f64> {
    let n = b.len() - 2;
    let mut out = vec![0.0; n + 2];
    let mut edge = 0;
    for i in 0..n {
        for j in i + 1..n {
            let p = c[edge][pol] * (b[i] - b[j]);
            out[i] += p;
            out[j] -= p;
            edge += 1;
        }
    }
    for i in 0..n {
        for j in n..n + 2 {
            let p = c[edge][pol] * (b[i] - b[j]);
            out[i] += p;
            out[j] -= p;
            edge += 1;
        }
    }
    let p = c[edge][pol] * (b[n] - b[n + 1]);
    out[n] += p;
    out[n + 1] -= p;
    out
}

#[test]
fn g1_ordered_stack_matches_independent_dense_transport_and_retains_tiny_escape() {
    for n in [1, 2, 3, 5, 16] {
        let ports = (0..n)
            .map(|i| {
                let t = 0.1 + 0.015 * i as f64;
                [
                    scalar(0.12 + 0.02 * i as f64, t),
                    scalar(0.43 - 0.01 * i as f64, t),
                ]
            })
            .collect::<Vec<_>>();
        let c = exchange_fractions(&ports, || Ok(())).unwrap();
        let b = (0..n + 2)
            .map(|i| 0.25 + (i as f64 * 1.7).sin().powi(2) * 2.0)
            .collect::<Vec<_>>();
        for (actual, reference) in network_balance(&c, &b, 0)
            .into_iter()
            .zip(dense_balance(&ports, &b, 0))
        {
            close(actual, reference, 4e-12);
        }
    }
    // R rounds to exactly one and every naive A*A product underflows. Divide
    // the two gap balance equations by epsilon: 6*g0-2*g1=s0 and
    // -2*g0+6*g1=s1. Their inverse gives these finite leading coefficients.
    let tiny = 1e-200;
    let c = exchange_fractions(&[[scalar(tiny, 2.0 * tiny); 2]; 3], || Ok(())).unwrap();
    let expected = [
        0.25, 0.0625, 0.25, 1.375, 0.125, 0.5, 0.5, 0.125, 1.375, 0.25,
    ];
    for (actual, expected) in c.into_iter().zip(expected) {
        close(actual[0], expected * tiny, 6e-15);
    }
}

#[test]
fn g1_coherent_three_body_spectrum_closes_all_terminals_and_mirrors_geometry() {
    let bounds = (4000.0, 20_000.0);
    let films = [
        [
            film(bounds, 2.7, 0.4, 0.42e-6, 1),
            film(bounds, 1.3, 0.02, 1.6e-6, 2),
        ],
        [
            film(bounds, 1.8, 0.09, 1.0e-6, 3),
            film(bounds, 3.5, 0.6, 0.25e-6, 4),
        ],
        [
            film(bounds, 2.1, 0.2, 0.7e-6, 5),
            film(bounds, 1.4, 0.03, 1.3e-6, 6),
        ],
    ];
    let bodies =
        core::array::from_fn::<_, 3, _>(|i| emitter(&films[i], bounds, i == 1, 10 + i as u8));
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 12,
        cosine_bins: 14,
    };
    let area = 0.025;
    let stack = prepare(&bodies, area, q);
    let temperatures = [390.0, 280.0, 440.0, 510.0, 220.0];
    let rates = stack
        .evaluate(
            &temperatures[..3],
            Some(temperatures[3]),
            Some(temperatures[4]),
            || Ok(()),
        )
        .unwrap();
    let mut expected = [0.0; 5];
    let width = (bounds.1 - bounds.0) / q.wavelength_bins as f64;
    for w in 0..q.wavelength_bins {
        let lambda = bounds.0 + (w as f64 + 0.5) * width;
        let b = temperatures.map(|t| planck(t, lambda));
        for a in 0..q.cosine_bins {
            let mu = (a as f64 + 0.5) / q.cosine_bins as f64;
            let ports = bodies.each_ref().map(|body| {
                [false, true].map(|opposite| {
                    body.reciprocal_vacuum_port_power(lambda, mu, opposite)
                        .unwrap()
                        .unwrap()
                })
            });
            for pol in 0..2 {
                for (out, direct) in expected.iter_mut().zip(dense_balance(&ports, &b, pol)) {
                    *out += area * PI * mu * width / q.cosine_bins as f64 * direct;
                }
            }
        }
    }
    for (actual, reference) in rates
        .body_outward_power_w()
        .iter()
        .copied()
        .chain(rates.bath_outward_power_w())
        .zip(expected)
    {
        close(actual, reference, 4e-12);
    }
    assert_eq!(stack.optical_evaluations(), 6 * 12 * 14);
    assert_eq!(stack.channel_spectral_nodes(), 10 * 12);
    assert_eq!(
        stack.optical_identities(),
        &bodies.each_ref().map(|b| b.optical_identity())
    );
    assert_eq!(stack.wavelength_bounds_nm(), bounds);
    assert_eq!(stack.area_m2(), area);
    assert_eq!(stack.quadrature(), q);
    assert_eq!(stack.body_count(), 3);
    let same = stack
        .evaluate(&[300.0; 3], Some(300.0), Some(300.0), || Ok(()))
        .unwrap();
    assert_eq!(same.body_outward_power_w(), [0.0; 3]);
    assert_eq!(same.bath_outward_power_w(), [0.0; 2]);
    assert!(
        same.body_pairs()
            .iter()
            .all(|p| p.exchange.outward_power_w == 0.0 && p.exchange.secant_w_k > 0.0)
    );
    let mirrored = core::array::from_fn::<_, 3, _>(|i| {
        emitter(&films[2 - i], bounds, 2 - i != 1, 20 + i as u8)
    });
    let mirrored = prepare(&mirrored, area, q)
        .evaluate(&[440.0, 280.0, 390.0], Some(220.0), Some(510.0), || Ok(()))
        .unwrap();
    for i in 0..3 {
        close(
            mirrored.body_outward_power_w()[i],
            rates.body_outward_power_w()[2 - i],
            5e-12,
        );
    }
    close(
        mirrored.bath_transfer.outward_power_w,
        -rates.bath_transfer.outward_power_w,
        5e-12,
    );
}

#[test]
fn g3_one_and_two_body_limits_reduce_to_existing_surface_and_pair_owners() {
    let bounds = (4000.0, 20_000.0);
    let first = [
        film(bounds, 2.5, 0.3, 0.6e-6, 30),
        film(bounds, 1.3, 0.05, 1.2e-6, 31),
    ];
    let second = [film(bounds, 1.8, 0.2, 0.8e-6, 32)];
    let left = emitter(&first, bounds, false, 33);
    let reversed = emitter(&first, bounds, true, 34);
    let right = emitter(&second, bounds, false, 35);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 8,
        cosine_bins: 10,
    };
    let ordinary = CoatingExchangeBudget {
        max_spectral_nodes: 8,
        max_optical_evaluations: 80,
    };
    let one = prepare(&[left.clone()], 0.02, q)
        .evaluate(&[400.0], Some(290.0), Some(310.0), || Ok(()))
        .unwrap();
    assert!(one.body_pairs().is_empty());
    for (i, body) in [&left, &reversed].into_iter().enumerate() {
        let expected = CoatingExchangeResponse::from_emitter(body, 0.02, q, ordinary, || Ok(()))
            .unwrap()
            .evaluate(400.0, Some([290.0, 310.0][i]), || Ok(()))
            .unwrap();
        close(
            one.body_baths()[0][i].outward_power_w,
            expected.outward_power_w,
            3e-12,
        );
    }
    let pair = CoatingTransmittingPairResponse::from_emitters(
        &reversed,
        &right,
        0.02,
        q,
        CoatingTransmittingPairBudget {
            max_channel_spectral_nodes: 48,
            max_optical_evaluations: 320,
        },
        || Ok(()),
    )
    .unwrap()
    .evaluate(400.0, 350.0, Some(290.0), Some(310.0), || Ok(()))
    .unwrap();
    let two = prepare(&[left, right], 0.02, q)
        .evaluate(&[400.0, 350.0], Some(290.0), Some(310.0), || Ok(()))
        .unwrap();
    close(
        two.body_pairs()[0].exchange.outward_power_w,
        pair.pair.outward_power_w,
        3e-12,
    );
    close(
        two.bath_transfer.outward_power_w,
        pair.bath_transfer.outward_power_w,
        3e-12,
    );
    for (actual, expected) in two
        .body_baths()
        .iter()
        .flatten()
        .zip(pair.first_baths.into_iter().chain(pair.second_baths))
    {
        close(actual.outward_power_w, expected.outward_power_w, 3e-12);
        close(actual.secant_w_k, expected.secant_w_k, 3e-12);
    }
}

#[test]
fn g1_opaque_shields_disconnect_nonadjacent_bodies_and_transparent_bodies_only_pass_baths() {
    let ports = [
        [scalar(0.2, 0.0), scalar(0.3, 0.0)],
        [scalar(0.4, 0.0), scalar(0.5, 0.0)],
        [scalar(0.6, 0.0), scalar(0.7, 0.0)],
    ];
    let c = exchange_fractions(&ports, || Ok(())).unwrap();
    let expected = [
        0.3 * 0.4 / (0.3 + 0.4 - 0.3 * 0.4),
        0.0,
        0.5 * 0.6 / (0.5 + 0.6 - 0.5 * 0.6),
        0.2,
        0.0,
        0.0,
        0.0,
        0.0,
        0.7,
        0.0,
    ];
    for (actual, expected) in c.into_iter().zip(expected) {
        close(actual[0], expected, 3e-15);
    }
    let mut trapped = ports;
    trapped[0][1] = scalar(0.0, 0.0);
    trapped[1][0] = scalar(0.0, 0.0);
    let trapped = exchange_fractions(&trapped, || Ok(())).unwrap();
    assert_eq!(trapped[0], [0.0; 2]);
    assert!(trapped[2][0] > 0.0);
    assert_eq!(trapped[3], [0.2; 2]);
    let bounds = (4000.0, 80_000.0);
    let body = emitter(&[], bounds, false, 40);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 256,
        cosine_bins: 1,
    };
    let stack = prepare(&vec![body; MAX_RADIATIVE_STACK_BODIES], 0.01, q);
    let rates = stack
        .evaluate(
            &[300.0; MAX_RADIATIVE_STACK_BODIES],
            Some(320.0),
            Some(280.0),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(
        rates.body_outward_power_w(),
        [0.0; MAX_RADIATIVE_STACK_BODIES]
    );
    assert!(rates.body_pairs().iter().all(|p| p.exchange == ZERO_RATE));
    assert!(rates.body_baths().iter().flatten().all(|p| *p == ZERO_RATE));
    assert!(rates.bath_transfer.outward_power_w > 0.0);
    assert_eq!(
        rates.bath_outward_power_w(),
        [
            rates.bath_transfer.outward_power_w,
            -rates.bath_transfer.outward_power_w
        ]
    );
    let dark = stack
        .evaluate(&[300.0; MAX_RADIATIVE_STACK_BODIES], None, None, || Ok(()))
        .unwrap();
    assert_eq!(dark.bath_transfer, ZERO_RATE);
}

#[test]
fn g1_weak_physical_absorption_and_planck_scaling_remain_representable() {
    let bounds = (380.0, 780.0);
    let films = [film(bounds, 1e-120, 100.0, 5e-6, 50)];
    let body = emitter(&films, bounds, false, 51);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 1,
        cosine_bins: 1,
    };
    let stack = prepare(&vec![body.clone(); 3], 1e-202, q);
    let rates = stack
        .evaluate(&[1e300, 300.0, 300.0], None, None, || Ok(()))
        .unwrap();
    let absorption = body.reciprocal_absorbed_fraction(580.0, 0.5).unwrap();
    let slope = 2.0 * 299_792_458.0 * 1.380_649e-23 * 1e27 / 580.0_f64.powi(4);
    let expected =
        (slope * 1e300) * (absorption[0] + absorption[1]) * (0.5 * 1e-202) * (PI * 0.5 * 400.0);
    close(
        rates.body_pairs()[0].exchange.outward_power_w,
        expected,
        4e-12,
    );
    assert!(rates.body_pairs()[0].exchange.secant_w_k > 0.0);
    assert_eq!(rates.body_pairs()[1].exchange, ZERO_RATE);
}

#[test]
fn g4_stack_admission_checks_budgets_domains_and_cancellation_without_partial_results() {
    let bounds = (4000.0, 20_000.0);
    let body = emitter(&[], bounds, false, 60);
    let q = CoatingExchangeQuadrature {
        wavelength_bins: 2,
        cosine_bins: 3,
    };
    for b in [
        CoatingRadiativeStackBudget {
            max_bodies: 2,
            ..budget(3, q)
        },
        CoatingRadiativeStackBudget {
            max_channel_spectral_nodes: 19,
            ..budget(3, q)
        },
        CoatingRadiativeStackBudget {
            max_optical_evaluations: 35,
            ..budget(3, q)
        },
    ] {
        let mut calls = 0;
        let response =
            CoatingRadiativeStackResponse::from_emitters(&vec![body.clone(); 3], 0.1, q, b, || {
                calls += 1;
                Ok(())
            });
        assert!(matches!(response, Err(CoatingExchangeError::Budget)));
        assert_eq!(calls, 0);
    }
    for bodies in [vec![], vec![body.clone(); MAX_RADIATIVE_STACK_BODIES + 1]] {
        assert!(matches!(
            CoatingRadiativeStackResponse::from_emitters(
                &bodies,
                0.1,
                q,
                budget(bodies.len(), q),
                || Ok(())
            ),
            Err(CoatingExchangeError::Budget)
        ));
    }
    let disjoint = emitter(&[], (30_000.0, 40_000.0), false, 61);
    assert!(matches!(
        CoatingRadiativeStackResponse::from_emitters(
            &[body.clone(), disjoint],
            0.1,
            q,
            budget(2, q),
            || Ok(())
        ),
        Err(CoatingExchangeError::WavelengthDomain)
    ));
    let nonvacuum = CoatingThermalEmission::try_transmitting(
        &[],
        &[],
        medium(bounds, 1.4),
        false,
        ContentHash([62; 32]),
    )
    .unwrap();
    assert!(matches!(
        CoatingRadiativeStackResponse::from_emitters(
            &[body.clone(), nonvacuum],
            0.1,
            q,
            budget(2, q),
            || Ok(())
        ),
        Err(CoatingExchangeError::VacuumTransmittingRequired)
    ));
    let mut calls = 0;
    let response = CoatingRadiativeStackResponse::from_emitters(
        &vec![body.clone(); 3],
        0.1,
        q,
        budget(3, q),
        || {
            calls += 1;
            if calls == 15 {
                Err(CoatingExchangeError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(response, Err(CoatingExchangeError::Cancelled)));
    assert_eq!(calls, 15);
    let stack = prepare(&vec![body; 3], 0.1, q);
    assert_eq!(
        stack.evaluate(&[300.0; 2], None, None, || Ok(())),
        Err(CoatingExchangeError::TemperatureChannels)
    );
    for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert_eq!(
            stack.evaluate(&[300.0, invalid, 300.0], None, None, || Ok(())),
            Err(CoatingExchangeError::Thermal(
                ThermalSurfaceError::Temperature
            ))
        );
        assert_eq!(
            stack.evaluate(&[300.0; 3], Some(invalid), None, || Ok(())),
            Err(CoatingExchangeError::Thermal(
                ThermalSurfaceError::Temperature
            ))
        );
    }
    assert_eq!(
        stack.evaluate(&[300.0; 3], None, None, || Err(
            CoatingExchangeError::Cancelled
        )),
        Err(CoatingExchangeError::Cancelled)
    );
}
