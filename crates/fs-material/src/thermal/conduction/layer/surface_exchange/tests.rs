//! Independent nonlinear two-surface balances, support, and source feedback.

use super::*;
use std::cell::Cell;

const AREA: f64 = 0.02;
const SIGMA: f64 = 5.670_374_419e-8;
const DISTANCES: [f64; 2] = [0.08, 0.05];
const ALPHA: [f64; 2] = [0.04, 0.02];
const FULL: [[f64; 2]; 2] = [[300.0, 700.0]; 2];

fn controls() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 64,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn curves() -> [ConductivityCurve; 2] {
    std::array::from_fn(|side| {
        ConductivityCurve::new(
            [300.0, 440.0, 600.0, 700.0]
                .map(|t| [t, ALPHA[side] * t])
                .to_vec(),
        )
        .unwrap()
    })
}

// Deliberately asymmetric surface coordinates: swapping hot/cold must not
// replace the original A/B inputs. Both film and emissivity depend on both
// ACTUAL material faces, not either bulk temperature.
fn coefficients([a, b]: [f64; 2]) -> [f64; 2] {
    [
        10.0 + 0.015 * a + 0.025 * b,
        0.4 + 0.0002 * (a - 450.0) + 0.0001 * (b - 450.0),
    ]
}

fn nonlinear(surface: [f64; 2]) -> SurfaceExchangeState<[f64; 2]> {
    let [a, b] = surface;
    let [htc, eps] = coefficients(surface);
    let delta = a - b;
    let fourth_factor = SIGMA * (a + b) * (a * a + b * b);
    let conductance = htc + eps * fourth_factor;
    SurfaceExchangeState {
        flux_a_to_b_w_m2: conductance * delta,
        conductance_w_m2_k: conductance,
        surface_derivatives_w_m2_k: [
            htc + 0.015 * delta + 4.0 * eps * SIGMA * a.powi(3) + 0.0002 * fourth_factor * delta,
            -htc + 0.025 * delta - 4.0 * eps * SIGMA * b.powi(3) + 0.0001 * fourth_factor * delta,
        ],
        source: surface,
    }
}

fn fixed(surface: [f64; 2], conductance: f64) -> SurfaceExchangeState<[f64; 2]> {
    SurfaceExchangeState {
        flux_a_to_b_w_m2: conductance * (surface[0] - surface[1]),
        conductance_w_m2_k: conductance,
        surface_derivatives_w_m2_k: [conductance, -conductance],
        source: surface,
    }
}

// Independent quadratic material primitives and square-root inverses, with
// direct fourth powers and scalar bisection. No production integral, inverse,
// interface solve, or derivative participates in this physical reference.
fn reference(temperatures: [f64; 2]) -> (f64, [f64; 2]) {
    let [a, b] = temperatures;
    let contact = (a * a - b * b) / (2.0 * (DISTANCES[0] / ALPHA[0] + DISTANCES[1] / ALPHA[1]));
    let surfaces = |flux: f64| {
        [
            (a * a - 2.0 * DISTANCES[0] * flux / ALPHA[0]).sqrt(),
            (b * b + 2.0 * DISTANCES[1] * flux / ALPHA[1]).sqrt(),
        ]
    };
    let (mut lo, mut hi) = (contact.min(0.0), contact.max(0.0));
    for _ in 0..100 {
        let flux = 0.5 * (lo + hi);
        let surface = surfaces(flux);
        let [htc, eps] = coefficients(surface);
        let source = htc * (surface[0] - surface[1])
            + eps * SIGMA * (surface[0].powi(4) - surface[1].powi(4));
        if source > flux {
            lo = flux;
        } else {
            hi = flux;
        }
    }
    let flux = 0.5 * (lo + hi);
    (AREA * flux, surfaces(flux))
}

#[test]
fn g1_actual_surface_law_matches_independent_integrals_and_signed_endpoint_order() {
    let curves = curves();
    let model = TwoLayerConduction::new([&curves[0], &curves[1]], DISTANCES, AREA, 0.0).unwrap();
    for temperatures in [[650.0, 330.0], [350.0, 640.0]] {
        let mut calls = 0;
        let response = model
            .evaluate_surface_exchange(
                temperatures,
                FULL,
                controls(),
                |surface| {
                    calls += 1;
                    Ok(nonlinear(surface))
                },
                || Ok(()),
            )
            .unwrap();
        let (power, surfaces) = reference(temperatures);
        near(response.power_a_to_b_w, power, 2e-8);
        for side in 0..2 {
            near(response.interface_temperature_k[side], surfaces[side], 2e-9);
        }
        near(
            response.conductance_w_k * (temperatures[0] - temperatures[1]),
            power,
            2e-8,
        );
        near(response.layer_heat_residual_w, 0.0, 2e-8);
        near(response.surface_heat_residual_w, 0.0, 2e-8);
        assert_eq!(response.source, response.interface_temperature_k);
        assert_eq!(calls, response.iterations);
        assert!(calls <= controls().max_iterations);

        let partials = nonlinear(surfaces).surface_derivatives_w_m2_k;
        let denominator = 1.0 + partials[0] * DISTANCES[0] / (ALPHA[0] * surfaces[0])
            - partials[1] * DISTANCES[1] / (ALPHA[1] * surfaces[1]);
        near(
            response.source_flux_derivative_scale_m2,
            AREA / denominator,
            2e-12,
        );
        for side in 0..2 {
            let analytic =
                AREA * partials[side] * temperatures[side] / (surfaces[side] * denominator);
            near(response.endpoint_derivative_w_k[side], analytic, 2e-10);
            let mut plus = temperatures;
            let mut minus = temperatures;
            plus[side] += 0.001;
            minus[side] -= 0.001;
            near(
                response.endpoint_derivative_w_k[side],
                (reference(plus).0 - reference(minus).0) / 0.002,
                3e-7,
            );
        }
    }
}

#[test]
fn g1_disjoint_surface_support_retains_both_flux_bounds_and_exact_supported_endpoints() {
    let a = ConductivityCurve::new(vec![[500.0, 10.0], [700.0, 10.0]]).unwrap();
    let b = ConductivityCurve::new(vec![[250.0, 20.0], [400.0, 20.0]]).unwrap();
    for reversed in [false, true] {
        let (model, bulk, expected, sign) = if reversed {
            (
                TwoLayerConduction::new([&b, &a], [0.2, 0.1], AREA, 0.0).unwrap(),
                [300.0, 600.0],
                [350.0, 550.0],
                -1.0,
            )
        } else {
            (
                TwoLayerConduction::new([&a, &b], [0.1, 0.2], AREA, 0.0).unwrap(),
                [600.0, 300.0],
                [550.0, 350.0],
                1.0,
            )
        };
        // Interior root, exact upper/lower supported flux, and a singleton.
        // None of these source intervals contains either original bulk state.
        for original in [
            [[540.0, 560.0], [340.0, 360.0]],
            [[550.0, 560.0], [340.0, 350.0]],
            [[540.0, 550.0], [350.0, 360.0]],
            [[550.0, 550.0], [350.0, 350.0]],
        ] {
            let bounds = if reversed {
                [original[1], original[0]]
            } else {
                original
            };
            let mut calls = 0;
            let response = model
                .evaluate_surface_exchange(
                    bulk,
                    bounds,
                    controls(),
                    |surface| {
                        calls += 1;
                        for side in 0..2 {
                            assert!((bounds[side][0]..=bounds[side][1]).contains(&surface[side]));
                        }
                        Ok(fixed(surface, 25.0))
                    },
                    || Ok(()),
                )
                .unwrap();
            near(response.power_a_to_b_w, sign * 100.0, 1e-9);
            assert_eq!(response.interface_temperature_k, expected);
            assert_eq!(response.source, expected);
            assert_eq!(response.iterations, calls);
        }
    }

    let model = TwoLayerConduction::new([&a, &b], [0.1, 0.2], AREA, 0.0).unwrap();
    let bounds = [[560.0, 570.0], [330.0, 340.0]];
    let mut calls = 0;
    assert!(
        model
            .evaluate_surface_exchange(
                [600.0, 300.0],
                bounds,
                controls(),
                |surface| {
                    calls += 1;
                    for side in 0..2 {
                        assert!((bounds[side][0]..=bounds[side][1]).contains(&surface[side]));
                    }
                    Ok(fixed(surface, 25.0))
                },
                || Ok(()),
            )
            .is_err()
    );
    assert!(calls > 0 && calls <= controls().max_iterations);
    // These two otherwise supported face intervals cannot carry the SAME
    // conserved flux, so refuse before querying any source at all.
    let mut impossible_calls = 0;
    assert!(
        model
            .evaluate_surface_exchange(
                [600.0, 300.0],
                [[540.0, 550.0], [330.0, 340.0]],
                controls(),
                |surface| {
                    impossible_calls += 1;
                    Ok(fixed(surface, 25.0))
                },
                || Ok(()),
            )
            .is_err()
    );
    assert_eq!(impossible_calls, 0);
}

#[test]
fn g1_passive_signed_feedback_and_zero_limits_preserve_the_actual_differential() {
    let curve = ConductivityCurve::new(vec![[300.0, 10.0], [700.0, 10.0]]).unwrap();
    let model = TwoLayerConduction::new([&curve; 2], [0.1; 2], AREA, 0.0).unwrap();
    let feedback = |surface: [f64; 2], slope: f64| {
        let conductance = 25.0 + slope * (surface[0] - 550.0);
        SurfaceExchangeState {
            surface_derivatives_w_m2_k: [
                conductance + slope * (surface[0] - surface[1]),
                -conductance,
            ],
            ..fixed(surface, conductance)
        }
    };
    let response = model
        .evaluate_surface_exchange(
            [600.0, 300.0],
            [[300.0, 650.0]; 2],
            controls(),
            |surface| Ok(feedback(surface, -0.2)),
            || Ok(()),
        )
        .unwrap();
    near(response.power_a_to_b_w, 100.0, 2e-8);
    near(response.interface_temperature_k[0], 550.0, 2e-9);
    near(response.interface_temperature_k[1], 350.0, 2e-9);
    near(
        response.endpoint_derivative_w_k[0],
        AREA * -15.0 / 1.1,
        2e-10,
    );
    near(
        response.endpoint_derivative_w_k[1],
        AREA * -25.0 / 1.1,
        2e-10,
    );
    near(response.source_flux_derivative_scale_m2, AREA / 1.1, 2e-12);
    assert!(response.power_a_to_b_w > 0.0 && response.endpoint_derivative_w_k[0] < 0.0);
    // The same positive conductance at this exact root can have unstable
    // feedback: D=1-175*.01+25*.01=-.5. Passive signs alone cannot admit it.
    assert!(matches!(
        model.evaluate_surface_exchange(
            [600.0, 300.0],
            [[550.0, 550.0], [350.0, 350.0]],
            controls(),
            |surface| Ok(feedback(surface, -1.0)),
            || Ok(()),
        ),
        Err(DepositionError::InvalidInput(_))
    ));

    let insulated = model
        .evaluate_surface_exchange(
            [600.0, 300.0],
            FULL,
            controls(),
            |surface| Ok(fixed(surface, 0.0)),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(insulated.power_a_to_b_w, 0.0);
    assert_eq!(insulated.conductance_w_k, 0.0);
    assert_eq!(insulated.interface_temperature_k, [600.0, 300.0]);
    assert_eq!(insulated.endpoint_derivative_w_k, [0.0; 2]);
    assert_eq!(insulated.source, [600.0, 300.0]);
    let equilibrium = model
        .evaluate_surface_exchange(
            [500.0; 2],
            FULL,
            controls(),
            |surface| Ok(fixed(surface, 25.0)),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(equilibrium.power_a_to_b_w, 0.0);
    assert_eq!(equilibrium.interface_temperature_k, [500.0; 2]);
    let conductance = AREA / (0.01 + 0.04 + 0.01);
    near(equilibrium.conductance_w_k, conductance, 2e-12);
    near(equilibrium.endpoint_derivative_w_k[0], conductance, 2e-12);
    near(equilibrium.endpoint_derivative_w_k[1], -conductance, 2e-12);
    near(
        equilibrium.source_flux_derivative_scale_m2,
        AREA / 1.5,
        2e-12,
    );
}

#[test]
fn g4_original_source_errors_cancellation_and_callback_budget_are_retained() {
    let curves = curves();
    let model = TwoLayerConduction::new([&curves[0], &curves[1]], DISTANCES, AREA, 0.0).unwrap();
    assert!(matches!(
        model.evaluate_surface_exchange(
            [650.0, 330.0],
            FULL,
            controls(),
            |_| Err::<SurfaceExchangeState<()>, _>(DepositionError::InvalidInput(
                "original surface source refusal"
            )),
            || Ok(()),
        ),
        Err(DepositionError::InvalidInput(
            "original surface source refusal"
        ))
    ));
    let queries = Cell::new(0);
    let cancelled = model.evaluate_surface_exchange(
        [650.0, 330.0],
        FULL,
        controls(),
        |surface| {
            queries.set(queries.get() + 1);
            Ok(nonlinear(surface))
        },
        || {
            if queries.get() > 0 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(cancelled, Err(DepositionError::Cancelled)));
    assert_eq!(queries.get(), 1);

    let mut calls = 0;
    assert!(
        model
            .evaluate_surface_exchange(
                [650.0, 330.0],
                FULL,
                LayerConductionConfig {
                    max_iterations: 1,
                    ..controls()
                },
                |surface| {
                    calls += 1;
                    Ok(nonlinear(surface))
                },
                || Ok(()),
            )
            .is_err()
    );
    assert_eq!(calls, 1);
}
