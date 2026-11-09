use super::*;

// Independent symmetric two-port scattering eigenvalues. This subtraction is
// deliberately confined to moderate-loss TEST oracles, never the product path.
fn scattering_loss(epsilon: C64, mu: f64, cosine: f64, depth: f64, te: bool) -> [f64; 2] {
    let q = (epsilon.scale(mu) - C64::from_re((1.0 - cosine) * (1.0 + cosine))).sqrt();
    let p = C64::new(det::cos(q.re * depth), -det::sin(q.re * depth)).scale(det::exp(q.im * depth));
    let (a, b) = if te {
        (C64::from_re(mu * cosine), q)
    } else {
        (q, epsilon.scale(cosine))
    };
    let even =
        (a * (C64::ONE + p) - b * (C64::ONE - p)) / (a * (C64::ONE + p) + b * (C64::ONE - p));
    let odd = (a * (C64::ONE - p) - b * (C64::ONE + p)) / (a * (C64::ONE - p) + b * (C64::ONE + p));
    [1.0 - even.norm_sq(), 1.0 - odd.norm_sq()]
}

fn direction(
    epsilon: C64,
    mu: f64,
    cosine: f64,
    depth: f64,
    depsilon: C64,
    dcosine: f64,
    ddepth: f64,
    te: bool,
) -> SlabLossDirection {
    let q = (epsilon.scale(mu) - C64::from_re((1.0 - cosine) * (1.0 + cosine))).sqrt();
    SlabLossDirection {
        normal: q,
        normal_scaled_derivative: (depsilon.scale(mu) + C64::from_re(2.0 * cosine * dcosine))
            / q.scale(2.0),
        permittivity_ratio: epsilon,
        permittivity_ratio_scaled_derivative: depsilon,
        permeability_ratio: mu,
        cosine,
        cosine_scaled_derivative: dcosine,
        depth,
        depth_scaled_derivative: ddepth,
        phase_parts: Dd::from_f64(q.re) * Dd::from_f64(depth),
        loss: scattering_loss(epsilon, mu, cosine, depth, te),
        transverse_electric: te,
    }
}

fn check_scattering(input: SlabLossDirection, step: f64, tolerance: f64) {
    let sample = |t: f64| {
        scattering_loss(
            input.permittivity_ratio + input.permittivity_ratio_scaled_derivative.scale(t),
            input.permeability_ratio,
            input.cosine + input.cosine_scaled_derivative * t,
            input.depth + input.depth_scaled_derivative * t,
            input.transverse_electric,
        )
    };
    let values = [
        sample(-2.0 * step),
        sample(-step),
        sample(step),
        sample(2.0 * step),
    ];
    let actual = scaled_loss_derivatives(input).unwrap();
    for i in 0..2 {
        let expected =
            (values[0][i] - 8.0 * values[1][i] + 8.0 * values[2][i] - values[3][i]) / (12.0 * step);
        assert!(
            (actual[i] - expected).abs() <= tolerance * (1.0 + expected.abs()),
            "mode={i}, actual={}, expected={expected}",
            actual[i]
        );
    }
}

#[test]
fn original_even_odd_loss_tangents_match_independent_scattering_eigenvalues() {
    for te in [true, false] {
        for depth in [0.03, 0.43, 4.2] {
            check_scattering(
                direction(
                    C64::new(2.2, -0.35),
                    0.8,
                    0.68,
                    depth,
                    C64::new(0.2, -0.09),
                    0.11,
                    1.3 * depth,
                    te,
                ),
                2e-4,
                3e-9,
            );
        }
    }
}

#[test]
fn very_thin_tiny_loss_keeps_odd_cubic_work_and_material_direction() {
    let gamma = 1e-180;
    let depth = 1e-6;
    let material_slope = 0.7;
    let x = depth * depth;
    // Independent first order in gamma at a matched lossless base. Errors
    // O(gamma^2) are many orders below f64; no 1-|S|^2 subtraction is used.
    let odd_phase = depth * x * (1.0 / 6.0 - x * (1.0 / 120.0 - x / 5040.0));
    let even_phase = depth + det::sin(depth);
    let loss = [gamma * even_phase, gamma * odd_phase];
    let expected = [
        gamma * (material_slope * even_phase + depth * (1.0 + det::cos(depth))),
        gamma * (material_slope * odd_phase + 2.0 * depth * det::sin(0.5 * depth).powi(2)),
    ];
    for te in [true, false] {
        let mut input = direction(
            C64::new(1.0, -gamma),
            1.0,
            1.0,
            depth,
            C64::new(0.0, -material_slope * gamma),
            0.0,
            depth,
            te,
        );
        input.loss = loss;
        let actual = scaled_loss_derivatives(input).unwrap();
        for i in 0..2 {
            assert!(actual[i] > 0.0);
            assert!(
                (actual[i] / expected[i] - 1.0).abs() < 2e-12,
                "mode={i}, actual={}, expected={}",
                actual[i],
                expected[i]
            );
        }
    }
}

#[test]
fn nearly_critical_lossy_film_retains_geometric_and_constitutive_directions() {
    for te in [true, false] {
        check_scattering(
            direction(
                C64::new(0.75, -1e-6),
                1.0,
                0.5,
                0.25,
                C64::new(0.08, -3e-7),
                0.13,
                0.4,
                te,
            ),
            1e-3,
            2e-10,
        );
    }
}

#[test]
fn zero_width_is_exact_and_opaque_or_unrepresentable_jets_refuse() {
    let mut input = direction(
        C64::new(2.2, -0.35),
        0.8,
        0.68,
        0.0,
        C64::ZERO,
        0.0,
        0.0,
        true,
    );
    input.loss = [0.0; 2];
    assert_eq!(scaled_loss_derivatives(input), Some([0.0; 2]));
    input.depth_scaled_derivative = 1.0;
    assert!(scaled_loss_derivatives(input).is_none());
    let mut input = direction(
        C64::new(2.2, -0.35),
        0.8,
        0.68,
        10_000.0,
        C64::ZERO,
        0.0,
        10_000.0,
        true,
    );
    assert!(scaled_loss_derivatives(input).is_none());
    input.depth = 0.4;
    input.normal_scaled_derivative = C64::new(f64::INFINITY, 0.0);
    assert!(scaled_loss_derivatives(input).is_none());
    input.normal_scaled_derivative = C64::ZERO;
    input.normal = C64::ZERO;
    assert!(scaled_loss_derivatives(input).is_none());
}
