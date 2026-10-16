//! Positive equal-area quadrature for the existing flat-board displacement.
//! Refinement only samples that same P1 or edge-cubic field more densely; it
//! never changes structural nodes, modal coordinates, bridge maps or forces.

/// Uniform barycentric subdivision. Level zero deliberately returns the
/// original coordinates directly, preserving its floating-point path.
pub(super) fn triangle_weights(levels: usize) -> Vec<[f64; 3]> {
    let rule = [
        [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
    ];
    if levels == 0 {
        return rule.to_vec();
    }
    let mut cells = vec![[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]];
    for _ in 0..levels {
        let mut children = Vec::with_capacity(4 * cells.len());
        for [a, b, c] in cells {
            let midpoint = |x: [f64; 3], y: [f64; 3]| std::array::from_fn(|i| 0.5 * (x[i] + y[i]));
            let ab = midpoint(a, b);
            let bc = midpoint(b, c);
            let ca = midpoint(c, a);
            children.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
        }
        cells = children;
    }
    cells
        .iter()
        .flat_map(|cell| {
            rule.map(|weights| {
                std::array::from_fn(|i| (0..3).map(|j| weights[j] * cell[j][i]).sum())
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_plate::{edge_cubic_transverse_mean_shape, edge_cubic_transverse_shape};
    use std::f64::consts::TAU;

    #[test]
    fn refined_acoustic_rule_integrates_independent_barycentric_moments() {
        for level in 0..=3 {
            let points = triangle_weights(level);
            let mean = |f: fn([f64; 3]) -> f64| {
                points.iter().map(|&p| f(p)).sum::<f64>() / points.len() as f64
            };
            assert_eq!(points.len(), 3 * 4usize.pow(level as u32));
            assert!(
                points
                    .iter()
                    .all(|p| p.iter().all(|&x| x > 0.0)
                        && (p.iter().sum::<f64>() - 1.0).abs() < 1e-14)
            );
            assert!((mean(|p| p[0]) - 1.0 / 3.0).abs() < 1e-14);
            assert!((mean(|p| p[0] * p[0]) - 1.0 / 6.0).abs() < 1e-14);
            assert!((mean(|p| p[0] * p[1]) - 1.0 / 12.0).abs() < 1e-14);
        }
    }

    #[test]
    fn refinement_preserves_every_signed_edge_cubic_area_moment() {
        // Every local displacement AND slope coordinate is checked against
        // the owner's analytic Bernstein integral, including skew triangles.
        // A nodal P1 substitution would lose these nonzero slope moments.
        for (x, y) in [
            ([0.0, 0.2, 0.02], [0.0, 0.01, 0.18]),
            ([-0.3, 0.8, -0.1], [-0.2, 0.1, 0.4]),
        ] {
            let exact = edge_cubic_transverse_mean_shape(&x, &y);
            assert!(
                exact
                    .iter()
                    .enumerate()
                    .any(|(i, v)| i % 3 != 0 && v.abs() > 1e-5)
            );
            for level in 0..=3 {
                let points = triangle_weights(level);
                let mut actual = [0.0; 9];
                for weights in &points {
                    for (sum, value) in actual
                        .iter_mut()
                        .zip(edge_cubic_transverse_shape(&x, &y, *weights))
                    {
                        *sum += value / points.len() as f64;
                    }
                }
                for local in 0..9 {
                    assert!(
                        (actual[local] - exact[local]).abs() < 2e-14,
                        "level {level}, coordinate {local}: {} vs {}",
                        actual[local],
                        exact[local]
                    );
                }
            }
        }
    }

    // Independently prescribed quadratic normal acceleration [m/s²]. The
    // edge-cubic field reproduces it from vertex values and physical slopes;
    // replacing that field by P1 changes both its shape and its radiation.
    fn acceleration(x: f64, y: f64) -> (f64, f64, f64) {
        (
            (x - 0.09).powi(2) - 0.7 * (y - 0.08).powi(2) + 0.04 * (x - y),
            2.0 * (x - 0.09) + 0.04,
            -1.4 * (y - 0.08) - 0.04,
        )
    }
    fn kernel(x: f64, y: f64, receiver: [f64; 3], frequency: f64) -> [f64; 2] {
        let radius =
            ((x - receiver[0]).powi(2) + (y - receiver[1]).powi(2) + receiver[2] * receiver[2])
                .sqrt();
        let phase = -TAU * frequency * radius / 343.0;
        let gain = 1.2 / (TAU * radius);
        [gain * phase.cos(), gain * phase.sin()]
    }
    fn rule_pressure(
        x: [f64; 3],
        y: [f64; 3],
        level: usize,
        receiver: [f64; 3],
        frequency: f64,
    ) -> [f64; 2] {
        let points = triangle_weights(level);
        let area = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]));
        let mut dofs = [0.0; 9];
        for node in 0..3 {
            let (a, ax, ay) = acceleration(x[node], y[node]);
            dofs[3 * node..3 * node + 3].copy_from_slice(&[a, ax, ay]);
        }
        let mut pressure = [0.0; 2];
        for weights in &points {
            let px = (0..3).map(|i| weights[i] * x[i]).sum();
            let py = (0..3).map(|i| weights[i] * y[i]).sum();
            let a = edge_cubic_transverse_shape(&x, &y, *weights)
                .iter()
                .zip(dofs)
                .map(|(shape, value)| shape * value)
                .sum::<f64>();
            for (sum, value) in pressure.iter_mut().zip(kernel(px, py, receiver, frequency)) {
                *sum += area / points.len() as f64 * a * value;
            }
        }
        pressure
    }
    fn reference_pressure(
        x: [f64; 3],
        y: [f64; 3],
        intervals: usize,
        receiver: [f64; 3],
        frequency: f64,
    ) -> [f64; 2] {
        // Independent tensor Simpson integration through a Duffy map. It
        // directly evaluates the polynomial and 1/r * exp(-ikr), using none
        // of the production triangle rule, mode shapes or delay-filter taps.
        let jacobian = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
        let weight = |i| {
            if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            }
        };
        let mut pressure = [0.0; 2];
        for i in 0..=intervals {
            let u = i as f64 / intervals as f64;
            for j in 0..=intervals {
                let v = j as f64 / intervals as f64;
                let weights = [(1.0 - u) * (1.0 - v), u, (1.0 - u) * v];
                let px = (0..3).map(|k| weights[k] * x[k]).sum();
                let py = (0..3).map(|k| weights[k] * y[k]).sum();
                let scale = weight(i) * weight(j) * (1.0 - u) * jacobian
                    / (9.0 * (intervals * intervals) as f64)
                    * acceleration(px, py).0;
                for (sum, value) in pressure.iter_mut().zip(kernel(px, py, receiver, frequency)) {
                    *sum += scale * value;
                }
            }
        }
        pressure
    }
    fn distance_squared(a: [f64; 2], b: [f64; 2]) -> f64 {
        (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)
    }

    #[test]
    fn refined_cubic_rayleigh_transfer_converges_to_an_independent_surface_integral() {
        let (x, y) = ([0.0, 0.2, 0.02], [0.0, 0.01, 0.18]);
        let mut error = [0.0; 4];
        let (mut reference_error, mut signal) = (0.0, 0.0);
        for receiver in [[0.4, -0.1, 0.5], [-0.2, 0.3, 0.4]] {
            for frequency in [1500.0, 4000.0, 8000.0, 12000.0] {
                let reference = reference_pressure(x, y, 512, receiver, frequency);
                let coarser = reference_pressure(x, y, 256, receiver, frequency);
                reference_error += distance_squared(reference, coarser);
                signal += distance_squared(reference, [0.0; 2]);
                for (level, total) in error.iter_mut().enumerate() {
                    *total += distance_squared(
                        rule_pressure(x, y, level, receiver, frequency),
                        reference,
                    );
                }
            }
        }
        assert!(signal > 0.0);
        let reference_rms = (reference_error / signal).sqrt();
        let rms = error.map(|value| (value / signal).sqrt());
        assert!(
            reference_rms < 2e-6,
            "independent integration did not converge: {reference_rms}"
        );
        assert!(
            rms[0] > 0.03,
            "reference must expose unresolved high-band spatial phase: {rms:?}"
        );
        assert!(
            rms[3] < 0.02 && rms[3] < rms[0] / 20.0,
            "cubic Rayleigh complex RMS error did not improve: {rms:?}"
        );
        eprintln!(
            "cubic Rayleigh spatial relative RMS by level: {rms:?}; reference refinement {reference_rms:e}"
        );
    }
}
