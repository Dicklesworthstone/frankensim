//! G1: independent quadrature and spatial derivatives of the outgoing field.
use super::*;

fn gauss(count: usize) -> Vec<(f64, f64)> {
    (0..count)
        .map(|index| {
            let mut x = (std::f64::consts::PI * (index as f64 + 0.75) / (count as f64 + 0.5)).cos();
            let mut derivative = 0.0;
            for _ in 0..32 {
                let (mut a, mut b) = (1.0, x);
                for degree in 2..=count {
                    let next =
                        ((2 * degree - 1) as f64 * x * b - (degree - 1) as f64 * a) / degree as f64;
                    a = b;
                    b = next;
                }
                derivative = count as f64 * (x * b - a) / (x * x - 1.0);
                let update = b / derivative;
                x -= update;
                if update.abs() < 2e-16 {
                    break;
                }
            }
            (
                0.5 * (x + 1.0),
                1.0 / ((1.0 - x * x) * derivative * derivative),
            )
        })
        .collect()
}

// Explicit G and source-normal derivative on z=0. This does not call the
// production wave, kernel, derivative, subdivision or quadrature functions.
fn layers(delta: Point, k: C64) -> (C64, C64) {
    let r = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
    let phase = k.re * r;
    let g = C64::new(phase.cos(), phase.sin())
        .scale((-k.im * r).exp() / (4.0 * std::f64::consts::PI * r));
    let d = g * (C64::ONE - C64::new(-k.im, k.re).scale(r)).scale(delta[2] / (r * r));
    (g, d)
}

fn reference_pressure(x: Point, k: C64, p: C64, q: C64, nodes: &[(f64, f64)]) -> C64 {
    let mut result = C64::ZERO;
    for &(a, wa) in nodes {
        for &(b, wb) in nodes {
            let y = [a * (1.0 - b), a * b, 0.0];
            let (s, d) = layers(std::array::from_fn(|axis| x[axis] - y[axis]), k);
            result = result + (d * p - s * q).scale(wa * wb * a);
        }
    }
    result
}

fn reference_gradient(x: Point, k: C64, p: C64, q: C64, nodes: &[(f64, f64)]) -> [C64; 3] {
    let step = 1e-5;
    std::array::from_fn(|axis| {
        let mut positive = x;
        let mut negative = x;
        positive[axis] += step;
        negative[axis] -= step;
        (reference_pressure(positive, k, p, q, nodes)
            - reference_pressure(negative, k, p, q, nodes))
        .scale(0.5 / step)
    })
}

fn relative_vector(actual: [C64; 3], expected: [C64; 3]) -> f64 {
    (actual
        .iter()
        .zip(expected)
        .map(|(&a, b)| (a - b).norm_sq())
        .sum::<f64>()
        / expected.iter().map(|v| v.norm_sq()).sum::<f64>())
    .sqrt()
}

#[test]
fn complex_green_layer_gradients_match_independent_spatial_differences() {
    for k in [C64::new(4.2, 1.3), C64::new(12.0, 3.0)] {
        let wave = DampedWaveNumber {
            k,
            phase_origin_m: 0.0,
        };
        let delta = [0.13, -0.21, 0.07];
        let sample = kernel(wave, delta, [0.0, 0.0, 1.0]).unwrap();
        let mut integral = Integral::default();
        integral.accumulate(&sample, 1.0, true).unwrap();
        let step = 1e-6;
        for axis in 0..3 {
            let mut a = delta;
            let mut b = delta;
            a[axis] += step;
            b[axis] -= step;
            let (sa, da) = layers(a, k);
            let (sb, db) = layers(b, k);
            let ds = (sa - sb).scale(0.5 / step);
            let dd = (da - db).scale(0.5 / step);
            assert!((integral.gradient.s[axis] - ds).abs() < 2e-8 * ds.abs());
            assert!((integral.gradient.d[axis] - dd).abs() < 2e-8 * dd.abs());
        }
    }
}

#[test]
fn adaptive_absorbing_triangle_matches_independent_refined_pressure_and_velocity() {
    let triangle = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let normal = [0.0, 0.0, 1.0];
    let pressure = C64::new(0.8, -0.4);
    let rho_omega = 415.3;
    let q = C64::new(0.002, 0.001) * C64::new(0.0, rho_omega);
    let inverse = C64::new(0.0, -1.0 / rho_omega);
    let coarse = gauss(96);
    let refined = gauss(128);
    for k in [C64::new(4.2, 1.3), C64::new(12.0, 3.0)] {
        for point in [
            [0.31, 0.23, 0.08],
            [0.48, -0.035, 0.045],
            [0.21, 0.34, -0.065],
        ] {
            let expected = reference_pressure(point, k, pressure, q, &refined);
            let reference_v = reference_gradient(point, k, pressure, q, &refined)
                .map(|gradient| gradient * inverse);
            assert!(
                (reference_pressure(point, k, pressure, q, &coarse) - expected).abs()
                    < 1e-7 * expected.abs()
            );
            assert!(
                relative_vector(
                    reference_gradient(point, k, pressure, q, &coarse)
                        .map(|gradient| gradient * inverse),
                    reference_v,
                ) < 1e-7
            );
            let mut work = Work {
                options: Options {
                    relative_tolerance: 1e-9,
                    ..Options::default()
                },
                evaluations: 0,
            };
            let row = integrate_fields_cancellable(
                DampedWaveNumber {
                    k,
                    phase_origin_m: 0.0,
                },
                point,
                triangle,
                normal,
                0,
                &mut work,
                true,
                &mut || false,
            )
            .unwrap();
            let actual = row.d * pressure - row.s * q;
            let velocity = std::array::from_fn(|axis| {
                (row.gradient.d[axis] * pressure - row.gradient.s[axis] * q) * inverse
            });
            assert!((actual - expected).abs() < 3e-7 * expected.abs());
            assert!(relative_vector(velocity, reference_v) < 3e-7);
            assert!(work.evaluations > 80);
            assert!(row.error_d.is_finite() && row.error_s.is_finite());
        }
    }
}
