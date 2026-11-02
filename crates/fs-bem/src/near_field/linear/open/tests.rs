//! G1/G3: original P1 recovery and an independent absorbing triangle observer.
use super::*;

const TRIANGLE: Triangle = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
const NORMAL: Point = [0.0, 0.0, 1.0];
const JUMP: [C64; 3] = [C64::new(0.7, -0.4), C64::new(-0.2, 0.5), C64::new(1.1, 0.3)];

fn prepare(k: C64, points: &[Point]) -> OpenRows {
    OpenRows::new(
        &[TRIANGLE],
        &[NORMAL],
        DampedWaveNumber {
            k,
            phase_origin_m: 0.0,
        },
        423.7,
        points,
        1e-7,
        Options {
            relative_tolerance: 1e-9,
            ..Options::default()
        },
        required_row_bytes(1, points.len()).unwrap(),
        &mut || false,
    )
    .unwrap()
}

fn gauss(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut derivative = 0.0;
            for _ in 0..32 {
                let (mut a, mut b) = (1.0, x);
                for degree in 2..=n {
                    let next =
                        ((2 * degree - 1) as f64 * x * b - (degree - 1) as f64 * a) / degree as f64;
                    a = b;
                    b = next;
                }
                derivative = n as f64 * (x * b - a) / (x * x - 1.0);
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

// Explicit double-layer potential on z=0 with independently generated Gauss
// points. No production Green function, moment, subdivision or derivative.
fn reference(x: Point, k: C64, rule: &[(f64, f64)]) -> C64 {
    let mut result = C64::ZERO;
    for &(a, wa) in rule {
        for &(b, wb) in rule {
            let y = [a * (1.0 - b), a * b, 0.0];
            let r = (0..3)
                .map(|axis| (x[axis] - y[axis]).powi(2))
                .sum::<f64>()
                .sqrt();
            let phase = k.re * r;
            let green = C64::new(phase.cos(), phase.sin())
                .scale((-k.im * r).exp() / (4.0 * std::f64::consts::PI * r));
            let kernel = green * C64::new(1.0 / r + k.im, -k.re).scale(x[2] / r);
            let p = JUMP[0].scale(1.0 - a) + JUMP[1].scale(a * (1.0 - b)) + JUMP[2].scale(a * b);
            result = result + (kernel * p).scale(wa * wb * a);
        }
    }
    result
}

fn derivative(x: Point, k: C64, rule: &[(f64, f64)]) -> [C64; 3] {
    let h = 4e-5;
    core::array::from_fn(|axis| {
        let value = |step: f64| {
            let mut point = x;
            point[axis] += step * h;
            reference(point, k, rule)
        };
        (-value(2.0) + value(1.0).scale(8.0) - value(-1.0).scale(8.0) + value(-2.0))
            .scale(1.0 / (12.0 * h))
    })
}

#[test]
fn zero_absorption_recovers_original_open_p1_pressure_on_both_sides() {
    let points = [[0.23, 0.29, 0.14], [0.23, 0.29, -0.14]];
    let prepared = prepare(C64::from_re(4.2), &points);
    let actual = prepared.evaluate(|_| JUMP, &mut || false).unwrap();
    let original = crate::near_field::linear::open_pressure(
        4.2,
        &[TRIANGLE],
        &[NORMAL],
        &[JUMP],
        &points,
        1e-7,
        Options {
            relative_tolerance: 1e-9,
            ..Options::default()
        },
        || false,
    )
    .unwrap();
    for (&a, b) in actual.scalar.pressure.iter().zip(original) {
        assert!((a - b).abs() < 3e-8 * b.abs());
    }
    assert!((actual.scalar.pressure[0] + actual.scalar.pressure[1]).abs() < 1e-13);
    let repeated = prepared.evaluate(|_| JUMP, &mut || false).unwrap();
    assert_eq!(actual.scalar.pressure, repeated.scalar.pressure);
    assert_eq!(actual.particle_velocity_m_s, repeated.particle_velocity_m_s);
}

#[test]
fn absorbing_p1_pressure_and_particle_velocity_match_independent_quadrature() {
    let k = C64::new(4.2, 0.8);
    let points = [[0.23, 0.29, 0.1], [0.7, -0.09, -0.11]];
    let coarse = gauss(64);
    let refined = gauss(96);
    let prepared = prepare(k, &points);
    let actual = prepared.evaluate(|_| JUMP, &mut || false).unwrap();
    assert!(prepared.kernel_evaluations() >= 80 * points.len());
    for (receiver, &point) in points.iter().enumerate() {
        let expected = reference(point, k, &refined);
        assert!((reference(point, k, &coarse) - expected).abs() < 1e-8 * expected.abs());
        assert!((actual.scalar.pressure[receiver] - expected).abs() < 3e-7 * expected.abs());
        let gradient = derivative(point, k, &refined);
        let velocity = gradient.map(|value| value * C64::new(0.0, -1.0 / 423.7));
        let error = (0..3)
            .map(|axis| (actual.particle_velocity_m_s[receiver][axis] - velocity[axis]).norm_sq())
            .sum::<f64>();
        let norm_sq = velocity.iter().map(|v| v.norm_sq()).sum::<f64>();
        assert!((error / norm_sq).sqrt() < 3e-7);
        assert!(actual.scalar.quadrature_error_estimate_pa[receiver] >= 0.0);
        assert!(
            actual.quadrature_error_estimate_m_s[receiver]
                .iter()
                .all(|v| *v >= 0.0)
        );
    }
}

#[test]
fn complete_receiver_bank_resources_geometry_and_cancellation_are_admitted() {
    let wave = DampedWaveNumber {
        k: C64::new(4.2, 0.8),
        phase_origin_m: 0.0,
    };
    let points = [[0.23, 0.29, 0.1]];
    let bytes = required_row_bytes(1, 1).unwrap();
    let attempt =
        |locations: &[Point], options: Options, cap: usize, cancelled: &mut dyn FnMut() -> bool| {
            OpenRows::new(
                &[TRIANGLE],
                &[NORMAL],
                wave,
                423.7,
                locations,
                1e-7,
                options,
                cap,
                &mut || cancelled(),
            )
        };
    assert!(attempt(&points, Options::default(), bytes - 1, &mut || false).is_err());
    assert!(attempt(&[[0.2, 0.2, 0.0]], Options::default(), bytes, &mut || false).is_err());
    assert!(
        attempt(
            &points,
            Options {
                maximum_depth: 0,
                ..Options::default()
            },
            bytes,
            &mut || false
        )
        .is_err()
    );
    assert!(
        attempt(
            &[points[0], points[0]],
            Options {
                maximum_kernel_evaluations: 80,
                ..Options::default()
            },
            2 * bytes,
            &mut || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        attempt(&points, Options::default(), bytes, &mut || {
            polls += 1;
            polls == 6
        }),
        Err(HelmholtzError::Cancelled)
    ));
    let empty = attempt(
        &[],
        Options::default(),
        required_row_bytes(1, 0).unwrap(),
        &mut || false,
    )
    .unwrap();
    assert!(
        empty
            .evaluate(|_| JUMP, &mut || false)
            .unwrap()
            .scalar
            .pressure
            .is_empty()
    );
    assert!(
        empty
            .evaluate(|_| [C64::new(f64::NAN, 0.0); 3], &mut || false)
            .is_err()
    );
    assert!(matches!(
        empty.evaluate(|_| JUMP, &mut || true),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(required_row_bytes(usize::MAX, 2).is_none());
}
