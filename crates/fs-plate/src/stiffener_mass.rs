//! Translational inertia of the existing two-node Hermite stiffener field.
//! The coordinate order is `(w1, dw/ds1, w2, dw/ds2)`, with both slopes
//! measured along the same directed segment. This is the exact integral of
//! `rho*A*w_dot(s)^2`, not an added rotary or eccentric axial inertia term.

/// Transverse translational mass law for the existing Euler–Bernoulli beams.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum StiffenerMass {
    /// Half the segment mass at each endpoint displacement; the original law.
    #[default]
    Lumped,
    /// Exact mass integral of the cubic Hermite displacement used by bending.
    ConsistentHermite,
}

/// Exact `integral rho*A*N(s)^T*N(s) ds`. Zero line density admits a
/// deliberately massless brace, as in the existing stiffness-isolation tests.
pub(crate) fn hermite_mass(length: f64, line_density: f64) -> Option<[[f64; 4]; 4]> {
    if !length.is_finite() || length <= 0.0 || !line_density.is_finite() || line_density < 0.0 {
        return None;
    }
    let l = length;
    let l2 = l * l;
    let scale = line_density * l / 420.0;
    let coefficients = [
        [156.0, 22.0 * l, 54.0, -13.0 * l],
        [22.0 * l, 4.0 * l2, 13.0 * l, -3.0 * l2],
        [54.0, 13.0 * l, 156.0, -22.0 * l],
        [-13.0 * l, -3.0 * l2, -22.0 * l, 4.0 * l2],
    ];
    let mass = coefficients.map(|row| row.map(|value| scale * value));
    mass.iter()
        .flatten()
        .all(|value| value.is_finite())
        .then_some(mass)
}

#[cfg(test)]
mod tests {
    use super::hermite_mass;

    fn bilinear(m: &[[f64; 4]; 4], a: [f64; 4], b: [f64; 4]) -> f64 {
        (0..4)
            .flat_map(|i| (0..4).map(move |j| a[i] * m[i][j] * b[j]))
            .sum()
    }

    // Independent polynomial endpoints, slopes and continuum product integral.
    // These do not evaluate the element shape functions or repeat its matrix.
    fn endpoints(p: [f64; 4], length: f64) -> [f64; 4] {
        let value = p[0] + length * (p[1] + length * (p[2] + length * p[3]));
        let slope = p[1] + length * (2.0 * p[2] + 3.0 * length * p[3]);
        [p[0], p[1], value, slope]
    }

    fn integral(a: [f64; 4], b: [f64; 4], length: f64) -> f64 {
        (0..4)
            .flat_map(|i| {
                (0..4).map(move |j| {
                    a[i] * b[j] * length.powi((i + j + 1) as i32) / (i + j + 1) as f64
                })
            })
            .sum()
    }

    #[test]
    fn consistent_hermite_mass_integrates_every_cubic_velocity_product() {
        let fields = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.7, -1.1, 0.4, 0.9],
        ];
        for length in [0.03, 0.7, 2.0] {
            let density = 3.7;
            let m = hermite_mass(length, density).unwrap();
            for a in fields {
                for b in fields {
                    let actual = bilinear(&m, endpoints(a, length), endpoints(b, length));
                    let expected = density * integral(a, b, length);
                    assert!((actual - expected).abs() < 3e-13 * expected.abs().max(1.0));
                }
                assert!(bilinear(&m, endpoints(a, length), endpoints(a, length)) > 0.0);
            }
            assert!(
                (bilinear(&m, [1.0, 0.0, 1.0, 0.0], [1.0, 0.0, 1.0, 0.0]) - density * length).abs()
                    < 2e-14
            );
        }
    }

    #[test]
    fn consistent_hermite_mass_preserves_energy_when_beam_direction_is_reversed() {
        let q = [0.3, 0.8, -1.2, 0.4];
        let reversed = [q[2], -q[3], q[0], -q[1]];
        let m = hermite_mass(0.73, 4.2).unwrap();
        assert!((bilinear(&m, q, q) - bilinear(&m, reversed, reversed)).abs() < 1e-14);
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(m[i][j], m[j][i]);
            }
        }
    }

    #[test]
    fn consistent_hermite_mass_preserves_polynomial_energy_under_subdivision() {
        let p = [0.4, -0.2, 0.7, -0.1];
        let length = 1.7;
        let density = 2.3;
        let expected = density * integral(p, p, length);
        for segments in [1, 2, 7, 16] {
            let h = length / f64::from(segments);
            let m = hermite_mass(h, density).unwrap();
            let mut actual = 0.0;
            for segment in 0..segments {
                let x = f64::from(segment) * h;
                let shifted = [
                    p[0] + x * (p[1] + x * (p[2] + x * p[3])),
                    p[1] + x * (2.0 * p[2] + 3.0 * x * p[3]),
                    p[2] + 3.0 * x * p[3],
                    p[3],
                ];
                let q = endpoints(shifted, h);
                actual += bilinear(&m, q, q);
            }
            assert!((actual - expected).abs() < 2e-13 * expected);
        }
    }

    #[test]
    fn consistent_hermite_mass_refuses_nonphysical_and_overflowing_inputs() {
        for (length, density) in [
            (0.0, 1.0),
            (-1.0, 1.0),
            (1.0, -1.0),
            (f64::NAN, 1.0),
            (1.0, f64::INFINITY),
            (1e200, 1e200),
        ] {
            assert!(hermite_mass(length, density).is_none());
        }
        assert_eq!(hermite_mass(0.5, 0.0).unwrap(), [[0.0; 4]; 4]);
    }
}
