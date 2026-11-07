//! Exact radial transmission problem for bonded isotropic annuli with common nu.
//!
//! With x=r/R, axial warping is R^3 (x^3/4 + a*x + b/x) cos(theta)
//! per unit curvature gradient. Common nu makes the in-plane Poisson field
//! continuous. Displacement and G*(v'+nu*x^2/2) are continuous at interfaces;
//! the latter vanishes on the free bore/exterior. The exact annular harmonic
//! Dirichlet-to-Neumann map is a symmetric tridiagonal system. It uses the
//! shared adjacent-pivot solve; there is no radial mesh or dense section solve.

use fs_la::tridiagonal::TridiagonalWorkspace;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct LayerCoefficients {
    pub a: f64,
    pub b: f64,
    pub stress_scale: f64,
    pub peak: f64,
}

impl LayerCoefficients {
    /// Radial/tangential stress eigenvalues multiplying Q/A at r/R.
    pub fn factors(self, radius: f64, nu: f64) -> [f64; 2] {
        let s = radius * radius;
        let inverse = if radius == 0.0 {
            0.0
        } else {
            (self.b / radius) / radius
        };
        [
            self.stress_scale * ((0.75 + 0.5 * nu) * s + self.a - inverse),
            self.stress_scale * ((0.25 - 0.5 * nu) * s + self.a + inverse),
        ]
    }
}

pub(super) fn prepare(radii: &[f64], g: &[f64], nu: f64) -> Option<Vec<LayerCoefficients>> {
    let n = g.len();
    if n == 0
        || radii.len() != n + 1
        || radii[n] != 1.0
        || !nu.is_finite()
        || nu <= -1.0
        || nu >= 0.5
        || g.iter().any(|x| !x.is_finite() || *x <= 0.0 || *x > 1.0)
        || radii
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0 || (*x > 0.0 && x * x == 0.0))
        || radii.windows(2).any(|x| x[0] >= x[1])
    {
        return None;
    }
    let first = usize::from(radii[0] == 0.0);
    let dimension = n + 1 - first;
    let mut diagonal = vec![0.0; dimension];
    let mut off = vec![0.0; dimension - 1];
    let mut rhs = vec![0.0; dimension];
    let mut moment = 0.0;
    for i in 0..n {
        let (a, b) = (radii[i], radii[i + 1]);
        let d = (b - a) * (b + a);
        moment += g[i] * d * (b * b + a * a);
        let right = i + 1 - first;
        if a == 0.0 {
            diagonal[right] += g[i];
            rhs[right] -= 0.5 * g[i] * (1.0 + nu) * b.powi(3);
        } else {
            let left = right - 1;
            let diag = g[i] * ((a * a + b * b) / d);
            diagonal[left] += diag;
            diagonal[right] += diag;
            off[left] = -2.0 * g[i] * (a * b / d);
            rhs[left] += 0.5 * g[i] * a * ((1.0 + nu) * a * a - b * b);
            rhs[right] += 0.5 * g[i] * b * (a * a - (1.0 + nu) * b * b);
        }
    }
    if !moment.is_finite() || moment <= 0.0 {
        return None;
    }
    let mut v = vec![0.0; dimension];
    TridiagonalWorkspace::new(dimension)
        .ok()?
        .solve_into(&off, &diagonal, &off, &rhs, &mut v, || false)
        .ok()?;
    // Check the original equations, including both free contours and every
    // interface. A poorly resolved stiffness contrast is refused, never clipped.
    for i in 0..dimension {
        let center = diagonal[i] * v[i];
        let left = if i == 0 { 0.0 } else { off[i - 1] * v[i - 1] };
        let right = if i + 1 == dimension {
            0.0
        } else {
            off[i] * v[i + 1]
        };
        let residual = (center + left + right - rhs[i]).abs();
        let scale = center.abs() + left.abs() + right.abs() + rhs[i].abs();
        if !residual.is_finite()
            || !scale.is_finite()
            || residual > 128.0 * f64::EPSILON * dimension as f64 * scale
        {
            return None;
        }
    }
    let area = (1.0 - radii[0]) * (1.0 + radii[0]);
    let mut layers = Vec::with_capacity(n);
    for i in 0..n {
        let (inner, outer) = (radii[i], radii[i + 1]);
        let z_outer = v[i + 1 - first] - outer.powi(3) / 4.0;
        let (a, b) = if inner == 0.0 {
            (z_outer / outer, 0.0)
        } else {
            let z_inner = v[i - first] - inner.powi(3) / 4.0;
            let d = (outer - inner) * (outer + inner);
            (
                (outer * z_outer - inner * z_inner) / d,
                inner * outer * ((outer * z_inner - inner * z_outer) / d),
            )
        };
        let stress_scale = -(2.0 * g[i]) * (area / ((1.0 + nu) * moment));
        let mut layer = LayerCoefficients {
            a,
            b,
            stress_scale,
            peak: 0.0,
        };
        for radius in [inner, outer] {
            for value in layer.factors(radius, nu) {
                layer.peak = layer.peak.max(value.abs());
            }
        }
        // For c*s+a+d/s, the sole interior extremum is s=sqrt(d/c).
        // Include both eigenvalue extrema, even when an auxetic layer reverses
        // a component. This bounds the entire physical annulus, not one point.
        for (c, d) in [(0.75 + 0.5 * nu, -b), (0.25 - 0.5 * nu, b)] {
            if d > 0.0 {
                let radius = (d / c).sqrt().sqrt();
                if radius > inner && radius < outer {
                    for value in layer.factors(radius, nu) {
                        layer.peak = layer.peak.max(value.abs());
                    }
                }
            }
        }
        if [a, b, stress_scale, layer.peak]
            .iter()
            .any(|x| !x.is_finite())
            || layer.peak <= 0.0
            || stress_scale == 0.0
        {
            return None;
        }
        layers.push(layer);
    }
    Some(layers)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64) {
        assert!(
            (a - b).abs() < 3e-11 * a.abs().max(b.abs()).max(1.0),
            "{a} != {b}"
        );
    }

    #[test]
    fn homogeneous_solid_and_tube_recover_the_closed_form_without_a_radial_mesh() {
        for bore in [0.0, 0.2, 0.8] {
            for nu in [-0.75, 0.0, 0.25, 0.49] {
                let radii = [
                    bore,
                    bore + 0.3 * (1.0 - bore),
                    bore + 0.7 * (1.0 - bore),
                    1.0,
                ];
                let layers = prepare(&radii, &[1.0; 3], nu).unwrap();
                let c = 0.5 / ((1.0 + nu) * (1.0 + bore * bore));
                for (i, layer) in layers.iter().enumerate() {
                    for t in [0.0, 0.3, 0.8, 1.0] {
                        let r = radii[i] + t * (radii[i + 1] - radii[i]);
                        let q = if r == 0.0 { 0.0 } else { (bore / r).powi(2) };
                        let expected = [
                            c * (3.0 + 2.0 * nu) * (1.0 - r * r) * (1.0 - q),
                            c * ((3.0 + 2.0 * nu) * (1.0 + bore * bore + q)
                                - (1.0 - 2.0 * nu) * r * r),
                        ];
                        let actual = layer.factors(r, nu);
                        near(actual[0], expected[0]);
                        near(actual[1], expected[1]);
                    }
                }
            }
        }
    }

    #[test]
    fn heterogeneous_layers_match_displacement_traction_and_original_resultant() {
        let radii = [0.12, 0.3, 0.72, 1.0];
        let g = [0.02, 1.0, 0.15];
        for nu in [-0.7, 0.0, 0.3, 0.49] {
            let layers = prepare(&radii, &g, nu).unwrap();
            near(layers[0].factors(radii[0], nu)[0], 0.0);
            near(layers[2].factors(1.0, nu)[0], 0.0);
            for i in 1..3 {
                let r = radii[i];
                let l = layers[i - 1];
                let h = layers[i];
                near(l.a * r + l.b / r, h.a * r + h.b / r);
                near(l.factors(r, nu)[0], h.factors(r, nu)[0]);
                near(l.factors(r, nu)[1] / g[i - 1], h.factors(r, nu)[1] / g[i]);
            }
            // Exact independent polar resultant integral. The Laurent terms
            // cancel in (radial+tangential)*r, leaving a cubic polynomial.
            let mut resultant = 0.0;
            for i in 0..3 {
                let l = layers[i];
                let (a, b) = (radii[i], radii[i + 1]);
                resultant +=
                    l.stress_scale * ((b.powi(4) - a.powi(4)) / 4.0 + l.a * (b * b - a * a));
                for j in 0..101 {
                    let r = a + (b - a) * j as f64 / 100.0;
                    assert!(
                        l.factors(r, nu)
                            .iter()
                            .all(|x| x.abs() <= l.peak * (1.0 + 1e-13))
                    );
                }
            }
            near(resultant / (1.0 - radii[0] * radii[0]), 1.0);
            // Replacing an original layer by identical sublayers changes no field.
            let split = prepare(&[0.12, 0.3, 0.5, 0.72, 1.0], &[0.02, 1.0, 1.0, 0.15], nu).unwrap();
            for (i, k, r) in [(0, 0, 0.2), (1, 1, 0.4), (1, 2, 0.6), (2, 3, 0.8)] {
                for (a, b) in layers[i]
                    .factors(r, nu)
                    .into_iter()
                    .zip(split[k].factors(r, nu))
                {
                    near(a, b);
                }
            }
        }
    }
}
