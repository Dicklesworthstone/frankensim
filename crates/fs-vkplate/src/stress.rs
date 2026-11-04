//! Precompute the native Airy solve in an explicitly supplied displacement space.
//!
//! Retained, unreferenced drum experiment: the fine-grid 32/64-mode screen
//! passed potential/gradient parity but lost paired kernel cost. Those spaces
//! also omitted most high-band source modes. This is not a production fast path.
//!
//! This compiles the actual grid, rest curvature and Airy boundary operator. It
//! does not select modes, approximate the stress space, or install a time stepper.
//! Contact, damping, mass and retarded pressure ports remain the caller's work;
//! potential/gradient agreement alone is not a sounded-bandwidth certificate.

use super::{Body, dot};

/// Sum-of-squares representation of the native shell potential at `w = Phi a`.
///
/// For features `y = [a_i, a_i a_j (i <= j)]`, native compatibility is `q = S y`.
/// With the native factor `A = L L'`, QR of `L^-1 S` gives `R`, so
/// `V = h^2/(2 E t) |R y|^2`. Linear features retain the curved rest tangent.
/// No column or stress channel is discarded by a rank or magnitude threshold.
/// Evaluation reuses scratch storage and makes no allocations.
#[derive(Clone, Debug)]
pub struct ShellStress {
    coordinates: usize,
    rows: usize,
    scale: f64,
    r: Vec<f64>,
    pairs: Vec<(usize, usize)>,
    features: Vec<f64>,
    feature_force: Vec<f64>,
    gradient: Vec<f64>,
}

impl ShellStress {
    /// Compile finite, caller-authored displacement shapes in the body's node order.
    /// Shapes need not be independent or mass-orthonormal for this energy pullback.
    /// This snapshot is invalidated by changing the source operators or rest shape.
    ///
    /// Refuses non-shell bodies, malformed shapes, nonfinite arithmetic and more
    /// than 4,194,304 coefficient entries (32 MiB per dense work array). That is a
    /// resource limit, not permission to drop modes until a model fits.
    pub fn compile(body: &Body, basis: &[Vec<f64>]) -> Result<Self, String> {
        let vk = body
            .vk
            .as_ref()
            .ok_or("stress compilation requires a native von Karman body")?;
        let m = basis.len();
        let p = m
            .checked_add(3)
            .and_then(|k| m.checked_mul(k))
            .map(|k| k / 2)
            .ok_or("stress feature count overflow")?;
        let n = vk.nb.len();
        let rows = n.min(p);
        const MAX_ENTRIES: usize = 4_194_304;
        if m == 0
            || n == 0
            || basis
                .iter()
                .any(|v| v.len() != body.n || v.iter().any(|x| !x.is_finite()))
        {
            return Err(
                "stress compilation requires nonempty finite shapes in the source node order"
                    .into(),
            );
        }
        if n.checked_mul(p).is_none_or(|k| k > MAX_ENTRIES)
            || rows.checked_mul(p).is_none_or(|k| k > MAX_ENTRIES)
        {
            return Err(format!(
                "native stress compilation needs {n} x {p} coefficients; limit {MAX_ENTRIES}; no modes were dropped"
            ));
        }
        let h2 = body.h * body.h;
        let scale = h2 / vk.eh;
        if !scale.is_finite() || scale <= 0.0 {
            return Err("native shell scale must be finite and positive".into());
        }
        let curvature: Vec<Vec<[f64; 3]>> = basis
            .iter()
            .map(|shape| {
                let at = |node: Option<usize>| node.map_or(0.0, |i| shape[i]);
                vk.nb
                    .iter()
                    .map(|s| {
                        [
                            (at(s[1]) - 2.0 * at(s[0]) + at(s[2])) / h2,
                            (at(s[3]) - 2.0 * at(s[0]) + at(s[4])) / h2,
                            (at(s[5]) - at(s[6]) - at(s[7]) + at(s[8])) / (4.0 * h2),
                        ]
                    })
                    .collect()
            })
            .collect();
        let pairs: Vec<_> = (0..m).flat_map(|i| (i..m).map(move |j| (i, j))).collect();
        // Column major: each compatibility column traverses the unchanged
        // native forward factor. The incumbent fused RHS+solve is untouched.
        let mut z = vec![0.0; n * p];
        for feature in 0..p {
            let column = &mut z[feature * n..(feature + 1) * n];
            for c in 0..n {
                let value = if feature < m {
                    let [xx, yy, xy] = curvature[feature][c];
                    let (zx, zy, zxy) = vk.rest_curv[c];
                    xx * zy + yy * zx - 2.0 * xy * zxy
                } else {
                    let (i, j) = pairs[feature - m];
                    let [xx, yy, xy] = curvature[i][c];
                    let [jx, jy, jxy] = curvature[j][c];
                    if i == j {
                        xx * yy - xy * xy
                    } else {
                        xx * jy + jx * yy - 2.0 * xy * jxy
                    }
                };
                column[c] = -vk.eh * value;
            }
            let band = &vk.airy;
            for i in 0..n {
                let lo = i.saturating_sub(band.b);
                let start = i * (band.b + 1) + lo + band.b - i;
                let row = &band.l[start..i * (band.b + 1) + band.b + 1];
                column[i] = (column[i] - dot(&row[..i - lo], &column[lo..i])) / row[i - lo];
            }
        }
        if z.iter().any(|x| !x.is_finite()) {
            return Err("nonfinite native compatibility coefficients".into());
        }
        // Thin Householder QR preserves the whole column image, including
        // dependent or exactly zero columns. No Gram matrix or eigenvalue
        // clipping can introduce negative storage or erase the rest tangent.
        let mut reflector = vec![0.0; n];
        for k in 0..rows {
            let column = &z[k * n + k..(k + 1) * n];
            let max = column.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
            if max == 0.0 {
                continue;
            }
            let v = &mut reflector[..n - k];
            for (v, x) in v.iter_mut().zip(column) {
                *v = *x / max;
            }
            let norm = dot(v, v).sqrt();
            let sign = if v[0] >= 0.0 { 1.0 } else { -1.0 };
            let diagonal = -sign * max * norm;
            v[0] += sign * norm;
            let beta = 2.0 / dot(v, v);
            for j in k + 1..p {
                let column = &mut z[j * n + k..(j + 1) * n];
                let projection = beta * dot(v, column);
                for (x, v) in column.iter_mut().zip(v.iter()) {
                    *x -= projection * v;
                }
            }
            z[k * n + k] = diagonal;
            z[k * n + k + 1..(k + 1) * n].fill(0.0);
        }
        let mut r = vec![0.0; rows * p];
        for i in 0..rows {
            for j in i..p {
                r[i * p + j] = z[j * n + i];
            }
        }
        if r.iter().any(|x| !x.is_finite()) {
            return Err("nonfinite native stress factor".into());
        }
        Ok(Self {
            coordinates: m,
            rows,
            scale,
            r,
            pairs,
            features: vec![0.0; p],
            feature_force: vec![0.0; p],
            gradient: vec![0.0; m],
        })
    }

    pub fn coordinates(&self) -> usize {
        self.coordinates
    }

    /// Number of scalar features; useful for reporting the actual evaluation cost.
    pub fn features(&self) -> usize {
        self.features.len()
    }

    /// Number of retained sum-of-squares rows. This is not a selected modal rank.
    pub fn stress_rows(&self) -> usize {
        self.rows
    }

    /// Evaluate native nonlinear potential and its coordinate gradient.
    /// Floating-point association differs from the nodal solver; compare with
    /// a declared tolerance. A refusal leaves the caller's gradient unchanged.
    pub fn evaluate(&mut self, a: &[f64], gradient: &mut [f64]) -> Result<f64, String> {
        let m = self.coordinates;
        let p = self.features.len();
        if a.len() != m || gradient.len() != m || a.iter().any(|x| !x.is_finite()) {
            return Err(
                "stress evaluation requires finite coordinates and an equally sized gradient"
                    .into(),
            );
        }
        self.features[..m].copy_from_slice(a);
        for (k, &(i, j)) in self.pairs.iter().enumerate() {
            self.features[m + k] = a[i] * a[j];
        }
        if self.features.iter().any(|x| !x.is_finite()) {
            return Err("stress features overflowed".into());
        }
        self.feature_force.fill(0.0);
        let mut energy = 0.0;
        for i in 0..self.rows {
            let row = &self.r[i * p + i..(i + 1) * p];
            let u = dot(row, &self.features[i..]);
            energy += u * u;
            for (force, coefficient) in self.feature_force[i..].iter_mut().zip(row) {
                *force += coefficient * u;
            }
        }
        self.gradient.copy_from_slice(&self.feature_force[..m]);
        for (k, &(i, j)) in self.pairs.iter().enumerate() {
            let force = self.feature_force[m + k];
            self.gradient[i] += a[j] * force;
            self.gradient[j] += a[i] * force;
        }
        let energy = 0.5 * self.scale * energy;
        for g in &mut self.gradient {
            *g *= self.scale;
        }
        if !energy.is_finite() || self.gradient.iter().any(|x| !x.is_finite()) {
            return Err("native stress evaluation overflowed".into());
        }
        gradient.copy_from_slice(&self.gradient);
        Ok(energy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Edge, Nonlinear, PlateSpec};

    fn shell(curved: bool) -> Body {
        Body::plate(
            "stress source",
            &PlateSpec {
                radius: 3.5,
                hole: 0.0,
                grid: 1.0,
                edge: Edge::Free,
                thickness: &|_| 1.0 / 128.0,
                rest_height: &|r| {
                    if !curved || r < 0.5 {
                        0.0
                    } else if r < 1.1 {
                        0.5
                    } else {
                        1.0
                    }
                },
                young: 128.0,
                poisson: 0.3,
                density: 1000.0,
                tension: 0.0,
                bending: true,
                nonlinear: Nonlinear::VonKarman,
                t60: [(100.0, 1e12), (1000.0, 1e12)],
                support: None,
                muffle: None,
            },
        )
        .unwrap()
    }

    fn compare(body: &mut Body, basis: &[Vec<f64>], a: &[f64], compiled: &mut ShellStress) -> f64 {
        for k in 0..body.n {
            body.w[k] = basis.iter().zip(a).map(|(v, a)| v[k] * a).sum();
        }
        let expected = body.nonlinear_gradient();
        let gradient: Vec<_> = basis.iter().map(|v| dot(v, &body.g)).collect();
        let mut actual = vec![0.0; a.len()];
        let energy = compiled.evaluate(a, &mut actual).unwrap();
        let gnorm = dot(&gradient, &gradient).sqrt().max(1e-300);
        let error = actual
            .iter()
            .zip(&gradient)
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(
            (expected - energy).abs() <= 1e-11 * expected.abs().max(1e-300),
            "potential {expected:e} vs {energy:e}"
        );
        assert!(error <= 1e-11 * gnorm, "gradient {error:e} / {gnorm:e}");
        energy
    }

    #[test]
    fn complete_source_basis_preserves_flat_and_curved_energy_and_gradient() {
        for curved in [false, true] {
            let mut body = shell(curved);
            let basis: Vec<Vec<f64>> = (0..body.n)
                .map(|i| {
                    (0..body.n)
                        .map(|j| if i == j { 1.0 } else { 0.0 })
                        .collect()
                })
                .collect();
            let mut compiled = ShellStress::compile(&body, &basis).unwrap();
            assert_eq!(
                compare(&mut body, &basis, &vec![0.0; basis.len()], &mut compiled),
                0.0
            );
            for amplitude in [0.2, -0.01, 2.0_f64.powi(-55), -2.0_f64.powi(-55)] {
                let a: Vec<_> = (0..body.n)
                    .map(|i| amplitude * ((i * 17 % 23) as f64 - 11.0))
                    .collect();
                assert!(compare(&mut body, &basis, &a, &mut compiled) > 0.0);
            }
        }
    }

    #[test]
    fn mixed_dependent_and_zero_shapes_keep_the_native_pullback() {
        let mut body = shell(true);
        let a: Vec<_> = body.xy.iter().map(|p| p[0] * p[0] + p[1]).collect();
        let b: Vec<_> = body.xy.iter().map(|p| p[0] * p[1]).collect();
        let basis = vec![a.clone(), b, a, vec![0.0; body.n]];
        let mut compiled = ShellStress::compile(&body, &basis).unwrap();
        for amplitudes in [[0.02, -0.07, 0.01, 4.0], [-0.1, 0.02, 0.1, 0.0]] {
            assert!(compare(&mut body, &basis, &amplitudes, &mut compiled) > 0.0);
        }
        let mut zero = ShellStress::compile(&body, &[vec![0.0; body.n]; 1]).unwrap();
        let mut gradient = [1.0];
        assert_eq!(zero.evaluate(&[4.0], &mut gradient).unwrap(), 0.0);
        assert_eq!(gradient, [0.0]);
    }

    #[test]
    fn invalid_sources_and_nonfinite_inputs_refuse_without_overwriting_output() {
        assert!(ShellStress::compile(&Body::point("point", 1.0), &[vec![1.0]]).is_err());
        let body = shell(true);
        for basis in [
            vec![],
            vec![vec![0.0]],
            vec![vec![f64::NAN; body.n]],
            // Even one native core needs more than MAX_ENTRIES features.
            vec![vec![0.0; body.n]; 3000],
        ] {
            assert!(ShellStress::compile(&body, &basis).is_err());
        }
        let mut compiled = ShellStress::compile(&body, &[vec![1.0; body.n]]).unwrap();
        let mut gradient = [73.0];
        for a in [vec![], vec![f64::INFINITY], vec![f64::MAX]] {
            assert!(compiled.evaluate(&a, &mut gradient).is_err());
            assert_eq!(gradient, [73.0]);
        }
        assert!(compiled.evaluate(&[0.0], &mut []).is_err());
        assert_eq!(compiled.evaluate(&[0.0], &mut gradient).unwrap(), 0.0);
        assert_eq!(gradient, [0.0]);

        // Finite monomials can still overflow the squared stress or force.
        // A nonzero source shape reaches that later refusal independently.
        let mut bump = vec![0.0; body.n];
        bump[body.node_at(0.0, 0.0).unwrap()] = 1.0;
        let mut compiled = ShellStress::compile(&body, &[bump]).unwrap();
        gradient[0] = 73.0;
        assert!(compiled.evaluate(&[1e100], &mut gradient).is_err());
        assert_eq!(gradient, [73.0]);
        assert_eq!(compiled.evaluate(&[0.0], &mut gradient).unwrap(), 0.0);
        assert_eq!(gradient, [0.0]);
    }
}
