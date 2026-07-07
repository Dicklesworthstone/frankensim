//! Free-form deformation (plan §7.6): a trivariate Bernstein control
//! lattice warps any chart's ambient space (Sederberg–Parry FFD). The
//! warp is LINEAR in the control displacements, so the Jacobian action is
//! the exact basis contraction — no approximation anywhere.
//!
//! Points outside the lattice box pass through unchanged (classic FFD
//! semantics). Continuity across the box boundary is the CALLER's
//! control-point responsibility (keep boundary-layer controls at zero for
//! C⁰ blends) — documented, not silently smoothed.

use crate::{Parameterization, XformError};
use fs_geom::{Point3, Vec3};

/// A trivariate Bernstein FFD lattice over an axis-aligned box.
/// θ layout: `3·nx·ny·nz` displacements, node-major
/// (`node = (i·ny + j)·nz + k`, then x/y/z).
#[derive(Debug, Clone)]
pub struct FfdLattice {
    /// Box minimum corner.
    pub origin: Point3,
    /// Box edge lengths (all > 0).
    pub size: Vec3,
    /// Control counts per axis (each ≥ 2).
    pub counts: [usize; 3],
}

impl FfdLattice {
    /// Local coordinates in [0,1]³, or `None` outside the box.
    fn local(&self, x: Point3) -> Option<[f64; 3]> {
        let d = x.delta_from(self.origin);
        let t = [d.x / self.size.x, d.y / self.size.y, d.z / self.size.z];
        if t.iter().all(|&v| (0.0..=1.0).contains(&v)) {
            Some(t)
        } else {
            None
        }
    }

    fn check_theta(&self, theta: &[f64]) -> Result<(), XformError> {
        let expected = self.dof();
        if theta.len() != expected {
            return Err(XformError::DofMismatch {
                expected,
                got: theta.len(),
            });
        }
        Ok(())
    }

    /// The displacement field D(x) = Σ B(t)·θ and (optionally) its
    /// spatial gradient rows.
    fn displacement(&self, theta: &[f64], t: [f64; 3]) -> (Vec3, [[f64; 3]; 3]) {
        let (nx, ny, nz) = (self.counts[0], self.counts[1], self.counts[2]);
        let bx = bernstein(nx - 1, t[0]);
        let by = bernstein(ny - 1, t[1]);
        let bz = bernstein(nz - 1, t[2]);
        let dbx = bernstein_derivative(nx - 1, t[0]);
        let dby = bernstein_derivative(ny - 1, t[1]);
        let dbz = bernstein_derivative(nz - 1, t[2]);
        let mut disp = [0.0f64; 3];
        let mut grad = [[0.0f64; 3]; 3]; // grad[component][axis]
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let w = bx[i] * by[j] * bz[k];
                    let wg = [
                        dbx[i] * by[j] * bz[k] / self.size.x,
                        bx[i] * dby[j] * bz[k] / self.size.y,
                        bx[i] * by[j] * dbz[k] / self.size.z,
                    ];
                    let node = (i * ny + j) * nz + k;
                    for c in 0..3 {
                        let p = theta[3 * node + c];
                        disp[c] += w * p;
                        for a in 0..3 {
                            grad[c][a] += wg[a] * p;
                        }
                    }
                }
            }
        }
        (Vec3::new(disp[0], disp[1], disp[2]), grad)
    }
}

impl Parameterization for FfdLattice {
    fn dof(&self) -> usize {
        3 * self.counts[0] * self.counts[1] * self.counts[2]
    }

    fn apply(&self, theta: &[f64], x: Point3) -> Result<Point3, XformError> {
        self.check_theta(theta)?;
        match self.local(x) {
            None => Ok(x),
            Some(t) => {
                let (d, _) = self.displacement(theta, t);
                Ok(x.offset(d))
            }
        }
    }

    fn jacobian_action(
        &self,
        theta: &[f64],
        dtheta: &[f64],
        x: Point3,
    ) -> Result<Vec3, XformError> {
        self.check_theta(theta)?;
        self.check_theta(dtheta)?;
        match self.local(x) {
            None => Ok(Vec3::new(0.0, 0.0, 0.0)),
            Some(t) => {
                // Linear in θ: the action is the basis contraction with δθ.
                let (d, _) = self.displacement(dtheta, t);
                Ok(d)
            }
        }
    }

    fn spatial_jacobian(&self, theta: &[f64], x: Point3) -> Result<[[f64; 3]; 3], XformError> {
        self.check_theta(theta)?;
        let mut jac = [[0.0f64; 3]; 3];
        for (r, row) in jac.iter_mut().enumerate() {
            row[r] = 1.0; // identity
        }
        if let Some(t) = self.local(x) {
            let (_, grad) = self.displacement(theta, t);
            for c in 0..3 {
                for a in 0..3 {
                    jac[c][a] += grad[c][a];
                }
            }
        }
        Ok(jac)
    }
}

/// Bernstein basis values `B_{i,n}(t)` for i = 0..=n (iterative de
/// Casteljau-style construction; no factorials, no overflow).
#[must_use]
pub fn bernstein(n: usize, t: f64) -> Vec<f64> {
    let mut b = vec![0.0f64; n + 1];
    b[0] = 1.0;
    for d in 1..=n {
        // Build degree d from degree d−1, in place, right to left.
        b[d] = t * b[d - 1];
        for i in (1..d).rev() {
            b[i] = t * b[i - 1] + (1.0 - t) * b[i];
        }
        b[0] *= 1.0 - t;
    }
    b
}

/// Bernstein basis derivatives `dB_{i,n}/dt = n·(B_{i−1,n−1} − B_{i,n−1})`.
#[must_use]
pub fn bernstein_derivative(n: usize, t: f64) -> Vec<f64> {
    if n == 0 {
        return vec![0.0];
    }
    let lower = bernstein(n - 1, t);
    let nf = n as f64;
    (0..=n)
        .map(|i| {
            let left = if i > 0 { lower[i - 1] } else { 0.0 };
            let right = if i <= n - 1 { lower[i] } else { 0.0 };
            nf * (left - right)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bernstein_partitions_unity_and_derivative_sums_to_zero() {
        for &t in &[0.0, 0.25, 0.5, 0.9, 1.0] {
            for n in 0..6 {
                let b = bernstein(n, t);
                let sum: f64 = b.iter().sum();
                assert!(
                    (sum - 1.0).abs() < 1e-12,
                    "partition of unity at n={n}, t={t}"
                );
                let d = bernstein_derivative(n, t);
                let dsum: f64 = d.iter().sum();
                assert!(dsum.abs() < 1e-12, "derivative sum at n={n}, t={t}");
            }
        }
    }

    #[test]
    fn identity_lattice_is_the_identity_map() {
        let ffd = FfdLattice {
            origin: Point3::new(0.0, 0.0, 0.0),
            size: Vec3::new(1.0, 1.0, 1.0),
            counts: [3, 3, 3],
        };
        let theta = vec![0.0; ffd.dof()];
        let x = Point3::new(0.3, 0.7, 0.5);
        let y = ffd.apply(&theta, x).unwrap();
        for (a, b) in [(y.x, x.x), (y.y, x.y), (y.z, x.z)] {
            assert!((a - b).abs() < 1e-15, "identity drifted: {a} vs {b}");
        }
        let j = ffd.spatial_jacobian(&theta, x).unwrap();
        assert!((crate::det3(&j) - 1.0).abs() < 1e-12);
        // Outside the box: pass-through, zero velocity.
        let out = Point3::new(5.0, 5.0, 5.0);
        assert!((ffd.apply(&theta, out).unwrap().x - 5.0).abs() < 1e-15);
        let mut dtheta = vec![0.0; ffd.dof()];
        dtheta[0] = 1.0;
        let v = ffd.jacobian_action(&theta, &dtheta, out).unwrap();
        assert!(
            v.x.abs() + v.y.abs() + v.z.abs() == 0.0_f64.abs(),
            "outside is exactly zero"
        );
    }
}
