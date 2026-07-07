//! Level-set velocity parameterization (plan §7.6): θ = normal-velocity
//! DOFs on a grid, active inside a narrow band of the interface — THE
//! topology-optimization workhorse (Appendix C's
//! `xform.level-set-velocity vessel :band 12mm :dof 4096`), plus a
//! working first-order upwind advection step so the lever demonstrably
//! DRIVES a level set (the bead's acceptance criterion; WENO/octree
//! narrow bands are topo-levelset's).

use crate::XformError;
use fs_geom::Point3;

/// A scalar velocity field on a regular grid: trilinear in space, LINEAR
/// in θ (one DOF per node), masked to the narrow band `|φ(x)| ≤ band`.
#[derive(Debug, Clone)]
pub struct VelocityBand {
    /// Grid minimum corner.
    pub origin: Point3,
    /// Node spacing (> 0).
    pub spacing: f64,
    /// Node counts per axis (each ≥ 2).
    pub dims: [usize; 3],
    /// Narrow-band half-width in world units.
    pub band: f64,
}

impl VelocityBand {
    /// DOF count (one scalar normal velocity per node).
    #[must_use]
    pub fn dof(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    fn check_theta(&self, theta: &[f64]) -> Result<(), XformError> {
        if theta.len() != self.dof() {
            return Err(XformError::DofMismatch {
                expected: self.dof(),
                got: theta.len(),
            });
        }
        Ok(())
    }

    fn node(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.dims[1] + j) * self.dims[2] + k
    }

    /// Trilinear velocity at `x`, masked by the band around `phi_at_x`
    /// (pass the SDF value at the same point; outside the band the
    /// velocity — and every Jacobian entry — is exactly zero).
    ///
    /// # Errors
    /// [`XformError::DofMismatch`] on a wrong-length θ.
    pub fn velocity(&self, theta: &[f64], x: Point3, phi_at_x: f64) -> Result<f64, XformError> {
        self.check_theta(theta)?;
        if phi_at_x.abs() > self.band {
            return Ok(0.0);
        }
        let d = x.delta_from(self.origin);
        let g = [d.x / self.spacing, d.y / self.spacing, d.z / self.spacing];
        // Clamp into the grid (velocity extends constantly at the border).
        let cell = |v: f64, n: usize| -> (usize, f64) {
            let max = (n - 2) as f64;
            let c = v.clamp(0.0, max);
            let i = c.floor().min(max) as usize;
            (i, (c - i as f64).clamp(0.0, 1.0))
        };
        let (i, fx) = cell(g[0], self.dims[0]);
        let (j, fy) = cell(g[1], self.dims[1]);
        let (k, fz) = cell(g[2], self.dims[2]);
        let mut acc = 0.0;
        for (di, wi) in [(0usize, 1.0 - fx), (1, fx)] {
            for (dj, wj) in [(0usize, 1.0 - fy), (1, fy)] {
                for (dk, wk) in [(0usize, 1.0 - fz), (1, fz)] {
                    acc += wi * wj * wk * theta[self.node(i + di, j + dj, k + dk)];
                }
            }
        }
        Ok(acc)
    }

    /// The Jacobian action `∂v/∂θ · δθ` — identical basis contraction
    /// (linear lever), same band mask.
    ///
    /// # Errors
    /// [`XformError::DofMismatch`] on a wrong-length δθ.
    pub fn jacobian_action(
        &self,
        dtheta: &[f64],
        x: Point3,
        phi_at_x: f64,
    ) -> Result<f64, XformError> {
        self.velocity(dtheta, x, phi_at_x)
    }
}

/// One first-order Godunov upwind advection step of a node-sampled SDF
/// under the normal-speed field `v`: `φ ← φ − dt·v·|∇φ|` with upwinded
/// one-sided differences (Rouy–Tourin). Interior nodes only (the boundary
/// layer keeps its values; callers keep interfaces away from the grid
/// edge — documented).
pub fn advect_sdf(
    phi: &mut [f64],
    dims: [usize; 3],
    spacing: f64,
    v: &dyn Fn(usize, usize, usize) -> f64,
    dt: f64,
) {
    let node = |i: usize, j: usize, k: usize| (i * dims[1] + j) * dims[2] + k;
    let old = phi.to_vec();
    for i in 1..dims[0] - 1 {
        for j in 1..dims[1] - 1 {
            for k in 1..dims[2] - 1 {
                let c = old[node(i, j, k)];
                let dxm = (c - old[node(i - 1, j, k)]) / spacing;
                let dxp = (old[node(i + 1, j, k)] - c) / spacing;
                let dym = (c - old[node(i, j - 1, k)]) / spacing;
                let dyp = (old[node(i, j + 1, k)] - c) / spacing;
                let dzm = (c - old[node(i, j, k - 1)]) / spacing;
                let dzp = (old[node(i, j, k + 1)] - c) / spacing;
                let speed = v(i, j, k);
                // Godunov switches for outward (v>0) vs inward motion.
                let grad2 = if speed > 0.0 {
                    dxm.max(0.0).powi(2)
                        + dxp.min(0.0).powi(2)
                        + dym.max(0.0).powi(2)
                        + dyp.min(0.0).powi(2)
                        + dzm.max(0.0).powi(2)
                        + dzp.min(0.0).powi(2)
                } else {
                    dxp.max(0.0).powi(2)
                        + dxm.min(0.0).powi(2)
                        + dyp.max(0.0).powi(2)
                        + dym.min(0.0).powi(2)
                        + dzp.max(0.0).powi(2)
                        + dzm.min(0.0).powi(2)
                };
                phi[node(i, j, k)] = c - dt * speed * grad2.sqrt();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn band_mask_zeroes_velocity_and_action() {
        let band = VelocityBand {
            origin: Point3::new(0.0, 0.0, 0.0),
            spacing: 1.0,
            dims: [4, 4, 4],
            band: 0.5,
        };
        let theta = vec![2.0; band.dof()];
        let x = Point3::new(1.5, 1.5, 1.5);
        assert!(
            (band.velocity(&theta, x, 0.4).unwrap() - 2.0).abs() < 1e-15,
            "inside band"
        );
        assert_eq!(
            band.velocity(&theta, x, 0.9).unwrap().to_bits(),
            0u64,
            "outside band"
        );
        assert_eq!(
            band.jacobian_action(&theta, x, 0.9).unwrap().to_bits(),
            0u64
        );
        assert!(matches!(
            band.velocity(&[1.0], x, 0.0),
            Err(XformError::DofMismatch {
                expected: 64,
                got: 1
            })
        ));
    }
}
