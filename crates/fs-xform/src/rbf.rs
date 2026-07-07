//! Compactly supported RBF morphs (plan §7.6): scattered-handle
//! deformation with Wendland-C2 kernels — locality control by
//! construction (support radius), frame-equivariant because the kernel is
//! radial (the G3 law). Linear in θ: exact Jacobian actions.

use crate::{Parameterization, XformError};
use fs_geom::{Point3, Vec3};

/// Wendland C2 kernel `φ(q) = (1−q)⁴(4q+1)` on q ∈ [0,1], 0 beyond.
#[must_use]
pub fn wendland_c2(q: f64) -> f64 {
    if q >= 1.0 {
        return 0.0;
    }
    let omq = 1.0 - q;
    let omq2 = omq * omq;
    omq2 * omq2 * (4.0 * q + 1.0)
}

/// Derivative `φ'(q) = −20·q·(1−q)³` on q ∈ [0,1], 0 beyond.
#[must_use]
pub fn wendland_c2_derivative(q: f64) -> f64 {
    if q >= 1.0 {
        return 0.0;
    }
    let omq = 1.0 - q;
    -20.0 * q * omq * omq * omq
}

/// A scattered-handle morph: `T(x) = x + Σ φ(|x−cᵢ|/r)·θᵢ`.
/// θ layout: 3 displacements per handle, handle-major.
#[derive(Debug, Clone)]
pub struct RbfMorph {
    /// Handle centers.
    pub centers: Vec<Point3>,
    /// Common support radius (> 0).
    pub radius: f64,
}

impl RbfMorph {
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

    fn weighted_sum(&self, theta: &[f64], x: Point3) -> Vec3 {
        let mut acc = [0.0f64; 3];
        for (i, &c) in self.centers.iter().enumerate() {
            let q = x.delta_from(c).norm() / self.radius;
            let w = wendland_c2(q);
            if w != 0.0 {
                for a in 0..3 {
                    acc[a] += w * theta[3 * i + a];
                }
            }
        }
        Vec3::new(acc[0], acc[1], acc[2])
    }
}

impl Parameterization for RbfMorph {
    fn dof(&self) -> usize {
        3 * self.centers.len()
    }

    fn apply(&self, theta: &[f64], x: Point3) -> Result<Point3, XformError> {
        self.check_theta(theta)?;
        Ok(x.offset(self.weighted_sum(theta, x)))
    }

    fn jacobian_action(
        &self,
        theta: &[f64],
        dtheta: &[f64],
        x: Point3,
    ) -> Result<Vec3, XformError> {
        self.check_theta(theta)?;
        self.check_theta(dtheta)?;
        // Linear in θ.
        Ok(self.weighted_sum(dtheta, x))
    }

    fn spatial_jacobian(&self, theta: &[f64], x: Point3) -> Result<[[f64; 3]; 3], XformError> {
        self.check_theta(theta)?;
        let mut jac = [[0.0f64; 3]; 3];
        for (r, row) in jac.iter_mut().enumerate() {
            row[r] = 1.0;
        }
        for (i, &c) in self.centers.iter().enumerate() {
            let d = x.delta_from(c);
            let dist = d.norm();
            let q = dist / self.radius;
            if q >= 1.0 || dist < 1e-300 {
                // φ'(0) = 0: the kernel is flat at its center — smooth.
                continue;
            }
            let dphi = wendland_c2_derivative(q) / (self.radius * dist);
            let g = [dphi * d.x, dphi * d.y, dphi * d.z];
            for comp in 0..3 {
                let p = theta[3 * i + comp];
                for a in 0..3 {
                    jac[comp][a] += g[a] * p;
                }
            }
        }
        Ok(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kernel_shape_and_compact_support() {
        assert!((wendland_c2(0.0) - 1.0).abs() < 1e-15);
        assert!(wendland_c2(1.0).to_bits() == 0.0_f64.to_bits());
        assert!(wendland_c2(2.0).to_bits() == 0.0_f64.to_bits());
        assert!(
            wendland_c2_derivative(0.0).to_bits() == 0.0_f64.to_bits(),
            "flat at the center (C2)"
        );
        assert!(wendland_c2_derivative(1.5).to_bits() == 0.0_f64.to_bits());
        // Derivative matches a central difference mid-support.
        let (q, h) = (0.4, 1e-6);
        let fd = (wendland_c2(q + h) - wendland_c2(q - h)) / (2.0 * h);
        assert!((wendland_c2_derivative(q) - fd).abs() < 1e-8);
    }

    #[test]
    fn displacement_vanishes_outside_support() {
        let morph = RbfMorph {
            centers: vec![Point3::new(0.0, 0.0, 0.0)],
            radius: 1.0,
        };
        let theta = vec![1.0, 2.0, 3.0];
        let far = Point3::new(5.0, 0.0, 0.0);
        let y = morph.apply(&theta, far).unwrap();
        for (a, b) in [(y.x, 5.0), (y.y, 0.0), (y.z, 0.0)] {
            assert!((a - b).abs() < 1e-15, "outside support must pass through");
        }
        let near = Point3::new(0.25, 0.0, 0.0);
        let z = morph.apply(&theta, near).unwrap();
        assert!(z.x > 0.25 && z.y > 0.0, "inside support the handle pulls");
    }
}
