//! Coupling-layer 1-D ODE embedding of the fs-tribo Stribeck rung.
//!
//! `fs_tribo::FrictionLaw::evaluate` returns zero traction at rest
//! (it never invents a stick reaction). An explicit modal stepper
//! needs a continuous force, so this type is the same regularized
//! ramp `FrictionLaw::regularized_traction_1d` exposes, with the
//! driven-body sign used by a coupling port (`+` when the driver is
//! faster). A bow, a brake, and a fault are the same law.
//!
//! The coefficients match that Stribeck rung
//! (`μ_k + (μ_s−μ_k) exp(−(v/v₀)²)`, viscous term 0). A direct
//! `fs-tribo` dependency waits on `Cargo.toml`.

/// Regularized friction coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StribeckFriction {
    /// Static coefficient.
    pub mu_static: f64,
    /// Dynamic coefficient (`<= mu_static`).
    pub mu_dynamic: f64,
    /// Velocity scale of the stiction ramp and Stribeck decay [m/s].
    pub stiction_m_s: f64,
}

impl StribeckFriction {
    /// Tangential force [N] on the driven body for relative velocity
    /// `v_rel = v_driver − v_driven` and normal `n`.
    #[must_use]
    pub fn traction(self, v_rel: f64, normal_n: f64) -> f64 {
        let a = v_rel.abs();
        let v0 = self.stiction_m_s;
        if !(v0 > 0.0) || !v_rel.is_finite() || !normal_n.is_finite() {
            return 0.0;
        }
        let mu = if a < v0 {
            self.mu_static * a / v0
        } else {
            let decay = (-(a / v0) * (a / v0)).exp();
            self.mu_dynamic + (self.mu_static - self.mu_dynamic) * decay
        };
        mu * normal_n * v_rel.signum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn odd_in_velocity_and_falls_to_kinetic() {
        let law = StribeckFriction {
            mu_static: 0.8,
            mu_dynamic: 0.3,
            stiction_m_s: 0.05,
        };
        assert!((law.traction(0.01, 2.0) + law.traction(-0.01, 2.0)).abs() < 1.0e-15);
        let slope = law.traction(0.01, 1.0) / 0.01;
        assert!((slope - 0.8 / 0.05).abs() < 1.0e-12);
        assert!((law.traction(10.0, 1.0) - 0.3).abs() < 0.02);
    }
}
