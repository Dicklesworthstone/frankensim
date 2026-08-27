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
    /// Validate friction coefficients.
    ///
    /// # Errors
    /// Refuses non-finite or negative coefficients and a non-positive
    /// stiction scale (fresh-eyes review bead 9svup: the former
    /// silent-zero fallbacks in [`StribeckFriction::traction`] masked
    /// misconfiguration as physics).
    pub fn try_new(
        mu_static: f64,
        mu_dynamic: f64,
        stiction_m_s: f64,
    ) -> Result<Self, &'static str> {
        if !mu_static.is_finite() || mu_static < 0.0 {
            return Err("mu_static must be finite and non-negative");
        }
        if !mu_dynamic.is_finite() || mu_dynamic < 0.0 || mu_dynamic > mu_static {
            return Err("mu_dynamic must be finite, non-negative, and <= mu_static");
        }
        if !stiction_m_s.is_finite() || stiction_m_s <= 0.0 {
            return Err("stiction scale must be finite and strictly positive");
        }
        Ok(Self {
            mu_static,
            mu_dynamic,
            stiction_m_s,
        })
    }

    /// Tangential force [N] on the driven body for relative velocity
    /// `v_rel = v_driver − v_driven` and normal `n`.
    ///
    /// # Errors
    /// Refuses non-finite inputs instead of silently returning zero
    /// traction (bead 9svup).
    pub fn traction(self, v_rel: f64, normal_n: f64) -> Result<f64, &'static str> {
        let a = v_rel.abs();
        let v0 = self.stiction_m_s;
        if !v_rel.is_finite() || !normal_n.is_finite() {
            return Err("velocity and normal force must be finite");
        }
        debug_assert!(v0 > 0.0 && v0.is_finite(), "construct via try_new");
        let mu = if a < v0 {
            self.mu_static * a / v0
        } else {
            let decay = (-(a / v0) * (a / v0)).exp();
            self.mu_dynamic + (self.mu_static - self.mu_dynamic) * decay
        };
        Ok(mu * normal_n * v_rel.signum())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn odd_in_velocity_and_falls_to_kinetic() {
        let law = StribeckFriction::try_new(0.8, 0.3, 0.05).unwrap();
        let f_p = law.traction(0.01, 2.0).unwrap();
        let f_m = law.traction(-0.01, 2.0).unwrap();
        assert!((f_p + f_m).abs() < 1.0e-15);
        let slope = law.traction(0.01, 1.0).unwrap() / 0.01;
        assert!((slope - 0.8 / 0.05).abs() < 1.0e-12);
        assert!((law.traction(10.0, 1.0).unwrap() - 0.3).abs() < 0.02);
    }

    /// Regression (bead 9svup): non-finite inputs REFUSE instead of
    /// publishing zero traction; misconfigured coefficients refuse at
    /// construction.
    #[test]
    fn nonfinite_inputs_refuse_as_typed_errors() {
        let law = StribeckFriction::try_new(0.8, 0.3, 0.05).unwrap();
        assert!(law.traction(f64::NAN, 1.0).is_err());
        assert!(law.traction(f64::INFINITY, 1.0).is_err());
        assert!(law.traction(0.01, f64::NAN).is_err());
        assert_eq!(law.traction(0.01, 0.0).unwrap(), 0.0);
        assert!(StribeckFriction::try_new(0.3, 0.8, 0.05).is_err());
        assert!(StribeckFriction::try_new(f64::NAN, 0.3, 0.05).is_err());
        assert!(StribeckFriction::try_new(0.8, 0.3, 0.0).is_err());
    }
}
