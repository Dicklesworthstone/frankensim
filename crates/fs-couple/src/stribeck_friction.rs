//! Coupling-layer 1-D ODE embedding of the fs-tribo Stribeck rung.
//!
//! `fs_tribo::FrictionLaw::evaluate` returns zero traction at rest
//! (it never invents a stick reaction). An explicit modal stepper
//! needs a continuous force, so this type is the same regularized
//! ramp `FrictionLaw::regularized_traction_1d` exposes, with the
//! driven-body sign used by a coupling port (`+` when the driver is
//! faster). A bow, a brake, and a fault are the same law.
//!
//! Evaluation delegates to that shared owner, including its deterministic
//! exponential and input/overflow checks. Coefficients remain caller supplied;
//! this adapter does not admit a sourced interface or invent a stick reaction.

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
    /// Refuses invalid coefficients (including direct struct construction),
    /// negative/non-finite normal loads, non-finite velocities, and
    /// unrepresentable results through the shared friction owner.
    pub fn traction(self, v_rel: f64, normal_n: f64) -> Result<f64, &'static str> {
        fs_tribo::FrictionLaw::Stribeck {
            static_mu: self.mu_static,
            kinetic_mu: self.mu_dynamic,
            characteristic_speed: self.stiction_m_s,
            viscous_per_speed: 0.0,
        }
        // fs-tribo takes body-minus-driver velocity. Its opposing traction
        // then acts in the driver's direction, as this coupling port requires.
        .regularized_traction_1d(-v_rel, normal_n, self.stiction_m_s)
        .map_err(|error| match error {
            fs_tribo::TriboError::InvalidInput { field } => field,
            _ => "shared friction evaluation refused",
        })
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
        assert_eq!(law.traction(0.01, 0.0).unwrap().to_bits(), 0.0f64.to_bits());
        assert!(StribeckFriction::try_new(0.3, 0.8, 0.05).is_err());
        assert!(StribeckFriction::try_new(f64::NAN, 0.3, 0.05).is_err());
        assert!(StribeckFriction::try_new(0.8, 0.3, 0.0).is_err());
    }

    /// G1/G3: equal-and-opposite contact forces remove relative kinetic
    /// energy; the independent formula also checks the port's sign and SI load.
    #[test]
    fn shared_friction_matches_reference_and_dissipates_relative_work() {
        let law = StribeckFriction::try_new(0.8, 0.3, 0.05).unwrap();
        let owner = fs_tribo::FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.3,
            characteristic_speed: 0.05,
            viscous_per_speed: 0.0,
        };
        for velocity in [-10.0_f64, -0.1, -0.05, -0.01, 0.0, 0.01, 0.05, 0.1, 10.0] {
            for load in [0.0, 0.1, 2.0, 100.0] {
                let force = law.traction(velocity, load).unwrap();
                assert_eq!(
                    force.to_bits(),
                    owner
                        .regularized_traction_1d(-velocity, load, 0.05)
                        .unwrap()
                        .to_bits()
                );
                let ratio = velocity.abs() / 0.05;
                let mu = if ratio < 1.0 {
                    0.8 * ratio
                } else {
                    0.3 + 0.5 * (-ratio * ratio).exp()
                };
                let reference = mu * load * velocity.signum();
                assert!((force - reference).abs() <= 1.0e-12 * load.max(1.0));
                let body_velocity = 0.7;
                let driver_velocity = body_velocity + velocity;
                let pair_power = force * body_velocity - force * driver_velocity;
                assert!(pair_power <= 1.0e-12);
                assert!((pair_power + force * velocity).abs() <= 1.0e-12);
            }
        }
    }

    #[test]
    fn shared_friction_refuses_forged_coefficients_negative_load_and_overflow() {
        let valid = StribeckFriction::try_new(0.8, 0.3, 0.05).unwrap();
        for invalid in [
            StribeckFriction {
                mu_static: -0.8,
                ..valid
            },
            StribeckFriction {
                mu_dynamic: 0.9,
                ..valid
            },
            StribeckFriction {
                mu_dynamic: f64::NAN,
                ..valid
            },
            StribeckFriction {
                stiction_m_s: 0.0,
                ..valid
            },
        ] {
            assert!(invalid.traction(0.01, 2.0).is_err());
        }
        assert_eq!(valid.traction(0.01, -1.0), Err("normal_force"));
        let large = StribeckFriction::try_new(4.0, 4.0, 0.05).unwrap();
        assert_eq!(large.traction(0.1, f64::MAX), Err("regularized_traction"));
    }
}
