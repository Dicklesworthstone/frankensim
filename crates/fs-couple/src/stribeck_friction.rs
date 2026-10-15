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

/// A force/velocity pair satisfying a discrete collocated contact equation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImplicitFrictionStep {
    /// Traction on the driven body [N].
    pub traction_n: f64,
    /// Driver-minus-body velocity after the held-force step [m/s].
    pub relative_velocity_m_s: f64,
    /// Signed velocity residual of the constitutive/kinematic solve [m/s].
    pub residual_m_s: f64,
}

/// A contact refusal never substitutes an explicit force or a different law.
#[derive(Debug, Clone, PartialEq)]
pub enum ImplicitFrictionError<E> {
    /// Invalid input or an unrepresentable intermediate value.
    InvalidInput(&'static str),
    /// The caller's constitutive owner refused an evaluation.
    Constitutive(E),
    /// The supplied sliding law would add relative mechanical energy.
    NonDissipative,
    /// A discontinuity or unresolved root prevented a consistent step.
    NoConvergence,
}

/// Solve `v = v_free - mobility * F(v)` for an odd, dissipative sliding law.
///
/// `mobility` is the *discrete* collocated endpoint velocity response to one
/// newton held over this step [m/(N s)], not continuous-time inverse mass.
/// Both quantities must come from the same transition that the caller commits.
/// A string, a brake and a sliding mechanical port use the same scalar solve;
/// the constitutive callback retains ownership of materials and state validity.
///
/// The bracket does not allow traction to reverse the slip direction. Static
/// reaction/capture remains the caller's separate constraint solve. A
/// velocity-weakening law can have multiple roots; this deterministic bracket
/// selects one and does not certify uniqueness or continuous-time passivity.
/// In particular, an interval shrinking onto a constitutive jump is NOT
/// convergence: the velocity residual must also pass. There is no allocation.
///
/// # Errors
/// Invalid response data, constitutive refusal, active traction, or failure to
/// resolve a root to floating-point tolerance within 96 bisections.
pub fn solve_implicit_friction<E>(
    free_relative_velocity_m_s: f64,
    mobility_m_per_n_s: f64,
    mut traction: impl FnMut(f64) -> Result<f64, E>,
) -> Result<ImplicitFrictionStep, ImplicitFrictionError<E>> {
    if !free_relative_velocity_m_s.is_finite() {
        return Err(ImplicitFrictionError::InvalidInput("free_relative_velocity"));
    }
    if !mobility_m_per_n_s.is_finite() || mobility_m_per_n_s < 0.0 {
        return Err(ImplicitFrictionError::InvalidInput("discrete_mobility"));
    }
    let sign = if free_relative_velocity_m_s < 0.0 { -1.0 } else { 1.0 };
    let speed = free_relative_velocity_m_s.abs();
    let mut evaluate = |relative_speed: f64| {
        let force = traction(sign * relative_speed)
            .map_err(ImplicitFrictionError::Constitutive)?;
        if !force.is_finite() {
            return Err(ImplicitFrictionError::InvalidInput("traction"));
        }
        if force * sign < 0.0 || (relative_speed == 0.0 && force != 0.0) {
            return Err(ImplicitFrictionError::NonDissipative);
        }
        let velocity_change = mobility_m_per_n_s * force;
        if !velocity_change.is_finite() {
            return Err(ImplicitFrictionError::InvalidInput("contact_velocity_change"));
        }
        let endpoint = free_relative_velocity_m_s - velocity_change;
        let residual = sign * relative_speed - endpoint;
        if !endpoint.is_finite() || !residual.is_finite() {
            return Err(ImplicitFrictionError::InvalidInput("contact_residual"));
        }
        Ok(ImplicitFrictionStep {
            traction_n: force,
            relative_velocity_m_s: endpoint,
            residual_m_s: residual,
        })
    };
    // Admit the zero-slip endpoint too: this helper is not a static-reaction
    // solver, and cannot bracket a law that invents a force at rest.
    let at_zero = evaluate(0.0)?;
    if speed == 0.0 {
        return Ok(at_zero);
    }
    let at_free = evaluate(speed)?;
    if mobility_m_per_n_s == 0.0 || at_free.traction_n == 0.0 {
        return Ok(at_free);
    }
    let tolerance = (64.0 * f64::EPSILON * speed).max(f64::MIN_POSITIVE);
    let (mut lo, mut hi) = (0.0, speed);
    for _ in 0..96 {
        let mid = lo + 0.5 * (hi - lo);
        let candidate = evaluate(mid)?;
        if candidate.residual_m_s.abs() <= tolerance
            && candidate.relative_velocity_m_s * sign >= 0.0
        {
            return Ok(candidate);
        }
        if mid == lo || mid == hi {
            break;
        }
        if candidate.residual_m_s * sign > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Err(ImplicitFrictionError::NoConvergence)
}

#[cfg(test)]
mod implicit_contact_tests {
    use super::*;

    #[test]
    fn implicit_linear_contact_matches_closed_form_across_mechanical_scales() {
        for free in [-3.0, -0.01, 0.0, 0.01, 3.0] {
            for mobility in [0.0, 1.0e-6, 0.01, 1.0, 100.0] {
                let damping = 4.0;
                let step = solve_implicit_friction(free, mobility,
                    |v| Ok::<_, ()>(damping * v)).unwrap();
                let exact = free / (1.0 + mobility * damping);
                let tol = 2.0e-12 * free.abs().max(1.0e-20);
                assert!((step.relative_velocity_m_s - exact).abs() < tol);
                assert!((step.traction_n - damping * exact).abs() < 4.0 * tol);
                assert!(step.traction_n * step.relative_velocity_m_s >= 0.0);
            }
        }
    }

    #[test]
    fn implicit_material_law_retains_constitutive_ownership_and_odd_symmetry() {
        let law = StribeckFriction::try_new(0.8, 0.3, 0.05).unwrap();
        for free in [0.01, 0.2, 2.0] {
            for mobility in [1.0e-5, 0.01, 0.2] {
                let solve = |sign| solve_implicit_friction(sign * free, mobility,
                    |v| law.traction(v, 0.7)).unwrap();
                let a = solve(1.0);
                let b = solve(-1.0);
                assert_eq!(a.traction_n.to_bits(), (-b.traction_n).to_bits());
                let owned = law.traction(a.relative_velocity_m_s, 0.7).unwrap();
                assert!((owned - a.traction_n).abs() < 1.0e-11);
                assert!((free - mobility * a.traction_n
                    - a.relative_velocity_m_s).abs() < 1.0e-14);
            }
        }
    }

    #[test]
    fn implicit_contact_prevents_explicit_velocity_reversal() {
        // An explicit viscous force gives v_next = 1 - 10 * 2 = -19 m/s.
        let step = solve_implicit_friction(1.0, 10.0,
            |v| Ok::<_, ()>(2.0 * v)).unwrap();
        assert!((step.relative_velocity_m_s - 1.0 / 21.0).abs() < 1.0e-12);
        assert!(step.traction_n * step.relative_velocity_m_s > 0.0);
    }

    #[test]
    fn implicit_contact_refuses_active_laws_bad_response_and_constitutive_jumps() {
        assert!(matches!(solve_implicit_friction(1.0, -1.0,
            |v| Ok::<_, ()>(v)), Err(ImplicitFrictionError::InvalidInput(_))));
        assert!(matches!(solve_implicit_friction(f64::NAN, 1.0,
            |v| Ok::<_, ()>(v)), Err(ImplicitFrictionError::InvalidInput(_))));
        assert!(matches!(solve_implicit_friction(1.0, 1.0,
            |v| Ok::<_, ()>(-v)), Err(ImplicitFrictionError::NonDissipative)));
        assert_eq!(solve_implicit_friction(1.0, 1.0,
            |_| Err::<f64, _>("outside material domain")),
            Err(ImplicitFrictionError::Constitutive("outside material domain")));
        // There is no root: left limit is -0.5 and right limit is +0.5.
        assert_eq!(solve_implicit_friction(1.0, 1.0,
            |v| Ok::<_, ()>(if v < 0.5 { 0.0 } else { 1.0 })),
            Err(ImplicitFrictionError::NoConvergence));
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
