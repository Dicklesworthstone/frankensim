//! Frictionless Newton restitution for one active normal contact.
//!
//! This is an instantaneous, solver-independent scalar law. The caller owns
//! geometry, contact selection, the normal Jacobian, and the inverse effective
//! mass `w = J M^-1 J^T`. Negative relative normal velocity means approach;
//! a positive normal impulse increases that velocity. Configuration must
//! already be admissible: this law never repairs penetration or advances time.
//!
//! A sequence of these single-contact solves is not a globally admissible
//! multicontact impact law. Friction, simultaneous impact selection, persistent
//! contact, gap stabilization, and partitioning loss into heat or damage require
//! their own owners. The reported dissipation is not automatically thermal heat.

/// Admitted scalar geometry, relative motion, and inertia of one contact.
/// Use one consistent unit system: with SI inputs, gap is in m, speed in m/s,
/// and inverse effective mass in 1/kg.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalImpactInput {
    /// Finite nonnegative signed separation. Any penetration is refused,
    /// including one smaller than the activation tolerance.
    pub gap: f64,
    /// Finite nonnegative gap at or below which the contact is active.
    /// This admits proximity; it is not a penetration-repair allowance.
    pub activation_tolerance: f64,
    /// Relative contact-point velocity along the separating normal.
    /// Negative values approach; zero and positive values need no impulse.
    pub relative_normal_speed: f64,
    /// Positive finite inverse effective mass `J M^-1 J^T`.
    pub inverse_effective_mass: f64,
    /// Newton coefficient of restitution in the closed interval `[0, 1]`.
    pub restitution: f64,
}

/// Scalar response for an admitted single contact. With SI inputs, impulse
/// is in N s, speed in m/s, and dissipated energy in J.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct NormalImpactResult {
    /// Nonnegative impulse along the separating normal; zero when inactive
    /// or nonapproaching. The caller applies action/reaction through its Jacobian.
    pub normal_impulse: f64,
    /// Newton target `-e*v_before` for an impact, or unchanged speed otherwise.
    /// The caller must check its actual state update against this target.
    pub relative_normal_speed_after: f64,
    /// Nonnegative loss `-(v_before + v_target)/2 * normal_impulse`.
    /// No contact elastic storage or heat-partition model is implied.
    pub dissipated_energy: f64,
    /// Whether the supplied gap is within the activation tolerance.
    pub active: bool,
    /// Whether this active contact was approaching and received an impulse.
    pub impacted: bool,
    /// Independently reconstructed `w*impulse + v_before - v_target`, using
    /// fused multiply-add to avoid overflow in a canceling velocity increment.
    pub impulse_velocity_residual: f64,
}

/// An invalid scalar contact state or a response that cannot be represented
/// faithfully by finite floating-point values. No state is mutated on refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalImpactError {
    /// Gap is negative or nonfinite; configuration repair is not this law.
    InvalidGap,
    /// Activation tolerance is negative or nonfinite.
    InvalidActivationTolerance,
    /// Relative normal speed is nonfinite.
    InvalidRelativeNormalSpeed,
    /// Inverse effective mass is nonpositive or nonfinite.
    InvalidInverseEffectiveMass,
    /// Restitution is nonfinite or outside `[0, 1]`.
    InvalidRestitution,
    /// A required impulse, velocity, or loss is not finitely representable.
    UnrepresentableResponse,
    /// The representable impulse cannot reproduce the restitution target to
    /// floating-point accuracy, for example after severe subnormal quantization.
    UnresolvedVelocityResponse,
}

impl core::fmt::Display for NormalImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidGap => {
                "normal impact requires a finite nonnegative gap; penetration is not repaired"
            }
            Self::InvalidActivationTolerance => {
                "normal-impact activation tolerance must be finite and nonnegative"
            }
            Self::InvalidRelativeNormalSpeed => "normal-impact relative speed must be finite",
            Self::InvalidInverseEffectiveMass => {
                "normal-impact inverse effective mass must be finite and positive"
            }
            Self::InvalidRestitution => "normal-impact restitution must lie in [0, 1]",
            Self::UnrepresentableResponse => "normal-impact response is not finitely representable",
            Self::UnresolvedVelocityResponse => {
                "normal-impact impulse cannot resolve the restitution velocity"
            }
        })
    }
}

impl core::error::Error for NormalImpactError {}

/// Solve one frictionless normal impact without changing configuration.
///
/// For an active approaching contact, `p = -(1+e)*v/w`, the Newton target is
/// `v_after = -e*v`, and loss is `(1-e)*(-v)*p/2`. The active normal law has
/// `p >= 0`, `v_after + e*v >= 0`, and their complementarity product is zero.
/// Inactive contacts retain their velocity; no velocity complementarity claim
/// is made for a contact outside the declared activation set.
///
/// Impulse and loss evaluation reorder only when a direct intermediate would
/// overflow. An overflowing final result, an unrepresentable positive impulse,
/// or a reconstructed velocity residual above `64*EPSILON*abs(v)` refuses.
/// Exactly elastic impacts avoid evaluating a potentially overflowing energy
/// product whose loss coefficient is zero. The routine allocates nothing.
///
/// # Errors
/// [`NormalImpactError`] for invalid inputs or an unresolved finite response.
pub fn solve_normal_impact(
    input: NormalImpactInput,
) -> Result<NormalImpactResult, NormalImpactError> {
    if !input.gap.is_finite() || input.gap < 0.0 {
        return Err(NormalImpactError::InvalidGap);
    }
    if !input.activation_tolerance.is_finite() || input.activation_tolerance < 0.0 {
        return Err(NormalImpactError::InvalidActivationTolerance);
    }
    if !input.relative_normal_speed.is_finite() {
        return Err(NormalImpactError::InvalidRelativeNormalSpeed);
    }
    if !input.inverse_effective_mass.is_finite() || input.inverse_effective_mass <= 0.0 {
        return Err(NormalImpactError::InvalidInverseEffectiveMass);
    }
    if !input.restitution.is_finite() || !(0.0..=1.0).contains(&input.restitution) {
        return Err(NormalImpactError::InvalidRestitution);
    }
    let active = input.gap <= input.activation_tolerance;
    if !active || input.relative_normal_speed >= 0.0 {
        return Ok(NormalImpactResult {
            normal_impulse: 0.0,
            relative_normal_speed_after: input.relative_normal_speed,
            dissipated_energy: 0.0,
            active,
            impacted: false,
            impulse_velocity_residual: 0.0,
        });
    }
    let approach_speed = -input.relative_normal_speed;
    let factor = 1.0 + input.restitution;
    let numerator = factor * approach_speed;
    let normal_impulse = if numerator.is_finite() {
        numerator / input.inverse_effective_mass
    } else {
        // Here approach_speed is large, so divide-first cannot suffer the
        // avoidable underflow possible for tiny speeds and finite numerator.
        (approach_speed / input.inverse_effective_mass) * factor
    };
    if !(normal_impulse.is_finite() && normal_impulse > 0.0) {
        return Err(NormalImpactError::UnrepresentableResponse);
    }
    let relative_normal_speed_after = input.restitution * approach_speed;
    let actual_speed = input
        .inverse_effective_mass
        .mul_add(normal_impulse, input.relative_normal_speed);
    let impulse_velocity_residual = actual_speed - relative_normal_speed_after;
    if !actual_speed.is_finite() || !impulse_velocity_residual.is_finite() {
        return Err(NormalImpactError::UnrepresentableResponse);
    }
    if impulse_velocity_residual.abs() / approach_speed > 64.0 * f64::EPSILON {
        return Err(NormalImpactError::UnresolvedVelocityResponse);
    }
    let coefficient = 0.5 * (1.0 - input.restitution);
    let dissipated_energy = if coefficient == 0.0 {
        0.0
    } else {
        let product = approach_speed * normal_impulse;
        if product.is_finite() {
            coefficient * product
        } else {
            // Product overflow implies both factors exceed one. For e<1,
            // coefficient >= 2^-54, so this scaling cannot underflow.
            (coefficient * approach_speed) * normal_impulse
        }
    };
    if !dissipated_energy.is_finite() {
        return Err(NormalImpactError::UnrepresentableResponse);
    }
    Ok(NormalImpactResult {
        normal_impulse,
        relative_normal_speed_after,
        dissipated_energy,
        active,
        impacted: true,
        impulse_velocity_residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input() -> NormalImpactInput {
        NormalImpactInput {
            gap: 0.0,
            activation_tolerance: 1e-9,
            relative_normal_speed: -4.0,
            inverse_effective_mass: 0.5,
            restitution: 0.5,
        }
    }

    #[test]
    fn normal_impulse_restitution_and_energy_match_effective_mass_oracle() {
        for restitution in [0.0, 0.3, 0.5, 1.0] {
            let request = NormalImpactInput {
                restitution,
                ..input()
            };
            let result = solve_normal_impact(request).unwrap();
            assert!(result.active && result.impacted);
            assert!((result.normal_impulse - 8.0 * (1.0 + restitution)).abs() < 2e-15);
            assert_eq!(result.relative_normal_speed_after, 4.0 * restitution);
            let actual_speed = request.relative_normal_speed
                + request.inverse_effective_mass * result.normal_impulse;
            let before_energy =
                request.relative_normal_speed.powi(2) / (2.0 * request.inverse_effective_mass);
            let after_energy = actual_speed.powi(2) / (2.0 * request.inverse_effective_mass);
            assert!((before_energy - after_energy - result.dissipated_energy).abs() < 5e-15);
            let slack = actual_speed + restitution * request.relative_normal_speed;
            assert!(slack.abs() < 2e-15);
            assert!((result.normal_impulse * slack).abs() < 2e-14);
            assert!(result.dissipated_energy >= 0.0);
        }
    }

    #[test]
    fn activation_boundary_and_separating_contacts_receive_no_attraction() {
        let mut request = input();
        request.gap = request.activation_tolerance;
        assert!(solve_normal_impact(request).unwrap().impacted);
        request.gap = request.activation_tolerance.next_up();
        let inactive = solve_normal_impact(request).unwrap();
        assert!(!inactive.active && !inactive.impacted);
        assert_eq!(inactive.normal_impulse, 0.0);
        assert_eq!(
            inactive.relative_normal_speed_after,
            request.relative_normal_speed
        );
        request.gap = 0.0;
        for speed in [-0.0, 0.0, 0.2, f64::MAX] {
            request.relative_normal_speed = speed;
            let separating = solve_normal_impact(request).unwrap();
            assert!(separating.active && !separating.impacted);
            assert_eq!(separating.normal_impulse, 0.0);
            assert_eq!(
                separating.relative_normal_speed_after.to_bits(),
                speed.to_bits()
            );
            assert_eq!(separating.dissipated_energy, 0.0);
        }
    }

    #[test]
    fn penetrating_and_malformed_inputs_are_refused_without_a_repair_tolerance() {
        for value in [-f64::MIN_POSITIVE, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                solve_normal_impact(NormalImpactInput {
                    gap: value,
                    ..input()
                }),
                Err(NormalImpactError::InvalidGap)
            );
            assert_eq!(
                solve_normal_impact(NormalImpactInput {
                    activation_tolerance: value,
                    ..input()
                }),
                Err(NormalImpactError::InvalidActivationTolerance)
            );
        }
        for value in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                solve_normal_impact(NormalImpactInput {
                    inverse_effective_mass: value,
                    ..input()
                }),
                Err(NormalImpactError::InvalidInverseEffectiveMass)
            );
        }
        for value in [f64::NEG_INFINITY, f64::INFINITY, f64::NAN] {
            assert_eq!(
                solve_normal_impact(NormalImpactInput {
                    relative_normal_speed: value,
                    ..input()
                }),
                Err(NormalImpactError::InvalidRelativeNormalSpeed)
            );
        }
        for value in [-f64::EPSILON, 1.0 + f64::EPSILON, f64::NAN, f64::INFINITY] {
            assert_eq!(
                solve_normal_impact(NormalImpactInput {
                    restitution: value,
                    ..input()
                }),
                Err(NormalImpactError::InvalidRestitution)
            );
        }
    }

    #[test]
    fn mass_and_speed_scaling_preserve_the_single_contact_law() {
        let baseline = solve_normal_impact(input()).unwrap();
        for mass_scale in [0.125, 1.0, 32.0] {
            for speed_scale in [0.25, 1.0, 16.0] {
                let scaled = solve_normal_impact(NormalImpactInput {
                    relative_normal_speed: input().relative_normal_speed * speed_scale,
                    inverse_effective_mass: input().inverse_effective_mass / mass_scale,
                    ..input()
                })
                .unwrap();
                assert_eq!(
                    scaled.normal_impulse,
                    baseline.normal_impulse * mass_scale * speed_scale
                );
                assert_eq!(
                    scaled.relative_normal_speed_after,
                    baseline.relative_normal_speed_after * speed_scale
                );
                assert_eq!(
                    scaled.dissipated_energy,
                    baseline.dissipated_energy * mass_scale * speed_scale * speed_scale
                );
            }
        }
    }

    #[test]
    fn finite_extreme_responses_avoid_intermediate_overflow_and_underflow() {
        let speed = 0.75 * f64::MAX;
        let elastic = solve_normal_impact(NormalImpactInput {
            relative_normal_speed: -speed,
            inverse_effective_mass: 2.0,
            restitution: 1.0,
            ..input()
        })
        .unwrap();
        assert_eq!(elastic.normal_impulse, speed);
        assert_eq!(elastic.relative_normal_speed_after, speed);
        assert_eq!(elastic.dissipated_energy, 0.0);
        let smallest = f64::from_bits(1);
        let tiny = solve_normal_impact(NormalImpactInput {
            relative_normal_speed: -smallest,
            inverse_effective_mass: 2.0,
            restitution: 1.0,
            ..input()
        })
        .unwrap();
        assert_eq!(tiny.normal_impulse, smallest);
        assert_eq!(tiny.relative_normal_speed_after, smallest);
        let restitution = 1.0 - 2.0_f64.powi(-20);
        let finite_loss = solve_normal_impact(NormalImpactInput {
            relative_normal_speed: -1e170,
            inverse_effective_mass: 1e30,
            restitution,
            ..input()
        })
        .unwrap();
        let expected = (0.5 * (1.0 - restitution) * (1.0 + restitution) * 1e140) * 1e170;
        assert!((finite_loss.dissipated_energy / expected - 1.0).abs() < 8e-16);
    }

    #[test]
    fn unrepresentable_or_unresolved_responses_refuse() {
        for request in [
            NormalImpactInput {
                relative_normal_speed: -f64::MAX,
                inverse_effective_mass: 0.5,
                ..input()
            },
            NormalImpactInput {
                relative_normal_speed: -1e200,
                inverse_effective_mass: 1.0,
                restitution: 0.0,
                ..input()
            },
            NormalImpactInput {
                relative_normal_speed: -f64::from_bits(1),
                inverse_effective_mass: f64::MAX,
                ..input()
            },
        ] {
            assert_eq!(
                solve_normal_impact(request),
                Err(NormalImpactError::UnrepresentableResponse)
            );
        }
        assert_eq!(
            solve_normal_impact(NormalImpactInput {
                relative_normal_speed: -1e-15,
                inverse_effective_mass: f64::MAX,
                restitution: 0.0,
                ..input()
            }),
            Err(NormalImpactError::UnresolvedVelocityResponse)
        );
    }
}
