//! Ideal Coulomb impulses with an isotropic, decoupled tangential response.
//!
//! The normal impulse is supplied by its existing contact owner and is never
//! modified here. The tangent response is `v_after = v_before + w_t*p_t`.
//! With one common static/kinetic coefficient, the exact maximum-dissipation
//! solution is the projection of `-v_before/w_t` onto the disk
//! `|p_t| <= mu*p_normal`. This is a Coulomb tangent disk at a fixed normal
//! impulse, not an associated normal/tangent cone law.
//!
//! Isotropic sphere/plane and sphere/sphere contact supply this response when
//! the inverse effective mass includes spin. Geometry, the tangent basis,
//! inertia, moving-surface work, and timing belong to the consuming mechanics.
//! General anisotropic or normal/tangent-coupled contacts are outside this rung.
//! This instantaneous rigid law is distinct from finite-patch partial slip and
//! held-effort constitutive integration. No heat partition is inferred.

use fs_math::det;
use fs_tribo::TriboError;
/// The original constitutive owner and branch labels; no new coefficient law.
pub use fs_tribo::{FrictionLaw, FrictionRegime};

const RESPONSE_TOLERANCE: f64 = 256.0 * f64::EPSILON;

/// Work-conjugate inputs in any one consistent unit system. With SI inputs,
/// normal impulse is N s, tangent velocity is m/s, and inverse mass is 1/kg.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoulombImpulseInput {
    /// Nonnegative finite normal impulse already admitted by normal mechanics.
    pub normal_impulse: f64,
    /// Relative contact-point velocity in an orthonormal tangent basis.
    pub relative_tangent_velocity: [f64; 2],
    /// Positive scalar inverse effective mass, including both translation and
    /// rotation. For one sphere against a plane, `1/m + radius^2/I`.
    pub inverse_tangent_mass: f64,
}

/// Accepted ideal-Coulomb response. With SI inputs, impulse is N s, velocity
/// is m/s, and loss/work residual is J. Numerical residuals remain explicit.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct CoulombImpulseResult {
    /// Tangent impulse opposing the initial relative slip.
    pub tangent_impulse: [f64; 2],
    /// Actual `w_t.mul_add(p_t, v_before)` in each tangent coordinate.
    pub relative_tangent_velocity_after: [f64; 2],
    /// Arrest within capacity, or continuing sliding at the friction limit.
    pub regime: FrictionRegime,
    /// Nonnegative kinetic loss of the tangential relative mode.
    pub dissipated_energy: f64,
    /// Actual endpoint velocity minus the ideal radial solution.
    pub impulse_velocity_residual: [f64; 2],
    /// Scalar-mode kinetic loss plus independently evaluated midpoint impulse work.
    pub work_residual: f64,
    /// `mu * normal_impulse`.
    pub friction_capacity: f64,
    /// Nonnegative excess of the measured impulse magnitude over capacity.
    pub friction_limit_residual: f64,
}

/// Refusal leaves all geometry and mechanical state owned by the caller intact.
#[derive(Debug, Clone, PartialEq)]
pub enum CoulombImpulseError {
    /// The original tribology owner refused the constitutive parameters.
    Constitutive(TriboError),
    /// Only the declared constant Coulomb rung is admitted here.
    UnsupportedLaw,
    /// Static and kinetic coefficients differ; this owner has one friction
    /// disk so its impulse and continuous slip-arrest interpretations agree.
    UnequalCoefficients,
    /// A kinematic or inertial input is inadmissible.
    InvalidInput(&'static str),
    /// A required intermediate or result is not finite.
    NonFinite,
    /// Floating-point resolution could not close velocity, capacity, or work.
    UnresolvedResponse(&'static str),
}

impl core::fmt::Display for CoulombImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "isotropic Coulomb impulse: {self:?}")
    }
}
impl std::error::Error for CoulombImpulseError {}

/// Solve the two-dimensional ideal-Coulomb maximum-dissipation impulse.
///
/// The original [`FrictionLaw::Coulomb`] must have equal static and kinetic
/// coefficients. Constitutive validation and coefficient queries remain in
/// `fs-tribo`. The candidate sticks when the stopping impulse fits the supplied
/// friction disk, otherwise it slides with magnitude `mu*normal_impulse`.
///
/// The same law can consume an impact impulse or a resolved support impulse
/// over time. The caller must resolve the actual slip-arrest time when evolving
/// configuration; this function does not advance positions or repair contact.
/// Cancellation is unnecessary for this fixed-work, allocation-free operation.
pub fn solve_isotropic_coulomb_impulse(
    input: CoulombImpulseInput,
    law: &FrictionLaw,
) -> Result<CoulombImpulseResult, CoulombImpulseError> {
    law.validate().map_err(CoulombImpulseError::Constitutive)?;
    if !matches!(law, FrictionLaw::Coulomb { .. }) {
        return Err(CoulombImpulseError::UnsupportedLaw);
    }
    let mu = law
        .static_coefficient()
        .map_err(CoulombImpulseError::Constitutive)?;
    let kinetic = law
        .zero_speed_kinetic_coefficient()
        .map_err(CoulombImpulseError::Constitutive)?;
    if mu != kinetic {
        return Err(CoulombImpulseError::UnequalCoefficients);
    }
    if !input.normal_impulse.is_finite() || input.normal_impulse < 0.0 {
        return Err(CoulombImpulseError::InvalidInput("normal impulse"));
    }
    let w = input.inverse_tangent_mass;
    if !w.is_finite() || w <= 0.0 {
        return Err(CoulombImpulseError::InvalidInput("inverse tangential mass"));
    }
    let velocity = input.relative_tangent_velocity;
    if !velocity.into_iter().all(f64::is_finite) {
        return Err(CoulombImpulseError::InvalidInput("tangent velocity"));
    }
    let speed = norm(velocity)?;
    let capacity = finite(mu * input.normal_impulse)?;
    if speed == 0.0 || capacity == 0.0 {
        return Ok(CoulombImpulseResult {
            tangent_impulse: [0.0; 2],
            relative_tangent_velocity_after: velocity,
            regime: if speed == 0.0 {
                FrictionRegime::Sticking
            } else {
                FrictionRegime::Sliding
            },
            dissipated_energy: 0.0,
            impulse_velocity_residual: [0.0; 2],
            work_residual: 0.0,
            friction_capacity: capacity,
            friction_limit_residual: 0.0,
        });
    }
    // Compare against stopping capacity without forming an overflowing w*cap
    // or an unrepresentable stopping impulse that a sliding branch never needs.
    let stopping = if w > 1.0 && capacity > f64::MAX / w {
        true
    } else {
        speed <= w * capacity
    };
    let impulse = if stopping {
        velocity.map(|v| -v / w)
    } else {
        velocity.map(|v| -capacity * (v / speed))
    };
    for value in impulse {
        finite(value)?;
    }
    let magnitude = norm(impulse)?;
    let decrement = finite(w * magnitude)?;
    let target_speed = if stopping {
        0.0
    } else {
        (-w).mul_add(capacity, speed)
    };
    let target = velocity.map(|v| (target_speed / speed) * v);
    let after = core::array::from_fn(|i| w.mul_add(impulse[i], velocity[i]));
    let velocity_residual = core::array::from_fn(|i| after[i] - target[i]);
    for value in after {
        finite(value)?;
    }
    if norm(velocity_residual)? / speed > RESPONSE_TOLERANCE {
        return Err(CoulombImpulseError::UnresolvedResponse("endpoint velocity"));
    }
    let limit_residual = (magnitude - capacity).max(0.0);
    if magnitude > 0.0 && limit_residual / magnitude > RESPONSE_TOLERANCE {
        return Err(CoulombImpulseError::UnresolvedResponse("friction capacity"));
    }
    let dissipated_energy = finite(magnitude * (speed - 0.5 * decrement))?;
    let midpoint = core::array::from_fn::<_, 2, _>(|i| 0.5 * velocity[i] + 0.5 * after[i]);
    let work = finite(impulse[0].mul_add(midpoint[0], impulse[1] * midpoint[1]))?;
    if dissipated_energy < 0.0 || work > 0.0 {
        return Err(CoulombImpulseError::UnresolvedResponse(
            "nonpassive impulse",
        ));
    }
    let work_residual = finite(dissipated_energy + work)?;
    let energy_scale = dissipated_energy.max(work.abs());
    if energy_scale > 0.0 && work_residual.abs() / energy_scale > RESPONSE_TOLERANCE {
        return Err(CoulombImpulseError::UnresolvedResponse("midpoint work"));
    }
    Ok(CoulombImpulseResult {
        tangent_impulse: impulse,
        relative_tangent_velocity_after: after,
        regime: if stopping {
            FrictionRegime::Sticking
        } else {
            FrictionRegime::Sliding
        },
        dissipated_energy,
        impulse_velocity_residual: velocity_residual,
        work_residual,
        friction_capacity: capacity,
        friction_limit_residual: limit_residual,
    })
}

fn finite(value: f64) -> Result<f64, CoulombImpulseError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CoulombImpulseError::NonFinite)
    }
}
fn norm(vector: [f64; 2]) -> Result<f64, CoulombImpulseError> {
    finite(det::hypot(vector[0], vector[1]))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn law(mu: f64) -> FrictionLaw {
        FrictionLaw::Coulomb {
            static_mu: mu,
            kinetic_mu: mu,
        }
    }
    fn input() -> CoulombImpulseInput {
        CoulombImpulseInput {
            normal_impulse: 4.0,
            relative_tangent_velocity: [1.0, 0.0],
            inverse_tangent_mass: 1.75,
        }
    }
    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 2e-13 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn sphere_spin_inertia_gives_exact_sticking_and_sliding_impulses() {
        // m=2, R=.5, I=.2: inverse tangent mass=1/m+R^2/I=1.75.
        // An initially spinning sphere has unit contact slip even at zero
        // translational speed. The solver must use contact-point kinematics.
        let sticking = solve_isotropic_coulomb_impulse(input(), &law(0.2)).unwrap();
        assert_eq!(sticking.regime, FrictionRegime::Sticking);
        near(sticking.tangent_impulse[0], -4.0 / 7.0);
        near(sticking.relative_tangent_velocity_after[0], 0.0);
        near(sticking.dissipated_energy, 2.0 / 7.0);
        let sliding = solve_isotropic_coulomb_impulse(input(), &law(0.05)).unwrap();
        assert_eq!(sliding.regime, FrictionRegime::Sliding);
        near(sliding.tangent_impulse[0], -0.2);
        near(sliding.relative_tangent_velocity_after[0], 0.65);
        near(sliding.dissipated_energy, 0.165);
    }

    #[test]
    fn two_axis_projection_is_maximum_dissipation_and_closes_kinetic_work() {
        let request = CoulombImpulseInput {
            relative_tangent_velocity: [3.0, -4.0],
            inverse_tangent_mass: 2.0,
            normal_impulse: 2.0,
        };
        for mu in [0.25, 1.25, 2.0] {
            let result = solve_isotropic_coulomb_impulse(request, &law(mu)).unwrap();
            let objective = |p: [f64; 2]| p[0] * 3.0 - p[1] * 4.0 + p[0] * p[0] + p[1] * p[1];
            let accepted = objective(result.tangent_impulse);
            for direction in [
                [1.0, 0.0],
                [-1.0, 0.0],
                [0.0, 1.0],
                [0.0, -1.0],
                [0.6, -0.8],
            ] {
                for fraction in [0.0, 0.25, 0.5, 1.0] {
                    let candidate = direction.map(|v| v * fraction * mu * 2.0);
                    assert!(accepted <= objective(candidate) + 2e-13);
                }
            }
            let after = result.relative_tangent_velocity_after;
            let kinetic_change = (after[0] * after[0] + after[1] * after[1] - 25.0) / 4.0;
            near(kinetic_change, -result.dissipated_energy);
            assert!(result.work_residual.abs() < 1e-13);
            assert!(result.friction_limit_residual < 1e-13);
        }
    }

    #[test]
    fn tangent_rotation_and_unit_scaling_preserve_the_response() {
        let baseline = solve_isotropic_coulomb_impulse(input(), &law(0.05)).unwrap();
        let rotate = |v: [f64; 2]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1]];
        let rotated = solve_isotropic_coulomb_impulse(
            CoulombImpulseInput {
                relative_tangent_velocity: rotate(input().relative_tangent_velocity),
                ..input()
            },
            &law(0.05),
        )
        .unwrap();
        for i in 0..2 {
            near(
                rotated.tangent_impulse[i],
                rotate(baseline.tangent_impulse)[i],
            );
            near(
                rotated.relative_tangent_velocity_after[i],
                rotate(baseline.relative_tangent_velocity_after)[i],
            );
        }
        near(rotated.dissipated_energy, baseline.dissipated_energy);
        let scaled = solve_isotropic_coulomb_impulse(
            CoulombImpulseInput {
                normal_impulse: input().normal_impulse * 32.0,
                relative_tangent_velocity: input().relative_tangent_velocity.map(|v| v * 8.0),
                inverse_tangent_mass: input().inverse_tangent_mass / 4.0,
            },
            &law(0.05),
        )
        .unwrap();
        near(
            scaled.tangent_impulse[0],
            baseline.tangent_impulse[0] * 32.0,
        );
        near(scaled.dissipated_energy, baseline.dissipated_energy * 256.0);
    }

    #[test]
    fn zero_load_zero_friction_and_rest_have_no_invented_reaction() {
        for request in [
            CoulombImpulseInput {
                normal_impulse: 0.0,
                ..input()
            },
            CoulombImpulseInput {
                relative_tangent_velocity: [0.0; 2],
                ..input()
            },
        ] {
            let result = solve_isotropic_coulomb_impulse(request, &law(0.2)).unwrap();
            assert_eq!(result.tangent_impulse, [0.0; 2]);
            assert_eq!(
                result.relative_tangent_velocity_after,
                request.relative_tangent_velocity
            );
            assert_eq!(result.dissipated_energy, 0.0);
        }
        let frictionless = solve_isotropic_coulomb_impulse(input(), &law(0.0)).unwrap();
        assert_eq!(frictionless.tangent_impulse, [0.0; 2]);
        assert_eq!(frictionless.relative_tangent_velocity_after, [1.0, 0.0]);
    }

    #[test]
    fn constitutive_model_and_kinematic_refusals_are_explicit() {
        let unequal = FrictionLaw::Coulomb {
            static_mu: 0.3,
            kinetic_mu: 0.2,
        };
        assert_eq!(
            solve_isotropic_coulomb_impulse(input(), &unequal),
            Err(CoulombImpulseError::UnequalCoefficients)
        );
        let stribeck = FrictionLaw::Stribeck {
            static_mu: 0.3,
            kinetic_mu: 0.2,
            characteristic_speed: 1.0,
            viscous_per_speed: 0.0,
        };
        assert_eq!(
            solve_isotropic_coulomb_impulse(input(), &stribeck),
            Err(CoulombImpulseError::UnsupportedLaw)
        );
        assert!(matches!(
            solve_isotropic_coulomb_impulse(input(), &law(-0.1)),
            Err(CoulombImpulseError::Constitutive(_))
        ));
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                solve_isotropic_coulomb_impulse(
                    CoulombImpulseInput {
                        inverse_tangent_mass: value,
                        ..input()
                    },
                    &law(0.2)
                )
                .is_err()
            );
        }
        for value in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                solve_isotropic_coulomb_impulse(
                    CoulombImpulseInput {
                        normal_impulse: value,
                        ..input()
                    },
                    &law(0.2)
                )
                .is_err()
            );
        }
        assert!(
            solve_isotropic_coulomb_impulse(
                CoulombImpulseInput {
                    relative_tangent_velocity: [0.0, f64::NAN],
                    ..input()
                },
                &law(0.2)
            )
            .is_err()
        );
    }

    #[test]
    fn scale_safe_norm_and_unrepresentable_responses_do_not_publish_false_balance() {
        let large = solve_isotropic_coulomb_impulse(
            CoulombImpulseInput {
                normal_impulse: 1e150,
                relative_tangent_velocity: [1e150, -1e150],
                inverse_tangent_mass: 2.0,
            },
            &law(1.0),
        )
        .unwrap();
        assert_eq!(large.regime, FrictionRegime::Sticking);
        near(large.dissipated_energy / 1e300, 0.5);
        assert_eq!(
            solve_isotropic_coulomb_impulse(
                CoulombImpulseInput {
                    normal_impulse: f64::MAX,
                    ..input()
                },
                &law(2.0)
            ),
            Err(CoulombImpulseError::NonFinite)
        );
        assert!(
            solve_isotropic_coulomb_impulse(
                CoulombImpulseInput {
                    normal_impulse: 1.0,
                    relative_tangent_velocity: [f64::from_bits(1), 0.0],
                    inverse_tangent_mass: f64::MAX,
                },
                &law(1.0)
            )
            .is_err()
        );
    }
}
