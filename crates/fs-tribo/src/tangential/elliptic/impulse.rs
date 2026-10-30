//! Instantaneous maximum-dissipation impulses on a declared tangent ellipse.
//!
//! The normal impulse is already owned by normal mechanics. The tangent
//! problem minimizes p.dot(v) + p.dot(W*p)/2 inside that fixed ellipse.
//! The original stateless scalar/matrix projection supplies the solution;
//! no time interval, physical slip increment or thermal state is invented.
//!
//! Only one explicit ellipse is admitted: static and kinetic principal
//! coefficients must agree. This conditional tangent solve does not solve a
//! normal/tangent-coupled restitution law. Its consumer must supply the proper
//! normal owner and independently apply/audit the complete body impulse.

use super::compliance::{SymmetricResponse, resolve_matrix};
use super::*;

/// Physical impulse boundary; normal impulse is nonnegative and W is SPD.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticImpulseInput {
    pub normal_impulse_ns: f64,
    pub velocity_before_m_s: [f64; 2],
    pub mobility_per_kg: [[f64; 2]; 2],
}

/// Bounded work and independent physical acceptance tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticImpulseControls {
    /// At most 128 scalar root iterations, with no allocations.
    pub maximum_iterations: usize,
    pub maximum_impulse_ns: f64,
    pub impulse_absolute_tolerance_ns: f64,
    pub impulse_relative_tolerance: f64,
    pub velocity_tolerance_m_s: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticImpulseReport {
    /// Actual impulse on the positive-relative-velocity body; opposite to the
    /// positive traction convention used by the finite-displacement law.
    pub impulse_ns: [f64; 2],
    /// Actual v_before + W*p; a floating arrest residual is not projected away.
    pub velocity_after_m_s: [f64; 2],
    pub regime: FrictionRegime,
    /// Nonnegative mechanical impact loss; no thermal partition is inferred.
    pub dissipated_j: f64,
    /// Independently reconstructed in the inverse-mobility mass metric.
    pub kinetic_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub impulse_velocity_residual_m_s: [f64; 2],
    /// Ellipse excess scaled into impulse units, including forbidden zero axes.
    pub capacity_residual_ns: f64,
    /// Actual reaction versus the ellipse support at admitted endpoint slip.
    /// An explicitly arrested active component uses its admitted zero target.
    pub stationarity_residual_ns: f64,
    pub iterations: usize,
}

#[derive(Debug)]
pub enum EllipticImpulseError {
    Constitutive(TriboError),
    UnequalCoefficients,
    Invalid(&'static str),
    Unresolved(&'static str),
    Cancelled,
}
impl From<TriboError> for EllipticImpulseError {
    fn from(error: TriboError) -> Self {
        Self::Constitutive(error)
    }
}
impl core::fmt::Display for EllipticImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "directional Coulomb impulse: {self:?}")
    }
}
impl std::error::Error for EllipticImpulseError {}

impl EllipticCoulomb {
    /// Solve an instantaneous tangent impulse using the retained material axis.
    ///
    /// The borrowed history is checked for identity and never advanced. In
    /// particular impact loss is midpoint impulse work, not finite-slip heat:
    /// arrest can dissipate energy while instantaneous physical travel is zero.
    ///
    /// Cancellation is polled before and after the allocation-free projection
    /// and again after admission. Its longest unpolled unit is at most 128
    /// bounded two-dimensional scalar-root iterations.
    pub fn solve_impulse(
        &self,
        state: &EllipticCoulombState,
        input: &EllipticImpulseInput,
        controls: EllipticImpulseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EllipticImpulseReport, EllipticImpulseError> {
        poll(cancelled)?;
        self.check(state)?;
        if self.static_mu != self.kinetic_mu {
            return Err(EllipticImpulseError::UnequalCoefficients);
        }
        if !(1..=128).contains(&controls.maximum_iterations)
            || [
                controls.maximum_impulse_ns,
                controls.impulse_absolute_tolerance_ns,
                controls.impulse_relative_tolerance,
                controls.velocity_tolerance_m_s,
                controls.energy_tolerance_j,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0)
            || controls.impulse_relative_tolerance >= 1.0
        {
            return Err(EllipticImpulseError::Invalid(
                "positive finite physical tolerances and at most 128 iterations required",
            ));
        }
        if !input.normal_impulse_ns.is_finite()
            || input.normal_impulse_ns < 0.0
            || input
                .velocity_before_m_s
                .into_iter()
                .any(|x| !x.is_finite())
        {
            return Err(EllipticImpulseError::Invalid(
                "finite velocity and nonnegative normal impulse required",
            ));
        }
        let response = SymmetricResponse::new(input.mobility_per_kg)?;
        if response.scale == 0.0 || response.determinant <= 0.0 {
            return Err(EllipticImpulseError::Invalid(
                "strictly positive definite tangent mobility required",
            ));
        }
        let caps = self.kinetic_mu.map(|mu| mu * input.normal_impulse_ns);
        for capacity in caps {
            checked(capacity)?;
        }
        let free = to_principal(state.axis, input.velocity_before_m_s)?;
        let solve = ProjectionControls {
            maximum_iterations: controls.maximum_iterations,
            absolute_tolerance: controls.impulse_absolute_tolerance_ns,
            relative_tolerance: controls.impulse_relative_tolerance,
        };
        let m = response.matrix;
        let (reaction, endpoint, regime, iterations) = if m[0][1] == 0.0 && m[0][0] == m[1][1] {
            resolve(
                free,
                caps,
                caps,
                m[0][0],
                input.normal_impulse_ns > 0.0,
                solve,
            )?
        } else {
            let (c, determinant) = response.material(state.axis)?;
            resolve_matrix(
                free,
                caps,
                caps,
                c,
                determinant,
                response,
                state.axis,
                input.velocity_before_m_s,
                input.normal_impulse_ns > 0.0,
                solve,
            )?
        };
        poll(cancelled)?;
        let impulse_ns = from_principal(state.axis, reaction)?.map(|x| -x);
        let impulse_magnitude = norm(impulse_ns)?;
        if impulse_magnitude > controls.maximum_impulse_ns {
            return Err(EllipticImpulseError::Invalid(
                "maximum tangent impulse exceeded",
            ));
        }
        let before = input.velocity_before_m_s;
        let after = [0, 1]
            .map(|i| m[i][0].mul_add(impulse_ns[0], m[i][1].mul_add(impulse_ns[1], before[i])));
        for x in after {
            checked(x)?;
        }
        let target = from_principal(state.axis, endpoint)?;
        let velocity_residual = [0, 1].map(|i| after[i] - target[i]);
        if norm(velocity_residual)? > controls.velocity_tolerance_m_s {
            return Err(EllipticImpulseError::Unresolved(
                "impulse endpoint velocity",
            ));
        }

        // Reconstruct the physical reaction from the returned world-frame
        // impulse, including rotation residuals on a degenerate ellipse.
        let actual_reaction = to_principal(state.axis, impulse_ns.map(|x| -x))?;
        let capacity_residual_ns = capacity_residual(actual_reaction, caps)?;
        let tolerance = solve.tolerance(impulse_magnitude, impulse_magnitude)?;
        if capacity_residual_ns > tolerance {
            return Err(EllipticImpulseError::Unresolved("impulse ellipse capacity"));
        }
        let mut admitted_endpoint = to_principal(state.axis, after)?;
        for i in 0..2 {
            // Only a zero target from the projection can use the endpoint
            // tolerance. Never infer arrest from an arbitrary small speed.
            if endpoint[i] == 0.0 {
                if admitted_endpoint[i].abs() > controls.velocity_tolerance_m_s {
                    return Err(EllipticImpulseError::Unresolved(
                        "arrested impulse component",
                    ));
                }
                admitted_endpoint[i] = 0.0;
            }
        }
        let active_slip = [0, 1].map(|i| {
            if caps[i] == 0.0 {
                0.0
            } else {
                admitted_endpoint[i]
            }
        });
        let stationarity_residual_ns = if active_slip == [0.0; 2] {
            0.0
        } else {
            let expected = support(caps, admitted_endpoint)?;
            norm([
                actual_reaction[0] - expected[0],
                actual_reaction[1] - expected[1],
            ])?
        };
        if stationarity_residual_ns > tolerance {
            return Err(EllipticImpulseError::Unresolved(
                "final-slip maximum dissipation",
            ));
        }

        let midpoint = [0, 1].map(|i| 0.5 * before[i] + 0.5 * after[i]);
        let dissipated_j = checked(-dot(impulse_ns, midpoint))?;
        // Independent kinetic change: measured endpoint difference, original
        // SPD inverse mobility, and midpoint velocity. This avoids subtracting
        // two much larger kinetic energies in a small-friction impact.
        let delta_velocity = [after[0] - before[0], after[1] - before[1]];
        let g = delta_velocity.map(|v| v / response.scale);
        let [a, b, d] = response.normalized;
        let momentum_change = [
            (d * g[0] - b * g[1]) / response.determinant,
            (a * g[1] - b * g[0]) / response.determinant,
        ];
        for x in momentum_change {
            checked(x)?;
        }
        let kinetic_energy_change_j = checked(dot(midpoint, momentum_change))?;
        let energy_residual_j = checked(kinetic_energy_change_j + dissipated_j)?;
        if dissipated_j < 0.0 || energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(EllipticImpulseError::Unresolved(
                "passive impulse energy balance",
            ));
        }
        poll(cancelled)?;
        Ok(EllipticImpulseReport {
            impulse_ns,
            velocity_after_m_s: after,
            regime: if endpoint == [0.0; 2] {
                FrictionRegime::Sticking
            } else {
                regime
            },
            dissipated_j,
            kinetic_energy_change_j,
            energy_residual_j,
            impulse_velocity_residual_m_s: velocity_residual,
            capacity_residual_ns,
            stationarity_residual_ns,
            iterations,
        })
    }
}

fn capacity_residual(reaction: [f64; 2], caps: [f64; 2]) -> Result<f64, EllipticImpulseError> {
    let mut scaled = [0.0; 2];
    let mut forbidden = [0.0; 2];
    for i in 0..2 {
        if caps[i] == 0.0 {
            forbidden[i] = reaction[i];
        } else {
            scaled[i] = reaction[i] / caps[i];
        }
    }
    let excess = checked((norm(scaled)? - 1.0).max(0.0) * caps[0].max(caps[1]))?;
    Ok(excess.max(norm(forbidden)?))
}
fn checked(x: f64) -> Result<f64, EllipticImpulseError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(EllipticImpulseError::Invalid(
            "unrepresentable impulse response",
        ))
    }
}
fn norm(v: [f64; 2]) -> Result<f64, EllipticImpulseError> {
    checked(det::hypot(v[0], v[1]))
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[0], a[1] * b[1])
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EllipticImpulseError> {
    if cancelled() {
        Err(EllipticImpulseError::Cancelled)
    } else {
        Ok(())
    }
}
