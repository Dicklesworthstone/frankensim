//! Exact transient derivatives of the original Maxwell ramp, including heat.
//!
//! The logarithmic parameter is `ln(tau)` at fixed physical step duration.
//! At fixed strain increment, differentiating `ln(dt)` instead has the opposite
//! sign. These are local constitutive derivatives: a caller differentiating a
//! history must propagate the initial branch-strain derivative through every
//! accepted step.

use super::PronyRamp;

/// Original Maxwell ramp response and its exact local Jacobian.
///
/// Rows are endpoint branch strain, interval-mean branch strain, and physical
/// heat per unit branch modulus. Columns are initial branch strain `q0`, total
/// strain increment `delta`, and `ln(tau)` at fixed step duration. Multiply
/// mean/endpoint strain by the original modulus for stresses, and heat by the
/// original modulus for energy density. Work is mean stress times `delta`;
/// differentiating work must also differentiate that final multiplication.
///
/// The last column holds the initial history fixed. It is not the derivative
/// of a complete experiment unless the retained history is propagated with
/// [`Self::push_forward`] or [`Self::pull_back`]. No thermal evolution, source
/// extrapolation, fractional-fit derivative, or material-identification claim
/// follows from this local constitutive derivative.
#[derive(Debug, Clone, Copy)]
pub struct PronyRampLinearization {
    /// The bit-identical tuple returned by the original [`PronyRamp::apply`].
    pub response: (f64, f64, f64),
    /// Jacobian in the row/column order documented above.
    pub jacobian: [[f64; 3]; 3],
}

impl PronyRampLinearization {
    /// Propagate `[dq0, ddelta, dln(tau)]` through the exact ramp.
    #[must_use]
    #[inline]
    pub fn push_forward(&self, direction: [f64; 3]) -> [f64; 3] {
        self.jacobian
            .map(|row| row[0] * direction[0] + row[1] * direction[1] + row[2] * direction[2])
    }

    /// Pull back seeds for endpoint strain, interval mean, and physical heat.
    ///
    /// The result is the cotangent for `[q0, delta, ln(tau)]`. Accumulate the
    /// last entry over a branch's retained history; feed the first entry into
    /// the preceding step's endpoint-strain seed. This uses the transpose of
    /// the same Jacobian as [`Self::push_forward`], without an adjoint solve,
    /// allocation, finite differences, or a second constitutive history.
    #[must_use]
    #[inline]
    pub fn pull_back(&self, seed: [f64; 3]) -> [f64; 3] {
        core::array::from_fn(|column| {
            self.jacobian[0][column] * seed[0]
                + self.jacobian[1][column] * seed[1]
                + self.jacobian[2][column] * seed[2]
        })
    }
}

/// Reusable exact ramp coefficients with opt-in relaxation-time derivatives.
///
/// Prepare once per original branch and time step, then reuse for every
/// independent strain coordinate. This owns the original [`PronyRamp`]; its
/// primal result is unchanged. Derivatives require only scalar arithmetic and
/// cancellation-free small-rate series at preparation. Application performs no
/// allocation, exponentials, or finite differences. Forward-only consumers of
/// [`PronyRamp`] pay no additional storage or arithmetic cost.
///
/// The caller admits `r = dt/tau >= 0`, finite strains, and representable
/// responses and derivatives, as for [`PronyRamp`]. The exact elastic `r = 0`
/// and fully relaxed `r = +inf` limits have zero log-time derivatives.
#[derive(Debug, Clone, Copy)]
pub struct PronyRampSensitivity {
    ramp: PronyRamp,
    decay_log_tau: f64,
    ramp_log_tau: f64,
    mean_log_tau: f64,
    heat_shift_log_tau: f64,
    heat_scale_log_tau: f64,
    heat_remainder_log_tau: f64,
}

impl PronyRampSensitivity {
    /// Prepare the original ramp and its derivative with respect to `ln(tau)`
    /// at fixed physical duration, including both zero and infinite `r`.
    #[must_use]
    pub fn new(r: f64) -> Self {
        let ramp = PronyRamp::new(r);
        if r == 0.0 || r == f64::INFINITY {
            return Self {
                ramp,
                decay_log_tau: 0.0,
                ramp_log_tau: 0.0,
                mean_log_tau: 0.0,
                heat_shift_log_tau: 0.0,
                heat_scale_log_tau: 0.0,
                heat_remainder_log_tau: 0.0,
            };
        }

        // With r = dt/tau, d/d ln(tau) = -r d/dr. Reuse the SAME
        // exponentials and positive heat factorization as the original ramp.
        let decay_log_tau = r * ramp.decay;
        let heat_scale_log_tau = -ramp.decay * decay_log_tau;
        let (ramp_log_tau, mean_log_tau, heat_shift_log_tau, heat_remainder_log_tau) = if r < 0.1 {
            // b-a and c-(b-a)/r cancel precisely where weak damping must
            // remain differentiable. These analytic series retain their
            // leading O(r) terms even when exp(-r) rounds to one.
            let ramp_log_tau = r
                * (0.5
                    + r * (-1.0 / 3.0
                        + r * (1.0 / 8.0
                            + r * (-1.0 / 30.0
                                + r * (1.0 / 144.0
                                    + r * (-1.0 / 840.0
                                        + r * (1.0 / 5760.0
                                            + r * (-1.0 / 45360.0 + r / 403200.0))))))));
            let mean_log_tau = r
                * (1.0 / 6.0
                    + r * (-1.0 / 12.0
                        + r * (1.0 / 40.0
                            + r * (-1.0 / 180.0
                                + r * (1.0 / 1008.0
                                    + r * (-1.0 / 6720.0 + r * (1.0 / 51840.0 - r / 453600.0)))))));
            let r2 = r * r;
            let heat_shift_log_tau = r2
                * (1.0 / 12.0
                    + r2 * (-1.0 / 60.0
                        + r2 * (17.0 / 6720.0
                            + r2 * (-31.0 / 90720.0 + r2 * (691.0 / 15966720.0)))));
            let heat_remainder_log_tau = r
                * (-1.0 / 12.0
                    + r2 * (1.0 / 40.0
                        + r2 * (-17.0 / 4032.0
                            + r2 * (31.0 / 51840.0 - r2 * (691.0 / 8870400.0)))));
            (
                ramp_log_tau,
                mean_log_tau,
                heat_shift_log_tau,
                heat_remainder_log_tau,
            )
        } else {
            let ramp_log_tau = ramp.ramp - ramp.decay;
            let denominator = 1.0 + ramp.decay;
            let heat_shift_log_tau =
                ramp.heat_shift - 2.0 * ramp.decay / (denominator * denominator);
            (
                ramp_log_tau,
                ramp.mean_ramp - ramp_log_tau / r,
                heat_shift_log_tau,
                ramp.heat_remainder - 2.0 * heat_shift_log_tau / r,
            )
        };
        Self {
            ramp,
            decay_log_tau,
            ramp_log_tau,
            mean_log_tau,
            heat_shift_log_tau,
            heat_scale_log_tau,
            heat_remainder_log_tau,
        }
    }

    /// The original immutable preparation, for reuse without derivatives.
    #[must_use]
    pub const fn ramp(&self) -> &PronyRamp {
        &self.ramp
    }

    /// Original endpoint, mean and physical heat, with their exact Jacobian.
    ///
    /// Heat is differentiated through its independent positive quadratic
    /// integral, never through a subtraction of work and stored energy.
    #[must_use]
    #[inline]
    pub fn apply(&self, q0: f64, delta: f64) -> PronyRampLinearization {
        let shifted = q0 + self.ramp.heat_shift * delta;
        let heat_q0 = 2.0 * (self.ramp.heat_scale * shifted);
        let heat_delta = heat_q0 * self.ramp.heat_shift + 2.0 * (self.ramp.heat_remainder * delta);
        let heat_log_tau = self.heat_scale_log_tau * shifted * shifted
            + (heat_q0 * self.heat_shift_log_tau) * delta
            + (self.heat_remainder_log_tau * delta) * delta;
        PronyRampLinearization {
            response: self.ramp.apply(q0, delta),
            jacobian: [
                [
                    self.ramp.decay,
                    self.ramp.ramp,
                    self.decay_log_tau * q0 + self.ramp_log_tau * delta,
                ],
                [
                    self.ramp.ramp,
                    self.ramp.mean_ramp,
                    self.ramp_log_tau * q0 + self.mean_log_tau * delta,
                ],
                [heat_q0, heat_delta, heat_log_tau],
            ],
        }
    }
}
