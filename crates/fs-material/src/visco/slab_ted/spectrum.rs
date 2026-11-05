//! Positive thermal diffusion modes for the existing Maxwell runtime owners.
//!
//! The omitted complex response is bounded analytically using the convex
//! midpoint integral test on the odd thermal spectrum, not a frequency grid.
//! Scalar coefficient rounding is enclosed separately. The rank-one source
//! factor, its projection and runtime arithmetic are outside this scalar bound.

use super::{SlabThermoelasticBending, ViscoError};
use crate::visco::tensor::MAX_PRONY_BRANCHES;

/// One retained insulated-slab diffusion mode, in increasing thermal order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabPronyTerm {
    weight: f64,
    relaxation_time_s: f64,
}

impl SlabPronyTerm {
    /// Positive share of the full thermal bending relaxation tensor.
    /// Retained weights are never renormalized to replace the omitted tail.
    #[must_use]
    pub const fn weight(&self) -> f64 {
        self.weight
    }

    /// Positive relaxation time [s].
    #[must_use]
    pub const fn relaxation_time_s(&self) -> f64 {
        self.relaxation_time_s
    }
}

/// An admitted finite positive spectrum and its normalized complex-error bound.
///
/// For every signed frequency whose absolute value lies in `band_rad_s`,
/// `|H_infinite - sum_j w_j i omega tau_j/(1+i omega tau_j)| <= error_bound`.
/// The scalar functions treat the stored branch coefficients as exact inputs.
/// Error is normalized by the full relaxation tensor, NOT by its possibly
/// zero loss at a particular frequency. For curvature kappa the corresponding
/// quadratic-form error is bounded by `error_bound * (f dot kappa)^2`.
///
/// The finite branch law is causal and passive at all frequencies. Its
/// approximation guarantee belongs only to the declared band. Finite spectra
/// miss some adiabatic stiffness: the tail is disclosed, never reassigned.
/// This object adds no state/history owner; consume its factors and times with
/// the existing Maxwell/modal relaxation runtime.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabPronySpectrum {
    terms: Vec<SlabPronyTerm>,
    band_rad_s: (f64, f64),
    normalized_tail_bound: f64,
    normalized_coefficient_bound: f64,
    normalized_error_bound: f64,
}

impl SlabPronySpectrum {
    /// Retained positive modes, starting with the longest relaxation time.
    #[must_use]
    pub fn terms(&self) -> &[SlabPronyTerm] {
        &self.terms
    }

    /// Closed angular-frequency band [rad/s], interpreted by absolute value.
    #[must_use]
    pub const fn band_rad_s(&self) -> (f64, f64) {
        self.band_rad_s
    }

    /// Outward upper bound on the normalized complex response error.
    #[must_use]
    pub const fn normalized_error_bound(&self) -> f64 {
        self.normalized_error_bound
    }

    /// Analytic omitted-mode contribution, with outward scalar arithmetic.
    #[must_use]
    pub const fn normalized_tail_bound(&self) -> f64 {
        self.normalized_tail_bound
    }

    /// Contribution from rounded retained weights and relaxation times.
    /// Excludes source tensor, projection and runtime evaluation rounding.
    #[must_use]
    pub const fn normalized_coefficient_bound(&self) -> f64 {
        self.normalized_coefficient_bound
    }
}

impl SlabThermoelasticBending {
    /// Retain the smallest thermal prefix meeting an explicit complex-error gate.
    ///
    /// `tolerance` is dimensionless absolute error per full bending relaxation
    /// strength. The band has `0 <= lo < hi < infinity`. The work allowance is
    /// `1 <= max_terms <= MAX_PRONY_BRANCHES` (64), so preparation is bounded.
    ///
    /// With N retained odd modes, convexity gives
    /// `sum_{j=N}^infinity (2j+1)^-p <= (2N)^(1-p)/(2(p-1))`.
    /// Since `|ix/(1+ix)| <= min(1,|x|)`, the omitted response is bounded by
    /// `(96/pi^4) min((2N)^-3/6, omega_hi*tau_0*(2N)^-5/10)`.
    /// This bound covers storage AND loss without sampling or a loss floor.
    ///
    /// The calculation encloses mathematical pi between adjacent binary64
    /// neighbors and propagates outward arithmetic for positive products and
    /// quotients. For rounded times, the exact identity
    /// `sup_omega |F(omega*tau)-F(omega*tau_hat)| =
    /// |tau-tau_hat|/(tau+tau_hat)` supplies the additional coefficient bound.
    /// No tensor/projection or time-integration error is included here.
    ///
    /// # Errors
    /// Invalid band/tolerance/count, unrepresentable modes, or an exhausted
    /// term allowance that does not meet the declared truncation tolerance.
    pub fn lower_to_prony(
        &self,
        band_rad_s: (f64, f64),
        tolerance: f64,
        max_terms: usize,
    ) -> Result<SlabPronySpectrum, ViscoError> {
        let (lo, hi) = band_rad_s;
        if !lo.is_finite()
            || lo < 0.0
            || !hi.is_finite()
            || hi <= lo
            || !tolerance.is_finite()
            || tolerance <= 0.0
            || max_terms == 0
            || max_terms > MAX_PRONY_BRANCHES
        {
            return Err(ViscoError::Parameters {
                what: "slab spectrum requires a finite nonnegative band, positive tolerance and 1..=64 terms",
            });
        }
        let pi = core::f64::consts::PI;
        let pi2_lo = down_mul(pi.next_down(), pi.next_down());
        let pi2_hi = up_mul(pi.next_up(), pi.next_up());
        let coefficient_lo = down_div(96.0, up_mul(pi2_hi, pi2_hi));
        let coefficient_hi = up_div(96.0, down_mul(pi2_lo, pi2_lo));
        let coefficient = 96.0 / ((pi * pi) * (pi * pi));
        if !coefficient_lo.is_finite() || !(coefficient_lo > 0.0) || !coefficient_hi.is_finite() {
            return Err(ViscoError::Parameters {
                what: "thermal diffusion coefficient is unrepresentable",
            });
        }
        let tau0 = self.first_relaxation_time_s();
        let mut terms = Vec::with_capacity(max_terms);
        let mut rounding = 0.0;
        let mut error = f64::INFINITY;
        for index in 0..max_terms {
            // Integers and their squares are exact at the admitted cap.
            let odd = (2 * index + 1) as f64;
            let square = odd * odd;
            let fourth = square * square;
            let weight = coefficient / fourth;
            let weight_lo = down_div(coefficient_lo, fourth);
            let weight_hi = up_div(coefficient_hi, fourth);
            let tau = tau0 / square;
            let tau_lo = down_div(tau0, square);
            let tau_hi = up_div(tau0, square);
            if !(weight.is_finite()
                && weight > 0.0
                && tau.is_finite()
                && tau > 0.0
                && tau_lo > 0.0
                && tau_hi.is_finite())
            {
                return Err(ViscoError::Parameters {
                    what: "retained thermal diffusion weight or time is unrepresentable",
                });
            }
            let weight_error = (weight - weight_lo)
                .next_up()
                .max((weight_hi - weight).next_up());
            let time_error = (tau - tau_lo).next_up().max((tau_hi - tau).next_up());
            // Divide by tau first to avoid overflowing tau+tau_lo. Both
            // quotients are enclosed; tau_lo/tau >=0 and <=1.
            let relative_time_error = up_div(
                up_div(time_error, tau),
                (1.0 + down_div(tau_lo, tau)).next_down(),
            );
            rounding = up_add(
                rounding,
                up_add(weight_error, up_mul(weight_hi, relative_time_error)),
            );
            terms.push(SlabPronyTerm {
                weight,
                relaxation_time_s: tau,
            });

            let n2 = (2 * terms.len()) as f64;
            let inverse = up_div(1.0, n2);
            let cube = up_mul(up_mul(inverse, inverse), inverse);
            let fifth = up_mul(up_mul(cube, inverse), inverse);
            let uniform = up_div(up_mul(coefficient_hi, cube), 6.0);
            let band = up_div(
                up_mul(up_mul(up_mul(coefficient_hi, fifth), tau0), hi),
                10.0,
            );
            let tail = uniform.min(band);
            error = up_add(tail, rounding);
            if error <= tolerance {
                return Ok(SlabPronySpectrum {
                    terms,
                    band_rad_s,
                    normalized_tail_bound: tail,
                    normalized_coefficient_bound: rounding,
                    normalized_error_bound: error,
                });
            }
        }
        Err(ViscoError::TruncationTolerance {
            achieved: error,
            required: tolerance,
            terms: max_terms,
        })
    }
}

// Directed operations on nonnegative finite inputs for the spectrum bounds.
// Positive overflow is a conservative infinity; a subnormal rounded to zero
// is widened upward so no unrepresented tail becomes a zero certificate.
fn up_add(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        return b;
    }
    if b == 0.0 {
        return a;
    }
    (a + b).next_up()
}
fn up_mul(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        return 0.0;
    }
    (a * b).next_up()
}
fn down_mul(a: f64, b: f64) -> f64 {
    (a * b).next_down().max(0.0)
}
fn up_div(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        return 0.0;
    }
    (a / b).next_up()
}
fn down_div(a: f64, b: f64) -> f64 {
    (a / b).next_down().max(0.0)
}
