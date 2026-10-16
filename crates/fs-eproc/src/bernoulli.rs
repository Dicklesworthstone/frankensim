//! Bernoulli likelihood-mixture confidence sequences, not posterior intervals.
//!
//! With s successes and f failures, integrate the sequence likelihood against
//! a FIXED Beta(1/2, 1/2) mixing distribution. For each null probability p,
//! E(p) = B(s+1/2,f+1/2) / (B(1/2,1/2) * p^s * (1-p)^f).
//! Under Bernoulli observations with fixed conditional success probability p,
//! this is a nonnegative test martingale (a supermartingale at p=0 or p=1).
//! Ville's inequality gives simultaneous coverage by {p: E(p) < 1/alpha}.
//! This is the conjugate-mixture construction of Howard et al. (2021),
//! "Time-uniform, nonparametric, nonasymptotic confidence sequences", section 3.2.
//!
//! The mixing distribution is an inference tuning choice, NOT an assertion
//! that the unknown physical probability is random. It cannot be selected after
//! inspecting data. No iid-point interpretation of QMC, continuous scores,
//! missing/selected outcomes, model error or outward-rounded proof is supplied.
use core::fmt;
use fs_math::det;

/// Numerical reference envelope, matching the bounded product-UQ sample cap.
/// A caller may stop sooner. This is not an assertion about unlimited horizons.
pub const MAX_BERNOULLI_SAMPLES: u64 = 1_000_000;

/// Invalid inference input or exhaustion; none changes an accepted process.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BernoulliCsError {
    /// Alpha must be finite and strictly inside (0,1).
    InvalidAlpha,
    /// A queried null probability must be finite and in [0,1].
    InvalidProbability,
    /// The explicitly documented numerical sample envelope was exhausted.
    SampleLimit,
    /// Nonfinite or inconsistent likelihood arithmetic; no interval is returned.
    NumericalRange,
}
impl fmt::Display for BernoulliCsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidAlpha => "Bernoulli confidence alpha must be finite and inside (0,1)",
            Self::InvalidProbability => "Bernoulli null probability must be finite and in [0,1]",
            Self::SampleLimit => "Bernoulli mixture sample envelope exhausted",
            Self::NumericalRange => "Bernoulli mixture arithmetic cannot produce a finite interval",
        })
    }
}
impl std::error::Error for BernoulliCsError {}

/// Numerical likelihood-mixture confidence bounds. They are generally asymmetric
/// about the empirical mean; neither field is a radius or a posterior quantile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BernoulliInterval {
    /// Empirical success frequency, not the posterior mean.
    pub mean: f64,
    /// Lower probability bound.
    pub lo: f64,
    /// Upper probability bound.
    pub hi: f64,
    /// Complete binary observations consumed.
    pub n: u64,
}

/// Constant-storage, clone-resumable Beta(1/2,1/2) likelihood mixture.
///
/// Updates consume `bool`, so fractional outcomes cannot accidentally inherit
/// Bernoulli validity. Alpha and the prior are fixed at construction. Replaying
/// the same ordered observations preserves the numerical state; exact real
/// arithmetic is exchangeable, but bitwise permutation invariance is not claimed.
#[derive(Clone, Debug, PartialEq)]
pub struct BernoulliMixtureCs {
    log_threshold: f64,
    n: u64,
    successes: u64,
    log_marginal: f64,
    correction: f64,
}
impl BernoulliMixtureCs {
    /// Admit alpha before any observations. No reciprocal is formed, so small
    /// positive alpha does not overflow 1/alpha. Logs use fs-math strict kernels.
    pub fn new(alpha: f64) -> Result<Self, BernoulliCsError> {
        if !alpha.is_finite() || alpha <= 0.0 || alpha >= 1.0 {
            return Err(BernoulliCsError::InvalidAlpha);
        }
        let log_threshold = -det::ln(alpha);
        if !log_threshold.is_finite() || log_threshold <= 0.0 {
            return Err(BernoulliCsError::NumericalRange);
        }
        Ok(Self { log_threshold, n: 0, successes: 0, log_marginal: 0.0, correction: 0.0 })
    }

    /// Consume a completed binary observation. The beta predictive recurrence
    /// needs one logarithm and a compensated sum, not factorials or gamma fits.
    /// All fallible arithmetic is staged before mutating the accepted state.
    pub fn observe(&mut self, success: bool) -> Result<(), BernoulliCsError> {
        if self.n >= MAX_BERNOULLI_SAMPLES { return Err(BernoulliCsError::SampleLimit); }
        let count = if success { self.successes } else { self.n - self.successes };
        let increment = det::ln((count as f64 + 0.5) / (self.n as f64 + 1.0));
        let adjusted = increment - self.correction;
        let next = self.log_marginal + adjusted;
        let correction = (next - self.log_marginal) - adjusted;
        if !next.is_finite() || !correction.is_finite() || next > 0.0 {
            return Err(BernoulliCsError::NumericalRange);
        }
        self.log_marginal = next;
        self.correction = correction;
        self.n += 1;
        self.successes += u64::from(success);
        Ok(())
    }

    /// Log e-value for a fixed null. Impossible null endpoints return +infinity;
    /// before data every admitted null has log e-value zero. No state changes.
    pub fn log_e_value(&self, probability: f64) -> Result<f64, BernoulliCsError> {
        if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
            return Err(BernoulliCsError::InvalidProbability);
        }
        let failures = self.n - self.successes;
        if (probability == 0.0 && self.successes != 0) || (probability == 1.0 && failures != 0) {
            return Ok(f64::INFINITY);
        }
        let mut value = self.log_marginal;
        if self.successes != 0 { value -= self.successes as f64 * det::ln(probability); }
        if failures != 0 { value -= failures as f64 * det::ln(1.0 - probability); }
        if !value.is_finite() { return Err(BernoulliCsError::NumericalRange); }
        Ok(value)
    }

    /// Complete numerical confidence interval; None before data, never zero
    /// uncertainty merely because every observation agreed. Invert each monotone
    /// branch around s/n with at most 80 bisections and retain the OUTER endpoints
    /// of the computed brackets. Early precision exhaustion only widens them.
    /// Mathematical time-uniform coverage requires the stated sampling model;
    /// deterministic floating-point logs/updates are not interval certificates.
    pub fn interval(&self) -> Result<Option<BernoulliInterval>, BernoulliCsError> {
        if self.n == 0 { return Ok(None); }
        let mean = self.successes as f64 / self.n as f64;
        if self.log_e_value(mean)? >= self.log_threshold { return Err(BernoulliCsError::NumericalRange); }
        let mut lo = 0.0;
        let mut hi = 1.0;
        if self.successes != 0 { lo = self.outer_bound(0.0, mean)?; }
        if self.successes != self.n { hi = self.outer_bound(1.0, mean)?; }
        Ok(Some(BernoulliInterval { mean, lo, hi, n: self.n }))
    }

    fn outer_bound(&self, mut outside: f64, mut inside: f64) -> Result<f64, BernoulliCsError> {
        for _ in 0..80 {
            let mid = f64::midpoint(outside, inside);
            if mid == outside || mid == inside { break; }
            if self.log_e_value(mid)? >= self.log_threshold { outside = mid; }
            else { inside = mid; }
        }
        Ok(outside)
    }

    /// Completed observations, including both successes and failures.
    #[must_use]
    pub const fn len(&self) -> u64 { self.n }
    /// True before the first complete observation.
    #[must_use]
    pub const fn is_empty(&self) -> bool { self.n == 0 }
    /// Completed successes; failures are len() - successes().
    #[must_use]
    pub const fn successes(&self) -> u64 { self.successes }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numerical_envelope_refuses_without_mutating_the_checkpoint() {
        let mut state = BernoulliMixtureCs::new(0.05).unwrap();
        state.n = MAX_BERNOULLI_SAMPLES;
        let before = state.clone();
        assert_eq!(state.observe(false), Err(BernoulliCsError::SampleLimit));
        assert_eq!(state, before);
    }
}
