//! Sampling-confidence assessment of a fixed compliance event.
//!
//! Reuses fs-eproc's Gaussian-mixture confidence sequence on the indicators
//! `QoI <= threshold`, rather than treating the descriptive Monte Carlo standard
//! error as an optional-stopping bound. No physical-model error is bounded here.

use core::fmt;

use fs_eproc::GaussianMixtureCs;

use crate::{AnytimeEstimate, UqExecution, UqStatus};

/// An execution cannot support the requested compliance-confidence assessment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UqComplianceError {
    /// Alpha must be finite, strictly inside (0,1), with a finite reciprocal.
    InvalidAlpha,
    /// The requested radius must be finite and nonnegative.
    InvalidHalfWidth,
    /// The admitted plan did not declare a compliance ceiling.
    MissingThreshold,
    /// A failed run's accepted prefix must not masquerade as an uncensored sample.
    RefusedExecution,
    /// The underlying numerical confidence-sequence calculation was non-finite.
    NumericalRange,
}

impl fmt::Display for UqComplianceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::InvalidAlpha => "alpha must be inside (0,1) with a finite reciprocal",
            Self::InvalidHalfWidth => "target half-width must be finite and nonnegative",
            Self::MissingThreshold => "the UQ plan must declare a compliance threshold",
            Self::RefusedExecution => "a refused UQ execution cannot report compliance confidence",
            Self::NumericalRange => "compliance confidence exceeds the finite numerical range",
        })
    }
}

impl std::error::Error for UqComplianceError {}

impl UqExecution {
    /// Assess `P(QoI <= threshold)` using EVERY accepted observation so far.
    /// Returns `None` before the first observation; no evaluator is invoked.
    ///
    /// Indicators lie in [0,1], so the existing Gaussian-mixture confidence
    /// sequence uses sigma=1/2 and rho=1. Under its fixed-conditional-mean
    /// sampling assumptions, its mathematical interval is valid uniformly over
    /// sample counts at the supplied alpha. Choose alpha before inspecting the
    /// run and keep the model and threshold fixed. Searching across seeds, models,
    /// or confidence levels needs separate multiplicity control. The numerical
    /// implementation is not an outward-rounded interval certificate.
    ///
    /// `converged` means the UNCLIPPED confidence radius is at most the requested
    /// `target_half_width`; clipping to [0,1] never creates artificial precision.
    /// Calling this between `advance` chunks supports stopping decisions without
    /// confusing empirical standard error with a confidence sequence. This does
    /// not change execution status: `Complete` still means the original sample
    /// budget was exhausted, not that a statistical target was reached.
    ///
    /// Replays the retained indicators in O(n) time and O(1) extra storage, so
    /// checkpoint/restored and uninterrupted prefixes give the same assessment.
    /// It bounds sampling uncertainty for the DECLARED model only, not numerical,
    /// geometric, material-card, or physical-model error; evidence stays Estimated.
    ///
    /// # Errors
    /// Refuses invalid confidence parameters, missing compliance thresholds,
    /// failed executions, or a non-finite confidence-sequence calculation.
    pub fn assess_compliance(
        &self,
        alpha: f64,
        target_half_width: f64,
    ) -> Result<Option<AnytimeEstimate>, UqComplianceError> {
        if !alpha.is_finite() || alpha <= 0.0 || alpha >= 1.0 || !(1.0 / alpha).is_finite() {
            return Err(UqComplianceError::InvalidAlpha);
        }
        if !target_half_width.is_finite() || target_half_width < 0.0 {
            return Err(UqComplianceError::InvalidHalfWidth);
        }
        let threshold = self
            .plan
            .compliance_threshold
            .ok_or(UqComplianceError::MissingThreshold)?;
        if self.status == UqStatus::Refused || self.failure.is_some() {
            return Err(UqComplianceError::RefusedExecution);
        }
        let mut confidence = GaussianMixtureCs::new(0.5, 1.0, alpha);
        for &value in &self.values {
            confidence.observe(if value <= threshold { 1.0 } else { 0.0 });
        }
        let Some((center, radius)) = confidence.interval() else {
            return Ok(None);
        };
        if !center.is_finite() || !radius.is_finite() || radius < 0.0 {
            return Err(UqComplianceError::NumericalRange);
        }
        Ok(Some(AnytimeEstimate {
            mean: center,
            lo: (center - radius).max(0.0),
            hi: (center + radius).min(1.0),
            n: confidence.len(),
            converged: radius <= target_half_width,
        }))
    }
}
