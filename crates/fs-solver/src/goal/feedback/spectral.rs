//! The existing checked spectral inverse closes the complete feedback equation.

use fs_sparse::Csr;
use super::{FeedbackResidualLimits, FeedbackResidualReport, GoalResidualError,
    enclose_feedback, enclose_feedback_inner, schur};
use super::super::inverse::spectral::{
    SpectralError, SpectralInverse, SpectralInverseLimits, SpectralPreparation,
    SpectralStop, prepare_spectral_inverse,
};

/// Retained preparation diagnostics, never detached inverse authority.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralFeedbackDiagnostics {
    stop: SpectralStop,
    shift: Option<f64>,
    defect_upper: Option<f64>,
    coercivity_lower: Option<f64>,
    work_entries: usize,
    peak_storage_entries: usize,
    shift_attempts: usize,
}
impl SpectralFeedbackDiagnostics {
    /// Why sparse preparation succeeded or stopped without a certificate.
    #[must_use]
    pub const fn stop(&self) -> SpectralStop { self.stop }
    /// Independently checked shift, absent when no certificate was produced.
    #[must_use]
    pub const fn shift(&self) -> Option<f64> { self.shift }
    /// Full original-matrix Gram defect, not just the symmetric proposal's defect.
    #[must_use]
    pub const fn defect_upper(&self) -> Option<f64> { self.defect_upper }
    /// Proved positive coercivity, including stored assembly asymmetry.
    #[must_use]
    pub const fn coercivity_lower(&self) -> Option<f64> { self.coercivity_lower }
    /// Preparation visits including the initial diagonal proposal scan.
    /// Zero for an already checked, caller-owned spectral inverse.
    #[must_use]
    pub const fn work_entries(&self) -> usize { self.work_entries }
    /// Peak logical preparation storage; not allocator bytes or peak RSS.
    /// Zero when no new preparation was requested.
    #[must_use]
    pub const fn peak_storage_entries(&self) -> usize { self.peak_storage_entries }
    /// Shifts tried during this preparation; zero when reusing a certificate.
    #[must_use]
    pub const fn shift_attempts(&self) -> usize { self.shift_attempts }

    fn prepared(preparation: &SpectralPreparation, proposal_work: usize) -> Self {
        Self {
            stop: preparation.stop,
            shift: preparation.certificate.as_ref().map(SpectralInverse::shift),
            defect_upper: preparation.certificate.as_ref().map(SpectralInverse::defect_upper),
            coercivity_lower: preparation.certificate.as_ref().map(SpectralInverse::coercivity_lower),
            // The wrapper reserved proposal_work before passing the remaining
            // finite budget to preparation, so this sum cannot overflow.
            work_entries: proposal_work + preparation.work_entries,
            peak_storage_entries: preparation.peak_storage_entries,
            shift_attempts: preparation.shift_attempts,
        }
    }
}

/// Enclose the full coupled error using an already checked spectral inverse.
///
/// All residuals and response-column defects use the exact operator OWNED by
/// `certificate`. There is no separate matrix argument or untrusted scalar
/// inverse norm. The complete state-contraction or port-Schur check must still
/// succeed; a solid certificate alone never certifies the coupled equation.
/// No primal, response, factorization or eigenvalue solve is repeated.
///
/// # Errors
/// Ordinary residual input, resource, arithmetic and cancellation refusals.
#[allow(clippy::too_many_arguments)]
pub fn enclose_affine_feedback_error_with_spectral_inverse(
    certificate: &SpectralInverse, rhs: &[f64], primal: &[f64],
    injection: &Csr, feedback: &Csr, offset: &[f64],
    responses: Option<&[Vec<f64>]>, scaling: Option<&[f64]>,
    limits: FeedbackResidualLimits, mut checkpoint: impl FnMut() -> bool,
) -> Result<FeedbackResidualReport, GoalResidualError> {
    let (mut report, used) = enclose_feedback_inner(
        certificate.matrix(), rhs, primal, injection, feedback, offset, responses,
        scaling, None, Some(certificate), limits, &mut checkpoint,
    )?;
    report.solid_spectral = Some(SpectralFeedbackDiagnostics {
        stop: SpectralStop::Certified, shift: Some(certificate.shift()),
        defect_upper: Some(certificate.defect_upper()),
        coercivity_lower: Some(certificate.coercivity_lower()), work_entries: 0,
        peak_storage_entries: 0, shift_attempts: 0,
    });
    schur::finish_report(certificate.matrix(), feedback, responses, limits, used, report, checkpoint)
}

/// Assess `(A-B C)x=b+B d` with a bounded existing spectral-inverse fallback.
///
/// Comparison dominance runs first. When it already proves the solid inverse,
/// its existing result and Schur fallback are unchanged. Otherwise the minimum
/// diagonal divided by n proposes a shift; the existing sparse Gram verifier
/// independently checks every shift against the ORIGINAL stored matrix.
///
/// The first residual pass, proposal scan, ALL preparation retries, the second
/// residual/response pass and optional Schur check share one verification-work
/// allowance. Independent logical-storage limits include the owned CSR copy.
/// Failed preparation retains its paid-work/stop diagnostics, never zero error.
/// All supplied response errors remain in the full coupled bound.
///
/// # Errors
/// Existing input, allocation and cancellation refusals. A sparse work/storage
/// stop returns an explicit no-bound report, matching `SpectralPreparation`.
#[allow(clippy::too_many_arguments)]
pub fn enclose_affine_feedback_error_with_spectral(
    matrix: &Csr, rhs: &[f64], primal: &[f64],
    injection: &Csr, feedback: &Csr, offset: &[f64],
    responses: Option<&[Vec<f64>]>, scaling: Option<&[f64]>,
    limits: FeedbackResidualLimits, mut spectral: SpectralInverseLimits,
    mut checkpoint: impl FnMut() -> bool,
) -> Result<FeedbackResidualReport, SpectralError> {
    let (mut original, first_work) = enclose_feedback(
        matrix, rhs, primal, injection, feedback, offset, responses, scaling,
        None, limits, &mut checkpoint,
    )?;
    if original.solid_inverse_infinity_upper().is_some() {
        return schur::finish_report(
            matrix, feedback, responses, limits, first_work, original, checkpoint,
        ).map_err(Into::into);
    }
    let n = matrix.nrows();
    let required = first_work.checked_mul(2).and_then(|work| work.checked_add(n))
        .ok_or(GoalResidualError::Allocation)?;
    let remaining = limits.max_verification_entries.checked_sub(required).ok_or(
        GoalResidualError::Limit { field: "spectral feedback verification entries",
            required, allowed: limits.max_verification_entries },
    )?;
    // This is only a proposal, not an assumed spectral lower bound.
    let mut diagonal = f64::INFINITY;
    for i in 0..n {
        if i % 512 == 0 && !checkpoint() { return Err(GoalResidualError::Cancelled.into()); }
        diagonal = diagonal.min(matrix.get(i, i));
    }
    if !checkpoint() { return Err(GoalResidualError::Cancelled.into()); }
    let shift = diagonal / n as f64;
    if !(shift.is_finite() && shift > 0.0) {
        original.solid_spectral = Some(SpectralFeedbackDiagnostics {
            stop: SpectralStop::NotEstablished, shift: None, defect_upper: None,
            coercivity_lower: None, work_entries: n, peak_storage_entries: 0, shift_attempts: 0,
        });
        return Ok(original);
    }
    spectral.system.max_rows = spectral.system.max_rows.min(limits.solid.max_rows);
    spectral.system.max_nonzeros = spectral.system.max_nonzeros.min(limits.solid.max_nonzeros);
    spectral.max_work_entries = spectral.max_work_entries.min(remaining);
    let preparation = prepare_spectral_inverse(matrix, shift, spectral, &mut checkpoint)?;
    let diagnostics = SpectralFeedbackDiagnostics::prepared(&preparation, n);
    let Some(certificate) = preparation.certificate else {
        original.solid_spectral = Some(diagnostics);
        return Ok(original);
    };
    let before_second = first_work.checked_add(diagnostics.work_entries)
        .ok_or(GoalResidualError::Allocation)?;
    let mut second_limits = limits;
    second_limits.max_verification_entries = limits.max_verification_entries
        .checked_sub(before_second).ok_or(GoalResidualError::Allocation)?;
    let (mut report, second_work) = enclose_feedback_inner(
        certificate.matrix(), rhs, primal, injection, feedback, offset, responses,
        scaling, None, Some(&certificate), second_limits, &mut checkpoint,
    )?;
    report.solid_spectral = Some(diagnostics);
    let used = before_second.checked_add(second_work).ok_or(GoalResidualError::Allocation)?;
    schur::finish_report(certificate.matrix(), feedback, responses, limits, used, report, checkpoint)
        .map_err(Into::into)
}
