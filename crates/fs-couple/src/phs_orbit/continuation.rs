//! Original-system frequency response through folds, using fs-orbit's owner.

use super::{
    InputKind, PeriodicPortInput, PhsOrbitAudit, PhsOrbitError, PhsOrbitProblem,
    audit_forced_orbit, bad, forced_problem,
};
use fs_orbit::{ContinuableProblem, ContinuationBudget, ContinuationPoint, continue_branch};
use fs_phs::PortHamiltonian;

/// One accepted nonlinear frequency-response point with original physical work.
#[derive(Debug, Clone)]
pub struct PhsFrequencyResponsePoint {
    /// Native continuation result. `lambda` and `orbit.omega` both hold angular
    /// frequency [rad/s]; coefficients retain the original pHS state ordering.
    /// Residual, AFT sampling defect and arclength step are retained unchanged.
    pub point: ContinuationPoint,
    /// Original Hamiltonian, port supply and resistance loss integrated over
    /// this orbit at the accepted frequency, including every memory coordinate.
    pub audit: PhsOrbitAudit,
}

impl ContinuableProblem for PhsOrbitProblem<'_> {
    fn set_parameter(&mut self, lambda: f64) {
        if let InputKind::Forced { omega_rad_s, .. } = &mut self.input {
            *omega_rad_s = lambda;
        }
    }
}

/// Follow a forced nonlinear resonance branch through frequency folds.
///
/// Starts at `input.omega_rad_s()` and follows fs-orbit's pseudo-arclength
/// branch toward `omega_end_rad_s`. The frequency may reverse at a fold; it is
/// not a monotone frequency grid. The final accepted point may pass the end
/// frequency. The authored DC values, complex harmonic amplitudes and relative
/// phases stay fixed while their shared fundamental changes. Every source
/// evaluation uses the ORIGINAL storage, full dissipation and memory states.
///
/// Each point retains fs-orbit's residual and AFT sampling evidence, plus an
/// independent physical cycle audit. This is a fixed retained-harmonic solve:
/// callers must repeat with larger `budget.hb.harmonics` and `audit_samples` to
/// establish resolution. Small balance or power residuals do not prove harmonic
/// completeness, stability, or the size of a physical hysteresis loop. In
/// particular unstable periodic branches are returned without a Floquet label.
///
/// Work is bounded by the native continuation/Newton/AFT budgets and the
/// 32..=65536 audit samples per accepted point. Like fs-orbit, this synchronous
/// API has no cooperative cancellation hook. No time integration is performed.
///
/// # Errors
/// Invalid input, unresolved authored forcing, original-source refusal, or any
/// native continuation refusal returns an error without a partial branch.
pub fn continue_forced_frequency_response(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    omega_end_rad_s: f64,
    guess_amplitude: f64,
    budget: &ContinuationBudget,
    audit_samples: usize,
) -> Result<Vec<PhsFrequencyResponsePoint>, PhsOrbitError> {
    if !omega_end_rad_s.is_finite()
        || omega_end_rad_s <= 0.0
        || omega_end_rad_s == input.omega_rad_s()
    {
        return Err(bad(
            "positive finite distinct continuation endpoint required",
        ));
    }
    if !guess_amplitude.is_finite() || guess_amplitude <= 0.0 {
        return Err(bad("positive finite first-state seed amplitude required"));
    }
    if budget.max_steps == 0
        || !budget.initial_step.is_finite()
        || budget.initial_step <= 0.0
        || !budget.min_step.is_finite()
        || budget.min_step <= 0.0
        || budget.min_step > budget.initial_step
    {
        return Err(bad(
            "positive finite continuation steps and nonzero step budget required",
        ));
    }
    if budget.hb.harmonics == 0 || input.highest_nonzero_harmonic() > budget.hb.harmonics {
        return Err(bad(
            "nonzero forcing harmonic exceeds the retained HB truncation",
        ));
    }
    if !(32..=65_536).contains(&audit_samples) {
        return Err(bad("cycle audit needs 32..=65536 samples"));
    }
    let mut problem = forced_problem(system, input)?;
    let result = continue_branch(
        &mut problem,
        |omega| omega,
        input.omega_rad_s(),
        omega_end_rad_s,
        guess_amplitude,
        budget,
    );
    if problem.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    let points = result.map_err(PhsOrbitError::Orbit)?;
    let mut audit_input = input.clone();
    let mut response = Vec::with_capacity(points.len());
    for point in points {
        if !point.lambda.is_finite()
            || point.lambda <= 0.0
            || point.lambda.to_bits() != point.orbit.omega.to_bits()
        {
            return Err(bad(
                "continuation returned an inconsistent forcing frequency",
            ));
        }
        audit_input.omega_rad_s = point.lambda;
        let audit = audit_forced_orbit(system, &audit_input, &point.orbit, audit_samples)?;
        response.push(PhsFrequencyResponsePoint { point, audit });
    }
    Ok(response)
}

#[cfg(test)]
mod tests;
