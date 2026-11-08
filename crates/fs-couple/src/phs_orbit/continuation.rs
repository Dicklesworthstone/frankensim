//! Original-system frequency response through folds, using fs-orbit's owner.

use super::{
    InputKind, PeriodicPortInput, PhsOrbitAudit, PhsOrbitError, PhsOrbitProblem,
    audit_forced_orbit, bad, forced_problem,
};
use fs_orbit::{
    ContinuableProblem, ContinuationBudget, ContinuationPoint, ContinuationScales,
    ForcedParameterJacobian, ForcedParameterPoint, HbScales, ScaledContinuationPoint,
    continue_branch, continue_branch_scaled_seeded,
};
use fs_phs::PortHamiltonian;

mod response;
pub use response::{
    PhsForcedBranchStabilityControls, PhsForcedPortBranchControls, PhsForcedPortBranchError,
    PhsForcedPortBranchPoint, continue_forced_port_response_analytic_scaled_seeded,
    continue_forced_port_response_scaled_seeded,
};

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

/// A physically scaled forced branch point and original-system cycle work.
#[derive(Debug, Clone)]
pub struct PhsScaledFrequencyResponsePoint {
    /// Physical Fourier waveform/frequency and normalized continuation evidence.
    pub point: ScaledContinuationPoint,
    /// Original Hamiltonian, signed supply and all resistance/memory loss.
    pub audit: PhsOrbitAudit,
}

impl ContinuableProblem for PhsOrbitProblem<'_> {
    fn set_parameter(&mut self, lambda: f64) {
        if let InputKind::Forced { omega_rad_s, .. } = &mut self.input {
            *omega_rad_s = lambda;
        }
    }

    fn forced_parameter_jacobian(&self) -> Option<ForcedParameterJacobian<Self>> {
        matches!(self.input, InputKind::Forced { .. }).then_some(Self::frequency_parameter_jacobian)
    }
}

impl PhsOrbitProblem<'_> {
    /// At fixed phase and Fourier coefficients the entire original vector
    /// field is frequency independent: u(theta) keeps its authored complex
    /// coefficients, and every material-memory variable remains an ODE state.
    /// Only the original sI port changes, hence dR_n/domega = i*n*X_n.
    fn frequency_parameter_jacobian(
        &self,
        point: ForcedParameterPoint<'_>,
        out: &mut [f64],
    ) -> bool {
        let InputKind::Forced { omega_rad_s, .. } = &self.input else {
            return false;
        };
        // These private adapters are only used with omega_of(lambda)=lambda.
        // Do not silently reuse this tangent under another clock/source law.
        if !omega_rad_s.is_finite()
            || *omega_rad_s <= 0.0
            || point.parameter.to_bits() != omega_rad_s.to_bits()
            || point.angular_frequency.to_bits() != omega_rad_s.to_bits()
        {
            return false;
        }
        let d = self.system.state_dim();
        let Some(expected) = point
            .harmonics
            .checked_mul(2)
            .and_then(|orders| orders.checked_add(1))
            .and_then(|orders| d.checked_mul(orders))
        else {
            return false;
        };
        if point.coefficients.len() != expected
            || out.len() != expected
            || point.coefficients.iter().any(|value| !value.is_finite())
        {
            return false;
        }
        out.fill(0.0);
        for harmonic in 1..=point.harmonics {
            let real = d + 2 * (harmonic - 1) * d;
            let imaginary = real + d;
            for state in 0..d {
                out[real + state] = -(harmonic as f64) * point.coefficients[imaginary + state];
                out[imaginary + state] = harmonic as f64 * point.coefficients[real + state];
            }
        }
        out.iter().all(|value| value.is_finite())
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
/// Its continuation frequency column is analytic in this fixed-phase chart:
/// the complete original ODE gives `dR_n/domega = i*n*X_n`, including every
/// memory coordinate, without perturbing or replaying the physical source.
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
    admit_response(
        input,
        omega_end_rad_s,
        Some(guess_amplitude),
        budget,
        audit_samples,
    )?;
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
        let audit = audit_point(system, &mut audit_input, &point, audit_samples)?;
        response.push(PhsFrequencyResponsePoint { point, audit });
    }
    Ok(response)
}

/// Continue a full physical pHS waveform with independently declared units.
///
/// `seed` contains every original state and retained harmonic in canonical
/// physical packing, including material memory and forced relative phases.
/// The original equations are solved again at the starting frequency. State
/// and equation scales are those of the original pHS ODE; the frequency-column
/// and arclength parameter use `scales.angular_frequency`. Step budgets measure
/// dimensionless arclength. No source amplitude, Hamiltonian or loss is rescaled.
///
/// Each accepted point retains the native normalized HB/AFT/arclength defects
/// plus the same physical cycle audit as the existing frequency-response path.
/// Fixed harmonic truncation, possible Newton branch changes, numerical fold
/// resolution and lack of a stability guarantee remain disclosed limitations.
///
/// # Errors
/// The original response admission/source refusals and native scaled
/// continuation errors. No incomplete branch is returned on refusal.
#[allow(clippy::too_many_arguments)]
pub fn continue_forced_frequency_response_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    omega_end_rad_s: f64,
    seed: Vec<f64>,
    budget: &ContinuationBudget,
    scales: &HbScales,
    audit_samples: usize,
) -> Result<Vec<PhsScaledFrequencyResponsePoint>, PhsOrbitError> {
    admit_response(input, omega_end_rad_s, None, budget, audit_samples)?;
    let mut problem = forced_problem(system, input)?;
    let result = continue_branch_scaled_seeded(
        &mut problem,
        |omega| omega,
        input.omega_rad_s(),
        omega_end_rad_s,
        seed,
        budget,
        &ContinuationScales {
            hb: scales.clone(),
            parameter: scales.angular_frequency,
        },
    );
    if problem.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    let points = result.map_err(PhsOrbitError::Orbit)?;
    let mut audit_input = input.clone();
    let mut response = Vec::with_capacity(points.len());
    for point in points {
        let audit = audit_point(system, &mut audit_input, &point.point, audit_samples)?;
        response.push(PhsScaledFrequencyResponsePoint { point, audit });
    }
    Ok(response)
}

fn admit_response(
    input: &PeriodicPortInput,
    omega_end_rad_s: f64,
    guess_amplitude: Option<f64>,
    budget: &ContinuationBudget,
    audit_samples: usize,
) -> Result<(), PhsOrbitError> {
    if !omega_end_rad_s.is_finite()
        || omega_end_rad_s <= 0.0
        || omega_end_rad_s == input.omega_rad_s()
    {
        return Err(bad(
            "positive finite distinct continuation endpoint required",
        ));
    }
    if guess_amplitude.is_some_and(|value| !value.is_finite() || value <= 0.0) {
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
    Ok(())
}

fn audit_point(
    system: &PortHamiltonian,
    input: &mut PeriodicPortInput,
    point: &ContinuationPoint,
    samples: usize,
) -> Result<PhsOrbitAudit, PhsOrbitError> {
    if !point.lambda.is_finite()
        || point.lambda <= 0.0
        || point.lambda.to_bits() != point.orbit.omega.to_bits()
    {
        return Err(bad(
            "continuation returned an inconsistent forcing frequency",
        ));
    }
    input.omega_rad_s = point.lambda;
    audit_forced_orbit(system, input, &point.orbit, samples)
}

#[cfg(test)]
mod tests;
