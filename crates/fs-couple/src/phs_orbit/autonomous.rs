//! Self-excited periodic regimes driven through original physical pHS ports.

use super::{
    InputKind, PhsOrbitAudit, PhsOrbitError, PhsOrbitProblem, audit_orbit_with_input, bad,
};
use core::cell::{Cell, RefCell};
use fs_orbit::{
    AdaptiveHbOrbit, HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, solve_hb, solve_hb_adaptive,
};
use fs_phs::PortHamiltonian;

/// Original physical-port feedback families on the shared continuation owner.
pub mod continuation;

mod shooting;
pub use shooting::{
    PhsAutonomousControls, PhsAutonomousLevel, PhsAutonomousStability,
    analyze_autonomous_feedback_stability,
};

/// Caller-authored, instantaneous autonomous feedback on physical pHS ports.
///
/// `flows` is the original power-conjugate output G^T grad H. Fill the complete
/// physical effort vector `out`, so supplied power is exactly out dot flows.
/// The law may supply or absorb energy; component passivity does not establish
/// passivity of this externally driven loop. No material or excitation law is
/// inferred here. Return false for a source-domain refusal; unwritten/nonfinite
/// outputs also refuse. Calls must be deterministic and read-only: controller,
/// transport or constitutive memory belongs in explicit original pHS states,
/// never hidden callback history. There is no externally prescribed clock.
pub trait AutonomousPortFeedback {
    /// Exactly the number of original physical input/output ports.
    fn port_count(&self) -> usize;
    /// Evaluate all efforts from the full physical state and conjugate flows.
    fn input(&self, state: &[f64], flows: &[f64], out: &mut [f64]) -> bool;
}

pub(super) fn feedback_problem<'a>(
    system: &'a PortHamiltonian,
    feedback: &'a dyn AutonomousPortFeedback,
) -> Result<PhsOrbitProblem<'a>, PhsOrbitError> {
    if system.state_dim() == 0 || system.port_dim() != feedback.port_count() {
        return Err(bad(
            "autonomous feedback must cover every original physical port",
        ));
    }
    Ok(PhsOrbitProblem {
        system,
        input: InputKind::Feedback {
            law: feedback,
            flows: RefCell::new(vec![0.0; system.port_dim()]),
        },
        // The complete feedback derivative is not declared by this law.
        tangent: super::tangent::StorageTangent::Unavailable,
        effort: RefCell::new(vec![0.0; system.state_dim()]),
        ports: RefCell::new(vec![0.0; system.port_dim()]),
        source_refused: Cell::new(false),
    })
}

fn admit_hb(omega: f64, amplitude: f64, budget: &HbBudget) -> Result<(), PhsOrbitError> {
    if !omega.is_finite()
        || omega <= 0.0
        || !amplitude.is_finite()
        || amplitude <= 0.0
        || budget.harmonics == 0
        || budget.max_newton == 0
        || !budget.tolerance.is_finite()
        || budget.tolerance <= 0.0
    {
        return Err(bad(
            "positive finite autonomous seeds and nonzero HB budgets required",
        ));
    }
    Ok(())
}

/// Solve a self-excited periodic regime using fs-orbit's original autonomous HB
/// owner. Frequency and amplitude emerge from the complete original ODE and
/// the explicit feedback law; no negative resistance is inserted into the pHS.
/// A constant passive bias does not automatically sustain a dissipative cycle.
///
/// # Errors
/// Invalid seeds/budgets/ports, original source refusal, or native HB failure.
pub fn solve_autonomous_feedback_harmonic_balance(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    omega_guess: f64,
    guess_amplitude: f64,
    budget: &HbBudget,
) -> Result<HbOrbit, PhsOrbitError> {
    admit_hb(omega_guess, guess_amplitude, budget)?;
    let problem = feedback_problem(system, feedback)?;
    let result = solve_hb(
        &problem,
        HbAnchor::Autonomous { omega_guess },
        guess_amplitude,
        budget,
    );
    if problem.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

/// Autonomous feedback with the native retained/omitted-harmonic and successive
/// waveform/frequency refinement checks. Original source work remains physical;
/// sampled Fourier resolution is not an infinite-spectrum or existence proof.
///
/// # Errors
/// Invalid settings/ports, source refusal, or native solve/refinement failure.
pub fn solve_autonomous_feedback_harmonic_balance_adaptive(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    omega_guess: f64,
    guess_amplitude: f64,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    admit_hb(omega_guess, guess_amplitude, budget)?;
    let problem = feedback_problem(system, feedback)?;
    let result = solve_hb_adaptive(
        &problem,
        HbAnchor::Autonomous { omega_guess },
        guess_amplitude,
        budget,
        refinement,
    );
    if problem.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

/// Original Hamiltonian, supplied feedback work and native resistance loss on
/// the supplied HB waveform. This is the same physical quadrature used for
/// forced pHS response, with the authored law evaluated at every sample.
///
/// # Errors
/// Invalid ports/waveform/quadrature, source refusal, or nonfinite physical work.
pub fn audit_autonomous_feedback_orbit(
    system: &PortHamiltonian,
    feedback: &dyn AutonomousPortFeedback,
    candidate: &HbOrbit,
    samples: usize,
) -> Result<PhsOrbitAudit, PhsOrbitError> {
    if feedback.port_count() != system.port_dim() {
        return Err(bad(
            "autonomous feedback must cover every original physical port",
        ));
    }
    audit_orbit_with_input(system, candidate, samples, |_, state, flows, out| {
        feedback.input(state, flows, out)
    })
}

#[cfg(test)]
mod tests;
