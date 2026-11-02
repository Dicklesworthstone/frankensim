//! Continuous, possibly inharmonic forcing of one complete coupled modal owner.
//!
//! Springs, dashpots, material histories and reciprocal acoustic coordinates
//! evolve together. The cached source propagator evaluates actual force work
//! and each original positive loss channel over the same interval.

use super::exact::tones::PreparedJointTones;
use super::{CoupledModalFrame, CoupledModalSystem, ModalCouplingError, poll};
use fs_exec::CancelGate;

/// Explicit cold dense resources for the complete state and forcing generator.
/// Existing mechanical/setup/pressure/energy limits remain binding.
pub use fs_phs::linear_exact::DenseExponentialBudget as CoupledToneBudget;

/// Owning fixed-clock runtime for a finite sum of continuous generalized forces.
///
/// Each force row is in the original flattened modal order and contains every
/// declared tone: forces[mode * tone_count + tone], in N/sqrt(kg). The coefficients describe
/// Fc*cos(Omega*t) + Fs*sin(Omega*t) relative to the START of this interval.
/// Frequencies need not be harmonic or FFT bins. Exactly equal frequencies are
/// combined without merging nearby tones; their interference contributes to
/// work and every heat channel.
///
/// Construction retains every original modal, material and acoustic coordinate.
/// The complete source pencil is assembled by the existing exact owner; no
/// independent-component approximation or steady-state substitution is used.
/// Real instantaneous pressure-per-velocity maps remain valid. A nonzero
/// imaginary natural-mode pressure map refuses, because its quadrature is not
/// a broadband observer of the actual forced trajectory.
///
/// The runtime owns the source, preventing a cached propagator from outliving
/// a coefficient/material/radiation change. Access is read-only; consume the
/// runtime with into_network before changing the source and preparing again.
/// The cached physical propagation allocates no matrices and evaluates no
/// exponentials during a step; the original transactional staging and port
/// reporting retain their bounded vector allocations.
pub struct CoupledModalTones {
    network: CoupledModalSystem,
    prepared: PreparedJointTones,
}

impl CoupledModalTones {
    /// Prepare one complete fixed operator at the network's existing clock.
    /// Accepted histories are retained, including a previously started source.
    /// Harmonic-only construction and unresolved tones refuse; no mode is
    /// removed to satisfy the explicit dense dimension/byte/work allowance.
    pub fn new(
        network: CoupledModalSystem,
        angular_frequencies_rad_s: &[f64],
        budget: CoupledToneBudget,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        let prepared = PreparedJointTones::new(&network, angular_frequencies_rad_s, budget, gate)?;
        Ok(Self { network, prepared })
    }

    /// The original accepted mechanical/material/radiation owner.
    #[must_use]
    pub const fn network(&self) -> &CoupledModalSystem {
        &self.network
    }

    /// Recover the accepted original owner without changing its state or clock.
    #[must_use]
    pub fn into_network(self) -> CoupledModalSystem {
        self.network
    }

    /// Original caller tone order, including repeated or zero frequencies.
    #[must_use]
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        self.prepared.angular_frequencies_rad_s()
    }

    #[must_use]
    pub fn tone_count(&self) -> usize {
        self.angular_frequencies_rad_s().len()
    }

    /// Rotate one source's interval-local phase through exactly one sample.
    /// Advance phase only after its physical step has been accepted.
    pub fn advance_force_phase(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), ModalCouplingError> {
        self.prepared.advance_force_phase(tone, cosine, sine)
    }

    /// Advance all original coupled states under the complete continuous load.
    /// A refusal leaves accepted mechanics, histories, frame, clock and ports
    /// unchanged. Scratch may change, as in the original held-force owner.
    pub fn step(
        &mut self,
        cosine_forces: &[f64],
        sine_forces: &[f64],
    ) -> Result<&CoupledModalFrame, ModalCouplingError> {
        self.step_inner(cosine_forces, sine_forces, None)
    }

    /// Cancellable equivalent of step, sharing the original commit boundary.
    pub fn step_under_gate(
        &mut self,
        cosine_forces: &[f64],
        sine_forces: &[f64],
        gate: &CancelGate,
    ) -> Result<&CoupledModalFrame, ModalCouplingError> {
        self.step_inner(cosine_forces, sine_forces, Some(gate))
    }

    fn step_inner(
        &mut self,
        cosine_forces: &[f64],
        sine_forces: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<&CoupledModalFrame, ModalCouplingError> {
        self.prepared
            .stage(&mut self.network, cosine_forces, sine_forces, gate)?;
        self.network.validate_staged_material_points(gate)?;
        poll(gate)?;
        self.network.publish_staged();
        Ok(&self.network.frame)
    }
}

#[cfg(test)]
mod tests;
