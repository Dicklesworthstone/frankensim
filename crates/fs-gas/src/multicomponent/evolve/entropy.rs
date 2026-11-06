//! Native first-order inviscid evolution with original-mixture entropy admission.

use super::*;
use crate::multicomponent::{MixtureEntropyControl, MixtureEntropyReport};

impl MixtureEuler1d {
    /// Advance one native SSP-RK2 step with an original-source entropy-stable
    /// interior flux and measured entropy admission of the actual final state.
    ///
    /// Species, momentum and referenced total energy share a flux derived from
    /// the original ideal-mixture calorics, including temperature-dependent
    /// heat capacity. The native CFL, intermediate-state checks, bounded retry
    /// loop, conservative balances and clock publication remain authoritative.
    /// Final-state temperatures already inverted by the native candidate loop
    /// are reused by the same entropy account as molecular transport.
    ///
    /// This entry point requires first-order space, fixed constant-area cells,
    /// and either periodic boundaries or two stationary reflecting walls.
    /// The original reflecting-wave closure supplies wall momentum; stationary
    /// walls transport no species or energy and carry no external entropy.
    /// Every active species must be positive on both sides of each interior
    /// face, or absent on both. No trace seeding, state clipping, or switch to
    /// the engineering flux is performed.
    ///
    /// The returned duration may be shorter than the maximum because of the
    /// native CFL or entropy/admissibility retries. The entropy allowance applies
    /// to this complete accepted step, including each species continuity check;
    /// it is not a cumulative allowance across separate calls. Refusal and
    /// cancellation leave the entire accepted field and clock unchanged.
    ///
    /// Semidiscrete entropy dissipation does not imply a positive finite-step
    /// invariant domain or exact moving-contact pressure equilibrium. The actual
    /// finite-step checks remain necessary. MUSCL, moving/variable-area geometry,
    /// open advective boundaries, forcing, chemistry and molecular composition
    /// need their own complete entropy treatment and are not enabled here.
    pub fn step_entropy_consistent(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(StepReport, MixtureEntropyReport), Error> {
        poll(&mut cancelled)?;
        if controls.reconstruction != Reconstruction::FirstOrder {
            return Err(Error::EntropyDomain(
                "entropy-consistent flow requires first-order reconstruction",
            ));
        }
        if self.mesh.uniform_cross_section().is_none() {
            return Err(Error::EntropyDomain(
                "entropy-consistent flow requires fixed constant-area geometry",
            ));
        }
        if !matches!(self.boundaries,
            Boundaries::Periodic
            | Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            }
        ) {
            return Err(Error::EntropyDomain(
                "entropy-consistent flow requires periodic or stationary reflecting ends",
            ));
        }
        let mut account = EntropyTracker::new(self, entropy, &mut cancelled)?;
        let report = self.advance_inner::<true>(
            maximum_duration, controls, [None; 2], None, None,
            Some(&mut account), &mut cancelled,
        )?;
        Ok((report, account.report()))
    }
}

#[cfg(test)]
mod tests;
