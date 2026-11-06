//! Native first-order inviscid evolution with original-mixture entropy admission.

use super::*;
use crate::multicomponent::{ExteriorMixtureEntropyReport, MixtureEntropyControl, MixtureEntropyReport};

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
    /// Advance one entropy-admitted SSP-RK2 step with prescribed exterior gas
    /// states at either physical end.
    ///
    /// Each `Some` replaces that end's reflecting or transmissive boundary with
    /// the original-mixture entropy-stable interface flux. The exterior state is
    /// fixed throughout both stages and every retry, and may have nonzero axial
    /// velocity. A `None` end must be a stationary reflecting wall. Two absent
    /// exterior states execute the closed entropy step exactly.
    ///
    /// These are exterior Riemann traces, not stagnation conditions for the
    /// compact-reservoir characteristic closure in `step_with_reservoirs`.
    /// Borrowed cells must use this gas owner's canonical species order and
    /// species-referenced total energy. Each active species must have common
    /// positive support across every face; common absence is allowed. No trace
    /// seeding, clipping, or change of flux model is performed.
    ///
    /// The boundary supply is `-v_b · transfer + sign * dt * area * psi_b`,
    /// where `v_b` is the entropy variable of the ORIGINAL exterior calorics and
    /// `psi_b = u_b * sum_i(R * rho_i / M_i)`. It includes species chemical
    /// efforts, momentum, kinetic energy and the advective entropy potential.
    /// Transfers use the actual accepted SSP stage fluxes. The direct final
    /// source-state change minus this supply must satisfy one absolute entropy
    /// allowance; each species' change must match its own boundary transfer.
    /// Only entropy from the separately bounded numerical continuity residual
    /// is removed, preserving physical entropy transported by each species.
    ///
    /// First-order space and fixed constant-area geometry are required. Native
    /// CFL, stage admissibility, retries and final cancellation remain in force.
    /// A rejected attempt or cancelled call publishes no field or clock change.
    /// The exterior states are prescribed and are not depleted; finite exterior
    /// storage, moving geometry, source composition and a positivity theorem
    /// require their own coupled treatment.
    pub fn step_entropy_consistent_with_exterior_states(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        entropy: MixtureEntropyControl,
        exterior_states: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(StepReport, ExteriorMixtureEntropyReport), Error> {
        if exterior_states.iter().all(Option::is_none) {
            let (step, closed) =
                self.step_entropy_consistent(maximum_duration, controls, entropy, cancelled)?;
            return Ok((step, ExteriorMixtureEntropyReport {
                source_entropy_change_j_k: closed.source_entropy_change_j_k,
                boundary_entropy_input_j_k: [0.0; 2],
                entropy_generation_j_k: closed.entropy_change_j_k,
                source_entropy_generation_j_k: closed.source_entropy_change_j_k,
                reference_entropy_residual_j_k: closed.reference_entropy_residual_j_k,
                minimum_step_entropy_generation_j_k: closed.minimum_step_entropy_change_j_k,
                maximum_relative_species_residual: closed.maximum_relative_species_residual,
                rejected_entropy_steps: closed.rejected_entropy_steps,
            }));
        }
        poll(&mut cancelled)?;
        if controls.reconstruction != Reconstruction::FirstOrder {
            return Err(Error::EntropyDomain(
                "entropy-consistent exterior flow requires first-order reconstruction",
            ));
        }
        if self.mesh.uniform_cross_section().is_none() {
            return Err(Error::EntropyDomain(
                "entropy-consistent exterior flow requires fixed constant-area geometry",
            ));
        }
        let Boundaries::Ends { left, right } = self.boundaries else {
            return Err(Error::EntropyDomain(
                "prescribed exterior states require physical ends",
            ));
        };
        let mut admitted = [None; 2];
        for (end, boundary) in [left, right].into_iter().enumerate() {
            if !matches!(boundary, Boundary::Reflecting | Boundary::Transmissive) {
                return Err(Error::EntropyDomain(
                    "exterior entropy flow requires reflecting or transmissive end declarations",
                ));
            }
            if let Some(state) = exterior_states[end] {
                poll(&mut cancelled)?;
                admitted[end] = Some((state, self.gas.primitive(state, &mut cancelled)?));
            } else if boundary != Boundary::Reflecting {
                return Err(Error::EntropyDomain(
                    "an open advective end requires its exterior thermodynamic state",
                ));
            }
        }
        let mut account =
            EntropyTracker::new_with_exterior_states(self, entropy, admitted, &mut cancelled)?;
        let report = self.advance_inner::<true>(
            maximum_duration, controls, admitted, None, None,
            Some(&mut account), &mut cancelled,
        )?;
        Ok((report, account.exterior_report()))
    }

}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod exterior_tests;
