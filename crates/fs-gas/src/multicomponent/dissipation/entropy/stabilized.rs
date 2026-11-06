//! Source-conjugate RKL2 evolution with actual macrostep entropy admission.
//!
//! The native Legendre recurrence owns all stages and conservative transfers.
//! Its internal stage states are algebraic auxiliaries, not sequential physical
//! substeps, so entropy monotonicity is checked only at accepted macrostep
//! endpoints and over the complete requested interval. Positivity, caloric,
//! constitutive and parabolic admission remain active at every native stage.

use super::*;
use crate::multicomponent::dissipation::{
    StabilizedDissipationControl, StabilizedDissipationReport,
};

impl<L: TransportLaw> DissipationOperator<L> {
    /// Advance a closed molecular field using the native RKL2 recurrence and
    /// source-conjugate faces, with entropy and species checks before publication.
    ///
    /// This combines extended diffusion timesteps with the same original
    /// caloric, Maxwell--Stefan, Fourier, viscous and Soret/Dufour owners as
    /// `advance_entropy_consistent_for`. The native stage/attempt/RHS budgets
    /// remain authoritative. Every physical macrostep endpoint and the direct
    /// original-to-final state must pass one absolute entropy allowance.
    ///
    /// Final-stage temperatures already admitted by the native RHS are reused.
    /// The entropy check adds no stage constitutive evaluations, no extra EOS
    /// inversion, and no additional per-stage allocation. A refused endpoint
    /// leaves the cached anchor temperatures intact for the original bounded
    /// halving retry. The entire caller field and clock remain unchanged on
    /// cancellation, exhausted work, or final-interval refusal.
    ///
    /// Common species support and stationary zero-heat/species ends are required.
    /// The RKL2 negative-real stability interval is not a nonlinear positivity,
    /// entropy, or temporal-error theorem. Internal Legendre stages retain
    /// native admissibility checks but are not required to increase entropy.
    /// Reservoir, moving, advective and reactive evolution are separate paths.
    pub fn advance_stabilized_entropy_consistent_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        entropy: MolecularEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(StabilizedDissipationReport, MolecularEntropyReport), DissipationError> {
        controls.validate()?;
        let mut account = EntropyTracker::new(flow, entropy, &mut cancelled)?;
        let report = crate::multicomponent::dissipation::stabilized::advance_entropy_consistent(
            self, flow, duration, controls, &mut account, &mut cancelled,
        )?;
        Ok((report, account.report))
    }
}

impl EntropyTracker {
    /// Measure a complete RKL2 candidate against its unchanged physical anchor.
    /// Rejected endpoint temperatures must not become the next retry's anchor.
    pub(in crate::multicomponent::dissipation) fn check_stabilized_candidate(
        &mut self,
        gas: &MixtureGas,
        volumes: &[f64],
        old: &[f64],
        new: &[f64],
        primitive: &[Primitive],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        let n = gas.species_count();
        let m = n + 2;
        self.trial_entropy = Dd::ZERO;
        self.trial_mass.fill(Dd::ZERO);
        self.trial_relative = 0.0;
        let check = (|| {
            for (i, (previous, candidate)) in old
                .chunks_exact(m)
                .zip(new.chunks_exact(m))
                .enumerate()
            {
                poll(cancelled)?;
                let increment = cell_change(
                    gas, &self.reference_density, self.reference_temperature,
                    previous, candidate, self.current_temperature[i],
                    primitive[i].temperature, cancelled,
                )?;
                let volume = Dd::from_f64(volumes[i]);
                self.trial_entropy = self.trial_entropy + volume * increment;
                for k in 0..n {
                    self.trial_mass[k] = self.trial_mass[k]
                        + volume * (Dd::from_f64(candidate[k]) - Dd::from_f64(previous[k]));
                }
            }
            self.trial_relative = self.inventory(n, false)?;
            let change = scalar(self.trial_entropy)?;
            if change < -self.controls.absolute_tolerance_j_k {
                return Err(DissipationError::EntropyDecrease {
                    change_j_k: change,
                    tolerance_j_k: self.controls.absolute_tolerance_j_k,
                });
            }
            Ok(())
        })();
        if matches!(&check, Err(
            DissipationError::EntropyDecrease { .. }
                | DissipationError::SpeciesInventoryDrift { .. }
        )) {
            self.report.rejected_entropy_steps += 1;
        }
        check
    }

    pub(in crate::multicomponent::dissipation) fn accept_stabilized_candidate(
        &mut self,
        n: usize,
        primitive: &[Primitive],
    ) -> Result<(), DissipationError> {
        self.accept_attempt(n)?;
        for (temperature, state) in self.current_temperature.iter_mut().zip(primitive) {
            *temperature = state.temperature;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
