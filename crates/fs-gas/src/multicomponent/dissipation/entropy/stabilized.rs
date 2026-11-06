//! Source-conjugate RKL2 evolution with actual macrostep entropy admission.
//!
//! The native Legendre recurrence owns all stages and conservative transfers.
//! Its internal stage states are algebraic auxiliaries, not sequential physical
//! substeps, so entropy monotonicity is checked only at accepted macrostep
//! endpoints and over the complete requested interval. Positivity, caloric,
//! constitutive and parabolic admission remain active at every native stage.

use super::*;
use crate::multicomponent::dissipation::{
    Rkl2Scheme, StabilizedDissipationControl, StabilizedDissipationReport,
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<(StabilizedDissipationReport, MolecularEntropyReport), DissipationError> {
        self.advance_stabilized_entropy_consistent_budget(
            flow, duration, controls, entropy, &mut 0, cancelled,
        )
    }

    /// Original closed RKL2/entropy owner with a caller-retained RHS counter.
    /// Coupled retries must not obtain a fresh molecular work budget.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent) fn advance_stabilized_entropy_consistent_budget(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        entropy: MolecularEntropyControl,
        evaluations: &mut usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(StabilizedDissipationReport, MolecularEntropyReport), DissipationError> {
        controls.validate()?;
        let mut account = EntropyTracker::new(flow, entropy, &mut cancelled)?;
        let report = crate::multicomponent::dissipation::stabilized::advance_entropy_consistent_budget(
            self, flow, duration, controls, [None; 2], evaluations,
            &mut account, &mut cancelled,
        )?;
        Ok((report, account.report))
    }

    /// Advance source-conjugate RKL2 transport with stationary thermal/species baths.
    ///
    /// Bath entropy and species continuity use the SAME signed native RKL2
    /// boundary transfers as the conservative field. Relative inventory
    /// tolerances use absolute throughput from the actual tableau weights,
    /// including the negative anchor weight and every contributing stage/end.
    /// Cancellation between stages or ends cannot hide this throughput.
    ///
    /// Each actual endpoint and the complete requested interval must pass
    /// state-minus-bath entropy admission before the field publishes. Source
    /// gauges, stationary bath chemical efforts, common species support and
    /// original caloric/constitutive domains remain those of the SSP entry point
    /// `advance_entropy_consistent_with_reservoirs_for`. Time-varying or moving
    /// baths and advective/reactive exchange are not included.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_stabilized_entropy_consistent_with_reservoirs_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        entropy: MolecularEntropyControl,
        reservoirs: [Option<super::super::super::Cell<'_>>; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(StabilizedDissipationReport, ReservoirMolecularEntropyReport), DissipationError> {
        let (report, measured, _) = self.advance_stabilized_reservoir_entropy(
            flow, duration, controls, entropy, reservoirs, None, cancelled,
        )?;
        Ok((report, measured))
    }

    /// Retain the original RHS budget across molecular halves and abandoned
    /// coupled attempts, while exporting actual accepted end transfers and
    /// absolute signed-tableau throughput for the complete source account.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent) fn advance_stabilized_entropy_consistent_with_reservoir_inventory_budget(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        entropy: MolecularEntropyControl,
        reservoirs: [Option<super::super::super::Cell<'_>>; 2],
        evaluations: &mut usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(
        StabilizedDissipationReport,
        ReservoirMolecularEntropyReport,
        ExteriorTransportInventory,
    ), DissipationError> {
        self.advance_stabilized_reservoir_entropy(
            flow, duration, controls, entropy, reservoirs, Some(evaluations), cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_stabilized_reservoir_entropy(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        entropy: MolecularEntropyControl,
        reservoirs: [Option<super::super::super::Cell<'_>>; 2],
        evaluations: Option<&mut usize>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(
        StabilizedDissipationReport,
        ReservoirMolecularEntropyReport,
        ExteriorTransportInventory,
    ), DissipationError> {
        controls.validate()?;
        entropy.validate()?;
        let reservoirs = if reservoirs.iter().all(Option::is_none) {
            [None; 2]
        } else {
            crate::multicomponent::dissipation::reservoir::admit(
                &self.gas, self.boundaries, reservoirs, &mut cancelled,
            )?
        };
        let mut account = EntropyTracker::new(flow, entropy, &mut cancelled)?;
        account.boundary = Some(ReservoirEntropyAccount::new(
            &self.gas, account.reference_temperature, &account.reference_density,
            &account.reference_entropy, reservoirs, &mut cancelled,
        )?);
        let report = if let Some(evaluations) = evaluations {
            crate::multicomponent::dissipation::stabilized::advance_entropy_consistent_budget(
                self, flow, duration, controls, reservoirs, evaluations,
                &mut account, &mut cancelled,
            )?
        } else {
            // Preserve the standalone public owner's original local counter.
            crate::multicomponent::dissipation::stabilized::advance_entropy_consistent(
                self, flow, duration, controls, reservoirs, &mut account, &mut cancelled,
            )?
        };
        // Native publication and all fallible gates have completed. Only copy
        // the actual signed end transfers and admitted absolute stage throughput.
        let boundary = account.boundary.as_ref()
            .expect("native RKL2 reservoir interval retained its admitted bath account");
        let inventory = ExteriorTransportInventory::from_molecular_report(
            &report.transport, boundary.total.throughput,
        );
        Ok((report, account.completed_reservoir_report
            .expect("native RKL2 final admission retained reservoir entropy"), inventory))
    }
}

impl EntropyTracker {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent::dissipation) fn begin_stabilized_quadrature(
        &mut self,
        scheme: Rkl2Scheme,
        weights: &mut [f64; 256],
        anchor: &crate::multicomponent::dissipation::evolve::Stage,
        dt: f64,
        areas: [f64; 2],
        n: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        if let Some(boundary) = &mut self.boundary {
            scheme.rhs_weights(&mut weights[..scheme.stage_count()])?;
            boundary.begin_quadrature();
            boundary.observe_quadrature_stage(anchor, dt * weights[0], areas, n, cancelled)?;
        }
        Ok(())
    }

    pub(in crate::multicomponent::dissipation) fn observe_stabilized_quadrature(
        &mut self,
        stage: &crate::multicomponent::dissipation::evolve::Stage,
        weighted_dt: f64,
        areas: [f64; 2],
        n: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        if let Some(boundary) = &mut self.boundary {
            boundary.observe_quadrature_stage(stage, weighted_dt, areas, n, cancelled)?;
        }
        Ok(())
    }

    /// Measure a complete RKL2 candidate against its unchanged physical anchor.
    /// Rejected endpoint temperatures must not become the next retry's anchor.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent::dissipation) fn check_stabilized_candidate(
        &mut self,
        gas: &MixtureGas,
        volumes: &[f64],
        old: &[f64],
        new: &[f64],
        primitive: &[Primitive],
        transfers: &[[f64; crate::multicomponent::dissipation::MAX_COMPONENTS]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        let n = gas.species_count();
        let m = n + 2;
        self.trial_entropy = Dd::ZERO;
        self.trial_mass.fill(Dd::ZERO);
        self.trial_relative = 0.0;
        // Native RKL2 reports retain the molecular error surface; the common
        // state/bath account now returns the original gas owner's errors.
        let check: Result<(), DissipationError> = (|| {
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
            if let Some(boundary) = &mut self.boundary {
                boundary.finish_quadrature(transfers, n, cancelled)?;
            }
            self.trial_relative = self.inventory(n, false)?;
            let input = self.boundary.as_ref().map_or(Dd::ZERO, |b| b.trial.input);
            let change = scalar(self.trial_entropy - input)?;
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
#[cfg(test)]
mod reservoir_tests;
