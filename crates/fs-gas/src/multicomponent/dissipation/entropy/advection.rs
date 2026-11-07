//! Prescribed exterior entropy on the native inviscid SSP-RK2 transfers.
//!
//! Each exterior effort uses the same fixed entropy basis as the state account.
//! The actual species, momentum and referenced energy transfers, together with
//! the advective entropy potential, supply the open-system entropy input.
//! No molecular heat-only bath formula or characteristic-wave entropy is inferred.

use super::reservoir::EntropyExchange;
use super::{Dd, EntropyTracker, FlowError, MAX_SPECIES, MixtureEuler1d, MixtureGas,
    MolecularEntropyControl, Primitive, poll, positive, scalar};
use crate::multicomponent::evolve::Stage;
use crate::multicomponent::flux::ExteriorEntropyEffort;
use crate::multicomponent::{Cell, MAX_COMPONENTS};

/// Actual state-minus-exterior entropy of one accepted native inviscid step.
/// This finite-step admission does not establish a general invariant domain.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ExteriorMixtureEntropyReport {
    /// Original thermochemical source-state entropy change [J/K].
    pub source_entropy_change_j_k: f64,
    /// Entropy supplied to the gas at the left and right prescribed exteriors
    /// [J/K], using their actual species, momentum and referenced-energy
    /// transfers and advective entropy potentials. Positive means into the gas.
    /// Network receipts sum each end class in retained branch order; internal
    /// finite-plenum ports and periodic faces are excluded from both entries.
    pub boundary_entropy_input_j_k: [f64; 2],
    /// State change minus exterior input, with only the separately bounded
    /// numerical species-continuity reference defect removed [J/K].
    pub entropy_generation_j_k: f64,
    /// Same state-minus-input value INCLUDING the continuity reference defect.
    pub source_entropy_generation_j_k: f64,
    /// Source-reference entropy of (species change - boundary transfer) [J/K].
    pub reference_entropy_residual_j_k: f64,
    /// The admitted complete-step generation [J/K].
    pub minimum_step_entropy_generation_j_k: f64,
    /// Largest individual species residual divided by original inventory plus
    /// absolute throughput at both actual SSP stages and both ends.
    pub maximum_relative_species_residual: f64,
    /// Native attempts refused by entropy or species-continuity admission.
    pub rejected_entropy_steps: usize,
}

/// Actual accepted native end transfers, independent of the entropy basis.
/// Constructed only from native admitted advection stages or complete molecular
/// reports and accepted-stage throughput. Coupled owners contract the transfers
/// in their own immutable original-state basis.
#[derive(Clone, Copy)]
pub(in crate::multicomponent) struct ExteriorTransportInventory {
    species_count: usize,
    transfers: [[f64; MAX_COMPONENTS]; 2],
    signed_measure: [f64; 2],
    throughput: [Dd; MAX_SPECIES],
}

impl ExteriorTransportInventory {
    /// Copy a complete molecular interval after its native state, end-transfer
    /// and bath-entropy gates have succeeded. The report owns the actual signed
    /// conservative transfers; throughput retains every accepted SSP or signed
    /// RKL stage's absolute species contribution. Neither is reconstructed from
    /// a rounded entropy report or a net species change.
    ///
    /// This infallible adapter performs only copies after native publication:
    /// all source axes, finite transfers and throughput were already admitted.
    /// Molecular transport has NO advective entropy-potential contribution.
    pub(super) fn from_molecular_report(
        report: &super::DissipationReport,
        throughput: [Dd; MAX_SPECIES],
    ) -> Self {
        let n = report.before.species_count;
        let mut transfers = [[0.0; MAX_COMPONENTS]; 2];
        for (end, total) in [
            &report.left_boundary_change,
            &report.right_boundary_change,
        ].into_iter().enumerate() {
            transfers[end][..n].copy_from_slice(&total.species_masses[..n]);
            transfers[end][n] = total.momentum;
            transfers[end][n + 1] = total.energy;
        }
        Self {
            species_count: n,
            transfers,
            signed_measure: [0.0; 2],
            throughput,
        }
    }
}

/// Fixed efforts and the last completely admitted native transport inventory.
pub(super) struct ExteriorEntropyAccount {
    ends: [Option<ExteriorEntropyEffort>; 2],
    accepted_inventory: Option<ExteriorTransportInventory>,
}

impl ExteriorEntropyAccount {
    /// Equal source reference coordinates alone do not identify the prescribed
    /// exterior states. Accepted transfers must use the same cached efforts.
    pub(super) fn same_efforts(&self, other: &Self) -> bool {
        self.ends.iter().zip(&other.ends).all(|(left, right)| match (left, right) {
            (None, None) => true,
            (Some(left), Some(right)) => left.same_effort(right),
            _ => false,
        })
    }

    /// A complete interval keeps the same prescribed exterior states and
    /// immutable entropy basis. Copy their cached efforts for a fresh local
    /// step without repeating source intervals/logarithms or retaining an old
    /// accepted transport inventory.
    pub(super) fn clone_efforts(&self) -> Self {
        Self {
            ends: self.ends.clone(),
            accepted_inventory: None,
        }
    }

    pub(super) fn new(
        gas: &MixtureGas,
        reference_temperature: f64,
        reference_density: &[f64; MAX_SPECIES],
        exterior: [Option<(Cell<'_>, Primitive)>; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FlowError> {
        let mut ends = [None, None];
        for (end, trace) in exterior.into_iter().enumerate() {
            if let Some((state, primitive)) = trace {
                ends[end] = Some(ExteriorEntropyEffort::new(
                    gas, state, primitive, reference_temperature, reference_density,
                    cancelled,
                )?);
            }
        }
        Ok(Self { ends, accepted_inventory: None })
    }

    #[allow(clippy::too_many_arguments)]
    fn stage_pair(
        &self,
        gas: &MixtureGas,
        a: &Stage,
        b: &Stage,
        dt: f64,
        areas: [f64; 2],
        accounted_ends: [bool; 2],
        reference_entropy: &[f64; MAX_SPECIES],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(EntropyExchange, ExteriorTransportInventory), FlowError> {
        if !positive(dt) || !areas.into_iter().all(positive) {
            return Err(FlowError::EntropyDomain("exterior entropy quadrature measure"));
        }
        let n = gas.species_count();
        let m = n + 2;
        let mut inventory = ExteriorTransportInventory {
            species_count: n,
            transfers: [[0.0; MAX_COMPONENTS]; 2],
            signed_measure: [0.0; 2],
            throughput: [Dd::ZERO; MAX_SPECIES],
        };
        for end in 0..2 {
            poll(cancelled)?;
            if !accounted_ends[end] {
                if self.ends[end].is_some() {
                    return Err(FlowError::EntropyDomain(
                        "prescribed exterior cannot be an internal network port",
                    ));
                }
                continue;
            }
            let (first, second, sign) = if end == 0 {
                (&a.left, &b.left, 1.0)
            } else {
                (&a.right, &b.right, -1.0)
            };
            for k in 0..m {
                // Match the native accepted per-end report arithmetic exactly.
                let value = sign * 0.5 * dt * areas[end] * (first[k] + second[k]);
                if !value.is_finite() {
                    return Err(FlowError::EntropyDomain("non-finite exterior transfer"));
                }
                inventory.transfers[end][k] = value;
            }
            inventory.signed_measure[end] = sign * dt * areas[end];
            match &self.ends[end] {
                Some(effort) => {
                    // Opposite stage fluxes cannot hide a zero-activity refusal.
                    effort.admit_flux(&first[..m], cancelled)?;
                    effort.admit_flux(&second[..m], cancelled)?;
                }
                None => {
                    // A stationary reflecting end transfers only pressure impulse.
                    if first[..n].iter().chain(&second[..n]).any(|x| *x != 0.0)
                        || first[n + 1] != 0.0 || second[n + 1] != 0.0
                    {
                        return Err(FlowError::EntropyDomain(
                            "unowned advective species or energy transfer",
                        ));
                    }
                }
            }
            let weight = Dd::from_f64(0.5 * dt * areas[end]);
            for k in 0..n {
                poll(cancelled)?;
                inventory.throughput[k] = inventory.throughput[k]
                    + weight * (Dd::from_f64(first[k].abs())
                        + Dd::from_f64(second[k].abs()));
                scalar(inventory.throughput[k])?;
            }
        }
        let exchange = self.exchange_from_inventory(
            gas, &inventory, reference_entropy, cancelled,
        )?;
        Ok((exchange, inventory))
    }

    /// Recontract accepted physical transfers in this account's fixed source
    /// basis. No entropy report from a different split-stage basis is summed.
    pub(super) fn exchange_from_inventory(
        &self,
        gas: &MixtureGas,
        inventory: &ExteriorTransportInventory,
        reference_entropy: &[f64; MAX_SPECIES],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EntropyExchange, FlowError> {
        let n = gas.species_count();
        if inventory.species_count != n {
            return Err(FlowError::EntropyDomain("exterior transfer species axis"));
        }
        let mut exchange = EntropyExchange {
            throughput: inventory.throughput,
            ..EntropyExchange::default()
        };
        for end in 0..2 {
            poll(cancelled)?;
            let transfers = &inventory.transfers[end];
            let input = match &self.ends[end] {
                Some(effort) => effort.entropy_input(
                    &transfers[..n + 2], inventory.signed_measure[end], cancelled,
                )?,
                None => {
                    if transfers[..n].iter().any(|value| *value != 0.0)
                        || transfers[n + 1] != 0.0
                    {
                        return Err(FlowError::EntropyDomain(
                            "unowned accepted exterior species or energy transfer",
                        ));
                    }
                    Dd::ZERO
                }
            };
            let mut reference = Dd::ZERO;
            for (k, model) in gas.source().models().iter().enumerate() {
                poll(cancelled)?;
                let mass = Dd::from_f64(transfers[k]);
                exchange.species[k] = exchange.species[k] + mass;
                reference = reference
                    + mass / Dd::from_f64(model.molar_mass().value())
                        * Dd::from_f64(reference_entropy[k]);
                scalar(exchange.species[k])?;
                scalar(exchange.throughput[k])?;
            }
            exchange.input = exchange.input + input;
            exchange.source_input[end] = input + reference;
            scalar(exchange.source_input[end])?;
        }
        scalar(exchange.input)?;
        Ok(exchange)
    }

}

impl EntropyTracker {
    /// The caller has admitted the exterior states through this flow's original
    /// source, and retains those same fixed traces at both native RK stages.
    pub(in crate::multicomponent) fn new_with_exterior_states(
        flow: &MixtureEuler1d,
        controls: MolecularEntropyControl,
        exterior: [Option<(Cell<'_>, Primitive)>; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FlowError> {
        let mut account = Self::new(flow, controls, cancelled)?;
        if exterior.iter().any(Option::is_some) {
            account.prepare_exterior_pairs(flow.gas(), 1, [exterior], cancelled)?;
        }
        Ok(account)
    }

    /// Cache each admitted physical branch end in the SAME complete-storage
    /// basis. A pair may contain only stationary walls; the native reducer must
    /// still reject unowned species/energy there. Topology supplies the mask
    /// excluding internal finite-volume ports and periodic faces at admission.
    pub(in crate::multicomponent) fn prepare_exterior_pairs<'a>(
        &mut self,
        gas: &MixtureGas,
        count: usize,
        exterior: impl IntoIterator<Item = [Option<(Cell<'a>, Primitive)>; 2]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        poll(cancelled)?;
        if self.boundary.is_some() || self.accepted != 0 || !self.exterior.is_empty()
            || !self.matches_source(gas)
        {
            return Err(FlowError::EntropyDomain("exterior entropy requires a fresh source account"));
        }
        self.exterior.try_reserve_exact(count).map_err(|_| FlowError::Allocation)?;
        for pair in exterior {
            poll(cancelled)?;
            if self.exterior.len() == count {
                return Err(FlowError::EntropyDomain("exterior entropy branch count"));
            }
            self.exterior.push(ExteriorEntropyAccount::new(
                gas, self.reference_temperature, &self.reference_density, pair, cancelled,
            )?);
        }
        if self.exterior.len() != count {
            return Err(FlowError::EntropyDomain("exterior entropy branch count"));
        }
        Ok(())
    }

    /// One complete native candidate uses the original endpoint reducer and
    /// the actual two-stage exterior transfers before any field publication.
    /// Closed transport retains its original report arithmetic and admission.
    pub(in crate::multicomponent) fn admit_transport_step(
        &mut self,
        gas: &MixtureGas,
        a: &Stage,
        b: &Stage,
        dt: f64,
        areas: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        let Some(exterior) = self.exterior.first() else {
            return self.admit_single_step(gas, cancelled);
        };
        if self.boundary.is_some() || self.accepted != 0 {
            return Err(FlowError::EntropyDomain(
                "single-step exterior admission requires a fresh inviscid account",
            ));
        }
        let completion = (|| {
            let (exchange, inventory) = exterior.stage_pair(
                gas, a, b, dt, areas, [true; 2], &self.reference_entropy, cancelled,
            )?;
            Ok((self.exterior_step_report(gas, &exchange, cancelled)?, inventory))
        })();
        match completion {
            Ok((report, inventory)) => {
                self.exterior[0].accepted_inventory = Some(inventory);
                self.completed_exterior_report = Some(report);
                self.accepted = 1;
                Ok(())
            }
            Err(error) => {
                if matches!(
                    &error,
                    FlowError::EntropyDecrease { .. } | FlowError::SpeciesInventoryDrift { .. }
                ) {
                    self.report.rejected_entropy_steps += 1;
                }
                Err(error)
            }
        }
    }

    /// Admit all actual native branch-stage transfers against one network state
    /// increment. Each pair's mask selects only the physical accounting boundary;
    /// finite-plenum ports cancel inside the state inventory, never as supplies.
    /// No branch-local entropy allowance is added to the global allowance.
    pub(in crate::multicomponent) fn admit_network_transport_step<'a>(
        &mut self,
        gas: &MixtureGas,
        dt: f64,
        stages: impl IntoIterator<Item = (&'a Stage, &'a Stage, [f64; 2], [bool; 2])>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        if self.exterior.is_empty() {
            return self.admit_single_step(gas, cancelled);
        }
        if self.boundary.is_some() || self.accepted != 0 {
            return Err(FlowError::EntropyDomain(
                "network exterior admission requires a fresh inviscid account",
            ));
        }
        let completion = (|| {
            let mut total = EntropyExchange::default();
            let mut count = 0;
            for (index, (a, b, areas, mask)) in stages.into_iter().enumerate() {
                poll(cancelled)?;
                let exterior = self.exterior.get_mut(index).ok_or(FlowError::EntropyDomain(
                    "network exterior entropy stage count",
                ))?;
                let (exchange, inventory) = exterior.stage_pair(
                    gas, a, b, dt, areas, mask, &self.reference_entropy, cancelled,
                )?;
                // These physical SSP transfers become readable by a composing
                // owner only after the complete network account is admitted.
                exterior.accepted_inventory = Some(inventory);
                total.input = total.input + exchange.input;
                scalar(total.input)?;
                for end in 0..2 {
                    total.source_input[end] = total.source_input[end] + exchange.source_input[end];
                    scalar(total.source_input[end])?;
                }
                for k in 0..gas.species_count() {
                    poll(cancelled)?;
                    total.species[k] = total.species[k] + exchange.species[k];
                    total.throughput[k] = total.throughput[k] + exchange.throughput[k];
                    scalar(total.species[k])?;
                    scalar(total.throughput[k])?;
                }
                count += 1;
            }
            if count != self.exterior.len() {
                return Err(FlowError::EntropyDomain("network exterior entropy stage count"));
            }
            self.exterior_step_report(gas, &total, cancelled)
        })();
        match completion {
            Ok(report) => {
                self.completed_exterior_report = Some(report);
                self.accepted = 1;
                Ok(())
            }
            Err(error) => {
                for exterior in &mut self.exterior {
                    exterior.accepted_inventory = None;
                }
                if matches!(
                    &error,
                    FlowError::EntropyDecrease { .. } | FlowError::SpeciesInventoryDrift { .. }
                ) {
                    self.report.rejected_entropy_steps += 1;
                }
                Err(error)
            }
        }
    }

    /// Recontract a completely admitted native network step in THIS original
    /// source basis. Prescribed efforts differ between physical ends, so their
    /// conservative inventories cannot be combined before this contraction.
    /// Internal ports were masked by the native stage owner, not inferred here.
    pub(super) fn exchange_from_network_transport(
        &self,
        gas: &MixtureGas,
        transport: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EntropyExchange, FlowError> {
        poll(cancelled)?;
        if !self.matches_source(gas) || !transport.matches_source(gas)
            || self.exterior.is_empty()
            || self.exterior.len() != transport.exterior.len()
            || self.initial_temperature.len() != transport.initial_temperature.len()
            || transport.accepted != 1 || transport.completed_exterior_report.is_none()
        {
            return Err(FlowError::EntropyDomain(
                "reacting network requires an admitted matching exterior transport account",
            ));
        }
        let mut total = EntropyExchange::default();
        for (original, accepted) in self.exterior.iter().zip(&transport.exterior) {
            poll(cancelled)?;
            let inventory = accepted.accepted_inventory.as_ref().ok_or(
                FlowError::EntropyDomain("missing accepted network exterior transfers"),
            )?;
            let exchange = original.exchange_from_inventory(
                gas, inventory, &self.reference_entropy, cancelled,
            )?;
            total.input = total.input + exchange.input;
            scalar(total.input)?;
            for end in 0..2 {
                total.source_input[end] = total.source_input[end] + exchange.source_input[end];
                scalar(total.source_input[end])?;
            }
            for k in 0..gas.species_count() {
                poll(cancelled)?;
                total.species[k] = total.species[k] + exchange.species[k];
                total.throughput[k] = total.throughput[k] + exchange.throughput[k];
                scalar(total.species[k])?;
                scalar(total.throughput[k])?;
            }
        }
        Ok(total)
    }

    fn exterior_step_report(
        &self,
        gas: &MixtureGas,
        exchange: &EntropyExchange,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExteriorMixtureEntropyReport, FlowError> {
        // Subtract actual boundary transfer before applying source reference
        // constants; physical entropy carried by entering/leaving species stays.
        let relative = self.inventory_with_reaction(
            gas.species_count(), true,
            Some((&exchange.species, &exchange.throughput)),
        )?;
        let generation = self.trial_entropy - exchange.input;
        let change = scalar(generation)?;
        if change < -self.controls.absolute_tolerance_j_k {
            return Err(FlowError::EntropyDecrease {
                change_j_k: change,
                tolerance_j_k: self.controls.absolute_tolerance_j_k,
            });
        }
        let (reference, reference_residual, _) =
            self.reference_changes(gas, Some(&exchange.species), cancelled)?;
        Ok(ExteriorMixtureEntropyReport {
            source_entropy_change_j_k: scalar(self.trial_entropy + reference)?,
            boundary_entropy_input_j_k: [
                scalar(exchange.source_input[0])?, scalar(exchange.source_input[1])?,
            ],
            entropy_generation_j_k: change,
            source_entropy_generation_j_k: scalar(generation + reference_residual)?,
            reference_entropy_residual_j_k: scalar(reference_residual)?,
            minimum_step_entropy_generation_j_k: change,
            maximum_relative_species_residual: relative,
            rejected_entropy_steps: self.report.rejected_entropy_steps,
        })
    }

    /// Only a fully admitted native SSP candidate can supply a coupled owner.
    pub(in crate::multicomponent) fn exterior_transport_inventory(
        &self,
    ) -> ExteriorTransportInventory {
        self.exterior.first()
            .and_then(|account| account.accepted_inventory)
            .expect("accepted exterior step retained native conservative transfers")
    }

    /// All source, balance and report arithmetic already completed before the
    /// native accepted-field swap. The closed case preserves the original gate.
    pub(in crate::multicomponent) fn exterior_report(&self) -> ExteriorMixtureEntropyReport {
        if self.exterior.is_empty() {
            let report = self.report;
            return ExteriorMixtureEntropyReport {
                source_entropy_change_j_k: report.source_entropy_change_j_k,
                boundary_entropy_input_j_k: [0.0; 2],
                entropy_generation_j_k: report.entropy_change_j_k,
                source_entropy_generation_j_k: report.source_entropy_change_j_k,
                reference_entropy_residual_j_k: report.reference_entropy_residual_j_k,
                minimum_step_entropy_generation_j_k: report.minimum_step_entropy_change_j_k,
                maximum_relative_species_residual: report.maximum_relative_species_residual,
                rejected_entropy_steps: report.rejected_entropy_steps,
            };
        }
        self.completed_exterior_report
            .expect("accepted inviscid step retained complete exterior entropy admission")
    }
}
