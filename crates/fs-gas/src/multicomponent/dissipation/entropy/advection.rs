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
/// Only this module constructs the inventory after checking BOTH SSP stages;
/// coupled owners can contract it in their own immutable original-state basis.
#[derive(Clone, Copy)]
pub(in crate::multicomponent) struct ExteriorTransportInventory {
    species_count: usize,
    transfers: [[f64; MAX_COMPONENTS]; 2],
    signed_measure: [f64; 2],
    throughput: [Dd; MAX_SPECIES],
}

/// Fixed efforts and the last completely admitted native transport inventory.
pub(super) struct ExteriorEntropyAccount {
    ends: [Option<ExteriorEntropyEffort>; 2],
    accepted_inventory: Option<ExteriorTransportInventory>,
}

impl ExteriorEntropyAccount {
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
            account.exterior.try_reserve_exact(1).map_err(|_| FlowError::Allocation)?;
            account.exterior.push(ExteriorEntropyAccount::new(
                flow.gas(), account.reference_temperature, &account.reference_density,
                exterior, cancelled,
            )?);
        }
        Ok(account)
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
                gas, a, b, dt, areas, &self.reference_entropy, cancelled,
            )?;
            // These shared reducers accept an additional physical species
            // source independently of the state-entropy basis. Here that source
            // is the actual exterior transfer, not a chemical reaction.
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
            Ok((ExteriorMixtureEntropyReport {
                source_entropy_change_j_k: scalar(self.trial_entropy + reference)?,
                boundary_entropy_input_j_k: [
                    scalar(exchange.source_input[0])?,
                    scalar(exchange.source_input[1])?,
                ],
                entropy_generation_j_k: change,
                source_entropy_generation_j_k: scalar(generation + reference_residual)?,
                reference_entropy_residual_j_k: scalar(reference_residual)?,
                minimum_step_entropy_generation_j_k: change,
                maximum_relative_species_residual: relative,
                rejected_entropy_steps: self.report.rejected_entropy_steps,
            }, inventory))
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
