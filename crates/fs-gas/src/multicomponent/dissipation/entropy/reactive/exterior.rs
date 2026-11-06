//! Complete reacting state entropy with prescribed exterior transport.
//!
//! The immutable pre-chemistry basis is S=F+sum(alpha_i*m_i/M_i). Actual
//! chemistry changes Q_i and actual exterior transfers B_i are independent
//! physical sources. Admit Delta F + sum(alpha_i*Q_i/M_i) - I_normalized,
//! where I_normalized is reconstructed from the accepted native transport
//! inventory in THIS original basis. The numerical continuity defect is
//! Delta m_i-Q_i-B_i; neither physical source is removed with that defect.

use super::{
    ExteriorEntropyAccount, FlowError, MixtureEuler1d, MolecularEntropyControl,
    ReactiveEntropyAccount, poll, scalar,
};
use super::super::ExteriorTransportInventory;
use crate::multicomponent::Cell;

/// Direct original-to-final entropy of one complete open reacting step.
/// Component entropy reports and their allowances are not added into this gate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReactiveExteriorEntropyReport {
    /// Original thermochemical source-state entropy change [J/K], including
    /// thermal, mixing and real chemical conversion contributions.
    pub source_entropy_change_j_k: f64,
    /// Entropy supplied at the left and right fixed exterior states [J/K].
    /// Derived from the actual accepted transport stage/end transfers, including
    /// momentum/kinetic efforts and the advective entropy potential.
    pub boundary_entropy_input_j_k: [f64; 2],
    /// Complete state change minus exterior input, with only numerical
    /// species-source continuity reference defects removed [J/K]. This is admitted.
    pub entropy_generation_j_k: f64,
    /// Same generation INCLUDING the separately reported continuity defect [J/K].
    pub source_entropy_generation_j_k: f64,
    /// Sum_i alpha_i*(Delta m_i-Q_i-B_i)/M_i [J/K], in the immutable original basis.
    pub reference_entropy_residual_j_k: f64,
    /// Largest species residual divided by original inventory plus absolute
    /// actual chemical and boundary-stage throughput. The zero scale requires
    /// an exactly zero residual.
    pub maximum_relative_species_residual: f64,
    /// Whole-composition attempts refused by entropy or species-source admission.
    /// Internal chemistry and middle transport retries keep their original counters.
    pub rejected_entropy_steps: usize,
}

impl ReactiveEntropyAccount {
    /// Retain one original state/source account and cache exterior efforts in its
    /// immutable pre-chemistry basis. No second flow field or endpoint owner exists.
    pub(in crate::multicomponent) fn new_with_exterior_states(
        initial: &MixtureEuler1d,
        controls: MolecularEntropyControl,
        exterior: [Option<Cell<'_>>; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FlowError> {
        poll(cancelled)?;
        if exterior.iter().all(Option::is_none) {
            return Err(FlowError::EntropyDomain(
                "open reacting entropy requires a prescribed exterior state",
            ));
        }
        let mut account = Self::new(initial, controls, cancelled)?;
        let gas = initial.gas();
        let mut admitted = [None, None];
        for (end, trace) in exterior.into_iter().enumerate() {
            if let Some(state) = trace {
                let primitive = gas.decode_cell(state, None, cancelled)?;
                admitted[end] = Some((state, primitive));
            }
        }
        account.exterior.try_reserve_exact(1).map_err(|_| FlowError::Allocation)?;
        account.exterior.push(ExteriorEntropyAccount::new(
            gas, account.state.reference_temperature, &account.state.reference_density,
            admitted, cancelled,
        )?);
        Ok(account)
    }

    /// Admit both actual chemical halves and the actual middle transport against
    /// the ORIGINAL complete-step state. Rounded component entropy reports and
    /// the middle step's reference basis never enter this reduction.
    pub(in crate::multicomponent) fn check_candidate_with_exterior(
        &mut self,
        initial: &MixtureEuler1d,
        candidate: &MixtureEuler1d,
        temperature_bounds: Option<[f64; 2]>,
        transport: &ExteriorTransportInventory,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReactiveExteriorEntropyReport, FlowError> {
        if candidate.gas() != initial.gas()
            || initial.mesh != candidate.mesh
            || initial.boundaries != candidate.boundaries
        {
            return Err(FlowError::EntropyDomain("reacting entropy candidate owner or dimensions"));
        }
        self.reset_endpoint();
        let check = self.check_exterior_candidate_inner(
            &initial.states, candidate, temperature_bounds, transport, cancelled,
        );
        if matches!(&check, Err(
            FlowError::EntropyDecrease { .. } | FlowError::SpeciesInventoryDrift { .. }
        )) {
            self.rejected_entropy_steps += 1;
        }
        check
    }

    fn check_exterior_candidate_inner(
        &mut self,
        initial: &[f64],
        candidate: &MixtureEuler1d,
        temperature_bounds: Option<[f64; 2]>,
        transport: &ExteriorTransportInventory,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReactiveExteriorEntropyReport, FlowError> {
        self.collect_candidate_rows(initial, candidate, temperature_bounds, cancelled)?;
        let gas = candidate.gas();
        let exterior = self.exterior.first().ok_or(FlowError::EntropyDomain(
            "missing original reacting exterior entropy efforts",
        ))?;
        let exchange = exterior.exchange_from_inventory(
            gas, transport, &self.state.reference_entropy, cancelled,
        )?;

        let mut physical_change = self.chemical_change;
        let mut throughput = self.chemical_throughput;
        for k in 0..self.species_count {
            poll(cancelled)?;
            physical_change[k] = physical_change[k] + exchange.species[k];
            throughput[k] = throughput[k] + exchange.throughput[k];
            scalar(physical_change[k])?;
            scalar(throughput[k])?;
        }
        let relative = self.state.inventory_with_reaction(
            self.species_count, true, Some((&physical_change, &throughput)),
        )?;
        let (reference, reference_residual, chemical_reference) =
            self.state.reference_changes_with_transfer(
                gas, Some(&self.chemical_change), Some(&exchange.species), cancelled,
            )?;
        // Retain the true source chemical entropy. Boundary input uses the
        // original basis, so only the numerical source-continuity reference
        // defect is absent from this admitted state-minus-input quantity.
        let generation = self.state.trial_entropy + chemical_reference - exchange.input;
        let change = scalar(generation)?;
        if change < -self.state.controls.absolute_tolerance_j_k {
            return Err(FlowError::EntropyDecrease {
                change_j_k: change,
                tolerance_j_k: self.state.controls.absolute_tolerance_j_k,
            });
        }
        Ok(ReactiveExteriorEntropyReport {
            source_entropy_change_j_k: scalar(self.state.trial_entropy + reference)?,
            boundary_entropy_input_j_k: [
                scalar(exchange.source_input[0])?,
                scalar(exchange.source_input[1])?,
            ],
            entropy_generation_j_k: change,
            source_entropy_generation_j_k: scalar(generation + reference_residual)?,
            reference_entropy_residual_j_k: scalar(reference_residual)?,
            maximum_relative_species_residual: relative,
            rejected_entropy_steps: self.rejected_entropy_steps,
        })
    }
}
