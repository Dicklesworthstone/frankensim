//! Complete reacting state entropy with prescribed exterior transport.
//!
//! The immutable pre-composition basis is S=F+sum(alpha_i*m_i/M_i). Actual
//! chemistry changes Q_i and actual advective/molecular exterior transfers B_i are independent
//! physical sources. Admit Delta F + sum(alpha_i*Q_i/M_i) - I_normalized,
//! where I_normalized is reconstructed from the accepted native transport
//! inventory in THIS original basis. The numerical continuity defect is
//! Delta m_i-Q_i-B_i; neither physical source is removed with that defect.

use super::{
    Dd, EntropyTracker, FlowError, MAX_SPECIES, MixtureEuler1d,
    MolecularEntropyControl, ReactiveEntropyAccount, poll, scalar,
};
use super::super::{ExteriorTransportInventory, MixtureGas};
use super::super::advection::ExteriorEntropyAccount;
use super::super::reservoir::EntropyExchange;
use crate::multicomponent::{Boundaries, Cell};

/// Direct original-to-final entropy of one complete open reacting step.
/// Component entropy reports and their allowances are not added into this gate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReactiveExteriorEntropyReport {
    /// Original thermochemical source-state entropy change [J/K], including
    /// thermal, mixing and real chemical conversion contributions.
    pub source_entropy_change_j_k: f64,
    /// Entropy supplied at the left and right fixed exterior states [J/K].
    /// Derived from all actual accepted advective and molecular stage/end transfers,
    /// including momentum/kinetic efforts and the advective entropy potential.
    /// Molecular transfer contributes no advective entropy potential.
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
    /// Internal chemistry, molecular and inviscid retries keep their original counters.
    pub rejected_entropy_steps: usize,
}

/// One optional allocation holds immutable original efforts, actual physical
/// transfers and a receipt only after complete endpoint admission. In a step,
/// boundary_exchange contains molecular halves until the final advective merge.
/// In a window, it contains full exchanges from accepted steps only.
pub(super) struct ReactiveExteriorAccount {
    efforts: ExteriorEntropyAccount,
    scope_boundaries: Boundaries,
    boundary_exchange: EntropyExchange,
    accepted_exchange: Option<EntropyExchange>,
}

impl ReactiveExteriorAccount {
    pub(super) fn begin_attempt(&mut self) {
        self.boundary_exchange = EntropyExchange::default();
        self.accepted_exchange = None;
    }

    pub(super) fn invalidate_receipt(&mut self) {
        self.accepted_exchange = None;
    }

    /// A local copy keeps accumulation transactional on finite-arithmetic or
    /// cancellation refusal. Only physical transfer accounts share this basis;
    /// separately normalized component entropy reports never enter this sum.
    fn merge_exchange(
        &self,
        incoming: &EntropyExchange,
        species_count: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EntropyExchange, FlowError> {
        let mut combined = self.boundary_exchange;
        combined.input = combined.input + incoming.input;
        scalar(combined.input)?;
        for end in 0..2 {
            poll(cancelled)?;
            combined.source_input[end] = combined.source_input[end] + incoming.source_input[end];
            scalar(combined.source_input[end])?;
        }
        for k in 0..species_count {
            poll(cancelled)?;
            combined.species[k] = combined.species[k] + incoming.species[k];
            combined.throughput[k] = combined.throughput[k] + incoming.throughput[k];
            scalar(combined.species[k])?;
            scalar(combined.throughput[k])?;
        }
        poll(cancelled)?;
        Ok(combined)
    }
}

impl ReactiveEntropyAccount {
    /// Retain one original state/source account and cache exterior efforts in its
    /// immutable pre-composition basis, before either molecular or chemical half.
    /// No second flow field or endpoint owner exists.
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
        account.exterior.push(ReactiveExteriorAccount {
            efforts: ExteriorEntropyAccount::new(
                gas, account.state.reference_temperature, &account.state.reference_density,
                admitted, cancelled,
            )?,
            scope_boundaries: initial.boundaries,
            boundary_exchange: EntropyExchange::default(),
            accepted_exchange: None,
        });
        Ok(account)
    }

    /// Create a fresh LOCAL retry anchor without changing the immutable window
    /// basis. Original exterior efforts are copied, not re-evaluated. Local
    /// temperatures, masses and all physical transfer arrays start afresh.
    pub(in crate::multicomponent) fn new_step_in_basis(
        initial: &MixtureEuler1d,
        original: &Self,
        controls: MolecularEntropyControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FlowError> {
        poll(cancelled)?;
        let original_exterior = original.exterior.first().ok_or(FlowError::EntropyDomain(
            "an open entropy step requires its original exterior account",
        ))?;
        if initial.boundaries != original_exterior.scope_boundaries
            || initial.cell_count() != original.state.initial_temperature.len()
            || initial.mesh.uniform_cross_section().is_none()
        {
            return Err(FlowError::EntropyDomain("reacting interval exterior scope changed"));
        }
        let state = EntropyTracker::new_with_basis(
            initial, controls, Some(&original.state), cancelled,
        )?;
        let mut exterior = Vec::new();
        exterior.try_reserve_exact(1).map_err(|_| FlowError::Allocation)?;
        exterior.push(ReactiveExteriorAccount {
            efforts: original_exterior.efforts.clone_efforts(),
            scope_boundaries: original_exterior.scope_boundaries,
            boundary_exchange: EntropyExchange::default(),
            accepted_exchange: None,
        });
        Ok(Self {
            state,
            species_count: original.species_count,
            chemical_change: [Dd::ZERO; MAX_SPECIES],
            chemical_throughput: [Dd::ZERO; MAX_SPECIES],
            rejected_entropy_steps: 0,
            exterior,
        })
    }

    pub(in crate::multicomponent) fn has_exterior_states(&self) -> bool {
        !self.exterior.is_empty()
    }

    /// Merge only the most recently admitted complete step in this identical
    /// source, scope, reference basis and exterior effort chart. Both chemistry
    /// and boundary data finish on local copies before either becomes visible.
    pub(in crate::multicomponent) fn merge_exterior_step(
        &mut self,
        step: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        poll(cancelled)?;
        if self.species_count != step.species_count || !self.state.same_basis(&step.state) {
            return Err(FlowError::EntropyDomain("reacting interval entropy basis mismatch"));
        }
        let original = self.exterior.first().ok_or(FlowError::EntropyDomain(
            "missing original reacting exterior account",
        ))?;
        let accepted = step.exterior.first().ok_or(FlowError::EntropyDomain(
            "missing accepted reacting exterior account",
        ))?;
        if original.scope_boundaries != accepted.scope_boundaries
            || !original.efforts.same_efforts(&accepted.efforts)
        {
            return Err(FlowError::EntropyDomain("reacting interval exterior efforts changed"));
        }
        let exchange = accepted.accepted_exchange.as_ref().ok_or(FlowError::EntropyDomain(
            "reacting interval requires an admitted complete exterior exchange",
        ))?;
        let chemical = self.merged_inventory(&step.accepted_inventory(), cancelled)?;
        let boundary = original.merge_exchange(exchange, self.species_count, cancelled)?;
        poll(cancelled)?;
        self.chemical_change = chemical.chemical_change;
        self.chemical_throughput = chemical.chemical_throughput;
        self.rejected_entropy_steps = chemical.rejected_entropy_steps;
        let original = self.exterior.first_mut().expect("original exterior account retained");
        original.boundary_exchange = boundary;
        original.accepted_exchange = None;
        Ok(())
    }

    /// Recontract one completely admitted molecular half in the ORIGINAL
    /// pre-split entropy basis. Its opaque inventory carries the actual signed
    /// end transfers and accepted absolute stage throughput, including RKL2's
    /// signed weights; it has no advective time-area potential contribution.
    /// Cancellation/refusal leaves all previously observed halves intact.
    pub(in crate::multicomponent) fn observe_molecular_inventory(
        &mut self,
        gas: &MixtureGas,
        inventory: &ExteriorTransportInventory,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FlowError> {
        poll(cancelled)?;
        if gas.species_count() != self.species_count || !self.state.matches_source(gas) {
            return Err(FlowError::EntropyDomain("reacting molecular entropy source species axis"));
        }
        let exterior = self.exterior.first_mut().ok_or(FlowError::EntropyDomain(
            "missing original reacting exterior entropy efforts",
        ))?;
        let incoming = exterior.efforts.exchange_from_inventory(
            gas, inventory, &self.state.reference_entropy, cancelled,
        )?;
        let combined = exterior.merge_exchange(&incoming, self.species_count, cancelled)?;
        exterior.boundary_exchange = combined;
        exterior.accepted_exchange = None;
        Ok(())
    }

    /// Admit both actual chemical/molecular halves and the middle transport
    /// against the ORIGINAL complete-step state. Rounded component entropy
    /// reports and the middle step's reference basis never enter this reduction.
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
        if let Some(exterior) = self.exterior.first_mut() {
            exterior.invalidate_receipt();
        }
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

    /// Admit a complete window directly against its native conserved rollback
    /// rows and the accepted physical transfers already merged into this account.
    /// There is ONE whole-window allowance, independent of local step allowances.
    pub(in crate::multicomponent) fn check_candidate_exterior_rows(
        &mut self,
        initial: &[f64],
        candidate: &MixtureEuler1d,
        temperature_bounds: Option<[f64; 2]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReactiveExteriorEntropyReport, FlowError> {
        self.reset_endpoint();
        if let Some(exterior) = self.exterior.first_mut() {
            exterior.invalidate_receipt();
        }
        let check = (|| {
            self.collect_candidate_rows(initial, candidate, temperature_bounds, cancelled)?;
            let exterior = self.exterior.first().ok_or(FlowError::EntropyDomain(
                "missing original reacting exterior account",
            ))?;
            if candidate.boundaries != exterior.scope_boundaries {
                return Err(FlowError::EntropyDomain("reacting interval exterior scope changed"));
            }
            let exchange = exterior.boundary_exchange;
            let report = self.exterior_report_from_exchange(candidate.gas(), &exchange, cancelled)?;
            poll(cancelled)?;
            self.exterior.first_mut().expect("original exterior account retained")
                .accepted_exchange = Some(exchange);
            Ok(report)
        })();
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
        let advective = exterior.efforts.exchange_from_inventory(
            gas, transport, &self.state.reference_entropy, cancelled,
        )?;
        let exchange = exterior.merge_exchange(&advective, self.species_count, cancelled)?;
        if candidate.boundaries != exterior.scope_boundaries {
            return Err(FlowError::EntropyDomain("reacting exterior scope changed"));
        }
        let report = self.exterior_report_from_exchange(gas, &exchange, cancelled)?;
        poll(cancelled)?;
        self.exterior.first_mut().expect("original exterior account retained")
            .accepted_exchange = Some(exchange);
        Ok(report)
    }

    /// The step and full window share the exact state/source/boundary reduction.
    /// Every actual entropy increment and tolerance is evaluated here once.
    fn exterior_report_from_exchange(
        &self,
        gas: &MixtureGas,
        exchange: &EntropyExchange,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ReactiveExteriorEntropyReport, FlowError> {
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
