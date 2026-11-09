//! Actual finite-rate mineral phases and their original speciated aqueous view.
//!
//! The native mineral owner retains the selected source reaction, complete
//! analytical stocks, independently represented dissolved/solid amounts and
//! explicit surface-rate policy. It solves the same acid/base and charge
//! equations at each actual dissolved inventory. This adapter neither computes
//! another reaction rate nor projects the crystal stock to equilibrium.
//!
//! Solvent-volume changes hold the actual phases fixed and re-speciate only
//! the aqueous chemistry on the original explicit activity basis. Published
//! osmotic particles are H+, OH-, every actual family form and spectator once;
//! neither solid formula units nor an inferred dissociation multiplier enter.
//! Temperature and pressure remain those of the original ideal dilute source.
//! Heat, reaction water, species transport and additional complexes or solids
//! are not supplied by this material handoff.

use super::{
    AcidBaseSpeciation, AqueousSolution, Result, SolventActivityBasis, SpeciesId, SpectatorIon,
    admit_layout, poll, publish_aqueous,
};
use fs_qty::{Amount, Time, Volume};
use fs_thermochem::stoichiometry::thermodynamics::acid_base::{
    AcidBaseSpecies, FamilySpeciation, SpeciationOptions,
    precipitation::kinetics::{
        MineralKineticControlsV1, MineralKineticOutcomeV1, MineralKineticStateV1, MineralKineticsV1,
    },
};
use std::sync::Arc;

/// One original finite mineral state and its dissolved-only material view.
/// The native state owns all analytical stocks and source bindings; this type
/// stores no second component balance or separate crystal inventory.
#[derive(Clone, Debug)]
pub struct KineticElectrolyteSolubilitySolution {
    basis: SolventActivityBasis,
    state: MineralKineticStateV1,
    species: Arc<[SpeciesId]>,
    solution: AqueousSolution,
}

/// Accepted material endpoint and the original kinetic integration report.
#[derive(Clone, Debug)]
pub struct KineticElectrolyteSolubilityStep {
    solution: KineticElectrolyteSolubilitySolution,
    kinetics: MineralKineticOutcomeV1,
}
impl KineticElectrolyteSolubilityStep {
    pub const fn solution(&self) -> &KineticElectrolyteSolubilitySolution {
        &self.solution
    }
    pub const fn kinetics(&self) -> &MineralKineticOutcomeV1 {
        &self.kinetics
    }
    pub fn into_parts(self) -> (KineticElectrolyteSolubilitySolution, MineralKineticOutcomeV1) {
        (self.solution, self.kinetics)
    }
}

impl KineticElectrolyteSolubilitySolution {
    /// Bind explicitly supplied actual phases to the original rate/source owner.
    ///
    /// `dissolved_cation` and every `dissolved_family_totals` entry describe
    /// aqueous component stocks ONLY. `solid` counts actual formula units;
    /// the original reaction stoichiometry adds their contribution to complete
    /// analytical totals. A supersaturated or undersaturated starting state
    /// is preserved. Acid/base partition and charge balance are solved at
    /// these actual amounts, without dissolving or precipitating the solid.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        model: MineralKineticsV1,
        basis: SolventActivityBasis,
        solvent_volume: Volume,
        dissolved_cation: Amount,
        solid: Amount,
        dissolved_family_totals: Vec<(String, Amount)>,
        spectators: Vec<SpectatorIon>,
        options: SpeciationOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        admit_layout(model.source().aqueous_model(), spectators.len())?;
        let amount = basis.amount_at_unit_activity(solvent_volume)?;
        let state = model.prepare(
            dissolved_cation,
            solid,
            dissolved_family_totals,
            spectators,
            amount,
            options,
            &mut cancelled,
        )?;
        Self::publish(basis, solvent_volume, state, None, &mut cancelled)
    }

    /// Adopt an already prepared native phase state without another chemistry
    /// solve. Its activity-reference amount must exactly match the supplied
    /// physical volume and explicit reference density; neither is inferred.
    pub fn from_state(
        basis: SolventActivityBasis,
        solvent_volume: Volume,
        state: MineralKineticStateV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::publish(basis, solvent_volume, state, None, &mut cancelled)
    }

    pub const fn basis(&self) -> SolventActivityBasis {
        self.basis
    }
    pub const fn state(&self) -> &MineralKineticStateV1 {
        &self.state
    }
    pub fn model(&self) -> &MineralKineticsV1 {
        self.state.model()
    }
    /// Actual pH, original acid/base source, dissolved family forms and charge
    /// diagnostics. Complete analytical totals are exposed separately below.
    pub fn aqueous(&self) -> &AcidBaseSpeciation {
        self.state.aqueous()
    }
    pub const fn solution(&self) -> &AqueousSolution {
        &self.solution
    }
    /// Retained numerical buffers that are copied with this material state,
    /// excluding the inline value. This includes the actual aqueous report and
    /// published particle inventory, with their retained identifier buffers.
    ///
    /// Immutable source models, analytical stocks and the species axis are
    /// shared through `Arc` and excluded: they already exist before numerical
    /// workspace preparation. Transient allocations inside a native chemistry
    /// solve are also excluded. This is a per-state numerical storage estimate,
    /// not a bound on total or peak solver memory.
    pub fn workspace_heap_bytes(&self) -> usize {
        let aqueous = self.aqueous();
        let mut bytes = self.solution.heap_bytes()
            + aqueous.families.capacity() * core::mem::size_of::<FamilySpeciation>()
            + aqueous.spectators.capacity() * core::mem::size_of::<SpectatorIon>();
        for family in &aqueous.families {
            bytes += family.id.capacity()
                + family.forms.capacity() * core::mem::size_of::<AcidBaseSpecies>();
            for form in &family.forms {
                bytes += form.species.heap_bytes();
            }
        }
        for spectator in &aqueous.spectators {
            bytes += spectator.species.heap_bytes();
        }
        bytes
    }
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }
    pub fn solvent(&self) -> &SpeciesId {
        self.aqueous().model().water().water()
    }
    pub fn solid_amount(&self) -> Amount {
        self.state.solid_amount()
    }
    pub fn dissolved_cation(&self) -> Amount {
        self.state.dissolved_cation()
    }
    /// Actual dissolved total of the salt's selected family, when present.
    /// Other dissolved families remain available through `aqueous().families`.
    pub fn dissolved_family(&self) -> Option<Amount> {
        self.state.dissolved_family()
    }
    /// Complete cation stock, including the solid's original coefficient.
    pub fn component_total(&self) -> Amount {
        self.state.component_total()
    }
    /// Complete analytical family totals, including their solid contribution.
    pub fn family_totals(&self) -> &[(String, Amount)] {
        self.state.family_totals()
    }
    /// Original fixed spectator declarations, excluding the solved cation.
    pub fn spectators(&self) -> &[SpectatorIon] {
        self.state.spectators()
    }

    /// Pure-solvent volume change at the source temperature and pressure.
    /// Both actual component phases remain fixed. Only aqueous acid/base
    /// partition and charge balance are re-solved at the new activity basis.
    pub fn with_solvent_volume_m3(&self, solvent_volume_m3: f64) -> Result<Self> {
        self.with_solvent_volume_m3_cancelled(solvent_volume_m3, || false)
    }
    pub fn with_solvent_volume_m3_cancelled(
        &self,
        solvent_volume_m3: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let volume = Volume::new(solvent_volume_m3);
        let amount = self.basis.amount_at_unit_activity(volume)?;
        let state = self.model().respeciate_at_basis(&self.state, amount, &mut cancelled)?;
        Self::publish(
            self.basis,
            volume,
            state,
            Some(Arc::clone(&self.species)),
            &mut cancelled,
        )
    }

    /// Apply explicitly supplied aqueous dosing/removal while retaining the
    /// actual crystal stock and every original chemical identity and source.
    ///
    /// Amounts follow `family_totals()` and `spectators()` ordering. Family
    /// amounts here are DISSOLVED pools, unlike the complete analytical getter;
    /// the same native reaction adds the unchanged solid's stoichiometric
    /// contribution to form new complete totals. Spectator charges remain
    /// unchanged. Counterions must be supplied through already declared ions;
    /// none are invented. Only the aqueous acid/base partition is solved now;
    /// later `advance` calls evolve any resulting supersaturation at finite rate.
    pub fn with_dissolved_analytical_totals(
        &self,
        dissolved_cation: Amount,
        dissolved_family_amounts: &[Amount],
        spectator_amounts: &[Amount],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if dissolved_family_amounts.len() != self.family_totals().len()
            || spectator_amounts.len() != self.spectators().len()
        {
            return Err(super::ElectrolyteSolubilityError::Invalid(
                "aqueous dosing must retain original family and spectator ordering",
            ));
        }
        let families = self.family_totals().iter()
            .zip(dissolved_family_amounts)
            .map(|((id, _), &amount)| (id.clone(), amount))
            .collect();
        let spectators = self.spectators().iter()
            .zip(spectator_amounts)
            .map(|(ion, &amount)| SpectatorIon {
                species: ion.species.clone(), charge: ion.charge, amount,
            })
            .collect();
        let state = self.model().prepare(
            dissolved_cation,
            self.solid_amount(),
            families,
            spectators,
            self.state.amount_at_unit_activity(),
            self.state.options(),
            &mut cancelled,
        )?;
        Self::publish(
            self.basis,
            Volume::new(self.solution.solvent_volume_m3()),
            state,
            Some(Arc::clone(&self.species)),
            &mut cancelled,
        )
    }

    /// Evolve the actual mineral at fixed physical solvent volume. The native
    /// activity-driven rate and adaptive controller own the complete step.
    /// A refusal, including late publication cancellation, changes no input.
    pub fn advance(
        &self,
        duration: Time,
        controls: MineralKineticControlsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<KineticElectrolyteSolubilityStep> {
        poll(&mut cancelled)?;
        let kinetics = self.model().advance(&self.state, duration, controls, &mut cancelled)?;
        let solution = Self::publish(
            self.basis,
            Volume::new(self.solution.solvent_volume_m3()),
            kinetics.final_state().clone(),
            Some(Arc::clone(&self.species)),
            &mut cancelled,
        )?;
        Ok(KineticElectrolyteSolubilityStep { solution, kinetics })
    }

    fn publish(
        basis: SolventActivityBasis,
        volume: Volume,
        state: MineralKineticStateV1,
        species_axis: Option<Arc<[SpeciesId]>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if state.amount_at_unit_activity() != basis.amount_at_unit_activity(volume)? {
            return Err(super::ElectrolyteSolubilityError::Invalid(
                "native mineral activity basis must match the original physical solvent volume",
            ));
        }
        let source = state.model().source();
        admit_layout(source.aqueous_model(), state.spectators().len())?;
        let (species, solution) = publish_aqueous(
            basis,
            volume,
            state.aqueous(),
            source.solid(),
            source.standard().temperature(),
            species_axis,
            cancelled,
        )?;
        Ok(Self { basis, state, species, solution })
    }
}

#[cfg(test)]
mod tests;
