//! Finite copper-sulfide solubility coupled to the original metal/ligand owner.
//!
//! Every volume or analytical-stock update invokes the original simultaneous
//! Cu/S partition, metal association, ligand and charge solve. Complete input
//! stocks remain distinct from the actual dissolved report and selected solid.
//! No apparent neutral solubility or fixed-pH approximation replaces that solve.
//!
//! The material view counts each actual dissolved particle once. Charged metal
//! members already occur in the aqueous report; neutral metal members are added
//! from the original ligand result. Solid formula units and bound ligand units
//! are not additional particles. Temperature, pressure, the explicit activity
//! basis and the selected pH domain remain those of the original ideal model.
//! The dilute osmotic view retains the water-ion background described by the
//! parent module. It supplies no reaction heat, kinetics, ionic transport,
//! omitted complexes, competing solids or concentrated-electrolyte correction.

use super::{
    AcidBaseSpeciation, AqueousSolution, ElectrolyteSolubilityError, MAX_OSMOTIC_SOLUTES, Result,
    SolventActivityBasis, SpeciesId, SpectatorIon, poll, publish_aqueous_with_additional_particles,
};
use fs_qty::{Amount, Volume};
use fs_thermochem::stoichiometry::thermodynamics::acid_base::complexation::{
    LigandSpeciation,
    precipitation::{CopperSulfideOptions, CopperSulfidePrecipitation, CopperSulfideSpeciation},
};
use std::sync::Arc;

#[derive(Clone, Debug)]
struct AnalyticalStocks {
    ligand_total: Amount,
    metal_totals: Box<[(String, Amount)]>,
    other_acid_totals: Box<[(String, Amount)]>,
    spectators: Box<[SpectatorIon]>,
}
impl AnalyticalStocks {
    fn new(
        ligand_total: Amount,
        mut metal_totals: Vec<(String, Amount)>,
        mut other_acid_totals: Vec<(String, Amount)>,
        mut spectators: Vec<SpectatorIon>,
    ) -> Self {
        metal_totals.sort_by(|a, b| a.0.cmp(&b.0));
        other_acid_totals.sort_by(|a, b| a.0.cmp(&b.0));
        spectators.sort_by(|a, b| a.species.cmp(&b.species));
        Self {
            ligand_total,
            metal_totals: metal_totals.into_boxed_slice(),
            other_acid_totals: other_acid_totals.into_boxed_slice(),
            spectators: spectators.into_boxed_slice(),
        }
    }
}

/// Original finite Cu/S/CN equilibrium and its complete dissolved material view.
///
/// Cloned or volume-updated values retain the complete analytical stocks and
/// dissolved species axis. All phase and aqueous chemistry belongs to the
/// retained native owner; no separate material phase calculation is performed.
#[derive(Clone, Debug)]
pub struct ComplexedElectrolyteSolubilitySolution {
    basis: SolventActivityBasis,
    options: CopperSulfideOptions,
    stocks: Arc<AnalyticalStocks>,
    equilibrium: CopperSulfideSpeciation,
    species: Arc<[SpeciesId]>,
    solution: AqueousSolution,
}
impl ComplexedElectrolyteSolubilitySolution {
    /// Bind one selected Cu2S phase and its complete original ligand chemistry.
    ///
    /// `ligand_total` includes free, protonated and metal-bound ligand. Copper
    /// and sulfide analytical totals include respectively two and one component
    /// units per solid formula unit. Supply every original metal and other acid
    /// family exactly once, including zero totals. Counterions remain explicit.
    #[allow(clippy::too_many_arguments)]
    pub fn from_copper_sulfide(
        model: CopperSulfidePrecipitation,
        basis: SolventActivityBasis,
        solvent_volume: Volume,
        ligand_total: Amount,
        metal_totals: Vec<(String, Amount)>,
        other_acid_totals: Vec<(String, Amount)>,
        spectators: Vec<SpectatorIon>,
        options: CopperSulfideOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        admit_layout(
            &model,
            metal_totals.len(),
            other_acid_totals.len(),
            spectators.len(),
        )?;
        let amount = basis.amount_at_unit_activity(solvent_volume)?;
        let stocks = Arc::new(AnalyticalStocks::new(
            ligand_total,
            metal_totals,
            other_acid_totals,
            spectators,
        ));
        let equilibrium = model.speciate(
            stocks.ligand_total,
            &stocks.metal_totals,
            &stocks.other_acid_totals,
            &stocks.spectators,
            amount,
            options,
            &mut cancelled,
        )?;
        Self::from_equilibrium(
            basis,
            options,
            solvent_volume,
            stocks,
            equilibrium,
            None,
            &mut cancelled,
        )
    }

    /// Actual dissolved particles, including neutral complexes, with no solid.
    pub const fn solution(&self) -> &AqueousSolution {
        &self.solution
    }
    /// H+, OH-, acid forms, aqueous spectators, then neutral metal members.
    /// Every identity occurs exactly once and remains fixed across updates.
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }
    pub fn solvent(&self) -> &SpeciesId {
        self.aqueous().model().water().water()
    }
    pub const fn basis(&self) -> SolventActivityBasis {
        self.basis
    }
    pub const fn options(&self) -> CopperSulfideOptions {
        self.options
    }
    /// Original selected phase, all source receipts and simultaneous diagnostics.
    pub const fn equilibrium(&self) -> &CopperSulfideSpeciation {
        &self.equilibrium
    }
    pub fn model(&self) -> &CopperSulfidePrecipitation {
        self.equilibrium.model()
    }
    /// Original coupled metal distributions and free/protonated/bound ligand.
    pub const fn ligand(&self) -> &LigandSpeciation {
        &self.equilibrium.ligand
    }
    /// Original shared pH, acid families and charged metal members.
    /// Neutral metal members are additionally available through `ligand()`.
    pub const fn aqueous(&self) -> &AcidBaseSpeciation {
        &self.equilibrium.ligand.aqueous
    }
    pub const fn solid_amount(&self) -> Amount {
        self.equilibrium.solid_amount
    }
    pub const fn dissolved_copper(&self) -> Amount {
        self.equilibrium.dissolved_copper
    }
    pub const fn dissolved_sulfide(&self) -> Amount {
        self.equilibrium.dissolved_sulfide
    }
    /// Complete free + protonated + metal-bound ligand inventory.
    pub fn ligand_total(&self) -> Amount {
        self.stocks.ligand_total
    }
    /// Complete metal components in canonical family order, including solid Cu.
    pub fn metal_totals(&self) -> &[(String, Amount)] {
        &self.stocks.metal_totals
    }
    /// Complete non-ligand acid families, including the sulfide in solid Cu2S.
    pub fn other_acid_totals(&self) -> &[(String, Amount)] {
        &self.stocks.other_acid_totals
    }
    /// Original fixed spectators only, excluding dynamically speciated metals.
    pub fn spectators(&self) -> &[SpectatorIon] {
        &self.stocks.spectators
    }

    /// Re-equilibrate solvent exchange with every analytical inventory retained.
    /// Source activity amount scales through the original explicit volume basis.
    pub fn with_solvent_volume_m3(&self, solvent_volume_m3: f64) -> Result<Self> {
        self.with_solvent_volume_m3_cancelled(solvent_volume_m3, || false)
    }
    pub fn with_solvent_volume_m3_cancelled(
        &self,
        solvent_volume_m3: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        self.resolve(
            Volume::new(solvent_volume_m3),
            Arc::clone(&self.stocks),
            &mut cancelled,
        )
    }

    /// Re-equilibrate caller-declared analytical dosing or withdrawal.
    ///
    /// Arrays follow `metal_totals()`, `other_acid_totals()` and `spectators()`.
    /// Identities, spectator charges and original sources cannot be replaced.
    /// New metal/acid totals include retained solid components; ligand counts
    /// all uncomplexed and bound forms. No counterion or reagent is inferred.
    pub fn with_analytical_totals(
        &self,
        ligand_total: Amount,
        metal_totals: &[Amount],
        other_acid_totals: &[Amount],
        spectator_amounts: &[Amount],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if metal_totals.len() != self.stocks.metal_totals.len()
            || other_acid_totals.len() != self.stocks.other_acid_totals.len()
            || spectator_amounts.len() != self.stocks.spectators.len()
        {
            return Err(ElectrolyteSolubilityError::Invalid(
                "updated amounts must retain original metal, acid and spectator layouts",
            ));
        }
        let stocks = Arc::new(AnalyticalStocks {
            ligand_total,
            metal_totals: self
                .stocks
                .metal_totals
                .iter()
                .zip(metal_totals)
                .map(|((id, _), &amount)| (id.clone(), amount))
                .collect(),
            other_acid_totals: self
                .stocks
                .other_acid_totals
                .iter()
                .zip(other_acid_totals)
                .map(|((id, _), &amount)| (id.clone(), amount))
                .collect(),
            spectators: self
                .stocks
                .spectators
                .iter()
                .zip(spectator_amounts)
                .map(|(ion, &amount)| SpectatorIon {
                    species: ion.species.clone(),
                    charge: ion.charge,
                    amount,
                })
                .collect(),
        });
        self.resolve(
            Volume::new(self.solution.solvent_volume_m3()),
            stocks,
            &mut cancelled,
        )
    }

    fn resolve(
        &self,
        volume: Volume,
        stocks: Arc<AnalyticalStocks>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let amount = self.basis.amount_at_unit_activity(volume)?;
        let equilibrium = self.model().speciate(
            stocks.ligand_total,
            &stocks.metal_totals,
            &stocks.other_acid_totals,
            &stocks.spectators,
            amount,
            self.options,
            &mut *cancelled,
        )?;
        Self::from_equilibrium(
            self.basis,
            self.options,
            volume,
            stocks,
            equilibrium,
            Some(Arc::clone(&self.species)),
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_equilibrium(
        basis: SolventActivityBasis,
        options: CopperSulfideOptions,
        volume: Volume,
        stocks: Arc<AnalyticalStocks>,
        equilibrium: CopperSulfideSpeciation,
        species_axis: Option<Arc<[SpeciesId]>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let neutral: Vec<_> = equilibrium
            .ligand
            .metals
            .iter()
            .flat_map(|metal| &metal.members)
            .filter(|member| member.charge == 0)
            .map(|member| (&member.species, member.amount))
            .collect();
        let (species, solution) = publish_aqueous_with_additional_particles(
            basis,
            volume,
            &equilibrium.ligand.aqueous,
            &[equilibrium.model().solid()],
            &neutral,
            equilibrium.model().standard().temperature(),
            species_axis,
            cancelled,
        )?;
        Ok(Self {
            basis,
            options,
            stocks,
            equilibrium,
            species,
            solution,
        })
    }
}

fn admit_layout(
    model: &CopperSulfidePrecipitation,
    metals: usize,
    other_acids: usize,
    spectators: usize,
) -> Result<()> {
    let ligand = model.ligand_model();
    let aqueous = ligand.aqueous_model();
    if metals != ligand.families().len() || other_acids != aqueous.families().len() - 1 {
        return Err(ElectrolyteSolubilityError::Invalid(
            "supply every original metal and non-ligand acid family exactly once",
        ));
    }
    let count = 2
        + aqueous
            .families()
            .iter()
            .map(|family| family.species().len())
            .sum::<usize>()
        + ligand
            .families()
            .iter()
            .map(|family| family.species().len())
            .sum::<usize>();
    if count.saturating_add(spectators) > MAX_OSMOTIC_SOLUTES {
        return Err(ElectrolyteSolubilityError::Invalid(
            "complexed material requires at most 64 original dissolved particles",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
