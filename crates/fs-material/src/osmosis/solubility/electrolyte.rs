//! Source-backed mineral solubility with the actual speciated aqueous particles.
//!
//! The original hydrolysis/salt owner solves the selected finite solid, pH,
//! every acid/base family and charge balance together. This adapter retains
//! that complete result and the original TOTAL analytical inventories. Its
//! shared ideal aqueous view counts only actual dissolved species: each H+,
//! OH-, family member and spectator once. Neither undissolved formula units
//! nor a second inferred dissociation multiplier become osmotic particles.
//!
//! Physical solvent volume and the amount at unit activity are distinct. An
//! explicit fixed reference-amount density converts volume to the source's
//! activity basis; molality data require the caller's solvent-mass conversion.
//! The original fixed-temperature, fixed-pressure, ideal dilute chemistry uses
//! unit water activity. The aqueous view supplies the leading dilute
//! particle pressure/activity approximation, not a concentrated electrolyte
//! activity model. It retains the pure-water H+/OH- background: its reference
//! osmolarity is `2 * sqrt(K_w) * reference_amount_per_m3`. Thus empty analytical
//! stocks still have this small positive particle pressure and an activity
//! slightly below one; this view does not exactly normalize to pure solvent.
//! Reaction water, solute/crystal partial volumes, heat and competing solids
//! remain outside this adapter. Explicit finite-rate crystal inventories use
//! [`kinetics`]; selected Cu/S/CN metal association uses [`complexation`]. Both
//! retain their original chemical owners and share aqueous publication.

pub mod complexation;
pub mod kinetics;

use std::sync::Arc;

use fs_qty::{Amount, AmountConcentration, Volume};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        ConservedReaction,
        thermodynamics::{
            StandardReactionGibbs,
            acid_base::{
                AcidBaseError, AcidBaseSpeciation, IdealAcidBaseSystem, SpectatorIon,
                complexation::LigandSpeciationError,
                precipitation::{
                    HydrolysisPrecipitation, PrecipitationOptions, PrecipitationSpeciation,
                    PrecipitationStatus,
                    salt::{SaltPrecipitation, SaltSpeciation},
                },
            },
        },
    },
};

use super::super::{AqueousSolution, IdealDiluteSolution, MAX_OSMOTIC_SOLUTES, OsmosisError};

/// Original chemistry or dissolved-material refusal; no partial state is published.
#[derive(Clone, Debug, PartialEq)]
pub enum ElectrolyteSolubilityError {
    Speciation(AcidBaseError),
    Complexation(LigandSpeciationError),
    Aqueous(OsmosisError),
    Invalid(&'static str),
}
impl core::fmt::Display for ElectrolyteSolubilityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electrolyte solubility material: {self:?}")
    }
}
impl std::error::Error for ElectrolyteSolubilityError {}
impl From<AcidBaseError> for ElectrolyteSolubilityError {
    fn from(error: AcidBaseError) -> Self {
        Self::Speciation(error)
    }
}
impl From<OsmosisError> for ElectrolyteSolubilityError {
    fn from(error: OsmosisError) -> Self {
        Self::Aqueous(error)
    }
}
impl From<LigandSpeciationError> for ElectrolyteSolubilityError {
    fn from(error: LigandSpeciationError) -> Self {
        Self::Complexation(error)
    }
}
type Result<T> = core::result::Result<T, ElectrolyteSolubilityError>;

/// Explicit incompressible solvent and source-activity volume convention.
///
/// `reference_amount_per_m3 * V` is the common amount n_ref at unit activity.
/// It is not the solute concentration or a salt dissociation factor. For a
/// source using unit molality, supply solvent density times 1 mol/kg; this
/// conversion cannot be inferred from a species label or from molar volume alone.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolventActivityBasis {
    solvent_molar_volume_m3_mol: f64,
    reference_amount_per_m3: AmountConcentration,
}
impl SolventActivityBasis {
    pub fn new(
        solvent_molar_volume_m3_mol: f64,
        reference_amount_per_m3: AmountConcentration,
    ) -> Result<Self> {
        if !positive(solvent_molar_volume_m3_mol) || !positive(reference_amount_per_m3.value()) {
            return Err(ElectrolyteSolubilityError::Invalid(
                "positive finite solvent molar volume and activity-reference density required",
            ));
        }
        Ok(Self {
            solvent_molar_volume_m3_mol,
            reference_amount_per_m3,
        })
    }
    pub const fn solvent_molar_volume_m3_mol(&self) -> f64 {
        self.solvent_molar_volume_m3_mol
    }
    pub const fn reference_amount_per_m3(&self) -> AmountConcentration {
        self.reference_amount_per_m3
    }
    /// Original reference density applied to an actual finite solvent volume.
    pub fn amount_at_unit_activity(&self, solvent_volume: Volume) -> Result<Amount> {
        let amount = self.reference_amount_per_m3.value() * solvent_volume.value();
        if !positive(solvent_volume.value()) || !positive(amount) {
            return Err(ElectrolyteSolubilityError::Invalid(
                "positive finite solvent volume and reference amount must be representable",
            ));
        }
        Ok(Amount::new(amount))
    }
}

/// Original selected-solid result, including source, pH and balance diagnostics.
#[derive(Clone, Debug)]
pub enum ElectrolytePhaseEquilibrium {
    Hydrolysis(PrecipitationSpeciation),
    Salt(SaltSpeciation),
}
impl ElectrolytePhaseEquilibrium {
    pub const fn aqueous(&self) -> &AcidBaseSpeciation {
        match self {
            Self::Hydrolysis(result) => &result.aqueous,
            Self::Salt(result) => &result.aqueous,
        }
    }
    pub fn reaction(&self) -> &ConservedReaction {
        match self {
            Self::Hydrolysis(result) => result.model().reaction(),
            Self::Salt(result) => result.model().reaction(),
        }
    }
    pub fn standard(&self) -> &StandardReactionGibbs {
        match self {
            Self::Hydrolysis(result) => result.model().standard(),
            Self::Salt(result) => result.model().standard(),
        }
    }
    pub fn solid(&self) -> &SpeciesId {
        match self {
            Self::Hydrolysis(result) => result.model().solid(),
            Self::Salt(result) => result.model().solid(),
        }
    }
    pub const fn solid_amount(&self) -> Amount {
        match self {
            Self::Hydrolysis(result) => result.solid_amount,
            Self::Salt(result) => result.solid_amount,
        }
    }
    /// Complete cation-component stock, including its stoichiometric solid share.
    pub const fn component_total(&self) -> Amount {
        match self {
            Self::Hydrolysis(result) => result.component_total,
            Self::Salt(result) => result.cation_total,
        }
    }
    pub const fn status(&self) -> PrecipitationStatus {
        match self {
            Self::Hydrolysis(result) => result.status,
            Self::Salt(result) => result.status,
        }
    }
}

#[derive(Clone, Debug)]
struct AnalyticalStocks {
    component_total: Amount,
    family_totals: Box<[(String, Amount)]>,
    spectators: Box<[SpectatorIon]>,
}
impl AnalyticalStocks {
    fn new(
        component_total: Amount,
        mut family_totals: Vec<(String, Amount)>,
        mut spectators: Vec<SpectatorIon>,
    ) -> Self {
        family_totals.sort_by(|left, right| left.0.cmp(&right.0));
        spectators.sort_by(|left, right| left.species.cmp(&right.species));
        Self {
            component_total,
            family_totals: family_totals.into_boxed_slice(),
            spectators: spectators.into_boxed_slice(),
        }
    }
}

/// One finite mineral reservoir and its source-backed speciated aqueous view.
///
/// Original TOTAL family/cation stocks remain distinct from the dissolved-only
/// aqueous report. Immutable volume trials share those analytical stocks and
/// the dissolved species axis; all chemistry remains in the original L1 owner.
#[derive(Clone, Debug)]
pub struct ElectrolyteSolubilitySolution {
    basis: SolventActivityBasis,
    options: PrecipitationOptions,
    stocks: Arc<AnalyticalStocks>,
    equilibrium: ElectrolytePhaseEquilibrium,
    species: Arc<[SpeciesId]>,
    solution: AqueousSolution,
}
impl ElectrolyteSolubilitySolution {
    /// Bind one original hydrolyzing solid. Its fixed-solvent water convention
    /// is preserved; this adapter does not debit reaction water from the volume.
    #[allow(clippy::too_many_arguments)]
    pub fn from_hydrolysis(
        model: HydrolysisPrecipitation,
        basis: SolventActivityBasis,
        solvent_volume: Volume,
        component_total: Amount,
        family_totals: Vec<(String, Amount)>,
        spectators: Vec<SpectatorIon>,
        options: PrecipitationOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        admit_layout(model.aqueous_model(), spectators.len())?;
        let amount = basis.amount_at_unit_activity(solvent_volume)?;
        let stocks = Arc::new(AnalyticalStocks::new(
            component_total,
            family_totals,
            spectators,
        ));
        let result = model.speciate(
            stocks.component_total,
            &stocks.family_totals,
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
            ElectrolytePhaseEquilibrium::Hydrolysis(result),
            None,
            &mut cancelled,
        )
    }

    /// Bind the original finite salt and complete acid/base families. Analytical
    /// family/cation totals include their full stoichiometric solid contribution.
    #[allow(clippy::too_many_arguments)]
    pub fn from_salt(
        model: SaltPrecipitation,
        basis: SolventActivityBasis,
        solvent_volume: Volume,
        component_total: Amount,
        family_totals: Vec<(String, Amount)>,
        spectators: Vec<SpectatorIon>,
        options: PrecipitationOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        admit_layout(model.aqueous_model(), spectators.len())?;
        let amount = basis.amount_at_unit_activity(solvent_volume)?;
        let stocks = Arc::new(AnalyticalStocks::new(
            component_total,
            family_totals,
            spectators,
        ));
        let result = model.speciate(
            stocks.component_total,
            &stocks.family_totals,
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
            ElectrolytePhaseEquilibrium::Salt(result),
            None,
            &mut cancelled,
        )
    }

    /// Actual dissolved particle amounts; solid formula units are absent.
    /// The H+/OH- pure-water background is retained, so the ideal particle
    /// pressure/activity is a leading dilute approximation rather than an
    /// exactly normalized pure-solvent reference; see the module scope.
    pub const fn solution(&self) -> &AqueousSolution {
        &self.solution
    }
    /// H+, OH-, family forms, then spectators, each on its original identity.
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }
    pub fn solvent(&self) -> &SpeciesId {
        self.equilibrium.aqueous().model().water().water()
    }
    pub const fn basis(&self) -> SolventActivityBasis {
        self.basis
    }
    pub const fn options(&self) -> PrecipitationOptions {
        self.options
    }
    pub const fn equilibrium(&self) -> &ElectrolytePhaseEquilibrium {
        &self.equilibrium
    }
    pub fn component_total(&self) -> Amount {
        self.stocks.component_total
    }
    /// Complete analytical totals in canonical family order, including solid.
    pub fn family_totals(&self) -> &[(String, Amount)] {
        &self.stocks.family_totals
    }
    /// Original fixed spectators only, excluding the dynamically solved cation.
    pub fn spectators(&self) -> &[SpectatorIon] {
        &self.stocks.spectators
    }

    /// Re-equilibrate pure-solvent exchange with all analytical stocks retained.
    /// The physical volume also scales the original source-activity basis.
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

    /// Re-equilibrate explicitly changed analytical amounts in the SAME family
    /// and spectator order. Every original identity, spectator charge and source
    /// remains bound. The new cation/family totals include retained solid.
    pub fn with_analytical_totals(
        &self,
        component_total: Amount,
        family_totals: &[Amount],
        spectator_amounts: &[Amount],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if family_totals.len() != self.stocks.family_totals.len()
            || spectator_amounts.len() != self.stocks.spectators.len()
        {
            return Err(ElectrolyteSolubilityError::Invalid(
                "updated amounts must retain the original family and spectator layout",
            ));
        }
        let stocks = Arc::new(AnalyticalStocks {
            component_total,
            family_totals: self
                .stocks
                .family_totals
                .iter()
                .zip(family_totals)
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
        let equilibrium = match &self.equilibrium {
            ElectrolytePhaseEquilibrium::Hydrolysis(result) => {
                ElectrolytePhaseEquilibrium::Hydrolysis(result.model().speciate(
                    stocks.component_total,
                    &stocks.family_totals,
                    &stocks.spectators,
                    amount,
                    self.options,
                    &mut *cancelled,
                )?)
            }
            ElectrolytePhaseEquilibrium::Salt(result) => {
                ElectrolytePhaseEquilibrium::Salt(result.model().speciate(
                    stocks.component_total,
                    &stocks.family_totals,
                    &stocks.spectators,
                    amount,
                    self.options,
                    &mut *cancelled,
                )?)
            }
        };
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
        options: PrecipitationOptions,
        volume: Volume,
        stocks: Arc<AnalyticalStocks>,
        equilibrium: ElectrolytePhaseEquilibrium,
        species_axis: Option<Arc<[SpeciesId]>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let (species, solution) = publish_aqueous(
            basis,
            volume,
            equilibrium.aqueous(),
            equilibrium.solid(),
            equilibrium.standard().temperature(),
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

/// Publish exactly the original speciated particles for either a frozen finite
/// kinetic phase or the original equilibrium owner. No solid-phase calculation
/// or additional aqueous equilibrium is performed at this material boundary.
fn publish_aqueous(
    basis: SolventActivityBasis,
    volume: Volume,
    aqueous: &AcidBaseSpeciation,
    solid: &SpeciesId,
    temperature: fs_qty::Temperature,
    species_axis: Option<Arc<[SpeciesId]>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Arc<[SpeciesId]>, AqueousSolution)> {
    publish_aqueous_with_additional_particles(
        basis, volume, aqueous, &[solid], &[], temperature, species_axis, cancelled,
    )
}

/// Publish the original charged/acid-base state plus explicitly resolved
/// neutral particles. Metal complexes with nonzero charge are already present
/// in the native aqueous spectator axis and must not be appended again.
#[allow(clippy::too_many_arguments)]
fn publish_aqueous_with_additional_particles(
    basis: SolventActivityBasis,
    volume: Volume,
    aqueous: &AcidBaseSpeciation,
    solids: &[&SpeciesId],
    additional_particles: &[(&SpeciesId, Amount)],
    temperature: fs_qty::Temperature,
    species_axis: Option<Arc<[SpeciesId]>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Arc<[SpeciesId]>, AqueousSolution)> {
    poll(cancelled)?;
    let count = 2
        + aqueous
            .families
            .iter()
            .map(|family| family.forms.len())
            .sum::<usize>()
        + aqueous.spectators.len()
        + additional_particles.len();
    if count > MAX_OSMOTIC_SOLUTES
        || species_axis
            .as_ref()
            .is_some_and(|axis| axis.len() != count)
    {
        return Err(ElectrolyteSolubilityError::Invalid(
            "speciated material requires at most 64 particles in its original layout",
        ));
    }
    let mut amounts = Vec::with_capacity(count);
    let mut identities = if species_axis.is_none() {
        Vec::with_capacity(count)
    } else {
        Vec::new()
    };
    let water = aqueous.model().water();
    let mut append = |id: &SpeciesId, amount: Amount| -> Result<()> {
        poll(cancelled)?;
        if id == water.water() || solids.contains(&id) {
            return Err(ElectrolyteSolubilityError::Invalid(
                "solvent or solid cannot be counted as dissolved osmotic particles",
            ));
        }
        if let Some(axis) = &species_axis {
            if axis.get(amounts.len()) != Some(id) {
                return Err(ElectrolyteSolubilityError::Invalid(
                    "speciation changed the original dissolved species ordering",
                ));
            }
        } else if identities.contains(id) {
            return Err(ElectrolyteSolubilityError::Invalid(
                "a dissolved species cannot contribute osmotic particles twice",
            ));
        } else {
            identities.push(id.clone());
        }
        amounts.push(amount.value());
        Ok(())
    };
    append(water.proton(), aqueous.proton_amount)?;
    append(water.hydroxide(), aqueous.hydroxide_amount)?;
    for family in &aqueous.families {
        for form in &family.forms {
            append(&form.species, form.amount)?;
        }
    }
    for ion in &aqueous.spectators {
        append(&ion.species, ion.amount)?;
    }
    for &(id, amount) in additional_particles {
        append(id, amount)?;
    }
    let solution: AqueousSolution = IdealDiluteSolution::new(
        temperature.value(),
        basis.solvent_molar_volume_m3_mol,
        volume.value(),
        amounts,
    )?
    .into();
    solution.water_activity()?;
    poll(cancelled)?;
    Ok((species_axis.unwrap_or_else(|| identities.into()), solution))
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(AcidBaseError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn admit_layout(model: &IdealAcidBaseSystem, spectators: usize) -> Result<()> {
    let count = model
        .families()
        .iter()
        .try_fold(3usize, |n, family| n.checked_add(family.species().len()));
    if count
        .and_then(|n| n.checked_add(spectators))
        .is_none_or(|n| n > MAX_OSMOTIC_SOLUTES)
    {
        return Err(ElectrolyteSolubilityError::Invalid(
            "at most 64 dissolved species including proton, hydroxide and cation are supported",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
