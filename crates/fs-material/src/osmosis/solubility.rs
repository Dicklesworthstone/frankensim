//! Finite neutral solutes in equilibrium with their own pure solid phases.
//!
//! Each caller-selected solvent/solute pair supplies an apparent dilute
//! saturation concentration through the original L1 solubility law. The total
//! amount is retained separately from its dissolved and solid parts. Removing
//! solvent precipitates excess solute; adding solvent redissolves the finite
//! solid stock. Only the dissolved particles contribute to the shared ideal
//! dilute osmotic pressure and water activity.
//!
//! The independent pair approximation excludes electrolyte dissociation,
//! common-ion effects, chemical reactions, mixed crystals, nonideal activities,
//! solute partial volumes and crystal water. Equilibrium is instantaneous: this
//! owner supplies neither nucleation nor dissolution rates or heat feedback.
//! Numerical pair parameters do not establish measured material authority.

use std::sync::Arc;

use fs_qty::{Amount, AmountConcentration, Temperature, Volume};
use fs_thermochem::{
    SpeciesId,
    solution::solubility::{
        IdealSolubilityV1, SolubilityErrorV1, SolubilityPartitionV1, partition_at_saturation,
    },
};

use super::{AqueousSolution, IdealDiluteSolution, MAX_OSMOTIC_SOLUTES, OsmosisError};

/// Original pair-law or aqueous refusal; no partial state is published.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteSolubilityError {
    /// The original thermochemical saturation/partition law refused.
    Solubility(SolubilityErrorV1),
    /// The shared dissolved-particle thermodynamics refused.
    Aqueous(OsmosisError),
    /// Invalid species layout or unsupported composition.
    Invalid(&'static str),
}
impl core::fmt::Display for FiniteSolubilityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite solubility solution: {self:?}")
    }
}
impl std::error::Error for FiniteSolubilityError {}
impl From<SolubilityErrorV1> for FiniteSolubilityError {
    fn from(value: SolubilityErrorV1) -> Self {
        Self::Solubility(value)
    }
}
impl From<OsmosisError> for FiniteSolubilityError {
    fn from(value: OsmosisError) -> Self {
        Self::Aqueous(value)
    }
}
type Result<T> = core::result::Result<T, FiniteSolubilityError>;

/// One immutable pair law, total inventory and fixed-temperature saturation.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityComponent {
    law: IdealSolubilityV1,
    total: Amount,
    saturation: AmountConcentration,
}
impl SolubilityComponent {
    /// Original species pair, parameters and temperature validity.
    pub const fn law(&self) -> &IdealSolubilityV1 {
        &self.law
    }
    /// Conserved total of dissolved solute plus its solid phase [mol].
    pub const fn total(&self) -> Amount {
        self.total
    }
    /// Pair-law concentration evaluated at this solution's temperature [mol/m³].
    pub const fn saturation_concentration(&self) -> AmountConcentration {
        self.saturation
    }
}

/// A finite solution and every original solid stock, in explicit species order.
///
/// Pair laws, totals and fixed-temperature coefficients share one immutable
/// allocation across volume trials. Each bounded trial partitions at most 64
/// components without iterative chemistry or repeated temperature exponentials.
/// The aqueous view contains DISSOLVED amounts only and cannot replace this
/// owner when a consumer needs to preserve solid inventories.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteSolubilitySolution {
    solution: AqueousSolution,
    components: Arc<[SolubilityComponent]>,
    partitions: Box<[SolubilityPartitionV1]>,
}
impl FiniteSolubilitySolution {
    /// Admit 1..=64 independent neutral solutes sharing exactly one solvent.
    /// The supplied amounts are total component inventories, not concentrations
    /// and not necessarily already dissolved. One solute identity occurs once.
    /// Each dissolved neutral molecule counts as one osmotic particle.
    pub fn new(
        temperature: Temperature,
        solvent_molar_volume_m3_mol: f64,
        solvent_volume: Volume,
        components: Vec<(IdealSolubilityV1, Amount)>,
    ) -> Result<Self> {
        if components.is_empty() || components.len() > MAX_OSMOTIC_SOLUTES {
            return Err(FiniteSolubilityError::Invalid(
                "one through 64 independent soluble components required",
            ));
        }
        let solvent = components[0].0.solvent();
        for (index, (law, _)) in components.iter().enumerate() {
            if law.solvent() != solvent {
                return Err(FiniteSolubilityError::Invalid(
                    "all solubility pairs must use the same solvent identity",
                ));
            }
            if components[..index]
                .iter()
                .any(|(prior, _)| prior.solute() == law.solute())
            {
                return Err(FiniteSolubilityError::Invalid(
                    "each soluble component must occur exactly once",
                ));
            }
        }
        let prepared = components
            .into_iter()
            .map(|(law, total)| {
                let saturation = law.saturation_concentration(temperature)?;
                Ok(SolubilityComponent {
                    law,
                    total,
                    saturation,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Self::from_prepared(
            temperature,
            solvent_molar_volume_m3_mol,
            solvent_volume,
            prepared.into(),
        )
    }

    fn from_prepared(
        temperature: Temperature,
        solvent_molar_volume_m3_mol: f64,
        solvent_volume: Volume,
        components: Arc<[SolubilityComponent]>,
    ) -> Result<Self> {
        let partitions = components
            .iter()
            .map(|component| {
                partition_at_saturation(component.saturation, solvent_volume, component.total)
            })
            .collect::<core::result::Result<Vec<_>, _>>()?;
        let amounts = partitions
            .iter()
            .map(|partition| partition.dissolved.value())
            .collect();
        let solution: AqueousSolution = IdealDiluteSolution::new(
            temperature.value(),
            solvent_molar_volume_m3_mol,
            solvent_volume.value(),
            amounts,
        )?
        .into();
        solution.water_activity()?;
        Ok(Self {
            solution,
            components,
            partitions: partitions.into_boxed_slice(),
        })
    }

    /// Original common solvent identity, never guessed from a material label.
    pub fn solvent(&self) -> &SpeciesId {
        self.components[0].law.solvent()
    }
    /// Dissolved-only ideal thermodynamics at the accepted phase partition.
    pub const fn solution(&self) -> &AqueousSolution {
        &self.solution
    }
    /// Original laws and conserved component totals in caller-declared order.
    pub fn components(&self) -> &[SolubilityComponent] {
        &self.components
    }
    /// Dissolved and solid stocks, capacity and direct balances in that same order.
    pub fn phase_partitions(&self) -> &[SolubilityPartitionV1] {
        &self.partitions
    }

    /// Re-equilibrate after pure-solvent exchange, retaining each TOTAL amount.
    /// Concentrating an already saturated solution creates solid instead of
    /// incorrectly adding osmotic particles; dilution can exhaust that solid.
    pub fn with_solvent_volume_m3(&self, solvent_volume_m3: f64) -> Result<Self> {
        Self::from_prepared(
            Temperature::new(self.solution.temperature_k()),
            self.solution.solvent_molar_volume_m3_mol(),
            Volume::new(solvent_volume_m3),
            Arc::clone(&self.components),
        )
    }

    /// Re-equilibrate at a new externally prescribed temperature within every
    /// pair law's domain. The explicitly incompressible solvent molar volume is
    /// unchanged. Heat of dissolution is not applied to a thermal reservoir.
    pub fn with_temperature(&self, temperature: Temperature) -> Result<Self> {
        Self::new(
            temperature,
            self.solution.solvent_molar_volume_m3_mol(),
            Volume::new(self.solution.solvent_volume_m3()),
            self.components
                .iter()
                .map(|component| (component.law.clone(), component.total))
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests;
