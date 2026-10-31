//! Shared finite aqueous inventory with an explicit ideal or nonideal water law.
//!
//! Existing ideal solutions convert without allocation. A nonideal solution
//! retains the same inventory owner while sharing its immutable virial model
//! through `Arc`. State replacement checks the full candidate before publication.
//! Neither variant adds ionic transport, partial solute volumes, or caloric state.

use std::sync::Arc;

use super::nonideal::VirialSolutionModel;
use super::{IdealDiluteSolution, Result, finite, shape, thermal_scale, validate_inventory};

/// Finite independent-solute amounts and the thermodynamic law that interprets them.
#[derive(Clone, Debug, PartialEq)]
pub struct AqueousSolution {
    inventory: IdealDiluteSolution,
    model: Option<Arc<VirialSolutionModel>>,
}

impl From<IdealDiluteSolution> for AqueousSolution {
    fn from(inventory: IdealDiluteSolution) -> Self {
        Self {
            inventory,
            model: None,
        }
    }
}

impl AqueousSolution {
    /// Prepare a nonideal inventory using the model's own fixed temperature,
    /// solvent molar volume, species order and declared concentration support.
    pub fn nonideal(
        model: Arc<VirialSolutionModel>,
        solvent_volume_m3: f64,
        solute_moles: Vec<f64>,
    ) -> Result<Self> {
        model.evaluate_inventories(solvent_volume_m3, &solute_moles)?;
        let inventory = IdealDiluteSolution::new(
            model.temperature_k(),
            model.solvent_molar_volume_m3_mol(),
            solvent_volume_m3,
            solute_moles,
        )?;
        Ok(Self {
            inventory,
            model: Some(model),
        })
    }

    pub const fn is_ideal(&self) -> bool {
        self.model.is_none()
    }
    pub const fn temperature_k(&self) -> f64 {
        self.inventory.temperature_k()
    }
    pub const fn solvent_molar_volume_m3_mol(&self) -> f64 {
        self.inventory.solvent_molar_volume_m3_mol()
    }
    pub const fn solvent_volume_m3(&self) -> f64 {
        self.inventory.solvent_volume_m3()
    }
    pub fn solute_moles(&self) -> &[f64] {
        self.inventory.solute_moles()
    }
    pub fn nonideal_model(&self) -> Option<&VirialSolutionModel> {
        self.model.as_deref()
    }

    /// Conservative owned/shared heap estimate, excluding `size_of::<Self>()`.
    /// A shared model is counted in full for EACH solution, including its Arc
    /// control words and an alignment allowance, so callers may overestimate
    /// shared storage but cannot omit it from a preparation workspace budget.
    pub fn heap_bytes(&self) -> usize {
        core::mem::size_of_val(self.inventory.solute_moles())
            + self.model.as_ref().map_or(0, |model| {
                core::mem::size_of::<VirialSolutionModel>()
                    + 2 * core::mem::size_of::<usize>()
                    + core::mem::align_of::<VirialSolutionModel>()
                    + model.heap_bytes()
            })
    }

    pub fn concentrations_into(&self, output: &mut [f64]) -> Result<()> {
        self.inventory.concentrations_into(output)
    }

    /// Pressure from the actual selected model at the accepted inventory.
    pub fn osmotic_pressure_pa(&self) -> Result<f64> {
        self.osmotic_pressure_at(self.solvent_volume_m3(), self.solute_moles())
    }

    /// Water activity from the actual selected model at the accepted inventory.
    pub fn water_activity(&self) -> Result<f64> {
        match &self.model {
            None => self.inventory.water_activity(),
            Some(model) => Ok(model
                .evaluate_inventories(self.solvent_volume_m3(), self.solute_moles())?
                .water_activity),
        }
    }

    /// Allocation-free pressure evaluation of a candidate in the SAME species layout.
    /// The accepted solution is borrowed and never modified by an unsuccessful trial.
    pub fn osmotic_pressure_at(&self, solvent_volume_m3: f64, solute_moles: &[f64]) -> Result<f64> {
        shape(self.solute_moles().len(), solute_moles.len())?;
        match &self.model {
            Some(model) => Ok(model
                .evaluate_inventories(solvent_volume_m3, solute_moles)?
                .osmotic_pressure_pa),
            None => {
                let concentration =
                    validate_inventory(self.temperature_k(), solvent_volume_m3, solute_moles)?;
                finite(
                    thermal_scale(self.temperature_k())? * concentration,
                    "osmotic pressure is not representable",
                )
            }
        }
    }

    /// New owned inventory after pure solvent exchange. Solute amounts and the
    /// shared thermodynamic model are preserved; this owned operation allocates
    /// once for the copied amount array, unlike trial evaluation/replacement.
    pub fn with_solvent_volume_m3(&self, solvent_volume_m3: f64) -> Result<Self> {
        self.osmotic_pressure_at(solvent_volume_m3, self.solute_moles())?;
        Ok(Self {
            inventory: self.inventory.with_solvent_volume_m3(solvent_volume_m3)?,
            model: self.model.clone(),
        })
    }

    /// Atomically replace a complete finite inventory without allocating.
    /// Both the actual thermodynamic support and the original inventory checks
    /// must admit the candidate before the accepted state changes.
    pub fn replace_inventories(
        &mut self,
        solvent_volume_m3: f64,
        solute_moles: &[f64],
    ) -> Result<()> {
        self.osmotic_pressure_at(solvent_volume_m3, solute_moles)?;
        self.inventory
            .replace_inventories(solvent_volume_m3, solute_moles)
    }
}

#[cfg(test)]
mod tests;
