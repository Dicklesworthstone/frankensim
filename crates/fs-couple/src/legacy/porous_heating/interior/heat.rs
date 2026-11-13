//! Original thermal bindings and accepted-step accounts shared by frozen and
//! temperature-refreshed interior fields. Acoustic solves never occur here.

use super::{
    BoundCell, BoundWall, CellId, DepositionError, DepositionSimulation, EnergyAudit,
    InteriorCellHeating, InteriorHeatingEnergy, InteriorHeatingError, InteriorHeatingStep,
    InteriorWallHeating, finite,
};
use fs_math::dd::Dd;
use std::collections::BTreeMap;

/// Immutable, physically bound sources from one complete acoustic evaluation.
pub(super) struct InteriorHeatState<'t> {
    pub(super) cells: Vec<BoundCell<'t>>,
    pub(super) walls: Vec<BoundWall<'t>>,
    pub(super) sources_w: BTreeMap<CellId, f64>,
}

impl InteriorHeatState<'_> {
    pub(super) fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), InteriorHeatingError> {
        for binding in self
            .cells
            .iter()
            .map(|bound| &bound.binding)
            .chain(self.walls.iter().map(|bound| &bound.binding))
        {
            checkpoint()?;
            binding.validate(thermal)?;
        }
        Ok(())
    }

    pub(super) fn cells(&self) -> impl ExactSizeIterator<Item = &InteriorCellHeating> {
        self.cells.iter().map(|cell| &cell.heat)
    }

    pub(super) fn walls(&self) -> impl ExactSizeIterator<Item = &InteriorWallHeating> {
        self.walls.iter().map(|wall| &wall.heat)
    }

    pub(super) fn sources(&self) -> BTreeMap<CellId, f64> {
        self.sources_w.clone()
    }
}

/// Mutable energy history; cloning never shares partially accepted work.
#[derive(Clone)]
pub(super) struct InteriorHeatAccount {
    cell_energy: Vec<Dd>,
    wall_energy: Vec<Dd>,
    mechanisms: [Dd; 4],
}

impl InteriorHeatAccount {
    pub(super) fn new(heat: &InteriorHeatState<'_>) -> Self {
        Self {
            cell_energy: vec![Dd::ZERO; heat.cells.len()],
            wall_energy: vec![Dd::ZERO; heat.walls.len()],
            mechanisms: [Dd::ZERO; 4],
        }
    }

    fn validate_shape(&self, heat: &InteriorHeatState<'_>) -> Result<(), InteriorHeatingError> {
        if self.cell_energy.len() != heat.cells.len() || self.wall_energy.len() != heat.walls.len()
        {
            return Err(InteriorHeatingError::Invalid(
                "interior heat account needs the original volume and wall destinations",
            ));
        }
        Ok(())
    }

    /// Charge the accepted-start powers before any endpoint acoustic refresh.
    pub(super) fn record(
        &mut self,
        heat: &InteriorHeatState<'_>,
        dt: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), InteriorHeatingError> {
        self.validate_shape(heat)?;
        for (cell, energy) in heat.cells.iter().zip(&mut self.cell_energy) {
            checkpoint()?;
            *energy = *energy + Dd::from_f64(finite(cell.heat.absorbed_power_w() * dt)?);
            for (sum, power) in self.mechanisms.iter_mut().zip([
                cell.heat.density_dissipation_w,
                cell.heat.compressional_dissipation_w,
            ]) {
                *sum = *sum + Dd::from_f64(finite(power * dt)?);
            }
        }
        for (wall, energy) in heat.walls.iter().zip(&mut self.wall_energy) {
            checkpoint()?;
            *energy = *energy + Dd::from_f64(finite(wall.heat.absorbed_power_w() * dt)?);
            for (sum, power) in self.mechanisms[2..].iter_mut().zip([
                wall.heat.resistance_dissipation_w,
                wall.heat.relaxation_dissipation_w,
            ]) {
                *sum = *sum + Dd::from_f64(finite(power * dt)?);
            }
        }
        Ok(())
    }

    pub(super) fn energy(&self) -> Result<InteriorHeatingEnergy, InteriorHeatingError> {
        let deposited = self
            .cell_energy
            .iter()
            .chain(&self.wall_energy)
            .copied()
            .fold(Dd::ZERO, |a, b| a + b);
        Ok(InteriorHeatingEnergy {
            density_dissipation_j: finite(self.mechanisms[0].to_f64())?,
            compressional_dissipation_j: finite(self.mechanisms[1].to_f64())?,
            wall_resistance_dissipation_j: finite(self.mechanisms[2].to_f64())?,
            wall_relaxation_dissipation_j: finite(self.mechanisms[3].to_f64())?,
            deposited_j: finite(deposited.to_f64())?,
        })
    }

    pub(super) fn cell_energy_j(
        &self,
        heat: &InteriorHeatState<'_>,
    ) -> Result<Vec<(CellId, f64)>, InteriorHeatingError> {
        self.validate_shape(heat)?;
        heat.cells
            .iter()
            .zip(&self.cell_energy)
            .map(|(cell, energy)| Ok((cell.heat.cell, finite(energy.to_f64())?)))
            .collect()
    }

    pub(super) fn wall_energy_j(
        &self,
        heat: &InteriorHeatState<'_>,
    ) -> Result<Vec<(CellId, f64)>, InteriorHeatingError> {
        self.validate_shape(heat)?;
        heat.walls
            .iter()
            .zip(&self.wall_energy)
            .map(|(wall, energy)| Ok((wall.heat.cell, finite(energy.to_f64())?)))
            .collect()
    }

    pub(super) fn finish<'t>(
        &self,
        heat: &InteriorHeatState<'_>,
        thermal: DepositionSimulation<'t>,
        heat_report: EnergyAudit,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<InteriorHeatingStep<'t>, InteriorHeatingError> {
        heat.validate(&thermal, checkpoint)?;
        let acoustic = self.energy()?;
        let cell_energy_j = self.cell_energy_j(heat)?;
        let wall_energy_j = self.wall_energy_j(heat)?;
        checkpoint()?;
        Ok(InteriorHeatingStep {
            thermal,
            heat_report,
            acoustic,
            cell_energy_j,
            wall_energy_j,
        })
    }
}
