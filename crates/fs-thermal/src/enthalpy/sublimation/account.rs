//! One relative-enthalpy balance for split and coupled solid/vapor timelines.

use super::coupled::CoupledSublimationAudit;
use super::{SurfaceStep, ThermalSublimationError};
use crate::enthalpy::observe::EnergySum;
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, Event, finite};
use std::collections::BTreeMap;

#[derive(Default)]
pub(super) struct TimelineAccount {
    reference_h_j_kg: f64,
    external_heat: EnergySum,
    mass_transfer: EnergySum,
    latent: EnergySum,
    vapor: EnergySum,
    relative_vapor: EnergySum,
    recession: EnergySum,
    phase_steps: usize,
    last_step: Option<SurfaceStep>,
}

impl TimelineAccount {
    pub(super) fn new(
        initial: &DepositionSimulation<'_>,
        id: CellId,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        let mut reference = initial.cells.get(&id).map(|cell| cell.specific_h_j_kg);
        for (index, event) in initial.schedule[initial.cursor..].iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if let Event::Birth(cell) = &event.event {
                if cell.id == id && reference.is_none() {
                    reference = Some(cell.specific_h_j_kg);
                }
            }
        }
        let reference_h_j_kg = reference.ok_or(DepositionError::InvalidInput(
            "sublimating cell is neither active nor scheduled",
        ))?;
        finite(initial.mass_kg())?;
        finite(initial.energy_j())?;
        Ok(Self {
            reference_h_j_kg,
            ..Self::default()
        })
    }

    pub(super) fn heat(&mut self, joules: f64) -> Result<(), DepositionError> {
        self.external_heat.add(joules)
    }

    pub(super) fn phase(&mut self, step: SurfaceStep) -> Result<(), DepositionError> {
        self.mass_transfer.add(step.transferred_mass_kg)?;
        self.latent.add(step.latent_energy_j)?;
        self.vapor.add(step.vapor_energy_j)?;
        self.relative_vapor.add(finite(
            step.transferred_mass_kg
                * finite(step.solid.specific_enthalpy_j_kg - self.reference_h_j_kg)?,
        )?)?;
        self.relative_vapor.add(step.latent_energy_j)?;
        self.recession.add(step.recession_m)?;
        self.phase_steps = self
            .phase_steps
            .checked_add(1)
            .ok_or(DepositionError::Arithmetic)?;
        self.last_step = Some(step);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn finish(
        self,
        initial: &DepositionSimulation<'_>,
        next: &DepositionSimulation<'_>,
        substeps: usize,
        events: usize,
        energy_tolerance_j: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<CoupledSublimationAudit, ThermalSublimationError> {
        // Cancel each original/birth datum before reducing. A separate mass
        // balance prevents inventory roundoff from becoming artificial heat.
        let mut births = BTreeMap::new();
        let mut incoming_mass = EnergySum::default();
        let mut incoming_energy = EnergySum::default();
        for (index, event) in next.schedule[initial.cursor..next.cursor]
            .iter()
            .enumerate()
        {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if let Event::Birth(cell) = &event.event {
                births.insert(cell.id, (cell.mass_kg, cell.specific_h_j_kg));
                incoming_mass.add(cell.mass_kg)?;
                incoming_energy.add(finite(cell.mass_kg * cell.specific_h_j_kg)?)?;
            }
        }
        let incoming_mass_kg = incoming_mass.value()?;
        let transferred_mass_kg = self.mass_transfer.value()?;
        let external_heat_j = self.external_heat.value()?;
        let mut mass_balance = EnergySum::default();
        let mut storage = EnergySum::default();
        let mut energy_balance = EnergySum::default();
        for (index, (&id, cell)) in next.cells.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let original = initial.cells.get(&id);
            let (mass, h) = match original {
                Some(old) => (old.mass_kg, old.specific_h_j_kg),
                None => *births.get(&id).ok_or(DepositionError::Arithmetic)?,
            };
            let relative_energy = finite(cell.mass_kg * finite(cell.specific_h_j_kg - h)?)?;
            energy_balance.add(relative_energy)?;
            storage.add(relative_energy)?;
            storage.add(if original.is_some() {
                finite(finite(cell.mass_kg - mass)? * h)?
            } else {
                finite(cell.mass_kg * h)?
            })?;
            mass_balance.add(finite(cell.mass_kg - mass)?)?;
        }
        mass_balance.add(transferred_mass_kg)?;
        energy_balance.add(self.relative_vapor.value()?)?;
        energy_balance.add(-external_heat_j)?;
        let mass_balance_residual_kg = mass_balance.value()?;
        let scale = finite(initial.mass_kg())?
            .abs()
            .max(finite(next.mass_kg())?.abs())
            .max(incoming_mass_kg.abs())
            .max(transferred_mass_kg.abs());
        let terms = initial.cells.len() as f64
            + next.cells.len() as f64
            + events as f64
            + self.phase_steps as f64
            + 1.0;
        let mass_tolerance = finite(128.0 * f64::EPSILON * scale * terms)?;
        if mass_balance_residual_kg.abs() > mass_tolerance {
            return Err(ThermalSublimationError::MassBalance {
                residual_kg: mass_balance_residual_kg,
                tolerance_kg: mass_tolerance,
            });
        }
        let energy_residual_j = energy_balance.value()?;
        if energy_residual_j.abs() > energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(energy_residual_j).into());
        }
        finite(next.energy_j())?;
        Ok(CoupledSublimationAudit {
            incoming_mass_kg,
            incoming_energy_j: incoming_energy.value()?,
            external_heat_j,
            transferred_mass_kg,
            latent_energy_j: self.latent.value()?,
            vapor_energy_j: self.vapor.value()?,
            recession_m: self.recession.value()?,
            stored_energy_change_j: storage.value()?,
            mass_balance_residual_kg,
            energy_residual_j,
            substeps,
            events,
            phase_steps: self.phase_steps,
            last_step: self.last_step,
        })
    }
}
