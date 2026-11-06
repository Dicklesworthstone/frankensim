//! One physical open-mass account shared by explicit/implicit heat timelines.

use super::{CoupledEvaporationAudit, SurfaceStep, ThermalEvaporationError};
use crate::enthalpy::observe::EnergySum;
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, EnergyAudit, Event, finite};
use std::collections::BTreeMap;

#[derive(Default)]
pub(super) struct TimelineAccount {
    births: BTreeMap<CellId, (f64, f64)>,
    incoming_mass: EnergySum,
    incoming_energy: EnergySum,
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
    pub(super) fn birth(&mut self, id: CellId, mass: f64, h: f64) -> Result<(), DepositionError> {
        self.births.insert(id, (mass, h));
        self.incoming_mass.add(mass)?;
        self.incoming_energy.add(finite(mass * h)?)
    }

    pub(super) fn heat(&mut self, joules: f64) -> Result<(), DepositionError> {
        self.external_heat.add(joules)
    }

    /// Retain births and external heat from an accepted native heat interval.
    /// The split phase owner must use the same birth references as the coupled
    /// owner, even when the evaporating cell activates at an interval endpoint.
    pub(super) fn thermal_interval(
        &mut self,
        initial: &DepositionSimulation<'_>,
        next: &DepositionSimulation<'_>,
        heat: &EnergyAudit,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        for (index, event) in next.schedule[initial.cursor..next.cursor].iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if let Event::Birth(cell) = &event.event {
                self.birth(cell.id, cell.mass_kg, cell.specific_h_j_kg)?;
            }
        }
        self.heat(heat.external_heat_j)
    }

    pub(super) fn phase(
        &mut self,
        initial: &DepositionSimulation<'_>,
        id: CellId,
        step: SurfaceStep,
    ) -> Result<(), DepositionError> {
        let baseline_h = initial
            .cells
            .get(&id)
            .map(|cell| cell.specific_h_j_kg)
            .or_else(|| self.births.get(&id).map(|&(_, h)| h))
            .ok_or(DepositionError::Arithmetic)?;
        self.mass_transfer.add(step.transferred_mass_kg)?;
        self.latent.add(step.latent_energy_j)?;
        self.vapor.add(step.vapor_energy_j)?;
        self.relative_vapor.add(finite(
            step.transferred_mass_kg * finite(step.liquid.specific_enthalpy_j_kg - baseline_h)?
                + step.latent_energy_j,
        )?)?;
        self.recession.add(step.recession_m)?;
        self.phase_steps = self
            .phase_steps
            .checked_add(1)
            .ok_or(DepositionError::Arithmetic)?;
        self.last_step = Some(step);
        Ok(())
    }

    /// Reconstruct physical changes relative to each original/born cell.
    /// An explicit caller may use the original roundoff-scaled mass allowance;
    /// an implicit caller supplies its unchanged absolute kilogram tolerance.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn finish(
        self,
        initial: &DepositionSimulation<'_>,
        next: &DepositionSimulation<'_>,
        substeps: usize,
        events: usize,
        energy_tolerance_j: f64,
        mass_tolerance_kg: Option<f64>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<CoupledEvaporationAudit, ThermalEvaporationError> {
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
                None => *self.births.get(&id).ok_or(DepositionError::Arithmetic)?,
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
        let energy_residual_j = energy_balance.value()?;
        let incoming_mass_kg = self.incoming_mass.value()?;
        let mass_tolerance = if let Some(tolerance) = mass_tolerance_kg {
            tolerance
        } else {
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
            finite(128.0 * f64::EPSILON * scale * terms)?
        };
        if mass_balance_residual_kg.abs() > mass_tolerance {
            return Err(ThermalEvaporationError::MassBalance {
                residual_kg: mass_balance_residual_kg,
                tolerance_kg: mass_tolerance,
            });
        }
        if energy_residual_j.abs() > energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(energy_residual_j).into());
        }
        finite(next.mass_kg())?;
        finite(next.energy_j())?;
        Ok(CoupledEvaporationAudit {
            incoming_mass_kg,
            incoming_energy_j: self.incoming_energy.value()?,
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
