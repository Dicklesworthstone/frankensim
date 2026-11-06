//! Accepted interval and original birth-relative storage on the native graph.
use super::*;
use fs_thermal::enthalpy::{CellId, EnergyAudit, Event, TimedEvent};
use std::collections::BTreeMap;

#[derive(Default)]
pub(super) struct Account {
    events: usize,
    births: BTreeMap<CellId, (f64, f64)>,
    incoming_mass: Sum,
    incoming_energy: Sum,
    external_thermal: Sum,
    external_gas: Sum,
    wall_storage: Sum,
    film_storage: Sum,
    gas_storage: Sum,
    evaporated: Sum,
    condensed: Sum,
    transferred: Sum,
    latent: Sum,
    vapor_energy: Sum,
    sensible: Sum,
    phase_wall_heat: Sum,
    energy_residual: Sum,
    absolute_energy: Sum,
    absolute_mass: Sum,
}
impl Account {
    pub(super) fn add_events(
        &mut self,
        pending: &[TimedEvent],
        audit: &EnergyAudit,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitGasError> {
        if audit.substeps != 0 || audit.events > pending.len() {
            return Err(ImplicitGasError::Invalid(
                "invalid native wall-film event receipt",
            ));
        }
        self.events = count(self.events, audit.events)?;
        add(&mut self.incoming_mass, audit.incoming_mass_kg)?;
        add(&mut self.incoming_energy, audit.incoming_energy_j)?;
        add(&mut self.external_thermal, audit.external_heat_j)?;
        add(&mut self.energy_residual, audit.residual_j)?;
        add(&mut self.absolute_energy, audit.residual_j.abs())?;
        for event in pending.iter().take(audit.events) {
            checkpoint()?;
            if let Event::Birth(cell) = &event.event {
                if self
                    .births
                    .insert(cell.id, (cell.mass_kg, cell.specific_h_j_kg))
                    .is_some()
                {
                    return Err(ImplicitGasError::Invalid("duplicate native thermal birth"));
                }
            }
        }
        Ok(())
    }
    pub(super) fn add(&mut self, step: &DepositionWallFilmReport) -> Result<(), ImplicitGasError> {
        let phase = &step.phase;
        add(&mut self.external_thermal, step.external_thermal_heat_j)?;
        add(&mut self.external_gas, step.external_gas_heat_j)?;
        add(&mut self.wall_storage, step.thermal.stored_energy_change_j)?;
        add(&mut self.film_storage, phase.film_enthalpy_change_j)?;
        add(&mut self.gas_storage, phase.gas_internal_energy_change_j)?;
        add(&mut self.evaporated, phase.transferred_mass_kg.max(0.0))?;
        add(&mut self.condensed, (-phase.transferred_mass_kg).max(0.0))?;
        add(&mut self.transferred, phase.transferred_mass_kg)?;
        add(&mut self.latent, phase.latent_energy_j)?;
        add(&mut self.vapor_energy, phase.energy_into_gas_j)?;
        add(&mut self.sensible, step.sensible_heat_into_gas_j)?;
        add(
            &mut self.phase_wall_heat,
            step.wall_heat_into_phase_system_j,
        )?;
        add(&mut self.energy_residual, step.coupled_energy_residual_j)?;
        add(
            &mut self.absolute_energy,
            step.sum_absolute_energy_residual_j,
        )?;
        add(&mut self.absolute_mass, step.sum_absolute_mass_residual_kg)?;
        Ok(())
    }
    pub(super) fn gate(&self, window: WallFilmWindowSolve) -> Result<(), ImplicitGasError> {
        let energy = finite(self.absolute_energy.value())?;
        if energy > window.energy_tolerance_j {
            return Err(ImplicitGasError::CoupledEnergyBalance {
                residual_j: energy,
                tolerance_j: window.energy_tolerance_j,
            });
        }
        let mass = finite(self.absolute_mass.value())?;
        if mass > window.mass_tolerance_kg {
            return Err(ImplicitGasError::GasMassBalance {
                residual_kg: mass,
                tolerance_kg: window.mass_tolerance_kg,
            });
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn finish(
        mut self,
        initial: &DepositionSimulation<'_>,
        next: &DepositionSimulation<'_>,
        binding: DepositionWallFilmBinding<'_>,
        gas: &ClosedGasParcel,
        next_gas: &ClosedGasParcel,
        film: GasWallFilmInventory,
        external_gas_heat_j: f64,
        window: WallFilmWindowSolve,
        accepted: usize,
        attempted: usize,
        rejected: usize,
        meter: &WorkMeter,
        last_step: Option<DepositionWallFilmReport>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<DepositionWallFilmTimelineReport, ImplicitGasError> {
        checkpoint()?;
        if next.cells().len() != count(initial.cells().len(), self.births.len())?
            || attempted != count(accepted, rejected)?
        {
            return Err(ImplicitGasError::Invalid(
                "inconsistent native wall-film timeline inventory",
            ));
        }
        // Subtract each original caloric datum BEFORE multiplying and summing.
        // A newly born neighbor's formation/latent reference is incoming storage,
        // not heat supplied to the old wall or to the independent condensate.
        let mut stored = Sum::default();
        for (&id, cell) in next.cells() {
            checkpoint()?;
            let (mass, h0) = initial
                .cells()
                .get(&id)
                .map(|old| (old.mass_kg, old.specific_h_j_kg))
                .or_else(|| self.births.get(&id).copied())
                .ok_or(ImplicitGasError::Invalid(
                    "thermal endpoint has an unaccounted cell",
                ))?;
            if cell.mass_kg != mass {
                return Err(ImplicitGasError::Invalid(
                    "wall-film exchange changed native wall mass",
                ));
            }
            add(&mut stored, finite(mass * (cell.specific_h_j_kg - h0))?)?;
        }
        let telescoping = -difference(self.wall_storage, stored.value())?;
        let apportionment = difference(self.external_gas, external_gas_heat_j)?;
        add(&mut self.energy_residual, telescoping)?;
        add(&mut self.energy_residual, apportionment)?;
        add(&mut self.absolute_energy, telescoping.abs())?;
        add(&mut self.absolute_energy, apportionment.abs())?;
        self.gate(window)?;
        let energy = finite(self.energy_residual.value())?;
        if energy.abs() > window.energy_tolerance_j {
            return Err(ImplicitGasError::CoupledEnergyBalance {
                residual_j: energy,
                tolerance_j: window.energy_tolerance_j,
            });
        }
        let film_change = finite(film.mass_kg() - binding.film.mass_kg())?;
        let total_change = finite(next_gas.mass_kg() - gas.mass_kg())?;
        let species = &binding.boundary.binding.species;
        let volatile_change =
            finite(next_gas.species_mass_kg(species)? - gas.species_mass_kg(species)?)?;
        let mass_residual = finite(total_change + film_change)?;
        let volatile_residual = finite(volatile_change + film_change)?;
        let mut carriers = 0.0_f64;
        for source in gas.model().thermodynamics().components() {
            checkpoint()?;
            if source.species() != species {
                carriers = carriers.max(
                    finite(
                        next_gas.species_mass_kg(source.species())?
                            - gas.species_mass_kg(source.species())?,
                    )?
                    .abs(),
                );
            }
        }
        let maximum_mass = mass_residual
            .abs()
            .max(volatile_residual.abs())
            .max(carriers)
            .max(difference(self.transferred, volatile_change)?.abs())
            .max(difference(self.transferred, -film_change)?.abs())
            .max(self.absolute_mass.value());
        if maximum_mass > window.mass_tolerance_kg {
            return Err(ImplicitGasError::GasMassBalance {
                residual_kg: maximum_mass,
                tolerance_kg: window.mass_tolerance_kg,
            });
        }
        checkpoint()?;
        Ok(DepositionWallFilmTimelineReport {
            initial_time_s: initial.time_s(),
            final_time_s: next.time_s(),
            accepted_intervals: accepted,
            attempted_intervals: attempted,
            rejected_intervals: rejected,
            events: self.events,
            incoming_mass_kg: self.incoming_mass.value(),
            incoming_energy_j: self.incoming_energy.value(),
            external_thermal_heat_j: self.external_thermal.value(),
            external_gas_heat_j,
            wall_storage_change_j: stored.value(),
            film_enthalpy_change_j: self.film_storage.value(),
            gas_internal_energy_change_j: self.gas_storage.value(),
            evaporated_mass_kg: self.evaporated.value(),
            condensed_mass_kg: self.condensed.value(),
            latent_energy_j: self.latent.value(),
            vapor_energy_into_gas_j: self.vapor_energy.value(),
            sensible_heat_into_gas_j: self.sensible.value(),
            wall_heat_into_phase_system_j: self.phase_wall_heat.value(),
            gas_heat_apportionment_residual_j: apportionment,
            thermal_telescoping_residual_j: telescoping,
            combined_energy_residual_j: energy,
            sum_absolute_energy_residual_j: self.absolute_energy.value(),
            mass_residual_kg: mass_residual,
            volatile_mass_residual_kg: volatile_residual,
            maximum_carrier_mass_residual_kg: carriers,
            sum_absolute_mass_residual_kg: self.absolute_mass.value(),
            newton_iterations: meter.newton_iterations(),
            residual_evaluations: meter.residual_evaluations(),
            gas_iterations: meter.gas_iterations(),
            charged_gas_iterations: meter.charged_gas_iterations(),
            last_step,
        })
    }
}
