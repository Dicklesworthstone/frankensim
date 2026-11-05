//! Capture activation and its delayed heat on the shared thermal transaction.
//!
//! The original capture kernel retains the exponentially depleted parent and
//! the actual ages of newly formed radioactive products. Prompt irradiation
//! energy and product-decay energy enter the thermal owner once, through their
//! combined deposited channel. The original neutron grid, flux schedule,
//! temperature applicability and fixed thermal material mass remain authoritative.
//!
//! `KineticCellSource<CaptureActivationCellSource>` uses the existing neutron
//! kinetics bridge and its interval-mean absolute flux. The activation kernel
//! retains birth timing for that frozen flux; changing neutron amplitude still
//! requires timestep refinement because a mean flux does not preserve the exact
//! within-interval production history. No composition worth is inferred.

use fs_thermal::enthalpy::{CellId, MaterialId};

use super::source::{FluxDrivenHeatSource, NuclearHeatSource, sealed};
use super::{FissionThermalError, NeutronFluxSchedule, NuclearHeating, Sum, finite};
use crate::decay::activation::{
    CaptureActivationBudget, CaptureActivationError, CaptureActivationState, CaptureActivationStep,
};
use crate::irradiation::GroupFlux;
use crate::{EnergyChannels, FissionError};

/// One explicit irradiation parent and its capture-product inventory in a cell.
///
/// Capture identity and isotope data belong to the original activation state.
/// Flux changes use the existing right-continuous neutron-grid schedule. Initial
/// parent and product populations describe the binding time, including whatever
/// prior history the caller explicitly retained; no irradiation is replayed.
#[derive(Clone, Debug)]
pub struct CaptureActivationCellSource {
    /// Existing active thermal receiver.
    pub cell_id: CellId,
    /// Original thermal material identity.
    pub material_id: MaterialId,
    /// Accepted irradiation and radioactive-product histories.
    pub state: CaptureActivationState,
    /// Incident multigroup scalar flux on the parent's original energy grid,
    /// or the reference at neutron amplitude one when wrapped in `KineticCellSource`.
    pub fluxes: NeutronFluxSchedule,
}

/// The shared nuclear heat owner specialized for radiative capture activation.
pub type CaptureActivationHeating = NuclearHeating<CaptureActivationCellSource>;

impl sealed::Sealed for CaptureActivationCellSource {}

impl NuclearHeatSource for CaptureActivationCellSource {
    type Budget = CaptureActivationBudget;
    type Step = CaptureActivationStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }

    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.decay.energy_tolerance_j
    }

    fn nuclei_tolerance(budget: Self::Budget) -> Option<f64> {
        Some(budget.decay.nuclei_tolerance)
    }

    fn next_transition(&self, time_s: f64) -> Option<f64> {
        let changes = self.fluxes.changes();
        changes
            .get(changes.partition_point(|change| change.time_s <= time_s))
            .map(|change| change.time_s)
    }

    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError> {
        self.state.trial(
            0.0,
            self.fluxes.flux_at(time_s)?,
            temperature_k,
            budget,
            || Ok(()),
        )?;
        Ok(())
    }

    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        let model = self.state.parent().model();
        if !temperature_k.is_finite()
            || temperature_k < model.temperature_min_k()
            || temperature_k > model.temperature_max_k()
        {
            return Err(
                CaptureActivationError::Irradiation(FissionError::InvalidInput(
                    "thermal temperature is outside nuclear-data applicability",
                ))
                .into(),
            );
        }
        Ok(())
    }

    fn trial_source(
        &self,
        time_s: f64,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        self.trial_flux(
            self.fluxes.flux_at(time_s)?,
            duration_s,
            temperature_k,
            budget,
            checkpoint,
        )
    }

    fn released(step: &Self::Step) -> EnergyChannels {
        step.released
    }

    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.energy_roundoff_j
    }

    fn population_roundoff_nuclei(step: &Self::Step) -> f64 {
        step.population_roundoff_nuclei
    }

    fn decays(step: &Self::Step) -> f64 {
        crate::sum(&step.products.decay.decays)
    }

    fn produced_nuclei(step: &Self::Step) -> f64 {
        crate::sum(&step.products.produced_nuclei)
    }

    fn fissions(step: &Self::Step) -> f64 {
        step.irradiation.fissions
    }

    fn captures(step: &Self::Step) -> f64 {
        step.irradiation.captures
    }

    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        let mut channels = [Sum::default(); 3];
        for energy in [
            self.state.parent().cumulative_released(),
            self.state.products().cumulative_released(),
        ] {
            channels[0].add(energy.deposited_j)?;
            channels[1].add(energy.radiation_j)?;
            channels[2].add(energy.escaped_j)?;
        }
        let released = EnergyChannels {
            deposited_j: channels[0].value()?,
            radiation_j: channels[1].value()?,
            escaped_j: channels[2].value()?,
        };
        finite(released.total_j())?;
        Ok(released)
    }
}

impl FluxDrivenHeatSource for CaptureActivationCellSource {
    fn reference_flux(&self, time_s: f64) -> Result<&GroupFlux, FissionThermalError> {
        self.fluxes.flux_at(time_s)
    }

    fn max_neutron_groups(budget: Self::Budget) -> usize {
        budget.irradiation.max_neutron_groups
    }

    fn trial_flux(
        &self,
        flux: &GroupFlux,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        let mut checkpoint_error = None;
        let trial = self
            .state
            .trial(duration_s, flux, temperature_k, budget, || {
                checkpoint().map_err(|error| {
                    checkpoint_error = Some(error);
                    CaptureActivationError::Cancelled
                })
            });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let (state, mut step) = trial?;
        // Per-node counts remain in the original report. Charge any low part
        // lost when the common thermal audit exports one scalar count per cell.
        let mut defects = Sum::default();
        defects.add(step.population_roundoff_nuclei)?;
        for counts in [&step.products.decay.decays, &step.products.produced_nuclei] {
            let mut total = Sum::default();
            for &count in counts {
                total.add(count)?;
            }
            defects.add(total.export_roundoff()?)?;
        }
        step.population_roundoff_nuclei = defects.nonnegative_upper_bound()?;
        if step.population_roundoff_nuclei > budget.decay.nuclei_tolerance {
            return Err(FissionThermalError::NucleiBudget {
                defect_nuclei: step.population_roundoff_nuclei,
                tolerance_nuclei: budget.decay.nuclei_tolerance,
            });
        }
        Ok((
            Self {
                cell_id: self.cell_id,
                material_id: self.material_id,
                state,
                fluxes: self.fluxes.clone(),
            },
            step,
        ))
    }
}

#[cfg(test)]
mod tests;
