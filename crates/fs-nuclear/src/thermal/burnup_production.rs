//! Supplied nuclides and their subsequent reactions on the shared thermal clock.
//!
//! The original burnup source keeps its flux, material binding and composition
//! worth. Its inventory integrates production, decay, capture and fission over
//! the same accepted interval, with actual ages of newly supplied nuclei. Only
//! reaction energy enters heat; feed mass and incoming enthalpy are not inferred.
//!
//! This source also participates in the existing neutron-kinetics bridge. That
//! bridge uses its declared mean-flux approximation when neutron amplitude
//! changes, so coupled birth timing and heat still require timestep refinement.
//! Production schedules use the explicit absolute-time hook; an untimed flux
//! trial cannot select a changing production history.

use std::sync::Arc;

use super::source::{FluxDrivenHeatSource, NuclearHeatSource, sealed};
use super::{BurnupCellSource, BurnupHeatingBudget, FissionThermalError, NuclearHeating, Sum};
use crate::burnup::production::{BurnupProductionSchedule, BurnupProductionStep};
use crate::burnup::{BurnupError, BurnupInventory};
use crate::irradiation::GroupFlux;
use crate::{EnergyChannels, FissionError};
use fs_thermal::enthalpy::{CellId, MaterialId};

/// One original reaction-network heat source with an admitted production history.
/// Initial populations describe the binding time; previous production is never
/// replayed. The immutable schedule is shared across clones and must belong to
/// the exact inventory model. The original inventory owns preparation caching.
#[derive(Clone, Debug)]
pub struct BurnupProductionCellSource {
    source: BurnupCellSource,
    production: Arc<BurnupProductionSchedule>,
}

impl BurnupProductionCellSource {
    /// Bind the original source to an immutable schedule of the same physical model.
    pub fn new(
        source: BurnupCellSource,
        production: Arc<BurnupProductionSchedule>,
    ) -> Result<Self, BurnupError> {
        if !Arc::ptr_eq(source.inventory.model(), production.model()) {
            return Err(BurnupError::Invalid(
                "production schedule belongs to a different burnup model",
            ));
        }
        Ok(Self { source, production })
    }

    /// Original material, neutron-flux and composition-worth binding.
    pub fn source(&self) -> &BurnupCellSource {
        &self.source
    }

    /// Current accepted populations and their production/reaction accounts.
    pub fn inventory(&self) -> &BurnupInventory {
        &self.source.inventory
    }

    /// Immutable production history in the original nuclide order.
    pub fn production(&self) -> &Arc<BurnupProductionSchedule> {
        &self.production
    }

    fn trial_rates(
        &self,
        rates_nuclei_s: &[f64],
        flux: &GroupFlux,
        duration_s: f64,
        temperature_k: f64,
        budget: BurnupHeatingBudget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, BurnupProductionStep), FissionThermalError> {
        self.source.worth(budget, checkpoint)?;
        let (inventory, mut step) = self.source.inventory.trial_with_production(
            duration_s,
            flux,
            temperature_k,
            rates_nuclei_s,
            budget.source,
            &mut *checkpoint,
        )?;
        // The original burnup report already charges reaction-count reduction.
        // Charge the additional scalar production observation independently.
        let mut produced = Sum::default();
        for &count in &step.produced_nuclei {
            produced.add(count)?;
        }
        let mut defects = Sum::default();
        defects.add(step.burnup.population_roundoff_nuclei)?;
        defects.add(produced.export_roundoff()?)?;
        step.burnup.population_roundoff_nuclei = defects.nonnegative_upper_bound()?;
        if step.burnup.population_roundoff_nuclei > budget.source.nuclei_tolerance {
            return Err(FissionThermalError::NucleiBudget {
                defect_nuclei: step.burnup.population_roundoff_nuclei,
                tolerance_nuclei: budget.source.nuclei_tolerance,
            });
        }
        let source = BurnupCellSource {
            cell_id: self.source.cell_id,
            material_id: self.source.material_id,
            inventory,
            fluxes: self.source.fluxes.clone(),
            composition_feedback: self.source.composition_feedback.clone(),
        };
        // Incoming nuclei may change calibrated worth or leave its applicability
        // domain. Check the candidate before either heat or history can publish.
        source.worth(budget, checkpoint)?;
        Ok((
            Self {
                source,
                production: Arc::clone(&self.production),
            },
            step,
        ))
    }
}

/// Prescribed production and all subsequent reactions on the original heat owner.
pub type BurnupProductionHeating = NuclearHeating<BurnupProductionCellSource>;

impl sealed::Sealed for BurnupProductionCellSource {}
impl NuclearHeatSource for BurnupProductionCellSource {
    type Budget = BurnupHeatingBudget;
    type Step = BurnupProductionStep;

    fn cell_id(&self) -> CellId {
        self.source.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.source.material_id
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        BurnupCellSource::energy_tolerance_j(budget)
    }
    fn nuclei_tolerance(budget: Self::Budget) -> Option<f64> {
        BurnupCellSource::nuclei_tolerance(budget)
    }
    fn next_transition(&self, time_s: f64) -> Option<f64> {
        [
            self.source.next_transition(time_s),
            self.production.next_after(time_s),
        ]
        .into_iter()
        .flatten()
        .min_by(f64::total_cmp)
    }
    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError> {
        self.trial_flux_at(
            time_s,
            self.reference_flux(time_s)?,
            0.0,
            temperature_k,
            budget,
            &mut || Ok(()),
        )?;
        Ok(())
    }
    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        self.source.validate_temperature(temperature_k)
    }
    fn trial_source(
        &self,
        time_s: f64,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        self.trial_flux_at(
            time_s,
            self.reference_flux(time_s)?,
            duration_s,
            temperature_k,
            budget,
            checkpoint,
        )
    }
    fn released(step: &Self::Step) -> EnergyChannels {
        BurnupCellSource::released(&step.burnup)
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        BurnupCellSource::energy_roundoff_j(&step.burnup)
    }
    fn population_roundoff_nuclei(step: &Self::Step) -> f64 {
        BurnupCellSource::population_roundoff_nuclei(&step.burnup)
    }
    fn decays(step: &Self::Step) -> f64 {
        BurnupCellSource::decays(&step.burnup)
    }
    fn produced_nuclei(step: &Self::Step) -> f64 {
        crate::sum(&step.produced_nuclei)
    }
    fn fissions(step: &Self::Step) -> f64 {
        BurnupCellSource::fissions(&step.burnup)
    }
    fn captures(step: &Self::Step) -> f64 {
        BurnupCellSource::captures(&step.burnup)
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        self.source.cumulative_released()
    }
}

impl FluxDrivenHeatSource for BurnupProductionCellSource {
    fn reference_flux(&self, time_s: f64) -> Result<&GroupFlux, FissionThermalError> {
        self.source.reference_flux(time_s)
    }
    fn max_neutron_groups(budget: Self::Budget) -> usize {
        BurnupCellSource::max_neutron_groups(budget)
    }
    fn composition_reactivity(
        &self,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<f64, FissionThermalError> {
        self.source.worth(budget, checkpoint)
    }
    fn trial_flux(
        &self,
        flux: &GroupFlux,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        if !self.production.changes().is_empty() {
            return Err(FissionThermalError::Invalid(
                "scheduled production requires an absolute-time flux trial",
            ));
        }
        self.trial_rates(
            self.production.rates_at(0.0)?,
            flux,
            duration_s,
            temperature_k,
            budget,
            checkpoint,
        )
    }
    fn trial_flux_at(
        &self,
        time_s: f64,
        flux: &GroupFlux,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        let rates = self.production.rates_at(time_s)?;
        if self
            .production
            .next_after(time_s)
            .is_some_and(|next| duration_s > next - time_s)
        {
            return Err(FissionThermalError::Invalid(
                "production trial crosses a scheduled rate change",
            ));
        }
        self.trial_rates(rates, flux, duration_s, temperature_k, budget, checkpoint)
    }
}
