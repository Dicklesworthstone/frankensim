//! Scheduled production into radioactive inventories on the shared heat owner.
//!
//! Production rates are caller-authored expected nuclei per second. The original
//! production propagator integrates both newborn survivors and their decays over
//! the final thermal interval. It does not infer fission yields, capture products,
//! a neutron history or an energy carried into the cell by those incoming nuclei.
//! Only the supplied branch energy deposited during radioactive transformations
//! enters the existing thermal source account.

use std::sync::Arc;

use fs_la::decay::DecayLimits;
use fs_thermal::enthalpy::{CellId, MaterialId};

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NuclearHeating, Sum};
use crate::decay::production::{
    DecayProductionSchedule, DecayProductionStep, PreparedDecayProduction,
};
use crate::decay::{DecayBudget, DecayError, DecayInventory};
use crate::{EnergyChannels, FissionError};

#[derive(Clone, Debug)]
struct PreparedInterval {
    plan: Arc<PreparedDecayProduction>,
    limits: DecayLimits,
}

/// Explicit radioactive inventory and production history bound to one thermal cell.
///
/// The immutable schedule must belong to the exact inventory model. Initial
/// populations describe the binding time; earlier production is never replayed.
/// Preparation is cached for repeated durations and propagation limits, while
/// changed production rates are applied to the same immutable map. Thermal mass,
/// calorics and geometry remain with the original thermal material owner.
#[derive(Clone, Debug)]
pub struct DecayProductionCellSource {
    cell_id: CellId,
    material_id: MaterialId,
    inventory: DecayInventory,
    schedule: Arc<DecayProductionSchedule>,
    prepared: Option<PreparedInterval>,
}

impl DecayProductionCellSource {
    /// Bind one current inventory to an explicitly authored model-specific schedule.
    pub fn new(
        cell_id: CellId,
        material_id: MaterialId,
        inventory: DecayInventory,
        schedule: Arc<DecayProductionSchedule>,
    ) -> Result<Self, DecayError> {
        if !Arc::ptr_eq(inventory.model(), schedule.model()) {
            return Err(DecayError::Invalid(
                "production schedule belongs to a different radioactive model",
            ));
        }
        Ok(Self {
            cell_id,
            material_id,
            inventory,
            schedule,
            prepared: None,
        })
    }

    /// Original active thermal receiver.
    pub fn cell_id(&self) -> CellId {
        self.cell_id
    }

    /// Original thermal material identity.
    pub fn material_id(&self) -> MaterialId {
        self.material_id
    }

    /// Current accepted radioactive populations and their cumulative accounts.
    pub fn inventory(&self) -> &DecayInventory {
        &self.inventory
    }

    /// Immutable absolute-time production history, shared across clones.
    pub fn schedule(&self) -> &Arc<DecayProductionSchedule> {
        &self.schedule
    }
}

/// The existing nuclear/thermal transaction specialized for radioactive production.
pub type DecayProductionHeating = NuclearHeating<DecayProductionCellSource>;

impl sealed::Sealed for DecayProductionCellSource {}

impl NuclearHeatSource for DecayProductionCellSource {
    type Budget = DecayBudget;
    type Step = DecayProductionStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }

    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.energy_tolerance_j
    }

    fn nuclei_tolerance(budget: Self::Budget) -> Option<f64> {
        Some(budget.nuclei_tolerance)
    }

    fn next_transition(&self, time_s: f64) -> Option<f64> {
        self.schedule.next_after(time_s)
    }

    fn admit(&self, time_s: f64, _: f64, budget: Self::Budget) -> Result<(), FissionThermalError> {
        let prepared =
            self.inventory
                .model()
                .prepare_production(0.0, budget.propagation, || Ok(()))?;
        self.inventory.trial_production_prepared(
            &prepared,
            self.schedule.rates_at(time_s)?,
            budget,
            || Ok(()),
        )?;
        Ok(())
    }

    fn validate_temperature(&self, _: f64) -> Result<(), FissionThermalError> {
        // These admitted half-life/production laws are independent of temperature.
        // The original thermal law still admits its own material domain.
        Ok(())
    }

    fn trial_source(
        &self,
        time_s: f64,
        duration_s: f64,
        _: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        let rates = self.schedule.rates_at(time_s)?;
        let mut checkpoint_error = None;
        let trial = (|| -> Result<_, DecayError> {
            let mut check = || {
                checkpoint().map_err(|error| {
                    checkpoint_error = Some(error);
                    DecayError::Cancelled
                })
            };
            let prepared = match &self.prepared {
                Some(cached)
                    if cached.plan.duration_s() == duration_s
                        && cached.limits == budget.propagation =>
                {
                    cached.clone()
                }
                _ => PreparedInterval {
                    plan: Arc::new(self.inventory.model().prepare_production(
                        duration_s,
                        budget.propagation,
                        &mut check,
                    )?),
                    limits: budget.propagation,
                },
            };
            let (inventory, step) = self.inventory.trial_production_prepared(
                &prepared.plan,
                rates,
                budget,
                &mut check,
            )?;
            Ok((inventory, step, prepared))
        })();
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let (inventory, mut step, prepared) = trial?;
        // The shared audit exports one count per cell. Charge any low part lost
        // when per-nuclide counts are reduced, just as isotope energy reduction
        // is charged by the existing irradiation source.
        let mut defects = Sum::default();
        defects.add(step.decay.population_roundoff_nuclei)?;
        for counts in [&step.decay.decays, &step.produced_nuclei] {
            let mut total = Sum::default();
            for &count in counts {
                total.add(count)?;
            }
            defects.add(total.export_roundoff()?)?;
        }
        step.decay.population_roundoff_nuclei = defects.nonnegative_upper_bound()?;
        if step.decay.population_roundoff_nuclei > budget.nuclei_tolerance {
            return Err(FissionThermalError::NucleiBudget {
                defect_nuclei: step.decay.population_roundoff_nuclei,
                tolerance_nuclei: budget.nuclei_tolerance,
            });
        }
        Ok((
            Self {
                cell_id: self.cell_id,
                material_id: self.material_id,
                inventory,
                schedule: Arc::clone(&self.schedule),
                prepared: Some(prepared),
            },
            step,
        ))
    }

    fn released(step: &Self::Step) -> EnergyChannels {
        step.decay.released
    }

    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.decay.energy_roundoff_j
    }

    fn population_roundoff_nuclei(step: &Self::Step) -> f64 {
        step.decay.population_roundoff_nuclei
    }

    fn decays(step: &Self::Step) -> f64 {
        crate::sum(&step.decay.decays)
    }

    fn produced_nuclei(step: &Self::Step) -> f64 {
        crate::sum(&step.produced_nuclei)
    }

    fn fissions(_: &Self::Step) -> f64 {
        0.0
    }
    fn captures(_: &Self::Step) -> f64 {
        0.0
    }

    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.inventory.cumulative_released())
    }
}

#[cfg(test)]
mod tests;
