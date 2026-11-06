//! Product-forming neutron reactions, radioactive heat and composition feedback
//! on the original transactional thermal timeline.
//!
//! The model's explicit product inventories carry radioactive heat after
//! shutdown. Event energies must exclude future daughter-decay energy; this
//! source does not add effective fission decay-heat reservoirs a second time.
//! Material calorics and bulk mass retain their original thermal owner.

use super::source::{FluxDrivenHeatSource, NuclearHeatSource, sealed};
use super::{FissionThermalError, NeutronFluxSchedule, NuclearHeating, Sum};
use crate::burnup::{BurnupBudget, BurnupInventory, BurnupStep};
use crate::irradiation::GroupFlux;
use crate::kinetics::{CompositionFeedback, KineticsError};
use crate::{EnergyChannels, FissionError};
use fs_thermal::enthalpy::{CellId, MaterialId};

/// Reaction-network resources and optional composition-worth resources.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BurnupHeatingBudget {
    /// Original inventory evolution and energy admission. Its energy allowance
    /// also bounds the total nuclear/coupling defect across cells and substeps.
    pub source: BurnupBudget,
    /// Maximum independently supplied nuclide-worth fits for this cell.
    pub max_composition_feedback_terms: usize,
}

/// One complete physical reaction network sharing one spectrum and receiver.
#[derive(Clone, Debug)]
pub struct BurnupCellSource {
    /// Existing active thermal receiver.
    pub cell_id: CellId,
    /// Original caloric material identity.
    pub material_id: MaterialId,
    /// Fuel, capture products, fission products and radioactive daughters.
    pub inventory: BurnupInventory,
    /// Prescribed flux, or the reference flux at amplitude one when wrapped in
    /// `KineticCellSource<BurnupCellSource>`.
    pub fluxes: NeutronFluxSchedule,
    /// Explicit nuclide worths used by the neutron kinetics wrapper. With a
    /// prescribed flux they only enforce their declared applicability domains.
    pub composition_feedback: Vec<CompositionFeedback>,
}

/// Product and radioactive heating using the common nuclear/thermal owner.
pub type BurnupHeating = NuclearHeating<BurnupCellSource>;

impl BurnupCellSource {
    pub(super) fn worth(
        &self,
        budget: BurnupHeatingBudget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<f64, FissionThermalError> {
        if self.composition_feedback.len() > budget.max_composition_feedback_terms {
            return Err(KineticsError::Budget("composition-feedback terms").into());
        }
        let mut worth = Sum::default();
        for fit in &self.composition_feedback {
            checkpoint()?;
            let index = self
                .inventory
                .model()
                .nuclides()
                .iter()
                .position(|nuclide| nuclide.id == fit.nuclide())
                .ok_or(FissionThermalError::Invalid(
                    "composition feedback names an untracked nuclide",
                ))?;
            worth.add(fit.reactivity_shift(self.inventory.populations()[index])?)?;
        }
        worth.value()
    }
}

impl sealed::Sealed for BurnupCellSource {}
impl NuclearHeatSource for BurnupCellSource {
    type Budget = BurnupHeatingBudget;
    type Step = BurnupStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.source.energy_tolerance_j
    }
    fn nuclei_tolerance(budget: Self::Budget) -> Option<f64> {
        Some(budget.source.nuclei_tolerance)
    }
    fn next_transition(&self, time_s: f64) -> Option<f64> {
        self.fluxes.next_after(time_s)
    }

    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError> {
        // The same physical owner admits rates, data identity and all numerical
        // budgets even for a zero-duration/idle thermal binding.
        self.trial_flux(
            self.fluxes.flux_at(time_s)?,
            0.0,
            temperature_k,
            budget,
            &mut || Ok(()),
        )?;
        Ok(())
    }

    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        self.inventory.model().validate_temperature(temperature_k)?;
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
        step.decays
    }
    fn fissions(step: &Self::Step) -> f64 {
        step.fissions
    }
    fn captures(step: &Self::Step) -> f64 {
        step.captures
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.inventory.cumulative_released())
    }
}

impl FluxDrivenHeatSource for BurnupCellSource {
    fn reference_flux(&self, time_s: f64) -> Result<&GroupFlux, FissionThermalError> {
        self.fluxes.flux_at(time_s)
    }
    fn max_neutron_groups(budget: Self::Budget) -> usize {
        budget.source.limits.max_neutron_groups
    }
    fn composition_reactivity(
        &self,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<f64, FissionThermalError> {
        self.worth(budget, checkpoint)
    }
    fn trial_flux(
        &self,
        flux: &GroupFlux,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        self.worth(budget, checkpoint)?;
        let (inventory, report) = self.inventory.trial(
            duration_s,
            flux,
            temperature_k,
            budget.source,
            &mut *checkpoint,
        )?;
        let candidate = Self {
            cell_id: self.cell_id,
            material_id: self.material_id,
            inventory,
            fluxes: self.fluxes.clone(),
            composition_feedback: self.composition_feedback.clone(),
        };
        // Refuse a prospective product inventory outside the calibrated worth
        // domain BEFORE the common owner can publish its heat and source state.
        candidate.worth(budget, checkpoint)?;
        Ok((candidate, report))
    }
}
