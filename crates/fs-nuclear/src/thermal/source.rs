//! The shared transactional heat-owner seam for admitted nuclear source laws.

use super::{FissionCellSource, FissionThermalError};
use crate::irradiation::GroupFlux;
use crate::{EnergyChannels, FissionBudget, FissionError, FissionStep};
use fs_thermal::enthalpy::{CellId, MaterialId};

pub(super) mod sealed {
    pub trait Sealed {}
}

/// A cell-bound nuclear source using the shared thermal transaction.
///
/// This sealed interface keeps prescribed-event, neutron-irradiation and
/// radioactive-production laws on one binding, scheduling, heat-credit and
/// rollback implementation. Public source types carry their own physical
/// inventories and admitted budgets.
pub trait NuclearHeatSource: sealed::Sealed + Clone + core::fmt::Debug {
    /// Explicit work and numerical-accounting budget for this source law.
    type Budget: Copy + core::fmt::Debug;
    /// The source law's accepted interval report.
    type Step: Clone + core::fmt::Debug;

    /// Original thermal cell receiving only locally deposited energy.
    fn cell_id(&self) -> CellId;
    /// Original caloric material identity bound to the cell.
    fn material_id(&self) -> MaterialId;
    /// A spatial nuclear model can bind its number-density volume to the
    /// original thermal geometry. Local sources without geometry return None.
    fn expected_volume_m3(&self) -> Option<f64> {
        None
    }
    /// Advance-wide allowance for noncancelling nuclear and coupling defects.
    fn energy_tolerance_j(budget: Self::Budget) -> f64;
    /// Optional advance-wide allowance for noncancelling population defects.
    /// Existing source laws retain their own accounting when this is absent.
    fn nuclei_tolerance(_budget: Self::Budget) -> Option<f64> {
        None
    }
    /// Next authored discontinuity strictly after the accepted clock.
    fn next_transition(&self, time_s: f64) -> Option<f64>;
    /// Admit immutable source data, the initial environment and work limits.
    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError>;
    /// Refuse temperatures outside the declared nuclear-data domain.
    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError>;
    /// Prepare an independent source candidate over the final thermal interval.
    fn trial_source(
        &self,
        time_s: f64,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError>;
    /// Prepare all cell sources over the same final thermal interval. The
    /// default keeps independent local sources independent; a spatial source
    /// overrides this once for its shared neutron field and inventories.
    fn trial_sources(
        sources: &[&Self],
        time_s: f64,
        duration_s: f64,
        temperatures_k: &[f64],
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<(Self, Self::Step)>, FissionThermalError> {
        if sources.len() != temperatures_k.len() {
            return Err(FissionThermalError::Invalid(
                "one material temperature per nuclear source required",
            ));
        }
        sources
            .iter()
            .zip(temperatures_k)
            .map(|(source, &temperature)| {
                checkpoint()?;
                source.trial_source(time_s, duration_s, temperature, budget, checkpoint)
            })
            .collect()
    }
    /// Actual source release over an accepted interval.
    fn released(step: &Self::Step) -> EnergyChannels;
    /// Noncancelling numerical energy-account defects over that interval.
    fn energy_roundoff_j(step: &Self::Step) -> f64;
    /// Noncancelling population and production-account defects over the interval.
    fn population_roundoff_nuclei(_step: &Self::Step) -> f64 {
        0.0
    }
    /// Integrated radioactive transformations, including newly produced daughters.
    fn decays(_step: &Self::Step) -> f64 {
        0.0
    }
    /// Independently supplied nuclei integrated over the interval.
    fn produced_nuclei(_step: &Self::Step) -> f64 {
        0.0
    }
    /// Accepted fission events over that interval.
    fn fissions(step: &Self::Step) -> f64;
    /// Accepted competing captures over that interval.
    fn captures(step: &Self::Step) -> f64;
    /// The source's lifetime release, including all explicit initial reservoirs.
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError>;
}

/// An admitted irradiation law whose reference spectrum can be scaled by the
/// same lumped neutron kinetics owner. This is sealed through `NuclearHeatSource`.
/// Both independent parents and full product networks retain their own reaction
/// propagation; the shared bridge owns neutron normalization and thermal time.
pub trait FluxDrivenHeatSource: NuclearHeatSource {
    /// Authored absolute group flux at neutron population amplitude one.
    fn reference_flux(&self, time_s: f64) -> Result<&GroupFlux, FissionThermalError>;
    /// Maximum groups admitted by this source's numerical budget.
    fn max_neutron_groups(budget: Self::Budget) -> usize;
    /// Advance against the supplied absolute mean spectrum over this interval.
    fn trial_flux(
        &self,
        flux: &GroupFlux,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError>;
    /// Advance on the shared absolute clock when another source history, such
    /// as supplied nuclides, must be selected independently of neutron flux.
    /// Existing flux-only laws reuse their original trial unchanged. Coupling
    /// owners use this timed form so scheduled sources never guess a start time.
    fn trial_flux_at(
        &self,
        time_s: f64,
        flux: &GroupFlux,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(FissionThermalError::Invalid(
                "invalid absolute flux-trial time",
            ));
        }
        self.trial_flux(flux, duration_s, temperature_k, budget, checkpoint)
    }
    /// Caller-supplied composition worth at the beginning of the interval.
    /// The parent-only model has no inferred composition-to-reactivity law.
    fn composition_reactivity(
        &self,
        _budget: Self::Budget,
        _checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<f64, FissionThermalError> {
        Ok(0.0)
    }
}

impl sealed::Sealed for FissionCellSource {}
impl NuclearHeatSource for FissionCellSource {
    type Budget = FissionBudget;
    type Step = FissionStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.energy_tolerance_j
    }
    fn next_transition(&self, time_s: f64) -> Option<f64> {
        self.rates.next_after(time_s)
    }
    fn admit(&self, time_s: f64, _: f64, budget: Self::Budget) -> Result<(), FissionThermalError> {
        self.state
            .trial(0.0, self.rates.rate_at(time_s)?, budget, || Ok(()))?;
        Ok(())
    }
    fn validate_temperature(&self, _: f64) -> Result<(), FissionThermalError> {
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
        let (state, report) =
            self.state
                .trial(duration_s, self.rates.rate_at(time_s)?, budget, checkpoint)?;
        Ok((
            Self {
                cell_id: self.cell_id,
                material_id: self.material_id,
                state,
                rates: self.rates.clone(),
            },
            report,
        ))
    }
    fn released(step: &Self::Step) -> EnergyChannels {
        step.released
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.energy_roundoff_j
    }
    fn fissions(step: &Self::Step) -> f64 {
        step.fissions
    }
    fn captures(_: &Self::Step) -> f64 {
        0.0
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.state.cumulative_released())
    }
}
