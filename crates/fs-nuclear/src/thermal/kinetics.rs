//! Lumped neutron kinetics, material-temperature feedback and actual reactions.
//!
//! Each source owns one explicitly normalized, fixed-shape neutron region. The
//! native thermal cell temperature supplies its admitted feedback fits. The
//! point-kinetics population-time integral scales a reference multigroup flux;
//! that mean absolute flux advances the same finite parent-isotope inventories
//! used by the prescribed-irradiation owner. Only deposited joules enter heat.
//!
//! For fixed cross sections and spectral shape, integrated fluence determines
//! parent removal and competing reaction counts. The mean-flux replacement
//! does approximate the within-step time distribution of newly produced decay
//! heat. Temperature feedback is also frozen at the beginning of the native
//! thermal interval. Coupled histories therefore require timestep refinement;
//! analytic neutronics does not make this partitioned thermal solve exact.
//! No spatial neutron transport, neutron exchange between cells, burnup-driven
//! reactivity, or automatic fuel/isotope calibration is inferred.

use std::sync::Arc;

use super::source::{NuclearHeatSource, sealed};
use super::{
    FissionThermalError, IrradiationCellSource, IrradiationCellStep, IrradiationHeatingBudget,
    NuclearHeating, Sum, finite,
};
use crate::irradiation::GroupFlux;
use crate::kinetics::{
    KineticsBudget, KineticsError, KineticsStep, PointKineticsState, TemperatureFeedback,
};
use crate::{EnergyChannels, FissionError};
use fs_thermal::enthalpy::{CellId, MaterialId};

/// Constant controls over one authored interval, before thermal feedback.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KineticsControl {
    /// Reference/control reactivity rho, dimensionless and less than one.
    pub reactivity: f64,
    /// External neutron source in the kinetics population normalization per s.
    pub external_source_per_s: f64,
}

impl KineticsControl {
    fn validate(self) -> Result<(), FissionThermalError> {
        if !self.reactivity.is_finite()
            || self.reactivity >= 1.0
            || !self.external_source_per_s.is_finite()
            || self.external_source_per_s < 0.0
        {
            return Err(FissionThermalError::Invalid(
                "invalid neutron kinetics control",
            ));
        }
        Ok(())
    }
}

/// An absolute-time change in reference reactivity or external neutron source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KineticsControlChange {
    /// Strictly increasing positive absolute time [s].
    pub time_s: f64,
    /// Controls effective from this time, including at the transition itself.
    pub control: KineticsControl,
}

/// Immutable right-continuous kinetics controls, with binary-search lookup.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticsControlSchedule {
    initial: KineticsControl,
    changes: Arc<[KineticsControlChange]>,
}

impl KineticsControlSchedule {
    /// Validate every control and transition before retaining the schedule.
    pub fn new(
        initial: KineticsControl,
        changes: Vec<KineticsControlChange>,
    ) -> Result<Self, FissionThermalError> {
        initial.validate()?;
        let mut previous = 0.0;
        for change in &changes {
            change.control.validate()?;
            if !change.time_s.is_finite() || change.time_s <= previous {
                return Err(FissionThermalError::Invalid(
                    "kinetics controls need strictly increasing positive times",
                ));
            }
            previous = change.time_s;
        }
        Ok(Self {
            initial,
            changes: changes.into(),
        })
    }

    /// Controls at a finite nonnegative absolute time [s].
    pub fn control_at(&self, time_s: f64) -> Result<KineticsControl, FissionThermalError> {
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(FissionThermalError::Invalid(
                "invalid kinetics control time",
            ));
        }
        let index = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        Ok(if index == 0 {
            self.initial
        } else {
            self.changes[index - 1].control
        })
    }

    /// The complete authored schedule in increasing time order.
    pub fn changes(&self) -> &[KineticsControlChange] {
        &self.changes
    }

    fn next_after(&self, time_s: f64) -> Option<f64> {
        let index = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        self.changes.get(index).map(|change| change.time_s)
    }
}

/// Explicit work/accuracy limits for all owners in one coupled source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KineticHeatingBudget {
    /// Isotope/group limits and the common thermal nuclear-energy allowance.
    pub irradiation: IrradiationHeatingBudget,
    /// Population and integrated-population numerical estimates, separately
    /// budgeted from joule account closure and coupled-model approximation.
    pub kinetics: KineticsBudget,
    /// Maximum independently supplied temperature-feedback terms per cell.
    pub max_feedback_terms: usize,
}

/// One declared fixed-shape neutron region bound to one thermal receiver.
#[derive(Clone, Debug)]
pub struct KineticCellSource {
    /// Original isotope inventories and reference flux history. Here the
    /// authored flux is the absolute flux at neutron amplitude n=1; actual
    /// flux equals that reference times the explicitly normalized amplitude.
    pub irradiation: IrradiationCellSource,
    /// Current neutron and delayed-neutron precursor history for that region.
    pub kinetics: PointKineticsState,
    /// Independently authored control reactivity and external neutron source.
    pub controls: KineticsControlSchedule,
    /// Admitted feedback fits evaluated at this cell's original temperature.
    pub feedback: Vec<TemperatureFeedback>,
}

/// One accepted coupled interval, retaining both original physical reports.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticCellStep {
    /// Point-neutron solution and integrated neutron-population time.
    pub kinetics: KineticsStep,
    /// Actual finite-inventory fission, capture and energy accounts.
    pub irradiation: IrradiationCellStep,
    /// Reference/control reactivity plus the beginning-temperature feedback.
    pub applied_reactivity: f64,
    /// Integrated neutron-population time divided by the same interval length.
    pub mean_neutron_population: f64,
    /// Absolute roundoff in converting the integral to a mean and back [population*s].
    pub mean_population_roundoff_s: f64,
}

/// The shared transactional thermal owner with lumped neutron feedback.
pub type KineticHeating = NuclearHeating<KineticCellSource>;

impl KineticCellSource {
    fn reactivity_at(
        &self,
        control: KineticsControl,
        temperature_k: f64,
        budget: KineticHeatingBudget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<f64, FissionThermalError> {
        if self.feedback.len() > budget.max_feedback_terms {
            return Err(KineticsError::Budget("temperature-feedback terms").into());
        }
        let mut reactivity = Sum::default();
        reactivity.add(control.reactivity)?;
        for fit in &self.feedback {
            checkpoint()?;
            reactivity.add(fit.reactivity_shift(temperature_k)?)?;
        }
        reactivity.value()
    }

    fn mean_flux(
        &self,
        time_s: f64,
        amplitude: f64,
        budget: KineticHeatingBudget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<GroupFlux, FissionThermalError> {
        let reference = self.irradiation.fluxes.flux_at(time_s)?;
        if reference.groups().len() > budget.irradiation.source.max_neutron_groups {
            return Err(FissionError::Budget("neutron groups").into());
        }
        let mut values = Vec::with_capacity(reference.values_m2_s().len());
        for &value in reference.values_m2_s() {
            checkpoint()?;
            let scaled = finite(value * amplitude)?;
            if value > 0.0 && amplitude > 0.0 && scaled == 0.0 {
                return Err(FissionThermalError::Arithmetic);
            }
            values.push(scaled);
        }
        Ok(GroupFlux::new(Arc::clone(reference.groups()), values)?)
    }
}

impl sealed::Sealed for KineticCellSource {}
impl NuclearHeatSource for KineticCellSource {
    type Budget = KineticHeatingBudget;
    type Step = KineticCellStep;

    fn cell_id(&self) -> CellId {
        self.irradiation.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.irradiation.material_id
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        IrradiationCellSource::energy_tolerance_j(budget.irradiation)
    }
    fn next_transition(&self, time_s: f64) -> Option<f64> {
        [
            self.controls.next_after(time_s),
            self.irradiation.next_transition(time_s),
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
        self.irradiation
            .admit(time_s, temperature_k, budget.irradiation)?;
        let control = self.controls.control_at(time_s)?;
        let reactivity = self.reactivity_at(control, temperature_k, budget, &mut || Ok(()))?;
        self.kinetics.trial(
            0.0,
            reactivity,
            control.external_source_per_s,
            budget.kinetics,
            || Ok(()),
        )?;
        // Admit the actual initial flux, including its collapse against every
        // isotope, before a new or idle binding can be accepted.
        let flux = self.mean_flux(
            time_s,
            self.kinetics.neutron_population(),
            budget,
            &mut || Ok(()),
        )?;
        self.irradiation
            .trial_flux(
                &flux,
                0.0,
                temperature_k,
                budget.irradiation,
                &mut || Ok(()),
            )?;
        Ok(())
    }
    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        self.irradiation.validate_temperature(temperature_k)?;
        for fit in &self.feedback {
            fit.reactivity_shift(temperature_k)?;
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
        checkpoint()?;
        let control = self.controls.control_at(time_s)?;
        let reactivity = self.reactivity_at(control, temperature_k, budget, checkpoint)?;
        let mut checkpoint_error = None;
        let candidate = self.kinetics.trial(
            duration_s,
            reactivity,
            control.external_source_per_s,
            budget.kinetics,
            || {
                checkpoint().map_err(|error| {
                    checkpoint_error = Some(error);
                    KineticsError::Cancelled
                })
            },
        );
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let (kinetics, neutron_report) = candidate?;
        let mean = if duration_s == 0.0 {
            self.kinetics.neutron_population()
        } else {
            finite(neutron_report.integrated_population_s / duration_s)?
        };
        if mean == 0.0 && neutron_report.integrated_population_s > 0.0 {
            return Err(FissionThermalError::Arithmetic);
        }
        let mean_roundoff =
            finite((mean * duration_s - neutron_report.integrated_population_s).abs())?;
        if finite(neutron_report.integral_error_estimate_population_s + mean_roundoff)?
            > budget.kinetics.integral_tolerance_population_s
        {
            return Err(KineticsError::Budget("mean flux integral conversion").into());
        }
        let flux = self.mean_flux(time_s, mean, budget, checkpoint)?;
        let (irradiation, reaction_report) = self.irradiation.trial_flux(
            &flux,
            duration_s,
            temperature_k,
            budget.irradiation,
            checkpoint,
        )?;
        checkpoint()?;
        Ok((
            Self {
                irradiation,
                kinetics,
                controls: self.controls.clone(),
                feedback: self.feedback.clone(),
            },
            KineticCellStep {
                kinetics: neutron_report,
                irradiation: reaction_report,
                applied_reactivity: reactivity,
                mean_neutron_population: mean,
                mean_population_roundoff_s: mean_roundoff,
            },
        ))
    }
    fn released(step: &Self::Step) -> EnergyChannels {
        step.irradiation.released
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.irradiation.energy_roundoff_j
    }
    fn fissions(step: &Self::Step) -> f64 {
        step.irradiation.fissions
    }
    fn captures(step: &Self::Step) -> f64 {
        step.irradiation.captures
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        self.irradiation.cumulative_released()
    }
}
