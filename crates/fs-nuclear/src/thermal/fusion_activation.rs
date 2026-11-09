//! Finite activation and neutron kinetic deposition on the material heat clock.
//!
//! One complete cohort advances the original coupled neutron/nuclide owner
//! once over the accepted thermal interval. Only its combined deposited energy
//! heats material; the native kinetic, nuclear, radiation and leakage accounts
//! remain separately accessible. Material temperatures feed the original
//! microscopic capture tables, frozen over that interval.
//! Bound radioactive emitters share that finite inventory and continue without
//! external forcing. Their neutron kinetic energy is an internal transfer;
//! source audits compare only the external injection to the supplied history.

use std::{collections::BTreeSet, sync::Arc};

use super::source::{NuclearHeatSource, sealed};
use super::{
    FissionThermalError, FusionNeutronForcing, FusionNeutronPulse, NuclearHeating, Sum, finite,
    nonnegative,
};
use crate::burnup::BurnupStep;
use crate::diffusion::burnup::DiffusionBurnupError;
use crate::fusion_neutrons::activation::{
    FusionNeutronActivationBudget, FusionNeutronActivationError, FusionNeutronActivationState,
    FusionNeutronActivationStep,
};
use crate::fusion_neutrons::history::{FusionNeutronHistory, FusionNeutronHistoryLimits};
use crate::{EnergyChannels, FissionError};
use fs_qty::Time;
use fs_thermal::enthalpy::{CellId, MaterialId};

/// One original cell of a complete shared neutron and finite-material history.
#[derive(Clone, Debug)]
pub struct FusionNeutronActivationCellSource {
    cell_id: CellId,
    material_id: MaterialId,
    cell_index: usize,
    state: Arc<FusionNeutronActivationState>,
    forcing: FusionNeutronForcing,
    accepted_start_s: f64,
    accepted_duration_s: f64,
}

/// Original material calorics with finite activation, delayed and kinetic heat.
pub type FusionNeutronActivationHeating = NuclearHeating<FusionNeutronActivationCellSource>;

/// A cell view of the one shared accepted physical trial.
#[derive(Clone, Debug)]
pub struct FusionNeutronActivationCellStep {
    cell_index: usize,
    activation: Arc<FusionNeutronActivationStep>,
    forcing: FusionNeutronForcing,
    start_s: f64,
    source_number_residual: f64,
    source_energy_residual_j: f64,
    energy_roundoff_j: f64,
    population_roundoff_nuclei: f64,
}

impl FusionNeutronActivationCellStep {
    pub const fn cell_index(&self) -> usize {
        self.cell_index
    }
    /// Original coupled report, including separate kinetic and nuclear exports.
    pub fn activation(&self) -> &FusionNeutronActivationStep {
        &self.activation
    }
    /// Original finite-target, product and delayed-release report for this cell.
    pub fn burnup(&self) -> &BurnupStep {
        &self.activation.spatial.burnup[self.cell_index]
    }
    pub const fn forcing(&self) -> &FusionNeutronForcing {
        &self.forcing
    }
    pub const fn start(&self) -> Time {
        Time::new(self.start_s)
    }
    /// Global native injection minus this interval's original source count.
    pub const fn source_number_residual(&self) -> f64 {
        self.source_number_residual
    }
    /// Global native kinetic injection minus the original interval energy [J].
    pub const fn source_energy_residual_j(&self) -> f64 {
        self.source_energy_residual_j
    }
}

impl FusionNeutronActivationCellSource {
    /// Continue resident neutrons and finite radioactive emitters without any
    /// imposed neutron source. An explicit empty history reuses the same
    /// forcing, shared-cohort and native material transaction as supplied pulses.
    pub fn bind_without_external_source(
        state: Arc<FusionNeutronActivationState>,
        accepted_time: Time,
        cells: Vec<(CellId, MaterialId)>,
    ) -> Result<Vec<Self>, FissionThermalError> {
        let history = Arc::new(FusionNeutronHistory::project(
            Arc::clone(state.model()),
            Vec::new(),
            FusionNeutronHistoryLimits::default(),
        )?);
        Self::bind_with_history(state, accepted_time, cells, history)
    }

    /// Bind the supplied accepted material/neutron checkpoint to a finite pulse.
    /// Earlier irradiation and heat are never replayed. Every original neutron
    /// cell must map to one distinct active thermal cell of the original volume.
    pub fn bind(
        state: Arc<FusionNeutronActivationState>,
        accepted_time: Time,
        cells: Vec<(CellId, MaterialId)>,
        pulse: Arc<FusionNeutronPulse>,
    ) -> Result<Vec<Self>, FissionThermalError> {
        Self::bind_forcing(
            state,
            accepted_time,
            cells,
            FusionNeutronForcing::Pulse(pulse),
        )
    }

    /// Bind changing accepted emission in any subset of the spatial cells.
    /// Complete material coverage remains mandatory, including cells that only
    /// absorb or retain activated products after the external source shuts off.
    pub fn bind_with_history(
        state: Arc<FusionNeutronActivationState>,
        accepted_time: Time,
        cells: Vec<(CellId, MaterialId)>,
        history: Arc<FusionNeutronHistory>,
    ) -> Result<Vec<Self>, FissionThermalError> {
        Self::bind_forcing(
            state,
            accepted_time,
            cells,
            FusionNeutronForcing::History(history),
        )
    }

    fn bind_forcing(
        state: Arc<FusionNeutronActivationState>,
        accepted_time: Time,
        cells: Vec<(CellId, MaterialId)>,
        forcing: FusionNeutronForcing,
    ) -> Result<Vec<Self>, FissionThermalError> {
        if !nonnegative(accepted_time.value())
            || cells.len() != state.inventories().len()
            || **state.model() != **forcing.model()
        {
            return Err(FissionThermalError::Invalid(
                "matching activation source model, complete cell map and finite accepted time required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut sources = Vec::with_capacity(cells.len());
        for (cell_index, (cell_id, material_id)) in cells.into_iter().enumerate() {
            if !used.insert(cell_id) {
                return Err(FissionThermalError::Invalid(
                    "distinct activation thermal cells required",
                ));
            }
            sources.push(Self {
                cell_id,
                material_id,
                cell_index,
                state: Arc::clone(&state),
                forcing: forcing.clone(),
                accepted_start_s: accepted_time.value(),
                accepted_duration_s: 0.0,
            });
        }
        Ok(sources)
    }

    pub const fn state(&self) -> &Arc<FusionNeutronActivationState> {
        &self.state
    }
    pub const fn forcing(&self) -> &FusionNeutronForcing {
        &self.forcing
    }
    pub const fn cell_index(&self) -> usize {
        self.cell_index
    }

    fn admit_time(&self, time_s: f64) -> Result<(), FissionThermalError> {
        // The thermal endpoint can be the sum of start and duration or the
        // authored transition from which that duration was subtracted.
        if !nonnegative(time_s)
            || time_s < self.accepted_start_s
            || (time_s != self.accepted_start_s + self.accepted_duration_s
                && time_s - self.accepted_start_s != self.accepted_duration_s)
        {
            return Err(FissionThermalError::Invalid(
                "activation checkpoint and material must share the accepted clock",
            ));
        }
        Ok(())
    }
}

impl sealed::Sealed for FusionNeutronActivationCellSource {}
impl NuclearHeatSource for FusionNeutronActivationCellSource {
    type Budget = FusionNeutronActivationBudget;
    type Step = FusionNeutronActivationCellStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn expected_volume_m3(&self) -> Option<f64> {
        Some(self.state.model().diffusion().cell_volumes_m3()[self.cell_index])
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.energy_tolerance_j
    }
    fn nuclei_tolerance(budget: Self::Budget) -> Option<f64> {
        Some(crate::sum(&[
            budget.spatial.diffusion.population_tolerance,
            budget.spatial.burnup.nuclei_tolerance,
            budget.spatial.reaction_tolerance,
        ]))
    }
    fn next_transition(&self, time_s: f64) -> Option<f64> {
        self.forcing.next_after(time_s)
    }
    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError> {
        self.admit_time(time_s)?;
        budget.spatial.validate()?;
        if !nonnegative(budget.energy_tolerance_j)
            || !Self::nuclei_tolerance(budget).is_some_and(nonnegative)
        {
            return Err(FissionThermalError::Invalid(
                "finite nonnegative combined activation allowances required",
            ));
        }
        self.validate_temperature(temperature_k)?;
        if self.cell_index == 0 {
            self.state
                .spatial()
                .admit_transport(budget.spatial.burnup, &mut || Ok(()))?;
            self.state
                .neutrons()
                .trial_with_sources(
                    0.0,
                    self.forcing.rates_at(time_s)?,
                    &vec![0.0; self.state.neutrons().precursor_counts().len()],
                    budget.spatial.diffusion,
                    || Ok(()),
                )
                .map_err(DiffusionBurnupError::from)?;
        }
        let flux = self
            .state
            .neutrons()
            .group_flux(self.cell_index)
            .map_err(DiffusionBurnupError::from)?;
        self.state.inventories()[self.cell_index].trial(
            0.0,
            &flux,
            temperature_k,
            budget.spatial.burnup,
            || Ok(()),
        )?;
        Ok(())
    }
    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        self.state.inventories()[self.cell_index]
            .model()
            .validate_temperature(temperature_k)?;
        Ok(())
    }
    fn validate_sources<'source>(
        mut sources: impl ExactSizeIterator<Item = &'source Self>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionThermalError>,
    ) -> Result<(), FissionThermalError>
    where
        Self: 'source,
    {
        checkpoint()?;
        let count = sources.len();
        let first = sources
            .next()
            .ok_or(FissionThermalError::Invalid("empty activation cohort"))?;
        if count != first.state.inventories().len() {
            return Err(FissionThermalError::Invalid(
                "all activation cells must advance together",
            ));
        }
        let mut seen = vec![false; count];
        for source in std::iter::once(first).chain(sources) {
            checkpoint()?;
            if !Arc::ptr_eq(&first.state, &source.state)
                || !first.forcing.same_history(&source.forcing)
                || first.accepted_start_s != source.accepted_start_s
                || first.accepted_duration_s != source.accepted_duration_s
                || source.cell_index >= count
                || seen[source.cell_index]
            {
                return Err(FissionThermalError::Invalid(
                    "activation cells must share one complete accepted state, source and clock",
                ));
            }
            seen[source.cell_index] = true;
        }
        Ok(())
    }
    fn trial_source(
        &self,
        _: f64,
        _: f64,
        _: f64,
        _: Self::Budget,
        _: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        Err(FissionThermalError::Invalid(
            "fusion activation requires the complete cell collection",
        ))
    }
    fn trial_sources(
        sources: &[&Self],
        time_s: f64,
        duration_s: f64,
        temperatures_k: &[f64],
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<(Self, Self::Step)>, FissionThermalError> {
        Self::validate_sources(sources.iter().copied(), &mut || {
            checkpoint().map_err(Into::into)
        })?;
        let first = sources[0];
        first.admit_time(time_s)?;
        if temperatures_k.len() != sources.len()
            || !nonnegative(duration_s)
            || !(time_s + duration_s).is_finite()
            || (duration_s > 0.0 && time_s + duration_s <= time_s)
        {
            return Err(FissionThermalError::Invalid(
                "one temperature per activation cell and finite duration required",
            ));
        }
        if first
            .forcing
            .next_after(time_s)
            .is_some_and(|next| duration_s > next - time_s)
        {
            return Err(FissionThermalError::Invalid(
                "activation thermal interval crosses a fusion source transition",
            ));
        }
        let mut temperatures = vec![0.0; sources.len()];
        for (source, &temperature) in sources.iter().zip(temperatures_k) {
            checkpoint()?;
            source.validate_temperature(temperature)?;
            temperatures[source.cell_index] = temperature;
        }
        let (state, report) = first.state.trial_with_rates(
            duration_s,
            &temperatures,
            first.forcing.rates_at(time_s)?,
            budget,
            &mut *checkpoint,
        )?;
        let (expected_number, expected_energy) =
            first.forcing.interval_moments(time_s, duration_s)?;
        let source_number_residual = finite(report.external_injected_neutrons - expected_number)?;
        let source_energy_residual_j =
            finite(report.external_injected_kinetic_energy_j - expected_energy)?;
        let mut energy_defect = Sum::default();
        energy_defect.add(report.energy_roundoff_j)?;
        energy_defect.add(source_energy_residual_j.abs())?;
        let energy_roundoff_j = energy_defect.nonnegative_upper_bound()?;
        let mut neutron_defect = Sum::default();
        neutron_defect.add(report.spatial.neutrons.population_roundoff)?;
        neutron_defect.add(report.injection_number_roundoff)?;
        neutron_defect.add(source_number_residual.abs())?;
        if energy_roundoff_j > budget.energy_tolerance_j
            || neutron_defect.nonnegative_upper_bound()?
                > budget.spatial.diffusion.population_tolerance
        {
            return Err(FissionThermalError::Invalid(
                "original fusion source injection exceeds activation account allowance",
            ));
        }
        let state = Arc::new(state);
        let activation = Arc::new(report);
        let mut candidates = Vec::with_capacity(sources.len());
        for source in sources {
            checkpoint()?;
            let cell = source.cell_index;
            let mut population_defect = Sum::default();
            population_defect.add(activation.spatial.burnup[cell].population_roundoff_nuclei)?;
            population_defect.add(activation.spatial.reaction_defect_by_cell[cell])?;
            if cell == 0 {
                population_defect.add(neutron_defect.nonnegative_upper_bound()?)?;
            }
            let mut candidate = (*source).clone();
            candidate.state = Arc::clone(&state);
            candidate.accepted_start_s = time_s;
            candidate.accepted_duration_s = duration_s;
            candidates.push((
                candidate,
                FusionNeutronActivationCellStep {
                    cell_index: cell,
                    activation: Arc::clone(&activation),
                    forcing: source.forcing.clone(),
                    start_s: time_s,
                    source_number_residual,
                    source_energy_residual_j,
                    // This is one global account, regardless of thermal ID order.
                    energy_roundoff_j: if cell == 0 { energy_roundoff_j } else { 0.0 },
                    population_roundoff_nuclei: population_defect.nonnegative_upper_bound()?,
                },
            ));
        }
        checkpoint()?;
        Ok(candidates)
    }
    fn released(step: &Self::Step) -> EnergyChannels {
        step.activation.cell_released[step.cell_index]
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.energy_roundoff_j
    }
    fn population_roundoff_nuclei(step: &Self::Step) -> f64 {
        step.population_roundoff_nuclei
    }
    fn decays(step: &Self::Step) -> f64 {
        step.burnup().decays
    }
    fn produced_nuclei(step: &Self::Step) -> f64 {
        step.activation.spatial.produced_nuclei_by_cell[step.cell_index]
    }
    fn fissions(step: &Self::Step) -> f64 {
        step.burnup().fissions
    }
    fn captures(step: &Self::Step) -> f64 {
        // Finite-inventory events keep the native burnup counter. Fixed
        // background captures remain visible in the separate neutron report.
        step.burnup().captures
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        self.state
            .cell_released(self.cell_index)
            .ok_or(FissionThermalError::Invalid(
                "original activation cell absent from accepted state",
            ))
    }
}

impl From<FusionNeutronActivationError> for FissionThermalError {
    fn from(error: FusionNeutronActivationError) -> Self {
        match error {
            FusionNeutronActivationError::Source(error) => Self::Fission(error),
            FusionNeutronActivationError::Kinetic(error) => error.into(),
            FusionNeutronActivationError::Spatial(error) => error.into(),
            other => Self::FusionActivation(other),
        }
    }
}
