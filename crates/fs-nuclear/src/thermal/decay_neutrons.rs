//! Original radioactive chains, emitted-neutron transport and material heating.
//!
//! One complete spatial cohort advances once per accepted thermal interval.
//! The energy owner has already transferred emitted neutron energy out of the
//! radioactive radiation account. Only its disjoint deposited channel is heat;
//! neutron/precursor storage, later capture, outward radiation and leakage keep
//! their own accounts. Funded background fission also retains its actual event
//! count and delayed reserve; only accepted deposited energy heats material. The bound thermal law and physical mass are unchanged. Nuclear
//! coefficients remain authored and fixed; material temperature does not infer
//! a new half-life, neutron cross section or deposition fraction.

use std::{collections::BTreeSet, sync::Arc};

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NuclearHeating, Sum, finite, nonnegative};
use crate::diffusion::decay::DecayNeutronError;
use crate::diffusion::decay::energy::{
    DecayNeutronEnergyBudget, DecayNeutronEnergyError, DecayNeutronEnergyState,
    DecayNeutronEnergyStep,
};
use crate::diffusion::{DiffusionError, DiffusionReactionKind};
use crate::{EnergyChannels, FissionError};
use fs_qty::Time;
use fs_thermal::enthalpy::{CellId, MaterialId};

/// One original spatial cell in a complete, shared decay/neutron checkpoint.
#[derive(Clone, Debug)]
pub struct DecayNeutronCellSource {
    cell_id: CellId,
    material_id: MaterialId,
    cell_index: usize,
    state: Arc<DecayNeutronEnergyState>,
    // Native subtraction and addition may round differently when the thermal
    // owner snaps to an authored endpoint. Retain its original interval pair.
    accepted_start_s: f64,
    accepted_duration_s: f64,
}

/// Material heat from actual radioactive decay and its emitted neutron history.
pub type DecayNeutronHeating = NuclearHeating<DecayNeutronCellSource>;

/// A cell view of the one jointly accepted physical interval.
#[derive(Clone, Debug)]
pub struct DecayNeutronCellStep {
    cell_index: usize,
    spatial: Arc<DecayNeutronEnergyStep>,
    decays: f64,
    captures: f64,
    fissions: f64,
    count_roundoff: f64,
}

impl DecayNeutronCellStep {
    pub const fn cell_index(&self) -> usize {
        self.cell_index
    }
    /// Original count, kinetic-energy and radioactive-energy reports.
    pub fn spatial(&self) -> &DecayNeutronEnergyStep {
        &self.spatial
    }
}

impl DecayNeutronCellSource {
    /// Bind every neutron cell to a distinct original thermal cell/material.
    /// The supplied complete checkpoint represents `accepted_time`; neither
    /// earlier decay nor earlier heating is replayed. The thermal owner admits
    /// the same clock, reference volumes and immutable caloric material laws.
    pub fn bind(
        state: Arc<DecayNeutronEnergyState>,
        accepted_time: Time,
        cells: Vec<(CellId, MaterialId)>,
    ) -> Result<Vec<Self>, FissionThermalError> {
        if !nonnegative(accepted_time.value())
            || cells.len() != state.model().diffusion().cell_count()
        {
            return Err(FissionThermalError::Invalid(
                "complete decay-neutron thermal cell map and finite accepted time required",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut sources = Vec::with_capacity(cells.len());
        for (cell_index, (cell_id, material_id)) in cells.into_iter().enumerate() {
            if !seen.insert(cell_id) {
                return Err(FissionThermalError::Invalid(
                    "distinct decay-neutron thermal cells required",
                ));
            }
            sources.push(Self {
                cell_id,
                material_id,
                cell_index,
                state: Arc::clone(&state),
                accepted_start_s: accepted_time.value(),
                accepted_duration_s: 0.0,
            });
        }
        Ok(sources)
    }
    pub const fn state(&self) -> &Arc<DecayNeutronEnergyState> {
        &self.state
    }
    pub const fn cell_index(&self) -> usize {
        self.cell_index
    }
    fn admit_time(&self, time_s: f64) -> Result<(), FissionThermalError> {
        if !nonnegative(time_s)
            || time_s < self.accepted_start_s
            || (time_s != self.accepted_start_s + self.accepted_duration_s
                && time_s - self.accepted_start_s != self.accepted_duration_s)
        {
            return Err(FissionThermalError::Invalid(
                "decay-neutron checkpoint and material must share the accepted clock",
            ));
        }
        Ok(())
    }
}

impl sealed::Sealed for DecayNeutronCellSource {}
impl NuclearHeatSource for DecayNeutronCellSource {
    type Budget = DecayNeutronEnergyBudget;
    type Step = DecayNeutronCellStep;

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
            budget.spatial.decay.nuclei_tolerance,
        ]))
    }
    fn next_transition(&self, _: f64) -> Option<f64> {
        None
    }
    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError> {
        self.admit_time(time_s)?;
        self.validate_temperature(temperature_k)?;
        if self.cell_index == 0 {
            self.state.trial(0.0, budget, || Ok(()))?;
        }
        Ok(())
    }
    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(FissionThermalError::Invalid(
                "positive finite material temperature required",
            ));
        }
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
            .ok_or(FissionThermalError::Invalid("empty decay-neutron cohort"))?;
        if count != first.state.model().diffusion().cell_count() {
            return Err(FissionThermalError::Invalid(
                "all decay-neutron cells must advance together",
            ));
        }
        let mut seen = vec![false; count];
        for source in std::iter::once(first).chain(sources) {
            checkpoint()?;
            if !Arc::ptr_eq(&first.state, &source.state)
                || first.accepted_start_s != source.accepted_start_s
                || first.accepted_duration_s != source.accepted_duration_s
                || source.cell_index >= count
                || seen[source.cell_index]
            {
                return Err(FissionThermalError::Invalid(
                    "decay-neutron cells must share one complete state and clock",
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
            "decay-neutron evolution requires the complete cell collection",
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
                "one temperature per decay-neutron cell and finite duration required",
            ));
        }
        for (source, &temperature) in sources.iter().zip(temperatures_k) {
            checkpoint()?;
            source.validate_temperature(temperature)?;
        }
        let (state, report) = first.state.trial(duration_s, budget, &mut *checkpoint)?;
        let mut counts = vec![[Sum::default(); 3]; sources.len()];
        for (source, decay) in first
            .state
            .spatial()
            .sources()
            .iter()
            .zip(&report.spatial.decay)
        {
            for &count in &decay.decays {
                checkpoint()?;
                counts[source.cell][0].add(count)?;
            }
        }
        for (reaction, &count) in first
            .state
            .model()
            .diffusion()
            .reactions()
            .iter()
            .zip(&report.spatial.neutrons.reaction_counts)
        {
            checkpoint()?;
            match reaction.kind {
                DiffusionReactionKind::Capture => counts[reaction.cell][1].add(count)?,
                DiffusionReactionKind::Fission { .. } => counts[reaction.cell][2].add(count)?,
                _ => {}
            }
        }
        let mut defect = Sum::default();
        defect.add(report.spatial.neutrons.population_roundoff)?;
        defect.add(report.spatial.decay_population_roundoff)?;
        let mut count_roundoff = Vec::with_capacity(sources.len());
        for count in &counts {
            let value = finite(crate::sum(&[
                count[0].export_roundoff()?,
                count[1].export_roundoff()?,
                count[2].export_roundoff()?,
            ]))?;
            defect.add(value)?;
            count_roundoff.push(value);
        }
        let allowance = Self::nuclei_tolerance(budget).expect("combined population allowance");
        if defect.nonnegative_upper_bound()? > allowance {
            return Err(FissionThermalError::NucleiBudget {
                defect_nuclei: defect.nonnegative_upper_bound()?,
                tolerance_nuclei: allowance,
            });
        }
        let state = Arc::new(state);
        let spatial = Arc::new(report);
        let mut candidates = Vec::with_capacity(sources.len());
        for source in sources {
            checkpoint()?;
            let mut candidate = (*source).clone();
            candidate.state = Arc::clone(&state);
            candidate.accepted_start_s = time_s;
            candidate.accepted_duration_s = duration_s;
            candidates.push((
                candidate,
                DecayNeutronCellStep {
                    cell_index: source.cell_index,
                    spatial: Arc::clone(&spatial),
                    decays: counts[source.cell_index][0].value()?,
                    captures: counts[source.cell_index][1].value()?,
                    fissions: counts[source.cell_index][2].value()?,
                    count_roundoff: count_roundoff[source.cell_index],
                },
            ));
        }
        Ok(candidates)
    }
    fn released(step: &Self::Step) -> EnergyChannels {
        step.spatial.cell_released[step.cell_index]
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        // A shared account is charged once in original spatial order, even
        // when thermal cell IDs order receivers differently.
        if step.cell_index == 0 {
            step.spatial.energy_roundoff_j
        } else {
            0.0
        }
    }
    fn population_roundoff_nuclei(step: &Self::Step) -> f64 {
        if step.cell_index == 0 {
            crate::sum(&[
                step.spatial.spatial.neutrons.population_roundoff,
                step.spatial.spatial.decay_population_roundoff,
                step.count_roundoff,
            ])
        } else {
            step.count_roundoff
        }
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
        self.state
            .cell_released(self.cell_index)
            .ok_or(FissionThermalError::Invalid(
                "original decay-neutron cell absent from accepted state",
            ))
    }
}

impl From<DecayNeutronEnergyError> for FissionThermalError {
    fn from(error: DecayNeutronEnergyError) -> Self {
        match error {
            DecayNeutronEnergyError::Fission(error)
            | DecayNeutronEnergyError::Source(DecayNeutronError::Fission(error))
            | DecayNeutronEnergyError::Source(DecayNeutronError::Diffusion(
                DiffusionError::Fission(error),
            )) => Self::Fission(error),
            DecayNeutronEnergyError::Kinetic(error) => error.into(),
            DecayNeutronEnergyError::Source(DecayNeutronError::Decay(error)) => error.into(),
            other => Self::DecayNeutrons(other),
        }
    }
}

#[cfg(test)]
mod tests;
