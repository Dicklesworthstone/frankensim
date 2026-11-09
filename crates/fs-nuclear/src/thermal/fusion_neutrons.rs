//! Fusion-neutron deposition on the existing transactional thermal clock.
//!
//! One complete spatial cohort advances the original neutron owner once per
//! accepted thermal interval. Its local heat reaches the bound material cells;
//! neutron and precursor storage, outgoing radiation and leakage retain their
//! own accounts. Explicitly funded delayed fission follows the same accepted
//! native history; its future neutron energy becomes heat only upon deposition.
//! A single pulse or accepted-interval history retains its original channel
//! moments. Each source interval uses an explicit uniform/mean-energy projection.
//! This adapter supplies no new nuclear data or
//! temperature feedback to the authored neutron cross sections.

use std::{collections::BTreeSet, sync::Arc};

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NuclearHeating, finite, nonnegative};
use crate::diffusion::DiffusionReactionKind;
use crate::fusion_neutrons::history::FusionNeutronHistory;
use crate::fusion_neutrons::{
    FusionNeutronSource, NeutronKineticBudget, NeutronKineticModel, NeutronKineticState,
    NeutronKineticStep,
};
use crate::{EnergyChannels, FissionError};
use fs_qty::Time;
use fs_thermal::enthalpy::{CellId, MaterialId};

/// A single finite fusion-emission pulse on the absolute thermal clock.
/// The original projected source, duration, particle count and energy remain
/// immutable. Its rates apply on [start, end); zero rates apply outside it.
#[derive(Clone, Debug)]
pub struct FusionNeutronPulse {
    start_s: f64,
    end_s: f64,
    source: FusionNeutronSource,
    zero_rates_s: Vec<f64>,
}

impl FusionNeutronPulse {
    pub fn new(start: Time, source: FusionNeutronSource) -> Result<Self, FissionThermalError> {
        let start_s = start.value();
        let duration = source.emission().duration().value();
        let end_s = start_s + duration;
        if !nonnegative(start_s)
            || !end_s.is_finite()
            || (duration > 0.0 && end_s <= start_s)
            || ((end_s - start_s) - duration).abs() > 64.0 * f64::EPSILON * duration
        {
            return Err(FissionThermalError::Invalid(
                "finite pulse start and absolute clock preserving original duration to roundoff required",
            ));
        }
        Ok(Self {
            start_s,
            end_s,
            zero_rates_s: vec![0.0; source.neutron_rates_s().len()],
            source,
        })
    }

    pub const fn start(&self) -> Time {
        Time::new(self.start_s)
    }
    pub const fn end(&self) -> Time {
        Time::new(self.end_s)
    }
    /// The unchanged projection and its original channel-resolved emission.
    pub const fn source(&self) -> &FusionNeutronSource {
        &self.source
    }
    /// Active constant rates in original cell/group order. Calling after
    /// shutoff does not repeat the already emitted fusion interval.
    pub fn rates_at(&self, time: Time) -> Result<&[f64], FissionThermalError> {
        if !nonnegative(time.value()) {
            return Err(FissionThermalError::Invalid(
                "finite nonnegative pulse clock required",
            ));
        }
        Ok(if self.active_at(time.value()) {
            self.source.neutron_rates_s()
        } else {
            &self.zero_rates_s
        })
    }
    fn active_at(&self, time_s: f64) -> bool {
        time_s >= self.start_s && time_s < self.end_s
    }
    fn next_after(&self, time_s: f64) -> Option<f64> {
        if self.source.emission().total_particle_count() == 0.0 {
            None
        } else if time_s < self.start_s {
            Some(self.start_s)
        } else if time_s < self.end_s {
            Some(self.end_s)
        } else {
            None
        }
    }
}

/// Original immutable forcing retained by a complete neutron/thermal cohort.
/// Histories retain changing emission in every emitting cell, while Pulse is
/// the explicitly uniform single-interval source supported by the original API.
#[derive(Clone, Debug)]
pub enum FusionNeutronForcing {
    Pulse(Arc<FusionNeutronPulse>),
    History(Arc<FusionNeutronHistory>),
}

impl FusionNeutronForcing {
    pub(super) fn model(&self) -> &Arc<NeutronKineticModel> {
        match self {
            Self::Pulse(pulse) => pulse.source().model(),
            Self::History(history) => history.model(),
        }
    }
    pub(super) fn rates_at(&self, time_s: f64) -> Result<&[f64], FissionThermalError> {
        match self {
            Self::Pulse(pulse) => pulse.rates_at(Time::new(time_s)),
            Self::History(history) => Ok(history.rates_at(Time::new(time_s))?),
        }
    }
    pub(super) fn next_after(&self, time_s: f64) -> Option<f64> {
        match self {
            Self::Pulse(pulse) => pulse.next_after(time_s),
            Self::History(history) => history
                .next_after(Time::new(time_s))
                .map(|time| time.value()),
        }
    }
    pub(super) fn interval_moments(
        &self,
        time_s: f64,
        duration_s: f64,
    ) -> Result<(f64, f64), FissionThermalError> {
        match self {
            Self::Pulse(pulse) => {
                let fraction = if pulse.active_at(time_s) {
                    duration_s / pulse.source.emission().duration().value()
                } else {
                    0.0
                };
                Ok((
                    finite(pulse.source.emission().total_particle_count() * fraction)?,
                    finite(pulse.source.emission().total_kinetic_energy().value() * fraction)?,
                ))
            }
            Self::History(history) => {
                Ok(history.interval_moments(Time::new(time_s), Time::new(duration_s))?)
            }
        }
    }
    pub(super) fn same_history(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Pulse(a), Self::Pulse(b)) => Arc::ptr_eq(a, b),
            (Self::History(a), Self::History(b)) => Arc::ptr_eq(a, b),
            _ => false,
        }
    }
    pub fn pulse(&self) -> Option<&Arc<FusionNeutronPulse>> {
        if let Self::Pulse(pulse) = self {
            Some(pulse)
        } else {
            None
        }
    }
    pub fn history(&self) -> Option<&Arc<FusionNeutronHistory>> {
        if let Self::History(history) = self {
            Some(history)
        } else {
            None
        }
    }
}

/// A handle to one cell of a complete immutable neutron/thermal cohort.
#[derive(Clone, Debug)]
pub struct NeutronKineticCellSource {
    cell_id: CellId,
    material_id: MaterialId,
    cell_index: usize,
    state: Arc<NeutronKineticState>,
    forcing: FusionNeutronForcing,
    // Preserve the native start/duration pair. The thermal owner either adds
    // that duration or snaps to the authored endpoint from which it subtracted
    // the duration; those two floating-point paths need not agree bitwise.
    accepted_start_s: f64,
    accepted_duration_s: f64,
}

/// Existing conservative material thermal evolution with fusion-neutron heat.
pub type NeutronKineticHeating = NuclearHeating<NeutronKineticCellSource>;

/// One cell's view of the shared accepted neutron interval and original pulse.
#[derive(Clone, Debug)]
pub struct NeutronKineticCellStep {
    cell_index: usize,
    spatial: Arc<NeutronKineticStep>,
    forcing: FusionNeutronForcing,
    start_s: f64,
    captures: f64,
    capture_roundoff: f64,
    fissions: f64,
    fission_roundoff: f64,
    source_number_residual: f64,
    source_energy_residual_j: f64,
}

impl NeutronKineticCellStep {
    pub const fn cell_index(&self) -> usize {
        self.cell_index
    }
    pub fn spatial(&self) -> &NeutronKineticStep {
        &self.spatial
    }
    pub const fn forcing(&self) -> &FusionNeutronForcing {
        &self.forcing
    }
    pub fn pulse(&self) -> Option<&Arc<FusionNeutronPulse>> {
        self.forcing.pulse()
    }
    pub fn history(&self) -> Option<&Arc<FusionNeutronHistory>> {
        self.forcing.history()
    }
    pub const fn start(&self) -> Time {
        Time::new(self.start_s)
    }
    /// Global native injection minus the original pulse's interval share.
    pub const fn source_number_residual(&self) -> f64 {
        self.source_number_residual
    }
    /// Global native injected energy minus the original pulse's interval share.
    pub const fn source_energy_residual_j(&self) -> f64 {
        self.source_energy_residual_j
    }
}

impl NeutronKineticCellSource {
    /// Bind every original neutron cell to a distinct active thermal cell.
    /// The supplied neutron checkpoint represents `accepted_time`; no earlier
    /// part of the pulse is replayed. The thermal owner checks this time, cell
    /// volumes, masses and original material identities before it can advance.
    pub fn bind(
        state: Arc<NeutronKineticState>,
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
    /// Bind time-resolved, possibly concurrent emission in multiple cells.
    /// The original thermal clock clips every retained source transition; no
    /// earlier source history is replayed when binding an accepted checkpoint.
    pub fn bind_with_history(
        state: Arc<NeutronKineticState>,
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
        state: Arc<NeutronKineticState>,
        accepted_time: Time,
        cells: Vec<(CellId, MaterialId)>,
        forcing: FusionNeutronForcing,
    ) -> Result<Vec<Self>, FissionThermalError> {
        if !nonnegative(accepted_time.value())
            || cells.len() != state.model().diffusion().cell_count()
            || **state.model() != **forcing.model()
        {
            return Err(FissionThermalError::Invalid(
                "matching neutron source model, complete thermal cell map and finite accepted time required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut sources = Vec::with_capacity(cells.len());
        for (cell_index, (cell_id, material_id)) in cells.into_iter().enumerate() {
            if !used.insert(cell_id) {
                return Err(FissionThermalError::Invalid(
                    "distinct neutron thermal cells required",
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
    pub const fn state(&self) -> &Arc<NeutronKineticState> {
        &self.state
    }
    pub const fn forcing(&self) -> &FusionNeutronForcing {
        &self.forcing
    }
    pub fn pulse(&self) -> Option<&Arc<FusionNeutronPulse>> {
        self.forcing.pulse()
    }
    pub fn history(&self) -> Option<&Arc<FusionNeutronHistory>> {
        self.forcing.history()
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
                "neutron checkpoint and thermal interval must share the accepted clock",
            ));
        }
        Ok(())
    }
}

impl sealed::Sealed for NeutronKineticCellSource {}
impl NuclearHeatSource for NeutronKineticCellSource {
    type Budget = NeutronKineticBudget;
    type Step = NeutronKineticCellStep;

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
        Some(budget.diffusion.population_tolerance)
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
        self.validate_temperature(temperature_k)?;
        if self.cell_index == 0 {
            self.state
                .trial_with_rates(0.0, self.forcing.rates_at(time_s)?, budget, || Ok(()))?;
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
            .ok_or(FissionThermalError::Invalid("empty neutron cohort"))?;
        if count != first.state.model().diffusion().cell_count() {
            return Err(FissionThermalError::Invalid(
                "all neutron cells must advance together",
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
                    "neutron cells must share one complete accepted state, pulse and clock",
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
            "neutron evolution requires the complete cell collection",
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
                "one temperature per neutron cell and finite duration required",
            ));
        }
        if first
            .forcing
            .next_after(time_s)
            .is_some_and(|next| duration_s > next - time_s)
        {
            return Err(FissionThermalError::Invalid(
                "neutron thermal interval crosses a fusion pulse boundary",
            ));
        }
        for (source, &temperature) in sources.iter().zip(temperatures_k) {
            checkpoint()?;
            source.validate_temperature(temperature)?;
        }
        let (state, report) = first.state.trial_with_rates(
            duration_s,
            first.forcing.rates_at(time_s)?,
            budget,
            &mut *checkpoint,
        )?;
        let (expected_number, expected_energy) =
            first.forcing.interval_moments(time_s, duration_s)?;
        let source_number_residual = finite(report.diffusion.injected_neutrons - expected_number)?;
        let source_energy_residual_j = finite(report.injected_kinetic_energy_j - expected_energy)?;
        if finite(report.energy_roundoff_j + source_energy_residual_j.abs())?
            > budget.energy_tolerance_j
            || finite(report.diffusion.population_roundoff + source_number_residual.abs())?
                > budget.diffusion.population_tolerance
        {
            return Err(FissionThermalError::Invalid(
                "fusion pulse injection account exceeds neutron budget",
            ));
        }
        let mut captures = vec![super::Sum::default(); sources.len()];
        let mut fissions = vec![super::Sum::default(); sources.len()];
        for (reaction, &count) in first
            .state
            .model()
            .diffusion()
            .reactions()
            .iter()
            .zip(&report.diffusion.reaction_counts)
        {
            checkpoint()?;
            match reaction.kind {
                DiffusionReactionKind::Capture => captures[reaction.cell].add(count)?,
                DiffusionReactionKind::Fission { .. } => fissions[reaction.cell].add(count)?,
                _ => {}
            }
        }
        let mut population_defect = super::Sum::default();
        population_defect.add(report.diffusion.population_roundoff)?;
        population_defect.add(source_number_residual.abs())?;
        for account in captures.iter().chain(&fissions) {
            population_defect.add(account.export_roundoff()?)?;
        }
        if population_defect.nonnegative_upper_bound()? > budget.diffusion.population_tolerance {
            return Err(FissionThermalError::Invalid(
                "fusion pulse count and capture export exceed neutron budget",
            ));
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
                NeutronKineticCellStep {
                    cell_index: source.cell_index,
                    spatial: Arc::clone(&spatial),
                    forcing: source.forcing.clone(),
                    start_s: time_s,
                    captures: captures[source.cell_index].value()?,
                    capture_roundoff: captures[source.cell_index].export_roundoff()?,
                    fissions: fissions[source.cell_index].value()?,
                    fission_roundoff: fissions[source.cell_index].export_roundoff()?,
                    source_number_residual,
                    source_energy_residual_j,
                },
            ));
        }
        Ok(candidates)
    }
    fn released(step: &Self::Step) -> EnergyChannels {
        step.spatial.cell_released[step.cell_index]
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        // The underlying account is global: include it exactly once, even if
        // thermal cell IDs sort differently from the original neutron cells.
        if step.cell_index == 0 {
            step.spatial.energy_roundoff_j + step.source_energy_residual_j.abs()
        } else {
            0.0
        }
    }
    fn population_roundoff_nuclei(step: &Self::Step) -> f64 {
        if step.cell_index == 0 {
            crate::sum(&[
                step.spatial.diffusion.population_roundoff,
                step.source_number_residual.abs(),
                step.capture_roundoff,
                step.fission_roundoff,
            ])
        } else {
            crate::sum(&[step.capture_roundoff, step.fission_roundoff])
        }
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
                "original neutron cell absent from accepted state",
            ))
    }
}
