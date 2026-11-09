//! Time-resolved fusion forcing on the original neutron transport owner.
//!
//! Every retained interval preserves its original channel moments. Its source
//! is uniform only within that accepted interval. Different cells may emit
//! concurrently; intervals within one cell must not overlap. Preparing a sweep
//! of source changes avoids rescanning the burn history during transport.
//! Lookup is logarithmic and returns a borrowed original-order rate vector.
//! The neutron propagator, cross sections and spatial approximation are unchanged.

use std::sync::Arc;

use fs_qty::Time;
use fs_thermochem::fusion::neutrons::FusionNeutronIntervalV1;

use super::{
    FusionNeutronSource, NeutronKineticBudget, NeutronKineticError, NeutronKineticModel,
    NeutronKineticState, checked,
};
use crate::{Accumulator, EnergyChannels, FissionError, nonnegative};

/// One original burn interval and the neutron cell receiving its emitted source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FusionNeutronHistoryEntry {
    pub cell: usize,
    pub interval: FusionNeutronIntervalV1,
}

/// Construction work and storage bounds, independent of physical uncertainty.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FusionNeutronHistoryLimits {
    pub maximum_intervals: usize,
    /// Scalar rate slots in retained changes, zero rates and temporary projections.
    pub maximum_rate_values: usize,
}

impl Default for FusionNeutronHistoryLimits {
    fn default() -> Self {
        Self {
            maximum_intervals: 16_384,
            maximum_rate_values: 1_048_576,
        }
    }
}

#[derive(Clone, Debug)]
struct Change {
    time_s: f64,
    rates_s: Vec<f64>,
    original_particle_rate_s: f64,
    original_kinetic_power_w: f64,
}

#[derive(Clone, Copy)]
struct Event {
    time_s: f64,
    entry: usize,
    start: bool,
}

/// Immutable spatial source history. All times use the same supplied clock.
/// No prior emission is replayed when a caller resumes from a later checkpoint.
#[derive(Clone, Debug)]
pub struct FusionNeutronHistory {
    model: Arc<NeutronKineticModel>,
    entries: Arc<[FusionNeutronHistoryEntry]>,
    changes: Arc<[Change]>,
    zero_rates_s: Arc<[f64]>,
}

impl FusionNeutronHistory {
    pub fn project(
        model: Arc<NeutronKineticModel>,
        entries: Vec<FusionNeutronHistoryEntry>,
        limits: FusionNeutronHistoryLimits,
    ) -> Result<Self, NeutronKineticError> {
        Self::project_with_checkpoint(model, entries, limits, || Ok(()))
    }

    /// Retain the original intervals and project each through the existing
    /// count/energy-conserving mean-energy source owner. No spectrum is inferred.
    pub fn project_with_checkpoint(
        model: Arc<NeutronKineticModel>,
        entries: Vec<FusionNeutronHistoryEntry>,
        limits: FusionNeutronHistoryLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, NeutronKineticError> {
        checkpoint()?;
        let cells = model.diffusion().cell_count();
        let groups = model.group_energies_j().len();
        let nodes = cells * groups;
        if limits.maximum_intervals > 1_000_000 || entries.len() > limits.maximum_intervals {
            return Err(FissionError::Budget("fusion history interval count").into());
        }
        let projected_values = entries
            .len()
            .checked_mul(groups)
            .and_then(|n| n.checked_add(nodes))
            .ok_or(FissionError::Budget("fusion history rate storage"))?;
        if projected_values > limits.maximum_rate_values {
            return Err(FissionError::Budget("fusion history rate storage").into());
        }
        let mut order = Vec::new();
        order
            .try_reserve_exact(entries.len())
            .map_err(|_| FissionError::Budget("fusion history ordering"))?;
        for (index, entry) in entries.iter().enumerate() {
            checkpoint()?;
            if entry.cell >= cells {
                return Err(NeutronKineticError::Invalid(
                    "fusion history cell outside original geometry",
                ));
            }
            order.push(index);
        }
        order.sort_unstable_by(|&a, &b| {
            entries[a]
                .cell
                .cmp(&entries[b].cell)
                .then(
                    entries[a]
                        .interval
                        .start()
                        .value()
                        .total_cmp(&entries[b].interval.start().value()),
                )
                .then(a.cmp(&b))
        });
        for pair in order.windows(2) {
            checkpoint()?;
            let (a, b) = (entries[pair[0]], entries[pair[1]]);
            if a.cell == b.cell && b.interval.start().value() < a.interval.end().value() {
                return Err(NeutronKineticError::Invalid(
                    "overlapping fusion intervals in one cell",
                ));
            }
        }
        let mut projections = Vec::new();
        projections
            .try_reserve_exact(entries.len())
            .map_err(|_| FissionError::Budget("fusion history projections"))?;
        let mut events = Vec::new();
        let event_count = entries
            .len()
            .checked_mul(2)
            .ok_or(FissionError::Budget("fusion history events"))?;
        events
            .try_reserve_exact(event_count)
            .map_err(|_| FissionError::Budget("fusion history events"))?;
        for (index, entry) in entries.iter().enumerate() {
            checkpoint()?;
            let emission = *entry.interval.emission();
            let projected = FusionNeutronSource::project_interval_mean(
                Arc::clone(&model),
                emission,
                entry.cell,
            )?;
            let duration = emission.duration().value();
            let particle_rate = checked(emission.total_particle_count() / duration)?;
            let power = checked(emission.total_kinetic_energy().value() / duration)?;
            if (emission.total_particle_count() > 0.0 && particle_rate == 0.0)
                || (emission.total_kinetic_energy().value() > 0.0 && power == 0.0)
            {
                return Err(NeutronKineticError::Arithmetic(
                    "fusion history moment rate underflow",
                ));
            }
            projections.push((
                projected.neutron_rates_s()[entry.cell * groups..(entry.cell + 1) * groups]
                    .to_vec(),
                particle_rate,
                power,
            ));
            events.push(Event {
                time_s: entry.interval.start().value(),
                entry: index,
                start: true,
            });
            events.push(Event {
                time_s: entry.interval.end().value(),
                entry: index,
                start: false,
            });
        }
        // End events precede starts at a common boundary. Replacing each cell's
        // complete row avoids cancellation from subtracting a large old source.
        events.sort_unstable_by(|a, b| {
            a.time_s
                .total_cmp(&b.time_s)
                .then(a.start.cmp(&b.start))
                .then(a.entry.cmp(&b.entry))
        });
        let mut active = vec![None; cells];
        let mut rates = vec![0.0; nodes];
        let mut changes = Vec::new();
        let mut cursor = 0;
        while cursor < events.len() {
            checkpoint()?;
            let time_s = events[cursor].time_s;
            while cursor < events.len() && events[cursor].time_s == time_s {
                let event = events[cursor];
                let cell = entries[event.entry].cell;
                let row = &mut rates[cell * groups..(cell + 1) * groups];
                if event.start {
                    active[cell] = Some(event.entry);
                    row.copy_from_slice(&projections[event.entry].0);
                } else {
                    active[cell] = None;
                    row.fill(0.0);
                }
                cursor += 1;
            }
            let mut number = Accumulator::new(0.0);
            let mut energy = Accumulator::new(0.0);
            for &index in active.iter().flatten() {
                number.add(projections[index].1)?;
                energy.add(projections[index].2)?;
            }
            let retained_values = (changes.len() + 1)
                .checked_mul(nodes)
                .and_then(|n| n.checked_add(projected_values))
                .ok_or(FissionError::Budget("fusion history rate storage"))?;
            if retained_values > limits.maximum_rate_values {
                return Err(FissionError::Budget("fusion history rate storage").into());
            }
            changes
                .try_reserve(1)
                .map_err(|_| FissionError::Budget("fusion history changes"))?;
            changes.push(Change {
                time_s,
                rates_s: rates.clone(),
                original_particle_rate_s: number.value(),
                original_kinetic_power_w: energy.value(),
            });
        }
        checkpoint()?;
        Ok(Self {
            model,
            entries: entries.into(),
            changes: changes.into(),
            // Every interval has ended at the final event. Reuse that zero row
            // so preparation stays within its declared scalar-rate storage.
            zero_rates_s: rates.into(),
        })
    }

    pub fn model(&self) -> &Arc<NeutronKineticModel> {
        &self.model
    }
    pub fn entries(&self) -> &[FusionNeutronHistoryEntry] {
        &self.entries
    }

    fn change_at(&self, time_s: f64) -> Option<&Change> {
        let next = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        next.checked_sub(1).map(|index| &self.changes[index])
    }

    pub fn rates_at(&self, time: Time) -> Result<&[f64], NeutronKineticError> {
        if !nonnegative(time.value()) {
            return Err(NeutronKineticError::Invalid(
                "finite nonnegative fusion history time required",
            ));
        }
        Ok(self
            .change_at(time.value())
            .map_or(self.zero_rates_s.as_ref(), |change| {
                change.rates_s.as_slice()
            }))
    }

    /// First source transition strictly after the supplied time.
    pub fn next_after(&self, time: Time) -> Option<Time> {
        if !nonnegative(time.value()) {
            return None;
        }
        self.changes
            .get(
                self.changes
                    .partition_point(|change| change.time_s <= time.value()),
            )
            .map(|change| Time::new(change.time_s))
    }

    /// Original expected particle and kinetic-energy input over one constant
    /// source interval. Crossing a source transition is an explicit refusal.
    pub fn interval_moments(
        &self,
        start: Time,
        duration: Time,
    ) -> Result<(f64, f64), NeutronKineticError> {
        self.rates_at(start)?;
        if !nonnegative(duration.value())
            || !(start.value() + duration.value()).is_finite()
            || (duration.value() > 0.0 && start.value() + duration.value() <= start.value())
            || self
                .next_after(start)
                .is_some_and(|next| duration.value() > next.value() - start.value())
        {
            return Err(NeutronKineticError::Invalid(
                "finite fusion source interval without internal transition required",
            ));
        }
        let Some(change) = self.change_at(start.value()) else {
            return Ok((0.0, 0.0));
        };
        Ok((
            checked(change.original_particle_rate_s * duration.value())?,
            checked(change.original_kinetic_power_w * duration.value())?,
        ))
    }
}

/// Whole-request account limits and a bound on source-aligned neutron advances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NeutronKineticHistoryBudget {
    pub transport: NeutronKineticBudget,
    /// Finite positive transport interval ceiling [s], independent of burn
    /// recording intervals. Resolve fast interaction/residence times here.
    pub maximum_step_s: f64,
    pub maximum_substeps: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct NeutronKineticHistoryStep {
    pub substeps: usize,
    pub initial_kinetic_energy_j: f64,
    pub final_kinetic_energy_j: f64,
    pub initial_precursor_energy_j: f64,
    pub final_precursor_energy_j: f64,
    /// Reserved energy internally converted into delayed-neutron energy [J].
    pub precursor_conversion_energy_j: f64,
    /// Reserved energy exported by original material outflows, separate from
    /// released heat, radiation and neutron leakage [J].
    pub cell_exported_precursor_energy_j: Vec<f64>,
    pub exported_precursor_energy_j: f64,
    pub injected_neutrons: f64,
    pub injected_kinetic_energy_j: f64,
    /// Signed prompt Q per original material cell over this request [J].
    pub cell_nuclear_reaction_energy_j: Vec<f64>,
    pub nuclear_reaction_energy_j: f64,
    /// Additional future-neutron conversion reserved by induced fission [J].
    /// It supplements prompt Q and is not external source energy or local heat.
    pub cell_induced_precursor_energy_j: Vec<f64>,
    pub induced_precursor_energy_j: f64,
    pub cell_released: Vec<EnergyChannels>,
    pub released: EnergyChannels,
    pub fusion_number_residual: f64,
    pub fusion_energy_residual_j: f64,
    pub energy_balance_residual_j: f64,
    pub energy_roundoff_j: f64,
    pub population_roundoff: f64,
}

impl NeutronKineticState {
    /// Evolve from the supplied checkpoint clock, clipping at retained source
    /// transitions and the admitted transport step ceiling. Input state remains
    /// unchanged on any later refusal.
    /// This does not refine the original burn history or its within-bin spectrum.
    pub fn trial_with_history(
        &self,
        start: Time,
        end: Time,
        history: &FusionNeutronHistory,
        budget: NeutronKineticHistoryBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, NeutronKineticHistoryStep), NeutronKineticError> {
        checkpoint()?;
        if !nonnegative(start.value())
            || !nonnegative(end.value())
            || end.value() < start.value()
            || !crate::positive(budget.maximum_step_s)
            || budget.maximum_substeps == 0
            || budget.maximum_substeps > 1_000_000
            || **self.model() != **history.model()
        {
            return Err(NeutronKineticError::Invalid(
                "matching fusion history, ordered finite times and bounded work required",
            ));
        }
        // An idle call still admits the original count/energy/work budget.
        self.trial_with_rates(
            0.0,
            history.rates_at(start)?,
            budget.transport,
            &mut checkpoint,
        )?;
        let initial = self.kinetic_energy_j()?;
        let initial_precursor = self.precursor_energy_j()?;
        let mut candidate = self.clone();
        let mut time = start.value();
        let mut substeps = 0;
        let mut injected_number = Accumulator::new(0.0);
        let mut injected_energy = Accumulator::new(0.0);
        let mut number_residual = Accumulator::new(0.0);
        let mut energy_residual = Accumulator::new(0.0);
        let mut energy_defect = Accumulator::new(0.0);
        let mut count_defect = Accumulator::new(0.0);
        let mut cells = vec![[Accumulator::new(0.0); 3]; self.model().diffusion().cell_count()];
        let mut nuclear_cells = vec![Accumulator::new(0.0); cells.len()];
        let mut induced_cells = vec![Accumulator::new(0.0); cells.len()];
        let mut precursor_cells = vec![Accumulator::new(0.0); cells.len()];
        let mut precursor_conversion = Accumulator::new(0.0);
        while time < end.value() {
            checkpoint()?;
            if substeps >= budget.maximum_substeps {
                return Err(FissionError::Budget("fusion history neutron advances").into());
            }
            let stop = history
                .next_after(Time::new(time))
                .map_or(end.value(), |next| next.value().min(end.value()))
                .min(time + budget.maximum_step_s);
            if stop <= time {
                return Err(NeutronKineticError::Arithmetic(
                    "fusion history transport step cannot advance the supplied clock",
                ));
            }
            let duration = stop - time;
            let (expected_number, expected_energy) =
                history.interval_moments(Time::new(time), Time::new(duration))?;
            let (next, report) = candidate.trial_with_rates(
                duration,
                history.rates_at(Time::new(time))?,
                budget.transport,
                &mut checkpoint,
            )?;
            let nr = checked(report.diffusion.injected_neutrons - expected_number)?;
            let er = checked(report.injected_kinetic_energy_j - expected_energy)?;
            count_defect.add(report.diffusion.population_roundoff)?;
            count_defect.add(nr.abs())?;
            energy_defect.add(report.energy_roundoff_j)?;
            energy_defect.add(
                precursor_conversion
                    .add(report.precursor_conversion_energy_j)?
                    .abs(),
            )?;
            for (target, &value) in precursor_cells
                .iter_mut()
                .zip(&report.cell_exported_precursor_energy_j)
            {
                energy_defect.add(target.add(value)?.abs())?;
            }
            energy_defect.add(er.abs())?;
            energy_defect.add(injected_energy.add(report.injected_kinetic_energy_j)?.abs())?;
            count_defect.add(
                injected_number
                    .add(report.diffusion.injected_neutrons)?
                    .abs(),
            )?;
            number_residual.add(nr)?;
            energy_residual.add(er)?;
            for (target, value) in nuclear_cells
                .iter_mut()
                .zip(&report.cell_nuclear_reaction_energy_j)
            {
                energy_defect.add(target.add(*value)?.abs())?;
            }
            for (target, &value) in induced_cells
                .iter_mut()
                .zip(&report.cell_induced_precursor_energy_j)
            {
                energy_defect.add(target.add(value)?.abs())?;
            }
            for (target, released) in cells.iter_mut().zip(report.cell_released) {
                for (sum, value) in target.iter_mut().zip([
                    released.deposited_j,
                    released.radiation_j,
                    released.escaped_j,
                ]) {
                    energy_defect.add(sum.add(value)?.abs())?;
                }
            }
            if energy_defect.value() > budget.transport.energy_tolerance_j
                || count_defect.value() > budget.transport.diffusion.population_tolerance
            {
                return Err(NeutronKineticError::Balance(
                    "accumulated fusion history account exceeds allowance",
                ));
            }
            candidate = next;
            substeps += 1;
            time = stop;
        }
        let mut totals = [Accumulator::new(0.0); 3];
        let mut cell_released = Vec::with_capacity(cells.len());
        for cell in cells {
            for (total, value) in totals.iter_mut().zip(cell) {
                energy_defect.add(export_defect(value))?;
                energy_defect.add(total.add(value.value())?.abs())?;
            }
            cell_released.push(crate::channel_values(cell));
        }
        for total in totals {
            energy_defect.add(export_defect(total))?;
        }
        energy_defect.add(export_defect(injected_energy))?;
        count_defect.add(export_defect(injected_number))?;
        let mut nuclear_total = Accumulator::new(0.0);
        for &cell in &nuclear_cells {
            energy_defect.add(export_defect(cell))?;
            energy_defect.add(nuclear_total.add(cell.value())?.abs())?;
        }
        energy_defect.add(export_defect(nuclear_total))?;
        let mut induced_total = Accumulator::new(0.0);
        for &cell in &induced_cells {
            energy_defect.add(export_defect(cell))?;
            energy_defect.add(induced_total.add(cell.value())?.abs())?;
        }
        energy_defect.add(export_defect(induced_total))?;
        let mut precursor_exported = Accumulator::new(0.0);
        for &cell in &precursor_cells {
            energy_defect.add(export_defect(cell))?;
            energy_defect.add(precursor_exported.add(cell.value())?.abs())?;
        }
        energy_defect.add(export_defect(precursor_exported))?;
        energy_defect.add(export_defect(precursor_conversion))?;
        let released = crate::channel_values(totals);
        let final_energy = candidate.kinetic_energy_j()?;
        let final_precursor = candidate.precursor_energy_j()?;
        let residual = checked(crate::sum(&[
            final_energy,
            final_precursor,
            -initial,
            -initial_precursor,
            released.deposited_j,
            released.radiation_j,
            released.escaped_j,
            precursor_exported.value(),
            -injected_energy.value(),
            -nuclear_total.value(),
            -induced_total.value(),
        ]))?;
        energy_defect.add(residual.abs())?;
        if energy_defect.value() > budget.transport.energy_tolerance_j
            || count_defect.value() > budget.transport.diffusion.population_tolerance
        {
            return Err(NeutronKineticError::Balance(
                "fusion history export account exceeds allowance",
            ));
        }
        checkpoint()?;
        Ok((
            candidate,
            NeutronKineticHistoryStep {
                substeps,
                initial_kinetic_energy_j: initial,
                final_kinetic_energy_j: final_energy,
                initial_precursor_energy_j: initial_precursor,
                final_precursor_energy_j: final_precursor,
                precursor_conversion_energy_j: precursor_conversion.value(),
                cell_exported_precursor_energy_j: precursor_cells
                    .into_iter()
                    .map(Accumulator::value)
                    .collect(),
                exported_precursor_energy_j: precursor_exported.value(),
                injected_neutrons: injected_number.value(),
                injected_kinetic_energy_j: injected_energy.value(),
                cell_nuclear_reaction_energy_j: nuclear_cells
                    .into_iter()
                    .map(Accumulator::value)
                    .collect(),
                nuclear_reaction_energy_j: nuclear_total.value(),
                cell_induced_precursor_energy_j: induced_cells
                    .into_iter()
                    .map(Accumulator::value)
                    .collect(),
                induced_precursor_energy_j: induced_total.value(),
                cell_released,
                released,
                fusion_number_residual: number_residual.value(),
                fusion_energy_residual_j: energy_residual.value(),
                energy_balance_residual_j: residual,
                energy_roundoff_j: energy_defect.value(),
                population_roundoff: count_defect.value(),
            },
        ))
    }
}

fn export_defect(value: Accumulator) -> f64 {
    crate::sum(&[value.hi, value.lo, -value.value()]).abs()
}

#[cfg(test)]
mod tests;
