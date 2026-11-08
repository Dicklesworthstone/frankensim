//! Local time refinement of the original nuclear/thermal transaction.
//!
//! A full native interval is compared with two native half intervals. Only the
//! fine trajectory is published. This estimates frozen-temperature/mean-flux
//! coupling and explicit heat-transport error; it does not change their order,
//! replace their physical laws, or certify a global error bound. Native event
//! and CFL cuts are discovered before nuclear evolution, so identical hidden
//! subdivisions cannot masquerade as a successful refinement comparison.

use std::borrow::Cow;

use super::*;

mod source;
pub use source::AdaptiveNuclearHeatSource;

/// Componentwise local time-error scales, independent of accounting roundoff.
/// Relative scales use changes from the accepted state, not absolute Kelvin or
/// an arbitrary caloric datum. Each absolute tolerance must be strictly positive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NuclearTimeTolerances {
    /// Original specific stored-energy coordinate, J/kg (h or u as declared).
    pub specific_energy_j_kg: f64,
    pub temperature_k: f64,
    /// Each isotope/fuel population and interval reaction count, expected nuclei.
    pub nuclei: f64,
    /// Each delayed-energy reservoir, cell-local released channel and transported
    /// radiation coordinate (receiver absorption, emitted/absorbed/escaped), J.
    pub nuclear_energy_j: f64,
    /// Each neutron/precursor amplitude, in the point model's normalization.
    pub neutron_population: f64,
    /// Spatial neutron and precursor inventories and transfers, actual counts.
    pub neutron_count: f64,
    /// Finite value in [0,1), applied separately to each physical coordinate.
    pub relative: f64,
}

impl Default for NuclearTimeTolerances {
    fn default() -> Self {
        Self {
            specific_energy_j_kg: 1e-3,
            temperature_k: 1e-3,
            nuclei: 1e-6,
            nuclear_energy_j: 1e-6,
            neutron_population: 1e-8,
            neutron_count: 1e-6,
            relative: 1e-3,
        }
    }
}

impl NuclearTimeTolerances {
    fn validate(self) -> Result<(), NuclearAdaptiveError> {
        if [
            self.specific_energy_j_kg,
            self.temperature_k,
            self.nuclei,
            self.nuclear_energy_j,
            self.neutron_population,
            self.neutron_count,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
            || !self.relative.is_finite()
            || !(0.0..1.0).contains(&self.relative)
        {
            return Err(NuclearAdaptiveError::Invalid(
                "positive finite physical tolerances and relative tolerance in [0,1) required",
            ));
        }
        Ok(())
    }
}

/// Native numerical gates and additional whole-call time-refinement resources.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NuclearAdaptiveControls {
    /// Original maximum step and whole-call thermal accounting allowance.
    pub step: StepConfig,
    pub initial_dt_s: f64,
    /// Floor on a full comparison span; accepted physical half steps are shorter.
    /// A smaller final/event/CFL-clipped span may be tried, but never waived.
    pub minimum_dt_s: f64,
    /// Includes rejected comparisons and attempts cut by a native event or CFL.
    /// Each attempt performs at most three native physical intervals.
    pub maximum_trials: usize,
    /// Accepted physical half steps, also capped by the original thermal budget.
    pub maximum_substeps: usize,
    pub tolerances: NuclearTimeTolerances,
}

impl Default for NuclearAdaptiveControls {
    fn default() -> Self {
        Self {
            step: StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-8,
            },
            initial_dt_s: 1.0,
            minimum_dt_s: 1e-12,
            maximum_trials: 65_536,
            maximum_substeps: 131_072,
            tolerances: NuclearTimeTolerances::default(),
        }
    }
}

/// Maximum raw coarse/fine differences by physical kind. `normalized` is the
/// first-order fine-error estimate divided by the corresponding local scale;
/// its unequal-half correction is applied by the shared adaptive driver.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NuclearTimeError {
    pub specific_energy_difference_j_kg: f64,
    pub temperature_difference_k: f64,
    pub nuclei_difference: f64,
    pub nuclear_energy_difference_j: f64,
    pub neutron_population_difference: f64,
    pub neutron_count_difference: f64,
    pub normalized: f64,
}

impl NuclearTimeError {
    fn merge(&mut self, other: Self) {
        self.specific_energy_difference_j_kg = self
            .specific_energy_difference_j_kg
            .max(other.specific_energy_difference_j_kg);
        self.temperature_difference_k = self
            .temperature_difference_k
            .max(other.temperature_difference_k);
        self.nuclei_difference = self.nuclei_difference.max(other.nuclei_difference);
        self.nuclear_energy_difference_j = self
            .nuclear_energy_difference_j
            .max(other.nuclear_energy_difference_j);
        self.neutron_population_difference = self
            .neutron_population_difference
            .max(other.neutron_population_difference);
        self.neutron_count_difference = self
            .neutron_count_difference
            .max(other.neutron_count_difference);
        self.normalized = self.normalized.max(other.normalized);
    }
}

/// Accepted heat accounts and actual refinement work. Nuclear accounts remain
/// on the returned original `NuclearHeating::audit()` over the whole call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NuclearAdaptiveReport {
    pub thermal: EnergyAudit,
    pub trials: usize,
    pub comparisons: usize,
    pub accepted_comparisons: usize,
    pub accuracy_rejections: usize,
    pub native_interval_cuts: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<NuclearTimeError>,
    /// Explicit controller state for the next call; no hidden PI history.
    pub suggested_dt_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum NuclearAdaptiveError {
    Source(FissionThermalError),
    Invalid(&'static str),
    TrialBudget,
    SubstepBudget,
    TimeResolution {
        time_s: f64,
        span_s: f64,
    },
    TimeAccuracy {
        time_s: f64,
        span_s: f64,
        error: NuclearTimeError,
    },
}

impl From<FissionThermalError> for NuclearAdaptiveError {
    fn from(error: FissionThermalError) -> Self {
        Self::Source(error)
    }
}
impl From<DepositionError> for NuclearAdaptiveError {
    fn from(error: DepositionError) -> Self {
        Self::Source(error.into())
    }
}
impl core::fmt::Display for NuclearAdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive nuclear/thermal advance refused: {self:?}")
    }
}
impl std::error::Error for NuclearAdaptiveError {}

impl<S: AdaptiveNuclearHeatSource> NuclearHeating<S> {
    /// Refine the existing nuclear/thermal evolution with local time estimates.
    ///
    /// Every attempted trajectory uses the original heat/source transaction,
    /// including material domains, authored transitions, conservation gates and
    /// cancellation. Physical/domain/budget failures propagate unchanged; only
    /// excessive estimated time error requests smaller intervals. A failure at
    /// any point publishes neither a partial prefix nor comparison histories.
    ///
    /// State and release comparisons are componentwise: a stable total or zero
    /// heat cannot conceal isotope redistribution or neutron/precursor error.
    /// Radiation compares actual absorbed energy at every receiver and the
    /// emitted/absorbed/escaped accounts of the native comparison intervals.
    /// A small temperature change therefore cannot hide a large absorption error.
    /// Native numerical budgets are distinct from these temporal estimates.
    pub fn advance_adaptive_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: NuclearAdaptiveControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, NuclearAdaptiveReport), NuclearAdaptiveError> {
        checkpoint()?;
        controls.tolerances.validate()?;
        if !target_s.is_finite()
            || target_s < thermal.time_s()
            || !controls.step.max_dt_s.is_finite()
            || controls.step.max_dt_s <= 0.0
            || !nonnegative(controls.step.energy_tolerance_j)
            || !controls.minimum_dt_s.is_finite()
            || controls.minimum_dt_s <= 0.0
            || controls.minimum_dt_s > controls.step.max_dt_s
            || !controls.initial_dt_s.is_finite()
            || controls.initial_dt_s < controls.minimum_dt_s
            || controls.maximum_trials == 0
        {
            return Err(NuclearAdaptiveError::Invalid(
                "invalid target, comparison span, numerical allowance or trial budget",
            ));
        }
        self.validate(thermal, self.time_s, &mut checkpoint)?;
        let start_s = self.time_s;
        let duration_s = target_s - start_s;
        let maximum_substeps = controls.maximum_substeps.min(thermal.budget().max_substeps);
        // Begin exactly once. Later comparisons inherit the accepted native
        // compensated nuclear account, including its whole-window allowances.
        let (mut current, mut heating, idle) =
            thermal.advance_observed_limited(self, start_s, controls.step, 0, &mut checkpoint)?;
        let mut accounts = ThermalAccount::default();
        accounts.record(thermal, &current, idle, &mut checkpoint)?;
        let mut span = controls.initial_dt_s.min(controls.step.max_dt_s);
        let mut report = NuclearAdaptiveReport {
            thermal: idle,
            trials: 0,
            comparisons: 0,
            accepted_comparisons: 0,
            accuracy_rejections: 0,
            native_interval_cuts: 0,
            maximum_accepted_normalized_error: 0.0,
            last_error: None,
            suggested_dt_s: span,
        };
        let mut allocated_thermal_tolerance = Sum::default();
        while current.time_s() < target_s {
            checkpoint()?;
            if maximum_substeps.saturating_sub(accounts.substeps) < 2 {
                return Err(NuclearAdaptiveError::SubstepBudget);
            }
            if report.trials >= controls.maximum_trials {
                return Err(NuclearAdaptiveError::TrialBudget);
            }
            report.trials += 1;
            let time_s = current.time_s();
            let end_s = if span >= target_s - time_s {
                target_s
            } else {
                finite(time_s + span)?
            };
            if end_s <= time_s {
                return Err(NuclearAdaptiveError::TimeResolution {
                    time_s,
                    span_s: span,
                });
            }
            let actual_span = end_s - time_s;
            let coarse_tolerance = thermal_allowance(
                controls.step.energy_tolerance_j,
                allocated_thermal_tolerance,
                actual_span,
                duration_s,
                maximum_substeps,
            )?;
            let coarse =
                single_interval(&heating, &current, end_s, coarse_tolerance, &mut checkpoint);
            let (coarse_thermal, coarse_heating, _) = match coarse {
                Ok(value) => value,
                Err(IntervalError::Cut(next_s)) => {
                    // Original CFL, source and native event owners selected this
                    // boundary before evolving a nucleus or depositing heat.
                    if next_s <= time_s || next_s >= end_s {
                        return Err(FissionThermalError::Arithmetic.into());
                    }
                    span = next_s - time_s;
                    report.native_interval_cuts += 1;
                    continue;
                }
                Err(IntervalError::Source(error)) => return Err(error.into()),
            };
            let middle_s = time_s + actual_span * 0.5;
            if middle_s <= time_s || middle_s >= end_s {
                return Err(NuclearAdaptiveError::TimeResolution {
                    time_s,
                    span_s: actual_span,
                });
            }
            let mut proposed_allocation = allocated_thermal_tolerance;
            let first_tolerance = thermal_allowance(
                controls.step.energy_tolerance_j,
                proposed_allocation,
                middle_s - time_s,
                duration_s,
                maximum_substeps,
            )?;
            proposed_allocation.add(first_tolerance)?;
            let (middle_thermal, middle_heating, first_audit) = single_interval(
                &heating,
                &current,
                middle_s,
                first_tolerance,
                &mut checkpoint,
            )
            .map_err(interval_failure)?;
            let second_tolerance = thermal_allowance(
                controls.step.energy_tolerance_j,
                proposed_allocation,
                end_s - middle_s,
                duration_s,
                maximum_substeps,
            )?;
            proposed_allocation.add(second_tolerance)?;
            let (fine_thermal, fine_heating, second_audit) = single_interval(
                &middle_heating,
                &middle_thermal,
                end_s,
                second_tolerance,
                &mut checkpoint,
            )
            .map_err(interval_failure)?;
            let mut error = compare(
                &current,
                &coarse_thermal,
                &fine_thermal,
                &heating,
                &coarse_heating,
                &middle_heating,
                &fine_heating,
                actual_span,
                controls.tolerances,
                &mut checkpoint,
            )?;
            // For a first-order endpoint map, the fine leading error has
            // coefficient h1²+h2². The coarse/fine difference has h²-h1²-h2².
            // This also covers unequal representable halves at a large clock.
            let a = (middle_s - time_s) / actual_span;
            let b = (end_s - middle_s) / actual_span;
            let q = a * a + b * b;
            if !q.is_finite() || q >= 1.0 {
                return Err(NuclearAdaptiveError::TimeResolution {
                    time_s,
                    span_s: actual_span,
                });
            }
            error.normalized *= q / (1.0 - q);
            report.comparisons += 1;
            report.last_error = Some(error);
            if error.normalized > 1.0 {
                report.accuracy_rejections += 1;
                let smaller = actual_span * 0.5;
                if smaller < controls.minimum_dt_s || time_s + smaller <= time_s {
                    return Err(NuclearAdaptiveError::TimeAccuracy {
                        time_s,
                        span_s: actual_span,
                        error,
                    });
                }
                span = smaller;
                continue;
            }
            checkpoint()?;
            accounts.record(&current, &middle_thermal, first_audit, &mut checkpoint)?;
            accounts.record(
                &middle_thermal,
                &fine_thermal,
                second_audit,
                &mut checkpoint,
            )?;
            allocated_thermal_tolerance = proposed_allocation;
            current = fine_thermal;
            heating = fine_heating;
            report.accepted_comparisons += 1;
            report.maximum_accepted_normalized_error = report
                .maximum_accepted_normalized_error
                .max(error.normalized);
            // Dyadic changes are deterministic and let constant operators reuse
            // their original preparation caches across accepted half intervals.
            span = if error.normalized <= 0.25 {
                (actual_span * 2.0).min(controls.step.max_dt_s)
            } else {
                actual_span
            };
            span = span.max(controls.minimum_dt_s);
        }
        report.thermal = accounts.finish(
            thermal,
            &current,
            controls.step.energy_tolerance_j,
            &mut checkpoint,
        )?;
        report.suggested_dt_s = span;
        heating.validate(&current, heating.time_s, &mut checkpoint)?;
        heating.admit_audit(&heating.audit)?;
        checkpoint()?;
        Ok((current, heating, report))
    }
}

// The original explicit owner alone chooses events and conductive CFL. Reject a
// hidden subdivision before source evolution, then compare exactly that span.
#[derive(Clone)]
struct SingleInterval<'s, S: NuclearHeatSource> {
    source: Cow<'s, NuclearHeating<S>>,
    end_s: f64,
}

enum IntervalError {
    Source(FissionThermalError),
    Cut(f64),
}
impl From<DepositionError> for IntervalError {
    fn from(error: DepositionError) -> Self {
        Self::Source(error.into())
    }
}
impl From<FissionThermalError> for IntervalError {
    fn from(error: FissionThermalError) -> Self {
        Self::Source(error)
    }
}

impl<S: NuclearHeatSource> DepositionObserver for SingleInterval<'_, S> {
    type Error = IntervalError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.source.to_mut().begin().map_err(Into::into)
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.source
            .to_mut()
            .next_transition_time_s(thermal, checkpoint)
            .map_err(Into::into)
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        if end_s != self.end_s {
            return Err(IntervalError::Cut(end_s));
        }
        self.source
            .to_mut()
            .step_sources(thermal, start_s, duration_s, end_s, checkpoint)
            .map_err(Into::into)
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.source
            .to_mut()
            .event(event, thermal)
            .map_err(Into::into)
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.source
            .to_mut()
            .step(start_s, duration_s, old_temperatures_k, thermal, checkpoint)
            .map_err(Into::into)
    }
}

fn single_interval<'m, S: NuclearHeatSource>(
    source: &NuclearHeating<S>,
    thermal: &DepositionSimulation<'m>,
    end_s: f64,
    energy_tolerance_j: f64,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(DepositionSimulation<'m>, NuclearHeating<S>, EnergyAudit), IntervalError> {
    let observer = SingleInterval {
        source: Cow::Borrowed(source),
        end_s,
    };
    let (thermal, observer, audit) = thermal.advance_observed_continuing_limited(
        &observer,
        end_s,
        StepConfig {
            max_dt_s: end_s - thermal.time_s(),
            energy_tolerance_j,
        },
        1,
        &mut *checkpoint,
    )?;
    let source = observer.source.into_owned();
    source.validate(&thermal, source.time_s, checkpoint)?;
    source.admit_audit(&source.audit)?;
    Ok((thermal, source, audit))
}

fn interval_failure(error: IntervalError) -> NuclearAdaptiveError {
    match error {
        IntervalError::Source(error) => error.into(),
        // After the coarse interval has been admitted, its shorter halves must
        // not acquire another cut under these same fixed native coefficients.
        IntervalError::Cut(_) => NuclearAdaptiveError::Invalid(
            "native interval changed after coarse schedule/CFL admission",
        ),
    }
}

fn thermal_allowance(
    total: f64,
    allocated: Sum,
    span_s: f64,
    duration_s: f64,
    maximum_substeps: usize,
) -> Result<f64, FissionThermalError> {
    let remaining = (total - allocated.nonnegative_upper_bound()?).max(0.0);
    // Half the allowance follows physical duration; half is shared by the
    // bounded number of accepted native steps. The latter admits representable
    // storage roundoff even at a very short final/event-clipped interval.
    let per_step = (total * 0.5) / (maximum_substeps as f64).next_up();
    let proportional = total * 0.5 * (span_s / duration_s) + per_step;
    // Reserve numerical allowance only for accepted halves. The cap prevents
    // repeated calls from renewing the native unrepresented-heat allowance.
    finite(proportional.min(remaining))
}

#[allow(clippy::too_many_arguments)]
fn compare<S: AdaptiveNuclearHeatSource>(
    initial: &DepositionSimulation<'_>,
    coarse: &DepositionSimulation<'_>,
    fine: &DepositionSimulation<'_>,
    initial_source: &NuclearHeating<S>,
    coarse_source: &NuclearHeating<S>,
    middle_source: &NuclearHeating<S>,
    fine_source: &NuclearHeating<S>,
    duration_s: f64,
    tolerances: NuclearTimeTolerances,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<NuclearTimeError, FissionThermalError> {
    if coarse.cells().len() != fine.cells().len() {
        return Err(FissionThermalError::Invalid(
            "refinement changed cell identities",
        ));
    }
    let mut metric = source::Comparison::new(tolerances, checkpoint);
    for (&id, cell) in fine.cells() {
        let other = coarse
            .cells()
            .get(&id)
            .ok_or(FissionThermalError::InactiveCell(id))?;
        if cell.mass_kg != other.mass_kg || cell.material_id != other.material_id {
            return Err(FissionThermalError::BindingChanged(id));
        }
        // A newly born cell is activated exactly at this comparison's endpoint;
        // it has not aged on either path and supplies its own initial coordinate.
        let initial_h = initial
            .cells()
            .get(&id)
            .map_or(cell.specific_h_j_kg, |c| c.specific_h_j_kg);
        metric.specific_energy(initial_h, other.specific_h_j_kg, cell.specific_h_j_kg)?;
        let initial_t = if initial.cells().contains_key(&id) {
            initial.temperature_k(id)?
        } else {
            fine.temperature_k(id)?
        };
        metric.temperature(
            initial_t,
            coarse.temperature_k(id)?,
            fine.temperature_k(id)?,
        )?;
    }
    // single_interval admits exactly one native interval; event/CFL cuts are
    // retried before source evolution. The two accepted half-interval results
    // add once. Their continuing cumulative accounts already share a prefix
    // and must never be added or subtracted to recover these small releases.
    match (
        coarse_source.last_radiation_deposition(),
        middle_source.last_radiation_deposition(),
        fine_source.last_radiation_deposition(),
    ) {
        (None, None, None) => {}
        (Some(coarse), Some(first), Some(second)) => {
            metric.radiation(coarse, first, second)?;
        }
        _ => {
            return Err(FissionThermalError::Invalid(
                "refinement changed radiation response coverage",
            ));
        }
    }
    let mut error = metric.error;
    for (&id, bound) in &initial_source.sources {
        checkpoint()?;
        let coarse = &coarse_source.sources[&id];
        let middle = &middle_source.sources[&id];
        let fine = &fine_source.sources[&id];
        error.merge(bound.source.time_state_error(
            &coarse.source,
            &fine.source,
            tolerances,
            checkpoint,
        )?);
        error.merge(S::time_step_error(
            coarse
                .last_step
                .as_ref()
                .ok_or(FissionThermalError::Arithmetic)?,
            middle
                .last_step
                .as_ref()
                .ok_or(FissionThermalError::Arithmetic)?,
            fine.last_step
                .as_ref()
                .ok_or(FissionThermalError::Arithmetic)?,
            duration_s,
            tolerances,
            checkpoint,
        )?);
    }
    Ok(error)
}

#[derive(Default)]
struct ThermalAccount {
    external: Sum,
    incoming_mass: Sum,
    incoming_energy: Sum,
    absolute_residual: Sum,
    births: BTreeMap<CellId, f64>,
    substeps: usize,
    events: usize,
}

impl ThermalAccount {
    fn record(
        &mut self,
        before: &DepositionSimulation<'_>,
        after: &DepositionSimulation<'_>,
        audit: EnergyAudit,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), FissionThermalError> {
        self.external.add(audit.external_heat_j)?;
        self.incoming_mass.add(audit.incoming_mass_kg)?;
        self.incoming_energy.add(audit.incoming_energy_j)?;
        self.absolute_residual.add(audit.residual_j.abs())?;
        self.substeps = self
            .substeps
            .checked_add(audit.substeps)
            .ok_or(FissionThermalError::Arithmetic)?;
        self.events = self
            .events
            .checked_add(audit.events)
            .ok_or(FissionThermalError::Arithmetic)?;
        for (index, (&id, cell)) in after.cells().iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if !before.cells().contains_key(&id) {
                self.births.insert(id, cell.specific_h_j_kg);
            }
        }
        Ok(())
    }

    fn finish(
        &self,
        initial: &DepositionSimulation<'_>,
        final_state: &DepositionSimulation<'_>,
        tolerance_j: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<EnergyAudit, FissionThermalError> {
        let mut thermal_change = Sum::default();
        let mut stored_change = Sum::default();
        for (index, (&id, cell)) in final_state.cells().iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let old = initial.cells().get(&id);
            let baseline = old
                .map(|c| c.specific_h_j_kg)
                .or_else(|| self.births.get(&id).copied())
                .ok_or(FissionThermalError::Arithmetic)?;
            let change = finite(cell.mass_kg * finite(cell.specific_h_j_kg - baseline)?)?;
            thermal_change.add(change)?;
            stored_change.add(if old.is_some() {
                change
            } else {
                finite(cell.mass_kg * cell.specific_h_j_kg)?
            })?;
        }
        // Initial/final storage is evaluated once, relative to each actual birth
        // or initial coordinate. Do not subtract rounded total formation energies.
        thermal_change.add(-self.external.sum)?;
        thermal_change.add(-self.external.correction)?;
        let residual_j = thermal_change.value()?;
        if residual_j.abs() > tolerance_j
            || self.absolute_residual.nonnegative_upper_bound()? > tolerance_j
        {
            return Err(DepositionError::EnergyBalance(
                residual_j
                    .abs()
                    .max(self.absolute_residual.nonnegative_upper_bound()?),
            )
            .into());
        }
        Ok(EnergyAudit {
            incoming_mass_kg: self.incoming_mass.value()?,
            incoming_energy_j: self.incoming_energy.value()?,
            external_heat_j: self.external.value()?,
            stored_energy_change_j: stored_change.value()?,
            residual_j,
            substeps: self.substeps,
            events: self.events,
        })
    }
}
