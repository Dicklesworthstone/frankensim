//! Physical time accuracy for simultaneous heat and liquid/vapor exchange.
//!
//! Compare one backward-Euler interval with two actual half intervals. Only
//! the conservative fine trajectory is accepted. The original implicit heat
//! controller supplies order correction and bounded timestep factors; original
//! calorics and each interface's actual mass supply the error coordinates.

use super::super::account::TimelineAccount;
use super::super::{CoupledEvaporationAudit, LiquidBranch, ThermalEvaporationError};
use super::{
    EvaporatingCellBinding, EvaporationImplicitError, EvaporationNetworkImplicitControls,
    EvaporationNetworkImplicitReport,
};
use crate::enthalpy::backward_euler::BoundaryExchangeStepCoefficients;
use crate::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;
use crate::enthalpy::observe::EnergySum;
use crate::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, Event, finite, positive,
};
use crate::evaporation::{EvaporationLaw, LiquidSurfaceState, SurfaceStep};
use std::collections::{BTreeMap, BTreeSet};

/// Time-error allowances are distinct from nonlinear energy/mass residual gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveEvaporationControls {
    /// Whole-request numerical work and physical balance tolerances.
    pub solve: EvaporationNetworkImplicitControls,
    pub initial_dt_s: f64,
    /// Bounds on the full comparison span, not its two accepted half intervals.
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    pub maximum_comparisons: usize,
    pub maximum_reductions: usize,
    /// Original specific stored energy, including latent storage, for EVERY cell.
    pub enthalpy_absolute_tolerance_j_kg: f64,
    pub temperature_absolute_tolerance_k: f64,
    /// Actual inventory of EACH selected liquid, not net network mass.
    pub mass_absolute_tolerance_kg: f64,
    /// Scaled by changes from the accepted start, never arbitrary h or Kelvin datums.
    pub relative_tolerance: f64,
}
impl AdaptiveEvaporationControls {
    pub fn new(solve: EvaporationNetworkImplicitControls) -> Self {
        Self {
            solve,
            initial_dt_s: 1.0,
            minimum_dt_s: 1e-9,
            maximum_dt_s: 1.0,
            maximum_comparisons: 4096,
            maximum_reductions: 32,
            enthalpy_absolute_tolerance_j_kg: 1.0,
            temperature_absolute_tolerance_k: 1e-3,
            mass_absolute_tolerance_kg: 1e-9,
            relative_tolerance: 1e-3,
        }
    }
    fn validate(self) -> Result<(), AdaptiveEvaporationError> {
        self.solve.validate()?;
        if !positive(self.initial_dt_s)
            || !positive(self.minimum_dt_s)
            || !positive(self.maximum_dt_s)
            || self.minimum_dt_s > self.maximum_dt_s
            || self.initial_dt_s < self.minimum_dt_s
            || self.maximum_comparisons == 0
            || !positive(self.enthalpy_absolute_tolerance_j_kg)
            || !positive(self.temperature_absolute_tolerance_k)
            || !positive(self.mass_absolute_tolerance_kg)
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
        {
            return Err(AdaptiveEvaporationError::Invalid(
                "invalid adaptive evaporation controls",
            ));
        }
        Ok(())
    }
}

/// Coarse/fine discrepancies and their corrected local-error estimate.
/// This is not a certified bound on accumulated trajectory error.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EvaporationTimeError {
    pub maximum_enthalpy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_mass_difference_kg: f64,
    pub normalized: f64,
}

/// Signed, compensated phase transfers on one liquid inventory. Positive mass
/// and recession mean evaporation; condensation retains its negative sign.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EvaporationSurfaceTotals {
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    pub vapor_energy_j: f64,
    pub recession_m: f64,
    pub intervals: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveEvaporationReport {
    /// Accepted fine intervals and original events only. `phase_steps` counts
    /// individual interface applications; `last_step` is the last canonical id.
    /// Aggregate recession sums signed distances across interfaces; use
    /// `surface_totals[id].recession_m` for the distance at an individual cell.
    pub physical: CoupledEvaporationAudit,
    /// Accepted fine transfers per cell; speculative work never enters these sums.
    pub surface_totals: BTreeMap<CellId, EvaporationSurfaceTotals>,
    /// Complete per-interface receipts of the last accepted half interval.
    pub last_surfaces: BTreeMap<CellId, SurfaceStep>,
    /// Includes coarse trajectories and accuracy-rejected comparisons.
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    pub maximum_accepted_residual_norm_j: f64,
    pub absolute_energy_residual_j: f64,
    pub absolute_mass_rate_residual_kg: f64,
    pub comparisons: usize,
    pub accepted_comparisons: usize,
    pub accuracy_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<EvaporationTimeError>,
    /// Explicit next-step suggestion; no hidden controller history is retained.
    pub suggested_dt_s: f64,
}

#[derive(Debug)]
pub enum AdaptiveEvaporationError {
    Implicit(EvaporationImplicitError),
    Invalid(&'static str),
    ComparisonBudget,
    TimeAccuracy {
        time_s: f64,
        span_s: f64,
        error: EvaporationTimeError,
    },
}
impl From<EvaporationImplicitError> for AdaptiveEvaporationError {
    fn from(error: EvaporationImplicitError) -> Self {
        Self::Implicit(error)
    }
}
impl From<ThermalEvaporationError> for AdaptiveEvaporationError {
    fn from(error: ThermalEvaporationError) -> Self {
        Self::Implicit(error.into())
    }
}
impl From<DepositionError> for AdaptiveEvaporationError {
    fn from(error: DepositionError) -> Self {
        Self::Implicit(error.into())
    }
}
impl core::fmt::Display for AdaptiveEvaporationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive evaporation: {self:?}")
    }
}
impl std::error::Error for AdaptiveEvaporationError {}

#[derive(Default)]
struct Work {
    iterations: usize,
    evaluations: usize,
}

#[derive(Default)]
struct SurfaceAccumulator {
    mass: EnergySum,
    latent: EnergySum,
    vapor: EnergySum,
    recession: EnergySum,
    intervals: usize,
}
impl SurfaceAccumulator {
    fn add(&mut self, step: SurfaceStep) -> Result<(), DepositionError> {
        self.mass.add(step.transferred_mass_kg)?;
        self.latent.add(step.latent_energy_j)?;
        self.vapor.add(step.vapor_energy_j)?;
        self.recession.add(step.recession_m)?;
        self.intervals = add_count(self.intervals, 1)?;
        Ok(())
    }
    fn finish(self) -> Result<EvaporationSurfaceTotals, DepositionError> {
        Ok(EvaporationSurfaceTotals {
            transferred_mass_kg: self.mass.value()?,
            latent_energy_j: self.latent.value()?,
            vapor_energy_j: self.vapor.value()?,
            recession_m: self.recession.value()?,
            intervals: self.intervals,
        })
    }
}

impl<'m> DepositionSimulation<'m> {
    /// Advance one heat network and all selected liquid interfaces with local
    /// h/T/mass time-error control. An empty binding set retains heat-only BE.
    ///
    /// Start-time events apply before comparison. All three solves end before
    /// the next event changes topology or exposure; only the two fine intervals
    /// contribute physical heat, vapor, storage and payload. The final payload
    /// precedes endpoint events. No extrapolated mass or enthalpy is installed.
    /// A selected future birth activates only after its original event. Before
    /// then, existing cells retain ordinary heat transport and only active
    /// selected interfaces exchange vapor. Unscheduled cell identities refuse.
    ///
    /// Every comparison, including accuracy rejection, spends the same Newton
    /// and residual-evaluation budgets. A fine pair needs two native substep
    /// slots. Source, nonlinear, material and cancellation failures propagate;
    /// they do not trigger an unaccounted retry or select a different model.
    /// A final/event-clipped span below the floor still has to pass time accuracy.
    ///
    /// Callbacks must preserve physical history during speculative evaluations.
    /// Fixed reference geometry and prescribed films retain their original
    /// limits. Backward Euler remains first order; this estimator is local and
    /// requires independent validation near constitutive corners.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_adaptive_evaporating_network_to<L: EvaporationLaw + ?Sized, T>(
        &self,
        bindings: &[EvaporatingCellBinding<'_, L>],
        target_s: f64,
        controls: AdaptiveEvaporationControls,
        mut resolve: impl FnMut(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
        )
            -> Result<BoundaryExchangeStepCoefficients<T>, ThermalEvaporationError>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, Option<T>, AdaptiveEvaporationReport), AdaptiveEvaporationError> {
        checkpoint()?;
        controls.validate()?;
        if !target_s.is_finite() || target_s < self.time_s || bindings.len() > self.budget.max_cells
        {
            return Err(AdaptiveEvaporationError::Invalid(
                "invalid adaptive target or interface count",
            ));
        }
        let mut selected = BTreeMap::new();
        let mut active_bindings = Vec::with_capacity(bindings.len());
        let mut active_ids = BTreeSet::new();
        let mut unscheduled = BTreeSet::new();
        for (i, binding) in bindings.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if selected.insert(binding.id, *binding).is_some() {
                return Err(AdaptiveEvaporationError::Invalid(
                    "duplicate liquid interface inventory",
                ));
            }
            binding
                .film
                .validate()
                .map_err(ThermalEvaporationError::from)?;
            LiquidBranch::new(binding.id, binding.pure_liquid_enthalpy_range_j_kg)?.check(self)?;
            if self.cells.contains_key(&binding.id) {
                active_bindings.push(*binding);
                active_ids.insert(binding.id);
            } else {
                unscheduled.insert(binding.id);
            }
        }
        // Resolve future identities once, before any source callback. Actual
        // birth admission remains with the original event owner at its time;
        // a selected birth beyond the requested target stays pending.
        if !unscheduled.is_empty() {
            for (index, event) in self.schedule[self.cursor..].iter().enumerate() {
                if index % 256 == 0 {
                    checkpoint()?;
                }
                if let Event::Birth(cell) = &event.event {
                    unscheduled.remove(&cell.id);
                }
            }
            if !unscheduled.is_empty() {
                return Err(AdaptiveEvaporationError::Invalid(
                    "selected liquid interface is neither active nor a scheduled birth",
                ));
            }
        }
        finite(self.mass_kg())?;
        finite(self.energy_j())?;
        let mut next = self.clone();
        let mut native = EnergyAudit {
            incoming_mass_kg: 0.0,
            incoming_energy_j: 0.0,
            external_heat_j: 0.0,
            stored_energy_change_j: 0.0,
            residual_j: 0.0,
            substeps: 0,
            events: 0,
        };
        let mut account = TimelineAccount::default();
        let mut work = Work::default();
        let mut absolute_energy = EnergySum::default();
        let mut absolute_mass_rate = EnergySum::default();
        let mut surface_totals = BTreeMap::<CellId, SurfaceAccumulator>::new();
        let mut comparisons = 0usize;
        let mut accepted_comparisons = 0usize;
        let mut accuracy_rejections = 0usize;
        let mut reductions = 0usize;
        let mut maximum_accepted_residual_norm_j = 0.0_f64;
        let mut maximum_accepted_normalized_error = 0.0_f64;
        let mut last_error = None;
        let mut last_surfaces = BTreeMap::new();
        let mut last_payload = None;
        let mut span = controls.initial_dt_s.min(controls.maximum_dt_s);
        let order = ConductionTimeOrder::First;
        loop {
            checkpoint()?;
            while next.cursor < next.schedule.len()
                && next.schedule[next.cursor].time_s == next.time_s
            {
                checkpoint()?;
                let event = next.schedule[next.cursor].event.clone();
                let birth = match &event {
                    Event::Birth(cell) => Some((cell.id, cell.mass_kg, cell.specific_h_j_kg)),
                    _ => None,
                };
                next.apply(event, &mut native)?;
                next.cursor += 1;
                native.events = add_count(native.events, 1)?;
                if let Some((id, mass, h)) = birth {
                    account.birth(id, mass, h)?;
                    if let Some(binding) = selected.get(&id) {
                        active_bindings.push(*binding);
                        active_ids.insert(id);
                    }
                }
                for binding in &active_bindings {
                    LiquidBranch::new(binding.id, binding.pure_liquid_enthalpy_range_j_kg)?
                        .check(&next)?;
                }
                checkpoint()?;
            }
            if next.time_s == target_s {
                break;
            }
            let stop = next
                .schedule
                .get(next.cursor)
                .map_or(target_s, |event| event.time_s.min(target_s));
            if next.cells.is_empty() {
                next.time_s = stop;
                continue;
            }
            if next.budget.max_substeps.saturating_sub(native.substeps) < 2 {
                return Err(
                    DepositionError::Budget("adaptive evaporation fine half intervals").into(),
                );
            }
            if comparisons >= controls.maximum_comparisons {
                return Err(AdaptiveEvaporationError::ComparisonBudget);
            }
            let start = next.time_s;
            let end = if span >= stop - start {
                stop
            } else {
                finite(start + span)?
            };
            let middle = f64::midpoint(start, end);
            if !(start < middle && middle < end && end <= stop) {
                return Err(AdaptiveEvaporationError::Invalid(
                    "no representable adaptive midpoint",
                ));
            }
            let dt = finite(end - start)?;
            comparisons = add_count(comparisons, 1)?;
            let (coarse, _, _) = solve_interval(
                &next,
                &active_bindings,
                end,
                controls.solve,
                &mut work,
                &mut resolve,
                &mut checkpoint,
            )?;
            let (half, _, first) = solve_interval(
                &next,
                &active_bindings,
                middle,
                controls.solve,
                &mut work,
                &mut resolve,
                &mut checkpoint,
            )?;
            let (fine, payload, second) = solve_interval(
                &half,
                &active_bindings,
                end,
                controls.solve,
                &mut work,
                &mut resolve,
                &mut checkpoint,
            )?;
            let correction = order
                .fine_error_correction(middle - start, end - middle)
                .ok_or(DepositionError::Arithmetic)?;
            let error = time_error(
                &next,
                &coarse,
                &fine,
                &active_ids,
                controls,
                correction,
                &mut checkpoint,
            )?;
            if error.normalized <= 1.0 {
                for part in [&first, &second] {
                    account.heat(part.physical.external_heat_j)?;
                    absolute_energy.add(part.physical.energy_residual_j.abs())?;
                    for (&id, surface) in &part.physical.surfaces {
                        account.phase(self, id, *surface)?;
                        surface_totals.entry(id).or_default().add(*surface)?;
                        absolute_mass_rate.add(surface.mass_rate_residual_kg.abs())?;
                    }
                    maximum_accepted_residual_norm_j =
                        maximum_accepted_residual_norm_j.max(part.maximum_accepted_residual_norm_j);
                }
                let energy_defect = absolute_energy.value()?;
                if energy_defect > controls.solve.total_energy_tolerance_j {
                    return Err(DepositionError::EnergyBalance(energy_defect).into());
                }
                let mass_defect = absolute_mass_rate.value()?;
                if mass_defect > controls.solve.mass_tolerance_kg {
                    return Err(ThermalEvaporationError::MassBalance {
                        residual_kg: mass_defect,
                        tolerance_kg: controls.solve.mass_tolerance_kg,
                    }
                    .into());
                }
                native.substeps = add_count(native.substeps, 2)?;
                accepted_comparisons = add_count(accepted_comparisons, 1)?;
                maximum_accepted_normalized_error =
                    maximum_accepted_normalized_error.max(error.normalized);
                last_error = Some(error);
                last_surfaces = second.physical.surfaces;
                last_payload = Some(payload);
                let growth = order
                    .accepted_growth_factor(error.normalized)
                    .ok_or(DepositionError::Arithmetic)?;
                let suggested = if growth >= 1.0 {
                    dt.min(controls.maximum_dt_s / growth) * growth
                } else {
                    dt * growth
                };
                span = suggested.clamp(controls.minimum_dt_s, controls.maximum_dt_s);
                next = fine;
                reductions = 0;
                continue;
            }
            accuracy_rejections = add_count(accuracy_rejections, 1)?;
            let cause = AdaptiveEvaporationError::TimeAccuracy {
                time_s: start,
                span_s: dt,
                error,
            };
            if reductions >= controls.maximum_reductions || dt <= controls.minimum_dt_s {
                return Err(cause);
            }
            let factor = order
                .rejected_reduction_factor(error.normalized)
                .ok_or(DepositionError::Arithmetic)?;
            let wanted = (dt * factor).max(controls.minimum_dt_s);
            let mut smaller = finite(start + wanted)?;
            if smaller - start < controls.minimum_dt_s {
                smaller = smaller.next_up();
            }
            if !(start < smaller && smaller < end) {
                return Err(cause);
            }
            span = smaller - start;
            reductions = add_count(reductions, 1)?;
        }
        let physical = account.finish(
            self,
            &next,
            native.substeps,
            native.events,
            controls.solve.total_energy_tolerance_j,
            Some(controls.solve.mass_tolerance_kg),
            &mut checkpoint,
        )?;
        let surface_totals = surface_totals
            .into_iter()
            .enumerate()
            .map(|(index, (id, sum))| {
                if index % 256 == 0 {
                    checkpoint()?;
                }
                Ok((id, sum.finish()?))
            })
            .collect::<Result<BTreeMap<_, _>, DepositionError>>()?;
        let report = AdaptiveEvaporationReport {
            physical,
            surface_totals,
            last_surfaces,
            newton_iterations: work.iterations,
            residual_evaluations: work.evaluations,
            maximum_accepted_residual_norm_j,
            absolute_energy_residual_j: absolute_energy.value()?,
            absolute_mass_rate_residual_kg: absolute_mass_rate.value()?,
            comparisons,
            accepted_comparisons,
            accuracy_rejections,
            maximum_accepted_normalized_error,
            last_error,
            suggested_dt_s: span,
        };
        checkpoint()?;
        Ok((next, last_payload, report))
    }
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn solve_interval<'m, L: EvaporationLaw + ?Sized, T>(
    initial: &DepositionSimulation<'m>,
    bindings: &[EvaporatingCellBinding<'_, L>],
    target: f64,
    mut controls: EvaporationNetworkImplicitControls,
    work: &mut Work,
    resolve: &mut impl FnMut(
        &DepositionSimulation<'m>,
        &BTreeMap<CellId, LiquidSurfaceState>,
    )
        -> Result<BoundaryExchangeStepCoefficients<T>, ThermalEvaporationError>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<
    (
        DepositionSimulation<'m>,
        T,
        EvaporationNetworkImplicitReport,
    ),
    AdaptiveEvaporationError,
> {
    checkpoint()?;
    controls.maximum_iterations = controls.maximum_iterations.saturating_sub(work.iterations);
    controls.maximum_evaluations = controls
        .maximum_evaluations
        .saturating_sub(work.evaluations);
    if controls.maximum_iterations == 0 || controls.maximum_evaluations == 0 {
        return Err(DepositionError::Budget("adaptive evaporation whole-call solver work").into());
    }
    let (next, payload, report) = initial
        .solve_implicit_evaporating_network_step(bindings, target, controls, resolve, checkpoint)?;
    work.iterations = add_count(work.iterations, report.newton_iterations)?;
    work.evaluations = add_count(work.evaluations, report.residual_evaluations)?;
    Ok((next, payload, report))
}

fn time_error(
    old: &DepositionSimulation<'_>,
    coarse: &DepositionSimulation<'_>,
    fine: &DepositionSimulation<'_>,
    ids: &BTreeSet<CellId>,
    controls: AdaptiveEvaporationControls,
    correction: f64,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<EvaporationTimeError, DepositionError> {
    let mut error = EvaporationTimeError::default();
    for (index, (&id, cell)) in old.cells.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let (difference, ratio) = component_error(
            cell.specific_h_j_kg,
            coarse.cells[&id].specific_h_j_kg,
            fine.cells[&id].specific_h_j_kg,
            controls.enthalpy_absolute_tolerance_j_kg,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_enthalpy_difference_j_kg =
            error.maximum_enthalpy_difference_j_kg.max(difference);
        error.normalized = error.normalized.max(ratio);
        let (difference, ratio) = component_error(
            old.temperature_k(id)?,
            coarse.temperature_k(id)?,
            fine.temperature_k(id)?,
            controls.temperature_absolute_tolerance_k,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_temperature_difference_k =
            error.maximum_temperature_difference_k.max(difference);
        error.normalized = error.normalized.max(ratio);
        if ids.contains(&id) {
            let (difference, ratio) = component_error(
                cell.mass_kg,
                coarse.cells[&id].mass_kg,
                fine.cells[&id].mass_kg,
                controls.mass_absolute_tolerance_kg,
                controls.relative_tolerance,
                correction,
            )?;
            error.maximum_mass_difference_kg = error.maximum_mass_difference_kg.max(difference);
            error.normalized = error.normalized.max(ratio);
        }
    }
    checkpoint()?;
    Ok(error)
}

fn component_error(
    old: f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
    correction: f64,
) -> Result<(f64, f64), DepositionError> {
    let difference = finite(fine - coarse)?.abs();
    let change = finite(coarse - old)?.abs().max(finite(fine - old)?.abs());
    let allowance = finite(absolute + finite(relative * change)?)?;
    let ratio = (difference / allowance) * correction;
    if ratio.is_nan() || ratio < 0.0 || !positive(allowance) {
        return Err(DepositionError::Arithmetic);
    }
    Ok((difference, ratio))
}

fn add_count(value: usize, by: usize) -> Result<usize, DepositionError> {
    value
        .checked_add(by)
        .ok_or(DepositionError::Budget("adaptive evaporation work counter"))
}
