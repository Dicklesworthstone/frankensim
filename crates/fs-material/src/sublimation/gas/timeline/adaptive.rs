//! Adaptive integration of the COMPLETE heat/phase split, not phase alone.
//!
//! One full original heat-plus-molecular interval is compared with two half
//! intervals. Only the conservative fine endpoint is accepted. Every thermal
//! cell's temperature, mass and caloric increment, plus every catalog species
//! inventory and gas temperature, enters the error norm. Caloric increments
//! catch latent plateaus without making the norm depend on an energy datum.
//!
//! Trials end at scheduled events. Initial events are processed once; empty
//! pre-birth intervals jump directly to the next event. Coarse/rejected trials
//! contribute no physical history, but ALL their thermal substeps spend the
//! original timeline's single work budget. Coefficients and source laws stay
//! fixed. Inventory refusals alone permit retry; source-domain and caloric
//! failures are not silently reclassified as temporal error.
//!
//! This is first-order step doubling with a LOCAL estimate, not a rigorous
//! global error bound, positivity/entropy proof, or spatial convergence claim.

use super::{
    GasCellAudit, GasTimelineAudit, GasTimelineCandidate, GasTimelineError,
    GasTimelineSurface, HeatTotals, Sum, ThermalWork, advance_deposition, finite,
};
use crate::sublimation::gas::{
    GasSurfaceError, accommodation::ThermalAccommodation, chamber::GasChamberSolve,
};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation, Event, StepConfig};
use std::collections::BTreeMap;

/// Physical-state error scales and a bound on COMPLETE coupled trials.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasTimelineAdaptiveSolve {
    /// Absolute error allowance for each cell and gas-species mass [kg].
    pub absolute_mass_tolerance_kg: f64,
    /// Absolute error allowance for each thermal cell and gas temperature [K].
    pub absolute_temperature_tolerance_k: f64,
    /// Absolute caloric error scale [J/kg], including isothermal latent plateaus.
    pub absolute_specific_enthalpy_tolerance_j_kg: f64,
    /// Relative temperature/mass scale; caloric scaling uses only h increments,
    /// never absolute formation or reference enthalpy.
    pub relative_tolerance: f64,
    /// Minimum retried interval. Shorter event/final intervals may be tried once.
    pub minimum_dt_s: f64,
    /// Includes rejected trials; each trial evaluates at most three intervals.
    /// The physical maximum interval remains StepConfig::max_dt_s.
    pub maximum_trials: usize,
}

/// No partial accepted history escapes on any refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum GasTimelineAdaptiveError {
    /// Invalid controls, mismatched state identity or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original heat, molecular, source or conservation refusal.
    Timeline(GasTimelineError),
    /// All completed and rejected attempts count against this ceiling.
    TrialBudget { maximum_trials: usize },
    /// Requested retry fell below the physical minimum or cannot advance time.
    MinimumStep { requested_dt_s: f64, minimum_dt_s: f64 },
}
impl From<GasTimelineError> for GasTimelineAdaptiveError {
    fn from(error: GasTimelineError) -> Self { Self::Timeline(error) }
}
impl From<DepositionError> for GasTimelineAdaptiveError {
    fn from(error: DepositionError) -> Self { Self::Timeline(error.into()) }
}
impl core::fmt::Display for GasTimelineAdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive thermal/gas timeline: {self:?}")
    }
}
impl std::error::Error for GasTimelineAdaptiveError {}

/// Numerical work and independent complete-window conservation checks.
#[derive(Clone, Debug, PartialEq)]
pub struct GasTimelineAdaptiveAudit {
    /// Includes rejected full/fine attempts, excludes exact empty-time jumps.
    pub trials: usize,
    /// Error or inventory rejections, never source-domain refusals.
    pub rejected_trials: usize,
    /// Each accepted interval consists of two original heat/phase half steps.
    pub accepted_intervals: usize,
    /// Original driver calls, including admission, coarse, rejected and empty.
    pub interval_evaluations: usize,
    /// Thermal substeps actually computed, INCLUDING coarse/rejected work.
    pub thermal_substeps_evaluated: usize,
    /// Empty pre-birth advances clipped directly to the next event/target.
    pub dormant_advances: usize,
    /// Maximum local error ratio over accepted intervals; at most one.
    pub maximum_accepted_error_ratio: f64,
    /// Largest actual per-species mass defect over the WHOLE window [kg].
    pub maximum_species_mass_residual_kg: f64,
    /// Actual gas storage change minus the accepted integrated stream energy [J].
    pub gas_stream_energy_residual_j: f64,
    /// Actual solid storage change minus births and heat plus stream energy [J].
    pub solid_stream_energy_residual_j: f64,
}

/// Original physical endpoint with a whole-window audit, plus adaptive work.
#[derive(Clone)]
pub struct GasTimelineAdaptiveCandidate<'thermal, 'catalog> {
    endpoint: GasTimelineCandidate<'thermal, 'catalog>,
    audit: GasTimelineAdaptiveAudit,
}
impl<'thermal, 'catalog> GasTimelineAdaptiveCandidate<'thermal, 'catalog> {
    /// Its physical audit covers ALL accepted intervals, not merely the last.
    #[must_use]
    pub const fn endpoint(&self) -> &GasTimelineCandidate<'thermal, 'catalog> { &self.endpoint }
    /// Local-error and executed-work account; distinct from physical history.
    #[must_use]
    pub const fn audit(&self) -> &GasTimelineAdaptiveAudit { &self.audit }
    /// Publish both physical owners only after the enclosing consumer accepts.
    #[must_use]
    pub fn into_parts(self) -> (GasTimelineCandidate<'thermal, 'catalog>, GasTimelineAdaptiveAudit) {
        (self.endpoint, self.audit)
    }
}

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn increment(value: &mut usize) -> Result<(), DepositionError> {
    *value = value.checked_add(1).ok_or(DepositionError::Arithmetic)?;
    Ok(())
}
fn ratio(a: f64, b: f64, scale: f64) -> Result<f64, GasTimelineAdaptiveError> {
    let result = (a - b).abs() / scale;
    if !positive(scale) || !result.is_finite() {
        Err(GasTimelineAdaptiveError::Invalid("unrepresentable coupled error scale"))
    } else { Ok(result) }
}
fn relative_scale(a: f64, b: f64, initial: f64, absolute: f64, relative: f64) -> f64 {
    absolute + relative * a.abs().max(b.abs()).max(initial.abs())
}

// All cells participate, not just sublimating surfaces. For a new endpoint
// birth both trajectories have the identical incoming state, so it adds no
// pre-birth error. IDs/material bindings may never silently zip out of sync.
fn compare(
    before: &DepositionSimulation<'_>,
    original_gas: &ClosedGasParcel,
    coarse: &GasTimelineCandidate<'_, '_>,
    fine: &GasTimelineCandidate<'_, '_>,
    controls: GasTimelineAdaptiveSolve,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<f64, GasTimelineAdaptiveError> {
    let a = coarse.gas().state().state().temperature;
    let b = fine.gas().state().state().temperature;
    let mut maximum = ratio(a, b, relative_scale(a, b, original_gas.state().state().temperature,
        controls.absolute_temperature_tolerance_k, controls.relative_tolerance))?;
    if coarse.thermal().time_s() != fine.thermal().time_s()
        || coarse.thermal().cells().len() != fine.thermal().cells().len()
    {
        return Err(GasTimelineAdaptiveError::Invalid("coupled trial endpoints disagree on time or cell roster"));
    }
    for (&id, a) in coarse.thermal().cells() {
        checkpoint()?;
        let b = fine.thermal().cells().get(&id)
            .ok_or(GasTimelineAdaptiveError::Invalid("coupled trial cell identity changed"))?;
        if a.material_id != b.material_id || a.geometry != b.geometry {
            return Err(GasTimelineAdaptiveError::Invalid("coupled trial material or geometry changed"));
        }
        let initial = before.cells().get(&id).unwrap_or(a);
        maximum = maximum.max(ratio(a.mass_kg, b.mass_kg,
            relative_scale(a.mass_kg, b.mass_kg, initial.mass_kg,
                controls.absolute_mass_tolerance_kg, controls.relative_tolerance))?);
        let ta = coarse.thermal().temperature_at_specific_enthalpy_k(id, a.specific_h_j_kg)?;
        let tb = fine.thermal().temperature_at_specific_enthalpy_k(id, b.specific_h_j_kg)?;
        let t0 = coarse.thermal().temperature_at_specific_enthalpy_k(id, initial.specific_h_j_kg)?;
        maximum = maximum.max(ratio(ta, tb, relative_scale(ta, tb, t0,
            controls.absolute_temperature_tolerance_k, controls.relative_tolerance))?);
        let caloric_scale = controls.absolute_specific_enthalpy_tolerance_j_kg
            + controls.relative_tolerance * (a.specific_h_j_kg - initial.specific_h_j_kg).abs()
                .max((b.specific_h_j_kg - initial.specific_h_j_kg).abs());
        maximum = maximum.max(ratio(a.specific_h_j_kg, b.specific_h_j_kg, caloric_scale)?);
    }
    for species in coarse.catalog().thermodynamics().components() {
        checkpoint()?;
        let mass = |gas: &ClosedGasParcel| gas.species_mass_kg(species.species())
            .map_err(GasSurfaceError::from).map_err(GasTimelineError::from);
        let a = mass(coarse.gas())?;
        let b = mass(fine.gas())?;
        maximum = maximum.max(ratio(a, b, relative_scale(a, b, mass(original_gas)?,
            controls.absolute_mass_tolerance_kg, controls.relative_tolerance))?);
    }
    Ok(maximum)
}

#[derive(Default)]
struct CellTotals {
    outgoing: Sum, incoming: Sum, transferred: Sum, energy: Sum, reflected: Sum,
    recession: Sum, steps: usize,
}
impl CellTotals {
    fn add(&mut self, row: &GasCellAudit) -> Result<(), DepositionError> {
        self.outgoing.add(row.outgoing_mass_kg)?;
        self.incoming.add(row.incoming_mass_kg)?;
        self.transferred.add(row.transferred_mass_kg)?;
        self.energy.add(row.energy_into_gas_j)?;
        self.reflected.add(row.reflected_heat_into_gas_j)?;
        self.recession.add(row.recession_m)?;
        self.steps = self.steps.checked_add(row.phase_steps).ok_or(DepositionError::Arithmetic)?;
        Ok(())
    }
    fn report(self) -> Result<GasCellAudit, DepositionError> {
        Ok(GasCellAudit {
            outgoing_mass_kg: self.outgoing.value()?, incoming_mass_kg: self.incoming.value()?,
            transferred_mass_kg: self.transferred.value()?, energy_into_gas_j: self.energy.value()?,
            reflected_heat_into_gas_j: self.reflected.value()?, recession_m: self.recession.value()?,
            phase_steps: self.steps,
        })
    }
}
#[derive(Default)]
struct History {
    heat: HeatTotals, cells: BTreeMap<CellId, CellTotals>, outgoing: Sum, incoming: Sum,
    transferred: Sum, energy: Sum, reflected: Sum, birth_energy: Sum,
    solid_energy: Sum, gas_energy: Sum, solid_residual: Sum, gas_residual: Sum,
    phase_steps: usize, gas_iterations: usize,
}
impl History {
    fn add(&mut self, candidate: &GasTimelineCandidate<'_, '_>) -> Result<(), DepositionError> {
        let report = candidate.audit();
        self.heat.add(report.thermal)?;
        self.outgoing.add(report.outgoing_mass_kg)?;
        self.incoming.add(report.incoming_mass_kg)?;
        self.transferred.add(report.transferred_mass_kg)?;
        self.energy.add(report.energy_into_gas_j)?;
        self.reflected.add(report.reflected_heat_into_gas_j)?;
        self.birth_energy.add(report.incoming_energy_j)?;
        self.solid_energy.add(report.solid_energy_change_j)?;
        self.gas_energy.add(report.gas_energy_change_j)?;
        self.solid_residual.add(candidate.solid_stream_energy_residual_j)?;
        self.gas_residual.add(candidate.gas_stream_energy_residual_j)?;
        self.phase_steps = self.phase_steps.checked_add(report.phase_steps).ok_or(DepositionError::Arithmetic)?;
        self.gas_iterations = self.gas_iterations.checked_add(report.gas_iterations).ok_or(DepositionError::Arithmetic)?;
        for (&id, row) in &report.cells { self.cells.entry(id).or_default().add(row)?; }
        Ok(())
    }

    // Check actual whole-window inventories independently. Accepted original
    // caloric defects telescope; recomputing them from rounded absolute endpoint
    // energies would reintroduce each species' arbitrary formation datum.
    #[allow(clippy::too_many_arguments)]
    fn finish(
        self,
        original: &DepositionSimulation<'_>,
        original_gas: &ClosedGasParcel,
        surfaces: &[GasTimelineSurface<'_>],
        endpoint: &mut GasTimelineCandidate<'_, '_>,
        thermal_config: StepConfig,
        phase_config: GasChamberSolve,
        work: &mut GasTimelineAdaptiveAudit,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), GasTimelineAdaptiveError> {
        let mut births = BTreeMap::new();
        for event in original.pending_events() {
            checkpoint()?;
            if event.time_s > endpoint.thermal().time_s() { break; }
            if let Event::Birth(cell) = &event.event { births.entry(cell.id).or_insert(cell); }
        }
        let bindings: BTreeMap<_, _> = surfaces.iter().map(|surface| (surface.cell_id, surface.binding)).collect();
        let mut solid_mass = Sum::default();
        let mut species_mass: BTreeMap<_, Sum> = BTreeMap::new();
        for (&id, cell) in endpoint.thermal().cells() {
            checkpoint()?;
            let initial = original.cells().get(&id).or_else(|| births.get(&id).copied())
                .ok_or(GasTimelineAdaptiveError::Invalid("accepted cell has no original inventory or birth"))?;
            let dm = finite(cell.mass_kg - initial.mass_kg)?;
            solid_mass.add(dm)?;
            if let Some(binding) = bindings.get(&id) {
                species_mass.entry(binding.species.clone()).or_default().add(dm)?;
            } else if dm != 0.0 {
                return Err(GasTimelineAdaptiveError::Invalid("unselected thermal cell changed mass"));
            }
        }
        let gas_mass_change = finite(endpoint.gas().mass_kg() - original_gas.mass_kg())?;
        let gas_energy_change = self.gas_energy.value()?;
        let mass_tolerance = phase_config.surface.gas.mass_tolerance_kg;
        let mass_residual = finite(solid_mass.value()? + gas_mass_change)?;
        for species in endpoint.catalog().thermodynamics().components() {
            checkpoint()?;
            let id = species.species();
            let mass = |gas: &ClosedGasParcel| gas.species_mass_kg(id)
                .map_err(GasSurfaceError::from).map_err(GasTimelineError::from);
            let change = species_mass.get(id).map_or(Ok(0.0), Sum::value)?;
            let residual = finite(change + (mass(endpoint.gas())? - mass(original_gas)?))?;
            work.maximum_species_mass_residual_kg = work.maximum_species_mass_residual_kg.max(residual.abs());
        }
        let largest = mass_residual.abs().max(work.maximum_species_mass_residual_kg);
        if largest > mass_tolerance {
            return Err(GasTimelineError::MassBalance { residual_kg: largest, tolerance_kg: mass_tolerance }.into());
        }
        let heat = self.heat.audit()?;
        let stream_energy = self.energy.value()?;
        work.gas_stream_energy_residual_j = self.gas_residual.value()?;
        work.solid_stream_energy_residual_j = self.solid_residual.value()?;
        let energy_residual = finite(work.solid_stream_energy_residual_j + work.gas_stream_energy_residual_j)?;
        let energy_tolerance = thermal_config.energy_tolerance_j.max(phase_config.surface.gas.energy_tolerance_j);
        for residual in [energy_residual, work.gas_stream_energy_residual_j, work.solid_stream_energy_residual_j] {
            if residual.abs() > energy_tolerance {
                return Err(GasTimelineError::Exchange(GasSurfaceError::EnergyBalance {
                    residual_j: residual, tolerance_j: energy_tolerance,
                }).into());
            }
        }
        let mut cells = BTreeMap::new();
        for (id, total) in self.cells { checkpoint()?; cells.insert(id, total.report()?); }
        endpoint.audit = GasTimelineAudit {
            thermal: heat, cells, outgoing_mass_kg: self.outgoing.value()?,
            incoming_mass_kg: self.incoming.value()?, transferred_mass_kg: self.transferred.value()?,
            energy_into_gas_j: stream_energy, reflected_heat_into_gas_j: self.reflected.value()?,
            solid_energy_change_j: self.solid_energy.value()?,
            gas_energy_change_j: gas_energy_change, incoming_energy_j: self.birth_energy.value()?,
            mass_residual_kg: mass_residual, energy_residual_j: energy_residual,
            phase_steps: self.phase_steps, gas_iterations: self.gas_iterations,
            last_step: endpoint.audit.last_step.take(),
        };
        endpoint.gas_stream_energy_residual_j = work.gas_stream_energy_residual_j;
        endpoint.solid_stream_energy_residual_j = work.solid_stream_energy_residual_j;
        Ok(())
    }
}

fn smaller(requested: f64, old: f64, minimum: f64) -> Result<f64, GasTimelineAdaptiveError> {
    let dt = requested.min(old * 0.5);
    if !positive(dt) || dt < minimum {
        Err(GasTimelineAdaptiveError::MinimumStep { requested_dt_s: dt, minimum_dt_s: minimum })
    } else { Ok(dt) }
}

/// Adapt the entire heat/phase splitting error on the original event timeline.
/// Existing source laws, thermal CFL subcycling and gas reconstruction remain
/// the only physical owners. The fine history is accepted without extrapolation.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_adaptive<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    controls: GasTimelineAdaptiveSolve,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasTimelineAdaptiveCandidate<'thermal, 'catalog>, GasTimelineAdaptiveError> {
    advance(thermal, surfaces, None, gas, catalog, target_s, thermal_config, phase_config, controls, checkpoint)
}

/// The same coupled controller with explicit reflected heat on every selected
/// cell. Birth/contact timing, source support, budgets and rollback are shared.
/// Do not add the separately reported reflected subaccount to total heat again.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_adaptive_with_accommodation<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineSurface<'_>],
    accommodation: &BTreeMap<CellId, ThermalAccommodation>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    controls: GasTimelineAdaptiveSolve,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasTimelineAdaptiveCandidate<'thermal, 'catalog>, GasTimelineAdaptiveError> {
    advance(thermal, surfaces, Some(accommodation), gas, catalog, target_s, thermal_config, phase_config, controls, checkpoint)
}

#[allow(clippy::too_many_arguments)]
fn advance<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineSurface<'_>],
    accommodation: Option<&BTreeMap<CellId, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    controls: GasTimelineAdaptiveSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasTimelineAdaptiveCandidate<'thermal, 'catalog>, GasTimelineAdaptiveError> {
    checkpoint()?;
    if !target_s.is_finite() || target_s < thermal.time_s()
        || !positive(controls.absolute_mass_tolerance_kg)
        || !positive(controls.absolute_temperature_tolerance_k)
        || !positive(controls.absolute_specific_enthalpy_tolerance_j_kg)
        || !controls.relative_tolerance.is_finite() || controls.relative_tolerance < 0.0
        || !positive(controls.minimum_dt_s) || !positive(thermal_config.max_dt_s)
        || controls.minimum_dt_s > thermal_config.max_dt_s || controls.maximum_trials == 0
    {
        return Err(GasTimelineAdaptiveError::Invalid("finite target, positive error scales and bounded coupled-step controls required"));
    }
    let meter = ThermalWork::new(thermal.work_budget().max_substeps);
    let mut current = advance_deposition(thermal, surfaces, accommodation, gas, catalog,
        thermal.time_s(), thermal_config, phase_config, Some(&meter), &mut checkpoint)?;
    let mut history = History::default();
    history.add(&current)?;
    let mut work = GasTimelineAdaptiveAudit {
        trials: 0, rejected_trials: 0, accepted_intervals: 0, interval_evaluations: 1,
        thermal_substeps_evaluated: 0, dormant_advances: 0, maximum_accepted_error_ratio: 0.0,
        maximum_species_mass_residual_kg: 0.0, gas_stream_energy_residual_j: 0.0,
        solid_stream_energy_residual_j: 0.0,
    };
    let mut dt = thermal_config.max_dt_s;
    while current.thermal().time_s() < target_s {
        checkpoint()?;
        let start = current.thermal().time_s();
        let event_stop = current.thermal().pending_events().first().map_or(target_s,
            |event| event.time_s.min(target_s));
        if current.thermal().cells().is_empty() {
            increment(&mut work.interval_evaluations)?;
            current = advance_deposition(current.thermal(), surfaces, accommodation, current.gas(),
                catalog, event_stop, thermal_config, phase_config, Some(&meter), &mut checkpoint)?;
            history.add(&current)?;
            increment(&mut work.dormant_advances)?;
            continue;
        }
        if work.trials >= controls.maximum_trials {
            return Err(GasTimelineAdaptiveError::TrialBudget { maximum_trials: controls.maximum_trials });
        }
        let stop = (start + dt).min(event_stop);
        let span = stop - start;
        let middle = start + 0.5 * span;
        if !positive(span) || !(start < middle && middle < stop) {
            return Err(GasTimelineAdaptiveError::MinimumStep {
                requested_dt_s: span * 0.5, minimum_dt_s: controls.minimum_dt_s,
            });
        }
        increment(&mut work.trials)?;
        let mut evaluate = |state: &GasTimelineCandidate<'thermal, 'catalog>, end: f64| {
            increment(&mut work.interval_evaluations).map_err(GasTimelineError::from)?;
            advance_deposition(state.thermal(), surfaces, accommodation, state.gas(), catalog,
                end, StepConfig { max_dt_s: end - state.thermal().time_s(), ..thermal_config },
                phase_config, Some(&meter), &mut checkpoint)
        };
        let trial = (|| {
            let coarse = evaluate(&current, stop)?;
            let first = evaluate(&current, middle)?;
            let fine = evaluate(&first, stop)?;
            Ok::<_, GasTimelineError>((coarse, first, fine))
        })();
        let (coarse, first, fine) = match trial {
            Ok(value) => value,
            Err(GasTimelineError::Exchange(GasSurfaceError::InventoryLimit { maximum_duration_s })) => {
                increment(&mut work.rejected_trials)?;
                dt = smaller(0.8 * maximum_duration_s, span, controls.minimum_dt_s)?;
                continue;
            }
            Err(error) => return Err(error.into()),
        };
        let error = compare(current.thermal(), current.gas(), &coarse, &fine, controls, &mut checkpoint)?;
        if error > 1.0 {
            increment(&mut work.rejected_trials)?;
            dt = smaller(span * (0.9 / fs_math::det::sqrt(error)).max(0.1), span, controls.minimum_dt_s)?;
            continue;
        }
        history.add(&first)?;
        history.add(&fine)?;
        current = fine;
        increment(&mut work.accepted_intervals)?;
        work.maximum_accepted_error_ratio = work.maximum_accepted_error_ratio.max(error);
        let growth = if error == 0.0 { 2.0 } else { (0.9 / fs_math::det::sqrt(error)).clamp(0.5, 2.0) };
        dt = (span * growth).min(thermal_config.max_dt_s).max(controls.minimum_dt_s);
    }
    history.finish(thermal, gas, surfaces, &mut current, thermal_config, phase_config, &mut work, &mut checkpoint)?;
    work.thermal_substeps_evaluated = meter.used();
    checkpoint()?;
    Ok(GasTimelineAdaptiveCandidate { endpoint: current, audit: work })
}

#[cfg(test)]
#[path = "adaptive_tests.rs"]
mod tests;
