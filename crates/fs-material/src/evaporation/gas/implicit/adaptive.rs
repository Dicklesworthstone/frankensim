//! Adaptive physical time for the original simultaneous thermal/gas endpoint.
//!
//! Only two conservative fine intervals are accepted. The coarse comparison
//! and every rejected comparison spend the same numerical budgets but supply
//! no physical account or final source payload. The event-free and native
//! deposition timeline APIs share this comparison, work and conservation owner.

use super::super::chamber::Total;
use super::timeline::{Birth, Roster};
use super::*;
use fs_thermal::enthalpy::EnergyAudit;
use fs_thermal::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;
use fs_thermal::enthalpy::evaporation::implicit::AdaptiveEvaporationControls;
use fs_thermal::evaporation::SurfaceStep;
use fs_thermochem::SpeciesId;
use std::collections::BTreeSet;

/// Local time accuracy, whole-window conservation and shared numerical work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveImplicitGasControls {
    pub solve: ImplicitGasControls,
    pub initial_dt_s: f64,
    /// Bounds on the full comparison span, before its two fine intervals.
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    pub maximum_comparisons: usize,
    pub maximum_reductions: usize,
    /// Total original gas reconstruction iterations, including discarded work.
    /// Each inversion also retains the original per-solve iteration ceiling.
    pub maximum_gas_iterations: usize,
    pub enthalpy_absolute_tolerance_j_kg: f64,
    /// Applies independently to every thermal cell and the actual finite gas.
    pub temperature_absolute_tolerance_k: f64,
    /// Applies independently to each selected liquid and each gas species.
    pub mass_absolute_tolerance_kg: f64,
    /// Scales changes from the interval start, never absolute h/T/inventory.
    pub relative_tolerance: f64,
}
impl AdaptiveImplicitGasControls {
    pub fn new(solve: ImplicitGasControls) -> Self {
        let defaults = AdaptiveEvaporationControls::new(solve.thermal);
        Self {
            solve,
            initial_dt_s: defaults.initial_dt_s,
            minimum_dt_s: defaults.minimum_dt_s,
            maximum_dt_s: defaults.maximum_dt_s,
            maximum_comparisons: defaults.maximum_comparisons,
            maximum_reductions: defaults.maximum_reductions,
            maximum_gas_iterations: 100_000,
            enthalpy_absolute_tolerance_j_kg: defaults.enthalpy_absolute_tolerance_j_kg,
            temperature_absolute_tolerance_k: defaults.temperature_absolute_tolerance_k,
            mass_absolute_tolerance_kg: defaults.mass_absolute_tolerance_kg,
            relative_tolerance: defaults.relative_tolerance,
        }
    }

    fn validate(self, duration_s: f64) -> Result<(), AdaptiveImplicitGasError> {
        self.solve.validate(duration_s)?;
        if !positive(self.initial_dt_s)
            || !positive(self.minimum_dt_s)
            || !positive(self.maximum_dt_s)
            || self.minimum_dt_s > self.maximum_dt_s
            || self.initial_dt_s < self.minimum_dt_s
            || self.maximum_comparisons == 0
            || self.maximum_gas_iterations == 0
            || !positive(self.enthalpy_absolute_tolerance_j_kg)
            || !positive(self.temperature_absolute_tolerance_k)
            || !positive(self.mass_absolute_tolerance_kg)
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
        {
            return Err(AdaptiveImplicitGasError::Invalid(
                "invalid adaptive finite-gas controls",
            ));
        }
        Ok(())
    }
}

/// Measured coarse/fine discrepancies and their corrected local-error estimate.
/// This is not a certified bound on accumulated trajectory error.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ImplicitGasTimeError {
    pub maximum_enthalpy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_liquid_mass_difference_kg: f64,
    pub gas_temperature_difference_k: f64,
    pub maximum_species_mass_difference_kg: f64,
    pub normalized: f64,
}

/// Accepted signed transfers for one actual liquid cell in its original datum.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AdaptiveImplicitGasSurfaceTotals {
    pub transferred_mass_kg: f64,
    pub evaporated_mass_kg: f64,
    /// Positive magnitude of condensation, separately from signed net transfer.
    pub condensed_mass_kg: f64,
    /// Already included in vapor_energy_j; never add it a second time.
    pub latent_energy_j: f64,
    pub vapor_energy_j: f64,
    /// Also included in vapor_energy_j, retained before large datum additions.
    pub vapor_enthalpy_correction_j: f64,
    pub sensible_heat_into_liquid_j: f64,
    /// Signed physical distance for this surface only, without geometry remap.
    pub recession_m: f64,
    pub intervals: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveImplicitGasReport {
    /// Authored native events applied by the timeline; zero for an event-free window.
    pub events: usize,
    /// New native material, separate from internal liquid/gas phase exchange.
    pub incoming_mass_kg: f64,
    /// Birth energy in each cell's original thermal storage datum.
    pub incoming_energy_j: f64,
    /// Selected birth mass times its declared gas/liquid enthalpy offset.
    /// Add to incoming_energy_j only when expressing birth energy in the gas
    /// datum. The conservation gates retain source-relative increments instead.
    pub liquid_datum_incoming_energy_j: f64,
    pub external_heat_j: f64,
    /// Signed sums of accepted source-relative owner defects, not differences
    /// of absolute formation-energy totals. Coarse/rejected work is excluded.
    pub thermal_energy_residual_j: f64,
    pub gas_energy_residual_j: f64,
    pub energy_residual_j: f64,
    pub accumulated_mass_residual_kg: f64,
    /// Actual final total-gas mass change minus accepted phase mass transfer.
    pub gas_mass_residual_kg: f64,
    /// Signed sum of the independent volatile inventory defects below.
    /// Each species is gated separately, so opposing defects cannot cancel.
    pub vapor_mass_residual_kg: f64,
    /// Accepted thermal phase transfer for each declared volatile species.
    /// Future declarations with no active phase interval retain zero transfer.
    pub species_transferred_mass_kg: BTreeMap<SpeciesId, f64>,
    /// Actual final-minus-initial species inventory minus its accepted transfer.
    pub species_mass_residual_kg: BTreeMap<SpeciesId, f64>,
    /// Sum of native cell mass changes from their initial or actual birth
    /// inventories, plus the total-gas change. Incoming births cancel per cell.
    pub mass_residual_kg: f64,
    pub maximum_carrier_mass_residual_kg: f64,
    /// Per fine interval: max(|native thermal account|, sum of selected local
    /// |energy defects|) + |gas energy defect|. Unselected cells retain their
    /// original per-step native gates; their separate residuals are not exposed.
    pub absolute_energy_residual_j: f64,
    /// Per fine interval: max(|native mass account|, sum of selected local
    /// |represented-minus-requested transfers|) + original gas maximum defect.
    pub absolute_mass_residual_kg: f64,
    pub surface_totals: BTreeMap<CellId, AdaptiveImplicitGasSurfaceTotals>,
    pub accepted_intervals: usize,
    /// Includes all coarse and accuracy-rejected comparisons.
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    pub gas_iterations: usize,
    pub maximum_accepted_residual_norm_j: f64,
    pub comparisons: usize,
    pub accepted_comparisons: usize,
    pub accuracy_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<ImplicitGasTimeError>,
    pub suggested_dt_s: f64,
    /// Last accepted fine endpoint, including its complete vapor corrections.
    /// Timeline events at the target follow this endpoint and create no phase step.
    pub last_step: Option<ImplicitGasReport>,
}

#[derive(Debug)]
pub enum AdaptiveImplicitGasError {
    Implicit(ImplicitGasError),
    Invalid(&'static str),
    ComparisonBudget,
    TimeAccuracy {
        time_s: f64,
        span_s: f64,
        error: ImplicitGasTimeError,
    },
}
impl From<ImplicitGasError> for AdaptiveImplicitGasError {
    fn from(error: ImplicitGasError) -> Self {
        Self::Implicit(error)
    }
}
impl From<DepositionError> for AdaptiveImplicitGasError {
    fn from(error: DepositionError) -> Self {
        Self::Implicit(error.into())
    }
}
impl From<GasEvaporationError> for AdaptiveImplicitGasError {
    fn from(error: GasEvaporationError) -> Self {
        Self::Implicit(error.into())
    }
}
impl From<GasEnergyError> for AdaptiveImplicitGasError {
    fn from(error: GasEnergyError) -> Self {
        Self::Implicit(error.into())
    }
}
impl core::fmt::Display for AdaptiveImplicitGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive implicit finite-gas exchange: {self:?}")
    }
}
impl std::error::Error for AdaptiveImplicitGasError {}

/// The complete accepted fine trajectory, with no partial state on refusal.
pub struct AdaptiveImplicitGasCandidate<'m, 'catalog, T> {
    thermal: DepositionSimulation<'m>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    payload: Option<T>,
    report: AdaptiveImplicitGasReport,
}
impl<'m, 'catalog, T> AdaptiveImplicitGasCandidate<'m, 'catalog, T> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub fn payload(&self) -> Option<&T> {
        self.payload.as_ref()
    }
    pub const fn report(&self) -> &AdaptiveImplicitGasReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        Option<T>,
        AdaptiveImplicitGasReport,
    ) {
        (
            self.thermal,
            self.gas,
            self.catalog,
            self.payload,
            self.report,
        )
    }
}

/// Advance a positive, event-free window with joint h/T/liquid/species accuracy.
/// Pending events strictly before target refuse before any source callback;
/// events exactly at target remain pending. Every selected surface stays active.
/// With an empty roster, existing heat cells advance and gas remains unchanged.
///
/// Three original coupled solves form each comparison. Only its two fine
/// intervals contribute mass, heat and payload. No extrapolated state or
/// old prescribed-film timeline account is used. Source, domain, nonlinear
/// and cancellation errors propagate; only successful inaccurate comparisons
/// trigger a smaller span. Even a final span below the floor must pass accuracy.
///
/// All numerical work shares one request budget, while two original thermal
/// substep slots are reserved for each accepted fine pair. The original first-
/// order controller supplies correction and timestep factors. Callbacks must
/// preserve external accepted history throughout speculative evaluations.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_deposition_with_gas_adaptive<'m, 'catalog, T>(
    thermal: &DepositionSimulation<'m>,
    surfaces: &[ImplicitGasSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    controls: AdaptiveImplicitGasControls,
    resolve_heat: impl FnMut(
        &DepositionSimulation<'m>,
        &BTreeMap<CellId, LiquidSurfaceState>,
        &ClosedGasParcel,
    ) -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<AdaptiveImplicitGasCandidate<'m, 'catalog, T>, AdaptiveImplicitGasError> {
    advance_window(
        thermal,
        surfaces,
        gas,
        catalog,
        target_s,
        controls,
        false,
        resolve_heat,
        checkpoint,
    )
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(super) fn advance_window<'m, 'catalog, T>(
    thermal: &DepositionSimulation<'m>,
    surfaces: &[ImplicitGasSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    controls: AdaptiveImplicitGasControls,
    scheduled: bool,
    mut resolve_heat: impl FnMut(
        &DepositionSimulation<'m>,
        &BTreeMap<CellId, LiquidSurfaceState>,
        &ClosedGasParcel,
    ) -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<AdaptiveImplicitGasCandidate<'m, 'catalog, T>, AdaptiveImplicitGasError> {
    checkpoint()?;
    if !target_s.is_finite() || target_s <= thermal.time_s() {
        return Err(AdaptiveImplicitGasError::Invalid(
            "a finite positive window is required",
        ));
    }
    controls.validate(target_s - thermal.time_s())?;
    if !scheduled
        && thermal
            .pending_events()
            .first()
            .is_some_and(|event| event.time_s < target_s)
    {
        return Err(DepositionError::InvalidInput(
            "process pending events or split implicit step at the next event",
        )
        .into());
    }
    let mut ids = BTreeSet::new();
    for (index, surface) in surfaces.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        if (!scheduled && !thermal.cells().contains_key(&surface.cell_id))
            || !ids.insert(surface.cell_id)
        {
            return Err(AdaptiveImplicitGasError::Invalid(
                "unique active liquid cells required",
            ));
        }
    }
    let roster = if scheduled {
        Some(Roster::new(
            thermal,
            surfaces,
            gas,
            catalog,
            controls.solve,
            &mut checkpoint,
        )?)
    } else {
        None
    };
    let species_ids: BTreeSet<_> = catalog
        .thermodynamics()
        .components()
        .iter()
        .map(|species| species.species().clone())
        .collect();
    let volatile_ids: BTreeSet<_> = surfaces
        .iter()
        .map(|surface| surface.binding.species.clone())
        .collect();
    let mut next_thermal = thermal.clone();
    let mut next_gas = gas.clone();
    let mut active_surfaces = if let Some(roster) = &roster {
        roster.active(&next_thermal, &mut checkpoint)?
    } else {
        Vec::new()
    };
    let mut work = Work::default();
    let mut account = Account::default();
    let mut births = BTreeMap::new();
    let mut accepted_intervals = 0usize;
    let mut comparisons = 0usize;
    let mut accepted_comparisons = 0usize;
    let mut accuracy_rejections = 0usize;
    let mut reductions = 0usize;
    let mut maximum_accepted_normalized_error = 0.0_f64;
    let mut maximum_accepted_residual_norm_j = 0.0_f64;
    let mut last_error = None;
    let mut last_step = None;
    let mut last_payload = None;
    let mut span = controls.initial_dt_s.min(controls.maximum_dt_s);
    let order = ConductionTimeOrder::First;
    loop {
        checkpoint()?;
        if let Some(roster) = &roster {
            let now = next_thermal.time_s();
            let event = next_thermal
                .pending_events()
                .first()
                .map(|event| event.time_s);
            let stop = if event == Some(now) {
                Some(now)
            } else if next_thermal.cells().is_empty() && now < target_s {
                Some(event.map_or(target_s, |time| time.min(target_s)))
            } else {
                None
            };
            if let Some(stop) = stop {
                // The native owner applies authored events and empty-geometry
                // time only. Zero allowed substeps forbids a second heat solve.
                let (advanced, incoming, audit) = roster.advance_events(
                    &next_thermal,
                    &next_gas,
                    catalog,
                    stop,
                    controls.solve,
                    &mut checkpoint,
                )?;
                account.add_events(&incoming, &audit, &mut checkpoint)?;
                account.gate(controls)?;
                for (id, birth) in incoming {
                    if births.insert(id, birth).is_some() {
                        return Err(DepositionError::InvalidInput(
                            "duplicate native birth receipt",
                        )
                        .into());
                    }
                }
                next_thermal = advanced;
                active_surfaces = roster.active(&next_thermal, &mut checkpoint)?;
                if next_thermal.cells().is_empty() && next_thermal.time_s() < target_s {
                    continue;
                }
            }
        }
        if next_thermal.time_s() == target_s {
            break;
        }
        if thermal
            .work_budget()
            .max_substeps
            .saturating_sub(accepted_intervals)
            < 2
        {
            return Err(DepositionError::Budget("thermal substeps").into());
        }
        if comparisons >= controls.maximum_comparisons {
            return Err(AdaptiveImplicitGasError::ComparisonBudget);
        }
        let start = next_thermal.time_s();
        let stop = if scheduled {
            next_thermal
                .pending_events()
                .first()
                .map_or(target_s, |event| event.time_s.min(target_s))
        } else {
            target_s
        };
        let end = if span >= stop - start {
            stop
        } else {
            finite(start + span)?
        };
        let middle = f64::midpoint(start, end);
        if !(start < middle && middle < end && end <= target_s) {
            return Err(AdaptiveImplicitGasError::Invalid(
                "no representable adaptive gas midpoint",
            ));
        }
        let dt = finite(end - start)?;
        let surfaces = if scheduled {
            active_surfaces.as_slice()
        } else {
            surfaces
        };
        comparisons = add_count(comparisons, 1)?;
        let (coarse_thermal, coarse_gas, _, _, _) = solve_interval(
            &next_thermal,
            surfaces,
            &next_gas,
            catalog,
            end,
            controls,
            &mut work,
            &mut resolve_heat,
            &mut checkpoint,
        )?
        .into_parts();
        let (half_thermal, half_gas, _, _, first) = solve_interval(
            &next_thermal,
            surfaces,
            &next_gas,
            catalog,
            middle,
            controls,
            &mut work,
            &mut resolve_heat,
            &mut checkpoint,
        )?
        .into_parts();
        let (fine_thermal, fine_gas, _, payload, second) = solve_interval(
            &half_thermal,
            surfaces,
            &half_gas,
            catalog,
            end,
            controls,
            &mut work,
            &mut resolve_heat,
            &mut checkpoint,
        )?
        .into_parts();
        drop(half_thermal);
        drop(half_gas);
        let correction = order
            .fine_error_correction(middle - start, end - middle)
            .ok_or(DepositionError::Arithmetic)?;
        let error = time_error(
            &next_thermal,
            &coarse_thermal,
            &fine_thermal,
            &ids,
            &next_gas,
            &coarse_gas,
            &fine_gas,
            &species_ids,
            controls,
            correction,
            &mut checkpoint,
        )?;
        drop(coarse_thermal);
        drop(coarse_gas);
        if error.normalized <= 1.0 {
            for part in [&first, &second] {
                account.add(part, &mut checkpoint)?;
                maximum_accepted_residual_norm_j =
                    maximum_accepted_residual_norm_j.max(part.maximum_accepted_residual_norm_j);
            }
            account.gate(controls)?;
            accepted_intervals = add_count(accepted_intervals, 2)?;
            accepted_comparisons = add_count(accepted_comparisons, 1)?;
            maximum_accepted_normalized_error =
                maximum_accepted_normalized_error.max(error.normalized);
            last_error = Some(error);
            last_step = Some(second);
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
            next_thermal = fine_thermal;
            next_gas = fine_gas;
            reductions = 0;
            continue;
        }
        accuracy_rejections = add_count(accuracy_rejections, 1)?;
        let cause = AdaptiveImplicitGasError::TimeAccuracy {
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
    let mut species_transferred_mass_kg = BTreeMap::new();
    for id in &volatile_ids {
        checkpoint()?;
        let transferred = account
            .species_mass
            .get(id)
            .map(Total::value)
            .transpose()?
            .unwrap_or(0.0);
        species_transferred_mass_kg.insert(id.clone(), transferred);
    }
    let mass = final_mass_balance(
        thermal,
        &next_thermal,
        gas,
        &next_gas,
        &births,
        &species_ids,
        &species_transferred_mass_kg,
        account.phase_mass.value()?,
        mass_tolerance(controls),
        &mut checkpoint,
    )?;
    let mut surface_totals = BTreeMap::new();
    for (index, (id, sums)) in account.surfaces.into_iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        surface_totals.insert(id, sums.finish()?);
    }
    let report = AdaptiveImplicitGasReport {
        events: account.events,
        incoming_mass_kg: account.incoming_mass.value()?,
        incoming_energy_j: account.incoming_energy.value()?,
        liquid_datum_incoming_energy_j: account.incoming_offset.value()?,
        external_heat_j: account.external.value()?,
        thermal_energy_residual_j: account.thermal_energy.value()?,
        gas_energy_residual_j: account.gas_energy.value()?,
        energy_residual_j: account.coupled_energy.value()?,
        accumulated_mass_residual_kg: account.coupled_mass.value()?,
        gas_mass_residual_kg: mass.gas,
        vapor_mass_residual_kg: mass.vapor,
        species_transferred_mass_kg,
        species_mass_residual_kg: mass.species,
        mass_residual_kg: mass.coupled,
        maximum_carrier_mass_residual_kg: mass.carrier,
        absolute_energy_residual_j: account.absolute_energy.value()?,
        absolute_mass_residual_kg: account.absolute_mass.value()?,
        surface_totals,
        accepted_intervals,
        newton_iterations: work.iterations,
        residual_evaluations: work.evaluations,
        gas_iterations: work.gas_iterations,
        maximum_accepted_residual_norm_j,
        comparisons,
        accepted_comparisons,
        accuracy_rejections,
        maximum_accepted_normalized_error,
        last_error,
        suggested_dt_s: span,
        last_step,
    };
    checkpoint()?;
    Ok(AdaptiveImplicitGasCandidate {
        thermal: next_thermal,
        gas: next_gas,
        catalog,
        payload: last_payload,
        report,
    })
}

#[derive(Default)]
struct Work {
    iterations: usize,
    evaluations: usize,
    gas_iterations: usize,
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn solve_interval<'m, 'catalog, T>(
    thermal: &DepositionSimulation<'m>,
    surfaces: &[ImplicitGasSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    controls: AdaptiveImplicitGasControls,
    work: &mut Work,
    resolve: &mut impl FnMut(
        &DepositionSimulation<'m>,
        &BTreeMap<CellId, LiquidSurfaceState>,
        &ClosedGasParcel,
    ) -> Result<BoundaryExchangeStepCoefficients<T>, ImplicitGasError>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<ImplicitGasCandidate<'m, 'catalog, T>, AdaptiveImplicitGasError> {
    checkpoint()?;
    let mut solve = controls.solve;
    solve.thermal.maximum_iterations = solve
        .thermal
        .maximum_iterations
        .saturating_sub(work.iterations);
    solve.thermal.maximum_evaluations = solve
        .thermal
        .maximum_evaluations
        .saturating_sub(work.evaluations);
    if solve.thermal.maximum_iterations == 0 || solve.thermal.maximum_evaluations == 0 {
        return Err(DepositionError::Budget("adaptive finite-gas whole-window solver work").into());
    }
    let gas_iterations = controls
        .maximum_gas_iterations
        .saturating_sub(work.gas_iterations);
    if gas_iterations == 0 && !surfaces.is_empty() {
        return Err(
            DepositionError::Budget("adaptive finite-gas whole-window gas iterations").into(),
        );
    }
    if !surfaces.is_empty() {
        solve.phase.gas.gas.maximum_iterations =
            solve.phase.gas.gas.maximum_iterations.min(gas_iterations);
    }
    let next = solve_deposition_with_gas(
        thermal, surfaces, gas, catalog, target_s, solve, resolve, checkpoint,
    )?;
    work.iterations = add_count(work.iterations, next.report().newton_iterations)?;
    work.evaluations = add_count(work.evaluations, next.report().residual_evaluations)?;
    work.gas_iterations = add_count(work.gas_iterations, next.report().gas_iterations)?;
    Ok(next)
}

#[derive(Default)]
struct SurfaceTotals {
    mass: Total,
    evaporated: Total,
    condensed: Total,
    latent: Total,
    vapor: Total,
    correction: Total,
    sensible: Total,
    recession: Total,
    intervals: usize,
}
impl SurfaceTotals {
    fn add(
        &mut self,
        phase: &SurfaceStep,
        correction: f64,
        sensible: f64,
    ) -> Result<(), AdaptiveImplicitGasError> {
        self.mass.add(phase.transferred_mass_kg)?;
        self.evaporated.add(phase.transferred_mass_kg.max(0.0))?;
        self.condensed.add((-phase.transferred_mass_kg).max(0.0))?;
        self.latent.add(phase.latent_energy_j)?;
        self.vapor.add(phase.vapor_energy_j)?;
        self.correction.add(correction)?;
        self.sensible.add(sensible)?;
        self.recession.add(phase.recession_m)?;
        self.intervals = add_count(self.intervals, 1)?;
        Ok(())
    }
    fn finish(self) -> Result<AdaptiveImplicitGasSurfaceTotals, AdaptiveImplicitGasError> {
        Ok(AdaptiveImplicitGasSurfaceTotals {
            transferred_mass_kg: self.mass.value()?,
            evaporated_mass_kg: self.evaporated.value()?,
            condensed_mass_kg: self.condensed.value()?,
            latent_energy_j: self.latent.value()?,
            vapor_energy_j: self.vapor.value()?,
            vapor_enthalpy_correction_j: self.correction.value()?,
            sensible_heat_into_liquid_j: self.sensible.value()?,
            recession_m: self.recession.value()?,
            intervals: self.intervals,
        })
    }
}

#[derive(Default)]
struct Account {
    events: usize,
    incoming_mass: Total,
    incoming_energy: Total,
    incoming_offset: Total,
    surfaces: BTreeMap<CellId, SurfaceTotals>,
    phase_mass: Total,
    species_mass: BTreeMap<SpeciesId, Total>,
    external: Total,
    thermal_energy: Total,
    gas_energy: Total,
    coupled_energy: Total,
    coupled_mass: Total,
    absolute_energy: Total,
    absolute_mass: Total,
}
impl Account {
    fn add_events(
        &mut self,
        births: &BTreeMap<CellId, Birth>,
        audit: &EnergyAudit,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), AdaptiveImplicitGasError> {
        checkpoint()?;
        if audit.substeps != 0 {
            return Err(DepositionError::InvalidInput("event pass must not advance heat").into());
        }
        self.events = add_count(self.events, audit.events)?;
        self.incoming_mass.add(audit.incoming_mass_kg)?;
        self.incoming_energy.add(audit.incoming_energy_j)?;
        self.external.add(audit.external_heat_j)?;
        self.thermal_energy.add(audit.residual_j)?;
        self.coupled_energy.add(audit.residual_j)?;
        self.absolute_energy.add(audit.residual_j.abs())?;
        for (index, birth) in births.values().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            self.incoming_offset
                .add(finite(birth.mass_kg * birth.liquid_enthalpy_offset_j_kg)?)?;
        }
        Ok(())
    }

    fn add(
        &mut self,
        part: &ImplicitGasReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), AdaptiveImplicitGasError> {
        checkpoint()?;
        let native = &part.thermal.thermal;
        self.external.add(part.external_heat_j)?;
        self.thermal_energy.add(native.energy_residual_j)?;
        self.gas_energy.add(part.gas_energy_residual_j)?;
        self.coupled_energy.add(part.energy_residual_j)?;
        self.coupled_mass.add(part.mass_residual_kg)?;
        let mut local_energy = Total::default();
        let mut local_mass = Total::default();
        for (index, (&id, phase)) in native.surfaces.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            self.surfaces.entry(id).or_default().add(
                phase,
                part.thermal.vapor_enthalpy_correction_j[&id],
                part.sensible_heat_into_liquid_j[&id],
            )?;
            self.phase_mass.add(phase.transferred_mass_kg)?;
            local_energy.add(phase.energy_residual_j.abs())?;
            local_mass.add(phase.mass_rate_residual_kg.abs())?;
        }
        for (id, &mass) in &part.species_transferred_mass_kg {
            checkpoint()?;
            self.species_mass.entry(id.clone()).or_default().add(mass)?;
        }
        self.absolute_energy
            .add(native.energy_residual_j.abs().max(local_energy.value()?))?;
        self.absolute_energy.add(part.gas_energy_residual_j.abs())?;
        self.absolute_mass.add(
            native
                .mass_balance_residual_kg
                .abs()
                .max(local_mass.value()?),
        )?;
        self.absolute_mass.add(part.gas_mass_residual_kg)?;
        Ok(())
    }
    fn gate(&self, controls: AdaptiveImplicitGasControls) -> Result<(), AdaptiveImplicitGasError> {
        let residual_j = self.absolute_energy.value()?;
        if residual_j > controls.solve.thermal.total_energy_tolerance_j {
            return Err(ImplicitGasError::CoupledEnergyBalance {
                residual_j,
                tolerance_j: controls.solve.thermal.total_energy_tolerance_j,
            }
            .into());
        }
        let residual_kg = self.absolute_mass.value()?;
        let tolerance_kg = mass_tolerance(controls);
        if residual_kg > tolerance_kg {
            return Err(ImplicitGasError::GasMassBalance {
                residual_kg,
                tolerance_kg,
            }
            .into());
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn time_error(
    old: &DepositionSimulation<'_>,
    coarse: &DepositionSimulation<'_>,
    fine: &DepositionSimulation<'_>,
    selected: &BTreeSet<CellId>,
    old_gas: &ClosedGasParcel,
    coarse_gas: &ClosedGasParcel,
    fine_gas: &ClosedGasParcel,
    species: &BTreeSet<SpeciesId>,
    controls: AdaptiveImplicitGasControls,
    correction: f64,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<ImplicitGasTimeError, AdaptiveImplicitGasError> {
    let mut error = ImplicitGasTimeError::default();
    for (index, (&id, cell)) in old.cells().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let (difference, ratio) = component_error(
            cell.specific_h_j_kg,
            coarse.cells()[&id].specific_h_j_kg,
            fine.cells()[&id].specific_h_j_kg,
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
        if selected.contains(&id) {
            let (difference, ratio) = component_error(
                cell.mass_kg,
                coarse.cells()[&id].mass_kg,
                fine.cells()[&id].mass_kg,
                controls.mass_absolute_tolerance_kg,
                controls.relative_tolerance,
                correction,
            )?;
            error.maximum_liquid_mass_difference_kg =
                error.maximum_liquid_mass_difference_kg.max(difference);
            error.normalized = error.normalized.max(ratio);
        }
    }
    let (difference, ratio) = component_error(
        old_gas.state().state().temperature,
        coarse_gas.state().state().temperature,
        fine_gas.state().state().temperature,
        controls.temperature_absolute_tolerance_k,
        controls.relative_tolerance,
        correction,
    )?;
    error.gas_temperature_difference_k = difference;
    error.normalized = error.normalized.max(ratio);
    let old_inventory = old_gas.species_inventory_kg()?;
    let coarse_inventory = coarse_gas.species_inventory_kg()?;
    let fine_inventory = fine_gas.species_inventory_kg()?;
    for (index, id) in species.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let (difference, ratio) = component_error(
            old_inventory.get(id).copied().unwrap_or(0.0),
            coarse_inventory.get(id).copied().unwrap_or(0.0),
            fine_inventory.get(id).copied().unwrap_or(0.0),
            controls.mass_absolute_tolerance_kg,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_species_mass_difference_kg =
            error.maximum_species_mass_difference_kg.max(difference);
        error.normalized = error.normalized.max(ratio);
    }
    checkpoint()?;
    Ok(error)
}

struct MassBalance {
    gas: f64,
    vapor: f64,
    species: BTreeMap<SpeciesId, f64>,
    coupled: f64,
    carrier: f64,
}

#[allow(clippy::too_many_arguments)]
fn final_mass_balance(
    initial: &DepositionSimulation<'_>,
    final_thermal: &DepositionSimulation<'_>,
    initial_gas: &ClosedGasParcel,
    final_gas: &ClosedGasParcel,
    births: &BTreeMap<CellId, Birth>,
    species: &BTreeSet<SpeciesId>,
    transfers: &BTreeMap<SpeciesId, f64>,
    transfer: f64,
    tolerance_kg: f64,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<MassBalance, AdaptiveImplicitGasError> {
    let gas_change = finite(final_gas.mass_kg() - initial_gas.mass_kg())?;
    let mut coupled = Total::default();
    for (index, (id, cell)) in final_thermal.cells().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let initial_mass = initial
            .cells()
            .get(id)
            .map(|cell| cell.mass_kg)
            .or_else(|| births.get(id).map(|birth| birth.mass_kg))
            .ok_or(DepositionError::InvalidInput(
                "missing native cell mass baseline",
            ))?;
        coupled.add(finite(cell.mass_kg - initial_mass)?)?;
    }
    coupled.add(gas_change)?;
    let old = initial_gas.species_inventory_kg()?;
    let next = final_gas.species_inventory_kg()?;
    let mut carrier = 0.0_f64;
    let mut vapor = Total::default();
    let mut species_residuals = BTreeMap::new();
    for (index, id) in species.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let change =
            finite(next.get(id).copied().unwrap_or(0.0) - old.get(id).copied().unwrap_or(0.0))?;
        if let Some(&mass) = transfers.get(id) {
            let residual_kg = finite(change - mass)?;
            if residual_kg.abs() > tolerance_kg {
                return Err(ImplicitGasError::GasMassBalance {
                    residual_kg,
                    tolerance_kg,
                }
                .into());
            }
            vapor.add(residual_kg)?;
            species_residuals.insert(id.clone(), residual_kg);
        } else {
            carrier = carrier.max(change.abs());
        }
    }
    let mass = MassBalance {
        gas: finite(gas_change - transfer)?,
        vapor: vapor.value()?,
        species: species_residuals,
        coupled: coupled.value()?,
        carrier,
    };
    for residual_kg in [mass.gas, mass.vapor, mass.coupled, mass.carrier] {
        if residual_kg.abs() > tolerance_kg {
            return Err(ImplicitGasError::GasMassBalance {
                residual_kg,
                tolerance_kg,
            }
            .into());
        }
    }
    checkpoint()?;
    Ok(mass)
}

fn component_error(
    old: f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
    correction: f64,
) -> Result<(f64, f64), AdaptiveImplicitGasError> {
    let difference = finite(fine - coarse)?.abs();
    let change = finite(coarse - old)?.abs().max(finite(fine - old)?.abs());
    let allowance = finite(absolute + finite(relative * change)?)?;
    let ratio = (difference / allowance) * correction;
    if ratio.is_nan() || ratio < 0.0 || !positive(allowance) {
        return Err(DepositionError::Arithmetic.into());
    }
    Ok((difference, ratio))
}
fn mass_tolerance(controls: AdaptiveImplicitGasControls) -> f64 {
    controls
        .solve
        .thermal
        .mass_tolerance_kg
        .min(controls.solve.phase.gas.mass_tolerance_kg)
}
fn finite(value: f64) -> Result<f64, AdaptiveImplicitGasError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
fn add_count(value: usize, by: usize) -> Result<usize, AdaptiveImplicitGasError> {
    value
        .checked_add(by)
        .ok_or_else(|| DepositionError::Budget("adaptive finite-gas work counter").into())
}
