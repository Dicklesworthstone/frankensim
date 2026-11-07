//! Adaptive heat, supported deposition and finite-gas evolution on one clock.
//!
//! Compare a complete original interval with two half intervals and accept only
//! the conservative fine history. Both substrate and coating caloric increments,
//! coating masses, all gas-species inventories and temperatures enter the norm.
//! This resolves deposition from a genuinely bare substrate without a seed or
//! using substrate mass as the error scale for a much smaller deposit.
//!
//! This is a first-order local temporal estimate, not a global error bound or
//! a spatial coating model. The original fixed-area, common-temperature and
//! molecular/film assumptions remain. Only inventory refusals are retryable;
//! unsupported source/caloric states retain their original refusal.

use super::super::{
    GasTimelineAudit, HeatTotals, Sum, ThermalWork,
    adaptive::{CellTotals, GasTimelineAdaptiveAudit, GasTimelineAdaptiveSolve},
    finite,
};
use super::{
    GasCoatingTimelineCandidate, GasCoatingTimelineError, GasCoatingTimelineSurface,
    GasSupportedChamberSolve, ThermalAccommodation, advance_coatings_metered,
};
use crate::sublimation::gas::GasSurfaceError;
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::enthalpy::sublimation::coating::CoatingSimulation;
use fs_thermal::enthalpy::{CellId, DepositionError, Event, StepConfig};
use std::collections::BTreeMap;

/// Original physical failures and bounded adaptive-work refusals.
#[derive(Clone, Debug, PartialEq)]
pub enum GasCoatingAdaptiveError {
    Invalid(&'static str),
    Timeline(GasCoatingTimelineError),
    TrialBudget {
        maximum_trials: usize,
    },
    MinimumStep {
        requested_dt_s: f64,
        minimum_dt_s: f64,
    },
}
impl From<GasCoatingTimelineError> for GasCoatingAdaptiveError {
    fn from(error: GasCoatingTimelineError) -> Self {
        Self::Timeline(error)
    }
}
impl From<DepositionError> for GasCoatingAdaptiveError {
    fn from(error: DepositionError) -> Self {
        Self::Timeline(error.into())
    }
}
impl From<GasSurfaceError> for GasCoatingAdaptiveError {
    fn from(error: GasSurfaceError) -> Self {
        Self::Timeline(error.into())
    }
}
impl core::fmt::Display for GasCoatingAdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive supported deposition: {self:?}")
    }
}
impl std::error::Error for GasCoatingAdaptiveError {}

/// Original complete physical endpoint and its local-error/work account.
#[derive(Clone)]
pub struct GasCoatingAdaptiveCandidate<'thermal, 'catalog> {
    endpoint: GasCoatingTimelineCandidate<'thermal, 'catalog>,
    audit: GasTimelineAdaptiveAudit,
}
impl<'thermal, 'catalog> GasCoatingAdaptiveCandidate<'thermal, 'catalog> {
    /// Physical accounts cover the whole accepted window, including births.
    #[must_use]
    pub const fn endpoint(&self) -> &GasCoatingTimelineCandidate<'thermal, 'catalog> {
        &self.endpoint
    }
    /// Rejected/coarse work is counted here, never in the physical heat history.
    #[must_use]
    pub const fn audit(&self) -> &GasTimelineAdaptiveAudit {
        &self.audit
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasCoatingTimelineCandidate<'thermal, 'catalog>,
        GasTimelineAdaptiveAudit,
    ) {
        (self.endpoint, self.audit)
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn increment(value: &mut usize, count: usize) -> Result<(), DepositionError> {
    *value = value
        .checked_add(count)
        .ok_or(DepositionError::Arithmetic)?;
    Ok(())
}
fn ratio(a: f64, b: f64, scale: f64) -> Result<f64, GasCoatingAdaptiveError> {
    let error = (a - b).abs() / scale;
    if !positive(scale) || !error.is_finite() {
        return Err(GasCoatingAdaptiveError::Invalid(
            "unrepresentable coating error scale",
        ));
    }
    Ok(error)
}
fn scale(a: f64, b: f64, initial: f64, absolute: f64, relative: f64) -> f64 {
    absolute + relative * a.abs().max(b.abs()).max(initial.abs())
}
fn caloric_ratio(
    a: f64,
    b: f64,
    initial: f64,
    controls: GasTimelineAdaptiveSolve,
) -> Result<f64, GasCoatingAdaptiveError> {
    ratio(
        a,
        b,
        controls.absolute_specific_enthalpy_tolerance_j_kg
            + controls.relative_tolerance * (a - initial).abs().max((b - initial).abs()),
    )
}

fn compare(
    before: &GasCoatingTimelineCandidate<'_, '_>,
    coarse: &GasCoatingTimelineCandidate<'_, '_>,
    fine: &GasCoatingTimelineCandidate<'_, '_>,
    controls: GasTimelineAdaptiveSolve,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<f64, GasCoatingAdaptiveError> {
    let a = coarse.gas().state().state().temperature;
    let b = fine.gas().state().state().temperature;
    let mut error = ratio(
        a,
        b,
        scale(
            a,
            b,
            before.gas().state().state().temperature,
            controls.absolute_temperature_tolerance_k,
            controls.relative_tolerance,
        ),
    )?;
    let original = before.coating();
    let left = coarse.coating();
    let right = fine.coating();
    if left.time_s() != right.time_s()
        || left.cells().len() != right.cells().len()
        || left.coatings().len() != right.coatings().len()
    {
        return Err(GasCoatingAdaptiveError::Invalid(
            "coating trials disagree on time or inventory roster",
        ));
    }
    for (&id, a) in left.cells() {
        checkpoint()?;
        let b = right
            .cells()
            .get(&id)
            .ok_or(GasCoatingAdaptiveError::Invalid(
                "coating trial substrate identity changed",
            ))?;
        let initial = original.cells().get(&id).unwrap_or(a);
        if a.material_id != b.material_id
            || a.geometry != b.geometry
            || a.mass_kg != b.mass_kg
            || a.mass_kg != initial.mass_kg
        {
            return Err(GasCoatingAdaptiveError::Invalid(
                "coating trial changed original substrate",
            ));
        }
        let ta = left.temperature_k(id)?;
        let tb = right.temperature_k(id)?;
        let t0 = left.temperature_at_specific_enthalpy_k(id, initial.specific_h_j_kg)?;
        error = error.max(ratio(
            ta,
            tb,
            scale(
                ta,
                tb,
                t0,
                controls.absolute_temperature_tolerance_k,
                controls.relative_tolerance,
            ),
        )?);
        error = error.max(caloric_ratio(
            a.specific_h_j_kg,
            b.specific_h_j_kg,
            initial.specific_h_j_kg,
            controls,
        )?);
    }
    for (&id, a) in left.coatings() {
        checkpoint()?;
        let b = right
            .coatings()
            .get(&id)
            .ok_or(GasCoatingAdaptiveError::Invalid(
                "coating trial deposit identity changed",
            ))?;
        // Endpoint births have no preceding exposure in either trajectory.
        let initial = original.coatings().get(&id).unwrap_or(a);
        error = error.max(ratio(
            a.mass_kg,
            b.mass_kg,
            scale(
                a.mass_kg,
                b.mass_kg,
                initial.mass_kg,
                controls.absolute_mass_tolerance_kg,
                controls.relative_tolerance,
            ),
        )?);
        if a.mass_kg > 0.0 || b.mass_kg > 0.0 {
            error = error.max(caloric_ratio(
                a.specific_enthalpy_j_kg,
                b.specific_enthalpy_j_kg,
                initial.specific_enthalpy_j_kg,
                controls,
            )?);
        }
    }
    for species in fine.catalog().thermodynamics().components() {
        checkpoint()?;
        let mass = |gas: &ClosedGasParcel| {
            gas.species_mass_kg(species.species())
                .map_err(GasSurfaceError::from)
        };
        let a = mass(coarse.gas())?;
        let b = mass(fine.gas())?;
        error = error.max(ratio(
            a,
            b,
            scale(
                a,
                b,
                mass(before.gas())?,
                controls.absolute_mass_tolerance_kg,
                controls.relative_tolerance,
            ),
        )?);
    }
    Ok(error)
}

#[derive(Default)]
struct History {
    heat: HeatTotals,
    cells: BTreeMap<CellId, CellTotals>,
    outgoing: Sum,
    incoming: Sum,
    transferred: Sum,
    energy: Sum,
    reflected: Sum,
    stored: Sum,
    gas_stored: Sum,
    birth_energy: Sum,
    solid_defect: Sum,
    gas_defect: Sum,
    phase_steps: usize,
    gas_iterations: usize,
}
impl History {
    fn add(
        &mut self,
        candidate: &GasCoatingTimelineCandidate<'_, '_>,
    ) -> Result<(), DepositionError> {
        let report = candidate.audit();
        self.heat.add(report.thermal)?;
        self.outgoing.add(report.outgoing_mass_kg)?;
        self.incoming.add(report.incoming_mass_kg)?;
        self.transferred.add(report.transferred_mass_kg)?;
        self.energy.add(report.energy_into_gas_j)?;
        self.reflected.add(report.reflected_heat_into_gas_j)?;
        self.stored.add(report.solid_energy_change_j)?;
        self.gas_stored.add(report.gas_energy_change_j)?;
        self.birth_energy.add(report.incoming_energy_j)?;
        self.solid_defect
            .add(candidate.solid_stream_energy_residual_j)?;
        self.gas_defect
            .add(candidate.gas_stream_energy_residual_j)?;
        increment(&mut self.phase_steps, report.phase_steps)?;
        increment(&mut self.gas_iterations, report.gas_iterations)?;
        for (&id, row) in &report.cells {
            self.cells.entry(id).or_default().add(row)?;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        self,
        original: &CoatingSimulation<'_>,
        gas: &ClosedGasParcel,
        surfaces: &[GasCoatingTimelineSurface<'_>],
        endpoint: &mut GasCoatingTimelineCandidate<'_, '_>,
        heat: StepConfig,
        controls: GasSupportedChamberSolve,
        work: &mut GasTimelineAdaptiveAudit,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), GasCoatingAdaptiveError> {
        let mut births = BTreeMap::new();
        for event in original.pending_events() {
            checkpoint()?;
            if event.time_s > endpoint.coating().time_s() {
                break;
            }
            if let Event::Birth(cell) = &event.event {
                births.insert(cell.id, cell);
            }
        }
        for (&id, cell) in endpoint.coating().cells() {
            checkpoint()?;
            let before = original
                .cells()
                .get(&id)
                .or_else(|| births.get(&id).copied())
                .ok_or(GasCoatingAdaptiveError::Invalid(
                    "substrate has no original inventory or birth",
                ))?;
            if cell.mass_kg != before.mass_kg {
                return Err(GasCoatingAdaptiveError::Invalid(
                    "deposition changed nonvolatile substrate mass",
                ));
            }
        }
        let mut cells = BTreeMap::new();
        for (id, total) in self.cells {
            checkpoint()?;
            cells.insert(id, total.report()?);
        }
        let mut deposit_change = Sum::default();
        let mut species_change: BTreeMap<_, Sum> = BTreeMap::new();
        let mass_tolerance = controls.chamber.surface.gas.mass_tolerance_kg;
        let bindings: BTreeMap<_, _> = surfaces
            .iter()
            .map(|surface| (surface.cell_id, surface.binding))
            .collect();
        for (id, binding) in bindings {
            checkpoint()?;
            let Some(after) = endpoint.coating().coatings().get(&id) else {
                continue;
            };
            let before = original
                .coatings()
                .get(&id)
                .map_or(original.coating_surfaces()[&id].initial_mass_kg, |state| {
                    state.mass_kg
                });
            let dm = finite(after.mass_kg - before)?;
            let local = finite(dm + cells[&id].transferred_mass_kg)?;
            if local.abs() > mass_tolerance {
                return Err(GasCoatingTimelineError::MassBalance {
                    residual_kg: local,
                    tolerance_kg: mass_tolerance,
                }
                .into());
            }
            deposit_change.add(dm)?;
            species_change
                .entry(binding.species.clone())
                .or_default()
                .add(dm)?;
        }
        let mass_residual =
            finite(deposit_change.value()? + (endpoint.gas().mass_kg() - gas.mass_kg()))?;
        for source in endpoint.catalog().thermodynamics().components() {
            checkpoint()?;
            let id = source.species();
            let condensed = species_change.get(id).map_or(Ok(0.0), Sum::value)?;
            let change = endpoint
                .gas()
                .species_mass_kg(id)
                .map_err(GasSurfaceError::from)?
                - gas.species_mass_kg(id).map_err(GasSurfaceError::from)?;
            work.maximum_species_mass_residual_kg = work
                .maximum_species_mass_residual_kg
                .max(finite(condensed + change)?.abs());
        }
        let largest = mass_residual
            .abs()
            .max(work.maximum_species_mass_residual_kg);
        if largest > mass_tolerance {
            return Err(GasCoatingTimelineError::MassBalance {
                residual_kg: largest,
                tolerance_kg: mass_tolerance,
            }
            .into());
        }
        work.solid_stream_energy_residual_j = self.solid_defect.value()?;
        work.gas_stream_energy_residual_j = self.gas_defect.value()?;
        let energy_residual =
            finite(work.solid_stream_energy_residual_j + work.gas_stream_energy_residual_j)?;
        let energy_tolerance = heat
            .energy_tolerance_j
            .max(controls.chamber.surface.gas.energy_tolerance_j);
        for residual in [
            energy_residual,
            work.solid_stream_energy_residual_j,
            work.gas_stream_energy_residual_j,
        ] {
            if residual.abs() > energy_tolerance {
                return Err(GasSurfaceError::EnergyBalance {
                    residual_j: residual,
                    tolerance_j: energy_tolerance,
                }
                .into());
            }
        }
        endpoint.audit = GasTimelineAudit {
            thermal: self.heat.audit()?,
            cells,
            condensed_phases: BTreeMap::new(),
            outgoing_mass_kg: self.outgoing.value()?,
            incoming_mass_kg: self.incoming.value()?,
            transferred_mass_kg: self.transferred.value()?,
            energy_into_gas_j: self.energy.value()?,
            reflected_heat_into_gas_j: self.reflected.value()?,
            solid_energy_change_j: self.stored.value()?,
            gas_energy_change_j: self.gas_stored.value()?,
            incoming_energy_j: self.birth_energy.value()?,
            mass_residual_kg: mass_residual,
            energy_residual_j: energy_residual,
            phase_steps: self.phase_steps,
            gas_iterations: self.gas_iterations,
            last_step: endpoint.audit.last_step.take(),
        };
        endpoint.solid_stream_energy_residual_j = work.solid_stream_energy_residual_j;
        endpoint.gas_stream_energy_residual_j = work.gas_stream_energy_residual_j;
        Ok(())
    }
}

fn smaller(requested: f64, old: f64, minimum: f64) -> Result<f64, GasCoatingAdaptiveError> {
    let dt = requested.min(0.5 * old);
    if !positive(dt) || dt < minimum {
        Err(GasCoatingAdaptiveError::MinimumStep {
            requested_dt_s: dt,
            minimum_dt_s: minimum,
        })
    } else {
        Ok(dt)
    }
}

/// Adapt the COMPLETE heat/coating/gas split, with exactly one gas inventory.
///
/// Controls are shared with the original bulk timeline. Deposits use their own
/// mass scale; coating and substrate enthalpies use datum-independent increment
/// scales. Every trial ends at the next scheduled event. Original thermal CFL
/// substeps, including completed work in rejected probes, spend one whole-call
/// allowance. Empty pre-birth spans jump directly to their next event.
///
/// Only the two conservative half steps of an accepted trial enter physical
/// history. Reflection and finite-film mass resistance use the existing kernels.
/// Source errors, exhausted budgets and cancellation return no partial state.
#[allow(clippy::too_many_arguments)]
pub fn advance_coatings_adaptive_with_gas<'thermal, 'catalog>(
    coating: &CoatingSimulation<'thermal>,
    surfaces: &[GasCoatingTimelineSurface<'_>],
    accommodation: Option<&BTreeMap<CellId, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasSupportedChamberSolve,
    controls: GasTimelineAdaptiveSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasCoatingAdaptiveCandidate<'thermal, 'catalog>, GasCoatingAdaptiveError> {
    checkpoint()?;
    if !target_s.is_finite()
        || target_s < coating.time_s()
        || !positive(controls.absolute_mass_tolerance_kg)
        || !positive(controls.absolute_temperature_tolerance_k)
        || !positive(controls.absolute_specific_enthalpy_tolerance_j_kg)
        || !controls.relative_tolerance.is_finite()
        || controls.relative_tolerance < 0.0
        || !positive(controls.minimum_dt_s)
        || !positive(thermal_config.max_dt_s)
        || controls.minimum_dt_s > thermal_config.max_dt_s
        || controls.maximum_trials == 0
    {
        return Err(GasCoatingAdaptiveError::Invalid(
            "finite target, positive error scales and bounded coating-step controls required",
        ));
    }
    let meter = ThermalWork::new(coating.work_budget().max_substeps);
    let mut current = advance_coatings_metered(
        coating,
        surfaces,
        accommodation,
        gas,
        catalog,
        coating.time_s(),
        thermal_config,
        phase_config,
        Some(&meter),
        &mut checkpoint,
    )?;
    let mut history = History::default();
    history.add(&current)?;
    let mut work = GasTimelineAdaptiveAudit {
        trials: 0,
        rejected_trials: 0,
        accepted_intervals: 0,
        interval_evaluations: 1,
        thermal_substeps_evaluated: 0,
        dormant_advances: 0,
        maximum_accepted_error_ratio: 0.0,
        maximum_species_mass_residual_kg: 0.0,
        gas_stream_energy_residual_j: 0.0,
        solid_stream_energy_residual_j: 0.0,
    };
    let mut dt = thermal_config.max_dt_s;
    while current.coating().time_s() < target_s {
        checkpoint()?;
        let start = current.coating().time_s();
        let event_stop = current
            .coating()
            .pending_events()
            .first()
            .map_or(target_s, |event| event.time_s.min(target_s));
        if current.coating().cells().is_empty() {
            increment(&mut work.interval_evaluations, 1)?;
            current = advance_coatings_metered(
                current.coating(),
                surfaces,
                accommodation,
                current.gas(),
                catalog,
                event_stop,
                thermal_config,
                phase_config,
                Some(&meter),
                &mut checkpoint,
            )?;
            history.add(&current)?;
            increment(&mut work.dormant_advances, 1)?;
            continue;
        }
        if work.trials >= controls.maximum_trials {
            return Err(GasCoatingAdaptiveError::TrialBudget {
                maximum_trials: controls.maximum_trials,
            });
        }
        let stop = (start + dt).min(event_stop);
        let span = stop - start;
        let middle = start + 0.5 * span;
        if !positive(span) || !(start < middle && middle < stop) {
            return Err(GasCoatingAdaptiveError::MinimumStep {
                requested_dt_s: 0.5 * span,
                minimum_dt_s: controls.minimum_dt_s,
            });
        }
        increment(&mut work.trials, 1)?;
        let mut evaluate = |state: &GasCoatingTimelineCandidate<'thermal, 'catalog>, end: f64| {
            increment(&mut work.interval_evaluations, 1)?;
            advance_coatings_metered(
                state.coating(),
                surfaces,
                accommodation,
                state.gas(),
                catalog,
                end,
                StepConfig {
                    max_dt_s: end - state.coating().time_s(),
                    ..thermal_config
                },
                phase_config,
                Some(&meter),
                &mut checkpoint,
            )
        };
        let trial = (|| {
            let coarse = evaluate(&current, stop)?;
            let first = evaluate(&current, middle)?;
            let fine = evaluate(&first, stop)?;
            Ok::<_, GasCoatingTimelineError>((coarse, first, fine))
        })();
        let (coarse, first, fine) = match trial {
            Ok(value) => value,
            Err(GasCoatingTimelineError::Exchange(GasSurfaceError::InventoryLimit {
                maximum_duration_s,
            })) => {
                increment(&mut work.rejected_trials, 1)?;
                dt = smaller(0.8 * maximum_duration_s, span, controls.minimum_dt_s)?;
                continue;
            }
            Err(error) => return Err(error.into()),
        };
        let error = compare(&current, &coarse, &fine, controls, &mut checkpoint)?;
        if error > 1.0 {
            increment(&mut work.rejected_trials, 1)?;
            dt = smaller(
                span * (0.9 / fs_math::det::sqrt(error)).max(0.1),
                span,
                controls.minimum_dt_s,
            )?;
            continue;
        }
        history.add(&first)?;
        history.add(&fine)?;
        current = fine;
        increment(&mut work.accepted_intervals, 1)?;
        work.maximum_accepted_error_ratio = work.maximum_accepted_error_ratio.max(error);
        let growth = if error == 0.0 {
            2.0
        } else {
            (0.9 / fs_math::det::sqrt(error)).clamp(0.5, 2.0)
        };
        dt = (span * growth)
            .min(thermal_config.max_dt_s)
            .max(controls.minimum_dt_s);
    }
    history.finish(
        coating,
        gas,
        surfaces,
        &mut current,
        thermal_config,
        phase_config,
        &mut work,
        &mut checkpoint,
    )?;
    work.thermal_substeps_evaluated = meter.used();
    checkpoint()?;
    Ok(GasCoatingAdaptiveCandidate {
        endpoint: current,
        audit: work,
    })
}
