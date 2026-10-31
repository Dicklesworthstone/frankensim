//! Several retained liquid films on one original thermal/gas timeline.
//!
//! Heat transport advances once per event-clipped interval. The original joint
//! film owner then exchanges every exposed liquid against the SAME gas snapshot
//! and reconstructs one finite gas endpoint. Cells retain their own temperatures,
//! masses, source charts and fixed caloric datums. Opposing evaporation and
//! condensation therefore cannot replenish the old vapor withdrawal allowance.
//!
//! This is a first-order heat/phase split with frozen film rates. Step-size
//! refusals retry the complete interval; rejected attempts spend thermal work
//! but publish no state or physical account. Reference geometry stays fixed.
//! Selected cells must contain positive liquid mass on an explicitly admitted
//! branch. Dry-wall nucleation, film extinction and moving-volume work are not
//! inferred. Sensible heat through walls remains in the original heat network.
//! The opt-in gas-film path solves reciprocal heat for all exposed films before
//! phase transfer; its explicit coefficients exclude vapor-carried energy.

use super::{
    Branch, GasEvaporationTimelineError, GasEvaporationTimelineSurface, add, add_heat, finite,
};
use crate::evaporation::gas::{
    GasEvaporationError, GasFilmHeatSolve,
    chamber::{
        GasLiquidChamberSolve, GasLiquidChamberStep, GasLiquidChamberSurface,
        GasLiquidChamberSurfaceStep, advance_chamber_with_gas,
        heat::{
            GasLiquidChamberHeatSolve, GasLiquidChamberHeatStep, GasLiquidChamberHeatSurface,
            GasLiquidChamberSurfaceHeatStep, advance_chamber_with_gas_and_heat,
        },
    },
};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig, evaporation::LiquidCellExchange,
};
use fs_thermal::evaporation::{EvaporationError, LiquidSurfaceState};
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};

/// An original thermal film with explicit reciprocal gas-film heat controls.
/// Active heat must use exactly the same admitted liquid enthalpy branch.
#[derive(Clone, Copy)]
pub struct GasLiquidTimelineHeatSurface<'a> {
    pub surface: GasEvaporationTimelineSurface<'a>,
    pub heat_transfer_coefficient_w_m2_k: f64,
    pub heat_controls: GasFilmHeatSolve,
}

/// Accepted signed exchange on one physical film over the complete window.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GasLiquidCellTimelineAudit {
    pub transferred_mass_kg: f64,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    pub latent_energy_j: f64,
    pub energy_into_gas_j: f64,
    /// Sensible gas-to-liquid heat, separate from vapor-carried energy [J].
    pub sensible_heat_into_liquid_j: f64,
    pub recession_m: f64,
    pub phase_steps: usize,
}

/// Actual whole-window accounts. Latent energy is part of carried vapor energy.
#[derive(Clone, Debug, PartialEq)]
pub struct GasLiquidChamberTimelineAudit {
    pub thermal: EnergyAudit,
    pub cells: BTreeMap<CellId, GasLiquidCellTimelineAudit>,
    pub transferred_mass_kg: f64,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    pub latent_energy_j: f64,
    pub energy_into_gas_j: f64,
    /// Sensible gas-to-liquid heat, separate from vapor-carried energy [J].
    pub sensible_heat_into_liquid_j: f64,
    /// All original cells, including birth energy and selected liquid offsets.
    pub thermal_energy_change_j: f64,
    pub gas_energy_change_j: f64,
    pub incoming_energy_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Source-relative caloric defects, without subtracting formation energies.
    pub energy_residual_j: f64,
    /// Sum of absolute accepted film/gas defects and interval thermal defects.
    /// One whole-window energy budget bounds this, so cancellation cannot hide
    /// errors or renew the allowance for each additional surface or interval.
    pub sum_absolute_energy_residual_j: f64,
    pub phase_steps: usize,
    pub gas_iterations: usize,
    pub heat_iterations: usize,
    pub liquid_heat_iterations: usize,
    pub attempted_thermal_substeps: usize,
    pub rejected_intervals: usize,
    pub last_step: Option<GasLiquidChamberStep>,
    pub last_heat: Option<GasLiquidChamberHeatStep>,
}

/// Original heat and finite-gas owners publish together, or neither does.
#[derive(Clone)]
pub struct GasLiquidChamberTimelineCandidate<'thermal, 'catalog> {
    thermal: DepositionSimulation<'thermal>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    audit: GasLiquidChamberTimelineAudit,
}
impl<'thermal, 'catalog> GasLiquidChamberTimelineCandidate<'thermal, 'catalog> {
    pub const fn thermal(&self) -> &DepositionSimulation<'thermal> {
        &self.thermal
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn audit(&self) -> &GasLiquidChamberTimelineAudit {
        &self.audit
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'thermal>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasLiquidChamberTimelineAudit,
    ) {
        (self.thermal, self.gas, self.catalog, self.audit)
    }
}

#[derive(Clone)]
struct Branches(BTreeMap<CellId, Branch>);
impl DepositionObserver for Branches {
    type Error = DepositionError;
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        for branch in self.0.values_mut() {
            branch.event(event, thermal)?;
        }
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        for branch in self.0.values() {
            checkpoint()?;
            branch.check(thermal)?;
        }
        Ok(())
    }
}

struct FilmInterval {
    surfaces: Vec<GasLiquidChamberSurfaceStep>,
    gas: ClosedGasParcel,
    step: GasLiquidChamberStep,
    gas_iterations: usize,
    heat: Option<GasLiquidChamberHeatStep>,
    heat_surfaces: Vec<GasLiquidChamberSurfaceHeatStep>,
}

/// Advance every selected retained liquid film with one heat clock and gas.
///
/// Stable cell identity canonicalizes input order. Exposure is sampled before
/// each interval: endpoint covering cannot erase prior exchange, and newborn
/// films receive no exchange from before their birth. The original thermal
/// work allowance includes discarded attempts. The larger original heat/gas
/// energy tolerance bounds the complete window, not each surface separately.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_gas<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasEvaporationTimelineSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasLiquidChamberSolve,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasLiquidChamberTimelineCandidate<'thermal, 'catalog>, GasEvaporationTimelineError> {
    advance_deposition(
        thermal,
        surfaces,
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        None,
        checkpoint,
    )
}

/// Original wall heat, simultaneous reciprocal gas-film heat, then joint phase
/// transfer on the SAME event-clipped clock. Refine the complete timestep for
/// the first-order splitting error. A rejected interval discards all three
/// proposals together and spends the original thermal work allowance.
///
/// Zero coefficients preserve the phase-only route exactly; unused heat-solve
/// controls are ignored. Active film branches must match their original thermal
/// admission, including before-birth and zero-duration requests. One window
/// energy allowance covers wall, sensible and phase defects without cancellation.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_gas_and_heat<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasLiquidTimelineHeatSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasLiquidChamberSolve,
    heat_config: GasLiquidChamberHeatSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasLiquidChamberTimelineCandidate<'thermal, 'catalog>, GasEvaporationTimelineError> {
    checkpoint()?;
    if phase_config.maximum_surfaces == 0 || surfaces.len() > phase_config.maximum_surfaces {
        return Err(
            DepositionError::InvalidInput("multi-film heat surface budget exceeded").into(),
        );
    }
    let mut original = Vec::with_capacity(surfaces.len());
    let mut coefficients = BTreeMap::new();
    let mut enabled = false;
    for input in surfaces {
        checkpoint()?;
        let coefficient = input.heat_transfer_coefficient_w_m2_k;
        if !coefficient.is_finite() || coefficient < 0.0 {
            return Err(DepositionError::InvalidInput(
                "finite nonnegative film heat coefficient required",
            )
            .into());
        }
        if coefficient > 0.0 {
            enabled = true;
            crate::evaporation::gas::heat::validate_heat_controls(
                coefficient,
                input.heat_controls,
            )?;
            if input.heat_controls.pure_liquid_enthalpy_range_j_kg
                != input.surface.pure_liquid_enthalpy_range_j_kg
            {
                return Err(DepositionError::InvalidInput(
                    "joint sensible heat must retain each original liquid branch",
                )
                .into());
            }
        }
        if coefficients
            .insert(input.surface.cell_id, (coefficient, input.heat_controls))
            .is_some()
        {
            return Err(
                DepositionError::InvalidInput("duplicate thermal film heat identity").into(),
            );
        }
        original.push(input.surface);
    }
    if enabled && !(1..=256).contains(&heat_config.maximum_iterations) {
        return Err(
            DepositionError::InvalidInput("joint film heat iterations must be in 1..=256").into(),
        );
    }
    advance_deposition(
        thermal,
        &original,
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        enabled.then_some((&coefficients, heat_config)),
        checkpoint,
    )
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance_deposition<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasEvaporationTimelineSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasLiquidChamberSolve,
    film_heat: Option<(
        &BTreeMap<CellId, (f64, GasFilmHeatSolve)>,
        GasLiquidChamberHeatSolve,
    )>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasLiquidChamberTimelineCandidate<'thermal, 'catalog>, GasEvaporationTimelineError> {
    checkpoint()?;
    if !target_s.is_finite()
        || target_s < thermal.time_s()
        || !thermal_config.max_dt_s.is_finite()
        || thermal_config.max_dt_s <= 0.0
        || !thermal_config.energy_tolerance_j.is_finite()
        || thermal_config.energy_tolerance_j < 0.0
        || phase_config.maximum_surfaces == 0
        || surfaces.len() > phase_config.maximum_surfaces
    {
        return Err(DepositionError::InvalidInput(
            "invalid multi-film interval or declared work/energy budget",
        )
        .into());
    }
    let mut known: BTreeSet<_> = thermal.cells().keys().copied().collect();
    let mut births = BTreeMap::new();
    for (index, event) in thermal.pending_events().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        if let Event::Birth(cell) = &event.event {
            known.insert(cell.id);
            if event.time_s <= target_s {
                births.entry(cell.id).or_insert(cell);
            }
        }
    }
    let mut ordered = BTreeMap::new();
    let mut branches = Branches(BTreeMap::new());
    for surface in surfaces {
        checkpoint()?;
        let id = surface.cell_id;
        let [lo, hi] = surface.pure_liquid_enthalpy_range_j_kg;
        if !known.contains(&id)
            || !lo.is_finite()
            || !hi.is_finite()
            || lo > hi
            || !surface.binding.liquid_enthalpy_offset_j_kg.is_finite()
            || ordered.insert(id, surface).is_some()
        {
            return Err(DepositionError::InvalidInput(
                "unknown/duplicate liquid cell or invalid caloric branch/datum",
            )
            .into());
        }
        crate::evaporation::gas::validate_controls(
            phase_config.surface,
            target_s - thermal.time_s(),
            surface.mass_transfer_coefficient_m_s,
        )?;
        let branch = Branch {
            id,
            range: [lo, hi],
            born_mass_kg: 0.0,
        };
        branch.check(thermal)?;
        branches.0.insert(id, branch);
    }
    let mut checkpoint_error = None;
    let admission = advance_chamber_with_gas(&[], gas, catalog, 0.0, phase_config, || {
        if let Err(error) = checkpoint() {
            checkpoint_error = Some(error);
            true
        } else {
            false
        }
    });
    if let Some(error) = checkpoint_error {
        return Err(error.into());
    }
    let (_, mut next_gas, _, _) = admission?.into_parts();
    let (mut next, mut branches, initial_heat) =
        thermal.advance_observed(&branches, thermal.time_s(), thermal_config, &mut checkpoint)?;
    let mut audit = GasLiquidChamberTimelineAudit {
        thermal: initial_heat,
        cells: ordered
            .keys()
            .map(|id| (*id, GasLiquidCellTimelineAudit::default()))
            .collect(),
        transferred_mass_kg: 0.0,
        evaporated_mass_kg: 0.0,
        condensed_mass_kg: 0.0,
        latent_energy_j: 0.0,
        energy_into_gas_j: 0.0,
        sensible_heat_into_liquid_j: 0.0,
        thermal_energy_change_j: 0.0,
        gas_energy_change_j: 0.0,
        incoming_energy_j: 0.0,
        mass_residual_kg: 0.0,
        maximum_species_mass_residual_kg: 0.0,
        energy_residual_j: initial_heat.residual_j,
        sum_absolute_energy_residual_j: initial_heat.residual_j.abs(),
        phase_steps: 0,
        gas_iterations: 0,
        heat_iterations: 0,
        liquid_heat_iterations: 0,
        attempted_thermal_substeps: initial_heat.substeps,
        rejected_intervals: 0,
        last_step: None,
        last_heat: None,
    };
    let mut epoch = next.time_s();
    let mut interval = 1usize;
    while next.time_s() < target_s {
        checkpoint()?;
        let start = next.time_s();
        let mut exposed = Vec::with_capacity(ordered.len());
        for surface in ordered.values() {
            checkpoint()?;
            if next.cells().contains_key(&surface.cell_id) {
                exposed.push((**surface, next.exposed_area_m2(surface.cell_id)?));
            }
        }
        let empty = next.cells().is_empty();
        let event_time = next.pending_events().first().map(|event| event.time_s);
        let mut stop = if empty {
            target_s
        } else {
            finite(epoch + (interval as f64) * thermal_config.max_dt_s)?.min(target_s)
        };
        if let Some(event_time) = event_time {
            stop = stop.min(event_time);
        }
        if !stop.is_finite() || stop <= start {
            return Err(DepositionError::Arithmetic.into());
        }
        let planned_stop = stop;
        loop {
            checkpoint()?;
            let remaining = next
                .work_budget()
                .max_substeps
                .checked_sub(audit.attempted_thermal_substeps)
                .ok_or(DepositionError::Budget("multi-film thermal substeps"))?;
            let (mut heated, observed, heat) = next.advance_observed_limited(
                &branches,
                stop,
                thermal_config,
                remaining,
                &mut checkpoint,
            )?;
            audit.attempted_thermal_substeps = audit
                .attempted_thermal_substeps
                .checked_add(heat.substeps)
                .ok_or(DepositionError::Arithmetic)?;
            if !exposed.is_empty() {
                let chart_error = RefCell::new(None);
                let calorics: Vec<_> = exposed
                    .iter()
                    .map(|(surface, _)| {
                        let state = &heated;
                        let errors = &chart_error;
                        let id = surface.cell_id;
                        let [lo, hi] = surface.pure_liquid_enthalpy_range_j_kg;
                        move |h| {
                            let value = if (lo..=hi).contains(&h) {
                                state.temperature_at_specific_enthalpy_k(id, h)
                            } else {
                                Err(DepositionError::MaterialDomain)
                            };
                            value.map_err(|error| {
                                *errors.borrow_mut() = Some(error);
                                EvaporationError::Invalid("original liquid caloric branch refused")
                            })
                        }
                    })
                    .collect();
                let inputs: Vec<_> = exposed
                    .iter()
                    .zip(&calorics)
                    .map(|((surface, area), caloric)| {
                        let cell = &heated.cells()[&surface.cell_id];
                        GasLiquidChamberSurface {
                            id: surface.cell_id.0,
                            law: surface.law,
                            binding: surface.binding,
                            liquid: LiquidSurfaceState {
                                mass_kg: cell.mass_kg,
                                specific_enthalpy_j_kg: cell.specific_h_j_kg,
                                area_m2: *area,
                            },
                            temperature_from_enthalpy: caloric,
                            mass_transfer_coefficient_m_s: surface.mass_transfer_coefficient_m_s,
                        }
                    })
                    .collect();
                let mut checkpoint_error = None;
                let cancelled = || {
                    if let Err(error) = checkpoint() {
                        checkpoint_error = Some(error);
                        true
                    } else {
                        false
                    }
                };
                let proposal = if let Some((coefficients, heat_config)) = film_heat {
                    let heat_inputs: Vec<_> = inputs
                        .iter()
                        .map(|surface| {
                            let (coefficient, heat_controls) = coefficients[&CellId(surface.id)];
                            GasLiquidChamberHeatSurface {
                                surface: *surface,
                                heat_transfer_coefficient_w_m2_k: coefficient,
                                heat_controls,
                            }
                        })
                        .collect();
                    advance_chamber_with_gas_and_heat(
                        &heat_inputs,
                        &next_gas,
                        catalog,
                        stop - start,
                        phase_config,
                        heat_config,
                        cancelled,
                    )
                    .map(|candidate| {
                        let gas_iterations = candidate.phase().gas_exchange().step().iterations;
                        let (surfaces, gas, _, step, heat, heat_surfaces) = candidate.into_parts();
                        FilmInterval {
                            surfaces,
                            gas,
                            step,
                            gas_iterations,
                            heat: Some(heat),
                            heat_surfaces,
                        }
                    })
                } else {
                    advance_chamber_with_gas(
                        &inputs,
                        &next_gas,
                        catalog,
                        stop - start,
                        phase_config,
                        cancelled,
                    )
                    .map(|candidate| {
                        let gas_iterations = candidate.gas_exchange().step().iterations;
                        let (surfaces, gas, _, step) = candidate.into_parts();
                        FilmInterval {
                            surfaces,
                            gas,
                            step,
                            gas_iterations,
                            heat: None,
                            heat_surfaces: Vec::new(),
                        }
                    })
                };
                drop(inputs);
                drop(calorics);
                if let Some(error) = checkpoint_error {
                    return Err(error.into());
                }
                let candidate = match proposal {
                    Ok(candidate) => candidate,
                    Err(error) => {
                        if let Some(original) = chart_error.into_inner() {
                            return Err(original.into());
                        }
                        if let GasEvaporationError::StepSizeLimit {
                            maximum_duration_s, ..
                        } = error
                        {
                            let reduced = maximum_duration_s.min(0.5 * (stop - start));
                            let retry_stop = start + reduced;
                            if maximum_duration_s.is_finite()
                                && reduced > 0.0
                                && retry_stop.is_finite()
                                && retry_stop > start
                                && retry_stop < stop
                                && heat.substeps > 0
                            {
                                audit.rejected_intervals = audit
                                    .rejected_intervals
                                    .checked_add(1)
                                    .ok_or(DepositionError::Arithmetic)?;
                                stop = retry_stop;
                                continue;
                            }
                        }
                        return Err(error.into());
                    }
                };
                if candidate.surfaces.len() != exposed.len() {
                    return Err(DepositionError::InvalidInput("joint film count changed").into());
                }
                let mut updates = Vec::with_capacity(exposed.len());
                for ((surface, area), endpoint) in exposed.iter().zip(&candidate.surfaces) {
                    checkpoint()?;
                    if endpoint.id != surface.cell_id.0 {
                        return Err(
                            DepositionError::InvalidInput("joint film identity changed").into()
                        );
                    }
                    updates.push(LiquidCellExchange {
                        expected: heated.cells()[&surface.cell_id].clone(),
                        interval_exposed_area_m2: *area,
                        replacement: endpoint.liquid,
                        pure_liquid_enthalpy_range_j_kg: surface.pure_liquid_enthalpy_range_j_kg,
                    });
                }
                heated.apply_liquid_exchanges(&updates, &mut checkpoint)?;
                audit.gas_iterations = audit
                    .gas_iterations
                    .checked_add(candidate.gas_iterations)
                    .ok_or(DepositionError::Arithmetic)?;
                let FilmInterval {
                    surfaces: endpoints,
                    gas: evolved_gas,
                    step,
                    heat: sensible_heat,
                    heat_surfaces,
                    ..
                } = candidate;
                if let Some(sensible) = sensible_heat {
                    add(
                        &mut audit.sensible_heat_into_liquid_j,
                        sensible.heat_into_liquid_j,
                    )?;
                    audit.heat_iterations = audit
                        .heat_iterations
                        .checked_add(sensible.iterations)
                        .ok_or(DepositionError::Arithmetic)?;
                    audit.liquid_heat_iterations = audit
                        .liquid_heat_iterations
                        .checked_add(sensible.liquid_iterations)
                        .ok_or(DepositionError::Arithmetic)?;
                    for receipt in heat_surfaces {
                        let account = audit.cells.get_mut(&CellId(receipt.id)).ok_or(
                            DepositionError::InvalidInput("unselected sensible-heat film endpoint"),
                        )?;
                        add(
                            &mut account.sensible_heat_into_liquid_j,
                            receipt.heat_into_liquid_j,
                        )?;
                    }
                    audit.last_heat = Some(sensible);
                }
                next_gas = evolved_gas;
                add(&mut audit.transferred_mass_kg, step.transferred_mass_kg)?;
                add(&mut audit.evaporated_mass_kg, step.evaporated_mass_kg)?;
                add(&mut audit.condensed_mass_kg, step.condensed_mass_kg)?;
                add(&mut audit.latent_energy_j, step.latent_energy_j)?;
                add(&mut audit.energy_into_gas_j, step.energy_into_gas_j)?;
                add(
                    &mut audit.thermal_energy_change_j,
                    step.liquid_enthalpy_change_j,
                )?;
                add(
                    &mut audit.gas_energy_change_j,
                    step.gas_internal_energy_change_j,
                )?;
                add(&mut audit.energy_residual_j, step.energy_residual_j)?;
                add(
                    &mut audit.sum_absolute_energy_residual_j,
                    step.sum_absolute_energy_residual_j,
                )?;
                audit.phase_steps = audit
                    .phase_steps
                    .checked_add(1)
                    .ok_or(DepositionError::Arithmetic)?;
                for endpoint in endpoints {
                    let account = audit.cells.get_mut(&CellId(endpoint.id)).ok_or(
                        DepositionError::InvalidInput("unselected joint film endpoint"),
                    )?;
                    add(
                        &mut account.transferred_mass_kg,
                        endpoint.transferred_mass_kg,
                    )?;
                    add(
                        &mut account.evaporated_mass_kg,
                        endpoint.transferred_mass_kg.max(0.0),
                    )?;
                    add(
                        &mut account.condensed_mass_kg,
                        (-endpoint.transferred_mass_kg).max(0.0),
                    )?;
                    add(&mut account.latent_energy_j, endpoint.latent_energy_j)?;
                    add(&mut account.energy_into_gas_j, endpoint.energy_into_gas_j)?;
                    add(&mut account.recession_m, endpoint.recession_m)?;
                    account.phase_steps = account
                        .phase_steps
                        .checked_add(1)
                        .ok_or(DepositionError::Arithmetic)?;
                }
                audit.last_step = Some(step);
            }
            next = heated;
            branches = observed;
            add_heat(&mut audit.thermal, heat)?;
            add(&mut audit.energy_residual_j, heat.residual_j)?;
            add(
                &mut audit.sum_absolute_energy_residual_j,
                heat.residual_j.abs(),
            )?;
            break;
        }
        if empty || event_time == Some(stop) {
            epoch = stop;
            interval = 1;
        } else if stop == planned_stop {
            interval = interval.checked_add(1).ok_or(DepositionError::Arithmetic)?;
        }
    }
    // Compare each actual retained cell with its original state/birth, before
    // summing. A massive unrelated wall must not hide tiny film mass changes.
    let mut liquid_mass_change = 0.0;
    let mut species_changes = BTreeMap::new();
    for (&id, cell) in next.cells() {
        checkpoint()?;
        let original = thermal
            .cells()
            .get(&id)
            .or_else(|| births.get(&id).copied())
            .ok_or(DepositionError::InvalidInput(
                "thermal endpoint has no original cell or birth",
            ))?;
        let change = finite(cell.mass_kg - original.mass_kg)?;
        add(&mut liquid_mass_change, change)?;
        if let Some(surface) = ordered.get(&id) {
            add(
                species_changes
                    .entry(surface.binding.species.clone())
                    .or_insert(0.0),
                change,
            )?;
        } else if change != 0.0 {
            return Err(
                DepositionError::InvalidInput("unselected thermal cell changed mass").into(),
            );
        }
    }
    audit.mass_residual_kg = finite(liquid_mass_change + (next_gas.mass_kg() - gas.mass_kg()))?;
    for source in catalog.thermodynamics().components() {
        checkpoint()?;
        let id = source.species();
        let residual = finite(
            species_changes.get(id).copied().unwrap_or(0.0)
                + (next_gas
                    .species_mass_kg(id)
                    .map_err(GasEvaporationError::from)?
                    - gas.species_mass_kg(id).map_err(GasEvaporationError::from)?),
        )?;
        audit.maximum_species_mass_residual_kg =
            audit.maximum_species_mass_residual_kg.max(residual.abs());
    }
    let mass_residual = audit
        .mass_residual_kg
        .abs()
        .max(audit.maximum_species_mass_residual_kg);
    if mass_residual > phase_config.surface.gas.mass_tolerance_kg {
        return Err(GasEvaporationTimelineError::MassBalance {
            residual_kg: mass_residual,
            tolerance_kg: phase_config.surface.gas.mass_tolerance_kg,
        });
    }
    let mut birth_offset = 0.0;
    for (&id, branch) in &branches.0 {
        add(
            &mut birth_offset,
            branch.born_mass_kg * ordered[&id].binding.liquid_enthalpy_offset_j_kg,
        )?;
    }
    add(
        &mut audit.thermal_energy_change_j,
        audit.thermal.stored_energy_change_j,
    )?;
    add(&mut audit.thermal_energy_change_j, birth_offset)?;
    audit.incoming_energy_j = finite(audit.thermal.incoming_energy_j + birth_offset)?;
    let tolerance = thermal_config
        .energy_tolerance_j
        .max(phase_config.surface.gas.energy_tolerance_j);
    for residual in [
        audit.energy_residual_j,
        audit.sum_absolute_energy_residual_j,
    ] {
        if residual.abs() > tolerance {
            return Err(GasEvaporationError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
    }
    checkpoint()?;
    Ok(GasLiquidChamberTimelineCandidate {
        thermal: next,
        gas: next_gas,
        catalog,
        audit,
    })
}
