//! Heat transport, evaporation and one finite atmosphere on the same clock.
//!
//! Original heat/contact/source laws act first, then the existing conservative
//! liquid/gas owner advances the same interval. The next interval sees its
//! actual vapor inventory, pressure and temperatures. Reference geometry stays
//! fixed. This first-order split requires full-timestep refinement. Explicit
//! inventory/temperature/equilibrium step limits retry from the same accepted
//! state; rejected attempts spend work but contribute no physical exchange.

use super::{
    GasEvaporationError, GasEvaporationSolve, GasEvaporationStep, GasLiquidBinding,
    advance_surface_with_gas,
};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig, ThermalStorageBasis, evaporation::LiquidCellExchange,
};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, LiquidSurfaceState};
use std::collections::BTreeMap;

/// Explicit original cell, constitutive identity and admitted pure-liquid branch.
#[derive(Clone, Copy)]
pub struct GasEvaporationTimelineSurface<'a> {
    pub cell_id: CellId,
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    pub pure_liquid_enthalpy_range_j_kg: [f64; 2],
    pub mass_transfer_coefficient_m_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum GasEvaporationTimelineError {
    Thermal(DepositionError),
    Exchange(GasEvaporationError),
    MassBalance { residual_kg: f64, tolerance_kg: f64 },
}
impl From<DepositionError> for GasEvaporationTimelineError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<GasEvaporationError> for GasEvaporationTimelineError {
    fn from(error: GasEvaporationError) -> Self {
        Self::Exchange(error)
    }
}
impl core::fmt::Display for GasEvaporationTimelineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal/finite-gas evaporation: {self:?}")
    }
}
impl std::error::Error for GasEvaporationTimelineError {}

/// Complete requested interval, with local source-relative defects accumulated
/// BEFORE adding large conserved formation references to diagnostic totals.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEvaporationTimelineAudit {
    pub thermal: EnergyAudit,
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    pub energy_into_gas_j: f64,
    /// All original thermal cells, including births and the liquid datum offset.
    pub thermal_energy_change_j: f64,
    pub gas_energy_change_j: f64,
    /// Birth energy in the same gas/liquid datum as the storage reports.
    pub incoming_energy_j: f64,
    pub mass_residual_kg: f64,
    pub energy_residual_j: f64,
    pub recession_m: f64,
    pub phase_steps: usize,
    /// Native thermal work, INCLUDING discarded heat/phase attempts.
    pub attempted_thermal_substeps: usize,
    /// Explicit step-limit refusals that shortened an attempted interval.
    pub rejected_intervals: usize,
    pub last_step: Option<GasEvaporationStep>,
}

/// Both original owners publish together, or neither does.
#[derive(Clone)]
pub struct GasEvaporationTimelineCandidate<'thermal, 'catalog> {
    thermal: DepositionSimulation<'thermal>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    audit: GasEvaporationTimelineAudit,
}
impl<'thermal, 'catalog> GasEvaporationTimelineCandidate<'thermal, 'catalog> {
    pub const fn thermal(&self) -> &DepositionSimulation<'thermal> {
        &self.thermal
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn audit(&self) -> &GasEvaporationTimelineAudit {
        &self.audit
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'thermal>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasEvaporationTimelineAudit,
    ) {
        (self.thermal, self.gas, self.catalog, self.audit)
    }
}

fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn add(a: &mut f64, b: f64) -> Result<(), DepositionError> {
    *a = finite(*a + b)?;
    Ok(())
}
fn add_heat(total: &mut EnergyAudit, part: EnergyAudit) -> Result<(), DepositionError> {
    add(&mut total.incoming_mass_kg, part.incoming_mass_kg)?;
    add(&mut total.incoming_energy_j, part.incoming_energy_j)?;
    add(&mut total.external_heat_j, part.external_heat_j)?;
    add(
        &mut total.stored_energy_change_j,
        part.stored_energy_change_j,
    )?;
    add(&mut total.residual_j, part.residual_j)?;
    total.substeps = total
        .substeps
        .checked_add(part.substeps)
        .ok_or(DepositionError::Arithmetic)?;
    total.events = total
        .events
        .checked_add(part.events)
        .ok_or(DepositionError::Arithmetic)?;
    Ok(())
}

#[derive(Clone, Copy)]
struct Branch {
    id: CellId,
    range: [f64; 2],
    born_mass_kg: f64,
}
impl Branch {
    fn check(&self, thermal: &DepositionSimulation<'_>) -> Result<(), DepositionError> {
        if let Some(cell) = thermal.cells().get(&self.id) {
            if thermal.material_storage_basis(cell.material_id)? != ThermalStorageBasis::Enthalpy {
                return Err(DepositionError::InvalidInput(
                    "finite-gas evaporation requires original liquid enthalpy storage",
                ));
            }
            if !(self.range[0]..=self.range[1]).contains(&cell.specific_h_j_kg) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        Ok(())
    }
}
impl DepositionObserver for Branch {
    type Error = DepositionError;
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.check(thermal)?;
        if let Event::Birth(cell) = event {
            if cell.id == self.id {
                add(&mut self.born_mass_kg, cell.mass_kg)?;
            }
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
        checkpoint()?;
        self.check(thermal)
    }
}

/// Propose one selected liquid cell, original heat network and finite gas.
///
/// Initial exposure applies through each event-clipped interval; endpoint
/// covering cannot rewrite it, and a newborn cell receives no earlier exchange.
/// The original thermal substep allowance is spent once over the whole call.
/// The gas owner's absolute mass budget bounds the COMPLETE mass residual.
/// The larger original absolute heat/gas energy budget bounds the COMPLETE
/// source-relative energy defect, with no multiplication by accepted step count.
/// Explicit phase step limits discard the complete attempted heat/gas state
/// and retry a strictly shorter interval. Rejected thermal substeps spend the
/// same original work allowance; unsupported sources/charts still refuse.
///
/// Heat to other cells remains in the original thermal network. The gas-volume
/// and film assumptions of advance_surface_with_gas remain unchanged; no
/// gas/liquid sensible-conduction, moving volume or mechanical law is inferred.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_gas<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surface: GasEvaporationTimelineSurface<'_>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasEvaporationSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<GasEvaporationTimelineCandidate<'thermal, 'catalog>, GasEvaporationTimelineError> {
    checkpoint()?;
    let id = surface.cell_id;
    let [lo, hi] = surface.pure_liquid_enthalpy_range_j_kg;
    if !target_s.is_finite()
        || target_s < thermal.time_s()
        || !thermal_config.max_dt_s.is_finite()
        || thermal_config.max_dt_s <= 0.0
        || !thermal_config.energy_tolerance_j.is_finite()
        || thermal_config.energy_tolerance_j < 0.0
        || !lo.is_finite()
        || !hi.is_finite()
        || lo > hi
        || !surface.mass_transfer_coefficient_m_s.is_finite()
        || surface.mass_transfer_coefficient_m_s < 0.0
        || !surface.binding.liquid_enthalpy_offset_j_kg.is_finite()
        || !phase_config.gas.mass_tolerance_kg.is_finite()
        || phase_config.gas.mass_tolerance_kg <= 0.0
        || !phase_config.gas.energy_tolerance_j.is_finite()
        || phase_config.gas.energy_tolerance_j <= 0.0
    {
        return Err(DepositionError::InvalidInput(
            "invalid finite-gas evaporation timeline, liquid branch or balance budget",
        )
        .into());
    }
    super::validate_controls(
        phase_config,
        target_s - thermal.time_s(),
        surface.mass_transfer_coefficient_m_s,
    )?;
    let mut checkpoint_error = None;
    let admission = gas.exchange_species_boundary_relative(catalog, &[], phase_config.gas, || {
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
    let (mut next_gas, _, _) = admission.map_err(GasEvaporationError::from)?.into_parts();
    let mut known = thermal.cells().contains_key(&id);
    for (index, event) in thermal.pending_events().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        known |= matches!(&event.event, Event::Birth(cell) if cell.id == id);
    }
    if !known {
        return Err(
            DepositionError::InvalidInput("liquid cell is neither active nor scheduled").into(),
        );
    }
    let branch = Branch {
        id,
        range: [lo, hi],
        born_mass_kg: 0.0,
    };
    branch.check(thermal)?;
    let (mut next, mut branch, heat) =
        thermal.advance_observed(&branch, thermal.time_s(), thermal_config, &mut checkpoint)?;
    let mut audit = GasEvaporationTimelineAudit {
        thermal: heat,
        transferred_mass_kg: 0.0,
        latent_energy_j: 0.0,
        energy_into_gas_j: 0.0,
        thermal_energy_change_j: 0.0,
        gas_energy_change_j: 0.0,
        incoming_energy_j: 0.0,
        mass_residual_kg: 0.0,
        energy_residual_j: 0.0,
        recession_m: 0.0,
        phase_steps: 0,
        attempted_thermal_substeps: heat.substeps,
        rejected_intervals: 0,
        last_step: None,
    };
    let mut epoch = next.time_s();
    let mut interval = 1usize;
    while next.time_s() < target_s {
        checkpoint()?;
        let start = next.time_s();
        let area = if next.cells().contains_key(&id) {
            Some(next.exposed_area_m2(id)?)
        } else {
            None
        };
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
                .ok_or(DepositionError::Budget(
                    "finite-gas evaporation thermal substeps",
                ))?;
            let (mut heated, observed, heat) = next.advance_observed_limited(
                &branch,
                stop,
                thermal_config,
                remaining,
                &mut checkpoint,
            )?;
            audit.attempted_thermal_substeps = audit
                .attempted_thermal_substeps
                .checked_add(heat.substeps)
                .ok_or(DepositionError::Arithmetic)?;
            if let Some(area) = area {
                let expected = heated.cells()[&id].clone();
                let chart_error = std::cell::RefCell::new(None);
                let mut checkpoint_error = None;
                let proposal = advance_surface_with_gas(
                    surface.law,
                    surface.binding,
                    LiquidSurfaceState {
                        mass_kg: expected.mass_kg,
                        specific_enthalpy_j_kg: expected.specific_h_j_kg,
                        area_m2: area,
                    },
                    &next_gas,
                    catalog,
                    |h| {
                        let value = if (lo..=hi).contains(&h) {
                            heated.temperature_at_specific_enthalpy_k(id, h)
                        } else {
                            Err(DepositionError::MaterialDomain)
                        };
                        value.map_err(|error| {
                            *chart_error.borrow_mut() = Some(error);
                            EvaporationError::Invalid("original pure-liquid caloric branch refused")
                        })
                    },
                    surface.mass_transfer_coefficient_m_s,
                    stop - start,
                    phase_config,
                    || {
                        if let Err(error) = checkpoint() {
                            checkpoint_error = Some(error);
                            true
                        } else {
                            false
                        }
                    },
                );
                if let Some(error) = checkpoint_error {
                    return Err(error.into());
                }
                let candidate = match proposal {
                    Ok(value) => value,
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
                            if !maximum_duration_s.is_finite()
                                || reduced <= 0.0
                                || !retry_stop.is_finite()
                                || retry_stop <= start
                                || retry_stop >= stop
                                || heat.substeps == 0
                            {
                                return Err(error.into());
                            }
                            audit.rejected_intervals = audit
                                .rejected_intervals
                                .checked_add(1)
                                .ok_or(DepositionError::Arithmetic)?;
                            stop = retry_stop;
                            // `next`, `branch`, gas, birth offsets and physical
                            // audits still describe the unchanged accepted start.
                            continue;
                        }
                        return Err(error.into());
                    }
                };
                let (liquid, evolved_gas, _, step) = candidate.into_parts();
                heated.apply_liquid_exchanges(
                    &[LiquidCellExchange {
                        expected,
                        interval_exposed_area_m2: area,
                        replacement: liquid,
                        pure_liquid_enthalpy_range_j_kg: [lo, hi],
                    }],
                    &mut checkpoint,
                )?;
                next_gas = evolved_gas;
                add(&mut audit.transferred_mass_kg, step.transferred_mass_kg)?;
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
                add(&mut audit.recession_m, step.recession_m)?;
                audit.phase_steps = audit
                    .phase_steps
                    .checked_add(1)
                    .ok_or(DepositionError::Arithmetic)?;
                audit.last_step = Some(step);
            }
            next = heated;
            branch = observed;
            add_heat(&mut audit.thermal, heat)?;
            break;
        }
        if empty || event_time == Some(stop) {
            epoch = stop;
            interval = 1;
        } else if stop == planned_stop {
            interval = interval.checked_add(1).ok_or(DepositionError::Arithmetic)?;
        }
    }
    audit.mass_residual_kg = finite(
        (next.mass_kg() - thermal.mass_kg()) - audit.thermal.incoming_mass_kg
            + (next_gas.mass_kg() - gas.mass_kg()),
    )?;
    if audit.mass_residual_kg.abs() > phase_config.gas.mass_tolerance_kg {
        return Err(GasEvaporationTimelineError::MassBalance {
            residual_kg: audit.mass_residual_kg,
            tolerance_kg: phase_config.gas.mass_tolerance_kg,
        });
    }
    let birth_offset = finite(branch.born_mass_kg * surface.binding.liquid_enthalpy_offset_j_kg)?;
    add(
        &mut audit.thermal_energy_change_j,
        audit.thermal.stored_energy_change_j,
    )?;
    add(&mut audit.thermal_energy_change_j, birth_offset)?;
    audit.incoming_energy_j = finite(audit.thermal.incoming_energy_j + birth_offset)?;
    add(&mut audit.energy_residual_j, audit.thermal.residual_j)?;
    let tolerance = thermal_config
        .energy_tolerance_j
        .max(phase_config.gas.energy_tolerance_j);
    if audit.energy_residual_j.abs() > tolerance {
        return Err(GasEvaporationError::EnergyBalance {
            residual_j: audit.energy_residual_j,
            tolerance_j: tolerance,
        }
        .into());
    }
    checkpoint()?;
    Ok(GasEvaporationTimelineCandidate {
        thermal: next,
        gas: next_gas,
        catalog,
        audit,
    })
}

#[cfg(test)]
mod tests;
