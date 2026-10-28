//! Several real solid cells exchanging one finite vapor inventory.
//!
//! All cells use the same volatile law and compatible solid enthalpy datum.
//! Each original caloric chart must remain on a pure-solid branch. One joint
//! chamber solve owns shared backpressure; no surface gets an independent
//! copy of the initial vapor inventory. Gas temperature remains prescribed.
//!
//! The complete advance uses first-order heat/phase splitting and the original
//! event, conduction and work-budget owners. Surface area is frozen at each
//! interval's start. Reference geometry is retained; per-cell recession does
//! not silently remap contact, flow or mechanical geometry.

use std::collections::{BTreeMap, BTreeSet};

use super::{ThermalSublimationError, add_thermal, sum};
use crate::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, Event, StepConfig,
};
use crate::sublimation::{
    SolidSurfaceState, SublimationError, SublimationLaw, SurfaceKinetics, VaporReservoir,
    chamber::{self as phase, ChamberSolverConfig, ChamberStep},
};

/// One whole-solid cell and its local surface kinetics. Selection is sorted by
/// persistent cell identity before solving, independently of caller order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChamberSurface {
    pub cell_id: CellId,
    pub kinetics: SurfaceKinetics,
}

/// Joint result and the corresponding stable cell identities.
#[derive(Clone, Debug, PartialEq)]
pub struct ChamberCellStep {
    /// Canonical cell order matching `exchange.surfaces` exactly.
    pub cell_ids: Vec<CellId>,
    /// One shared reservoir; every surface result repeats that SAME endpoint.
    pub exchange: ChamberStep,
}

/// Cumulative local phase exchange. Recession belongs to one surface; different
/// surfaces' recession lengths are never added into a meaningless global value.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CellPhaseAudit {
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    pub vapor_energy_j: f64,
    pub recession_m: f64,
    pub phase_steps: usize,
}

/// Conservation and actual work on the common accepted timeline.
#[derive(Clone, Debug, PartialEq)]
pub struct ChamberAudit {
    /// Original heat-operator audits, excluding the separate phase operators.
    pub thermal: EnergyAudit,
    /// Per-cell phase accounts, including zero accounts for future births.
    pub cells: BTreeMap<CellId, CellPhaseAudit>,
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    /// Complete sensible-plus-latent energy into the thermostatted vapor bath.
    /// Already includes latent energy; never add the latter a second time.
    pub vapor_energy_j: f64,
    pub stored_energy_change_j: f64,
    /// Actual change in solid PLUS vapor mass, minus incoming material births.
    pub mass_balance_residual_kg: f64,
    /// Sum of accepted joint-solver defects: vapor mass change minus transfer.
    /// The actual whole-window defect must agree with this to roundoff.
    pub phase_mass_residual_kg: f64,
    /// Solid energy change minus births minus external heat plus vapor energy.
    pub energy_residual_j: f64,
    pub phase_steps: usize,
    pub pressure_iterations: usize,
    pub last_step: Option<ChamberCellStep>,
}

fn finite(x: f64) -> Result<f64, DepositionError> {
    crate::enthalpy::finite(x)
}
fn finite_vapor_mass(reservoir: VaporReservoir) -> Result<f64, ThermalSublimationError> {
    match reservoir {
        VaporReservoir::Finite {
            vapor_mass_kg,
            volume_m3,
            temperature_k,
        } if crate::enthalpy::nonnegative(vapor_mass_kg)
            && crate::enthalpy::positive(volume_m3)
            && crate::enthalpy::positive(temperature_k) =>
        {
            Ok(vapor_mass_kg)
        }
        _ => Err(SublimationError::Invalid(
            "shared chamber needs finite nonnegative vapor mass and positive volume/temperature",
        )
        .into()),
    }
}

impl DepositionSimulation<'_> {
    fn chamber_selection(
        &self,
        surfaces: &[ChamberSurface],
        max_surfaces: usize,
        allow_scheduled: bool,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<ChamberSurface>, ThermalSublimationError> {
        checkpoint()?;
        if surfaces.is_empty() {
            return Err(
                DepositionError::InvalidInput("shared chamber needs selected solid cells").into(),
            );
        }
        if surfaces.len() > self.budget.max_cells || surfaces.len() > max_surfaces {
            return Err(DepositionError::Budget("shared chamber surface selection").into());
        }
        let mut known = BTreeSet::new();
        for (index, id) in self.cells.keys().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            known.insert(*id);
        }
        if allow_scheduled {
            for (index, event) in self.schedule[self.cursor..].iter().enumerate() {
                if index % 256 == 0 {
                    checkpoint()?;
                }
                if let Event::Birth(cell) = &event.event {
                    known.insert(cell.id);
                }
            }
        }
        let mut ordered = surfaces.to_vec();
        ordered.sort_unstable_by_key(|surface| surface.cell_id);
        for (index, surface) in ordered.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if !known.contains(&surface.cell_id)
                || (index > 0 && ordered[index - 1].cell_id == surface.cell_id)
            {
                return Err(DepositionError::InvalidInput(
                    "duplicate, inactive or unscheduled chamber cell",
                )
                .into());
            }
        }
        Ok(ordered)
    }

    // Capture START exposure once for all active selected surfaces. Future
    // births stay absent, and endpoint contact events cannot rewrite this batch.
    fn chamber_exposure(
        &self,
        ordered: &[ChamberSurface],
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<(ChamberSurface, f64)>, ThermalSublimationError> {
        let mut selected = Vec::with_capacity(ordered.len());
        for (index, surface) in ordered.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if self.cells.contains_key(&surface.cell_id) { selected.push(surface.cell_id); }
        }
        let areas = self.exposed_areas_m2(&selected, &mut *checkpoint)?;
        let mut exposed = Vec::with_capacity(selected.len());
        for (index, surface) in ordered.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if let Some(&area) = areas.get(&surface.cell_id) { exposed.push((*surface, area)); }
        }
        checkpoint()?;
        Ok(exposed)
    }

    /// Apply one simultaneous phase interval to actual active cells and ONE
    /// finite reservoir. Original caloric laws and current exposed areas are
    /// retained. Refusal leaves every cell and the reservoir unchanged.
    ///
    /// This phase operation does not advance the clock; use the complete
    /// `advance_chamber_sublimating_to` method for heat/phase time composition.
    /// Whole-cell solid inventory, a common enthalpy datum, fixed reference
    /// topology and a thermostatted gas remain explicit model requirements.
    #[allow(clippy::too_many_arguments)]
    pub fn sublimate_chamber<L: SublimationLaw + ?Sized>(
        &mut self,
        surfaces: &[ChamberSurface],
        law: &L,
        reservoir: &mut VaporReservoir,
        duration_s: f64,
        config: ChamberSolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ChamberCellStep, ThermalSublimationError> {
        finite_vapor_mass(*reservoir)?;
        let ordered =
            self.chamber_selection(surfaces, config.max_surfaces, false, &mut checkpoint)?;
        let exposed = self.chamber_exposure(&ordered, &mut checkpoint)?;
        self.sublimate_chamber_over_areas(
            &exposed,
            law,
            reservoir,
            duration_s,
            config,
            &mut checkpoint,
        )
    }

    fn sublimate_chamber_over_areas<L: SublimationLaw + ?Sized>(
        &mut self,
        exposed: &[(ChamberSurface, f64)],
        law: &L,
        reservoir: &mut VaporReservoir,
        duration_s: f64,
        config: ChamberSolverConfig,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ChamberCellStep, ThermalSublimationError> {
        let mut cell_ids = Vec::with_capacity(exposed.len());
        let mut inputs = Vec::with_capacity(exposed.len());
        let mut caloric = Vec::with_capacity(exposed.len());
        for (index, (surface, area)) in exposed.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let cell = self
                .cells
                .get(&surface.cell_id)
                .ok_or(DepositionError::InvalidInput(
                    "shared chamber phase step needs active solid cells",
                ))?;
            cell_ids.push(surface.cell_id);
            caloric.push(self.materials[&cell.material_id]);
            inputs.push(phase::ChamberSurface {
                solid: SolidSurfaceState {
                    mass_kg: cell.mass_kg,
                    specific_enthalpy_j_kg: cell.specific_h_j_kg,
                    area_m2: *area,
                },
                kinetics: surface.kinetics,
            });
        }
        let chart_error = std::cell::RefCell::new(None);
        let mut checkpoint_error = None;
        let result = phase::advance_chamber(
            law,
            &inputs,
            *reservoir,
            duration_s,
            config,
            |index, h| {
                caloric[index].temperature_k(h).map_err(|error| {
                    *chart_error.borrow_mut() = Some(error);
                    SublimationError::Invalid("thermal caloric chart refused chamber endpoint")
                })
            },
            || {
                checkpoint().map_err(|error| {
                    checkpoint_error = Some(error);
                    SublimationError::Cancelled
                })
            },
        );
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let exchange = match result {
            Ok(exchange) => exchange,
            Err(error) => {
                return Err(match chart_error.into_inner() {
                    Some(error) => error.into(),
                    None => error.into(),
                });
            }
        };
        if exchange.surfaces.len() != cell_ids.len() {
            return Err(DepositionError::InvalidInput(
                "chamber result does not match selected cell count",
            )
            .into());
        }
        for (index, (&id, step)) in cell_ids.iter().zip(&exchange.surfaces).enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if step.solid.mass_kg <= 0.0 {
                return Err(ThermalSublimationError::CellExhausted(id));
            }
            finite(step.solid.mass_kg)?;
            finite(step.solid.specific_enthalpy_j_kg)?;
        }
        checkpoint()?;
        // The bounded commit has no fallible work or callbacks. All endpoints
        // have passed; never publish an earlier surface before a later refusal.
        for (&id, step) in cell_ids.iter().zip(&exchange.surfaces) {
            let cell = self.cells.get_mut(&id).expect("validated chamber cell");
            cell.mass_kg = step.solid.mass_kg;
            cell.specific_h_j_kg = step.solid.specific_enthalpy_j_kg;
        }
        *reservoir = exchange.reservoir;
        Ok(ChamberCellStep { cell_ids, exchange })
    }

    /// Advance the existing heat graph and multiple selected surfaces sharing
    /// one finite vapor chamber as a single immutable candidate.
    ///
    /// Heat is followed by one JOINT phase solve at each bounded interval.
    /// Tighten max_dt_s for splitting accuracy. Intervals end exactly at events;
    /// start exposure is retained even if endpoint contacts change. A cell born
    /// at an endpoint first participates in the following interval. Per-cell
    /// recession is reported without remapping fixed reference geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_chamber_sublimating_to<L: SublimationLaw + ?Sized>(
        &self,
        surfaces: &[ChamberSurface],
        law: &L,
        mut reservoir: VaporReservoir,
        target_s: f64,
        thermal_config: StepConfig,
        phase_config: ChamberSolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, VaporReservoir, ChamberAudit), ThermalSublimationError> {
        checkpoint()?;
        if !target_s.is_finite()
            || target_s < self.time_s
            || !crate::enthalpy::positive(thermal_config.max_dt_s)
            || !crate::enthalpy::nonnegative(thermal_config.energy_tolerance_j)
        {
            return Err(DepositionError::InvalidInput(
                "invalid shared chamber timeline target, timestep or energy tolerance",
            )
            .into());
        }
        let initial_vapor_mass = finite_vapor_mass(reservoir)?;
        let ordered =
            self.chamber_selection(surfaces, phase_config.max_surfaces, true, &mut checkpoint)?;
        let initial_mass = finite(self.mass_kg())?;
        let initial_energy = finite(self.energy_j())?;
        let (mut next, initial_thermal) =
            self.advance_to(self.time_s, thermal_config, &mut checkpoint)?;
        let mut audit = ChamberAudit {
            thermal: initial_thermal,
            cells: ordered
                .iter()
                .map(|surface| (surface.cell_id, CellPhaseAudit::default()))
                .collect(),
            transferred_mass_kg: 0.0,
            latent_energy_j: 0.0,
            vapor_energy_j: 0.0,
            stored_energy_change_j: 0.0,
            mass_balance_residual_kg: 0.0,
            phase_mass_residual_kg: 0.0,
            energy_residual_j: 0.0,
            phase_steps: 0,
            pressure_iterations: 0,
            last_step: None,
        };
        while next.time_s < target_s {
            checkpoint()?;
            let start_s = next.time_s;
            let exposed = next.chamber_exposure(&ordered, &mut checkpoint)?;
            // No heat or phase work exists before the first active cell.
            // Jump to the next event/target instead of renewing empty intervals.
            let mut stop_s = if next.cells.is_empty() {
                target_s
            } else {
                (start_s + thermal_config.max_dt_s).min(target_s)
            };
            if let Some(event) = next.schedule.get(next.cursor) {
                stop_s = stop_s.min(event.time_s);
            }
            if !stop_s.is_finite() || stop_s <= start_s {
                return Err(DepositionError::Arithmetic.into());
            }
            let remaining = next
                .budget
                .max_substeps
                .checked_sub(audit.thermal.substeps)
                .ok_or(DepositionError::Budget("shared chamber thermal substeps"))?;
            let (heated, (), heat) = next.advance_observed_limited(
                &(),
                stop_s,
                thermal_config,
                remaining,
                &mut checkpoint,
            )?;
            next = heated;
            add_thermal(&mut audit.thermal, heat)?;
            if exposed.is_empty() {
                continue;
            }
            let step = next.sublimate_chamber_over_areas(
                &exposed,
                law,
                &mut reservoir,
                stop_s - start_s,
                phase_config,
                &mut checkpoint,
            )?;
            let exchange = &step.exchange;
            audit.transferred_mass_kg =
                sum(audit.transferred_mass_kg, exchange.transferred_mass_kg)?;
            audit.latent_energy_j = sum(audit.latent_energy_j, exchange.latent_energy_j)?;
            audit.vapor_energy_j = sum(audit.vapor_energy_j, exchange.vapor_energy_j)?;
            audit.phase_mass_residual_kg =
                sum(audit.phase_mass_residual_kg, exchange.mass_residual_kg)?;
            audit.phase_steps = audit
                .phase_steps
                .checked_add(1)
                .ok_or(DepositionError::Arithmetic)?;
            audit.pressure_iterations = audit
                .pressure_iterations
                .checked_add(exchange.pressure_iterations)
                .ok_or(DepositionError::Arithmetic)?;
            for (&id, surface) in step.cell_ids.iter().zip(&exchange.surfaces) {
                let account = audit.cells.get_mut(&id).expect("selected chamber account");
                account.transferred_mass_kg =
                    sum(account.transferred_mass_kg, surface.transferred_mass_kg)?;
                account.latent_energy_j = sum(account.latent_energy_j, surface.latent_energy_j)?;
                account.vapor_energy_j = sum(account.vapor_energy_j, surface.vapor_energy_j)?;
                account.recession_m = sum(account.recession_m, surface.recession_m)?;
                account.phase_steps = account
                    .phase_steps
                    .checked_add(1)
                    .ok_or(DepositionError::Arithmetic)?;
            }
            audit.last_step = Some(step);
        }
        let final_mass = finite(next.mass_kg())?;
        let final_vapor_mass = finite_vapor_mass(reservoir)?;
        let solid_residual = finite(
            final_mass - initial_mass - audit.thermal.incoming_mass_kg + audit.transferred_mass_kg,
        )?;
        audit.mass_balance_residual_kg = finite(
            final_mass - initial_mass - audit.thermal.incoming_mass_kg + final_vapor_mass
                - initial_vapor_mass,
        )?;
        let scale = initial_mass
            .abs()
            .max(final_mass.abs())
            .max(initial_vapor_mass.abs())
            .max(final_vapor_mass.abs())
            .max(audit.thermal.incoming_mass_kg.abs())
            .max(audit.transferred_mass_kg.abs());
        let terms = self.cells.len() as f64
            + next.cells.len() as f64
            + audit.thermal.events as f64
            + audit.phase_steps as f64
            + 1.0;
        let roundoff = finite(128.0 * f64::EPSILON * scale * terms)?;
        let unexplained_mass =
            finite(audit.mass_balance_residual_kg - audit.phase_mass_residual_kg)?;
        if solid_residual.abs() > roundoff {
            return Err(ThermalSublimationError::MassBalance {
                residual_kg: solid_residual,
                tolerance_kg: roundoff,
            });
        }
        if unexplained_mass.abs() > roundoff {
            return Err(ThermalSublimationError::MassBalance {
                residual_kg: unexplained_mass,
                tolerance_kg: roundoff,
            });
        }
        audit.stored_energy_change_j = finite(next.energy_j() - initial_energy)?;
        audit.energy_residual_j = finite(
            audit.stored_energy_change_j
                - audit.thermal.incoming_energy_j
                - audit.thermal.external_heat_j
                + audit.vapor_energy_j,
        )?;
        if audit.energy_residual_j.abs() > thermal_config.energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(audit.energy_residual_j).into());
        }
        checkpoint()?;
        Ok((next, reservoir, audit))
    }
}

#[cfg(test)]
#[path = "chamber/tests.rs"]
mod tests;
