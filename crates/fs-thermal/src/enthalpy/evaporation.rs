//! Continuum liquid evaporation and condensation on the original heat timeline.
//!
//! The original caloric law determines temperature, and the signed surface
//! transfer changes actual cell mass and enthalpy. The solver's h bracket is
//! the caller-admitted pure-liquid branch; births and EVERY heat substep are
//! checked against it. Temperature alone does not establish phase identity.
//! The prescribed gas film owns the complete carried vapor enthalpy. Latent
//! heat is its subaccount, never an additional source on the thermal timeline.
//!
//! Reference geometry, heat links and cell identities remain unchanged. This
//! is a cell-uniform, fixed-exposure model: recession is reported, not silently
//! turned into moving geometry. Exhaustion, freezing or boiling needs another
//! admitted phase/topology owner. Refine the full heat/phase split in time.

use super::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig, ThermalStorageBasis,
};
use crate::evaporation::{
    EvaporationError, EvaporationLaw, GasFilm, LiquidSurfaceState, SolverConfig, SurfaceStep,
    advance_surface,
};

/// A failed phase or thermal proposal leaves all accepted state unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalEvaporationError {
    Thermal(DepositionError),
    Surface(EvaporationError),
    /// Reference topology cannot retain a zero-mass thermal cell.
    CellExhausted(CellId),
    MassBalance {
        residual_kg: f64,
        tolerance_kg: f64,
    },
}
impl From<DepositionError> for ThermalEvaporationError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<EvaporationError> for ThermalEvaporationError {
    fn from(error: EvaporationError) -> Self {
        Self::Surface(error)
    }
}
impl core::fmt::Display for ThermalEvaporationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal evaporation refused: {self:?}")
    }
}
impl std::error::Error for ThermalEvaporationError {}

/// Whole-window accounts on the original caloric datum. These are numerical
/// balance checks, not estimates of time-splitting or physical model error.
#[derive(Clone, Debug, PartialEq)]
pub struct EvaporationAudit {
    pub thermal: EnergyAudit,
    /// Net liquid-to-vapor transfer [kg]; condensation is negative.
    pub transferred_mass_kg: f64,
    /// Signed latent-energy SUBACCOUNT of vapor_energy_j [J].
    pub latent_energy_j: f64,
    /// Complete carried enthalpy delivered to the prescribed film [J].
    pub vapor_energy_j: f64,
    /// Positive liquid recession [m]; reference geometry is retained.
    pub recession_m: f64,
    pub stored_energy_change_j: f64,
    pub mass_balance_residual_kg: f64,
    pub energy_residual_j: f64,
    pub phase_steps: usize,
    pub last_step: Option<SurfaceStep>,
}

/// An externally solved liquid endpoint bound to the original thermal cell.
/// The external owner must close the joint gas/liquid mass and energy account.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidCellExchange {
    /// Complete post-heat cell; stale mass, h, boundary or geometry refuses.
    pub expected: super::DepositedCell,
    /// Exposure during the preceding interval, before endpoint contact events.
    pub interval_exposed_area_m2: f64,
    pub replacement: LiquidSurfaceState,
    /// Caller-admitted pure-liquid interval on the original caloric datum.
    pub pure_liquid_enthalpy_range_j_kg: [f64; 2],
}

#[derive(Clone, Copy)]
struct LiquidBranch {
    id: CellId,
    range: [f64; 2],
}
impl LiquidBranch {
    fn new(id: CellId, range: [f64; 2]) -> Result<Self, DepositionError> {
        if !range.iter().all(|h| h.is_finite()) || range[0] > range[1] {
            return Err(DepositionError::InvalidInput(
                "finite ordered pure-liquid h bracket required",
            ));
        }
        Ok(Self { id, range })
    }
    fn check(&self, thermal: &DepositionSimulation<'_>) -> Result<(), DepositionError> {
        if let Some(cell) = thermal.cells.get(&self.id) {
            if thermal.materials[&cell.material_id].storage_basis() != ThermalStorageBasis::Enthalpy
            {
                return Err(DepositionError::InvalidInput(
                    "evaporation requires enthalpy storage and an explicit open-mass model",
                ));
            }
            if !(self.range[0]..=self.range[1]).contains(&cell.specific_h_j_kg) {
                return Err(DepositionError::MaterialDomain);
            }
        }
        Ok(())
    }
}
impl DepositionObserver for LiquidBranch {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check(thermal)
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &super::BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        self.check(thermal)
    }
}

fn sum(a: f64, b: f64) -> Result<f64, DepositionError> {
    super::finite(a + b)
}
fn add_heat(total: &mut EnergyAudit, part: EnergyAudit) -> Result<(), DepositionError> {
    total.incoming_mass_kg = sum(total.incoming_mass_kg, part.incoming_mass_kg)?;
    total.incoming_energy_j = sum(total.incoming_energy_j, part.incoming_energy_j)?;
    total.external_heat_j = sum(total.external_heat_j, part.external_heat_j)?;
    total.stored_energy_change_j = sum(total.stored_energy_change_j, part.stored_energy_change_j)?;
    total.residual_j = sum(total.residual_j, part.residual_j)?;
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

impl DepositionSimulation<'_> {
    /// Atomically admit a prepared batch from a finite-gas/interface owner.
    /// Identities, clocks, geometry, contacts and original calorics survive.
    /// Every endpoint and the final cancellation gate pass before any write.
    pub fn apply_liquid_exchanges(
        &mut self,
        exchanges: &[LiquidCellExchange],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        if exchanges.len() > self.budget.max_cells {
            return Err(DepositionError::Budget("external liquid exchange cells"));
        }
        let mut ids = std::collections::BTreeSet::new();
        for (index, exchange) in exchanges.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let id = exchange.expected.id;
            let cell = self.cells.get(&id).ok_or(DepositionError::InvalidInput(
                "external liquid exchange requires an active cell",
            ))?;
            if !ids.insert(id) || cell != &exchange.expected {
                return Err(DepositionError::InvalidInput(
                    "duplicate or stale external liquid exchange",
                ));
            }
            LiquidBranch::new(id, exchange.pure_liquid_enthalpy_range_j_kg)?.check(self)?;
            let [lo, hi] = exchange.pure_liquid_enthalpy_range_j_kg;
            let next = exchange.replacement;
            let area = exchange.interval_exposed_area_m2;
            if !(lo..=hi).contains(&next.specific_enthalpy_j_kg)
                || !super::positive(next.mass_kg)
                || !super::nonnegative(area)
                || area > cell.geometry.surface_area_m2
                || next.area_m2 != area
            {
                return Err(DepositionError::InvalidInput(
                    "external liquid endpoint leaves admitted branch, inventory or reference area",
                ));
            }
            let caloric = self.materials[&cell.material_id];
            for h in [cell.specific_h_j_kg, next.specific_enthalpy_j_kg] {
                if !super::positive(caloric.temperature_k(h)?) {
                    return Err(DepositionError::MaterialDomain);
                }
            }
            let mut candidate = cell.clone();
            candidate.mass_kg = next.mass_kg;
            candidate.specific_h_j_kg = next.specific_enthalpy_j_kg;
            caloric.validate_cell(&candidate)?;
        }
        checkpoint()?;
        for exchange in exchanges {
            let cell = self
                .cells
                .get_mut(&exchange.expected.id)
                .expect("all external liquid endpoints validated");
            cell.mass_kg = exchange.replacement.mass_kg;
            cell.specific_h_j_kg = exchange.replacement.specific_enthalpy_j_kg;
        }
        Ok(())
    }

    /// Advance signed evaporation on one actual active liquid cell without
    /// changing the clock. The caller owns admission of the pure-liquid h
    /// bracket and the prescribed, nonboiling continuum gas-film regime.
    /// Only exposed area exchanges mass. This does not alter convection or
    /// source heat, whose existing thermal owner remains independent.
    #[allow(clippy::too_many_arguments)]
    pub fn evaporate_cell<L: EvaporationLaw + ?Sized>(
        &mut self,
        id: CellId,
        law: &L,
        film: GasFilm,
        duration_s: f64,
        config: SolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<SurfaceStep, ThermalEvaporationError> {
        checkpoint()?;
        let area = self.exposed_area_m2(id)?;
        self.evaporate_cell_over_area(id, law, film, area, duration_s, config, &mut checkpoint)
    }

    #[allow(clippy::too_many_arguments)]
    fn evaporate_cell_over_area<L: EvaporationLaw + ?Sized>(
        &mut self,
        id: CellId,
        law: &L,
        film: GasFilm,
        area: f64,
        duration_s: f64,
        config: SolverConfig,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<SurfaceStep, ThermalEvaporationError> {
        LiquidBranch::new(id, config.specific_enthalpy_bracket_j_kg)?.check(self)?;
        let cell = self.cells.get(&id).ok_or(DepositionError::InvalidInput(
            "evaporation requires an active liquid cell",
        ))?;
        let caloric = self.materials[&cell.material_id];
        let chart_error = std::cell::RefCell::new(None);
        let mut checkpoint_error = None;
        let result = advance_surface(
            law,
            film,
            LiquidSurfaceState {
                mass_kg: cell.mass_kg,
                specific_enthalpy_j_kg: cell.specific_h_j_kg,
                area_m2: area,
            },
            duration_s,
            config,
            |h| {
                caloric.temperature_k(h).map_err(|error| {
                    *chart_error.borrow_mut() = Some(error);
                    EvaporationError::Invalid("original liquid caloric chart refused")
                })
            },
            || {
                checkpoint().map_err(|error| {
                    checkpoint_error = Some(error);
                    EvaporationError::Cancelled
                })
            },
        );
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let step = match result {
            Ok(step) => step,
            Err(error) => {
                if let Some(original) = chart_error.into_inner() {
                    return Err(original.into());
                }
                return Err(error.into());
            }
        };
        if !super::positive(step.liquid.mass_kg) {
            return Err(ThermalEvaporationError::CellExhausted(id));
        }
        checkpoint()?;
        // All fallible evaluation precedes the single physical-state commit.
        let cell = self
            .cells
            .get_mut(&id)
            .expect("validated active liquid cell");
        cell.mass_kg = step.liquid.mass_kg;
        cell.specific_h_j_kg = step.liquid.specific_enthalpy_j_kg;
        Ok(step)
    }

    /// Propose original heat transport and evaporation on one physical clock.
    ///
    /// Heat precedes the implicit phase operator over event-clipped intervals
    /// no longer than max_dt_s. The START exposed area applies throughout the
    /// interval: endpoint coverage never rewrites the preceding exposure, and
    /// cells cannot evaporate before birth. The original thermal substep budget
    /// is shared by the whole request, with no renewal between phase steps.
    ///
    /// This is a first-order split against one prescribed film. Its fixed
    /// humidity and temperature are boundary conditions, not a finite gas
    /// inventory. The returned vapor enthalpy includes latent energy once;
    /// a separate atmosphere owner must accept that complete account. A failed
    /// source, branch, budget, balance or cancellation discards the full advance.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_evaporating_to<L: EvaporationLaw + ?Sized>(
        &self,
        id: CellId,
        law: &L,
        film: GasFilm,
        target_s: f64,
        thermal_config: StepConfig,
        phase_config: SolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, EvaporationAudit), ThermalEvaporationError> {
        checkpoint()?;
        film.validate()?;
        if !target_s.is_finite()
            || target_s < self.time_s
            || !super::positive(thermal_config.max_dt_s)
            || !super::nonnegative(thermal_config.energy_tolerance_j)
        {
            return Err(DepositionError::InvalidInput(
                "invalid evaporation timeline target or accuracy budget",
            )
            .into());
        }
        let branch = LiquidBranch::new(id, phase_config.specific_enthalpy_bracket_j_kg)?;
        branch.check(self)?;
        let mut known = self.cells.contains_key(&id);
        for (index, event) in self.schedule[self.cursor..].iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            known |= matches!(&event.event, Event::Birth(cell) if cell.id == id);
        }
        if !known {
            return Err(DepositionError::InvalidInput(
                "evaporating cell is neither active nor scheduled",
            )
            .into());
        }
        let initial_mass = super::finite(self.mass_kg())?;
        let initial_energy = super::finite(self.energy_j())?;
        let (mut next, _, thermal) =
            self.advance_observed(&branch, self.time_s, thermal_config, &mut checkpoint)?;
        let mut audit = EvaporationAudit {
            thermal,
            transferred_mass_kg: 0.0,
            latent_energy_j: 0.0,
            vapor_energy_j: 0.0,
            recession_m: 0.0,
            stored_energy_change_j: 0.0,
            mass_balance_residual_kg: 0.0,
            energy_residual_j: 0.0,
            phase_steps: 0,
            last_step: None,
        };
        // Integer-indexed endpoints avoid a fictitious final phase interval
        // caused only by repeated addition of a decimal timestep.
        let mut epoch = next.time_s;
        let mut interval = 1usize;
        while next.time_s < target_s {
            checkpoint()?;
            let start = next.time_s;
            let area = if next.cells.contains_key(&id) {
                Some(next.exposed_area_m2(id)?)
            } else {
                None
            };
            let empty = next.cells.is_empty();
            let event_time = next.schedule.get(next.cursor).map(|event| event.time_s);
            let mut stop = if empty {
                target_s
            } else {
                super::finite(epoch + (interval as f64) * thermal_config.max_dt_s)?.min(target_s)
            };
            if let Some(event_time) = event_time {
                stop = stop.min(event_time);
            }
            if !stop.is_finite() || stop <= start {
                return Err(DepositionError::Arithmetic.into());
            }
            let remaining = next
                .budget
                .max_substeps
                .checked_sub(audit.thermal.substeps)
                .ok_or(DepositionError::Budget("evaporating thermal substeps"))?;
            let (heated, _, heat) = next.advance_observed_limited(
                &branch,
                stop,
                thermal_config,
                remaining,
                &mut checkpoint,
            )?;
            next = heated;
            add_heat(&mut audit.thermal, heat)?;
            if let Some(area) = area {
                let step = next.evaporate_cell_over_area(
                    id,
                    law,
                    film,
                    area,
                    stop - start,
                    phase_config,
                    &mut checkpoint,
                )?;
                audit.transferred_mass_kg =
                    sum(audit.transferred_mass_kg, step.transferred_mass_kg)?;
                audit.latent_energy_j = sum(audit.latent_energy_j, step.latent_energy_j)?;
                audit.vapor_energy_j = sum(audit.vapor_energy_j, step.vapor_energy_j)?;
                audit.recession_m = sum(audit.recession_m, step.recession_m)?;
                audit.phase_steps = audit
                    .phase_steps
                    .checked_add(1)
                    .ok_or(DepositionError::Arithmetic)?;
                audit.last_step = Some(step);
            }
            if empty || event_time == Some(stop) {
                epoch = stop;
                interval = 1;
            } else {
                interval = interval.checked_add(1).ok_or(DepositionError::Arithmetic)?;
            }
        }
        let final_mass = super::finite(next.mass_kg())?;
        audit.mass_balance_residual_kg = super::finite(
            final_mass - initial_mass - audit.thermal.incoming_mass_kg + audit.transferred_mass_kg,
        )?;
        let scale = initial_mass
            .abs()
            .max(final_mass.abs())
            .max(audit.thermal.incoming_mass_kg.abs())
            .max(audit.transferred_mass_kg.abs());
        let terms = self.cells.len() as f64
            + next.cells.len() as f64
            + audit.thermal.events as f64
            + audit.phase_steps as f64
            + 1.0;
        let tolerance = super::finite(128.0 * f64::EPSILON * scale * terms)?;
        if audit.mass_balance_residual_kg.abs() > tolerance {
            return Err(ThermalEvaporationError::MassBalance {
                residual_kg: audit.mass_balance_residual_kg,
                tolerance_kg: tolerance,
            });
        }
        audit.stored_energy_change_j = super::finite(next.energy_j() - initial_energy)?;
        audit.energy_residual_j = super::finite(
            audit.stored_energy_change_j
                - audit.thermal.incoming_energy_j
                - audit.thermal.external_heat_j
                + audit.vapor_energy_j,
        )?;
        if audit.energy_residual_j.abs() > thermal_config.energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(audit.energy_residual_j).into());
        }
        checkpoint()?;
        Ok((next, audit))
    }
}

#[cfg(test)]
#[path = "evaporation/tests.rs"]
mod tests;
