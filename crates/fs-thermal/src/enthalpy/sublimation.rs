//! Sublimation on the existing conservative thermal timeline.
//!
//! The cell must be a pure solid on its original caloric branch throughout the
//! advance. Temperature alone cannot identify a mixed solid/liquid plateau;
//! mixtures need an independently owned solid inventory and are not admitted by
//! this whole-cell model.
//!
//! The phase operator changes the actual cell mass and specific enthalpy. The
//! same accepted transfer supplies vapor energy and recession; latent heat is
//! never added again as a separate thermal source. Geometry remains the caller's
//! reference chart. Recession is reported, not silently converted into a mesh,
//! contact-area, conductivity, or mechanical remap. Exhaustion refuses until a
//! caller can supply that topology transition.
//!
//! `advance_sublimating_to` composes ordinary heat transport followed by the
//! implicit phase operator over bounded intervals (first-order Lie splitting).
//! Refine `StepConfig::max_dt_s` for coupling accuracy. An isothermal finite
//! vapor reservoir has a mass balance and an energy-transfer account, not a
//! solved gas temperature or flow field.

use super::{CellId, DepositionError, DepositionSimulation, EnergyAudit, Event, StepConfig};
use crate::sublimation::{
    SolidSurfaceState, SolverConfig, SublimationError, SublimationLaw, SurfaceKinetics,
    SurfaceStep, VaporReservoir, advance_surface,
};

/// Original thermal/constitutive failure, or an unsupported topology transition.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalSublimationError {
    /// Existing thermal timeline, caloric chart, budget, or cancellation refusal.
    Thermal(DepositionError),
    /// Shared surface law or nonlinear solve refused.
    Surface(SublimationError),
    /// A cell cannot disappear while its reference topology and links remain.
    CellExhausted(CellId),
    /// Signed solid-mass change did not match births minus vapor transfer.
    MassBalance { residual_kg: f64, tolerance_kg: f64 },
}
impl From<DepositionError> for ThermalSublimationError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl From<SublimationError> for ThermalSublimationError {
    fn from(value: SublimationError) -> Self {
        Self::Surface(value)
    }
}
impl core::fmt::Display for ThermalSublimationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal sublimation refused: {self:?}")
    }
}
impl std::error::Error for ThermalSublimationError {}

/// Accounts for one accepted heat-plus-phase advance. These check numerical
/// conservation; they are not estimates of splitting or spatial model error.
#[derive(Clone, Debug, PartialEq)]
pub struct SublimationAudit {
    /// Sum of the original heat-transport operator audits, including births.
    /// Its stored change excludes the separately reported phase operators.
    pub thermal: EnergyAudit,
    /// Net solid-to-vapor mass [kg]; negative means deposition.
    pub transferred_mass_kg: f64,
    /// Signed latent heat used by the phase operators [J].
    pub latent_energy_j: f64,
    /// Signed energy carried into the vapor [J], including solid sensible
    /// enthalpy on the SAME datum as the original thermal law. Already includes
    /// latent heat; do not add `latent_energy_j` again.
    pub vapor_energy_j: f64,
    /// Sum of normal recession increments [m]; negative means growth.
    /// Caller-owned reference geometry is unchanged.
    pub recession_m: f64,
    /// Actual change in the complete thermal timeline's stored energy [J].
    pub stored_energy_change_j: f64,
    /// Stored solid mass change minus births plus vapor transfer [kg].
    pub mass_balance_residual_kg: f64,
    /// Stored energy change minus births minus external heat plus vapor energy [J].
    pub energy_residual_j: f64,
    /// Number of phase operators, separate from ordinary conductive substeps.
    pub phase_steps: usize,
    /// Last phase result, or None when the cell was not yet born.
    pub last_step: Option<SurfaceStep>,
}

fn sum(a: f64, b: f64) -> Result<f64, DepositionError> {
    super::finite(a + b)
}
fn add_thermal(total: &mut EnergyAudit, part: EnergyAudit) -> Result<(), DepositionError> {
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
    /// Apply phase exchange to one actual active cell and its vapor reservoir.
    ///
    /// This phase operator does not advance the timeline's clock. Its duration
    /// is the same interval as a caller's separately integrated heat operator;
    /// use `advance_sublimating_to` for the complete bounded composition.
    /// The caller must supply a pure-solid cell and caloric branch throughout
    /// the solve; a melting plateau is not an all-solid inventory.
    /// The original material law converts enthalpy to temperature, and only the
    /// currently exposed area sublimates. A covered surface exchanges no mass.
    ///
    /// All evaluation, cancellation and extinction checks precede either write.
    /// A refusal leaves BOTH the thermal cell and reservoir unchanged. Stored
    /// reference geometry and contacts are retained; a changed mass must not be
    /// passed to a fixed-density/volume consumer without its explicit remap.
    #[allow(clippy::too_many_arguments)]
    pub fn sublimate_cell<L: SublimationLaw + ?Sized>(
        &mut self,
        id: CellId,
        law: &L,
        kinetics: SurfaceKinetics,
        reservoir: &mut VaporReservoir,
        duration_s: f64,
        config: SolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<SurfaceStep, ThermalSublimationError> {
        checkpoint()?;
        let area = self.exposed_area_m2(id)?;
        self.sublimate_cell_over_area(
            id,
            law,
            kinetics,
            reservoir,
            area,
            duration_s,
            config,
            &mut checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sublimate_cell_over_area<L: SublimationLaw + ?Sized>(
        &mut self,
        id: CellId,
        law: &L,
        kinetics: SurfaceKinetics,
        reservoir: &mut VaporReservoir,
        exposed_area_m2: f64,
        duration_s: f64,
        config: SolverConfig,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<SurfaceStep, ThermalSublimationError> {
        let cell = self.cells.get(&id).ok_or(DepositionError::InvalidInput(
            "sublimation requires an active cell",
        ))?;
        let caloric = self.materials[&cell.material_id];
        let solid = SolidSurfaceState {
            mass_kg: cell.mass_kg,
            specific_enthalpy_j_kg: cell.specific_h_j_kg,
            area_m2: exposed_area_m2,
        };
        // Preserve the original chart/checkpoint refusal in the adapter error.
        let chart_error = std::cell::RefCell::new(None);
        let mut checkpoint_error = None;
        let result = advance_surface(
            law,
            kinetics,
            solid,
            *reservoir,
            duration_s,
            config,
            |h| {
                caloric.temperature_k(h).map_err(|error| {
                    *chart_error.borrow_mut() = Some(error);
                    SublimationError::Invalid("thermal caloric chart refused phase endpoint")
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
        // A safeguarded solve may reject a trial outside the caloric domain
        // and still converge inside it. Only propagate a chart refusal when
        // the complete surface solve failed; successful brackets remain valid.
        let step = match result {
            Ok(step) => step,
            Err(error) => {
                if let Some(original) = chart_error.into_inner() {
                    return Err(original.into());
                }
                return Err(error.into());
            }
        };
        if step.solid.mass_kg <= 0.0 {
            return Err(ThermalSublimationError::CellExhausted(id));
        }
        checkpoint()?;
        // No fallible work follows: both owners commit the same accepted result.
        let cell = self
            .cells
            .get_mut(&id)
            .expect("validated active sublimating cell");
        cell.mass_kg = step.solid.mass_kg;
        cell.specific_h_j_kg = step.solid.specific_enthalpy_j_kg;
        *reservoir = step.reservoir;
        Ok(step)
    }

    /// Advance the original heat timeline and one surface's phase exchange as
    /// one transaction. The caller owns admission of a pure-solid caloric branch
    /// over the full heat/phase trajectory; mixed-phase cells need a phase mass
    /// partition that this whole-cell API does not infer.
    /// Conduction, authored heat, contact and scheduled events
    /// retain their existing owner. The phase solve uses the resulting enthalpy.
    ///
    /// Intervals end at every scheduled event and at most `max_dt_s` apart.
    /// The START exposed area applies throughout its preceding interval, so a
    /// contact/boundary change at the endpoint cannot retroactively cover it.
    /// A cell born at the endpoint receives no sublimation before its birth.
    /// The original total substep budget is shared across all inner advances.
    ///
    /// The supplied vapor reservoir belongs to this one surface. Reusing copied
    /// finite reservoirs for several surfaces does not model a shared chamber;
    /// that requires a joint mass/backpressure solve. The method solves neither
    /// gas temperature nor changing topology, and makes no moving-mesh claim.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sublimating_to<L: SublimationLaw + ?Sized>(
        &self,
        id: CellId,
        law: &L,
        kinetics: SurfaceKinetics,
        mut reservoir: VaporReservoir,
        target_s: f64,
        thermal_config: StepConfig,
        phase_config: SolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, VaporReservoir, SublimationAudit), ThermalSublimationError> {
        checkpoint()?;
        if !target_s.is_finite()
            || target_s < self.time_s
            || !super::positive(thermal_config.max_dt_s)
            || !super::nonnegative(thermal_config.energy_tolerance_j)
        {
            return Err(DepositionError::InvalidInput(
                "invalid sublimating timeline target, timestep or energy tolerance",
            )
            .into());
        }
        let mut known = self.cells.contains_key(&id);
        for (index, event) in self.schedule[self.cursor..].iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            known |= matches!(&event.event, Event::Birth(cell) if cell.id == id);
        }
        if !known {
            return Err(DepositionError::InvalidInput(
                "sublimating cell is neither active nor scheduled",
            )
            .into());
        }
        let initial_mass = super::finite(self.mass_kg())?;
        let initial_energy = super::finite(self.energy_j())?;
        // Activate exactly the ordinary same-time events before taking an area.
        let (mut next, initial_thermal) =
            self.advance_to(self.time_s, thermal_config, &mut checkpoint)?;
        let mut audit = SublimationAudit {
            thermal: initial_thermal,
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
        while next.time_s < target_s {
            checkpoint()?;
            let start_s = next.time_s;
            let start_area = if next.cells.contains_key(&id) {
                Some(next.exposed_area_m2(id)?)
            } else {
                None
            };
            let mut stop_s = (start_s + thermal_config.max_dt_s).min(target_s);
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
                .ok_or(DepositionError::Budget("sublimating thermal substeps"))?;
            let (heated, (), heat) = next.advance_observed_limited(
                &(),
                stop_s,
                thermal_config,
                remaining,
                &mut checkpoint,
            )?;
            next = heated;
            add_thermal(&mut audit.thermal, heat)?;
            if let Some(area) = start_area {
                let step = next.sublimate_cell_over_area(
                    id,
                    law,
                    kinetics,
                    &mut reservoir,
                    area,
                    stop_s - start_s,
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
        let terms = (self.cells.len() as f64)
            + (next.cells.len() as f64)
            + (audit.thermal.events as f64)
            + (audit.phase_steps as f64)
            + 1.0;
        let mass_tolerance = super::finite(128.0 * f64::EPSILON * scale * terms)?;
        if audit.mass_balance_residual_kg.abs() > mass_tolerance {
            return Err(ThermalSublimationError::MassBalance {
                residual_kg: audit.mass_balance_residual_kg,
                tolerance_kg: mass_tolerance,
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
        Ok((next, reservoir, audit))
    }
}

#[cfg(test)]
mod tests;
