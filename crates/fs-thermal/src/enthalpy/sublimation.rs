//! Sublimation on the existing conservative thermal timeline.
//!
//! Whole-cell operators require a pure solid on its original caloric branch
//! throughout the advance. Temperature alone cannot identify a mixed solid/liquid plateau;
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

/// Separate substrate/coating storage with heat transport, births and regrowth.
#[path = "sublimation/coating.rs"]
pub mod coating;

/// Multiple thermal cells sharing one finite vapor inventory.
#[path = "sublimation/chamber.rs"]
pub mod chamber;

/// Original thermal-owner admission for externally solved solid/gas exchange.
#[path = "sublimation/external.rs"]
pub mod external;

use super::observe::EnergySum;
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
    /// Evaluated from compensated per-cell relative enthalpies, independently of
    /// the mass-balance gate. Rounded aggregate energies can erase small transfers.
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
        if caloric.storage_basis() != super::ThermalStorageBasis::Enthalpy {
            return Err(DepositionError::InvalidInput(
                "sublimation requires enthalpy storage and an explicit open-mass model",
            )
            .into());
        }
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
        let mut candidate = cell.clone();
        candidate.mass_kg = step.solid.mass_kg;
        candidate.specific_h_j_kg = step.solid.specific_enthalpy_j_kg;
        caloric.validate_cell(&candidate)?;
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
        let mut phase_reference_h = self.cells.get(&id).map(|cell| cell.specific_h_j_kg);
        for (index, event) in self.schedule[self.cursor..].iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if let Event::Birth(cell) = &event.event {
                if cell.id == id && phase_reference_h.is_none() {
                    phase_reference_h = Some(cell.specific_h_j_kg);
                }
            }
        }
        let phase_reference_h = phase_reference_h.ok_or(DepositionError::InvalidInput(
            "sublimating cell is neither active nor scheduled",
        ))?;
        let initial_mass = super::finite(self.mass_kg())?;
        super::finite(self.energy_j())?;
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
        let mut phase_mass = EnergySum::default();
        let mut latent_energy = EnergySum::default();
        let mut vapor_energy = EnergySum::default();
        let mut relative_vapor = EnergySum::default();
        let mut recession = EnergySum::default();
        let mut external_heat = EnergySum::default();
        external_heat.add(audit.thermal.external_heat_j)?;
        // Indexed endpoints avoid an extra phase solve caused only by adding
        // decimal timesteps repeatedly. Native heat substeps keep their budget.
        let mut epoch_s = next.time_s;
        let mut interval = 1usize;
        while next.time_s < target_s {
            checkpoint()?;
            let start_s = next.time_s;
            let start_area = if next.cells.contains_key(&id) {
                Some(next.exposed_area_m2(id)?)
            } else {
                None
            };
            let empty = next.cells.is_empty();
            let event_time = next.schedule.get(next.cursor).map(|event| event.time_s);
            // No heat or phase work exists before the first active cell.
            // Jump to the next event/target instead of renewing empty intervals.
            let mut stop_s = if empty {
                target_s
            } else {
                super::finite(epoch_s + (interval as f64) * thermal_config.max_dt_s)?
                    .min(target_s)
            };
            if let Some(event_time) = event_time {
                stop_s = stop_s.min(event_time);
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
            external_heat.add(heat.external_heat_j)?;
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
                phase_mass.add(step.transferred_mass_kg)?;
                latent_energy.add(step.latent_energy_j)?;
                vapor_energy.add(step.vapor_energy_j)?;
                relative_vapor.add(super::finite(
                    step.transferred_mass_kg
                        * super::finite(step.solid.specific_enthalpy_j_kg - phase_reference_h)?,
                )?)?;
                relative_vapor.add(step.latent_energy_j)?;
                recession.add(step.recession_m)?;
                audit.phase_steps = audit
                    .phase_steps
                    .checked_add(1)
                    .ok_or(DepositionError::Arithmetic)?;
                audit.last_step = Some(step);
            }
            if empty || event_time == Some(stop_s) {
                epoch_s = stop_s;
                interval = 1;
            } else {
                interval = interval.checked_add(1).ok_or(DepositionError::Arithmetic)?;
            }
        }
        audit.transferred_mass_kg = phase_mass.value()?;
        audit.latent_energy_j = latent_energy.value()?;
        audit.vapor_energy_j = vapor_energy.value()?;
        audit.recession_m = recession.value()?;
        audit.thermal.external_heat_j = external_heat.value()?;

        // Cancel each cell's original or birth datum before reducing. Keeping
        // the mass balance separate avoids turning mass roundoff into a heat
        // source proportional to an arbitrary formation-enthalpy reference.
        let mut births = super::BTreeMap::new();
        let mut incoming_mass = EnergySum::default();
        let mut incoming_energy = EnergySum::default();
        for (index, event) in next.schedule[self.cursor..next.cursor].iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if let Event::Birth(cell) = &event.event {
                births.insert(cell.id, (cell.mass_kg, cell.specific_h_j_kg));
                incoming_mass.add(cell.mass_kg)?;
                incoming_energy.add(super::finite(cell.mass_kg * cell.specific_h_j_kg)?)?;
            }
        }
        audit.thermal.incoming_mass_kg = incoming_mass.value()?;
        audit.thermal.incoming_energy_j = incoming_energy.value()?;
        let mut mass_balance = EnergySum::default();
        let mut storage = EnergySum::default();
        let mut energy_balance = EnergySum::default();
        for (index, (&cell_id, cell)) in next.cells.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let original = self.cells.get(&cell_id);
            let (mass, h) = match original {
                Some(old) => (old.mass_kg, old.specific_h_j_kg),
                None => *births.get(&cell_id).ok_or(DepositionError::Arithmetic)?,
            };
            let relative_energy =
                super::finite(cell.mass_kg * super::finite(cell.specific_h_j_kg - h)?)?;
            energy_balance.add(relative_energy)?;
            storage.add(relative_energy)?;
            storage.add(if original.is_some() {
                super::finite(super::finite(cell.mass_kg - mass)? * h)?
            } else {
                super::finite(cell.mass_kg * h)?
            })?;
            mass_balance.add(super::finite(cell.mass_kg - mass)?)?;
        }
        mass_balance.add(audit.transferred_mass_kg)?;
        energy_balance.add(relative_vapor.value()?)?;
        energy_balance.add(-audit.thermal.external_heat_j)?;
        let final_mass = super::finite(next.mass_kg())?;
        audit.mass_balance_residual_kg = mass_balance.value()?;
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
        super::finite(next.energy_j())?;
        audit.stored_energy_change_j = storage.value()?;
        audit.energy_residual_j = energy_balance.value()?;
        if audit.energy_residual_j.abs() > thermal_config.energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(audit.energy_residual_j).into());
        }
        checkpoint()?;
        Ok((next, reservoir, audit))
    }
}

#[cfg(test)]
#[path = "sublimation/tests.rs"]
mod tests;
