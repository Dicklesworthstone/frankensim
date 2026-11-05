//! One finite vapor inventory for supported coatings on the native heat timeline.

use super::{
    CoatingAudit, CoatingSimulation, CoatingSurface, Ledger, accept_phase_step, admit_mass,
    supported_state,
};
use crate::enthalpy::sublimation::ThermalSublimationError;
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, StepConfig, finite};
use crate::sublimation::chamber::ChamberSolverConfig;
use crate::sublimation::coating::chamber::{SupportedChamberSurface, advance_supported_chamber};
use crate::sublimation::{SublimationError, SublimationLaw, SurfaceKinetics, VaporReservoir};
use std::collections::BTreeMap;

/// A persistent substrate admitted to heterogeneous deposition from one shared
/// volatile species. Zero inventory means an initially bare substrate. Every
/// surface receives the same original volatile and deposit enthalpy laws.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChamberCoatingSurface {
    pub cell_id: CellId,
    pub kinetics: SurfaceKinetics,
    pub initial_mass_kg: f64,
}

fn vapor_mass(reservoir: VaporReservoir) -> Result<f64, ThermalSublimationError> {
    match reservoir {
        VaporReservoir::Finite { vapor_mass_kg, volume_m3, temperature_k }
            if super::nonnegative(vapor_mass_kg)
                && super::positive(volume_m3)
                && super::positive(temperature_k) => Ok(vapor_mass_kg),
        _ => Err(SublimationError::Invalid(
            "shared coating chamber requires finite nonnegative vapor and positive volume/temperature",
        ).into()),
    }
}

impl<'a> CoatingSimulation<'a> {
    /// Connect active and scheduled substrates to ONE finite thermostatted gas.
    ///
    /// Substrate caloric datums remain independent. The deposit law and h(T)
    /// are shared, fixing the volatile species and its transported energy datum.
    /// Initial coating on a future substrate enters with that substrate's birth;
    /// its vapor supply is the chamber's then-current inventory. No nucleation,
    /// porosity, changing geometry or gas sensible-energy equation is inferred.
    pub fn new_shared_chamber(
        thermal: DepositionSimulation<'a>,
        surfaces: &[ChamberCoatingSurface],
        law: &'a dyn SublimationLaw,
        deposit_enthalpy: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
        reservoir: VaporReservoir,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalSublimationError> {
        checkpoint()?;
        vapor_mass(reservoir)?;
        if surfaces.len() > thermal.budget.max_cells {
            return Err(DepositionError::Budget("shared coating surfaces").into());
        }
        let mut selected = Vec::with_capacity(surfaces.len());
        for surface in surfaces {
            checkpoint()?;
            selected.push(CoatingSurface {
                cell_id: surface.cell_id,
                law,
                deposit_enthalpy,
                kinetics: surface.kinetics,
                initial_mass_kg: surface.initial_mass_kg,
                reservoir,
            });
        }
        let mut state = Self::new(thermal, &selected, &mut checkpoint)?;
        state.chamber = Some(reservoir);
        checkpoint()?;
        Ok(state)
    }

    /// The single shared inventory, including intervals with no active substrate.
    /// Per-coating reservoir observations repeat this state and must not be summed.
    pub fn shared_reservoir(&self) -> Option<VaporReservoir> {
        self.chamber
    }

    /// Replace the common external gas at a fully processed endpoint. This is an
    /// externally supplied inventory/thermostat change, not a phase-transfer step.
    pub fn set_shared_reservoir(
        &mut self,
        reservoir: VaporReservoir,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalSublimationError> {
        checkpoint()?;
        self.require_processed_events()?;
        if self.chamber.is_none() {
            return Err(
                DepositionError::InvalidInput("coatings have independent reservoirs").into(),
            );
        }
        vapor_mass(reservoir)?;
        // Complete all fallible work before publishing any of the observations.
        for _ in self.coatings.values() {
            checkpoint()?;
        }
        checkpoint()?;
        self.chamber = Some(reservoir);
        for coating in self.coatings.values_mut() {
            coating.reservoir = reservoir;
        }
        Ok(())
    }

    /// Advance original heat transport and simultaneous shared-vapor exchange.
    ///
    /// Heat then phase is first-order splitting on each native CFL substep.
    /// Explicit surface/pressure work budgets apply to every substep. The entire
    /// candidate, including late births and gas inventory, is atomic on refusal.
    pub fn advance_chamber_to(
        &self,
        target_s: f64,
        thermal_config: StepConfig,
        config: ChamberSolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, CoatingAudit), ThermalSublimationError> {
        checkpoint()?;
        let reservoir = self.chamber.ok_or(DepositionError::InvalidInput(
            "coatings have independent reservoirs",
        ))?;
        if self.surfaces.len() > config.max_surfaces {
            return Err(DepositionError::Budget("shared coating surface budget").into());
        }
        // Use the common numerical admission even before the first substrate is
        // born or when target==time. An empty interval cannot bypass bad controls.
        let mut refusal = None;
        let admission = advance_supported_chamber(
            &[],
            reservoir,
            0.0,
            config,
            |_, _| Err(SublimationError::Invalid("empty chamber has no substrate")),
            || {
                checkpoint().map_err(|error| {
                    refusal = Some(error);
                    SublimationError::Cancelled
                })
            },
        );
        if let Some(error) = refusal {
            return Err(error.into());
        }
        admission?;
        self.advance_impl(
            target_s,
            thermal_config,
            config.surface_solver,
            Some(config),
            checkpoint,
        )
    }

    pub(super) fn exchange_chamber(
        &mut self,
        areas: &BTreeMap<CellId, f64>,
        duration_s: f64,
        config: ChamberSolverConfig,
        sums: &mut Ledger,
        audit: &mut CoatingAudit,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalSublimationError> {
        let reservoir = self.chamber.ok_or(DepositionError::InvalidInput(
            "missing shared coating inventory",
        ))?;
        // BTreeMap cell identity fixes evaluation and summation order. Only
        // active cells participate; a future birth cannot consume vapor early.
        let mut ids = Vec::with_capacity(self.coatings.len());
        let mut surfaces = Vec::with_capacity(self.coatings.len());
        for (&id, coating) in &self.coatings {
            checkpoint()?;
            let model = &self.surfaces[&id];
            ids.push(id);
            surfaces.push(SupportedChamberSurface {
                state: supported_state(&self.thermal.cells[&id], *coating, areas[&id]),
                kinetics: model.kinetics,
                law: model.law,
                deposit_enthalpy: model.deposit_enthalpy,
            });
        }
        let mut refusal = None;
        let result = advance_supported_chamber(
            &surfaces,
            reservoir,
            duration_s,
            config,
            |index, h| {
                let material = self.thermal.cells[&ids[index]].material_id;
                self.thermal.materials[&material]
                    .temperature_k(h)
                    .map_err(|_| {
                        SublimationError::Invalid(
                            "original substrate calorics refused chamber trial",
                        )
                    })
            },
            || {
                checkpoint().map_err(|error| {
                    refusal = Some(error);
                    SublimationError::Cancelled
                })
            },
        );
        if let Some(error) = refusal {
            return Err(error.into());
        }
        let step = result?;
        let old = vapor_mass(reservoir)?;
        let new = vapor_mass(step.reservoir)?;
        let vapor_defect = finite((new - old) - step.transferred_mass_kg)?;
        admit_mass(old, new, step.transferred_mass_kg, vapor_defect)?;
        sums.mass_residual.add(vapor_defect)?;
        audit.pressure_iterations = audit
            .pressure_iterations
            .checked_add(step.pressure_iterations)
            .ok_or(DepositionError::Arithmetic)?;
        audit.absolute_nonlinear_mass_residual_kg = finite(
            audit.absolute_nonlinear_mass_residual_kg + step.nonlinear_mass_residual_kg.abs(),
        )?;
        for (id, surface_step) in ids.into_iter().zip(step.surfaces) {
            checkpoint()?;
            accept_phase_step(
                self.thermal.cells.get_mut(&id).unwrap(),
                self.coatings.get_mut(&id).unwrap(),
                surface_step,
                false,
                sums,
                audit,
            )?;
        }
        self.chamber = Some(step.reservoir);
        Ok(())
    }
}
