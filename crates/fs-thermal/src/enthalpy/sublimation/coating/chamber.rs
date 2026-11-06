//! One finite vapor inventory for supported coatings on the native heat timeline.

use super::{
    CoatingAudit, CoatingSimulation, CoatingState, CoatingSurface, Ledger, accept_phase_step,
    admit_mass, supported_state,
};
use crate::enthalpy::sublimation::ThermalSublimationError;
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, StepConfig, finite};
use crate::sublimation::chamber::{ChamberSolverConfig, compensated_add};
use crate::sublimation::coating::chamber::{
    SupportedChamberStep, SupportedChamberSurface, advance_supported_chamber,
    advance_supported_chamber_with_heat,
};
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

    /// Apply signed heat and simultaneous phase exchange to all active coated
    /// substrates sharing this chamber. Sparse heat entries use persistent cell
    /// identity; omitted cells receive zero heat and still exchange vapor.
    /// Duplicate, inactive, unknown or nonfinite entries refuse.
    ///
    /// The returned identities are sorted and match the surface receipts exactly.
    /// Original substrate/deposit calorics, inventory, exposure and one endpoint
    /// gas pressure participate in the same solve. Authored heat sources remain
    /// caller-owned; this operation neither integrates them nor advances time.
    /// Every cell, coating and the shared reservoir commit together.
    pub fn exchange_chamber_with_heat(
        &mut self,
        duration_s: f64,
        supplied_heat_j: &[(CellId, f64)],
        config: ChamberSolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Vec<CellId>, SupportedChamberStep), ThermalSublimationError> {
        checkpoint()?;
        self.require_processed_events()?;
        let reservoir = self.chamber.ok_or(DepositionError::InvalidInput(
            "coatings have independent reservoirs",
        ))?;
        let old_vapor = vapor_mass(reservoir)?;
        if self.coatings.len() > config.max_surfaces {
            return Err(DepositionError::Budget("shared coating surface budget").into());
        }
        if supplied_heat_j.len() > self.coatings.len() {
            return Err(DepositionError::InvalidInput("too many chamber heat entries").into());
        }
        let mut heat_by_cell = BTreeMap::new();
        for &(id, heat) in supplied_heat_j {
            checkpoint()?;
            if !heat.is_finite()
                || !self.coatings.contains_key(&id)
                || heat_by_cell.insert(id, heat).is_some()
            {
                return Err(DepositionError::InvalidInput(
                    "chamber heat needs unique active coated cells and finite values",
                )
                .into());
            }
        }
        let areas = self.exposed_areas(&mut checkpoint)?;
        let mut ids = Vec::with_capacity(self.coatings.len());
        let mut surfaces = Vec::with_capacity(self.coatings.len());
        let mut heat = Vec::with_capacity(self.coatings.len());
        for (&id, coating) in &self.coatings {
            checkpoint()?;
            let model = &self.surfaces[&id];
            ids.push(id);
            heat.push(heat_by_cell.get(&id).copied().unwrap_or(0.0));
            surfaces.push(SupportedChamberSurface {
                state: supported_state(&self.thermal.cells[&id], *coating, areas[&id]),
                kinetics: model.kinetics,
                law: model.law,
                deposit_enthalpy: model.deposit_enthalpy,
            });
        }
        let mut refusal = None;
        let chart_refusal = std::cell::RefCell::new(None);
        let result = advance_supported_chamber_with_heat(
            &surfaces,
            reservoir,
            duration_s,
            &heat,
            config,
            |index, h| {
                let material = self.thermal.cells[&ids[index]].material_id;
                self.thermal.materials[&material]
                    .temperature_k(h)
                    .map_err(|error| {
                        *chart_refusal.borrow_mut() = Some(error);
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
        let mut step = match result {
            Ok(step) => step,
            Err(error) => {
                return Err(match chart_refusal.into_inner() {
                    Some(error) => error.into(),
                    None => error.into(),
                });
            }
        };
        let new_vapor = vapor_mass(step.reservoir)?;
        admit_mass(
            old_vapor,
            new_vapor,
            step.transferred_mass_kg,
            step.mass_residual_kg,
        )?;
        let mut prepared = Vec::with_capacity(ids.len());
        let mut deposit_energy = 0.0;
        let mut deposit_correction = 0.0;
        let mut residual = 0.0;
        let mut residual_correction = 0.0;
        for (&id, surface) in ids.iter().zip(&mut step.surfaces) {
            checkpoint()?;
            let cell = &self.thermal.cells[&id];
            let coating = self.coatings[&id];
            let actual_dm = finite(surface.state.deposit_mass_kg - coating.mass_kg)?;
            let mass_defect = finite(actual_dm + surface.phase.transferred_mass_kg)?;
            admit_mass(
                coating.mass_kg,
                surface.state.deposit_mass_kg,
                surface.phase.transferred_mass_kg,
                mass_defect,
            )?;
            let deposit_h = surface.phase.solid.specific_enthalpy_j_kg;
            let deposit_sensible_j =
                finite(coating.mass_kg * (deposit_h - coating.specific_enthalpy_j_kg))?;
            // Keep the original represented-inventory correction. The scalar
            // residual already subtracts supplied heat exactly once.
            surface.deposit_energy_change_j = finite(deposit_sensible_j + actual_dm * deposit_h)?;
            surface.phase.energy_residual_j =
                finite(surface.phase.energy_residual_j + finite(mass_defect * deposit_h)?)?;
            let physical_scale_j = surface
                .phase
                .supplied_heat_j
                .abs()
                .max(surface.substrate_energy_change_j.abs())
                .max(deposit_sensible_j.abs())
                .max(surface.phase.latent_energy_j.abs());
            let tolerance_j = finite(
                config.surface_solver.absolute_energy_tolerance_j
                    + config.surface_solver.relative_tolerance * physical_scale_j,
            )?;
            if surface.phase.energy_residual_j.abs() > tolerance_j {
                return Err(DepositionError::EnergyBalance(surface.phase.energy_residual_j).into());
            }
            compensated_add(
                &mut deposit_energy,
                &mut deposit_correction,
                surface.deposit_energy_change_j,
            );
            compensated_add(
                &mut residual,
                &mut residual_correction,
                surface.phase.energy_residual_j,
            );
            let mut candidate = cell.clone();
            candidate.specific_h_j_kg = surface.state.substrate_specific_enthalpy_j_kg;
            self.thermal.materials[&cell.material_id].validate_cell(&candidate)?;
            finite(candidate.mass_kg * candidate.specific_h_j_kg)?;
            finite(surface.state.deposit_mass_kg * deposit_h)?;
            prepared.push((
                id,
                candidate.specific_h_j_kg,
                CoatingState {
                    mass_kg: surface.state.deposit_mass_kg,
                    specific_enthalpy_j_kg: deposit_h,
                    reservoir: step.reservoir,
                },
            ));
        }
        step.deposit_energy_change_j = finite(deposit_energy + deposit_correction)?;
        step.energy_residual_j = finite(residual + residual_correction)?;
        checkpoint()?;
        // No fallible work or callbacks after this joint publication gate.
        for (id, substrate_h, coating) in prepared {
            self.thermal.cells.get_mut(&id).unwrap().specific_h_j_kg = substrate_h;
            *self.coatings.get_mut(&id).unwrap() = coating;
        }
        self.chamber = Some(step.reservoir);
        Ok((ids, step))
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
