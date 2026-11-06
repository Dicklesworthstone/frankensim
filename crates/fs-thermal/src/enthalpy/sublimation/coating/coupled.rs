//! Authored heat and supported phase exchange on the original native substeps.
//!
//! Both materials consume heat at their common phase endpoint. The existing
//! timeline still owns events, exposure, ordinary cells, explicit rates and CFL.

use super::{
    AdvanceMode, CoatingAudit, CoatingSimulation, Ledger, accept_phase_step, supported_state,
};
use crate::enthalpy::sublimation::ThermalSublimationError;
use crate::enthalpy::{CellId, DepositionError, StepConfig, finite};
use crate::sublimation::chamber::ChamberSolverConfig;
use crate::sublimation::coating::SupportedSurfaceStep;
use crate::sublimation::coating::advance_supported_surface_with_heat;
use crate::sublimation::coating::chamber::advance_supported_chamber;
use crate::sublimation::{SolverConfig, SublimationError};
use std::collections::BTreeMap;

/// Complete substrate/coating/vapor accounts, with no heat-only intermediate.
/// Latent energy is already part of vapor energy. Residuals are numerical
/// admission checks, not estimates of time or thin-coating model error.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledCoatingAudit {
    pub incoming_mass_kg: f64,
    pub incoming_energy_j: f64,
    /// Authored source and thermal-bath heat; internal contacts cancel.
    pub external_heat_j: f64,
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    pub vapor_energy_j: f64,
    pub stored_energy_change_j: f64,
    pub mass_balance_residual_kg: f64,
    pub energy_residual_j: f64,
    pub absolute_energy_residual_j: f64,
    pub substeps: usize,
    pub events: usize,
    pub phase_steps: usize,
    pub pressure_iterations: usize,
    pub absolute_nonlinear_mass_residual_kg: f64,
    /// Original scalar receipts. Global residuals also include the native
    /// owner's correction for actual represented inventory changes.
    pub last_steps: BTreeMap<CellId, SupportedSurfaceStep>,
}

impl From<CoatingAudit> for CoupledCoatingAudit {
    fn from(audit: CoatingAudit) -> Self {
        Self {
            incoming_mass_kg: audit.thermal.incoming_mass_kg,
            incoming_energy_j: audit.thermal.incoming_energy_j,
            external_heat_j: audit.thermal.external_heat_j,
            transferred_mass_kg: audit.transferred_mass_kg,
            latent_energy_j: audit.latent_energy_j,
            vapor_energy_j: audit.vapor_energy_j,
            stored_energy_change_j: audit.stored_energy_change_j,
            mass_balance_residual_kg: audit.mass_balance_residual_kg,
            energy_residual_j: audit.energy_residual_j,
            absolute_energy_residual_j: audit.absolute_energy_residual_j,
            substeps: audit.thermal.substeps,
            events: audit.thermal.events,
            phase_steps: audit.phase_steps,
            pressure_iterations: audit.pressure_iterations,
            absolute_nonlinear_mass_residual_kg: audit.absolute_nonlinear_mass_residual_kg,
            last_steps: audit.last_steps,
        }
    }
}

impl CoatingSimulation<'_> {
    /// Advance original heat transport and independent-bath coating exchange.
    ///
    /// Each native interval supplies its original source, bath and contact heat
    /// directly to the supported phase endpoint. Substrate and coating sensible
    /// storage participate together with latent transfer; an uncoated cell keeps
    /// its native update. No heat-only coating state is formed.
    ///
    /// Original event ordering, exposure, CFL and whole-call work limits remain
    /// authoritative. Heat quadrature is explicit and this composition is first
    /// order. Refine time independently of the nonlinear energy tolerances.
    /// Failure or cancellation publishes neither timeline nor vapor inventory.
    pub fn advance_coupled_to(
        &self,
        target_s: f64,
        thermal_config: StepConfig,
        phase_config: SolverConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, CoupledCoatingAudit), ThermalSublimationError> {
        if self.chamber.is_some() {
            return Err(DepositionError::InvalidInput(
                "advance a shared coating chamber with advance_coupled_chamber_to",
            )
            .into());
        }
        self.advance_mode(
            target_s,
            thermal_config,
            phase_config,
            AdvanceMode::Coupled(None),
            checkpoint,
        )
        .map(|(next, audit)| (next, audit.into()))
    }

    /// Advance native heat and all active coatings against one shared vapor.
    ///
    /// Each original cell receives its own signed heat once. One common
    /// endpoint-pressure closure couples the supported phase solves, including
    /// unheated bare receivers. A future substrate joins the then-current gas
    /// only at birth. The existing event loop and explicit heat CFL are shared
    /// with the independent-bath and split advances; no gas-temperature or
    /// moving-geometry owner is inferred.
    pub fn advance_coupled_chamber_to(
        &self,
        target_s: f64,
        thermal_config: StepConfig,
        config: ChamberSolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, CoupledCoatingAudit), ThermalSublimationError> {
        self.advance_coupled_chamber_with_checkpoint(
            target_s,
            thermal_config,
            config,
            &mut checkpoint,
        )
    }

    fn advance_coupled_chamber_with_checkpoint(
        &self,
        target_s: f64,
        thermal_config: StepConfig,
        config: ChamberSolverConfig,
        mut checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, CoupledCoatingAudit), ThermalSublimationError> {
        checkpoint()?;
        let reservoir = self.chamber.ok_or(DepositionError::InvalidInput(
            "coatings have independent reservoirs",
        ))?;
        if self.surfaces.len() > config.max_surfaces {
            return Err(DepositionError::Budget("shared coating surface budget").into());
        }
        // Admit the original chamber controls even before any substrate is born.
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
        self.advance_mode(
            target_s,
            thermal_config,
            config.surface_solver,
            AdvanceMode::Coupled(Some(config)),
            checkpoint,
        )
        .map(|(next, audit)| (next, audit.into()))
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn exchange_coupled(
        &mut self,
        areas: &BTreeMap<CellId, f64>,
        duration_s: f64,
        supplied_heat_j: &BTreeMap<CellId, f64>,
        phase_config: SolverConfig,
        chamber_config: Option<ChamberSolverConfig>,
        sums: &mut Ledger,
        audit: &mut CoatingAudit,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalSublimationError> {
        if let Some(config) = chamber_config {
            return self.exchange_chamber(
                areas,
                duration_s,
                Some(supplied_heat_j),
                config,
                sums,
                audit,
                checkpoint,
            );
        }
        for (&id, coating) in &mut self.coatings {
            checkpoint()?;
            let heat = *supplied_heat_j
                .get(&id)
                .ok_or(DepositionError::Arithmetic)?;
            let model = &self.surfaces[&id];
            let cell = self.thermal.cells.get_mut(&id).unwrap();
            let caloric = self.thermal.materials[&cell.material_id];
            let chart_error = std::cell::RefCell::new(None);
            let mut refusal = None;
            let result = advance_supported_surface_with_heat(
                model.law,
                model.kinetics,
                supported_state(cell, *coating, areas[&id]),
                coating.reservoir,
                duration_s,
                heat,
                phase_config,
                |h| {
                    caloric.temperature_k(h).map_err(|error| {
                        *chart_error.borrow_mut() = Some(error);
                        SublimationError::Invalid(
                            "original substrate calorics refused coupled coating trial",
                        )
                    })
                },
                model.deposit_enthalpy,
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
            let step = match result {
                Ok(step) => step,
                Err(error) => {
                    if let Some(original) = chart_error.into_inner() {
                        return Err(original.into());
                    }
                    return Err(error.into());
                }
            };
            let mut candidate = cell.clone();
            candidate.specific_h_j_kg = step.state.substrate_specific_enthalpy_j_kg;
            caloric.validate_cell(&candidate)?;
            finite(candidate.mass_kg * candidate.specific_h_j_kg)?;
            finite(step.state.deposit_mass_kg * step.phase.solid.specific_enthalpy_j_kg)?;
            accept_phase_step(cell, coating, step, true, sums, audit)?;
        }
        Ok(())
    }
}
