//! Original main-inventory gas heat on the complete load's solved trajectory.
use super::*;
use crate::bernoulli_aperture::wall::material::{
    FiniteGasDgInternalHeat, prepare_finite_gas_internal_heat,
};
use crate::dg_port::DgPortTrial;

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct SourceHeat {
    pub(super) main: FiniteGasDgInternalHeat,
    pub(super) terminal_viscous_heat_j: f64,
}

impl SourceHeat {
    pub(crate) fn gas_heat(self) -> FiniteGasDgInternalHeat {
        FiniteGasDgInternalHeat {
            molecular_heat_j: self.main.molecular_heat_j,
            viscous_heat_j: self.main.viscous_heat_j + self.terminal_viscous_heat_j,
        }
    }

    pub(crate) fn total_j<E>(self) -> Result<f64, FiniteGasDgError<E>> {
        self.main.total_j()?;
        if !self.terminal_viscous_heat_j.is_finite() || self.terminal_viscous_heat_j < 0.0 {
            return Err(FiniteGasDgError::Invalid(
                "original terminal viscosity heat must be finite and nonnegative",
            ));
        }
        // Use the same grouping as the transport owner's two source accounts.
        self.gas_heat().total_j()
    }
}

/// One original complete acoustic step, retained across caloric/transport trials.
/// Only this private proposal advances; observing its heat adds no second step.
#[derive(Clone, Debug)]
pub(crate) struct PreparedFiniteGasDgPortHeat {
    marched: DgPortLoad,
    output: DgStepOutput,
    source_heat: SourceHeat,
}

impl PreparedFiniteGasDgPortHeat {
    /// Main-bore and main-passage viscosity share the original gas inventory.
    /// The public step report retains their individual physical contributions.
    pub(crate) fn heat(&self) -> FiniteGasDgInternalHeat {
        self.source_heat.gas_heat()
    }
}

impl DgPortTrial<'_> {
    pub(crate) fn finite_gas_heat<E, C: FnMut() -> bool>(
        &self,
        controls: FiniteGasDgSolve,
        budgets: FiniteGasDgHeatBudget,
        cancelled: &mut C,
    ) -> Result<SourceHeat, FiniteGasDgError<E>> {
        check(cancelled)?;
        let duct = &self.owner.duct;
        let has_main_source = duct.has_viscothermal_memory()
            || duct.molecular_relaxation().is_some()
            || duct.molecular_relaxation_field().is_some();
        let main = if has_main_source {
            let inputs = self
                .main_inputs
                .as_deref()
                .ok_or(FiniteGasDgError::Invalid(
                    "main gas heat requires the original solved held-flow DG trajectory",
                ))?;
            prepare_finite_gas_internal_heat(
                duct,
                &self.owner.state,
                &self.owner.propagator,
                inputs,
                controls,
                budgets,
                cancelled,
            )?
        } else {
            // No physical source: no integral preparation or budget is required.
            FiniteGasDgInternalHeat::default()
        };
        let main_heat_j = main.total_j()?;
        let main_loss_j = self.main_dissipated_j();
        if !main_loss_j.is_finite()
            || (main_heat_j > 0.0 && main_loss_j - main_heat_j < -controls.energy_tolerance_j)
        {
            return Err(FiniteGasDgError::Invalid(
                "main gas heat exceeds its own original native decrement",
            ));
        }
        let mut terminal_viscous_heat_j = 0.0;
        let mut correction = 0.0;
        for frame in &self.terminal_channel_frames {
            check(cancelled)?;
            // Match the original finite-inventory volume admission: a branch
            // and its outlet passage belong to that branch's independent gas.
            if self
                .owner
                .branches
                .iter()
                .any(|branch| branch.port() == frame.port)
            {
                continue;
            }
            let Some(step) = frame.step else {
                continue;
            };
            let heat = step.viscous_dissipation_j;
            if !heat.is_finite() || heat < 0.0 {
                return Err(FiniteGasDgError::Invalid(
                    "native terminal viscosity heat must be finite and nonnegative",
                ));
            }
            // Hydraulic pressure work also pays kinetic storage and BE loss.
            // Only the original positive viscous term returns to this gas.
            let corrected = heat - correction;
            let next = terminal_viscous_heat_j + corrected;
            correction = (next - terminal_viscous_heat_j) - corrected;
            terminal_viscous_heat_j = next;
        }
        let heat = SourceHeat {
            main,
            terminal_viscous_heat_j,
        };
        heat.total_j()?;
        check(cancelled)?;
        Ok(heat)
    }
}

impl DgPortLoad {
    pub(crate) fn prepare_finite_gas_source_heat<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        inlet_flow_m3_s: f64,
        controls: FiniteGasDgSolve,
        budgets: FiniteGasDgHeatBudget,
        cancelled: &mut C,
    ) -> Result<PreparedFiniteGasDgPortHeat, FiniteGasDgError<E>> {
        admit_initial(self, gas, controls, cancelled)?;
        let mut marched = self.clone();
        let trial = marched.preview_step(inlet_flow_m3_s).map_err(acoustic)?;
        let source_heat = trial.finite_gas_heat(controls, budgets, cancelled)?;
        let output = trial.accept();
        check(cancelled)?;
        Ok(PreparedFiniteGasDgPortHeat {
            marched,
            output,
            source_heat,
        })
    }

    /// Complete only the gas/material feedback around the same prepared step.
    /// The callback includes prepared heat in its actual caloric proposal and
    /// excludes that internal transfer from its returned external heat account.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn complete_finite_gas_source_heat<S, E, C: FnMut() -> bool>(
        &self,
        prepared: &PreparedFiniteGasDgPortHeat,
        gas: &ClosedGasParcel,
        controls: FiniteGasDgSolve,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<FiniteGasDgPortCandidate<S>, FiniteGasDgError<E>> {
        prepare_after_step_with_heat(
            self,
            &prepared.marched,
            prepared.output.clone(),
            gas,
            controls,
            prepared.source_heat,
            evolve,
            cancelled,
        )
    }

    /// Return physical viscous and molecular heat to the original finite gas.
    ///
    /// One native trial solves the original inlet, terminal and branch-junction
    /// held flows. The main bore's existing source-positive observers integrate
    /// its COMPLETE original state under those exact inputs; separate named
    /// budgets cover their sequential preparations, with at most 128 native
    /// main-bore states and no truncation. Finite channels contribute their
    /// original resolved viscous heat only when they belong to this main gas.
    /// Molecular memory requires the original molecular-diffusion refresh mode.
    ///
    /// The same gas caloric inversion and reciprocal parameter-work solve pay
    /// these internal transfers once. External heat stays separate. Numerical,
    /// thermal-wall, solid-material and independent-branch losses remain exports;
    /// rational terminal exchange remains signed port work. No heat is inferred
    /// from channel pressure work or from total acoustic energy loss.
    ///
    /// All original wave, wall, branch, transverse-channel and exterior histories
    /// survive. Any refusal leaves the load and gas unchanged. The original
    /// fixed-volume, well-mixed and frozen-initial-medium partition applies;
    /// active main-bore sources on a continuous-port backend require a separate
    /// trajectory integral and refuse here. Refine the original space/time/loss
    /// controls; this is numerical accounting, not measured acoustic validation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_viscous_molecular_heat_controlled<C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        inlet_flow_m3_s: f64,
        external_heat_j: f64,
        gas_solve: GasEnergySolve,
        controls: FiniteGasDgSolve,
        budgets: FiniteGasDgHeatBudget,
        cancelled: &mut C,
    ) -> Result<FiniteGasDgPortStep<GasCompositionStep>, FiniteGasDgError<GasEnergyError>> {
        let prepared = self.prepare_finite_gas_source_heat(
            gas,
            inlet_flow_m3_s,
            controls,
            budgets,
            cancelled,
        )?;
        let internal_heat_j = prepared.heat().total_j()?;
        let total_heat = if internal_heat_j == 0.0 {
            external_heat_j
        } else {
            external_heat_j + internal_heat_j
        };
        if !external_heat_j.is_finite() || !total_heat.is_finite() {
            return Err(FiniteGasDgError::Invalid(
                "external and original gas-loss heat must have a finite sum",
            ));
        }
        let candidate = self.complete_finite_gas_source_heat(
            &prepared,
            gas,
            controls,
            |old, _, work, cancel| {
                let mut background = heated_finite_gas(old, total_heat, work, gas_solve, cancel)?;
                background.supplied_heat_j = external_heat_j;
                Ok(background)
            },
            cancelled,
        )?;
        *self = candidate.load;
        *gas = candidate.gas;
        Ok(candidate.step)
    }
}
