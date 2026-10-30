//! Differential-pressure closure around the existing solid/material heat solve.
//!
//! A pressure trial does not freeze walls from another mass-flow state: it runs
//! the existing joint wall/enthalpy/coefficient solver, then evaluates momentum
//! resistance on that SAME returned fluid trajectory. The solid callback must
//! start every inner AND outer candidate from the same accepted history. Only
//! the pressure-closed candidate is published; no speculative solid step may
//! commit external state, debit a reservoir or write actuator commands.

use fs_convection::CorrelationRobinBoundary;
use fs_couple::iqn_ils::IqnIlsConfig;
use fs_exec::Cx;
use fs_qty::Pressure;

use crate::AirflowError;
use crate::conjugate::ConjugateConfig;
use super::{LiquidFeedbackBudget, LiquidHeatPath, LiquidResistanceSection, LiquidResistanceModel,
    PressureDrivenBudget, PressureDrivenError, PressureDrivenSolution,
    checkpoint, solve_pressure_root};
use super::super::solid::CoupledLiquidFeedbackSolution;
use super::super::super::LiquidSolidResponse;

/// Explicit local, thermal, pressure and total solid-work limits.
#[derive(Debug)]
pub struct PressureSolidControls {
    /// Outer pressure closure and explicit neglected-dissipation allowance.
    pub pressure: PressureDrivenBudget,
    /// Existing wall/reference and per-region watt gates, per flow trial.
    pub coupling: ConjugateConfig,
    /// Existing regional material/coefficient fixed-point gates.
    pub local: LiquidFeedbackBudget,
    /// Existing IQN-ILS policy, not a new acceleration implementation.
    pub acceleration: IqnIlsConfig,
    /// Shared limit across ALL outer pressure trials, not reset for each trial.
    pub max_total_solid_solves: usize,
}

/// Original solid payload, heat state and pressure state at one operating point.
#[derive(Debug)]
pub struct CoupledPressureSolution<'a, R> {
    /// Source-bound pressure closure around the actual existing coupled result.
    pub operating_point: PressureDrivenSolution<CoupledLiquidFeedbackSolution<'a, R>>,
    /// Actual solid callbacks across all accepted and discarded flow trials.
    pub total_solid_solves: usize,
}

impl<'a, 'card> LiquidHeatPath<'a, 'card> {
    /// Solve pressure, material-dependent flow, coolant enthalpy and solid
    /// temperature as one bounded nested problem. There is no second solid
    /// solver and no stale coefficient or wall state from a previous flow trial.
    ///
    /// Thermal convergence and regional watt checks are still owned by
    /// `solve_coupled_iqn`. Pressure closure adds another necessary acceptance
    /// gate; it does not upgrade source accuracy, prove uniqueness, or include
    /// the omitted pressure-work/viscous heating. The explicit power allowance
    /// applies to every evaluated trial. Source pressure remains fixed-context.
    ///
    /// # Errors
    /// All existing source, thermal, pressure and cancellation refusals; a
    /// malformed assignment/configuration; or exhausted shared solid work.
    pub fn solve_pressure_coupled_iqn<R, F>(
        &self, cx: &Cx<'_>, initial_walls_k: &[f64], resistances: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
        drive: Pressure, controls: &PressureSolidControls, mut solid: F,
    ) -> Result<CoupledPressureSolution<'a, R>, PressureDrivenError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        checkpoint(cx)?;
        controls.pressure.admit(drive)?;
        controls.local.admit()?;
        self.admit_resistances(resistances)?;
        if controls.max_total_solid_solves == 0 {
            return Err(PressureDrivenError::Invalid("total solid-solve budget must be positive"));
        }
        let mut total_solid_solves = 0;
        let operating_point = solve_pressure_root(cx, drive, controls.pressure, |mass| {
            let path = self.at_pressure_trial_mass(cx, mass)?;
            let mut exhausted = false;
            let thermal = path.solve_coupled_iqn(cx, initial_walls_k, &controls.coupling,
                controls.local, controls.acceleration, |cx, rows| {
                    if total_solid_solves == controls.max_total_solid_solves {
                        exhausted = true;
                        return Err(AirflowError::InvalidConjugateInput {
                            field: "shared pressure/thermal solid-solve budget", value_bits: total_solid_solves as u64,
                        });
                    }
                    total_solid_solves += 1;
                    solid(cx, rows)
                });
            if exhausted {
                return Err(PressureDrivenError::SolidSolveBudget {
                    performed: total_solid_solves, maximum: controls.max_total_solid_solves,
                });
            }
            let thermal = thermal?;
            let pressure = self.pressure_state(cx, resistances, &thermal.fluid, controls.pressure)?;
            Ok((thermal, pressure))
        })?;
        checkpoint(cx)?;
        Ok(CoupledPressureSolution { operating_point, total_solid_solves })
    }
}
