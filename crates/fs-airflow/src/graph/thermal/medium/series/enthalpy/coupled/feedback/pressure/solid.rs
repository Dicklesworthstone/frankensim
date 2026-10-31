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
    /// The opt-in dissipative solve requires this allowance to be zero.
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
        drive: Pressure, controls: &PressureSolidControls, solid: F,
    ) -> Result<CoupledPressureSolution<'a, R>, PressureDrivenError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        self.solve_pressure_coupled_iqn_inner(
            cx, initial_walls_k, resistances, drive, controls, false, solid,
        )
    }

    /// Close pressure, hydraulic heating and the actual solid temperature together.
    ///
    /// Every wall iterate uses `march_with_dissipation`: the same source-bound
    /// pressure loss heats the coolant, changes its viscosity and convection,
    /// and changes the wall-only Robin rows seen by the original solid solver.
    /// The accepted hydraulic state comes from the independent re-evaluation
    /// at that solid's returned walls, never a loss-free or stale trajectory.
    /// The original IQN iteration, regional wall-power gates and pressure root
    /// remain the sole numerical owners, with one shared solid-callback budget.
    ///
    /// Requires zero omitted-power allowance, constant source density, equal
    /// flow areas and the original dissipative march's equal-elevation, unchanged
    /// kinetic-energy and complete fluid-deposition assumptions. Pump inefficiency
    /// and heat retained by a porous skeleton are outside this model. Callbacks
    /// must start every trial from the same accepted solid history.
    ///
    /// # Errors
    /// Preserves all original source, geometry, thermal, pressure, cancellation
    /// and work refusals; a nonzero omitted-power allowance is invalid.
    pub fn solve_pressure_coupled_iqn_with_dissipation<R, F>(
        &self, cx: &Cx<'_>, initial_walls_k: &[f64], resistances: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
        drive: Pressure, controls: &PressureSolidControls, solid: F,
    ) -> Result<CoupledPressureSolution<'a, R>, PressureDrivenError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        self.solve_pressure_coupled_iqn_inner(
            cx, initial_walls_k, resistances, drive, controls, true, solid,
        )
    }

    fn solve_pressure_coupled_iqn_inner<R, F>(
        &self, cx: &Cx<'_>, initial_walls_k: &[f64], resistances: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
        drive: Pressure, controls: &PressureSolidControls, deposit_hydraulic_heat: bool, mut solid: F,
    ) -> Result<CoupledPressureSolution<'a, R>, PressureDrivenError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        checkpoint(cx)?;
        controls.pressure.admit(drive)?;
        controls.local.admit()?;
        self.admit_resistances(resistances)?;
        if deposit_hydraulic_heat && controls.pressure.maximum_neglected_hydraulic_power_w != 0.0 {
            return Err(PressureDrivenError::Invalid(
                "fully dissipative heat path requires zero omitted-power allowance",
            ));
        }
        if controls.max_total_solid_solves == 0 {
            return Err(PressureDrivenError::Invalid("total solid-solve budget must be positive"));
        }
        let mut total_solid_solves = 0;
        let operating_point = solve_pressure_root(cx, drive, controls.pressure, |mass| {
            let path = self.at_pressure_trial_mass(cx, mass)?;
            let mut exhausted = false;
            let thermal = path.solve_coupled_iqn_march(cx, initial_walls_k, &controls.coupling,
                controls.local, controls.acceleration, |cx, rows| {
                    if total_solid_solves == controls.max_total_solid_solves {
                        exhausted = true;
                        return Err(AirflowError::InvalidConjugateInput {
                            field: "shared pressure/thermal solid-solve budget", value_bits: total_solid_solves as u64,
                        });
                    }
                    total_solid_solves += 1;
                    solid(cx, rows)
                }, |cx, walls| {
                    if deposit_hydraulic_heat {
                        let march = path.march_with_dissipation(cx, walls, resistances, controls.local)?;
                        Ok::<_, PressureDrivenError>((march.thermal, Some(march.hydraulics)))
                    } else {
                        Ok((path.march(cx, walls, controls.local)?, None))
                    }
                });
            if exhausted {
                return Err(PressureDrivenError::SolidSolveBudget {
                    performed: total_solid_solves, maximum: controls.max_total_solid_solves,
                });
            }
            let (thermal, _, hydraulics) = thermal?;
            let pressure = match hydraulics {
                Some(hydraulics) => hydraulics,
                None => path.pressure_state(cx, resistances, &thermal.fluid, controls.pressure)?,
            };
            Ok((thermal, pressure))
        })?;
        checkpoint(cx)?;
        Ok(CoupledPressureSolution { operating_point, total_solid_solves })
    }
}
