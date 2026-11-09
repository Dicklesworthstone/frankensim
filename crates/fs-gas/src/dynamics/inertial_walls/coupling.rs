//! One mechanical solve for native gas fields of any spatial dimension.

use super::{
    InertialWallLoad, InertialWallsControl, InertialWallsStep, MechanicalExchange, Piston,
    PistonError, StepControl, Time, mechanical_report, positive, solve_with, validate_controls,
};

/// Original gas exchanges from one complete trial at a fixed duration.
///
/// Both walls use the increasing axial coordinate, ordered [left, right].
/// The boundary entries are signed transfers INTO THE GAS, not support loads
/// or outward fluxes. For example, positive pressure gives a positive left
/// gas impulse and a negative right gas impulse. Each work entry is that same
/// impulse paired with the corresponding trial wall speed. Spatial owners
/// must sum their actual stage-integrated face exchanges in this convention.
///
/// The gas owner retains admissibility, mass/species conservation and geometric
/// conservation. These extensive mechanical exchanges do not replace its
/// original report, source model or state acceptance. Interior source impulse
/// and work, if present, must be subtracted from the two change fields before
/// the wall-only balance is evaluated and retained separately in `gas`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallGasExchange<R> {
    /// The complete native gas report for this exact trial.
    pub gas: R,
    /// Actual axial gas momentum change, less independent body sources [N s].
    pub momentum_change_n_s: f64,
    /// Actual gas energy change, less independent body sources [J].
    pub energy_change_j: f64,
    /// Actual signed [left, right] impulses received by the gas [N s].
    pub boundary_impulse_n_s: [f64; 2],
    /// Actual signed [left, right] work received by the gas [J].
    pub boundary_work_j: [f64; 2],
    /// Represented endpoint-face travel along the increasing axis [m].
    pub displacement_m: [f64; 2],
}

impl InertialWallsControl {
    /// Admit the original wall data and step budget before preparing gas trials.
    ///
    /// This is the same admission used by the one-dimensional owners. Geometry,
    /// face area and source-dependent wave speeds remain with the spatial owner.
    /// Invalid inputs are rejected without invoking a gas proposal.
    pub fn validate(
        self,
        walls: [Piston; 2],
        maximum_duration: Time,
        gas_controls: StepControl,
        loads: [InertialWallLoad; 2],
    ) -> Result<(), PistonError> {
        validate_controls(walls, maximum_duration, gas_controls, loads, self)
    }
}

/// Solve both finite wall midpoint velocities using the original gas operator.
///
/// This exposes the same bounded two-variable Newton solve and independent
/// impulse/work gates used by the native one-dimensional fields. It does not
/// construct a reduced gas model. The callback must rerun its original spatial
/// operator from the SAME accepted state at the supplied fixed `duration`,
/// using the two proposed constant wall speeds. It must not publish a trial.
/// Native scratch can be reused: a successful result ALWAYS corresponds to
/// the most recent callback, including after derivative and line-search work.
///
/// The caller owns relative-wave/contraction CFL, retries of the whole interval,
/// geometry, transverse supports and final atomic publication. It must poll
/// cancellation in each gas proposal and supply a `terminal` predicate that
/// prevents cancellation, allocation and other non-retryable owner errors from
/// being swallowed by backtracking. The original error type is retained.
///
/// The returned wall velocities have representable finite kinetic energy.
/// No successful mechanical solve allocates or changes a wall, gas field,
/// geometry or clock. `gas_trials` counts every callback in this solve;
/// an outer retry owner may add counts from discarded attempts.
///
/// # Errors
/// Refuses invalid controls/area, nonfinite original exchanges, unrepresentable
/// wall storage, exhausted root work, failed independent balances or a callback
/// error. Numerical acceptance is not a general FSI stability certificate.
#[allow(clippy::too_many_arguments)]
pub fn solve_coupled<R, E: From<PistonError>>(
    walls: [Piston; 2],
    duration: Time,
    gas_controls: StepControl,
    loads: [InertialWallLoad; 2],
    controls: InertialWallsControl,
    area_m2: f64,
    mut trial: impl FnMut([f64; 2]) -> Result<WallGasExchange<R>, E>,
    terminal: impl Fn(&E) -> bool,
) -> Result<InertialWallsStep<R>, E> {
    controls.validate(walls, duration, gas_controls, loads)?;
    if !positive(area_m2)
        || (0..2).any(|i| {
            !(loads[i].external_force_n - [-1.0, 1.0][i] * area_m2 * loads[i].exterior_pressure_pa)
                .is_finite()
        })
    {
        return Err(PistonError::Invalid("invalid finite-wall area or external force").into());
    }
    let mut gas_trials = 0;
    let mut report = solve_with(
        walls.map(|wall| wall.velocity_m_s()),
        controls,
        |speed| {
            gas_trials += 1;
            let exchange = trial(speed)?;
            if [exchange.momentum_change_n_s, exchange.energy_change_j]
                .into_iter()
                .chain(exchange.boundary_impulse_n_s)
                .chain(exchange.boundary_work_j)
                .chain(exchange.displacement_m)
                .any(|value| !value.is_finite())
            {
                return Err(PistonError::Invalid("nonfinite original wall/gas exchange").into());
            }
            let report = mechanical_report(
                exchange.gas,
                MechanicalExchange {
                    duration_s: duration.value(),
                    momentum_change: exchange.momentum_change_n_s,
                    energy_change: exchange.energy_change_j,
                    boundary_impulse: exchange.boundary_impulse_n_s,
                    boundary_work: exchange.boundary_work_j,
                },
                walls,
                loads,
                speed,
                exchange.displacement_m,
                area_m2,
            )?;
            for (wall, velocity) in walls.iter().zip(report.velocity_after_m_s) {
                Piston::new(wall.mass_kg(), velocity)?;
            }
            Ok(report)
        },
        terminal,
    )?;
    report.gas_trials = gas_trials;
    Ok(report)
}

#[cfg(test)]
mod tests;
