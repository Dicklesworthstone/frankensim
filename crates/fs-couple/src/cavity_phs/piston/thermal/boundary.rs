//! A thermal port on the existing gas/momentum owner, without a second body store.
//! The caller supplies the other side's temperature for THIS nonlinear trial.
//! A finite network must debit the returned signed heat and accept both endpoints.

use super::{ConductanceSample, GasPiston, HeatExchangeSolve, PistonInput,
    PistonSolve, PistonStep, ResolvedThermochemicalGasState, ThermalPistonError,
    feedback, poll, positive};

/// Actual source and motion balances for one admitted boundary proposal.
#[derive(Clone, Debug, PartialEq)]
pub struct PistonBoundaryStep {
    /// Original coupled gas/momentum result, retaining its material midpoint.
    pub piston: PistonStep,
    /// Other side's temperature in this proposal [K]. Not another stored state.
    pub boundary_temperature_k: f64,
    /// Actual midpoint conductance [W/K]. Zero explicitly isolates the port.
    pub conductance_w_k: f64,
    /// Signed joules into the gas, to be debited ONCE from a finite counterpart.
    pub exchanged_heat_j: f64,
    /// Duration times the actual resolved boundary heat rate [J].
    pub evaluated_heat_j: f64,
    /// Applied minus evaluated boundary heat [J].
    pub source_residual_j: f64,
    /// Gas + kinetic + exterior work - force work - prescribed heat
    /// - EVALUATED boundary heat [J]. Includes the source residual.
    pub total_energy_residual_j: f64,
    /// All complete gas/momentum proposals, including rejected heat trials.
    pub exchange_trials: usize,
    /// All nested gas-material proposals in this boundary solve.
    pub material_trials: usize,
}

/// One unchanged input yields one tentative endpoint and its actual source.
#[derive(Clone, Debug, PartialEq)]
pub struct PistonBoundaryCandidate<S> {
    piston: GasPiston,
    source: S,
    step: PistonBoundaryStep,
}
impl<S> PistonBoundaryCandidate<S> {
    /// Proposed gas, geometry and momentum, never yet installed in the caller.
    pub const fn piston(&self) -> &GasPiston { &self.piston }
    /// Concrete interface evaluation at the accepted gas midpoint.
    pub const fn source(&self) -> &S { &self.source }
    /// Physical heat, source consistency and original motion account.
    pub const fn step(&self) -> &PistonBoundaryStep { &self.step }
    /// Publish with the other physical owners. Do not apply the heat again.
    pub fn into_parts(self) -> (GasPiston, S, PistonBoundaryStep) {
        (self.piston, self.source, self.step)
    }
}

impl GasPiston {
    /// Solve gas energy, momentum and an external thermal port simultaneously.
    /// Q = dt G(gas_mid) (T_boundary - T_gas_mid). The existing gas, piston and
    /// bounded scalar feedback owners retain ALL caloric and motion equations.
    /// `input.heat_j` is a separate prescribed source, not this port's heat.
    ///
    /// For a reservoir the supplied temperature is held during this step. For a
    /// finite counterpart, call this from its nonlinear residual at that trial's
    /// temperature, debit exactly `exchanged_heat_j`, and publish together.
    /// Do not advance a live source during callbacks or treat this candidate as
    /// proof that the other side conserved energy. The interface callback must
    /// be repeatable and retains its original concrete material result.
    ///
    /// Every trial begins at the SAME accepted gas and velocity. Failure,
    /// cancellation or source-domain refusal cannot partially update them.
    /// Work is bounded by (exchange.maximum_iterations+2) complete piston
    /// proposals with their existing nested limits. No uniqueness, monotonicity,
    /// time-error, phase-path or unconditional passivity claim is made.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_boundary_candidate<S, E>(
        &self, input: PistonInput, controls: PistonSolve, exchange: HeatExchangeSolve,
        boundary_temperature_k: f64,
        mut resolve: impl FnMut(&ResolvedThermochemicalGasState) -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PistonBoundaryCandidate<S>, ThermalPistonError<E>> {
        poll(&mut cancelled)?;
        if !positive(boundary_temperature_k) || !input.heat_j.is_finite() {
            return Err(ThermalPistonError::Invalid("positive boundary temperature and finite prescribed heat required"));
        }
        let mut material_trials = 0usize;
        let candidate = feedback::solve(exchange, |q, trials, cancelled| {
            let heat_j = input.heat_j + q;
            if !heat_j.is_finite() || (q != 0.0 && heat_j == input.heat_j) {
                return Err(ThermalPistonError::Invalid("boundary heat vanishes or overflows in gas input"));
            }
            let mut next = self.clone();
            let piston = next.step(PistonInput { heat_j, ..input }, controls, &mut *cancelled)?;
            material_trials += piston.material_trials;
            poll(cancelled)?;
            let sample = resolve(&piston.midpoint).map_err(ThermalPistonError::Interface)?;
            let g = sample.conductance_w_k;
            let delta = boundary_temperature_k - piston.midpoint.state().temperature;
            let power = g * delta;
            let evaluated = input.duration_s * power;
            let residual = q - evaluated;
            let total = piston.gas.stored_energy_change_j + piston.kinetic_energy_change_j
                + piston.exterior_pressure_work_j - piston.external_force_work_j
                - input.heat_j - evaluated;
            if !g.is_finite() || g < 0.0 || !power.is_finite() || !evaluated.is_finite()
                || !residual.is_finite() || !total.is_finite()
                || (g != 0.0 && delta != 0.0 && power == 0.0)
                || (power != 0.0 && evaluated == 0.0) {
                return Err(ThermalPistonError::Invalid("unrepresentable thermal boundary source or work"));
            }
            let step = PistonBoundaryStep { piston, boundary_temperature_k,
                conductance_w_k: g, exchanged_heat_j: q, evaluated_heat_j: evaluated,
                source_residual_j: residual, total_energy_residual_j: total,
                exchange_trials: trials, material_trials };
            Ok((PistonBoundaryCandidate { piston: next, source: sample.source, step }, residual))
        }, &mut cancelled)?;
        if candidate.step.total_energy_residual_j.abs() > exchange.energy_tolerance_j {
            return Err(ThermalPistonError::Invalid("evaluated boundary heat does not close gas and mechanical energy"));
        }
        poll(&mut cancelled)?;
        Ok(candidate)
    }
}

#[cfg(test)]
mod tests;
