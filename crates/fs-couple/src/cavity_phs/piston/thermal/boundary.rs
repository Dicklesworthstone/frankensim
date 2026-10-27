//! A thermal port on the existing gas/momentum owner, without a second body store.
//! The caller supplies the other side's temperature for THIS nonlinear trial.
//! A finite network must debit the returned signed heat and accept both endpoints.

/// Two thermal ports sharing one paired-gas momentum solve.
pub mod double_acting;
/// Finite circuit and thermal network feedback on this same boundary.
pub mod network;
/// Simultaneous signed heat contacts sharing one gas energy midpoint.
pub mod ports;
pub(crate) mod heat;

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
        let candidate = heat::solve(self, input, controls, exchange, |gas| {
            let sample = resolve(gas).map_err(ThermalPistonError::Interface)?;
            let evaluated = heat::conductive_heat(input.duration_s, sample.conductance_w_k,
                boundary_temperature_k, gas.state().temperature)?;
            Ok((evaluated, sample))
        }, cancelled)?;
        let (piston, sample, report) = candidate.into_parts();
        let step = PistonBoundaryStep { piston: report.piston, boundary_temperature_k,
            conductance_w_k: sample.conductance_w_k, exchanged_heat_j: report.exchanged_heat_j,
            evaluated_heat_j: report.evaluated_heat_j, source_residual_j: report.source_residual_j,
            total_energy_residual_j: report.total_energy_residual_j,
            exchange_trials: report.exchange_trials, material_trials: report.material_trials };
        Ok(PistonBoundaryCandidate { piston, source: sample.source, step })
    }
}

#[cfg(test)]
mod tests;
