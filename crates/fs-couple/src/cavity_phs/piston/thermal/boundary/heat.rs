//! One signed thermal-source solve shared by single and multiple gas contacts.
//! Gas caloric inversion, pressure work and momentum remain with GasPiston.
use super::*;

/// Gas and mechanical balances against the evaluated sum of boundary heat.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryHeatStep {
    /// Existing gas/momentum result and its source-bearing energy midpoint.
    pub piston: PistonStep,
    /// Total boundary heat actually supplied to the gas solve [J].
    pub exchanged_heat_j: f64,
    /// Sum of independently evaluated signed contact heats [J].
    pub evaluated_heat_j: f64,
    /// Applied gas heat minus evaluated contact heat [J].
    pub source_residual_j: f64,
    /// Gas + kinetic + exterior work - force work - prescribed gas heat
    /// - evaluated contact heat [J], including the source residual.
    pub total_energy_residual_j: f64,
    /// All complete gas/momentum proposals, including rejected source trials.
    pub exchange_trials: usize,
    /// All gas-material proposals across those trials.
    pub material_trials: usize,
}

/// One gas endpoint and the complete source evaluation at its common midpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryHeatCandidate<S> {
    piston: GasPiston,
    source: S,
    step: BoundaryHeatStep,
}
impl<S> BoundaryHeatCandidate<S> {
    /// Proposed gas, volume and momentum. The accepted input remains unchanged.
    pub const fn piston(&self) -> &GasPiston { &self.piston }
    /// Original per-contact sources; not a homogenized substitute material.
    pub const fn source(&self) -> &S { &self.source }
    /// Independently checked thermal and mechanical balance.
    pub const fn step(&self) -> &BoundaryHeatStep { &self.step }
    /// Publish with every finite counterpart; the heat is already in the gas.
    pub fn into_parts(self) -> (GasPiston, S, BoundaryHeatStep) {
        (self.piston, self.source, self.step)
    }
}

// The physical callback returns evaluated whole-step joules, not a second gas
// temperature update. Source refusals keep the original error classification.
pub(super) fn solve<S, E>(
    gas: &GasPiston, input: PistonInput, controls: PistonSolve, exchange: HeatExchangeSolve,
    mut resolve: impl FnMut(&ResolvedThermochemicalGasState) -> Result<(f64, S), ThermalPistonError<E>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BoundaryHeatCandidate<S>, ThermalPistonError<E>> {
    poll(&mut cancelled)?;
    if !input.heat_j.is_finite() {
        return Err(ThermalPistonError::Invalid("finite prescribed gas heat required"));
    }
    let mut material_trials = 0usize;
    let candidate = feedback::solve(exchange, |q, trials, cancelled| {
        let heat_j = input.heat_j + q;
        if !heat_j.is_finite() || (q != 0.0 && heat_j == input.heat_j) {
            return Err(ThermalPistonError::Invalid("boundary heat vanishes or overflows in gas input"));
        }
        let mut next = gas.clone();
        let piston = next.step(PistonInput { heat_j, ..input }, controls, &mut *cancelled)?;
        material_trials += piston.material_trials;
        poll(cancelled)?;
        let (evaluated, source) = resolve(&piston.midpoint)?;
        let residual = q - evaluated;
        let total = piston.gas.stored_energy_change_j + piston.kinetic_energy_change_j
            + piston.exterior_pressure_work_j - piston.external_force_work_j
            - input.heat_j - evaluated;
        if !evaluated.is_finite() || !residual.is_finite() || !total.is_finite() {
            return Err(ThermalPistonError::Invalid("unrepresentable thermal boundary source or work"));
        }
        let step = BoundaryHeatStep { piston, exchanged_heat_j: q, evaluated_heat_j: evaluated,
            source_residual_j: residual, total_energy_residual_j: total,
            exchange_trials: trials, material_trials };
        Ok((BoundaryHeatCandidate { piston: next, source, step }, residual))
    }, &mut cancelled)?;
    if candidate.step.total_energy_residual_j.abs() > exchange.energy_tolerance_j {
        return Err(ThermalPistonError::Invalid("evaluated boundary heat does not close gas and mechanical energy"));
    }
    poll(&mut cancelled)?;
    Ok(candidate)
}

pub(crate) fn conductive_heat<E>(dt: f64, g: f64, boundary_k: f64, gas_k: f64)
    -> Result<f64, ThermalPistonError<E>>
{
    let delta = boundary_k - gas_k;
    let power = g * delta;
    let heat = dt * power;
    if !g.is_finite() || g < 0.0 || !power.is_finite() || !heat.is_finite()
        || (g != 0.0 && delta != 0.0 && power == 0.0) || (power != 0.0 && heat == 0.0) {
        return Err(ThermalPistonError::Invalid("unrepresentable thermal boundary source or work"));
    }
    Ok(heat)
}
