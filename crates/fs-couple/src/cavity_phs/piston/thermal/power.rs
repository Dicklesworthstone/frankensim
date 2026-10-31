//! State-dependent body power coupled to heat exchange, phase and piston motion.
//!
//! Solve W = dt P(gas_mid, body_mid), where each proposal runs the existing
//! reciprocal thermal/momentum solve from the SAME accepted state. W is added
//! only to the body's internal source. Interface heat remains equal/opposite.
//! No constitutive power law, temperature increment or new time stepper lives here.

/// Source-resolved resistive power, using the existing conductor owner.
pub mod joule;

#[cfg(test)]
use super::ResolvedThermochemicalGasState;
use super::{
    ClosedGasParcel, ConductanceSample, EquilibriumPhaseState, GasPiston, HeatExchangeSolve,
    PistonGas, PistonInput, PistonSolve, ThermalPistonCandidate, ThermalPistonError,
    UniformEnthalpyStepInput, feedback, poll,
};

/// A signed instantaneous power and the original concrete source evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct PowerSample<S> {
    /// Watts INTO the finite thermal body. Negative values withdraw energy.
    pub power_w: f64,
    /// Original material/circuit/source result, without type erasure.
    pub source: S,
}

/// Actual source consistency and independently checked total energy balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PowerFeedbackStep {
    /// Whole-step source energy used by the accepted body/gas/motion solve [J].
    pub applied_energy_j: f64,
    /// Duration times the source power re-evaluated at that solve's midpoint [J].
    pub evaluated_energy_j: f64,
    /// Applied minus evaluated energy [J], never relative to formation energy.
    pub source_residual_j: f64,
    /// Gas + body + kinetic change + exterior work - external force work
    /// - both prescribed heats - evaluated source energy [J].
    pub total_energy_residual_j: f64,
    /// Complete outer source proposals, including zero and bracket evaluations.
    pub source_trials: usize,
    /// Total inner heat-exchange proposals over all outer trials.
    pub exchange_trials: usize,
    /// Total gas/momentum material proposals over all inner/outer trials.
    pub material_trials: usize,
}

/// Both physical endpoints and the source law at their common energy midpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct PoweredThermalCandidate<I, P, G: PistonGas = ClosedGasParcel> {
    coupled: ThermalPistonCandidate<I, G>,
    power: PowerSample<P>,
    step: PowerFeedbackStep,
}
impl<I, P, G: PistonGas> PoweredThermalCandidate<I, P, G> {
    /// The original reciprocal thermal/gas/momentum candidate, not a second state.
    pub const fn coupled(&self) -> &ThermalPistonCandidate<I, G> {
        &self.coupled
    }
    /// Exact source sample used in the physical consistency test.
    pub const fn power(&self) -> &PowerSample<P> { &self.power }
    /// Actual source and total-energy residuals and bounded work counts.
    pub const fn step(&self) -> &PowerFeedbackStep { &self.step }
    /// Publish both endpoints together. Source energy is ALREADY in their
    /// storage; neither the returned power nor interface heat is added again.
    pub fn into_parts(
        self,
    ) -> (
        ThermalPistonCandidate<I, G>,
        PowerSample<P>,
        PowerFeedbackStep,
    ) {
        (self.coupled, self.power, self.step)
    }
}

impl<G: PistonGas> GasPiston<G> {
    /// Couple a repeatable state-dependent body power to the existing thermal
    /// and mechanical equations. Both callbacks receive the ACTUAL gas and body
    /// energy midpoints; neither may debit a live external energy store.
    /// Their original concrete source values are retained independently.
    ///
    /// `body.internal_heat_j` and `input.heat_j` remain prescribed heats, in
    /// addition to this power. `source.heat_bracket_j` bounds the extra signed
    /// whole-step body energy and contains zero. `exchange` independently bounds
    /// internal transfer. Every queried state must be admitted by all owners.
    ///
    /// The same bounded secant/bisection owner solves both feedback equations.
    /// Work is bounded by (source.maximum_iterations+2) complete thermal solves,
    /// each with its own exchange and piston limits. Failure is not evidence of
    /// physical runaway. No uniqueness or unconditional stability is claimed.
    /// Cancellation and failures leave the gas, body, geometry and momentum
    /// unchanged; the source callback's external side effects are not reversible.
    #[allow(clippy::too_many_arguments)]
    pub fn powered_thermal_candidate<I, P, E>(
        &self,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        source: HeatExchangeSolve,
        mut conductance: impl FnMut(&G::State, EquilibriumPhaseState) -> Result<ConductanceSample<I>, E>,
        mut power: impl FnMut(&G::State, EquilibriumPhaseState) -> Result<PowerSample<P>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PoweredThermalCandidate<I, P, G>, ThermalPistonError<E>> {
        poll(&mut cancelled)?;
        feedback::admit(source)?;
        let mut exchange_trials = 0usize;
        let mut material_trials = 0usize;
        let candidate = feedback::solve(source, |energy, trials, cancelled| {
            let prescribed = body.internal_heat_j + energy;
            if !prescribed.is_finite() || (energy != 0.0 && prescribed == body.internal_heat_j) {
                return Err(ThermalPistonError::Invalid("body source energy disappears or overflows"));
            }
            let coupled = self.thermal_candidate(input,
                UniformEnthalpyStepInput { internal_heat_j: prescribed, ..body },
                controls, exchange, |gas, body| conductance(gas, body), &mut *cancelled)?;
            let report = coupled.step();
            exchange_trials += report.exchange_trials;
            // The callback above counts every exchange proposal's material work
            // through thermal_candidate's report, not only the accepted proposal.
            material_trials += report.material_trials;
            poll(cancelled)?;
            let sample = power(&report.piston.midpoint, report.thermal_midpoint)
                .map_err(ThermalPistonError::Interface)?;
            let evaluated = sample.power_w * input.duration_s;
            let residual = energy - evaluated;
            if !sample.power_w.is_finite() || !evaluated.is_finite() || !residual.is_finite()
                || (sample.power_w != 0.0 && evaluated == 0.0) {
                return Err(ThermalPistonError::Invalid("body power or integrated source is unrepresentable"));
            }
            let piston = &report.piston;
            let total = ((piston.gas.stored_energy_change_j + piston.kinetic_energy_change_j)
                + report.thermal_energy_change_j) + piston.exterior_pressure_work_j
                - piston.external_force_work_j - input.heat_j - body.internal_heat_j - evaluated;
            if !total.is_finite() {
                return Err(ThermalPistonError::Invalid("coupled source energy balance overflow"));
            }
            let step = PowerFeedbackStep { applied_energy_j: energy, evaluated_energy_j: evaluated,
                source_residual_j: residual, total_energy_residual_j: total,
                source_trials: trials, exchange_trials, material_trials };
            Ok((PoweredThermalCandidate { coupled, power: sample, step }, residual))
        }, &mut cancelled)?;
        if candidate.step.total_energy_residual_j.abs() > source.energy_tolerance_j {
            return Err(ThermalPistonError::Invalid("actual midpoint power does not close the combined joule balance"));
        }
        poll(&mut cancelled)?;
        Ok(candidate)
    }
}

#[cfg(test)]
mod tests;
