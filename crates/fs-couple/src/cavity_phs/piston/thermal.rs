//! Reciprocal finite-body heat exchange on the existing gas/piston step.
//!
//! The solid/liquid material owns h -> (T, phase); the gas owns u, EOS and
//! acoustic/transport response; `GasPiston::step` owns motion and pressure work.
//! One signed Q goes into the gas and OUT of the thermal body. Source laws are
//! evaluated on their actual energy midpoints, never averaged temperatures.
//! The thermal body uses the existing negligible-pV-work enthalpy storage model.
//! Its volume is fixed during this proposal; expansion/deformation is not inferred.

use super::{GasPiston, PistonError, PistonInput, PistonSolve, PistonStep};
use fs_material::phase::{EquilibriumPhaseState, PhaseStateError,
    UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

/// Whole-step exchange bracket and absolute residual/work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatExchangeSolve {
    /// Signed joules into the gas, containing zero. Every queried proposal must
    /// be supported by BOTH materials and the existing piston volume bracket.
    pub heat_bracket_j: [f64; 2],
    /// Maximum |Q - dt G(T_body_mid - T_gas_mid)| [J].
    pub source_tolerance_j: f64,
    /// Maximum absolute body and combined energy residual [J].
    pub energy_tolerance_j: f64,
    /// Interior safeguarded secant/bisection trials, in 1..=128. At most two
    /// extra trials check zero exchange and one bracket endpoint. Each trial
    /// has the separate, explicit work limits of `PistonSolve`.
    pub maximum_iterations: usize,
}

/// A conductance and the concrete source evaluation that supplied it.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductanceSample<S> {
    /// Nonnegative total conductance [W/K]; zero is explicit insulation.
    pub conductance_w_k: f64,
    /// Retained material/interface evaluation; never replaced with a name preset.
    pub source: S,
}

/// Numerical, material and interface failures retain their original owner.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalPistonError<E> {
    /// Invalid geometry, clocks, source energy or numerical controls.
    Invalid(&'static str),
    /// Original gas/mechanical failure.
    Piston(PistonError),
    /// Original enthalpy/phase chart refused a state.
    Thermal(PhaseStateError),
    /// Caller-owned interface material law refused the midpoint.
    Interface(E),
    /// Cancellation before returning an accepted joint candidate.
    Cancelled,
    /// No sign change in the caller's admitted exchange bracket.
    OutsideBracket { lower_residual_j: f64, upper_residual_j: f64 },
    /// Work exhausted or floating-point stagnation without a physical root.
    NotConverged { iterations: usize, residual_j: f64 },
}
impl<E> From<PistonError> for ThermalPistonError<E> {
    fn from(e: PistonError) -> Self {
        if e == PistonError::Cancelled { Self::Cancelled } else { Self::Piston(e) }
    }
}
impl<E> From<PhaseStateError> for ThermalPistonError<E> {
    fn from(e: PhaseStateError) -> Self { Self::Thermal(e) }
}
impl<E: core::fmt::Display> core::fmt::Display for ThermalPistonError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "gas/body thermal coupling: {s}"),
            Self::Piston(e) => write!(f, "gas/body thermal coupling: {e}"),
            Self::Thermal(e) => write!(f, "gas/body thermal material: {e}"),
            Self::Interface(e) => write!(f, "gas/body thermal interface: {e}"),
            Self::Cancelled => write!(f, "gas/body thermal coupling cancelled"),
            Self::OutsideBracket { lower_residual_j, upper_residual_j } =>
                write!(f, "heat exchange not bracketed: {lower_residual_j}, {upper_residual_j} J"),
            Self::NotConverged { iterations, residual_j } =>
                write!(f, "heat exchange did not converge in {iterations} trials: {residual_j} J"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ThermalPistonError<E> {}

/// Actual accepted balances; no independent entropy-production claim is implied.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalPistonStep {
    /// Original coupled gas/momentum result, including the gas midpoint sources.
    pub piston: PistonStep,
    /// Actual material state at midpoint BODY ENTHALPY, including liquid fraction.
    pub thermal_midpoint: EquilibriumPhaseState,
    /// Accepted total interface conductance [W/K].
    pub conductance_w_k: f64,
    /// Signed heat from the thermal body into the gas [J].
    pub exchanged_heat_j: f64,
    /// Q minus the independently re-evaluated midpoint heat law [J].
    pub source_residual_j: f64,
    /// Actual m_body*(h1-h0) [J].
    pub thermal_energy_change_j: f64,
    /// Actual body storage plus Q minus its internal heat source [J].
    pub thermal_energy_residual_j: f64,
    /// Gas internal + piston kinetic + body enthalpy + exterior work - applied
    /// force work - the TWO prescribed heat sources [J]. Q cancels internally.
    pub total_energy_residual_j: f64,
    /// Actual calls to the complete gas/piston and midpoint interface owners.
    pub exchange_trials: usize,
}

/// Both accepted physical states in one candidate. Neither input is mutated.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalPistonCandidate<S> {
    piston: GasPiston,
    thermal: UniformEnthalpyStep<S>,
    step: ThermalPistonStep,
}
impl<S> ThermalPistonCandidate<S> {
    /// Tentative gas, volume and momentum; old sources/history remain intact.
    pub const fn piston(&self) -> &GasPiston { &self.piston }
    /// Canonical thermal proposal with original concrete interface source.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> { &self.thermal }
    /// Independently checked coupled balances.
    pub const fn step(&self) -> &ThermalPistonStep { &self.step }
    /// Publish both states together in the caller's existing physical owner.
    pub fn into_parts(self) -> (GasPiston, UniformEnthalpyStep<S>, ThermalPistonStep) {
        (self.piston, self.thermal, self.step)
    }
}

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), ThermalPistonError<E>> {
    if cancelled() { Err(ThermalPistonError::Cancelled) } else { Ok(()) }
}

impl GasPiston {
    /// Simultaneously solve motion, gas energy, finite-body enthalpy and heat
    /// exchange. `input.heat_j` and `body.internal_heat_j` are independent
    /// prescribed sources; neither includes the inter-body exchange Q.
    ///
    /// `resolve(gas_midpoint, body_midpoint)` returns the actual nonnegative
    /// interface conductance and its original concrete source. The callback may
    /// use both materials, phase and caller-owned geometry/source selections.
    /// It must be repeatable; irreversible callback effects are not rolled back.
    ///
    /// The shared enthalpy carrier retains the curve and admitted initial state.
    /// Material names imply no heat capacity, density, conductivity or phase law.
    /// Thermal mass is conserved; geometry/density consistency, uniform-body
    /// validity and negligible body pV/deformation work remain caller assumptions.
    /// This body may be a separate chamber wall or thermal store, not necessarily
    /// the translating piston. A phase change does not silently change mechanics.
    ///
    /// Root finding uses actual source residuals, not bracket width or a frozen
    /// conductance. No monotonicity or unconditional stability is claimed for
    /// arbitrary state-dependent conductance. Midpoint may overshoot thermal
    /// equilibrium at large timesteps; callers must check temporal resolution.
    /// Any refusal/cancellation returns neither partial state. The caller must
    /// publish both candidate endpoints together, not debit a live reservoir in
    /// the callback or separately re-apply its returned boundary heat.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_candidate<S, E>(
        &self, input: PistonInput, body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve, exchange: HeatExchangeSolve,
        mut resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalPistonCandidate<S>, ThermalPistonError<E>> {
        poll(&mut cancelled)?;
        let [mut lo, mut hi] = exchange.heat_bracket_j;
        if !lo.is_finite() || !hi.is_finite() || lo >= hi || lo > 0.0 || hi < 0.0
            || !positive(exchange.source_tolerance_j) || !positive(exchange.energy_tolerance_j)
            || !(1..=128).contains(&exchange.maximum_iterations)
            || !positive(body.mass_kg) || !positive(body.volume_m3)
            || !positive(body.surface_area_m2) || !positive(body.duration_s)
            || body.duration_s != input.duration_s || !body.internal_heat_j.is_finite()
            || !input.heat_j.is_finite()
        { return Err(ThermalPistonError::Invalid("invalid heat bounds, clocks, body geometry/mass or controls")); }
        let h0 = body.initial.specific_enthalpy_j_kg();
        if body.curve.state_at_specific_enthalpy(h0)? != body.initial {
            return Err(ThermalPistonError::Invalid("initial thermal state belongs to a different source chart"));
        }
        let mut evaluate = |q: f64, trials| -> Result<ThermalPistonCandidate<S>, ThermalPistonError<E>> {
            poll(&mut cancelled)?;
            let body_heat = body.internal_heat_j - q;
            let dh = body_heat / body.mass_kg;
            let h1 = h0 + dh;
            let hm = h0 + 0.5*dh;
            let gas_heat = input.heat_j + q;
            if !body_heat.is_finite() || !dh.is_finite() || !h1.is_finite() || !hm.is_finite()
                || !gas_heat.is_finite() || (body_heat != 0.0 && (dh == 0.0 || h1 == h0))
                || (q != 0.0 && (body_heat == body.internal_heat_j || gas_heat == input.heat_j))
            { return Err(ThermalPistonError::Invalid("heat exchange disappears or overflows in conserved energy coordinates")); }
            let thermal_midpoint = body.curve.state_at_specific_enthalpy(hm)?;
            let endpoint = body.curve.state_at_specific_enthalpy(h1)?;
            let mut next = self.clone();
            let piston = next.step(PistonInput { heat_j: gas_heat, ..input }, controls, &mut cancelled)?;
            poll(&mut cancelled)?;
            let sample = resolve(&piston.midpoint, thermal_midpoint).map_err(ThermalPistonError::Interface)?;
            let g = sample.conductance_w_k;
            let delta = thermal_midpoint.temperature_k() - piston.midpoint.state().temperature;
            let power = g*delta;
            let heat = input.duration_s*power;
            let residual = q-heat;
            if !g.is_finite() || g < 0.0 || !power.is_finite() || !heat.is_finite()
                || !residual.is_finite() || (g != 0.0 && delta != 0.0 && power == 0.0)
                || (power != 0.0 && heat == 0.0)
            { return Err(ThermalPistonError::Invalid("interface heat law is negative, nonfinite or unrepresentable")); }
            let stored = body.mass_kg*(endpoint.specific_enthalpy_j_kg()-h0);
            let body_residual = (stored+q)-body.internal_heat_j;
            let total = ((piston.gas.stored_energy_change_j+piston.kinetic_energy_change_j)+stored)
                +piston.exterior_pressure_work_j-piston.external_force_work_j
                -input.heat_j-body.internal_heat_j;
            if !stored.is_finite() || !body_residual.is_finite() || !total.is_finite()
                || body_residual.abs() > exchange.energy_tolerance_j
                || total.abs() > exchange.energy_tolerance_j
            { return Err(ThermalPistonError::Invalid("finite-body and coupled energy balances do not close")); }
            Ok(ThermalPistonCandidate {
                piston: next,
                thermal: UniformEnthalpyStep { state: endpoint, external_heat_j: -q,
                    energy_residual_tolerance_j: exchange.energy_tolerance_j, report: sample.source },
                step: ThermalPistonStep { piston, thermal_midpoint, conductance_w_k: g,
                    exchanged_heat_j: q, source_residual_j: residual,
                    thermal_energy_change_j: stored, thermal_energy_residual_j: body_residual,
                    total_energy_residual_j: total, exchange_trials: trials },
            })
        };
        let first = evaluate(0.0, 1)?;
        let f0 = first.step.source_residual_j;
        if f0.abs() <= exchange.source_tolerance_j {
            drop(evaluate); poll(&mut cancelled)?; return Ok(first);
        }
        // Only sample the side of zero indicated by the actual midpoint source.
        // A non-monotone caller law may have other roots; none is inferred here.
        let (mut fl, mut fu);
        let edge = if f0 < 0.0 { lo = 0.0; hi } else { hi = 0.0; lo };
        let end = evaluate(edge, 2)?;
        let fe = end.step.source_residual_j;
        if fe.abs() <= exchange.source_tolerance_j {
            drop(evaluate); poll(&mut cancelled)?; return Ok(end);
        }
        if f0 < 0.0 { fl = f0; fu = fe; } else { fl = fe; fu = f0; }
        if fl > 0.0 || fu < 0.0 || lo == hi {
            return Err(ThermalPistonError::OutsideBracket { lower_residual_j: fl, upper_residual_j: fu });
        }
        let mut residual = f0;
        for iteration in 1..=exchange.maximum_iterations {
            let middle = 0.5*lo+0.5*hi;
            let fraction = -fl/(fu-fl);
            let q = if fraction.is_finite() && fraction > 0.1 && fraction < 0.9 {
                (1.0-fraction)*lo+fraction*hi
            } else { middle };
            if q <= lo || q >= hi {
                return Err(ThermalPistonError::NotConverged { iterations: iteration-1, residual_j: residual });
            }
            let candidate = evaluate(q, iteration+2)?;
            residual = candidate.step.source_residual_j;
            if residual.abs() <= exchange.source_tolerance_j {
                drop(evaluate); poll(&mut cancelled)?; return Ok(candidate);
            }
            if residual < 0.0 { lo = q; fl = residual; } else { hi = q; fu = residual; }
        }
        Err(ThermalPistonError::NotConverged { iterations: exchange.maximum_iterations, residual_j: residual })
    }
}

#[cfg(test)]
mod tests;
