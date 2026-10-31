//! Reciprocal finite-body heat exchange on the existing gas/piston step.
//!
//! The solid/liquid material owns h -> (T, phase); the gas owns u, EOS and
//! acoustic/transport response; `GasPiston::step` owns motion and pressure work.
//! One signed Q goes into the gas and OUT of the thermal body. Source laws are
//! evaluated on their actual energy midpoints, never averaged temperatures.
//! The thermal body uses the existing negligible-pV-work enthalpy storage model.
//! Its volume is fixed during this proposal; expansion/deformation is not inferred.

/// One finite enthalpy/phase body feeding both faces of the shared gas piston.
pub mod double_acting;
/// Source-bound gas/body conductivity without copied material coefficients.
pub mod layers;
/// State-dependent power on the same material, phase and mechanical step.
pub mod power;
/// Shared material/phase coordinates for thermal and electrical responses.
pub mod source;
pub(crate) mod feedback;
/// A thermal port for reservoirs and jointly solved finite network boundaries.
pub mod boundary;

use super::{GasPiston, PistonError, PistonGas, PistonGasState, PistonInput, PistonSolve, PistonStep};
use fs_material::phase::{EquilibriumPhaseState, PhaseStateError,
    UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// Whole-step source/exchange bracket and absolute residual/work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatExchangeSolve {
    /// Signed joules into the gas for an exchange, or into the body for the
    /// extra source in `powered_thermal_candidate`; the bracket contains zero.
    /// Every queried proposal must be supported by both materials and motion.
    pub heat_bracket_j: [f64; 2],
    /// Maximum |Q - dt G(T_body_mid - T_gas_mid)|, or |W - dt P| [J].
    pub source_tolerance_j: f64,
    /// Maximum absolute body and combined energy residual [J].
    pub energy_tolerance_j: f64,
    /// Interior safeguarded secant/bisection trials, in 1..=128. At most two
    /// extra trials check zero and one bracket endpoint. Every nested solve
    /// retains its separate explicit work limits.
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
    /// Caller-owned interface or power-source law refused the midpoint.
    Interface(E),
    /// Cancellation before returning an accepted joint candidate.
    Cancelled,
    /// No sign change in the caller's admitted source/exchange bracket.
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
            Self::Interface(e) => write!(f, "gas/body thermal interface/source: {e}"),
            Self::Cancelled => write!(f, "gas/body thermal coupling cancelled"),
            Self::OutsideBracket { lower_residual_j, upper_residual_j } =>
                write!(f, "thermal feedback not bracketed: {lower_residual_j}, {upper_residual_j} J"),
            Self::NotConverged { iterations, residual_j } =>
                write!(f, "thermal feedback did not converge in {iterations} trials: {residual_j} J"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ThermalPistonError<E> {}

/// Actual accepted balances; no independent entropy-production claim is implied.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalPistonStep<S = ResolvedThermochemicalGasState> {
    /// Original coupled gas/momentum result, including the gas midpoint sources.
    pub piston: PistonStep<S>,
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
    /// Gas material proposals summed across ALL exchange trials, not only the last.
    pub material_trials: usize,
}

/// Both accepted physical states in one candidate. Neither input is mutated.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalPistonCandidate<S, G: PistonGas = ClosedGasParcel> {
    piston: GasPiston<G>,
    thermal: UniformEnthalpyStep<S>,
    step: ThermalPistonStep<G::State>,
}
impl<S, G: PistonGas> ThermalPistonCandidate<S, G> {
    /// Tentative gas, volume and momentum; old sources/history remain intact.
    pub const fn piston(&self) -> &GasPiston<G> { &self.piston }
    /// Canonical thermal proposal with original concrete interface source.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> { &self.thermal }
    /// Independently checked coupled balances.
    pub const fn step(&self) -> &ThermalPistonStep<G::State> { &self.step }
    /// Publish both states together in the caller's existing physical owner.
    pub fn into_parts(self) -> (GasPiston<G>, UniformEnthalpyStep<S>, ThermalPistonStep<G::State>) {
        (self.piston, self.thermal, self.step)
    }
}

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), ThermalPistonError<E>> {
    if cancelled() { Err(ThermalPistonError::Cancelled) } else { Ok(()) }
}

impl<G: PistonGas> GasPiston<G> {
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
    /// Ideal and real-gas owners use this same reciprocal residual; the actual
    /// source-bearing midpoint carries its own pressure, density and transport.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_candidate<S, E>(
        &self, input: PistonInput, body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve, exchange: HeatExchangeSolve,
        mut resolve: impl FnMut(&G::State, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThermalPistonCandidate<S, G>, ThermalPistonError<E>> {
        poll(&mut cancelled)?;
        feedback::admit(exchange)?;
        if !positive(body.mass_kg) || !positive(body.volume_m3)
            || !positive(body.surface_area_m2) || !positive(body.duration_s)
            || body.duration_s != input.duration_s || !body.internal_heat_j.is_finite()
            || !input.heat_j.is_finite()
        { return Err(ThermalPistonError::Invalid("invalid heat bounds, clocks, body geometry/mass or controls")); }
        let h0 = body.initial.specific_enthalpy_j_kg();
        if body.curve.state_at_specific_enthalpy(h0)? != body.initial {
            return Err(ThermalPistonError::Invalid("initial thermal state belongs to a different source chart"));
        }
        let mut material_trials = 0usize;
        feedback::solve(exchange, |q, trials, cancelled| {
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
            let piston = next.step(PistonInput { heat_j: gas_heat, ..input }, controls, &mut *cancelled)?;
            material_trials += piston.material_trials;
            poll(cancelled)?;
            let sample = resolve(&piston.midpoint, thermal_midpoint).map_err(ThermalPistonError::Interface)?;
            let g = sample.conductance_w_k;
            let delta = thermal_midpoint.temperature_k() - piston.midpoint.temperature_k();
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
            Ok((ThermalPistonCandidate {
                piston: next,
                thermal: UniformEnthalpyStep { state: endpoint, external_heat_j: -q,
                    energy_residual_tolerance_j: exchange.energy_tolerance_j, report: sample.source },
                step: ThermalPistonStep { piston, thermal_midpoint, conductance_w_k: g,
                    exchanged_heat_j: q, source_residual_j: residual,
                    thermal_energy_change_j: stored, thermal_energy_residual_j: body_residual,
                    total_energy_residual_j: total, exchange_trials: trials, material_trials },
            }, residual))
        }, &mut cancelled)
    }
}

#[cfg(test)]
mod tests;
