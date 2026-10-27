//! Finite material heat, reaction and native elastic motion at one endpoint.
//!
//! The shared bounded heat-residual owner eliminates Q. Every Q trial restarts
//! the original coupled reaction/pressure/modal solve at the accepted state,
//! with U1=U0+Q_gas+Q-W and H1=H0+Q_body-Q. The original source conductance is
//! evaluated at both proposed endpoints. There is no chemistry-then-heat lag,
//! independent piston, heat-capacity approximation or duplicate reaction heat.

use super::coupled::{CoupledReactiveModalGasSolve, CoupledReactiveModalGasStep};
use super::{GasValveMotion, ReactiveModalGasError};
use crate::cavity_phs::piston::reactive::{GasReaction, GasReactionSolve};
use crate::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve, feedback};
use crate::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use fs_exec::Cx;
use fs_material::phase::{
    EquilibriumPhaseState, PhaseStateError, UniformEnthalpyStep, UniformEnthalpyStepInput,
};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

/// Each native solve retains its own work limits; no convergence is inferred
/// from bracket width. The heat law and gas chemistry use backward endpoints.
#[derive(Clone, Copy, Debug)]
pub struct CoupledReactiveModalThermalSolve {
    pub coupled: CoupledReactiveModalGasSolve,
    /// Q is positive from the finite body into the reacting gas.
    /// `source_tolerance_j` bounds |Q-dt G(T_body,1-T_gas,1)|.
    pub exchange: HeatExchangeSolve,
    /// Absolute allowance for each discrete thermal entropy inequality [J/K].
    /// The original chemical entropy allowance remains independently active.
    pub entropy_tolerance_j_k: f64,
}

#[derive(Debug)]
pub enum ReactiveModalThermalError<E> {
    Reactive(ReactiveModalGasError),
    Thermal(PhaseStateError),
    Interface(E),
    Invalid(&'static str),
    OutsideBracket { lower_residual_j: f64, upper_residual_j: f64 },
    NotConverged { iterations: usize, residual_j: f64 },
    Cancelled,
}
impl<E> From<ReactiveModalGasError> for ReactiveModalThermalError<E> {
    fn from(error: ReactiveModalGasError) -> Self {
        match error {
            ReactiveModalGasError::Cancelled => Self::Cancelled,
            other => Self::Reactive(other),
        }
    }
}
impl<E> From<PhaseStateError> for ReactiveModalThermalError<E> {
    fn from(error: PhaseStateError) -> Self { Self::Thermal(error) }
}
impl<E: core::fmt::Display> core::fmt::Display for ReactiveModalThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Reactive(error) => write!(f, "reacting modal heat: {error}"),
            Self::Thermal(error) => write!(f, "reacting modal thermal material: {error}"),
            Self::Interface(error) => write!(f, "reacting modal heat interface: {error}"),
            Self::Invalid(message) => write!(f, "reacting modal heat: {message}"),
            Self::OutsideBracket { lower_residual_j, upper_residual_j } => write!(f,
                "reacting modal heat not bracketed: {lower_residual_j}, {upper_residual_j} J"),
            Self::NotConverged { iterations, residual_j } => write!(f,
                "reacting modal heat did not converge in {iterations} iterations: {residual_j} J"),
            Self::Cancelled => write!(f, "reacting modal heat cancelled"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ReactiveModalThermalError<E> {}

/// One enclosing chemical, mechanical, caloric and interface account.
#[derive(Clone, Debug)]
pub struct CoupledReactiveModalThermalStep {
    pub coupled: CoupledReactiveModalGasStep,
    pub exchanged_heat_j: f64,
    pub conductance_w_k: f64,
    pub source_residual_j: f64,
    pub body_energy_change_j: f64,
    pub body_energy_residual_j: f64,
    /// dU_gas+dH_body+dE_modal+D_modal+W_exterior-Q_gas-Q_body [J].
    /// Interface Q and chemical formation energy are internal, never added again.
    pub total_energy_residual_j: f64,
    pub body_entropy_change_j_k: f64,
    /// dS_body-(Q_body-Q)/T_body,1; includes BE time-discretization entropy.
    pub body_entropy_residual_j_k: f64,
    /// Q*(1/T_gas,1-1/T_body,1), using the accepted transfer, not a second flux.
    pub interface_entropy_j_k: f64,
    /// dS_gas+dS_body-Q_gas/T_gas,1-Q_body/T_body,1 [J/K].
    /// This is a discrete inequality, not a continuous entropy-production rate.
    pub total_entropy_residual_j_k: f64,
    pub exchange_trials: usize,
    /// Native pressure/reaction endpoint trials summed across ALL heat probes.
    pub coupled_trials: usize,
}

/// Previewing does not advance either accepted owner. Publish all parts together.
#[derive(Clone, Debug)]
pub struct CoupledReactiveModalThermalCandidate<S> {
    motion: GasValveMotion,
    thermal: UniformEnthalpyStep<S>,
    step: CoupledReactiveModalThermalStep,
}
impl<S> CoupledReactiveModalThermalCandidate<S> {
    pub const fn motion(&self) -> &GasValveMotion { &self.motion }
    pub const fn thermal(&self) -> &UniformEnthalpyStep<S> { &self.thermal }
    pub const fn step(&self) -> &CoupledReactiveModalThermalStep { &self.step }
    pub fn into_parts(self) -> (GasValveMotion, UniformEnthalpyStep<S>, CoupledReactiveModalThermalStep) {
        (self.motion, self.thermal, self.step)
    }
}

fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn poll<E>(cx: &Cx<'_>) -> Result<(), ReactiveModalThermalError<E>> {
    cx.checkpoint().map_err(|_| ReactiveModalThermalError::Cancelled)
}
fn feedback_error<E>(error: feedback::ScalarError<ReactiveModalThermalError<E>>)
    -> ReactiveModalThermalError<E> {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(message) => ReactiveModalThermalError::Invalid(message),
        feedback::ScalarError::Cancelled => ReactiveModalThermalError::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } =>
            ReactiveModalThermalError::OutsideBracket { lower_residual_j: lower, upper_residual_j: upper },
        feedback::ScalarError::NotConverged { iterations, residual } =>
            ReactiveModalThermalError::NotConverged { iterations, residual_j: residual },
    }
}

impl GasValveMotion {
    /// Solve reciprocal finite-body heat, native motion and implicit chemistry.
    ///
    /// `resolve` sees the SAME final gas that determines rates, EOS and pressure
    /// work, and the actual body enthalpy/phase endpoint. It must be repeatable;
    /// no callback may debit a live reservoir. Each heat trial repeats the whole
    /// coupled sample from the same accepted history, never an already-reacted
    /// or already-moved trial. All owners advance by exactly one native sample.
    ///
    /// Both prescribed heat inputs exclude interface Q. Native structural
    /// damping leaves this enclosure and is not deposited again. The finite
    /// thermal body retains its source chart and negligible-pV-work enthalpy
    /// convention; uniform temperature and fixed geometry need justification.
    /// It is a separate thermal participant, not an implicit update of the
    /// elastic wall's source-fixed basis, thermal expansion or strength. Phase
    /// evolution here does not authorize a mechanical/fluid phase handoff.
    ///
    /// The method is first order: endpoint chemistry/conduction and native ZOH
    /// structural forcing. Native clock/admission gates remain active; refine
    /// the whole coupled step to assess temporal/splitting error. No uniqueness,
    /// flame front, shock or continuum combustion accuracy is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_thermal_candidate<S, E>(
        &self,
        reaction: &GasReaction,
        gas_heat_j: f64,
        exterior_pressure_pa: f64,
        body: UniformEnthalpyStepInput<'_>,
        chemistry: GasReactionSolve,
        controls: CoupledReactiveModalThermalSolve,
        mut resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        cx: &Cx<'_>,
    ) -> Result<CoupledReactiveModalThermalCandidate<S>, ReactiveModalThermalError<E>> {
        poll(cx)?;
        let invalid = ReactiveModalThermalError::Invalid;
        let exchange = controls.exchange;
        if !positive(body.mass_kg) || !positive(body.volume_m3)
            || !positive(body.surface_area_m2) || !positive(body.duration_s)
            || body.duration_s != self.motion().model().sample_period_s()
            || !body.internal_heat_j.is_finite() || !gas_heat_j.is_finite()
            || !positive(exchange.energy_tolerance_j)
            || !controls.entropy_tolerance_j_k.is_finite()
            || controls.entropy_tolerance_j_k < 0.0
        {
            return Err(invalid("finite heat, positive body geometry/mass, one clock and bounded energy/entropy required"));
        }
        let h0 = body.initial.specific_enthalpy_j_kg();
        if body.curve.state_at_specific_enthalpy(h0)? != body.initial {
            return Err(invalid("finite body initial state belongs to a different material chart"));
        }
        let mut coupled_trials = 0usize;
        let candidate = feedback::solve_scalar(
            exchange.heat_bracket_j, exchange.source_tolerance_j, exchange.maximum_iterations,
            |q, exchange_trials, _| {
                let body_heat = body.internal_heat_j - q;
                let gas_heat = gas_heat_j + q;
                let dh = body_heat / body.mass_kg;
                let h1 = h0 + dh;
                if !body_heat.is_finite() || !gas_heat.is_finite() || !dh.is_finite()
                    || !h1.is_finite() || (body_heat != 0.0 && (dh == 0.0 || h1 == h0))
                    || (q != 0.0 && (body_heat == body.internal_heat_j || gas_heat == gas_heat_j))
                {
                    return Err(invalid("reciprocal heat disappears or overflows in physical energy coordinates"));
                }
                let endpoint = body.curve.state_at_specific_enthalpy(h1)?;
                let (motion, coupled) = self.coupled_reactive_step_candidate(
                    reaction, gas_heat, exterior_pressure_pa, chemistry, controls.coupled, cx,
                )?.into_parts();
                coupled_trials = coupled_trials.checked_add(coupled.trials)
                    .ok_or_else(|| invalid("coupled heat trial count overflow"))?;
                poll(cx)?;
                let sample = resolve(motion.gas().state(), endpoint)
                    .map_err(ReactiveModalThermalError::Interface)?;
                let g = sample.conductance_w_k;
                let tg = motion.gas().state().state().temperature;
                let tb = endpoint.temperature_k();
                let delta = tb - tg;
                let power = g * delta;
                let heat = body.duration_s * power;
                let residual = q - heat;
                if !g.is_finite() || g < 0.0 || !power.is_finite() || !heat.is_finite()
                    || !residual.is_finite() || (g != 0.0 && delta != 0.0 && power == 0.0)
                    || (power != 0.0 && heat == 0.0)
                {
                    return Err(invalid("source conductance or endpoint heat flux is not nonnegative and representable"));
                }
                let body_energy = body.mass_kg * (endpoint.specific_enthalpy_j_kg() - h0);
                let body_residual = body_energy + q - body.internal_heat_j;
                let total = coupled.reaction.reaction.stored_energy_change_j
                    + body_energy + coupled.modal_energy_change_j
                    + coupled.motion.modal.viscous_dissipation_j
                    + coupled.exterior_pressure_work_j - gas_heat_j - body.internal_heat_j;
                if [body_energy, body_residual, total].iter().any(|x| !x.is_finite())
                    || body_residual.abs() > exchange.energy_tolerance_j
                    || total.abs() > exchange.energy_tolerance_j
                {
                    return Err(invalid("finite body, reaction and native modal energy do not close"));
                }
                let body_entropy = body.mass_kg * body.curve.specific_entropy_change_j_kg_k(h0, h1)?;
                let body_entropy_residual = body_entropy - body_heat / tb;
                let interface_entropy = q * ((tb - tg) / tb) / tg;
                let total_entropy = coupled.reaction.reaction.entropy_change_j_k
                    + body_entropy - gas_heat_j / tg - body.internal_heat_j / tb;
                if [body_entropy, body_entropy_residual, interface_entropy, total_entropy]
                    .iter().any(|x| !x.is_finite())
                    || body_entropy_residual < -controls.entropy_tolerance_j_k
                {
                    return Err(invalid("finite body caloric entropy is unrepresentable or decreases below supplied heat entropy"));
                }
                // A trial Q need not yet have the conductive direction. Only
                // the converged endpoint may be admitted by interface entropy.
                let step = CoupledReactiveModalThermalStep {
                    coupled, exchanged_heat_j: q, conductance_w_k: g,
                    source_residual_j: residual, body_energy_change_j: body_energy,
                    body_energy_residual_j: body_residual, total_energy_residual_j: total,
                    body_entropy_change_j_k: body_entropy,
                    body_entropy_residual_j_k: body_entropy_residual,
                    interface_entropy_j_k: interface_entropy,
                    total_entropy_residual_j_k: total_entropy,
                    exchange_trials, coupled_trials,
                };
                Ok((CoupledReactiveModalThermalCandidate {
                    motion,
                    thermal: UniformEnthalpyStep { state: endpoint, external_heat_j: -q,
                        energy_residual_tolerance_j: exchange.energy_tolerance_j,
                        report: sample.source },
                    step,
                }, residual))
            },
            &mut || cx.checkpoint().is_err(),
        ).map_err(feedback_error)?;
        if candidate.step.interface_entropy_j_k < -controls.entropy_tolerance_j_k
            || candidate.step.total_entropy_residual_j_k < -controls.entropy_tolerance_j_k
        {
            return Err(invalid("accepted material heat transfer or enclosing entropy inequality failed"));
        }
        poll(cx)?;
        Ok(candidate)
    }

    /// Use the original source-bound gas/body series conduction law at the
    /// SAME endpoints. Its body conductivity must belong to the enthalpy card;
    /// all source axes and the original reacted gas transport remain intact.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_thermal_layers_candidate(
        &self,
        reaction: &GasReaction,
        gas_heat_j: f64,
        exterior_pressure_pa: f64,
        body: UniformEnthalpyStepInput<'_>,
        chemistry: GasReactionSolve,
        controls: CoupledReactiveModalThermalSolve,
        layers: MaterialConductionLayers<'_>,
        cx: &Cx<'_>,
    ) -> Result<CoupledReactiveModalThermalCandidate<ResolvedMaterialStatePoint>,
        ReactiveModalThermalError<LayerConductionError>> {
        self.coupled_reactive_thermal_candidate(
            reaction, gas_heat_j, exterior_pressure_pa, body, chemistry, controls,
            |gas, thermal| layers.resolve(gas, thermal), cx,
        )
    }
}
