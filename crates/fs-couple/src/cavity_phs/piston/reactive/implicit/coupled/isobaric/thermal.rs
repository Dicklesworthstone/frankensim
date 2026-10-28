//! Endpoint heat feedback for constant-pressure homogeneous chemistry.
//!
//! One signed exchange Q drives H1=H0+Q_external+Q. Every proposal starts at
//! the same accepted gas; repeated heat trials never repeat physical burning.
//! The original implicit mechanism, material closure and scalar heat solver
//! remain the numerical owners. Contact heat is evaluated at the final gas.

use super::{GasReaction, GasReactionError, GasReactionSolve, ClosedGasParcel,
    IsobaricGasReactionStep};
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonError, feedback,
};
use crate::cavity_phs::piston::thermal::boundary::heat::conductive_heat;
use crate::cavity_phs::piston::thermal::boundary::ports::{ThermalPort, ThermalPortSample};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use std::collections::BTreeMap;

/// Heat-root limits plus a strict shared cap across ALL nested chemistry work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricHeatSolve {
    pub exchange: HeatExchangeSolve,
    /// Endpoint/Jacobian evaluations across accepted and discarded heat trials.
    /// Each trial also retains the original mechanism's smaller per-solve cap.
    /// Must lie in 1..=1_048_576. A cap is checked before the next evaluation.
    pub maximum_chemical_evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum IsobaricThermalError<E> {
    Invalid(&'static str),
    Reaction(GasReactionError),
    Thermal(ThermalPistonError<E>),
    ChemicalEvaluationLimit { evaluations: usize, maximum: usize },
    Cancelled,
}
impl<E> From<GasReactionError> for IsobaricThermalError<E> {
    fn from(error: GasReactionError) -> Self {
        match error {
            GasReactionError::Cancelled => Self::Cancelled,
            other => Self::Reaction(other),
        }
    }
}
impl<E> From<ThermalPistonError<E>> for IsobaricThermalError<E> {
    fn from(error: ThermalPistonError<E>) -> Self {
        match error {
            ThermalPistonError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for IsobaricThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "constant-pressure reacting heat exchange: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for IsobaricThermalError<E> {}

/// One chemical endpoint and its independently evaluated thermal boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricBoundaryStep {
    pub chemistry: IsobaricGasReactionStep,
    /// Signed heat INTO gas, already present in the candidate's enthalpy.
    pub exchanged_heat_j: f64,
    /// Sum of actual signed endpoint contact laws [J].
    pub evaluated_heat_j: f64,
    pub source_residual_j: f64,
    /// Delta H_gas - prescribed gas heat - evaluated contact heat [J].
    pub enthalpy_residual_j: f64,
    /// Delta U_gas + actual p Delta V - prescribed/evaluated heat [J].
    pub energy_residual_j: f64,
    /// Sum Q_law*(1/T_gas-1/T_contact) [J/K]. This is the endpoint contact
    /// law's entropy, not the entropy of an unspecified external reservoir.
    pub contact_entropy_j_k: f64,
    pub exchange_trials: usize,
    pub chemical_evaluations: usize,
}

/// Publish the gas and every finite counterpart together; never deposit Q again.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricBoundaryCandidate<S> {
    gas: ClosedGasParcel,
    source: S,
    step: IsobaricBoundaryStep,
}
impl<S> IsobaricBoundaryCandidate<S> {
    pub const fn gas(&self) -> &ClosedGasParcel { &self.gas }
    pub const fn source(&self) -> &S { &self.source }
    pub const fn step(&self) -> &IsobaricBoundaryStep { &self.step }
    pub fn into_parts(self) -> (ClosedGasParcel, S, IsobaricBoundaryStep) {
        (self.gas, self.source, self.step)
    }
}

fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), IsobaricThermalError<E>> {
    if cancelled() { Err(IsobaricThermalError::Cancelled) } else { Ok(()) }
}
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }

// Cancel reciprocal-temperature subtraction before rounding. Called only with
// an independently evaluated passive conduction law, never with guessed heat.
fn contact_entropy<E>(heat: f64, boundary: f64, gas: f64) -> Result<f64, IsobaricThermalError<E>> {
    let entropy = (heat / gas) * ((boundary - gas) / boundary);
    if !entropy.is_finite() || entropy < 0.0 || (heat != 0.0 && entropy == 0.0) {
        return Err(IsobaricThermalError::Invalid("contact entropy is unrepresentable"));
    }
    Ok(entropy)
}

struct HeatSample<S> {
    heat_j: f64,
    entropy_j_k: f64,
    source: S,
}

impl GasReaction {
    /// Close original endpoint chemistry against any finite list of thermal ports.
    ///
    /// Port temperatures belong to the enclosing thermal trial. All conductance
    /// callbacks see the SAME reacted endpoint; no frozen initial mixture or
    /// temperature is substituted. IDs establish deterministic canonical order.
    /// Opposing transfers remain separate even when their net gas heat is zero.
    /// A finite network must debit each returned contact exactly once and commit
    /// all histories together. Callbacks must be repeatable and side-effect free.
    ///
    /// This is first-order closed ideal-gas chemistry under a quasistatic
    /// pressure load, not an open-flow reactor or a whole-network passivity proof.
    #[allow(clippy::too_many_arguments)]
    pub fn isobaric_thermal_ports_candidate<S, E>(
        &self, gas: &ClosedGasParcel, duration_s: f64, external_heat_j: f64,
        chemistry: GasReactionSolve, controls: IsobaricHeatSolve,
        ports: &[ThermalPort], maximum_ports: usize,
        mut resolve: impl FnMut(ThermalPort, &ResolvedThermochemicalGasState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricBoundaryCandidate<Vec<ThermalPortSample<S>>>, IsobaricThermalError<E>> {
        poll(&mut cancelled)?;
        if ports.len() > maximum_ports {
            return Err(IsobaricThermalError::Invalid("isobaric thermal port budget exceeded"));
        }
        let mut ordered = BTreeMap::new();
        for port in ports {
            poll(&mut cancelled)?;
            if !positive(port.temperature_k) || ordered.insert(port.id, *port).is_some() {
                return Err(IsobaricThermalError::Invalid("unique port IDs and positive temperatures required"));
            }
        }
        solve(self, gas, duration_s, external_heat_j, chemistry, controls,
            |_, state, cancel| {
                let mut heat = Sum::default();
                let mut entropy = Sum::default();
                let mut samples = Vec::with_capacity(ordered.len());
                for &port in ordered.values() {
                    poll(cancel)?;
                    let sample = resolve(port, state).map_err(ThermalPistonError::Interface)?;
                    let heat_j = conductive_heat(duration_s, sample.conductance_w_k,
                        port.temperature_k, state.state().temperature)?;
                    heat.add(heat_j)?;
                    entropy.add(contact_entropy::<E>(heat_j, port.temperature_k, state.state().temperature)?)?;
                    samples.push(ThermalPortSample { port, conductance_w_k: sample.conductance_w_k,
                        heat_j, source: sample.source });
                }
                Ok(HeatSample { heat_j: heat.value, entropy_j_k: entropy.value, source: samples })
            }, &mut cancelled)
    }
}

// Fixed canonical reduction; keep small opposing contact contributions.
#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add<E>(&mut self, value: f64) -> Result<(), IsobaricThermalError<E>> {
        let corrected = value - self.correction;
        let next = self.value + corrected;
        self.correction = (next - self.value) - corrected;
        self.value = next;
        if !next.is_finite() || !self.correction.is_finite() {
            return Err(IsobaricThermalError::Invalid("summed thermal exchange is unrepresentable"));
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn solve<S, E, C: FnMut() -> bool>(
    law: &GasReaction, gas: &ClosedGasParcel, duration_s: f64, external_heat_j: f64,
    chemistry: GasReactionSolve, controls: IsobaricHeatSolve,
    mut resolve: impl FnMut(f64, &ResolvedThermochemicalGasState, &mut C)
        -> Result<HeatSample<S>, IsobaricThermalError<E>>,
    cancelled: &mut C,
) -> Result<IsobaricBoundaryCandidate<S>, IsobaricThermalError<E>> {
    use IsobaricThermalError::{Invalid, ChemicalEvaluationLimit};
    feedback::admit(controls.exchange)?;
    if !(1..=1_048_576).contains(&controls.maximum_chemical_evaluations) {
        return Err(Invalid("isobaric chemistry evaluation cap must be in 1..=1048576"));
    }
    let original = law.isobaric_policy(duration_s, external_heat_j, chemistry, cancelled)?;
    // Clone only the bounded numerical policy ONCE, not a mechanism/catalog per
    // heat trial. The original owner enforces the exact remaining shared budget.
    let mut policy = original.clone();
    let per_trial = policy.maximum_evaluations;
    let mut spent = 0usize;
    let exchange = controls.exchange;
    let result = feedback::solve_scalar(exchange.heat_bracket_j, exchange.source_tolerance_j,
        exchange.maximum_iterations, |q, trials, cancel| {
            let remaining = controls.maximum_chemical_evaluations - spent;
            if remaining == 0 {
                return Err(ChemicalEvaluationLimit { evaluations: spent,
                    maximum: controls.maximum_chemical_evaluations });
            }
            let heat = external_heat_j + q;
            if !heat.is_finite() || (q != 0.0 && heat == external_heat_j)
                || (external_heat_j != 0.0 && heat == q)
            {
                return Err(Invalid("prescribed or exchanged heat vanished or overflowed"));
            }
            policy.maximum_evaluations = per_trial.min(remaining);
            let proposal = law.isobaric_candidate_with_policy(gas, duration_s, heat,
                chemistry, &policy, cancel).map_err(|error| match error {
                    GasReactionError::ImplicitEvaluationLimit { evaluations, .. }
                        if remaining <= per_trial => ChemicalEvaluationLimit {
                            evaluations: spent + evaluations,
                            maximum: controls.maximum_chemical_evaluations,
                        },
                    other => other.into(),
                })?;
            let used = proposal.step().boundary.reaction.implicit_mechanism.as_ref()
                .ok_or(Invalid("isobaric endpoint lost its implicit mechanism report"))?.evaluations;
            spent = spent.checked_add(used).ok_or(Invalid("chemical work counter overflowed"))?;
            if spent > controls.maximum_chemical_evaluations {
                return Err(ChemicalEvaluationLimit { evaluations: spent,
                    maximum: controls.maximum_chemical_evaluations });
            }
            poll(cancel)?;
            let sample = resolve(q, proposal.gas().state(), cancel)?;
            let residual = q - sample.heat_j;
            let (gas, chemistry) = proposal.into_parts();
            let enthalpy_residual = (chemistry.enthalpy_change_j - external_heat_j) - sample.heat_j;
            let energy_residual = (chemistry.boundary.reaction.stored_energy_change_j
                + chemistry.boundary.boundary.work_by_gas_j - external_heat_j) - sample.heat_j;
            if [sample.heat_j, sample.entropy_j_k, residual, enthalpy_residual, energy_residual]
                .iter().any(|value| !value.is_finite()) || sample.entropy_j_k < 0.0
            {
                return Err(Invalid("isobaric boundary heat, work or entropy is unrepresentable"));
            }
            Ok((IsobaricBoundaryCandidate { gas, source: sample.source,
                step: IsobaricBoundaryStep { chemistry, exchanged_heat_j: q,
                    evaluated_heat_j: sample.heat_j, source_residual_j: residual,
                    enthalpy_residual_j: enthalpy_residual, energy_residual_j: energy_residual,
                    contact_entropy_j_k: sample.entropy_j_k, exchange_trials: trials,
                    chemical_evaluations: spent } }, residual))
        }, cancelled).map_err(|error| match error {
            feedback::ScalarError::Evaluation(error) => error,
            feedback::ScalarError::Cancelled => IsobaricThermalError::Cancelled,
            feedback::ScalarError::Invalid(message) => Invalid(message),
            feedback::ScalarError::OutsideBracket { lower, upper } =>
                IsobaricThermalError::Thermal(ThermalPistonError::OutsideBracket {
                    lower_residual_j: lower, upper_residual_j: upper }),
            feedback::ScalarError::NotConverged { iterations, residual } =>
                IsobaricThermalError::Thermal(ThermalPistonError::NotConverged {
                    iterations, residual_j: residual }),
        })?;
    if result.step.enthalpy_residual_j.abs() > exchange.energy_tolerance_j
        || result.step.energy_residual_j.abs() > exchange.energy_tolerance_j
    {
        return Err(Invalid("evaluated boundary heat did not close enthalpy and pressure work"));
    }
    poll(cancelled)?;
    Ok(result)
}
