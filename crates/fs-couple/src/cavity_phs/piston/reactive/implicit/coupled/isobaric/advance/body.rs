//! Adaptive constant-pressure chemistry and reciprocal finite-body heat.
//!
//! Every physical substep is the existing simultaneous chemical/thermal solve.
//! Gas and body remain one candidate; their actual source reports are retained
//! once, and discarded meshes never contribute heat, storage or entropy.

use super::*;
use super::super::thermal::{IsobaricThermalError,
    body::{IsobaricBodySolve, IsobaricBodyStep}};
use crate::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use crate::cavity_phs::piston::thermal::layers::{LayerConductionError, MaterialConductionLayers};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricBodyAdvanceSolve {
    /// Shared species/extents/gas-volume controls and total work limits.
    /// Its temperature tolerance applies to BOTH physical temperatures, and
    /// endpoint evaluations count attempted complete gas/body intervals.
    pub trajectory: IsobaricAdvanceSolve,
    /// Body specific enthalpy resolves latent plateaus even at unchanged T.
    /// Its relative scale uses h-h_initial, so a datum does not loosen accuracy.
    pub absolute_body_enthalpy_j_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum IsobaricBodyAdvanceError<E> {
    Trajectory(IsobaricAdvanceError),
    Thermal(IsobaricThermalError<E>),
}
impl<E> From<IsobaricAdvanceError> for IsobaricBodyAdvanceError<E> {
    fn from(error: IsobaricAdvanceError) -> Self { Self::Trajectory(error) }
}
impl<E> From<GasReactionError> for IsobaricBodyAdvanceError<E> {
    fn from(error: GasReactionError) -> Self { Self::Trajectory(error.into()) }
}
impl<E> From<IsobaricThermalError<E>> for IsobaricBodyAdvanceError<E> {
    fn from(error: IsobaricThermalError<E>) -> Self {
        match error {
            IsobaricThermalError::Cancelled => Self::Trajectory(IsobaricAdvanceError::Cancelled),
            other => Self::Thermal(other),
        }
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for IsobaricBodyAdvanceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive isobaric gas/body exchange: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for IsobaricBodyAdvanceError<E> {}

#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricBodyAdvanceSubstep<S> {
    pub start_s: f64,
    pub end_s: f64,
    pub prescribed_gas_heat_j: f64,
    pub prescribed_body_heat_j: f64,
    /// Actual accepted gas endpoint, for its prescribed heat-entropy account.
    pub gas_temperature_k: f64,
    /// Original continuing phase state and concrete interface source.
    pub thermal: UniformEnthalpyStep<S>,
    /// Complete same-endpoint chemistry, reciprocal heat and entropy report.
    pub step: IsobaricBodyStep,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricBodyAdvanceReport<S> {
    /// Original gas summary and the complete joint accepted history, retained
    /// once. Its supplied_heat_j includes prescribed PLUS reciprocal gas heat.
    pub trajectory: IsobaricAdvanceReport<IsobaricBodyAdvanceSubstep<S>>,
    pub prescribed_gas_heat_j: f64,
    pub prescribed_body_heat_j: f64,
    /// Actual Q into gas; -Q is already deposited into the returned body.
    pub exchanged_heat_j: f64,
    pub evaluated_contact_heat_j: f64,
    pub source_residual_j: f64,
    pub thermal_energy_change_j: f64,
    pub thermal_energy_residual_j: f64,
    pub total_enthalpy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub thermal_entropy_change_j_k: f64,
    pub thermal_entropy_residual_j_k: f64,
    pub applied_contact_entropy_j_k: f64,
    pub evaluated_contact_entropy_j_k: f64,
    /// Sum of BOTH prescribed heats divided by their actual local endpoint T.
    pub prescribed_heat_entropy_j_k: f64,
    pub total_entropy_residual_j_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricBodyAdvanceCandidate<S> {
    gas: ClosedGasParcel,
    thermal: UniformEnthalpyStep<IsobaricBodyAdvanceReport<S>>,
}
impl<S> IsobaricBodyAdvanceCandidate<S> {
    pub const fn gas(&self) -> &ClosedGasParcel { &self.gas }
    pub const fn thermal(&self) -> &UniformEnthalpyStep<IsobaricBodyAdvanceReport<S>> { &self.thermal }
    pub const fn report(&self) -> &IsobaricBodyAdvanceReport<S> { &self.thermal.report }
    pub fn into_parts(self) -> (ClosedGasParcel, UniformEnthalpyStep<IsobaricBodyAdvanceReport<S>>) {
        (self.gas, self.thermal)
    }
}

impl GasReaction {
    /// Advance simultaneous chemistry and finite-body exchange under one fixed
    /// pressure with error control over both temperatures and body enthalpy.
    ///
    /// Prescribed gas/body heats are uniform powers over this interval. Each
    /// original source heat bracket is scaled by subinterval duration, retaining
    /// the caller's signed exchange-power bounds. Source tolerances divide over
    /// accepted steps; whole-window energy and entropy allowances never grow.
    /// The same original phase chart handles sensible/latent traversal without
    /// effective heat capacity or a second latent-energy source.
    ///
    /// All trial chemistry consumes the global work budget, including failed
    /// thermal roots. Callbacks must be repeatable and must not alter accepted
    /// history. Publish gas and returned thermal state together; its reciprocal
    /// heat is already included. This retains the original uniform fixed-body,
    /// negligible-body-pV and homogeneous ideal-gas assumptions.
    #[allow(clippy::too_many_arguments)]
    pub fn isobaric_thermal_body_adaptive_candidate<S, E>(
        &self, gas: &ClosedGasParcel, duration_s: f64, external_heat_j: f64,
        body: UniformEnthalpyStepInput<'_>, chemistry: GasReactionSolve,
        controls: IsobaricBodySolve, accuracy: IsobaricBodyAdvanceSolve,
        mut resolve: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricBodyAdvanceCandidate<S>, IsobaricBodyAdvanceError<E>> {
        poll(&mut cancelled)?;
        accuracy.trajectory.validate()?;
        if !positive(accuracy.absolute_body_enthalpy_j_kg)
            || body.duration_s != duration_s || !body.internal_heat_j.is_finite()
        { return Err(invalid("positive body enthalpy accuracy and matching finite heat interval required").into()); }
        self.isobaric_policy(duration_s, external_heat_j, chemistry, &mut cancelled)?;
        let policy = accuracy.trajectory;
        let (mut count, mut meshes, mut refused) = (1, 0, 0);
        let mut work = Work::default();
        let mut coarse: Option<IsobaricBodyAdvanceCandidate<S>> = None;
        let mut last_ratio = None;
        loop {
            poll(&mut cancelled)?;
            meshes += 1;
            match march(self, gas, duration_s, external_heat_j, body, chemistry, controls,
                accuracy, count, &mut work, &mut resolve, &mut cancelled)
            {
                Ok(mut fine) => {
                    if let Some(previous) = coarse.as_ref() {
                        let ratio = compare(previous, &fine, body.initial, accuracy)?;
                        last_ratio = Some(ratio);
                        if ratio <= 1.0 {
                            let trajectory = &mut fine.thermal.report.trajectory;
                            trajectory.endpoint_evaluations = work.endpoints;
                            trajectory.chemical_evaluations = work.chemical;
                            trajectory.meshes_attempted = meshes;
                            trajectory.refused_meshes = refused;
                            trajectory.estimated_error_ratio = ratio;
                            poll(&mut cancelled)?;
                            return Ok(fine);
                        }
                    }
                    fine.thermal.report.trajectory.steps = Vec::new();
                    coarse = Some(fine);
                }
                Err(error) if retryable(&error) => {
                    coarse = None;
                    last_ratio = None;
                    refused += 1;
                }
                Err(error) => return Err(error),
            }
            if count == policy.maximum_substeps {
                return Err(IsobaricAdvanceError::RefinementLimit {
                    substeps: count, endpoint_evaluations: work.endpoints,
                    chemical_evaluations: work.chemical, last_error_ratio: last_ratio,
                }.into());
            }
            count *= 2;
        }
    }

    /// The original material/card-bound series-conduction adapter on every
    /// actual reacted gas and phase endpoint of the accepted trajectory.
    #[allow(clippy::too_many_arguments)]
    pub fn isobaric_thermal_layers_adaptive_candidate(
        &self, gas: &ClosedGasParcel, duration_s: f64, external_heat_j: f64,
        body: UniformEnthalpyStepInput<'_>, chemistry: GasReactionSolve,
        controls: IsobaricBodySolve, accuracy: IsobaricBodyAdvanceSolve,
        layers: MaterialConductionLayers<'_>, cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricBodyAdvanceCandidate<ResolvedMaterialStatePoint>,
        IsobaricBodyAdvanceError<LayerConductionError>>
    {
        if layers.area_m2 > body.surface_area_m2 {
            return Err(invalid("conducting patch exceeds body exposed area").into());
        }
        self.isobaric_thermal_body_adaptive_candidate(gas, duration_s, external_heat_j,
            body, chemistry, controls, accuracy, |gas, state| layers.resolve(gas, state), cancelled)
    }
}

fn retryable<E>(error: &IsobaricBodyAdvanceError<E>) -> bool {
    matches!(error,
        IsobaricBodyAdvanceError::Thermal(IsobaricThermalError::Reaction(
            GasReactionError::ConsumptionLimit { .. } | GasReactionError::EntropyDecrease { .. }
            | GasReactionError::ImplicitOutsideBracket { .. } | GasReactionError::ImplicitNotConverged { .. }))
        | IsobaricBodyAdvanceError::Trajectory(IsobaricAdvanceError::Reaction(
            GasReactionError::EntropyDecrease { .. }))
        | IsobaricBodyAdvanceError::Thermal(IsobaricThermalError::Thermal(
            ThermalPistonError::OutsideBracket { .. } | ThermalPistonError::NotConverged { .. })))
}

#[allow(clippy::too_many_arguments)]
fn march<S, E>(
    law: &GasReaction, initial: &ClosedGasParcel, duration: f64, external_heat: f64,
    body: UniformEnthalpyStepInput<'_>, chemistry: GasReactionSolve,
    controls: IsobaricBodySolve, accuracy: IsobaricBodyAdvanceSolve, count: usize,
    work: &mut Work,
    resolve: &mut impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
        -> Result<ConductanceSample<S>, E>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<IsobaricBodyAdvanceCandidate<S>, IsobaricBodyAdvanceError<E>> {
    let mut gas = initial.clone();
    let mut phase = body.initial;
    let mut steps = Vec::with_capacity(count);
    let (mut start, mut gas_given, mut body_given) = (0.0, 0.0, 0.0);
    let mut inner = controls;
    inner.heat.exchange.source_tolerance_j /= count as f64;
    inner.heat.exchange.heat_bracket_j = controls.heat.exchange.heat_bracket_j.map(|q| q / count as f64);
    let policy = accuracy.trajectory;
    for index in 0..count {
        poll(cancelled)?;
        if work.endpoints == policy.maximum_endpoint_evaluations {
            return Err(IsobaricAdvanceError::EndpointEvaluationLimit {
                evaluations: work.endpoints, maximum: policy.maximum_endpoint_evaluations,
            }.into());
        }
        let remaining = policy.maximum_chemical_evaluations - work.chemical;
        if remaining == 0 { return Err(chemical_limit(work, policy).into()); }
        let fraction = (index + 1) as f64 / count as f64;
        let (end, gas_end, body_end) = (duration * fraction, external_heat * fraction,
            body.internal_heat_j * fraction);
        let (dt, gas_heat, body_heat) = (end - start, gas_end - gas_given, body_end - body_given);
        if !positive(dt) || !gas_heat.is_finite() || !body_heat.is_finite()
            || (external_heat != 0.0 && gas_heat == 0.0)
            || (body.internal_heat_j != 0.0 && body_heat == 0.0)
        { return Err(invalid("gas/body subinterval or prescribed heat is unrepresentable").into()); }
        inner.heat.maximum_chemical_evaluations = controls.heat.maximum_chemical_evaluations.min(remaining);
        let spent = Cell::new(0);
        work.endpoints += 1;
        let proposal = law.isobaric_thermal_body_counted(&gas, dt, gas_heat,
            UniformEnthalpyStepInput { initial: phase, duration_s: dt,
                internal_heat_j: body_heat, ..body }, chemistry, inner,
            &mut *resolve, &mut *cancelled, &spent);
        work.chemical += spent.get();
        let proposal = proposal.map_err(|error| match error {
            IsobaricThermalError::ChemicalEvaluationLimit { .. }
                if remaining <= controls.heat.maximum_chemical_evaluations =>
                IsobaricBodyAdvanceError::Trajectory(chemical_limit(work, policy)),
            other => other.into(),
        })?;
        let (next, thermal, step) = proposal.into_parts();
        gas = next;
        phase = thermal.state;
        steps.push(IsobaricBodyAdvanceSubstep { start_s: start, end_s: end,
            prescribed_gas_heat_j: gas_heat, prescribed_body_heat_j: body_heat,
            gas_temperature_k: gas.state().state().temperature, thermal, step });
        (start, gas_given, body_given) = (end, gas_end, body_end);
    }
    summarize(law, initial, gas, duration, external_heat, body, chemistry, controls, steps, cancelled)
}

fn chemical_limit(work: &Work, policy: IsobaricAdvanceSolve) -> IsobaricAdvanceError {
    IsobaricAdvanceError::ChemicalEvaluationLimit { evaluations: work.chemical,
        maximum: policy.maximum_chemical_evaluations, endpoint_evaluations: work.endpoints }
}

#[allow(clippy::too_many_arguments)]
fn summarize<S, E>(
    law: &GasReaction, initial: &ClosedGasParcel, gas: ClosedGasParcel, duration: f64,
    external_heat: f64, body: UniformEnthalpyStepInput<'_>, chemistry: GasReactionSolve,
    controls: IsobaricBodySolve, steps: Vec<IsobaricBodyAdvanceSubstep<S>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<IsobaricBodyAdvanceCandidate<S>, IsobaricBodyAdvanceError<E>> {
    let (mut q, mut evaluated, mut source) = (Sum::default(), Sum::default(), Sum::default());
    let (mut stored, mut stored_residual, mut total_h, mut total_u) =
        (Sum::default(), Sum::default(), Sum::default(), Sum::default());
    let (mut entropy, mut body_heat_entropy, mut external_entropy) =
        (Sum::default(), Sum::default(), Sum::default());
    let (mut applied_contact, mut evaluated_contact, mut gas_given, mut body_given) =
        (Sum::default(), Sum::default(), Sum::default(), Sum::default());
    for interval in &steps {
        poll(cancelled)?;
        let step = &interval.step;
        let exchange = step.boundary.exchanged_heat_j;
        let tb = interval.thermal.state.temperature_k();
        if interval.thermal.external_heat_j.to_bits() != (-exchange).to_bits() {
            return Err(invalid("finite-body history lost reciprocal heat").into());
        }
        q.add(exchange)?;
        evaluated.add(step.boundary.evaluated_heat_j)?;
        source.add(step.boundary.source_residual_j)?;
        stored.add(step.thermal_energy_change_j)?;
        stored_residual.add(step.thermal_energy_residual_j)?;
        total_h.add(step.total_enthalpy_residual_j)?;
        total_u.add(step.total_energy_residual_j)?;
        entropy.add(step.thermal_entropy_change_j_k)?;
        body_heat_entropy.add((interval.prescribed_body_heat_j - exchange) / tb)?;
        external_entropy.add(interval.prescribed_gas_heat_j / interval.gas_temperature_k)?;
        external_entropy.add(interval.prescribed_body_heat_j / tb)?;
        applied_contact.add(step.applied_contact_entropy_j_k)?;
        evaluated_contact.add(step.boundary.contact_entropy_j_k)?;
        gas_given.add(interval.prescribed_gas_heat_j)?;
        body_given.add(interval.prescribed_body_heat_j)?;
    }
    let final_phase = steps.last().expect("positive admitted substep count").thermal.state;
    // The phase chart's h is the body's authoritative stored coordinate.
    // Independently compare its actual endpoint change with accepted history.
    let actual_body_storage = body.mass_kg * (final_phase.specific_enthalpy_j_kg()
        - body.initial.specific_enthalpy_j_kg());
    let (gas, trajectory) = super::summarize(law, initial, gas, duration,
        external_heat + q.value, chemistry, steps, |interval| &interval.step.boundary.chemistry, cancelled)?;
    let thermal_residual = (stored.value + q.value) - body.internal_heat_j;
    let enthalpy_residual = (trajectory.enthalpy_change_j + stored.value)
        - external_heat - body.internal_heat_j;
    let energy_residual = (trajectory.stored_energy_change_j + stored.value)
        + trajectory.work_by_gas_j - external_heat - body.internal_heat_j;
    let thermal_entropy_residual = entropy.value - body_heat_entropy.value;
    let combined_entropy = trajectory.entropy_change_j_k + entropy.value - external_entropy.value;
    let energy_tolerance = controls.heat.exchange.energy_tolerance_j;
    let entropy_tolerance = controls.entropy_tolerance_j_k;
    if [thermal_residual, enthalpy_residual, energy_residual, stored_residual.value,
        total_h.value, total_u.value, gas_given.value - external_heat,
        body_given.value - body.internal_heat_j, actual_body_storage - stored.value].into_iter()
        .any(|value| !value.is_finite() || value.abs() > energy_tolerance)
        || source.value.abs() > controls.heat.exchange.source_tolerance_j
        || (q.value - evaluated.value).abs() > controls.heat.exchange.source_tolerance_j
        || [thermal_entropy_residual, combined_entropy, applied_contact.value,
            evaluated_contact.value].into_iter()
            .any(|value| !value.is_finite() || value < -entropy_tolerance)
    {
        return Err(invalid("whole gas/body interval failed original heat, enthalpy, pressure-work or entropy budgets").into());
    }
    poll(cancelled)?;
    Ok(IsobaricBodyAdvanceCandidate { gas,
        thermal: UniformEnthalpyStep { state: final_phase, external_heat_j: -q.value,
            energy_residual_tolerance_j: energy_tolerance,
            report: IsobaricBodyAdvanceReport { trajectory,
                prescribed_gas_heat_j: external_heat, prescribed_body_heat_j: body.internal_heat_j,
                exchanged_heat_j: q.value, evaluated_contact_heat_j: evaluated.value,
                source_residual_j: source.value, thermal_energy_change_j: stored.value,
                thermal_energy_residual_j: thermal_residual,
                total_enthalpy_residual_j: enthalpy_residual, total_energy_residual_j: energy_residual,
                thermal_entropy_change_j_k: entropy.value,
                thermal_entropy_residual_j_k: thermal_entropy_residual,
                applied_contact_entropy_j_k: applied_contact.value,
                evaluated_contact_entropy_j_k: evaluated_contact.value,
                prescribed_heat_entropy_j_k: external_entropy.value,
                total_entropy_residual_j_k: combined_entropy,
            } },
    })
}

fn compare<S>(a: &IsobaricBodyAdvanceCandidate<S>, b: &IsobaricBodyAdvanceCandidate<S>,
    initial: EquilibriumPhaseState, policy: IsobaricBodyAdvanceSolve,
) -> Result<f64, IsobaricAdvanceError> {
    let mut error = super::compare(&a.gas, &a.thermal.report.trajectory,
        &b.gas, &b.thermal.report.trajectory, policy.trajectory)?;
    error = error.max(scaled(a.thermal.state.temperature_k(), b.thermal.state.temperature_k(),
        policy.trajectory.absolute_temperature_k, policy.trajectory.relative_tolerance)?);
    error = error.max(scaled(a.thermal.state.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg(),
        b.thermal.state.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg(),
        policy.absolute_body_enthalpy_j_kg, policy.trajectory.relative_tolerance)?);
    Ok(error)
}
