//! Time adaptation that cannot confuse a latent plateau with an unchanged state.

use super::{EnthalpyHistoryRun, EnthalpyStepConfig, EnthalpyStepSolution,
    HeterogeneousEnthalpyError, conduction, invalid, poll};
use super::super::events::{self, adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep,
    admit_proposal, component_error, select_step_with_recovery, temperature_error}};
use fs_exec::Cx;
use crate::transient::enthalpy::EnthalpyError;
use crate::transient::enthalpy::adaptive::combine_entropy;
use crate::transient::enthalpy::entropy::{EnthalpyEntropyBalance, EnthalpyEntropyError};

/// The accepted physical fine trajectory and its separate entropy accounts.
/// Coarse comparison and rejected trials contribute to none of these accounts.
#[derive(Debug, Clone)]
pub struct EntropyAdaptiveHistoryStep<T> {
    /// Original native half-step reports, accepted clock and adaptation work.
    pub physical: AdaptiveHistoryStep<T>,
    /// Native entropy account for each accepted half, in chronological order.
    pub half_entropy: [EnthalpyEntropyBalance; 2],
    /// Whole-pair account with exact original-to-final caloric storage change.
    /// Physical exchanges and signed numerical defects retain their identity.
    pub entropy: EnthalpyEntropyBalance,
}

impl<T> EntropyAdaptiveHistoryStep<T> {
    pub(super) fn from_halves(
        selected: AdaptiveHistoryStep<(T, EnthalpyEntropyBalance)>,
        entropy: EnthalpyEntropyBalance,
    ) -> Self {
        let AdaptiveHistoryStep {
            start_time_s, midpoint_time_s, end_time_s, first_half, second_half,
            normalized_error, rejected_attempts, endpoint_solves, suggested_step_s,
        } = selected;
        Self {
            physical: AdaptiveHistoryStep {
                start_time_s, midpoint_time_s, end_time_s,
                first_half: first_half.0, second_half: second_half.0,
                normalized_error, rejected_attempts, endpoint_solves, suggested_step_s,
            },
            half_entropy: [first_half.1, second_half.1],
            entropy,
        }
    }
}

/// Original history admission or native entropy refusal; neither publishes state.
#[derive(Debug, Clone, PartialEq)]
pub enum EntropyAdaptiveHistoryError {
    /// Schedule, material, nonlinear, energy, cancellation or temporal admission.
    History(HeterogeneousEnthalpyError),
    /// Native physical entropy account or its whole-pair allowance failed.
    Entropy(EnthalpyEntropyError),
}

impl From<HeterogeneousEnthalpyError> for EntropyAdaptiveHistoryError {
    fn from(error: HeterogeneousEnthalpyError) -> Self { Self::History(error) }
}
impl From<EnthalpyError> for EntropyAdaptiveHistoryError {
    fn from(error: EnthalpyError) -> Self { Self::History(error.into()) }
}
impl From<EnthalpyEntropyError> for EntropyAdaptiveHistoryError {
    fn from(error: EnthalpyEntropyError) -> Self { Self::Entropy(error) }
}
impl std::fmt::Display for EntropyAdaptiveHistoryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::History(error) => write!(f, "{error}"),
            Self::Entropy(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for EntropyAdaptiveHistoryError {}

pub(super) fn entropy_conduction(error: crate::ConductionError) -> EntropyAdaptiveHistoryError {
    conduction(error).into()
}

pub(super) fn admit_entropy_tolerance(tolerance: f64) -> Result<(), EntropyAdaptiveHistoryError> {
    // Reuse the native allowance admission before any candidate solve.
    EnthalpyEntropyBalance::default().admit(tolerance)?;
    Ok(())
}

pub(super) fn accepted_entropy(
    cx: &Cx<'_>,
    run: &EnthalpyHistoryRun<'_, '_, '_>,
    final_h: &[f64],
    halves: [EnthalpyEntropyBalance; 2],
    tolerance: f64,
) -> Result<EnthalpyEntropyBalance, EntropyAdaptiveHistoryError> {
    let storage = run.storage.thermal_entropy_change_j_k(cx, &run.specific_h, final_h)?;
    let entropy = combine_entropy(halves[0], halves[1], storage)?;
    entropy.admit(tolerance)?;
    Ok(entropy)
}

// Match the standalone native adaptive owner's physical recovery policy.
// Source support, material identity, invalid controls, budgets and cancellation
// are admissions, not reasons to silently seek a different physical problem.
fn recoverable_physical(error: &EnthalpyError) -> bool {
    matches!(error, EnthalpyError::NotConverged(_)
        | EnthalpyError::EndpointResidual { .. } | EnthalpyError::EnergyBalance { .. })
}

pub(super) fn recoverable_history(error: &HeterogeneousEnthalpyError) -> bool {
    matches!(error, HeterogeneousEnthalpyError::Enthalpy(error) if recoverable_physical(error))
}

pub(super) fn recoverable_entropy_history(error: &EntropyAdaptiveHistoryError) -> bool {
    match error {
        EntropyAdaptiveHistoryError::History(error) => recoverable_history(error),
        EntropyAdaptiveHistoryError::Entropy(EnthalpyEntropyError::Physical(error)) =>
            recoverable_physical(error),
        _ => false,
    }
}

/// Additional error coordinates for a primary-enthalpy trajectory.
/// All three gates (temperature, enthalpy, liquid fraction) must pass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseErrorTolerance {
    /// Positive absolute local specific-enthalpy allowance, J/kg.
    pub absolute_enthalpy_j_kg: f64,
    /// Relative allowance against the larger coarse/fine enthalpy CHANGE from
    /// history, not the arbitrary enthalpy reference offset. Must be in [0,1).
    pub relative_enthalpy: f64,
    /// Positive absolute allowance on liquid mass fraction, at most one.
    pub absolute_liquid_fraction: f64,
}
impl PhaseErrorTolerance {
    /// Validate every caloric tolerance before any candidate solve.
    pub fn validate(self) -> Result<(), HeterogeneousEnthalpyError> {
        if !self.absolute_enthalpy_j_kg.is_finite() || self.absolute_enthalpy_j_kg <= 0.0
            || !self.relative_enthalpy.is_finite() || !(0.0..1.0).contains(&self.relative_enthalpy)
            || !self.absolute_liquid_fraction.is_finite()
            || self.absolute_liquid_fraction <= 0.0 || self.absolute_liquid_fraction > 1.0
        { return Err(invalid("phase time-error control requires positive finite J/kg and fraction allowances and relative tolerance in [0,1)")); }
        Ok(())
    }
}

impl EnthalpyHistoryRun<'_, '_, '_> {
    /// Advance two original implicit half steps, accepting only a small estimated
    /// discrepancy in temperature, specific enthalpy AND liquid fraction.
    ///
    /// Each candidate uses the same retained phase charts, reference masses,
    /// k(T), phase multipliers and contacts. The event controller limits the
    /// interval first; held loads retain their exact interval-end left limits.
    /// No apparent capacity, extrapolated caloric state or phase reset is used.
    /// The accepted pair consumes two lifetime slots and publishes atomically.
    ///
    /// This is a local estimator, not a phase-transition-time certificate or a
    /// global error bound. Kinks in constitutive charts need not have smooth
    /// first-order asymptotics. Native nonlinear nonconvergence and endpoint
    /// residual/energy refusals bisect within the same attempt and step limits.
    /// Model/domain, cancellation and budget refusals propagate immediately.
    /// Exhausted recovery returns the original typed physical refusal.
    ///
    /// # Errors
    /// Invalid policy/time, cancellation, insufficient accepted slots, original
    /// solver refusal, or unmet temporal accuracy at the floor/attempt limit.
    /// Every failure preserves enthalpy, temperature, liquid fraction and time.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target_s: f64, proposed_step_s: f64,
        step: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
        -> Result<AdaptiveHistoryStep<EnthalpyStepSolution>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.accepted_steps).map_err(conduction)?;
        phase.validate()?;
        let start = self.time_s;
        let cap = admit_proposal(cx, self.accepted_steps, self.max_steps, proposed_step_s, adaptive).map_err(conduction)?;
        let end = events::next_endpoint(cx, self.problem.boundary, &self.problem.source,
            std::iter::empty(), start, target_s, cap).map_err(conduction)?;
        let solve = |old: &[f64], from: f64, to: f64| -> Result<EnthalpyStepSolution, HeterogeneousEnthalpyError> {
            let boundary = events::boundary_before(cx, self.problem.boundary, to).map_err(conduction)?;
            self.storage.admit_prescribed_boundary(cx, &boundary)?;
            let source = events::source_before(&self.problem.source, to, old.len()).map_err(conduction)?;
            Ok(self.storage.advance(cx, self.conduction_problem(&boundary, &source)?,
                self.problem.interfaces, old, to - from, step)?)
        };
        let selected = select_step_with_recovery(cx, start, end, adaptive, conduction,
            recoverable_history, |middle, end, solves| {
            *solves += 1;
            let coarse = solve(&self.specific_h, start, end)?;
            *solves += 1;
            let first = solve(&self.specific_h, start, middle)?;
            *solves += 1;
            let second = solve(&first.specific_enthalpy_j_kg, middle, end)?;
            let error = phase_error(cx, self, &coarse, &second, adaptive, phase)?;
            Ok((first, second, error))
        })?;
        self.publish(cx, selected.end_time_s, &selected.second_half)?;
        self.accepted_steps += 1; // Two reserved physical half steps, one atomic publication.
        Ok(selected)
    }
    /// Adapt using the original energy-only coarse comparison and entropy-admitted
    /// native fine steps. The complete accepted pair must also pass ONE entropy
    /// allowance before any state or accepted-step counter is published.
    ///
    /// Exact storage is queried between the original and final caloric states.
    /// Reservoir exchange, unclassified source entropy, physical generation and
    /// signed temporal/algebraic defects are accumulated only from accepted fine
    /// intervals. Rejected candidates never enter the returned entropy balance.
    ///
    /// Recoverable physical failures use the same bounded bisection as ordinary
    /// enthalpy histories; entropy and model/domain refusals still propagate
    /// immediately. All retries retain the original floor, event limits and
    /// work budget. Every failure preserves accepted history.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step_with_entropy(
        &mut self, cx: &Cx<'_>, target_s: f64, proposed_step_s: f64,
        step: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig,
        phase: PhaseErrorTolerance, entropy_tolerance_j_k: f64,
    ) -> Result<EntropyAdaptiveHistoryStep<EnthalpyStepSolution>, EntropyAdaptiveHistoryError> {
        poll(cx, self.accepted_steps).map_err(entropy_conduction)?;
        phase.validate()?;
        admit_entropy_tolerance(entropy_tolerance_j_k)?;
        let start = self.time_s;
        let cap = admit_proposal(cx, self.accepted_steps, self.max_steps,
            proposed_step_s, adaptive).map_err(entropy_conduction)?;
        let end = events::next_endpoint(cx, self.problem.boundary, &self.problem.source,
            std::iter::empty(), start, target_s, cap).map_err(entropy_conduction)?;
        let solve = |old: &[f64], from: f64, to: f64, entropy: bool|
            -> Result<(EnthalpyStepSolution, Option<EnthalpyEntropyBalance>), EntropyAdaptiveHistoryError>
        {
            let boundary = events::boundary_before(cx, self.problem.boundary, to)
                .map_err(entropy_conduction)?;
            self.storage.admit_prescribed_boundary(cx, &boundary)?;
            let source = events::source_before(&self.problem.source, to, old.len())
                .map_err(entropy_conduction)?;
            let problem = self.conduction_problem(&boundary, &source)?;
            if entropy {
                let (physical, balance) = self.storage.advance_with_entropy(cx, problem,
                    self.problem.interfaces, old, to - from, step, entropy_tolerance_j_k)?;
                Ok((physical, Some(balance)))
            } else {
                Ok((self.storage.advance(cx, problem, self.problem.interfaces,
                    old, to - from, step)?, None))
            }
        };
        let selected = select_step_with_recovery(cx, start, end, adaptive, entropy_conduction,
            recoverable_entropy_history, |middle, end, solves| -> Result<_, EntropyAdaptiveHistoryError> {
                *solves += 1;
                let coarse = solve(&self.specific_h, start, end, false)?.0;
                *solves += 1;
                let (first, first_entropy) = solve(&self.specific_h, start, middle, true)?;
                *solves += 1;
                let (second, second_entropy) = solve(&first.specific_enthalpy_j_kg, middle, end, true)?;
                let error = phase_error(cx, self, &coarse, &second, adaptive, phase)?;
                Ok(((first, first_entropy.expect("requested native fine-step entropy")),
                    (second, second_entropy.expect("requested native fine-step entropy")), error))
            })?;
        let entropy = accepted_entropy(cx, self, &selected.second_half.0.specific_enthalpy_j_kg,
            [selected.first_half.1, selected.second_half.1], entropy_tolerance_j_k)?;
        let selected = EntropyAdaptiveHistoryStep::from_halves(selected, entropy);
        self.publish(cx, selected.physical.end_time_s, &selected.physical.second_half)?;
        self.accepted_steps += 1;
        Ok(selected)
    }
}

pub(super) fn phase_error(cx: &Cx<'_>, run: &EnthalpyHistoryRun<'_, '_, '_>,
    coarse: &EnthalpyStepSolution, fine: &EnthalpyStepSolution,
    adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
    -> Result<f64, HeterogeneousEnthalpyError>
{
    phase_field_error(cx, run,
        PhaseFields { h: &coarse.specific_enthalpy_j_kg, temperature: &coarse.temperature,
            fraction: &coarse.liquid_mass_fraction },
        PhaseFields { h: &fine.specific_enthalpy_j_kg, temperature: &fine.temperature,
            fraction: &fine.liquid_mass_fraction }, adaptive, phase)
}

/// Borrowed caloric fields let every time method use the same error coordinates.
pub(super) struct PhaseFields<'a> {
    pub h: &'a [f64],
    pub temperature: &'a [f64],
    pub fraction: &'a [f64],
}

pub(super) fn phase_field_error(cx: &Cx<'_>, run: &EnthalpyHistoryRun<'_, '_, '_>,
    coarse: PhaseFields<'_>, fine: PhaseFields<'_>,
    adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
    -> Result<f64, HeterogeneousEnthalpyError>
{
    let n = run.specific_h.len();
    if [coarse.h.len(), fine.h.len(), coarse.fraction.len(), fine.fraction.len()]
        .iter().any(|&len| len != n)
    { return Err(invalid("adaptive caloric fields must preserve the original nodal shape")); }
    let mut error = temperature_error(cx, &run.temperature, &coarse.temperature,
        &fine.temperature, adaptive).map_err(conduction)?;
    for i in 0..n {
        if i % 512 == 0 { poll(cx, i).map_err(conduction)?; }
        error = error.max(component_error(run.specific_h[i], coarse.h[i],
            fine.h[i], phase.absolute_enthalpy_j_kg, phase.relative_enthalpy).map_err(conduction)?);
        error = error.max(component_error(run.liquid_mass_fraction[i], coarse.fraction[i],
            fine.fraction[i], phase.absolute_liquid_fraction, 0.0).map_err(conduction)?);
    }
    Ok(error)
}
