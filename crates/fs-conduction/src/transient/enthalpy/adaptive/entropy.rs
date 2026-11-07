//! Entropy on the accepted fine trajectory of the original adaptive controller.
//!
//! Coarse comparison steps remain energy-only. The two fine steps use native
//! final-operator entropy; only an accepted pair enters the retained account.
//! Storage is remeasured from the original enthalpies, including bonded films.
//! Double-double sums survive resume and keep cancellation between intervals.

use super::*;
use crate::transient::enthalpy::entropy::{
    EnthalpyEntropyBalance, EnthalpyEntropyError, admit_radiation,
};
use fs_math::dd::Dd;
use crate::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use crate::transient::enthalpy::sdirk::EnthalpySdirk2Input;

/// Original adaptive work/heat report and the whole accepted entropy history.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveEnthalpyEntropyReport {
    pub thermal: AdaptiveEnthalpyReport,
    pub entropy: EnthalpyEntropyBalance,
}

/// Resumable entropy-admitted counterpart of the SAME adaptive BE/SDIRK2 controller.
///
/// The physical history is available read-only so an unaudited step cannot
/// silently enter this account. Cloning retains controller state and compensated
/// entropy sums. A failure preserves all earlier accepted pairs; the failing
/// pair changes neither physical history nor entropy. Observer cancellation
/// after publication retains that accepted pair, as in the original controller.
#[derive(Debug, Clone)]
pub struct AdaptiveEnthalpyEntropyState<'s, 'm, 'c> {
    history: AdaptiveEnthalpyState<'s, 'm, 'c>,
    account: EntropyAccount,
}

impl<'m, 'c> EnthalpyBackwardEuler<'m, 'c> {
    /// Prepare the existing adaptive history with one absolute J/K allowance.
    ///
    /// The allowance applies to each native fine endpoint and to the COMPLETE
    /// accepted window. It is not multiplied by the accepted-step count.
    /// It bounds the measured entropy identity and signed admission checks,
    /// not spatial/time error or material uncertainty. Original chart, energy
    /// and adaptive error controls remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn adaptive_state_with_entropy<'s>(
        &'s self,
        cx: &Cx<'_>,
        initial_specific_h: &[f64],
        initial_time_s: f64,
        initial_step_s: f64,
        config: AdaptiveEnthalpyConfig,
        entropy_tolerance_j_k: f64,
    ) -> Result<AdaptiveEnthalpyEntropyState<'s, 'm, 'c>, EnthalpyEntropyError> {
        EnthalpyEntropyBalance::default().admit(entropy_tolerance_j_k)?;
        let history = self.adaptive_state(
            cx, initial_specific_h, initial_time_s, initial_step_s, config,
        )?;
        // Query the same caloric/film owners before accepting an entropy history.
        self.thermal_entropy_change_j_k(cx, initial_specific_h, initial_specific_h)?;
        Ok(AdaptiveEnthalpyEntropyState {
            history,
            account: EntropyAccount {
                sums: [Dd::ZERO; 10],
                balance: EnthalpyEntropyBalance::default(),
                tolerance_j_k: entropy_tolerance_j_k,
            },
        })
    }
}

impl<'s, 'm, 'c> AdaptiveEnthalpyEntropyState<'s, 'm, 'c> {
    /// Read original accepted fields, clock, controller and work counters.
    #[must_use]
    pub fn history(&self) -> &AdaptiveEnthalpyState<'s, 'm, 'c> {
        &self.history
    }

    /// Whole accepted-window entropy; rejected/comparison work is excluded.
    #[must_use]
    pub const fn entropy_balance(&self) -> EnthalpyEntropyBalance {
        self.account.balance
    }

    /// Advance SDIRK2 with native stage entropy and whole-history admission.
    ///
    /// Uses the same three-field controller and native stages as the energy-only
    /// `AdaptiveEnthalpyState::advance_to_sdirk2`. The coarse comparison remains
    /// energy-only; both fine intervals use `advance_sdirk2_with_entropy`, with
    /// their original stage temperatures, transport and reservoir quadrature.
    /// Only an accepted fine pair enters the retained entropy sums.
    ///
    /// The original absolute J/K allowance applies to each fine interval and
    /// once to the entire accepted window, without step-count inflation.
    /// Original-to-final caloric and film entropy is measured once. Signed
    /// temporal defects remain distinct from physical production.
    /// Backward Euler and SDIRK2 may alternate between event calls without
    /// resetting accepted enthalpy, heat or entropy. L-stability alone is not
    /// an entropy theorem; native entropy refusals propagate transactionally.
    pub fn advance_to_sdirk2(
        &mut self,
        cx: &Cx<'_>,
        input: EnthalpySdirk2Input<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyEntropyReport, EnthalpyEntropyError> {
        self.advance_to_sdirk2_observed(
            cx, input, end_time_s, max_trials, |_| ControlFlow::Continue(()),
        )
    }

    /// Observe the authentic fine SDIRK2 stages after physical state, heat and
    /// whole-window entropy are published together. Coarse or rejected stages
    /// never reach the observer. Cancellation after publication retains the
    /// complete accepted pair. Physical laws are fixed during this event call;
    /// stage material/contact/radiation responses retain their actual states.
    pub fn advance_to_sdirk2_observed(
        &mut self,
        cx: &Cx<'_>,
        input: EnthalpySdirk2Input<'_>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpySdirk2AcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyEntropyReport, EnthalpyEntropyError> {
        let thermal = self.history.advance_sdirk2_accounted(
            cx, input, end_time_s, max_trials, Some(&mut self.account), observe,
        )?;
        Ok(AdaptiveEnthalpyEntropyReport { thermal, entropy: self.account.balance })
    }

    /// Advance the original fixed-law spatial BE history to a declared event.
    ///
    /// Only the two fine endpoints are entropy-admitted; a coarse estimator is
    /// never mistaken for the accepted physical trajectory. Entropy refusals
    /// propagate unchanged. Original recoverable Newton/energy refusals retain
    /// their bounded step backoff. Source laws may change BETWEEN calls, while
    /// original masses, caloric charts and the accumulated entropy persist.
    pub fn advance_to(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyEntropyReport, EnthalpyEntropyError> {
        self.advance_to_observed(
            cx, problem, interfaces, end_time_s, max_trials, |_| ControlFlow::Continue(()),
        )
    }

    /// Observe the same accepted fine fields AFTER both accounts are published.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_observed(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyEntropyReport, EnthalpyEntropyError> {
        self.advance_to_with_radiation_observed(
            cx, problem, interfaces, EnthalpyRadiation::None, end_time_s, max_trials, observe,
        )
    }

    /// Include original fixed or material-state pointwise ambient or enclosure
    /// radiation. Internal reflection and explicit external black reservoirs
    /// retain their native heat accounts; photon entropy is not inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_radiation(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        radiation: EnthalpyRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyEntropyReport, EnthalpyEntropyError> {
        self.advance_to_with_radiation_observed(
            cx, problem, interfaces, radiation, end_time_s, max_trials,
            |_| ControlFlow::Continue(()),
        )
    }

    /// Observe only accepted radiative fine pairs with their original heat reports.
    /// Mean-patch radiation requires another entropy chart and refuses before
    /// any endpoint solve or accepted-history change.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to_with_radiation_observed(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        radiation: EnthalpyRadiation<'_>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpyAcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyEntropyReport, EnthalpyEntropyError> {
        admit_radiation(radiation)?;
        let unsupported = || EnthalpyError::InvalidInput(
            "entropy does not support mean-temperature radiation approximations",
        );
        let storage = self.history.storage;
        let config = self.history.config.endpoint;
        let tolerance = self.account.tolerance_j_k;
        let thermal = self.history.advance_with_account(
            cx, end_time_s, max_trials,
            |h, dt, fine| {
                if fine {
                    let (solution, entropy) = storage.advance_with_radiation_and_entropy(
                        cx, problem, interfaces, h, dt, config, radiation, tolerance,
                    )?;
                    Ok((solution, Some(entropy)))
                } else {
                    let solution = match radiation {
                        EnthalpyRadiation::None => storage.advance(
                            cx, problem, interfaces, h, dt, config,
                        )?,
                        EnthalpyRadiation::PointwiseAmbient(owner) =>
                            storage.advance_with_pointwise_ambient_radiation(
                                cx, problem, interfaces, h, dt, config, owner,
                            )?,
                        EnthalpyRadiation::MaterialAmbient(owner) =>
                            storage.advance_with_material_ambient_radiation(
                                cx, problem, interfaces, h, dt, config, owner,
                            )?,
                        EnthalpyRadiation::PointwiseEnclosure(owner) =>
                            storage.advance_with_pointwise_enclosure_radiation(
                                cx, problem, interfaces, h, dt, config, owner,
                            )?,
                        EnthalpyRadiation::MaterialEnclosure(owner) =>
                            storage.advance_with_material_enclosure_radiation(
                                cx, problem, interfaces, h, dt, config, owner,
                            )?,
                        _ => return Err(unsupported().into()),
                    };
                    Ok((solution, None))
                }
            },
            Some(&mut self.account),
            observe,
        )?;
        Ok(AdaptiveEnthalpyEntropyReport { thermal, entropy: self.account.balance })
    }
}

#[derive(Debug, Clone)]
pub(super) struct EntropyAccount {
    sums: [Dd; 10],
    balance: EnthalpyEntropyBalance,
    tolerance_j_k: f64,
}

impl EntropyAccount {
    pub(super) const fn tolerance_j_k(&self) -> f64 {
        self.tolerance_j_k
    }

    pub(super) fn extended(
        &self,
        storage: &EnthalpyBackwardEuler<'_, '_>,
        cx: &Cx<'_>,
        initial_h: &[f64],
        final_h: &[f64],
        first: Option<EnthalpyEntropyBalance>,
        second: Option<EnthalpyEntropyBalance>,
    ) -> Result<Self, EnthalpyEntropyError> {
        let missing = || EnthalpyError::InvalidInput("native fine-step entropy account is missing");
        let first = first.ok_or_else(missing)?;
        let second = second.ok_or_else(missing)?;
        let exact_storage = storage.thermal_entropy_change_j_k(cx, initial_h, final_h)?;
        let mut sums = self.sums;
        add(&mut sums, first)?;
        add(&mut sums, second)?;
        let balance = finish(
            sums, exact_storage,
            self.balance.has_unclassified_inputs
                || first.has_unclassified_inputs || second.has_unclassified_inputs,
        )?;
        balance.admit(self.tolerance_j_k)?;
        Ok(Self { sums, balance, tolerance_j_k: self.tolerance_j_k })
    }
}

/// Shared accepted-pair composition for scheduled BE and SDIRK2 histories.
/// Storage is supplied by the ORIGINAL initial-to-final caloric/film query.
pub(crate) fn combine_entropy(
    first: EnthalpyEntropyBalance,
    second: EnthalpyEntropyBalance,
    exact_storage_j_k: f64,
) -> Result<EnthalpyEntropyBalance, EnthalpyError> {
    let mut sums = [Dd::ZERO; 10];
    add(&mut sums, first)?;
    add(&mut sums, second)?;
    finish(sums, exact_storage_j_k, first.has_unclassified_inputs || second.has_unclassified_inputs)
}

fn add(sums: &mut [Dd; 10], value: EnthalpyEntropyBalance) -> Result<(), EnthalpyError> {
    let terms = [
        value.dirichlet_entropy_input_j_k,
        value.robin_entropy_input_j_k,
        value.radiation_entropy_input_j_k,
        value.conduction_generation_j_k,
        value.robin_generation_j_k,
        value.radiation_generation_j_k,
        value.source_at_node_j_k,
        value.neumann_at_node_j_k,
        value.temporal_defect_j_k,
        value.algebraic_residual_j_k,
    ];
    for (sum, term) in sums.iter_mut().zip(terms) {
        *sum = *sum + Dd::from_f64(finite(term)?);
        finite(sum.to_f64())?;
    }
    Ok(())
}

fn finish(
    sums: [Dd; 10], exact_storage: f64, has_unclassified_inputs: bool,
) -> Result<EnthalpyEntropyBalance, EnthalpyError> {
    let mut residual = Dd::from_f64(finite(exact_storage)?);
    for sum in sums { residual = residual - sum; }
    let [dirichlet, robin, radiation, conduction, robin_generation, radiation_generation,
        source, neumann, temporal, algebraic] = sums.map(Dd::to_f64);
    Ok(EnthalpyEntropyBalance {
        thermal_storage_change_j_k: exact_storage,
        dirichlet_entropy_input_j_k: dirichlet,
        robin_entropy_input_j_k: robin,
        radiation_entropy_input_j_k: radiation,
        conduction_generation_j_k: conduction,
        robin_generation_j_k: robin_generation,
        radiation_generation_j_k: radiation_generation,
        source_at_node_j_k: source,
        neumann_at_node_j_k: neumann,
        has_unclassified_inputs,
        temporal_defect_j_k: temporal,
        algebraic_residual_j_k: algebraic,
        balance_residual_j_k: finite(residual.to_f64())?,
    })
}
