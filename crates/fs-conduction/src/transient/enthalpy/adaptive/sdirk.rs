//! SDIRK2 steps through the existing enthalpy adaptive controller.

use super::*;
use super::entropy::EntropyAccount;
use crate::transient::enthalpy::entropy::{EnthalpyEntropyError, admit_radiation};
use crate::transient::enthalpy::sdirk::{EnthalpySdirk2Input, EnthalpySdirk2Step};

/// Two accepted fine SDIRK2 steps with their authentic stages and heat integrals.
///
/// Stage one is an internal tableau state, not an accepted time interval.
/// Only the two complete steps enter the trajectory's external heat account.
/// Comparison and rejected steps are never exposed. The references remain
/// valid only during observation; saving a trajectory is a caller decision.
#[derive(Debug, Clone, Copy)]
pub struct EnthalpySdirk2AcceptedTrial<'a> {
    pub accepted_trial: usize,
    pub start_time_s: f64,
    pub midpoint_time_s: f64,
    pub end_time_s: f64,
    pub first_half: &'a EnthalpySdirk2Step,
    pub second_half: &'a EnthalpySdirk2Step,
    /// Raw full/fine differences and the normalized fine-error estimate.
    /// This is not a bound on accumulated, spatial or material-model error.
    pub error: EnthalpyTimeError,
}

impl AdaptiveEnthalpyState<'_, '_, '_> {
    /// Last accepted SDIRK2 HALF-step, including its two native stages.
    ///
    /// The heat integrals and storage change cover that half step only; the
    /// adaptive report covers the whole accepted trajectory. Returns `None`
    /// initially or when the latest accepted step used backward Euler.
    #[must_use]
    pub fn last_sdirk2_step(&self) -> Option<&EnthalpySdirk2Step> {
        match self.endpoint.as_ref()? {
            AdaptiveStep::Sdirk2(step) => Some(step),
            AdaptiveStep::BackwardEuler(_) => None,
        }
    }

    /// Advance with second-order, L-stable SDIRK2 under fixed physical laws.
    ///
    /// Uses the SAME history, event clipping, three-field error controls,
    /// rejected-step backoff and transactional publication as backward Euler.
    /// Every trial compares one native full step with two native half steps;
    /// the fine state is accepted without extrapolating enthalpy or phase.
    /// The difference divided by three estimates the fine local error for
    /// equal halves in a smooth regime; rounded unequal halves retain their
    /// actual cubic weights. Step adaptation uses the third-order local error.
    ///
    /// Source and boundary laws are fixed during this call. Material responses,
    /// contacts and radiation still evaluate each stage's actual physical state.
    /// Resolve forcing jumps by separate event calls. Constitutive corners can
    /// reduce order; neither the estimator nor L-stability certifies phase
    /// monotonicity, entropy, global time accuracy or spatial resolution.
    ///
    /// `max_trials` bounds comparisons; each of the three native step attempts
    /// has at most two stages under `config.endpoint`'s original work gates.
    /// The method may change between calls without resetting physical history.
    pub fn advance_to_sdirk2(
        &mut self,
        cx: &Cx<'_>,
        input: EnthalpySdirk2Input<'_>,
        end_time_s: f64,
        max_trials: usize,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_to_sdirk2_observed(cx, input, end_time_s, max_trials, |_| {
            ControlFlow::Continue(())
        })
    }

    /// Observe only complete accepted fine pairs, after committing their
    /// physical state, integrated heat and controller. Stopping or requesting
    /// cancellation in the observer retains that pair for an exact resume.
    pub fn advance_to_sdirk2_observed(
        &mut self,
        cx: &Cx<'_>,
        input: EnthalpySdirk2Input<'_>,
        end_time_s: f64,
        max_trials: usize,
        observe: impl FnMut(EnthalpySdirk2AcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyError> {
        self.advance_sdirk2_accounted(cx, input, end_time_s, max_trials, None, observe)
            .map_err(|error| match error {
                EnthalpyEntropyError::Physical(error) => error,
                _ => EnthalpyError::InvalidInput(
                    "unexpected entropy gate on a physical-only trajectory",
                ),
            })
    }

    // Both public histories use the same native stages, estimator, heat and
    // observer path. Entropy adds no work to the coarse comparison or Newton
    // iterations; only the two fine steps enter a candidate physical account.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_sdirk2_accounted(
        &mut self,
        cx: &Cx<'_>,
        input: EnthalpySdirk2Input<'_>,
        end_time_s: f64,
        max_trials: usize,
        entropy: Option<&mut EntropyAccount>,
        mut observe: impl FnMut(EnthalpySdirk2AcceptedTrial<'_>) -> ControlFlow<()>,
    ) -> Result<AdaptiveEnthalpyReport, EnthalpyEntropyError> {
        let storage = self.storage;
        let config = self.config.endpoint;
        let tolerance = entropy.as_deref().map(EntropyAccount::tolerance_j_k);
        if tolerance.is_some() {
            admit_radiation(input.radiation)?;
        }
        self.advance_native_with_account(
            cx,
            end_time_s,
            max_trials,
            AdaptiveOrder::Sdirk2,
            |h, dt, fine| match (fine, tolerance) {
                (true, Some(tolerance)) => storage
                    .advance_sdirk2_with_entropy(cx, [input; 2], h, dt, config, tolerance)
                    .map(|(step, account)| (AdaptiveStep::Sdirk2(step), Some(account))),
                _ => Ok((
                    AdaptiveStep::Sdirk2(storage.advance_sdirk2(cx, [input; 2], h, dt, config)?),
                    None,
                )),
            },
            entropy,
            |trial| {
                let (AdaptiveStep::Sdirk2(first_half), AdaptiveStep::Sdirk2(second_half)) =
                    (trial.first_half, trial.second_half)
                else {
                    unreachable!("SDIRK2 advancement retains SDIRK2 reports")
                };
                observe(EnthalpySdirk2AcceptedTrial {
                    accepted_trial: trial.accepted_trial,
                    start_time_s: trial.start_time_s,
                    midpoint_time_s: trial.midpoint_time_s,
                    end_time_s: trial.end_time_s,
                    first_half,
                    second_half,
                    error: trial.error,
                })
            },
        )
    }
}
