//! Bounded automatic time integration with a first-order step-doubling estimate.
//!
//! A full sweep is compared with two half sweeps, and the two-half result is
//! accepted. This controls a local composition estimate, not a certified global
//! error. Both independently represented species participate in the estimate.
//! The entire requested interval is atomic, including cancellation and exhausted
//! work/accuracy budgets. No state or partial flux is published on refusal.

use super::{
    BinaryMiscibilityCandidate, BinaryMiscibilityError as Error, BinaryMiscibilityNetwork,
    BinaryMiscibilitySolve, BinaryMiscibilityState, Sum, finite, nonnegative, poll, positive,
};

/// Explicit temporal accuracy, step-size bounds and maximum attempted intervals.
/// A final remainder may be smaller than the minimum; it must still pass the
/// same estimator, root, conservation and energy requirements.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryMiscibilityTimeControl {
    pub initial_step_s: f64,
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    pub absolute_mole_fraction_tolerance: f64,
    pub relative_mole_fraction_tolerance: f64,
    /// Each attempt performs at most three bounded implicit face sweeps.
    pub maximum_attempts: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityAdvanceError {
    Solver(Error),
    WorkLimit {
        attempted_intervals: usize,
        accepted_intervals: usize,
        elapsed_s: f64,
    },
    TemporalAccuracy {
        elapsed_s: f64,
        duration_s: f64,
        error_ratio: f64,
    },
}
impl From<Error> for BinaryMiscibilityAdvanceError {
    fn from(error: Error) -> Self {
        Self::Solver(error)
    }
}
impl core::fmt::Display for BinaryMiscibilityAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive binary miscibility: {self:?}")
    }
}
impl std::error::Error for BinaryMiscibilityAdvanceError {}

/// Whole-window balances and work; energy allowances are not accumulated per
/// accepted interval. The original solver's single allowance applies again to
/// the complete interval and the sum of all accepted convex-split dissipations.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    /// Includes successful coarse sweeps and sweeps later rejected by accuracy.
    pub completed_sweeps: usize,
    /// Root evaluations in completed sweeps; a refused sweep may do extra work.
    pub completed_root_iterations: u64,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_next_step_s: f64,
    /// Integrated signed A transfer in canonical face order; B is opposite.
    pub transferred_a_mol: Vec<f64>,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub free_energy_before_j: f64,
    pub free_energy_after_j: f64,
    pub free_energy_change_j: f64,
    pub convex_split_dissipation_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityAdvanceCandidate {
    state: BinaryMiscibilityState,
    advance: BinaryMiscibilityAdvance,
}
impl BinaryMiscibilityAdvanceCandidate {
    pub const fn state(&self) -> &BinaryMiscibilityState {
        &self.state
    }
    pub const fn advance(&self) -> &BinaryMiscibilityAdvance {
        &self.advance
    }
    pub fn into_parts(self) -> (BinaryMiscibilityState, BinaryMiscibilityAdvance) {
        (self.state, self.advance)
    }
}

impl BinaryMiscibilityNetwork {
    /// Advance a closed isothermal field with automatic temporal refinement.
    /// Root/energy/balance refusals can retry at a smaller step. Invalid inputs,
    /// unrepresentable arithmetic and cancellation remain immediate refusals.
    /// All accepted fine sweeps retain the original homogeneous/gradient law.
    pub fn advance_candidate(
        &self,
        initial: &BinaryMiscibilityState,
        duration_s: f64,
        time: BinaryMiscibilityTimeControl,
        solve: BinaryMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityAdvanceCandidate, BinaryMiscibilityAdvanceError> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        if !positive(duration_s)
            || !positive(time.minimum_step_s)
            || !positive(time.initial_step_s)
            || !positive(time.maximum_step_s)
            || time.minimum_step_s > time.initial_step_s
            || time.initial_step_s > time.maximum_step_s
            || !positive(time.absolute_mole_fraction_tolerance)
            || time.absolute_mole_fraction_tolerance > 1.0
            || !nonnegative(time.relative_mole_fraction_tolerance)
            || time.relative_mole_fraction_tolerance > 1.0
            || !(1..=1_000_000).contains(&time.maximum_attempts)
        {
            return Err(Error::Invalid(
                "positive bounded adaptive time/accuracy/work controls required",
            )
            .into());
        }
        let target = finite(initial.elapsed_s + duration_s, "adaptive target time")?;
        if target <= initial.elapsed_s {
            return Err(Error::Arithmetic("adaptive time interval is unrepresentable").into());
        }
        let before = self.free_energy_j(initial, &mut cancelled)?;
        let mut state = initial.clone();
        let mut next_step = time.initial_step_s;
        let mut attempts = 0;
        let mut accepted = 0;
        let mut completed_sweeps = 0;
        let mut root_iterations = 0u64;
        let mut largest_error: f64 = 0.0;
        let mut smallest_step = f64::INFINITY;
        let mut largest_step: f64 = 0.0;
        let mut transfers = vec![Sum::default(); self.owner.faces.len()];
        let mut dissipation = Sum::default();
        while state.elapsed_s < target {
            poll(&mut cancelled)?;
            if attempts == time.maximum_attempts {
                return Err(BinaryMiscibilityAdvanceError::WorkLimit {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                    elapsed_s: state.elapsed_s,
                });
            }
            attempts += 1;
            let endpoint = trial_endpoint(state.elapsed_s, next_step, target);
            let step = endpoint - state.elapsed_s;
            if !positive(step) {
                return Err(
                    Error::Arithmetic("adaptive step cannot advance represented time").into(),
                );
            }
            let trial = (|| -> super::Result<_> {
                let mut sweep = |start: &BinaryMiscibilityState, dt| {
                    let candidate = self.step_candidate(start, dt, solve, &mut cancelled)?;
                    completed_sweeps += 1;
                    root_iterations += candidate.step.root_iterations as u64;
                    Ok::<BinaryMiscibilityCandidate, Error>(candidate)
                };
                let coarse = sweep(&state, step)?;
                let half = sweep(&state, 0.5 * step)?;
                // Form the second duration from its actual represented start,
                // so a long-running field lands on the requested endpoint.
                let fine = sweep(&half.state, endpoint - half.state.elapsed_s)?;
                if coarse.state.elapsed_s != endpoint || fine.state.elapsed_s != endpoint {
                    return Err(Error::Arithmetic(
                        "adaptive endpoint cannot be represented exactly",
                    ));
                }
                Ok((coarse, half, fine))
            })();
            let (coarse, half, fine) = match trial {
                Ok(trial) => trial,
                Err(error) => {
                    let retryable = matches!(
                        error,
                        Error::Root { .. } | Error::FreeEnergy { .. } | Error::Conservation { .. }
                    );
                    let smaller = (0.5 * step).max(time.minimum_step_s);
                    if !retryable || smaller >= step {
                        return Err(error.into());
                    }
                    next_step = smaller;
                    continue;
                }
            };
            let mut error_ratio: f64 = 0.0;
            for (cell, &capacity) in self.owner.capacities_mol.iter().enumerate() {
                poll(&mut cancelled)?;
                for species in 0..2 {
                    let difference = (fine.state.amounts[cell][species]
                        - coarse.state.amounts[cell][species])
                        .abs()
                        / capacity;
                    let scale =
                        state.fractions[cell][species].max(fine.state.fractions[cell][species]);
                    let allowance = time.absolute_mole_fraction_tolerance
                        + time.relative_mole_fraction_tolerance * scale;
                    error_ratio =
                        error_ratio.max(finite(difference / allowance, "adaptive temporal error")?);
                }
            }
            let factor = if error_ratio == 0.0 {
                2.0
            } else {
                (0.8 / error_ratio.sqrt()).clamp(0.2, 2.0)
            };
            if error_ratio > 1.0 {
                let smaller = (step * factor.min(0.8)).max(time.minimum_step_s);
                if smaller >= step {
                    return Err(BinaryMiscibilityAdvanceError::TemporalAccuracy {
                        elapsed_s: state.elapsed_s,
                        duration_s: step,
                        error_ratio,
                    });
                }
                next_step = smaller;
                continue;
            }
            for report in [&half.step, &fine.step] {
                dissipation.add(report.convex_split_dissipation_j)?;
                for (total, &transfer) in transfers.iter_mut().zip(&report.transferred_a_mol) {
                    poll(&mut cancelled)?;
                    total.add(transfer)?;
                }
            }
            accepted += 1;
            largest_error = largest_error.max(error_ratio);
            smallest_step = smallest_step.min(step);
            largest_step = largest_step.max(step);
            next_step = (step * factor).clamp(time.minimum_step_s, time.maximum_step_s);
            state = fine.state;
        }
        let mut balances = [Sum::default(); 2];
        let mut capacity_residual: f64 = 0.0;
        for (cell, (&old, &new)) in initial.amounts.iter().zip(&state.amounts).enumerate() {
            poll(&mut cancelled)?;
            for species in 0..2 {
                balances[species].add(finite(
                    new[species] - old[species],
                    "adaptive component balance",
                )?)?;
            }
            capacity_residual =
                capacity_residual.max((new[0] + new[1] - self.owner.capacities_mol[cell]).abs());
        }
        let balances = [balances[0].value()?, balances[1].value()?];
        if balances
            .iter()
            .any(|x| x.abs() > solve.amount_balance_tolerance_mol)
            || capacity_residual > solve.amount_balance_tolerance_mol
        {
            return Err(Error::Conservation {
                component_residual_mol: balances,
                capacity_residual_mol: capacity_residual,
            }
            .into());
        }
        let after = self.free_energy_j(&state, &mut cancelled)?;
        let change = self.free_energy_change(initial, &state, &mut cancelled)?;
        let dissipation = dissipation.value()?;
        let defect = finite(change + dissipation, "whole-window free-energy inequality")?;
        if change > solve.free_energy_tolerance_j || defect > solve.free_energy_tolerance_j {
            return Err(Error::FreeEnergy {
                change_j: change,
                dissipation_j: dissipation,
                tolerance_j: solve.free_energy_tolerance_j,
            }
            .into());
        }
        let transfers = transfers
            .into_iter()
            .map(Sum::value)
            .collect::<super::Result<Vec<_>>>()?;
        poll(&mut cancelled)?;
        Ok(BinaryMiscibilityAdvanceCandidate {
            advance: BinaryMiscibilityAdvance {
                duration_s: target - initial.elapsed_s,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: attempts - accepted,
                completed_sweeps,
                completed_root_iterations: root_iterations,
                maximum_accepted_error_ratio: largest_error,
                minimum_accepted_step_s: smallest_step,
                maximum_accepted_step_s: largest_step,
                suggested_next_step_s: next_step,
                transferred_a_mol: transfers,
                component_residual_mol: balances,
                maximum_capacity_residual_mol: capacity_residual,
                free_energy_before_j: before,
                free_energy_after_j: after,
                free_energy_change_j: change,
                convex_split_dissipation_j: dissipation,
            },
            state,
        })
    }
}

/// Endpoint of the next step-doubling trial: `elapsed_s + step_s`, or
/// `target_s` once the step reaches it. Each trial is compared with two half
/// steps, so the interval it leaves behind must itself be divisible. When the
/// midpoint of the leftover `[endpoint, target_s]` is not strictly inside it
/// (a remainder of one ulp, left by rounded sums of steps), no later trial
/// could take it, and this trial runs on to `target_s` instead. Every other
/// endpoint is the one the step alone gives.
pub(crate) fn trial_endpoint(elapsed_s: f64, step_s: f64, target_s: f64) -> f64 {
    let endpoint = if step_s >= target_s - elapsed_s {
        target_s
    } else {
        (elapsed_s + step_s).min(target_s)
    };
    let leftover_midpoint = endpoint + 0.5 * (target_s - endpoint);
    let divisible = endpoint < leftover_midpoint && leftover_midpoint < target_s;
    if endpoint < target_s && !divisible {
        target_s
    } else {
        endpoint
    }
}

#[cfg(test)]
#[path = "adaptive/tests.rs"]
mod tests;
