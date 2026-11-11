//! Adaptive prescribed stirring and demixing on the original closed graph.
//!
//! Compare a complete Lie step with two half steps and accept the latter. Each
//! trial uses ONE donor advection step: otherwise CFL subcycling can make the
//! coarse and fine paths identical, falsely reporting zero temporal error.
//! The estimate controls the coupled local composition error, not spatial
//! numerical mixing or a certified global error. The solution is first order.
//!
//! Only accepted stages contribute to integrated flux and energy accounts.
//! Stirring can increase free energy; diffusion retains its own inequality.
//! The original amount and diffusion-energy allowances apply once again to the
//! entire requested interval. No failed interval publishes a partial result.

use super::{BinaryMiscibilityTransportCandidate, BinaryMiscibilityTransportStep};
use crate::miscibility::adaptive::{
    BinaryMiscibilityAdvanceError, BinaryMiscibilityTimeControl, trial_endpoint,
};
use crate::miscibility::advection::{
    BinaryMiscibilityAdvection, BinaryMiscibilityAdvectionControls,
};
use crate::miscibility::{
    BinaryMiscibilityError as Error, BinaryMiscibilityNetwork, BinaryMiscibilitySolve,
    BinaryMiscibilityState, Result as SolverResult, Sum, finite, nonnegative, poll, positive,
};

/// Whole-interval physical transfers, energy changes and bounded solver work.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityTransportAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    /// Includes coarse and rejected trials that completed both operators.
    pub completed_transport_steps: usize,
    /// Work in completed transport steps; a refused step may do extra work.
    pub completed_advection_substeps: u64,
    pub completed_root_iterations: u64,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    /// CFL-capped suggestion. May be below the supplied minimum when only a
    /// shorter final remainder was admissible under that minimum.
    pub suggested_next_step_s: f64,
    /// Largest outgoing-volume Courant number of accepted advection stages.
    pub maximum_advection_courant: f64,
    /// Signed A/B transfers, in canonical face order, from accepted advection.
    pub advected_component_moles: Vec<[f64; 2]>,
    /// Signed A transfers from accepted diffusion; B transfers are opposite.
    pub diffused_a_mol: Vec<f64>,
    /// Sum of the two original operator transfer accounts in canonical order.
    pub transferred_component_moles: Vec<[f64; 2]>,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub free_energy_before_j: f64,
    pub free_energy_after_j: f64,
    /// Direct original-to-final increment, without subtracting bulk totals.
    pub free_energy_change_j: f64,
    /// Either sign is allowed. This is not a heat or mechanical-work account.
    pub advection_free_energy_change_j: f64,
    pub diffusion_free_energy_change_j: f64,
    pub diffusion_convex_split_dissipation_j: f64,
    /// Direct total increment minus the sum of accepted operator increments.
    pub stage_energy_residual_j: f64,
}

/// The complete requested interval succeeds or the original field is unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityTransportAdvanceCandidate {
    state: BinaryMiscibilityState,
    advance: BinaryMiscibilityTransportAdvance,
}
impl BinaryMiscibilityTransportAdvanceCandidate {
    pub const fn state(&self) -> &BinaryMiscibilityState {
        &self.state
    }
    pub const fn advance(&self) -> &BinaryMiscibilityTransportAdvance {
        &self.advance
    }
    pub fn into_parts(self) -> (BinaryMiscibilityState, BinaryMiscibilityTransportAdvance) {
        (self.state, self.advance)
    }
}

impl BinaryMiscibilityNetwork {
    /// Advance the full prescribed-flow/molecular evolution with local temporal
    /// error control on both species, preserving the original fixed volume frame.
    ///
    /// Each trial contains one donor step followed by one implicit face sweep.
    /// CFL limits therefore restrict the macro interval, even when the supplied
    /// advection controls allow subcycling. If the minimum time step prevents
    /// this restriction, the original advection work refusal is returned. A final
    /// remainder below the minimum remains allowed, with the same error checks.
    /// The prepared flow is frozen throughout the interval.
    ///
    /// Root, conservation, diffusion-energy and donor-CFL refusals may retry at
    /// smaller steps. Invalid inputs, unrepresentable arithmetic, cancellation
    /// and exhausted work/accuracy bounds publish no candidate. Work is at most
    /// three bounded transport trials per attempt; no time history is retained.
    /// O(cells + faces) payload fits the original graph's workspace admission.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_transport_candidate(
        &self,
        initial: &BinaryMiscibilityState,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        advection_controls: BinaryMiscibilityAdvectionControls,
        diffusion_controls: BinaryMiscibilitySolve,
        time: BinaryMiscibilityTimeControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityTransportAdvanceCandidate, BinaryMiscibilityAdvanceError> {
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
        // Validate the original controls before narrowing only their work cap.
        if !positive(advection_controls.maximum_courant)
            || advection_controls.maximum_courant > 1.0
            || !(1..=1_000_000).contains(&advection_controls.maximum_substeps)
            || !nonnegative(advection_controls.amount_balance_tolerance_mol)
        {
            return Err(Error::Invalid(
                "CFL in (0,1], bounded substeps and finite conservation allowance required",
            )
            .into());
        }
        let target = finite(
            initial.elapsed_s + duration_s,
            "adaptive transport target time",
        )?;
        if target <= initial.elapsed_s {
            return Err(Error::Arithmetic("adaptive transport interval is unrepresentable").into());
        }
        let maximum_trial = if flow.maximum_outgoing_rate_s() == 0.0 {
            time.maximum_step_s
        } else {
            // Leave a small roundoff margin at the donor CFL ceiling. The step
            // owner independently enforces one substep, including tiny scales.
            ((advection_controls.maximum_courant / flow.maximum_outgoing_rate_s())
                * (1.0 - 8.0 * f64::EPSILON))
                .min(time.maximum_step_s)
        };
        if !positive(maximum_trial) {
            return Err(Error::Arithmetic("single donor CFL interval is unrepresentable").into());
        }
        let advection = BinaryMiscibilityAdvectionControls {
            maximum_substeps: 1,
            ..advection_controls
        };
        let before = self.free_energy_j(initial, &mut cancelled)?;
        let mut state = initial.clone();
        let mut next_step = time
            .initial_step_s
            .min(maximum_trial)
            .max(time.minimum_step_s);
        let mut attempts = 0;
        let mut accepted = 0;
        let mut completed_steps = 0;
        let mut advection_substeps = 0u64;
        let mut root_iterations = 0u64;
        let mut largest_error: f64 = 0.0;
        let mut smallest_step = f64::INFINITY;
        let mut largest_step: f64 = 0.0;
        let mut maximum_courant: f64 = 0.0;
        let mut advected = vec![[Sum::default(); 2]; self.owner.faces.len()];
        let mut diffused = vec![Sum::default(); self.owner.faces.len()];
        let mut advection_energy = Sum::default();
        let mut diffusion_energy = Sum::default();
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
                return Err(Error::Arithmetic("adaptive transport cannot advance time").into());
            }
            let trial = (|| -> SolverResult<_> {
                let mut advance = |start: &BinaryMiscibilityState, dt| {
                    let result = self.transport_candidate(
                        start,
                        flow,
                        dt,
                        advection,
                        diffusion_controls,
                        &mut cancelled,
                    )?;
                    completed_steps += 1;
                    advection_substeps += result.step.advection.substeps as u64;
                    root_iterations += result.step.diffusion.root_iterations as u64;
                    Ok::<BinaryMiscibilityTransportCandidate, Error>(result)
                };
                let coarse = advance(&state, step)?;
                let half = advance(&state, 0.5 * step)?;
                let fine = advance(&half.state, endpoint - half.state.elapsed_s)?;
                if coarse.state.elapsed_s != endpoint || fine.state.elapsed_s != endpoint {
                    return Err(Error::Arithmetic(
                        "adaptive transport endpoint is unrepresentable",
                    ));
                }
                Ok((coarse, half, fine))
            })();
            let (coarse, half, fine) = match trial {
                Ok(trial) => trial,
                Err(error) => {
                    let retryable = matches!(
                        error,
                        Error::Root { .. }
                            | Error::FreeEnergy { .. }
                            | Error::Conservation { .. }
                            | Error::AdvectionWorkLimit { .. }
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
                    let mut difference = Sum::default();
                    difference.add(fine.state.amounts[cell][species])?;
                    difference.add(-coarse.state.amounts[cell][species])?;
                    difference.add(fine.state.amount_remainders[cell][species])?;
                    difference.add(-coarse.state.amount_remainders[cell][species])?;
                    let scale =
                        state.fractions[cell][species].max(fine.state.fractions[cell][species]);
                    let allowance = time.absolute_mole_fraction_tolerance
                        + time.relative_mole_fraction_tolerance * scale;
                    error_ratio = error_ratio.max(finite(
                        (difference.value()?.abs() / capacity) / allowance,
                        "adaptive transport temporal error",
                    )?);
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
                advection_energy.add(report.advection.free_energy_change_j)?;
                diffusion_energy.add(report.diffusion.free_energy_change_j)?;
                dissipation.add(report.diffusion.convex_split_dissipation_j)?;
                maximum_courant = maximum_courant.max(report.advection.maximum_courant);
                accumulate_transfers(report, &mut advected, &mut diffused, &mut cancelled)?;
            }
            accepted += 1;
            largest_error = largest_error.max(error_ratio);
            smallest_step = smallest_step.min(step);
            largest_step = largest_step.max(step);
            next_step = (step * factor).min(maximum_trial).max(time.minimum_step_s);
            state = fine.state;
        }
        let allowance = advection_controls
            .amount_balance_tolerance_mol
            .min(diffusion_controls.amount_balance_tolerance_mol);
        let mut balances = [Sum::default(); 2];
        let mut capacity_residual: f64 = 0.0;
        for (cell, (old, new)) in initial.amounts.iter().zip(&state.amounts).enumerate() {
            poll(&mut cancelled)?;
            for species in 0..2 {
                balances[species].add(finite(new[species] - old[species], "transport balance")?)?;
            }
            capacity_residual = capacity_residual.max(finite(
                (new[0] + new[1] - self.owner.capacities_mol[cell]).abs(),
                "transport capacity residual",
            )?);
        }
        let balances = [balances[0].value()?, balances[1].value()?];
        if balances.iter().any(|x| x.abs() > allowance) || capacity_residual > allowance {
            return Err(Error::Conservation {
                component_residual_mol: balances,
                capacity_residual_mol: capacity_residual,
            }
            .into());
        }
        let after = self.free_energy_j(&state, &mut cancelled)?;
        let change = self.free_energy_change(initial, &state, &mut cancelled)?;
        let advection_energy = advection_energy.value()?;
        let diffusion_energy = diffusion_energy.value()?;
        let dissipation = dissipation.value()?;
        let mut adjusted = Sum::default();
        adjusted.add(change)?;
        adjusted.add(-advection_energy)?;
        let adjusted_change = adjusted.value()?;
        let mut residual = adjusted;
        residual.add(-diffusion_energy)?;
        let stage_energy_residual_j = residual.value()?;
        // Audit accepted diffusion and the independently reconstructed total
        // after accounting for imposed advection. Never demand total F decrease.
        for diffusion_change in [diffusion_energy, adjusted_change] {
            let defect = finite(
                diffusion_change + dissipation,
                "transport diffusion inequality",
            )?;
            if diffusion_change > diffusion_controls.free_energy_tolerance_j
                || defect > diffusion_controls.free_energy_tolerance_j
            {
                return Err(Error::FreeEnergy {
                    change_j: diffusion_change,
                    dissipation_j: dissipation,
                    tolerance_j: diffusion_controls.free_energy_tolerance_j,
                }
                .into());
            }
        }
        let mut advected_component_moles = Vec::with_capacity(advected.len());
        let mut diffused_a_mol = Vec::with_capacity(diffused.len());
        let mut transferred_component_moles = Vec::with_capacity(diffused.len());
        for (advection, diffusion) in advected.into_iter().zip(diffused) {
            poll(&mut cancelled)?;
            let a = [advection[0].value()?, advection[1].value()?];
            let d = diffusion.value()?;
            let mut total = advection;
            total[0].add(diffusion.main)?;
            total[0].add(diffusion.correction)?;
            total[1].add(-diffusion.main)?;
            total[1].add(-diffusion.correction)?;
            advected_component_moles.push(a);
            diffused_a_mol.push(d);
            transferred_component_moles.push([total[0].value()?, total[1].value()?]);
        }
        poll(&mut cancelled)?;
        Ok(BinaryMiscibilityTransportAdvanceCandidate {
            state,
            advance: BinaryMiscibilityTransportAdvance {
                duration_s: target - initial.elapsed_s,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: attempts - accepted,
                completed_transport_steps: completed_steps,
                completed_advection_substeps: advection_substeps,
                completed_root_iterations: root_iterations,
                maximum_accepted_error_ratio: largest_error,
                minimum_accepted_step_s: smallest_step,
                maximum_accepted_step_s: largest_step,
                suggested_next_step_s: next_step.min(maximum_trial),
                maximum_advection_courant: maximum_courant,
                advected_component_moles,
                diffused_a_mol,
                transferred_component_moles,
                component_residual_mol: balances,
                maximum_capacity_residual_mol: capacity_residual,
                free_energy_before_j: before,
                free_energy_after_j: after,
                free_energy_change_j: change,
                advection_free_energy_change_j: advection_energy,
                diffusion_free_energy_change_j: diffusion_energy,
                diffusion_convex_split_dissipation_j: dissipation,
                stage_energy_residual_j,
            },
        })
    }
}

fn accumulate_transfers(
    report: &BinaryMiscibilityTransportStep,
    advected: &mut [[Sum; 2]],
    diffused: &mut [Sum],
    cancelled: &mut impl FnMut() -> bool,
) -> SolverResult<()> {
    for (face, (advection, diffusion)) in advected.iter_mut().zip(diffused).enumerate() {
        poll(cancelled)?;
        for (species, amount) in advection.iter_mut().enumerate() {
            amount.add(report.advection.transferred_component_moles[face][species])?;
        }
        diffusion.add(report.diffusion.transferred_a_mol[face])?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
