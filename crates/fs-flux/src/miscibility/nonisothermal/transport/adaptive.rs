//! Adaptive complete thermal stirring and molecular/heat diffusion.
//!
//! One full Lie step is compared with two half steps using BOTH original
//! component inventories and local temperatures. Each trial contains ONE donor
//! advection substep: allowing independent CFL subcycling would sometimes make
//! coarse and fine paths identical and silently erase their time-error estimate.
//! Only the accepted two-half-step path contributes to physical face receipts.
//! Original whole-window energy, entropy and amount budgets apply again before
//! publication. This estimates local temporal error, not a certified global or
//! spatial error bound; prescribed flow is frozen throughout the interval.

use super::*;
use crate::miscibility::adaptive::trial_endpoint;
use crate::miscibility::nonisothermal::adaptive::{
    BinaryNonisothermalAdvanceError, BinaryNonisothermalTimeControl, admit_time, error_ratio,
    factor, retryable,
};

/// Whole-window physical exchanges and bounded completed-trial work.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalTransportAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    /// Coarse and rejected trials count if they completed both operators.
    pub completed_transport_steps: usize,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    /// Capped at the single-donor CFL interval, including short final remainders.
    pub suggested_next_step_s: f64,
    pub maximum_advection_courant: f64,
    /// Work in completed trials; a refused trial may perform additional work.
    pub root_iterations: u64,
    pub caloric_iterations: u64,
    /// Integrated accepted transfers, oriented by the original canonical faces.
    pub transferred_component_moles: Vec<[f64; 2]>,
    pub transferred_energy_j: Vec<f64>,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

/// The full requested physical interval succeeds or no field is published.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalTransportAdvanceCandidate {
    field: BinaryNonisothermalField,
    advance: BinaryNonisothermalTransportAdvance,
}
impl BinaryNonisothermalTransportAdvanceCandidate {
    pub const fn field(&self) -> &BinaryNonisothermalField {
        &self.field
    }
    pub const fn advance(&self) -> &BinaryNonisothermalTransportAdvance {
        &self.advance
    }
    pub fn into_parts(self) -> (BinaryNonisothermalField, BinaryNonisothermalTransportAdvance) {
        (self.field, self.advance)
    }
}

impl BinaryNonisothermalField {
    /// Adapt the COMPLETE stirred heat/species operator, accepting two half
    /// steps with independent composition and temperature tolerances. Internal
    /// CFL subcycling cannot hide an advective temporal error. An incompatible
    /// minimum step, exhausted work, original-source refusal or cancellation
    /// leaves the original field unchanged, including every material binding.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_transport_candidate(
        &self,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        advection_controls: BinaryMiscibilityAdvectionControls,
        mut controls: BinaryNonisothermalSolve,
        time: BinaryNonisothermalTimeControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<
        BinaryNonisothermalTransportAdvanceCandidate,
        BinaryNonisothermalAdvanceError,
    > {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        admit_time(duration_s, time)?;
        flow.admit(self.network())
            .map_err(BinaryNonisothermalError::from)?;
        if !positive(advection_controls.maximum_courant)
            || advection_controls.maximum_courant > 1.0
            || !(1..=1_000_000).contains(&advection_controls.maximum_substeps)
            || !nonnegative(advection_controls.amount_balance_tolerance_mol)
        {
            return Err(BinaryNonisothermalError::Invalid(
                "CFL in (0,1], finite amount allowance and bounded advection work required",
            )
            .into());
        }
        controls.amount_balance_tolerance_mol = controls
            .amount_balance_tolerance_mol
            .min(advection_controls.amount_balance_tolerance_mol);
        let t = time.composition;
        let target = finite(self.elapsed_s() + duration_s, "coupled thermal target time")
            .map_err(BinaryNonisothermalError::from)?;
        if target <= self.elapsed_s() {
            return Err(BinaryNonisothermalError::Invalid(
                "coupled thermal interval cannot advance physical time",
            )
            .into());
        }
        let maximum_trial = if flow.maximum_outgoing_rate_s() == 0.0 {
            t.maximum_step_s
        } else {
            ((advection_controls.maximum_courant / flow.maximum_outgoing_rate_s())
                * (1.0 - 8.0 * f64::EPSILON))
                .min(t.maximum_step_s)
        };
        if !positive(maximum_trial) {
            return Err(BinaryNonisothermalError::Invalid(
                "single-donor thermal CFL interval is unrepresentable",
            )
            .into());
        }
        let advection = BinaryMiscibilityAdvectionControls {
            maximum_substeps: 1,
            ..advection_controls
        };
        let mut field = self.clone();
        let mut next_step = t.initial_step_s
            .min(maximum_trial)
            .max(t.minimum_step_s);
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut completed_steps = 0usize;
        let mut root_iterations = 0u64;
        let mut caloric_iterations = 0u64;
        let mut largest_error = 0.0_f64;
        let mut smallest_step = f64::INFINITY;
        let mut largest_step = 0.0_f64;
        let mut maximum_courant = 0.0_f64;
        let mut amounts = vec![[Sum::default(); 2]; self.owner.faces.len()];
        let mut energy = vec![Sum::default(); self.owner.faces.len()];
        while field.elapsed_s() < target {
            poll(&mut cancelled)?;
            if attempts == t.maximum_attempts {
                return Err(BinaryNonisothermalAdvanceError::WorkLimit {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                    elapsed_s: field.elapsed_s() - self.elapsed_s(),
                });
            }
            let endpoint = trial_endpoint(field.elapsed_s(), next_step, target);
            let dt = endpoint - field.elapsed_s();
            if !positive(dt) {
                return Err(BinaryNonisothermalError::Invalid(
                    "coupled thermal step cannot advance physical time",
                )
                .into());
            }
            attempts += 1;
            let trial = (|| -> Result<_> {
                let mut advance = |start: &BinaryNonisothermalField, duration| {
                    let result = start.transport_candidate(
                        flow,
                        duration,
                        advection,
                        controls,
                        &mut cancelled,
                    )?;
                    completed_steps += 1;
                    root_iterations += result.step().diffusion.root_iterations;
                    caloric_iterations += result.step().diffusion.caloric_iterations
                        + result.step().advection.caloric_iterations;
                    Ok::<_, BinaryNonisothermalError>(result)
                };
                let coarse = advance(&field, dt)?;
                let first = advance(&field, 0.5 * dt)?;
                let second = advance(first.field(), endpoint - first.field().elapsed_s())?;
                if coarse.field().elapsed_s() != endpoint
                    || second.field().elapsed_s() != endpoint
                {
                    return Err(BinaryNonisothermalError::Invalid(
                        "coupled thermal endpoint is unrepresentable",
                    ));
                }
                let ratio = error_ratio(coarse.field(), second.field(), time, &mut cancelled)?;
                Ok((first, second, ratio))
            })();
            let (first, second, ratio) = match trial {
                Ok(result) => result,
                Err(error) => {
                    let may_refine = retryable(&error)
                        || matches!(
                            error,
                            BinaryNonisothermalError::Species { .. }
                                | BinaryNonisothermalError::Spatial(
                                    BinaryMiscibilityError::AdvectionWorkLimit { .. }
                                )
                        );
                    let smaller = (0.5 * dt).max(t.minimum_step_s);
                    if !may_refine || smaller >= dt {
                        return Err(error.into());
                    }
                    next_step = smaller;
                    continue;
                }
            };
            if ratio > 1.0 {
                let smaller = (dt * factor(ratio)).max(t.minimum_step_s);
                if smaller >= dt {
                    return Err(BinaryNonisothermalAdvanceError::TemporalAccuracy {
                        elapsed_s: field.elapsed_s() - self.elapsed_s(),
                        duration_s: dt,
                        error_ratio: ratio,
                    });
                }
                next_step = smaller;
                continue;
            }
            for report in [first.step(), second.step()] {
                maximum_courant = maximum_courant.max(report.advection.maximum_courant);
                for face in 0..amounts.len() {
                    poll(&mut cancelled)?;
                    for species in 0..2 {
                        amounts[face][species]
                            .add(report.transferred_component_moles[face][species])
                            .map_err(BinaryNonisothermalError::from)?;
                    }
                    energy[face]
                        .add(report.transferred_energy_j[face])
                        .map_err(BinaryNonisothermalError::from)?;
                }
            }
            field = second.into_parts().0;
            accepted += 1;
            largest_error = largest_error.max(ratio);
            smallest_step = smallest_step.min(dt);
            largest_step = largest_step.max(dt);
            next_step = (dt * factor(ratio))
                .min(maximum_trial)
                .max(t.minimum_step_s);
        }
        let (
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
        ) = self.balances(&field, controls, &mut cancelled)?;
        let mut transferred_component_moles = Vec::with_capacity(amounts.len());
        let mut transferred_energy_j = Vec::with_capacity(energy.len());
        for face in 0..amounts.len() {
            poll(&mut cancelled)?;
            transferred_component_moles.push([
                amounts[face][0].value().map_err(BinaryNonisothermalError::from)?,
                amounts[face][1].value().map_err(BinaryNonisothermalError::from)?,
            ]);
            transferred_energy_j.push(
                energy[face].value().map_err(BinaryNonisothermalError::from)?,
            );
        }
        poll(&mut cancelled)?;
        Ok(BinaryNonisothermalTransportAdvanceCandidate {
            field,
            advance: BinaryNonisothermalTransportAdvance {
                duration_s: target - self.elapsed_s(),
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: attempts - accepted,
                completed_transport_steps: completed_steps,
                maximum_accepted_error_ratio: largest_error,
                minimum_accepted_step_s: smallest_step,
                maximum_accepted_step_s: largest_step,
                suggested_next_step_s: next_step.min(maximum_trial),
                maximum_advection_courant: maximum_courant,
                root_iterations,
                caloric_iterations,
                transferred_component_moles,
                transferred_energy_j,
                component_residual_mol,
                maximum_capacity_residual_mol,
                energy_residual_j,
                entropy_change_j_per_k,
            },
        })
    }
}

#[cfg(test)]
mod tests;
