//! Error control for the complete multicomponent stirring/mixing operator.
//!
//! Each coarse/fine trial uses one donor-volume step and full Maxwell–Stefan
//! diffusion over the SAME physical interval. The single-donor CFL cap prevents
//! independent subcycling from hiding temporal error.

use super::*;
use crate::miscibility::adaptive::trial_endpoint;
use super::super::{MulticomponentMiscibilityError, finite, nonnegative, positive};
use super::super::adaptive::{
    MulticomponentMiscibilityAdvanceError, MulticomponentMiscibilityTimeControl,
    admit_time, error_ratio, factor, retryable,
};

/// Accepted physical transfers and whole-window balances. This first-order
/// local temporal estimate is not a spatial-error certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityTransportAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    /// Completed coarse/rejected trials count. An incomplete refused trial may
    /// perform additional work within the explicit attempt bound.
    pub completed_transport_steps: usize,
    pub completed_advection_substeps: usize,
    pub completed_maxwell_stefan_evaluations: usize,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_next_step_s: f64,
    pub maximum_advection_courant: f64,
    /// Canonical face then species; only accepted fine paths contribute.
    pub advected_mol: Vec<Vec<f64>>,
    pub diffused_mol: Vec<Vec<f64>>,
    pub transferred_mol: Vec<Vec<f64>>,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub free_energy_change_j: f64,
    pub advection_free_energy_change_j: f64,
    pub diffusion_free_energy_change_j: f64,
    /// Direct whole-window change minus accepted stage changes. This is
    /// numerical accounting, never heat or mechanical stirring work.
    pub stage_energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityTransportAdvanceCandidate {
    field: MulticomponentMiscibilityField,
    advance: MulticomponentMiscibilityTransportAdvance,
}
impl MulticomponentMiscibilityTransportAdvanceCandidate {
    pub const fn field(&self) -> &MulticomponentMiscibilityField { &self.field }
    pub const fn advance(&self) -> &MulticomponentMiscibilityTransportAdvance { &self.advance }
    pub fn into_parts(self) -> (MulticomponentMiscibilityField, MulticomponentMiscibilityTransportAdvance) {
        (self.field, self.advance)
    }
}

impl MulticomponentMiscibilityField {
    /// Advance frozen prescribed flow and full molecular mixing atomically.
    /// Both generators enter one step-doubling estimate, including splitting
    /// error. Independent amounts and their low words participate. The stricter
    /// original stage allowances apply again to the complete physical interval.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_transport_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection_controls: MulticomponentMiscibilityAdvectionControls,
        mut controls: MulticomponentMiscibilitySolve,
        time: MulticomponentMiscibilityTimeControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<MulticomponentMiscibilityTransportAdvanceCandidate, MulticomponentMiscibilityAdvanceError> {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        admit_time(duration_s, time)?;
        flow.admit(self.network())?;
        if !positive(advection_controls.maximum_courant)
            || advection_controls.maximum_courant > 1.0
            || !(1..=1_000_000).contains(&advection_controls.maximum_substeps)
            || !nonnegative(advection_controls.amount_balance_tolerance_mol)
            || !nonnegative(advection_controls.free_energy_tolerance_j)
        {
            return Err(MulticomponentMiscibilityError::Invalid(
                "CFL in (0,1], finite allowances and bounded advection work required",
            ).into());
        }
        controls.amount_balance_tolerance_mol = controls.amount_balance_tolerance_mol
            .min(advection_controls.amount_balance_tolerance_mol);
        controls.free_energy_tolerance_j = controls.free_energy_tolerance_j
            .min(advection_controls.free_energy_tolerance_j);
        let target = finite(self.elapsed_s + duration_s, "stirred multicomponent target time")?;
        if target <= self.elapsed_s {
            return Err(MulticomponentMiscibilityError::Invalid(
                "stirred interval cannot advance physical time",
            ).into());
        }
        let maximum_trial = if flow.maximum_outgoing_rate_s() == 0.0 {
            time.maximum_step_s
        } else {
            ((advection_controls.maximum_courant / flow.maximum_outgoing_rate_s())
                * (1.0 - 8.0 * f64::EPSILON)).min(time.maximum_step_s)
        };
        if !positive(maximum_trial)
            || (time.minimum_step_s > maximum_trial && duration_s > maximum_trial)
        {
            return Err(MulticomponentMiscibilityError::Invalid(
                "minimum timestep is incompatible with the single-donor CFL interval",
            ).into());
        }
        let advection = MulticomponentMiscibilityAdvectionControls {
            maximum_substeps: 1, ..advection_controls
        };
        let n = self.network.species_count();
        let slots = self.network.faces().len() * n;
        let mut advected = vec![Sum::default(); slots];
        let mut diffused = vec![Sum::default(); slots];
        let mut advection_energy = Sum::default();
        let mut diffusion_energy = Sum::default();
        let mut field = self.clone();
        let mut next_step = time.initial_step_s.min(maximum_trial);
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut completed_steps = 0usize;
        let mut substeps = 0usize;
        let mut evaluations = 0usize;
        let mut maximum_error = 0.0_f64;
        let mut minimum_step = f64::INFINITY;
        let mut maximum_step = 0.0_f64;
        let mut maximum_courant = 0.0_f64;
        while field.elapsed_s < target {
            poll(&mut cancelled)?;
            if attempts == time.maximum_attempts {
                return Err(MulticomponentMiscibilityAdvanceError::WorkLimit {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                    elapsed_s: field.elapsed_s - self.elapsed_s,
                });
            }
            let remaining = target - field.elapsed_s;
            let endpoint = trial_endpoint(field.elapsed_s, next_step, target);
            let dt = endpoint - field.elapsed_s;
            let floor = time.minimum_step_s.min(remaining);
            if !positive(0.5 * dt) || field.elapsed_s + 0.5 * dt <= field.elapsed_s {
                return Err(MulticomponentMiscibilityError::Invalid(
                    "stirred adaptive trial cannot represent both half steps",
                ).into());
            }
            attempts += 1;
            let trial = (|| -> Result<_> {
                let mut step = |start: &MulticomponentMiscibilityField, duration| {
                    let result = start.transport_candidate(
                        flow, duration, advection, controls, &mut cancelled,
                    )?;
                    completed_steps = count(completed_steps, 1)?;
                    substeps = count(substeps, result.step().advection.substeps)?;
                    evaluations = count(evaluations, result.step().diffusion.maxwell_stefan_evaluations)?;
                    Ok::<_, MulticomponentMiscibilityError>(result)
                };
                let coarse = step(&field, dt)?;
                let first = step(&field, 0.5 * dt)?;
                let second = step(first.field(), endpoint - first.field().elapsed_s)?;
                if coarse.field().elapsed_s != endpoint || second.field().elapsed_s != endpoint {
                    return Err(MulticomponentMiscibilityError::Invalid(
                        "stirred adaptive endpoint is unrepresentable",
                    ));
                }
                let ratio = error_ratio(coarse.field(), second.field(), time, &mut cancelled)?;
                Ok((first, second, ratio))
            })();
            let (first, second, ratio) = match trial {
                Ok(result) => result,
                Err(error) if (retryable(&error)
                    || matches!(error, MulticomponentMiscibilityError::AdvectionWorkLimit { .. }))
                    && dt > floor => {
                    rejected += 1;
                    next_step = (0.5 * dt).max(floor);
                    continue;
                }
                Err(error) => return Err(error.into()),
            };
            if ratio > 1.0 {
                rejected += 1;
                if dt <= floor {
                    return Err(MulticomponentMiscibilityAdvanceError::TemporalAccuracy {
                        elapsed_s: field.elapsed_s - self.elapsed_s,
                        duration_s: dt,
                        error_ratio: ratio,
                    });
                }
                next_step = (dt * factor(ratio)).max(floor);
                continue;
            }
            for report in [first.step(), second.step()] {
                maximum_courant = maximum_courant.max(report.advection.maximum_courant);
                advection_energy.add(report.advection.free_energy_change_j)?;
                diffusion_energy.add(report.diffusion.free_energy_change_j)?;
                for face in 0..self.network.faces().len() {
                    poll(&mut cancelled)?;
                    for species in 0..n {
                        advected[face * n + species].add(report.advection.transferred_mol[face][species])?;
                        diffused[face * n + species].add(report.diffusion.transferred_mol[face][species])?;
                    }
                }
            }
            field = second.into_parts().0;
            accepted += 1;
            maximum_error = maximum_error.max(ratio);
            minimum_step = minimum_step.min(dt);
            maximum_step = maximum_step.max(dt);
            // Even an admitted short final interval must not suggest a next
            // step above the one-donor CFL bound.
            next_step = (dt * factor(ratio)).max(time.minimum_step_s).min(maximum_trial);
        }
        let (component_residual_mol, maximum_capacity_residual_mol, free_energy_change_j) =
            self.balances(&field, controls, &mut cancelled)?;
        let mut transferred_mol = Vec::with_capacity(self.network.faces().len());
        let mut advected_mol = Vec::with_capacity(transferred_mol.capacity());
        let mut diffused_mol = Vec::with_capacity(transferred_mol.capacity());
        for (a, d) in advected.chunks(n).zip(diffused.chunks(n)) {
            poll(&mut cancelled)?;
            let mut total = Vec::with_capacity(n);
            for (left, right) in a.iter().zip(d) {
                let mut sum = *left;
                sum.add(right.high)?;
                sum.add(right.low)?;
                total.push(sum.value()?);
            }
            transferred_mol.push(total);
            advected_mol.push(a.iter().map(Sum::value).collect::<Result<Vec<_>>>()?);
            diffused_mol.push(d.iter().map(Sum::value).collect::<Result<Vec<_>>>()?);
        }
        let advection_free_energy_change_j = advection_energy.value()?;
        let diffusion_free_energy_change_j = diffusion_energy.value()?;
        let mut residual = Sum::default();
        residual.add(free_energy_change_j)?;
        residual.add(-advection_free_energy_change_j)?;
        residual.add(-diffusion_free_energy_change_j)?;
        poll(&mut cancelled)?;
        Ok(MulticomponentMiscibilityTransportAdvanceCandidate {
            field,
            advance: MulticomponentMiscibilityTransportAdvance {
                duration_s, attempted_intervals: attempts, accepted_intervals: accepted,
                rejected_intervals: rejected, completed_transport_steps: completed_steps,
                completed_advection_substeps: substeps,
                completed_maxwell_stefan_evaluations: evaluations,
                maximum_accepted_error_ratio: maximum_error,
                minimum_accepted_step_s: minimum_step,
                maximum_accepted_step_s: maximum_step,
                suggested_next_step_s: next_step,
                maximum_advection_courant: maximum_courant,
                advected_mol, diffused_mol, transferred_mol,
                component_residual_mol, maximum_capacity_residual_mol,
                free_energy_change_j, advection_free_energy_change_j,
                diffusion_free_energy_change_j, stage_energy_residual_j: residual.value()?,
            },
        })
    }
}

fn count(current: usize, additional: usize) -> Result<usize> {
    current.checked_add(additional).ok_or(MulticomponentMiscibilityError::Arithmetic(
        "stirred multicomponent work count",
    ))
}

#[cfg(test)]
mod tests;
