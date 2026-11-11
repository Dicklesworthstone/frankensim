//! Temporal accuracy for complete stirring, molecular mixing and heat exchange.
//!
//! Each full/half trial contains both physical generators. One donor-CFL step
//! per trial prevents independently subcycled coarse/fine advection from
//! becoming the same trajectory and hiding its temporal error.

use super::super::adaptive::{admit_time, error_ratio, retryable};
use super::super::advection::MulticomponentNonisothermalAdvectionControls;
use super::super::*;
use crate::miscibility::adaptive::trial_endpoint;
use crate::miscibility::multicomponent::adaptive::factor;
use crate::miscibility::multicomponent::advection::MulticomponentMiscibilityAdvection;

pub use super::super::adaptive::{
    MulticomponentNonisothermalAdvanceError as MulticomponentNonisothermalTransportAdvanceError,
    MulticomponentNonisothermalTimeControl,
};

/// Accepted face receipts and independently checked complete-window balances.
/// Completed work counts include rejected trials; an incomplete refused stage
/// may do additional work within the explicit attempt/root/CFL limits.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalTransportAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub completed_transport_steps: usize,
    pub completed_advection_substeps: usize,
    pub completed_maxwell_stefan_evaluations: usize,
    pub completed_caloric_iterations: u64,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_next_step_s: f64,
    pub maximum_advection_courant: f64,
    /// Accumulated accepted donor resolution estimate, charged once against
    /// the complete window's energy allowance rather than once per interval.
    pub parcel_energy_resolution_bound_j: f64,
    pub advected_mol: Vec<Vec<f64>>,
    pub diffused_mol: Vec<Vec<f64>>,
    pub transferred_mol: Vec<Vec<f64>>,
    pub advected_energy_j: Vec<f64>,
    pub diffused_energy_j: Vec<f64>,
    pub transferred_energy_j: Vec<f64>,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalTransportAdvanceCandidate {
    field: MulticomponentNonisothermalField,
    advance: MulticomponentNonisothermalTransportAdvance,
}
impl MulticomponentNonisothermalTransportAdvanceCandidate {
    pub const fn field(&self) -> &MulticomponentNonisothermalField {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentNonisothermalTransportAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentNonisothermalField,
        MulticomponentNonisothermalTransportAdvance,
    ) {
        (self.field, self.advance)
    }
}

impl MulticomponentNonisothermalField {
    /// Advance one complete physical window with composition AND temperature
    /// error control, including the error of the operator splitting. The
    /// prescribed original flow and kinetic face data stay fixed; caloric
    /// recovery uses each actual stage's composition and local temperature.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_transport_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection_controls: MulticomponentNonisothermalAdvectionControls,
        mut controls: MulticomponentNonisothermalSolve,
        time: MulticomponentNonisothermalTimeControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<
        MulticomponentNonisothermalTransportAdvanceCandidate,
        MulticomponentNonisothermalTransportAdvanceError,
    > {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        admit_time(duration_s, time)?;
        flow.admit(self.network())?;
        if !positive(advection_controls.maximum_courant)
            || advection_controls.maximum_courant > 1.0
            || !(1..=1_000_000).contains(&advection_controls.maximum_substeps)
            || !nonnegative(advection_controls.amount_balance_tolerance_mol)
        {
            return Err(MulticomponentNonisothermalError::Invalid(
                "CFL in (0,1], bounded donor work and finite inventory allowance required",
            )
            .into());
        }
        controls.amount_balance_tolerance_mol = controls
            .amount_balance_tolerance_mol
            .min(advection_controls.amount_balance_tolerance_mol);
        let target = finite(
            self.elapsed_s() + duration_s,
            "thermal stirring target time",
        )?;
        if target <= self.elapsed_s() {
            return Err(MulticomponentNonisothermalError::Invalid(
                "stirred physical interval cannot advance time",
            )
            .into());
        }
        let maximum_trial = if flow.maximum_outgoing_rate_s() == 0.0 {
            time.composition.maximum_step_s
        } else {
            ((advection_controls.maximum_courant / flow.maximum_outgoing_rate_s())
                * (1.0 - 8.0 * f64::EPSILON))
                .min(time.composition.maximum_step_s)
        };
        if !positive(maximum_trial)
            || (time.composition.minimum_step_s > maximum_trial && duration_s > maximum_trial)
        {
            return Err(MulticomponentNonisothermalError::Invalid(
                "minimum thermal timestep conflicts with the single-donor CFL bound",
            )
            .into());
        }
        let advection = MulticomponentNonisothermalAdvectionControls {
            maximum_substeps: 1,
            ..advection_controls
        };
        let n = self.network().species_count();
        let faces = self.network().faces().len();
        let slots = faces * n;
        let mut advected = vec![Sum::default(); slots];
        let mut diffused = vec![Sum::default(); slots];
        let mut advected_heat = vec![Sum::default(); faces];
        let mut diffused_heat = vec![Sum::default(); faces];
        let mut field = self.clone();
        let mut next_step = time.composition.initial_step_s.min(maximum_trial);
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut completed_steps = 0usize;
        let mut substeps = 0usize;
        let mut evaluations = 0usize;
        let mut caloric_iterations = 0u64;
        let mut maximum_error = 0.0_f64;
        let mut minimum_step = f64::INFINITY;
        let mut maximum_step = 0.0_f64;
        let mut maximum_courant = 0.0_f64;
        let mut parcel_resolution = 0.0_f64;
        while field.elapsed_s() < target {
            poll(&mut cancelled)?;
            if attempts == time.composition.maximum_attempts {
                return Err(
                    MulticomponentNonisothermalTransportAdvanceError::WorkLimit {
                        attempted_intervals: attempts,
                        accepted_intervals: accepted,
                        elapsed_s: field.elapsed_s() - self.elapsed_s(),
                    },
                );
            }
            let remaining = target - field.elapsed_s();
            let endpoint = trial_endpoint(field.elapsed_s(), next_step, target);
            let dt = endpoint - field.elapsed_s();
            let floor = time.composition.minimum_step_s.min(remaining);
            if !positive(0.5 * dt) || field.elapsed_s() + 0.5 * dt <= field.elapsed_s() {
                return Err(MulticomponentNonisothermalError::Invalid(
                    "stirred thermal trial cannot represent both half intervals",
                )
                .into());
            }
            attempts += 1;
            let trial = (|| -> Result<_> {
                let mut step = |start: &MulticomponentNonisothermalField, duration| -> Result<_> {
                    let result = start.transport_candidate(
                        flow,
                        duration,
                        advection,
                        controls,
                        &mut cancelled,
                    )?;
                    completed_steps = count(completed_steps, 1)?;
                    substeps = count(substeps, result.step().advection.substeps)?;
                    evaluations = count(
                        evaluations,
                        result.step().diffusion.maxwell_stefan_evaluations,
                    )?;
                    caloric_iterations = caloric_iterations
                        .checked_add(result.step().advection.caloric_iterations)
                        .and_then(|v| v.checked_add(result.step().diffusion.caloric_iterations))
                        .ok_or(MulticomponentMiscibilityError::Arithmetic(
                            "coupled thermal caloric-work count",
                        ))?;
                    Ok(result)
                };
                let coarse = step(&field, dt)?;
                let first = step(&field, 0.5 * dt)?;
                let second = step(first.field(), endpoint - first.field().elapsed_s())?;
                if coarse.field().elapsed_s() != endpoint || second.field().elapsed_s() != endpoint
                {
                    return Err(MulticomponentNonisothermalError::Invalid(
                        "stirred thermal endpoint is unrepresentable",
                    ));
                }
                let ratio = error_ratio(coarse.field(), second.field(), time, &mut cancelled)?;
                Ok((first, second, ratio))
            })();
            let (first, second, ratio) = match trial {
                Ok(result) => result,
                Err(error)
                    if (retryable(&error)
                        || matches!(
                            error,
                            MulticomponentNonisothermalError::Spatial(
                                MulticomponentMiscibilityError::AdvectionWorkLimit { .. }
                            )
                        ))
                        && dt > floor =>
                {
                    rejected += 1;
                    next_step = (0.5 * dt).max(floor);
                    continue;
                }
                Err(error) => return Err(error.into()),
            };
            if ratio > 1.0 {
                rejected += 1;
                if dt <= floor {
                    return Err(
                        MulticomponentNonisothermalTransportAdvanceError::TemporalAccuracy {
                            elapsed_s: field.elapsed_s() - self.elapsed_s(),
                            duration_s: dt,
                            error_ratio: ratio,
                        },
                    );
                }
                next_step = (dt * factor(ratio)).max(floor);
                continue;
            }
            for report in [first.step(), second.step()] {
                maximum_courant = maximum_courant.max(report.advection.maximum_courant);
                if report.advection.parcel_energy_resolution_bound_j > 0.0 {
                    parcel_resolution = (parcel_resolution
                        + report.advection.parcel_energy_resolution_bound_j)
                        .next_up();
                }
                if parcel_resolution > 0.25 * controls.energy_tolerance_j {
                    return Err(MulticomponentNonisothermalError::UnresolvedParcelEnthalpy {
                        bound_j: parcel_resolution,
                        tolerance_j: 0.25 * controls.energy_tolerance_j,
                    }
                    .into());
                }
                for face in 0..faces {
                    poll(&mut cancelled)?;
                    for species in 0..n {
                        advected[face * n + species]
                            .add(report.advection.transferred_mol[face][species])?;
                        diffused[face * n + species]
                            .add(report.diffusion.transferred_mol[face][species])?;
                    }
                    advected_heat[face].add(report.advection.transferred_energy_j[face])?;
                    diffused_heat[face].add(report.diffusion.transferred_energy_j[face])?;
                }
            }
            field = second.into_parts().0;
            accepted += 1;
            maximum_error = maximum_error.max(ratio);
            minimum_step = minimum_step.min(dt);
            maximum_step = maximum_step.max(dt);
            // A short terminal remainder must also leave a CFL-valid next-step
            // suggestion when the user's ordinary minimum is larger than it.
            next_step = (dt * factor(ratio))
                .max(time.composition.minimum_step_s)
                .min(maximum_trial);
        }
        let (
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
        ) = self.balances(&field, controls, &mut cancelled)?;
        let mut advected_mol = Vec::with_capacity(faces);
        let mut diffused_mol = Vec::with_capacity(faces);
        let mut transferred_mol = Vec::with_capacity(faces);
        for (a, d) in advected.chunks(n).zip(diffused.chunks(n)) {
            poll(&mut cancelled)?;
            transferred_mol.push(
                a.iter()
                    .zip(d)
                    .map(|(a, d)| total(*a, *d))
                    .collect::<Result<Vec<_>>>()?,
            );
            advected_mol.push(
                a.iter()
                    .map(Sum::value)
                    .collect::<crate::miscibility::multicomponent::Result<Vec<_>>>()?,
            );
            diffused_mol.push(
                d.iter()
                    .map(Sum::value)
                    .collect::<crate::miscibility::multicomponent::Result<Vec<_>>>()?,
            );
        }
        let transferred_energy_j = advected_heat
            .iter()
            .zip(&diffused_heat)
            .map(|(a, d)| total(*a, *d))
            .collect::<Result<Vec<_>>>()?;
        let advected_energy_j = advected_heat
            .iter()
            .map(Sum::value)
            .collect::<crate::miscibility::multicomponent::Result<Vec<_>>>()?;
        let diffused_energy_j = diffused_heat
            .iter()
            .map(Sum::value)
            .collect::<crate::miscibility::multicomponent::Result<Vec<_>>>()?;
        poll(&mut cancelled)?;
        Ok(MulticomponentNonisothermalTransportAdvanceCandidate {
            field,
            advance: MulticomponentNonisothermalTransportAdvance {
                duration_s,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                completed_transport_steps: completed_steps,
                completed_advection_substeps: substeps,
                completed_maxwell_stefan_evaluations: evaluations,
                completed_caloric_iterations: caloric_iterations,
                maximum_accepted_error_ratio: maximum_error,
                minimum_accepted_step_s: minimum_step,
                maximum_accepted_step_s: maximum_step,
                suggested_next_step_s: next_step,
                maximum_advection_courant: maximum_courant,
                parcel_energy_resolution_bound_j: parcel_resolution,
                advected_mol,
                diffused_mol,
                transferred_mol,
                advected_energy_j,
                diffused_energy_j,
                transferred_energy_j,
                component_residual_mol,
                maximum_capacity_residual_mol,
                energy_residual_j,
                entropy_change_j_per_k,
            },
        })
    }
}

fn count(current: usize, added: usize) -> Result<usize> {
    current
        .checked_add(added)
        .ok_or(MulticomponentMiscibilityError::Arithmetic("coupled thermal work count").into())
}
fn total(mut a: Sum, b: Sum) -> Result<f64> {
    a.add(b.high)?;
    a.add(b.low)?;
    Ok(a.value()?)
}

#[cfg(test)]
mod tests;
