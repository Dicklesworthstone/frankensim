//! Coupled composition, Kelvin and interface-energy accuracy for demixing.
//!
//! One full explicit interval is compared with two half intervals. Every stage
//! evaluates the original thermal chemical forces and recovers bulk enthalpy
//! after half-edge energy allocation. Only accepted fine-path transfers enter
//! physical receipts. The complete interval rechecks original inventories,
//! total bulk-plus-interface energy, caloric precision and original entropy.

use super::*;
use crate::miscibility::adaptive::trial_endpoint;
use crate::miscibility::multicomponent::nonnegative;
use crate::miscibility::multicomponent::adaptive::factor;
use crate::miscibility::multicomponent::nonisothermal::adaptive::{
    MulticomponentNonisothermalTimeControl, admit_time, count_completed_step,
    count_evaluations, error_ratio as thermal_error_ratio, retryable,
};

pub use crate::miscibility::multicomponent::nonisothermal::adaptive::MulticomponentNonisothermalAdvanceError
    as MulticomponentNonisothermalPhaseFieldAdvanceError;

/// Original composition/Kelvin controls and one independent interface norm.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentNonisothermalPhaseFieldTimeControl {
    pub thermal: MulticomponentNonisothermalTimeControl,
    /// Allowance on the sum of absolute coarse/fine per-edge energy differences.
    pub absolute_interfacial_energy_tolerance_j: f64,
    /// Scale is the sum of max(coarse, fine) energies of each physical edge.
    pub relative_interfacial_energy_tolerance: f64,
}

/// Accepted transfers, actual completed work and complete-window physical laws.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalPhaseFieldAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub completed_steps: usize,
    /// Shared initial rates count once, including across rejected retries.
    /// Incomplete refused evaluations may perform additional bounded work.
    pub completed_maxwell_stefan_evaluations: usize,
    pub completed_caloric_iterations: u64,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_next_step_s: f64,
    /// Canonical face then species; excludes every coarse or rejected transfer.
    pub transferred_mol: Vec<Vec<f64>>,
    /// Paired total local-energy transfer, including original heat carriage.
    pub transferred_energy_j: Vec<f64>,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
    pub homogeneous_enthalpy_change_j: f64,
    pub interfacial_energy_change_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalPhaseFieldAdvanceCandidate {
    field: MulticomponentNonisothermalPhaseField,
    advance: MulticomponentNonisothermalPhaseFieldAdvance,
}
impl MulticomponentNonisothermalPhaseFieldAdvanceCandidate {
    pub const fn field(&self) -> &MulticomponentNonisothermalPhaseField {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentNonisothermalPhaseFieldAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentNonisothermalPhaseField,
        MulticomponentNonisothermalPhaseFieldAdvance,
    ) {
        (self.field, self.advance)
    }
}

impl MulticomponentNonisothermalPhaseField {
    /// Advance one closed thermal demixing interval atomically. Refine failed
    /// positivity, caloric bracket, energy, entropy and temporal trials within
    /// explicit work and time limits. Original source, precision and cancellation
    /// refusals remain hard failures. These local error estimates do not imply
    /// unconditional stability or a global trajectory/spatial error bound.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: MulticomponentNonisothermalPhaseFieldTimeControl,
        controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<
        MulticomponentNonisothermalPhaseFieldAdvanceCandidate,
        MulticomponentNonisothermalPhaseFieldAdvanceError,
    > {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        admit_time(duration_s, time.thermal)?;
        if !nonnegative(time.absolute_interfacial_energy_tolerance_j)
            || !nonnegative(time.relative_interfacial_energy_tolerance)
            || !positive(
                time.absolute_interfacial_energy_tolerance_j
                    + time.relative_interfacial_energy_tolerance,
            )
        {
            return Err(MulticomponentNonisothermalError::Invalid(
                "finite nonnegative and nonzero thermal interface time-error allowances required",
            )
            .into());
        }
        let target = finite(self.elapsed_s() + duration_s, "thermal phase-field target time")?;
        if target <= self.elapsed_s() {
            return Err(MulticomponentNonisothermalError::Invalid(
                "requested thermal phase-field interval does not advance time",
            )
            .into());
        }
        let t = time.thermal.composition;
        let n = self.network().species_count();
        let faces = self.network().faces().len();
        let mut field = self.clone();
        let mut transferred = vec![Sum::default(); faces * n];
        let mut energy_transferred = vec![Sum::default(); faces];
        let mut dt = t.initial_step_s.min(duration_s);
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut completed_steps = 0usize;
        let mut completed_evaluations = 0usize;
        let mut completed_caloric_iterations = 0u64;
        let mut cached_rates: Option<thermal_solve::Rates> = None;
        let mut maximum_ratio = 0.0_f64;
        let mut minimum_step = f64::INFINITY;
        let mut maximum_step = 0.0_f64;
        while field.elapsed_s() < target {
            poll(&mut cancelled)?;
            if attempts == t.maximum_attempts {
                return Err(MulticomponentNonisothermalPhaseFieldAdvanceError::WorkLimit {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                    elapsed_s: field.elapsed_s() - self.elapsed_s(),
                });
            }
            let remaining = target - field.elapsed_s();
            let endpoint = trial_endpoint(field.elapsed_s(), dt, target);
            dt = endpoint - field.elapsed_s();
            let floor = t.minimum_step_s.min(remaining);
            let half = 0.5 * dt;
            if !positive(half) || field.elapsed_s() + half <= field.elapsed_s() {
                return Err(MulticomponentNonisothermalError::Invalid(
                    "adaptive thermal phase-field interval cannot be represented",
                )
                .into());
            }
            attempts += 1;
            let trial = (|| -> Result<_> {
                if cached_rates.is_none() {
                    let rates = solve::rates(&field, controls, &mut cancelled)?;
                    count_evaluations(&rates, &mut completed_evaluations)?;
                    cached_rates = Some(rates);
                }
                let rates = cached_rates
                    .as_ref()
                    .expect("rates prepared for unchanged thermal phase field");
                let coarse =
                    solve::candidate_from_rates(&field, dt, rates, controls, &mut cancelled)?;
                count_completed_step(
                    coarse.step().caloric_iterations,
                    &mut completed_steps,
                    &mut completed_caloric_iterations,
                )?;
                let first =
                    solve::candidate_from_rates(&field, half, rates, controls, &mut cancelled)?;
                count_completed_step(
                    first.step().caloric_iterations,
                    &mut completed_steps,
                    &mut completed_caloric_iterations,
                )?;
                let midpoint_rates = solve::rates(first.field(), controls, &mut cancelled)?;
                count_evaluations(&midpoint_rates, &mut completed_evaluations)?;
                let second = solve::candidate_from_rates(
                    first.field(),
                    endpoint - first.field().elapsed_s(),
                    &midpoint_rates,
                    controls,
                    &mut cancelled,
                )?;
                count_completed_step(
                    second.step().caloric_iterations,
                    &mut completed_steps,
                    &mut completed_caloric_iterations,
                )?;
                if coarse.field().elapsed_s() != endpoint || second.field().elapsed_s() != endpoint {
                    return Err(MulticomponentNonisothermalError::Invalid(
                        "adaptive thermal phase-field endpoint is unrepresentable",
                    ));
                }
                let ratio = thermal_error_ratio(
                    &coarse.field().bulk,
                    &second.field().bulk,
                    time.thermal,
                    &mut cancelled,
                )?
                .max(interface_error_ratio(
                    coarse.field(),
                    second.field(),
                    time,
                    &mut cancelled,
                )?);
                Ok((first, second, ratio))
            })();
            let (first, second, ratio) = match trial {
                Ok(result) => result,
                Err(error) if retryable(&error) && dt > floor => {
                    rejected += 1;
                    dt = (0.5 * dt).max(floor);
                    continue;
                }
                Err(error) => return Err(error.into()),
            };
            if ratio > 1.0 {
                rejected += 1;
                if dt <= floor {
                    return Err(MulticomponentNonisothermalPhaseFieldAdvanceError::TemporalAccuracy {
                        elapsed_s: field.elapsed_s() - self.elapsed_s(),
                        duration_s: dt,
                        error_ratio: ratio,
                    });
                }
                dt = (dt * factor(ratio)).max(floor);
                continue;
            }
            for report in [first.step(), second.step()] {
                for (face, amounts) in report.transferred_mol.iter().enumerate() {
                    poll(&mut cancelled)?;
                    for (species, &amount) in amounts.iter().enumerate() {
                        transferred[face * n + species].add(amount)?;
                    }
                    energy_transferred[face].add(report.transferred_energy_j[face])?;
                }
            }
            field = second.into_parts().0;
            cached_rates = None;
            accepted += 1;
            maximum_ratio = maximum_ratio.max(ratio);
            minimum_step = minimum_step.min(dt);
            maximum_step = maximum_step.max(dt);
            dt = (dt * factor(ratio)).clamp(t.minimum_step_s, t.maximum_step_s);
        }
        let (
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
            homogeneous_enthalpy_change_j,
            interfacial_energy_change_j,
        ) = self.balances(&field, controls, &mut cancelled)?;
        let mut transferred_mol = Vec::with_capacity(faces);
        let mut transferred_energy_j = Vec::with_capacity(faces);
        for (amounts, energy) in transferred.chunks(n).zip(&energy_transferred) {
            poll(&mut cancelled)?;
            transferred_mol.push(
                amounts.iter().map(Sum::value)
                    .collect::<crate::miscibility::multicomponent::Result<Vec<_>>>()?,
            );
            transferred_energy_j.push(energy.value()?);
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentNonisothermalPhaseFieldAdvanceCandidate {
            field,
            advance: MulticomponentNonisothermalPhaseFieldAdvance {
                duration_s,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                completed_steps,
                completed_maxwell_stefan_evaluations: completed_evaluations,
                completed_caloric_iterations,
                maximum_accepted_error_ratio: maximum_ratio,
                minimum_accepted_step_s: minimum_step,
                maximum_accepted_step_s: maximum_step,
                suggested_next_step_s: dt,
                transferred_mol,
                transferred_energy_j,
                component_residual_mol,
                maximum_capacity_residual_mol,
                energy_residual_j,
                entropy_change_j_per_k,
                homogeneous_enthalpy_change_j,
                interfacial_energy_change_j,
            },
        })
    }
}

/// Whole-graph L1 edge-energy discrepancy; signed cancellation is forbidden.
fn interface_error_ratio(
    a: &MulticomponentNonisothermalPhaseField,
    b: &MulticomponentNonisothermalPhaseField,
    time: MulticomponentNonisothermalPhaseFieldTimeControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut discrepancy = Sum::default();
    let mut scale = Sum::default();
    for face in 0..a.network().faces().len() {
        poll(cancelled)?;
        discrepancy.add(solve::edge_energy_change_j(a, b, face)?.abs())?;
        scale.add(solve::edge_energy_j(a, face)?.max(solve::edge_energy_j(b, face)?))?;
    }
    let tolerance = finite(
        time.absolute_interfacial_energy_tolerance_j
            + time.relative_interfacial_energy_tolerance * scale.value()?,
        "thermal phase-field interface temporal allowance",
    )?;
    let discrepancy = discrepancy.value()?;
    Ok(if discrepancy == 0.0 {
        0.0
    } else if tolerance == 0.0 {
        f64::INFINITY
    } else {
        discrepancy / tolerance
    })
}

#[cfg(test)]
mod tests;
