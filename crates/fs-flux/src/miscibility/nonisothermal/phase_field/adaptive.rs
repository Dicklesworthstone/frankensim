//! Adaptive local-temperature demixing with an independent interface-energy norm.
//!
//! One frozen-temperature convex-split sweep is compared with two half sweeps.
//! Only the fine path contributes physical transfers. Composition, temperature
//! and the L1 difference of individual interface energies all control refinement;
//! compensating interface errors cannot disappear in a total-energy subtraction.
//! The complete returned interval rechecks the original conserved inventories,
//! homogeneous plus interfacial energy and bulk entropy. These are local time
//! error estimates, not spatial or global-trajectory accuracy certificates.

use super::*;
use crate::miscibility::adaptive::trial_endpoint;
use crate::miscibility::nonisothermal::adaptive::{
    BinaryNonisothermalAdvanceError, BinaryNonisothermalTimeControl, admit_time, error_ratio,
    factor, retryable,
};

/// Existing composition/Kelvin controls and one whole-graph interface allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryNonisothermalPhaseFieldTimeControl {
    pub thermal: BinaryNonisothermalTimeControl,
    /// Absolute allowance on the sum of absolute coarse/fine face-energy errors.
    pub absolute_interfacial_energy_tolerance_j: f64,
    /// Relative scale is the sum of max(coarse, fine) energies of each face.
    pub relative_interfacial_energy_tolerance: f64,
}

/// Accepted physical path and work from every completed trial sweep.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalPhaseFieldAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub completed_step_candidates: usize,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_next_step_s: f64,
    /// Completed coarse and rejected sweeps count; an incomplete sweep may do more work.
    pub root_iterations: u64,
    pub caloric_iterations: u64,
    /// Accepted face transfers, oriented by the original canonical graph.
    pub transferred_a_mol: Vec<f64>,
    pub transferred_energy_j: Vec<f64>,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub homogeneous_enthalpy_change_j: f64,
    pub interfacial_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalPhaseFieldAdvanceCandidate {
    field: BinaryNonisothermalPhaseField,
    advance: BinaryNonisothermalPhaseFieldAdvance,
}
impl BinaryNonisothermalPhaseFieldAdvanceCandidate {
    pub const fn field(&self) -> &BinaryNonisothermalPhaseField {
        &self.field
    }
    pub const fn advance(&self) -> &BinaryNonisothermalPhaseFieldAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        BinaryNonisothermalPhaseField,
        BinaryNonisothermalPhaseFieldAdvance,
    ) {
        (self.field, self.advance)
    }
}

impl BinaryNonisothermalPhaseField {
    /// Resolve a complete closed demixing/conduction interval. Refine failed
    /// thermal/entropy trials within explicit work and minimum-step limits;
    /// original source failures and cancellation immediately refuse. Failure
    /// leaves this field, its clock and every original binding unchanged.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: BinaryNonisothermalPhaseFieldTimeControl,
        controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<
        BinaryNonisothermalPhaseFieldAdvanceCandidate,
        BinaryNonisothermalAdvanceError,
    > {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        admit_time(duration_s, time.thermal)?;
        if !nonnegative(time.absolute_interfacial_energy_tolerance_j)
            || !nonnegative(time.relative_interfacial_energy_tolerance)
            || time.absolute_interfacial_energy_tolerance_j
                + time.relative_interfacial_energy_tolerance
                == 0.0
        {
            return Err(BinaryNonisothermalError::Invalid(
                "finite nonnegative and nonzero interface time-error tolerances required",
            )
            .into());
        }
        let t = time.thermal.composition;
        let target = finite(self.elapsed_s() + duration_s, "phase-field target time")
            .map_err(BinaryNonisothermalError::from)?;
        if target <= self.elapsed_s() {
            return Err(BinaryNonisothermalError::Invalid(
                "requested phase-field interval cannot advance time",
            )
            .into());
        }
        let mut field = self.clone();
        let mut next_step = t.initial_step_s;
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut completed = 0usize;
        let mut roots = 0u64;
        let mut caloric = 0u64;
        let mut maximum_error = 0.0_f64;
        let mut minimum_step = f64::INFINITY;
        let mut maximum_step = 0.0_f64;
        let count = self.network().owner.faces.len();
        let mut amounts = vec![Sum::default(); count];
        let mut energies = vec![Sum::default(); count];
        while field.elapsed_s() < target {
            poll(&mut cancelled)?;
            if attempts == t.maximum_attempts {
                return Err(BinaryNonisothermalAdvanceError::WorkLimit {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                    elapsed_s: field.elapsed_s() - self.elapsed_s(),
                });
            }
            let remaining = target - field.elapsed_s();
            let endpoint = trial_endpoint(field.elapsed_s(), next_step, target);
            let dt = endpoint - field.elapsed_s();
            if !positive(dt) {
                return Err(BinaryNonisothermalError::Invalid(
                    "phase-field trial cannot advance physical time",
                )
                .into());
            }
            let floor = t.minimum_step_s.min(remaining);
            attempts += 1;
            let trial = (|| -> Result<_> {
                let mut step = |start: &BinaryNonisothermalPhaseField, duration| {
                    let result = start.step_candidate(duration, controls, &mut cancelled)?;
                    completed += 1;
                    roots += result.step().root_iterations;
                    caloric += result.step().caloric_iterations;
                    Ok::<_, BinaryNonisothermalError>(result)
                };
                let coarse = step(&field, dt)?;
                let first = step(&field, 0.5 * dt)?;
                let second = step(first.field(), endpoint - first.field().elapsed_s())?;
                if coarse.field().elapsed_s() != endpoint || second.field().elapsed_s() != endpoint
                {
                    return Err(BinaryNonisothermalError::Invalid(
                        "phase-field trial endpoint is unrepresentable",
                    ));
                }
                let ratio = error_ratio(
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
                    next_step = (0.5 * dt).max(floor);
                    continue;
                }
                Err(error) => return Err(error.into()),
            };
            if ratio > 1.0 {
                if dt <= floor {
                    return Err(BinaryNonisothermalAdvanceError::TemporalAccuracy {
                        elapsed_s: field.elapsed_s() - self.elapsed_s(),
                        duration_s: dt,
                        error_ratio: ratio,
                    });
                }
                next_step = (dt * factor(ratio)).max(floor);
                continue;
            }
            for report in [first.step(), second.step()] {
                for face in 0..count {
                    poll(&mut cancelled)?;
                    amounts[face]
                        .add(report.transferred_a_mol[face])
                        .map_err(BinaryNonisothermalError::from)?;
                    energies[face]
                        .add(report.transferred_energy_j[face])
                        .map_err(BinaryNonisothermalError::from)?;
                }
            }
            field = second.into_parts().0;
            accepted += 1;
            maximum_error = maximum_error.max(ratio);
            minimum_step = minimum_step.min(dt);
            maximum_step = maximum_step.max(dt);
            next_step = (dt * factor(ratio)).clamp(t.minimum_step_s, t.maximum_step_s);
        }
        let (
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
            homogeneous_enthalpy_change_j,
            interfacial_energy_change_j,
        ) = self.balances(&field, controls, &mut cancelled)?;
        let mut transferred_a_mol = Vec::with_capacity(count);
        let mut transferred_energy_j = Vec::with_capacity(count);
        for face in 0..count {
            poll(&mut cancelled)?;
            transferred_a_mol.push(
                amounts[face]
                    .value()
                    .map_err(BinaryNonisothermalError::from)?,
            );
            transferred_energy_j.push(
                energies[face]
                    .value()
                    .map_err(BinaryNonisothermalError::from)?,
            );
        }
        Ok(BinaryNonisothermalPhaseFieldAdvanceCandidate {
            field,
            advance: BinaryNonisothermalPhaseFieldAdvance {
                duration_s,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: attempts - accepted,
                completed_step_candidates: completed,
                maximum_accepted_error_ratio: maximum_error,
                minimum_accepted_step_s: minimum_step,
                maximum_accepted_step_s: maximum_step,
                suggested_next_step_s: next_step,
                root_iterations: roots,
                caloric_iterations: caloric,
                transferred_a_mol,
                transferred_energy_j,
                component_residual_mol,
                maximum_capacity_residual_mol,
                homogeneous_enthalpy_change_j,
                interfacial_energy_change_j,
                energy_residual_j,
                entropy_change_j_per_k,
            },
        })
    }
}

fn interface_error_ratio(
    a: &BinaryNonisothermalPhaseField,
    b: &BinaryNonisothermalPhaseField,
    time: BinaryNonisothermalPhaseFieldTimeControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut discrepancy = Sum::default();
    let mut scale = Sum::default();
    for face in a.network().owner.faces.iter() {
        poll(cancelled)?;
        let [i, j] = face.physical.cells;
        let x = a.mole_fractions()[i][0] - a.mole_fractions()[j][0];
        let y = b.mole_fractions()[i][0] - b.mole_fractions()[j][0];
        let coefficient = 0.5 * face.gradient_j;
        discrepancy.add(((coefficient * (x - y)) * (x + y)).abs())?;
        let difference = x.abs().max(y.abs());
        scale.add((coefficient * difference) * difference)?;
    }
    let tolerance = time.absolute_interfacial_energy_tolerance_j
        + time.relative_interfacial_energy_tolerance * scale.value()?;
    if !tolerance.is_finite() {
        return Err(BinaryNonisothermalError::Invalid(
            "phase-field interface time-error tolerance overflow",
        ));
    }
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
