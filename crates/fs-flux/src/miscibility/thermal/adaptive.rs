//! Adaptive adiabatic miscibility with bounded step doubling.
//!
//! One full step is compared with two half steps. The fine result is accepted
//! only when both composition and temperature meet the requested local error.
//! Retryable thermodynamic refusals shrink the step. The complete requested
//! window is then audited once for species, energy and entropy so numerical
//! allowances do not grow with the number of accepted substeps.

use super::{
    BinaryMiscibilityThermalError as Error, BinaryMiscibilityThermalField,
    BinaryMiscibilityThermalSolve, Result as ThermalResult, total,
};
use crate::miscibility::{BinaryMiscibilityError, BinaryMiscibilitySolve, Sum, finite, nonnegative, poll, positive};
use crate::miscibility::adaptive::trial_endpoint;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryMiscibilityThermalTimeControl {
    pub initial_step_s: f64,
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    pub absolute_mole_fraction_tolerance: f64,
    pub relative_mole_fraction_tolerance: f64,
    pub absolute_temperature_tolerance_k: f64,
    pub relative_temperature_tolerance: f64,
    pub maximum_attempts: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityThermalAdvanceError {
    Thermal(Error),
    WorkLimit { attempted_intervals: usize, accepted_intervals: usize, elapsed_s: f64 },
    TemporalAccuracy { elapsed_s: f64, duration_s: f64, error_ratio: f64 },
}
impl From<Error> for BinaryMiscibilityThermalAdvanceError {
    fn from(error: Error) -> Self { Self::Thermal(error) }
}
impl core::fmt::Display for BinaryMiscibilityThermalAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive adiabatic miscibility: {self:?}")
    }
}
impl std::error::Error for BinaryMiscibilityThermalAdvanceError {}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityThermalAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_next_step_s: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityThermalAdvanceCandidate {
    field: BinaryMiscibilityThermalField,
    advance: BinaryMiscibilityThermalAdvance,
}
impl BinaryMiscibilityThermalAdvanceCandidate {
    pub const fn field(&self) -> &BinaryMiscibilityThermalField { &self.field }
    pub const fn advance(&self) -> &BinaryMiscibilityThermalAdvance { &self.advance }
    pub fn into_parts(self) -> (BinaryMiscibilityThermalField, BinaryMiscibilityThermalAdvance) {
        (self.field, self.advance)
    }
}

impl BinaryMiscibilityThermalField {
    pub fn advance_adiabatic_candidate(
        &self,
        duration_s: f64,
        time: BinaryMiscibilityThermalTimeControl,
        diffusion: BinaryMiscibilitySolve,
        thermal: BinaryMiscibilityThermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityThermalAdvanceCandidate, BinaryMiscibilityThermalAdvanceError> {
        poll(&mut cancelled).map_err(Error::from)?;
        self.admit_controls(thermal)?;
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
            || !positive(time.absolute_temperature_tolerance_k)
            || !nonnegative(time.relative_temperature_tolerance)
            || time.relative_temperature_tolerance > 1.0
            || !(1..=1_000_000).contains(&time.maximum_attempts)
        {
            return Err(Error::Invalid("positive bounded adaptive thermal controls required").into());
        }
        let initial_time = self.state.elapsed_s;
        let target = finite(initial_time + duration_s, "adaptive thermal target time").map_err(Error::from)?;
        if target <= initial_time {
            return Err(Error::Solver(BinaryMiscibilityError::Arithmetic("adaptive thermal interval is unrepresentable")).into());
        }

        let initial_temperature = self.network.model().temperature();
        let mut field = self.clone();
        let mut next_step = time.initial_step_s;
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut largest_error = 0.0_f64;
        let mut smallest_step = f64::INFINITY;
        let mut largest_step = 0.0_f64;

        while field.state.elapsed_s < target {
            poll(&mut cancelled).map_err(Error::from)?;
            if attempts == time.maximum_attempts {
                return Err(BinaryMiscibilityThermalAdvanceError::WorkLimit {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                    elapsed_s: field.state.elapsed_s,
                });
            }
            attempts += 1;
            let endpoint = trial_endpoint(field.state.elapsed_s, next_step, target);
            let step = endpoint - field.state.elapsed_s;
            if !positive(step) {
                return Err(Error::Solver(BinaryMiscibilityError::Arithmetic("adaptive thermal step cannot advance represented time")).into());
            }

            let trial = (|| -> ThermalResult<_> {
                let coarse = field.adiabatic_candidate(step, diffusion, thermal, &mut cancelled)?;
                let half = field.adiabatic_candidate(0.5 * step, diffusion, thermal, &mut cancelled)?;
                let second = endpoint - half.field.state.elapsed_s;
                let fine = half.field.adiabatic_candidate(second, diffusion, thermal, &mut cancelled)?;
                if coarse.field.state.elapsed_s != endpoint || fine.field.state.elapsed_s != endpoint {
                    return Err(Error::Solver(BinaryMiscibilityError::Arithmetic("adaptive thermal endpoint cannot be represented exactly")));
                }
                Ok((coarse, fine))
            })();

            let (coarse, fine) = match trial {
                Ok(value) => value,
                Err(error) => {
                    let retryable = matches!(
                        error,
                        Error::Solver(BinaryMiscibilityError::Root { .. })
                            | Error::Solver(BinaryMiscibilityError::FreeEnergy { .. })
                            | Error::Solver(BinaryMiscibilityError::Conservation { .. })
                            | Error::HeatOutsideBracket { .. }
                            | Error::IterationLimit
                            | Error::Energy { .. }
                            | Error::Entropy { .. }
                    );
                    let smaller = (0.5 * step).max(time.minimum_step_s);
                    if !retryable || smaller >= step { return Err(error.into()); }
                    next_step = smaller;
                    continue;
                }
            };

            let mut error_ratio = 0.0_f64;
            for (cell, &capacity) in self.network.owner.capacities_mol.iter().enumerate() {
                for species in 0..2 {
                    let difference = (fine.field.state.amounts[cell][species] - coarse.field.state.amounts[cell][species]).abs() / capacity;
                    let scale = field.state.fractions[cell][species].max(fine.field.state.fractions[cell][species]);
                    let allowance = time.absolute_mole_fraction_tolerance + time.relative_mole_fraction_tolerance * scale;
                    error_ratio = error_ratio.max(finite(difference / allowance, "adaptive thermal composition error").map_err(Error::from)?);
                }
            }
            let tc = coarse.field.network.model().temperature().value();
            let tf = fine.field.network.model().temperature().value();
            let scale = field.network.model().temperature().value().abs().max(tf.abs());
            let allowance = time.absolute_temperature_tolerance_k + time.relative_temperature_tolerance * scale;
            error_ratio = error_ratio.max(finite((tf - tc).abs() / allowance, "adaptive thermal temperature error").map_err(Error::from)?);

            let factor = if error_ratio == 0.0 { 2.0 } else { (0.8 / error_ratio.sqrt()).clamp(0.2, 2.0) };
            if error_ratio > 1.0 {
                let smaller = (step * factor.min(0.8)).max(time.minimum_step_s);
                if smaller >= step {
                    return Err(BinaryMiscibilityThermalAdvanceError::TemporalAccuracy { elapsed_s: field.state.elapsed_s, duration_s: step, error_ratio });
                }
                next_step = smaller;
                continue;
            }
            accepted += 1;
            largest_error = largest_error.max(error_ratio);
            smallest_step = smallest_step.min(step);
            largest_step = largest_step.max(step);
            next_step = (step * factor).clamp(time.minimum_step_s, time.maximum_step_s);
            field = fine.field;
        }

        let mut balances = [Sum::default(); 2];
        let mut capacity_residual = 0.0_f64;
        for (cell, (&old, &new)) in self.state.amounts.iter().zip(&field.state.amounts).enumerate() {
            for species in 0..2 { balances[species].add(finite(new[species] - old[species], "adaptive thermal component balance").map_err(Error::from)?).map_err(Error::from)?; }
            capacity_residual = capacity_residual.max((new[0] + new[1] - self.network.owner.capacities_mol[cell]).abs());
        }
        let balances = [balances[0].value().map_err(Error::from)?, balances[1].value().map_err(Error::from)?];
        if balances.iter().any(|x| x.abs() > diffusion.amount_balance_tolerance_mol) || capacity_residual > diffusion.amount_balance_tolerance_mol {
            return Err(Error::Solver(BinaryMiscibilityError::Conservation { component_residual_mol: balances, capacity_residual_mol: capacity_residual }).into());
        }

        let (mixing_h, mixing_s, interface) = self.composition_change(&field.state, &mut cancelled)?;
        let final_temperature = field.network.model().temperature();
        let sensible_h = self.sensible_change(initial_temperature, final_temperature, false, &mut cancelled)?;
        let sensible_s = self.sensible_change(initial_temperature, final_temperature, true, &mut cancelled)?;
        let energy_residual_j = total(&[mixing_h, interface, sensible_h])?;
        if energy_residual_j.abs() > thermal.energy_tolerance_j {
            return Err(Error::Energy { residual_j: energy_residual_j, tolerance_j: thermal.energy_tolerance_j }.into());
        }
        let entropy_change_j_per_k = total(&[mixing_s, sensible_s])?;
        if entropy_change_j_per_k < -thermal.entropy_tolerance_j_per_k {
            return Err(Error::Entropy { change_j_per_k: entropy_change_j_per_k, tolerance_j_per_k: thermal.entropy_tolerance_j_per_k }.into());
        }
        poll(&mut cancelled).map_err(Error::from)?;
        Ok(BinaryMiscibilityThermalAdvanceCandidate {
            field,
            advance: BinaryMiscibilityThermalAdvance {
                duration_s: target - initial_time,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: attempts - accepted,
                maximum_accepted_error_ratio: largest_error,
                minimum_accepted_step_s: smallest_step,
                maximum_accepted_step_s: largest_step,
                suggested_next_step_s: next_step,
                energy_residual_j,
                entropy_change_j_per_k,
            },
        })
    }
}
