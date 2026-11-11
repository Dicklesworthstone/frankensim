//! Adaptive full Maxwell–Stefan transport on the original liquid field.
//!
//! A full explicit step is compared with two half steps. Only the two-half
//! result and its face transfers enter the accepted trajectory. Every species
//! has its own amount and temporal error, including positive trace inventories
//! and compensated changes smaller than one ulp of the stored amount.

use super::*;
use crate::miscibility::adaptive::trial_endpoint;

/// First-order local composition accuracy and finite physical/work bounds.
/// A final remainder below the minimum step is admitted only if it passes the
/// same accuracy, positivity, conservation and free-energy checks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentMiscibilityTimeControl {
    pub initial_step_s: f64,
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    pub absolute_mole_fraction_tolerance: f64,
    pub relative_mole_fraction_tolerance: f64,
    /// At most two complete face sweeps per attempted interval. Unchanged
    /// initial rates are reused when an interval is retried at a smaller step.
    pub maximum_attempts: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentMiscibilityAdvanceError {
    Solver(MulticomponentMiscibilityError),
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
impl From<MulticomponentMiscibilityError> for MulticomponentMiscibilityAdvanceError {
    fn from(error: MulticomponentMiscibilityError) -> Self {
        Self::Solver(error)
    }
}
impl core::fmt::Display for MulticomponentMiscibilityAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive multicomponent liquid miscibility: {self:?}")
    }
}
impl std::error::Error for MulticomponentMiscibilityAdvanceError {}

/// Accepted physical transfers and complete-window numerical balances.
/// Constitutive dissipation rates are instantaneous diagnostics; the actual
/// whole-window free-energy change below is independently checked. An explicit
/// method does not inherit the implicit binary owner's dissipation inequality.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    /// Includes completed coarse/fine steps from rejected trials. An incomplete
    /// refused step may perform additional work, bounded by maximum_attempts.
    pub completed_steps: usize,
    /// Actual completed constitutive evaluations, including rates from rejected
    /// trials. Shared full/first-half rates are counted once, even across retries.
    pub completed_maxwell_stefan_evaluations: usize,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_next_step_s: f64,
    /// Signed transfer from canonical face cell 0 to cell 1, then species order.
    /// Rejected/coarse trials never contribute to these physical receipts.
    pub transferred_mol: Vec<Vec<f64>>,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub free_energy_change_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityAdvanceCandidate {
    field: MulticomponentMiscibilityField,
    advance: MulticomponentMiscibilityAdvance,
}
impl MulticomponentMiscibilityAdvanceCandidate {
    pub const fn field(&self) -> &MulticomponentMiscibilityField {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentMiscibilityAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentMiscibilityField,
        MulticomponentMiscibilityAdvance,
    ) {
        (self.field, self.advance)
    }
}

impl MulticomponentMiscibilityField {
    /// Advance one closed physical interval atomically. Positivity and actual
    /// free-energy refusals refine the step; unsupported thermodynamics,
    /// transport data, arithmetic and cancellation remain immediate refusals.
    /// One original conservation/energy allowance applies to the entire window.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: MulticomponentMiscibilityTimeControl,
        controls: MulticomponentMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<
        MulticomponentMiscibilityAdvanceCandidate,
        MulticomponentMiscibilityAdvanceError,
    > {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        admit_time(duration_s, time)?;
        let target = finite(self.elapsed_s + duration_s, "multicomponent target time")?;
        if target <= self.elapsed_s {
            return Err(MulticomponentMiscibilityError::Invalid(
                "requested physical interval does not advance time",
            )
            .into());
        }
        let n = self.network.species_count();
        let faces = self.network.faces().len();
        let mut field = self.clone();
        let mut transferred = vec![Sum::default(); faces * n];
        let mut dt = time.initial_step_s.min(duration_s);
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut completed_steps = 0usize;
        let mut completed_evaluations = 0usize;
        let mut cached_rates: Option<solve::Rates> = None;
        let mut maximum_ratio = 0.0_f64;
        let mut minimum_step = f64::INFINITY;
        let mut maximum_step = 0.0_f64;
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
            let endpoint = trial_endpoint(field.elapsed_s, dt, target);
            dt = endpoint - field.elapsed_s;
            let floor = time.minimum_step_s.min(remaining);
            let half = 0.5 * dt;
            if !positive(half) || field.elapsed_s + half <= field.elapsed_s {
                return Err(MulticomponentMiscibilityError::Invalid(
                    "adaptive interval cannot be represented",
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
                    .expect("rates prepared for unchanged field");
                let coarse =
                    solve::candidate_from_rates(&field, dt, rates, controls, &mut cancelled)?;
                count_step(&mut completed_steps)?;
                let first =
                    solve::candidate_from_rates(&field, half, rates, controls, &mut cancelled)?;
                count_step(&mut completed_steps)?;
                let midpoint_rates = solve::rates(first.field(), controls, &mut cancelled)?;
                count_evaluations(&midpoint_rates, &mut completed_evaluations)?;
                let second = solve::candidate_from_rates(
                    first.field(),
                    endpoint - first.field().elapsed_s,
                    &midpoint_rates,
                    controls,
                    &mut cancelled,
                )?;
                count_step(&mut completed_steps)?;
                if coarse.field().elapsed_s != endpoint || second.field().elapsed_s != endpoint {
                    return Err(MulticomponentMiscibilityError::Invalid(
                        "adaptive trial endpoint is unrepresentable",
                    ));
                }
                let ratio = error_ratio(coarse.field(), second.field(), time, &mut cancelled)?;
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
                    return Err(MulticomponentMiscibilityAdvanceError::TemporalAccuracy {
                        elapsed_s: field.elapsed_s - self.elapsed_s,
                        duration_s: dt,
                        error_ratio: ratio,
                    });
                }
                dt = (dt * factor(ratio)).max(floor);
                continue;
            }
            for report in [first.step(), second.step()] {
                for (face, amounts) in report.transferred_mol.iter().enumerate() {
                    for (species, &amount) in amounts.iter().enumerate() {
                        transferred[face * n + species].add(amount)?;
                    }
                }
            }
            field = second.into_parts().0;
            cached_rates = None;
            accepted += 1;
            maximum_ratio = maximum_ratio.max(ratio);
            minimum_step = minimum_step.min(dt);
            maximum_step = maximum_step.max(dt);
            dt = (dt * factor(ratio)).clamp(time.minimum_step_s, time.maximum_step_s);
        }
        let (component_residual_mol, maximum_capacity_residual_mol, free_energy_change_j) =
            self.balances(&field, controls, &mut cancelled)?;
        let transferred_mol = transferred
            .chunks(n)
            .map(|amounts| amounts.iter().map(|amount| amount.value()).collect())
            .collect::<Result<Vec<Vec<f64>>>>()?;
        poll(&mut cancelled)?;
        Ok(MulticomponentMiscibilityAdvanceCandidate {
            field,
            advance: MulticomponentMiscibilityAdvance {
                duration_s,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                completed_steps,
                completed_maxwell_stefan_evaluations: completed_evaluations,
                maximum_accepted_error_ratio: maximum_ratio,
                minimum_accepted_step_s: minimum_step,
                maximum_accepted_step_s: maximum_step,
                suggested_next_step_s: dt,
                transferred_mol,
                component_residual_mol,
                maximum_capacity_residual_mol,
                free_energy_change_j,
            },
        })
    }
}

pub(super) fn count_step(steps: &mut usize) -> Result<()> {
    *steps = steps
        .checked_add(1)
        .ok_or(MulticomponentMiscibilityError::Arithmetic(
            "adaptive completed-step count",
        ))?;
    Ok(())
}

pub(super) fn count_evaluations(rates: &solve::Rates, evaluations: &mut usize) -> Result<()> {
    *evaluations = evaluations.checked_add(rates.evaluations()).ok_or(
        MulticomponentMiscibilityError::Arithmetic("adaptive constitutive-work count"),
    )?;
    Ok(())
}

pub(super) fn admit_time(duration: f64, time: MulticomponentMiscibilityTimeControl) -> Result<()> {
    if !positive(duration)
        || !positive(time.initial_step_s)
        || !positive(time.minimum_step_s)
        || !positive(time.maximum_step_s)
        || time.initial_step_s < time.minimum_step_s
        || time.initial_step_s > time.maximum_step_s
        || time.maximum_step_s < time.minimum_step_s
        || time.maximum_attempts == 0
        || !nonnegative(time.absolute_mole_fraction_tolerance)
        || !nonnegative(time.relative_mole_fraction_tolerance)
        || !positive(time.absolute_mole_fraction_tolerance + time.relative_mole_fraction_tolerance)
    {
        return Err(MulticomponentMiscibilityError::Invalid(
            "finite ordered time steps, positive work and temporal accuracy required",
        ));
    }
    Ok(())
}

pub(super) fn retryable(error: &MulticomponentMiscibilityError) -> bool {
    matches!(
        error,
        MulticomponentMiscibilityError::Positivity { .. }
            | MulticomponentMiscibilityError::FreeEnergy { .. }
    )
}

pub(super) fn factor(error: f64) -> f64 {
    if error == 0.0 {
        2.0
    } else {
        (0.9 / fs_math::det::sqrt(error)).clamp(0.25, 2.0)
    }
}

pub(super) fn error_ratio(
    a: &MulticomponentMiscibilityField,
    b: &MulticomponentMiscibilityField,
    time: MulticomponentMiscibilityTimeControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let n = a.network.species_count();
    let mut worst = 0.0_f64;
    for (cell, &capacity) in a.network.capacities_mol().iter().enumerate() {
        poll(cancelled)?;
        for species in 0..n {
            let index = cell * n + species;
            let tolerance = finite(
                time.absolute_mole_fraction_tolerance
                    + time.relative_mole_fraction_tolerance
                        * a.fractions[index].abs().max(b.fractions[index].abs()),
                "multicomponent temporal tolerance",
            )?;
            let mut difference = Sum::default();
            difference.add(a.amounts[index])?;
            difference.add(-b.amounts[index])?;
            difference.add(a.amount_remainders[index])?;
            difference.add(-b.amount_remainders[index])?;
            let amount_error = difference.value()?.abs();
            // A relative trace allowance may underflow. Test zero before the
            // capacity division so representable inventory errors cannot vanish.
            let ratio = if amount_error == 0.0 {
                0.0
            } else if tolerance == 0.0 {
                f64::INFINITY
            } else {
                // Oppositely scaled denominators have a representable product.
                // Otherwise both divisions move in the same direction, so an
                // intermediate under/overflow cannot hide a finite large error.
                if (capacity >= 1.0) == (tolerance >= 1.0) {
                    (amount_error / capacity) / tolerance
                } else {
                    amount_error / (capacity * tolerance)
                }
            };
            if ratio.is_nan() {
                return Err(MulticomponentMiscibilityError::Arithmetic(
                    "undefined multicomponent temporal error estimate",
                ));
            }
            worst = worst.max(ratio);
        }
    }
    Ok(worst)
}

#[cfg(test)]
mod tests;
