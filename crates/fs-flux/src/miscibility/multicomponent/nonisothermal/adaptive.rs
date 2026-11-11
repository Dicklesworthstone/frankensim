//! Bounded adaptive species/heat evolution with reused original face rates.

use super::super::adaptive::{self as composition, MulticomponentMiscibilityTimeControl};
use super::*;
use crate::miscibility::adaptive::trial_endpoint;

/// Independent composition and Kelvin accuracy for the same complete step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentNonisothermalTimeControl {
    pub composition: MulticomponentMiscibilityTimeControl,
    pub absolute_temperature_tolerance_k: f64,
    pub relative_temperature_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentNonisothermalAdvanceError {
    Solver(MulticomponentNonisothermalError),
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
impl From<MulticomponentNonisothermalError> for MulticomponentNonisothermalAdvanceError {
    fn from(error: MulticomponentNonisothermalError) -> Self {
        Self::Solver(error)
    }
}
impl From<MulticomponentMiscibilityError> for MulticomponentNonisothermalAdvanceError {
    fn from(error: MulticomponentMiscibilityError) -> Self {
        Self::Solver(error.into())
    }
}
impl core::fmt::Display for MulticomponentNonisothermalAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "adaptive nonisothermal multicomponent miscibility: {self:?}"
        )
    }
}
impl std::error::Error for MulticomponentNonisothermalAdvanceError {}

/// Only accepted half-step transfers enter physical receipts. Work counts also
/// include successful stages of rejected trials. Incomplete refused stages may
/// perform additional work, bounded by the explicit attempt/root/face budgets.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub completed_steps: usize,
    pub completed_maxwell_stefan_evaluations: usize,
    pub completed_caloric_iterations: u64,
    pub maximum_accepted_error_ratio: f64,
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_next_step_s: f64,
    pub transferred_mol: Vec<Vec<f64>>,
    pub transferred_energy_j: Vec<f64>,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalAdvanceCandidate {
    field: MulticomponentNonisothermalField,
    advance: MulticomponentNonisothermalAdvance,
}
impl MulticomponentNonisothermalAdvanceCandidate {
    pub const fn field(&self) -> &MulticomponentNonisothermalField {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentNonisothermalAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentNonisothermalField,
        MulticomponentNonisothermalAdvance,
    ) {
        (self.field, self.advance)
    }
}

impl MulticomponentNonisothermalField {
    /// Compare one full coupled step with two half steps, accepting only the
    /// latter. All original species (including low words) and local temperatures
    /// participate in temporal error. Rates at the unchanged initial state are
    /// reused for coarse/first-half trials and rejected interval refinements.
    /// One original H/S and inventory allowance applies to the whole window.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: MulticomponentNonisothermalTimeControl,
        controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<
        MulticomponentNonisothermalAdvanceCandidate,
        MulticomponentNonisothermalAdvanceError,
    > {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        admit_time(duration_s, time)?;
        let target = finite(self.elapsed_s() + duration_s, "thermal target time")?;
        if target <= self.elapsed_s() {
            return Err(MulticomponentNonisothermalError::Invalid(
                "requested thermal interval does not advance time",
            )
            .into());
        }
        let n = self.network().species_count();
        let faces = self.network().faces().len();
        let mut field = self.clone();
        let mut transferred = vec![Sum::default(); faces * n];
        let mut energy_transferred = vec![Sum::default(); faces];
        let mut dt = time.composition.initial_step_s.min(duration_s);
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut completed_steps = 0usize;
        let mut completed_evaluations = 0usize;
        let mut completed_caloric_iterations = 0u64;
        let mut cached_rates: Option<solve::Rates> = None;
        let mut maximum_ratio = 0.0_f64;
        let mut minimum_step = f64::INFINITY;
        let mut maximum_step = 0.0_f64;
        while field.elapsed_s() < target {
            poll(&mut cancelled)?;
            if attempts == time.composition.maximum_attempts {
                return Err(MulticomponentNonisothermalAdvanceError::WorkLimit {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                    elapsed_s: field.elapsed_s() - self.elapsed_s(),
                });
            }
            let remaining = target - field.elapsed_s();
            let endpoint = trial_endpoint(field.elapsed_s(), dt, target);
            dt = endpoint - field.elapsed_s();
            let floor = time.composition.minimum_step_s.min(remaining);
            let half = 0.5 * dt;
            if !positive(half) || field.elapsed_s() + half <= field.elapsed_s() {
                return Err(MulticomponentNonisothermalError::Invalid(
                    "adaptive thermal interval cannot be represented",
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
                    .expect("rates prepared for unchanged thermal field");
                let coarse =
                    solve::candidate_from_rates(&field, dt, rates, controls, &mut cancelled)?;
                count_step(
                    coarse.step(),
                    &mut completed_steps,
                    &mut completed_caloric_iterations,
                )?;
                let first =
                    solve::candidate_from_rates(&field, half, rates, controls, &mut cancelled)?;
                count_step(
                    first.step(),
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
                count_step(
                    second.step(),
                    &mut completed_steps,
                    &mut completed_caloric_iterations,
                )?;
                if coarse.field().elapsed_s() != endpoint || second.field().elapsed_s() != endpoint
                {
                    return Err(MulticomponentNonisothermalError::Invalid(
                        "adaptive thermal endpoint is unrepresentable",
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
                    return Err(MulticomponentNonisothermalAdvanceError::TemporalAccuracy {
                        elapsed_s: field.elapsed_s() - self.elapsed_s(),
                        duration_s: dt,
                        error_ratio: ratio,
                    });
                }
                dt = (dt * composition::factor(ratio)).max(floor);
                continue;
            }
            for report in [first.step(), second.step()] {
                for (face, amounts) in report.transferred_mol.iter().enumerate() {
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
            dt = (dt * composition::factor(ratio)).clamp(
                time.composition.minimum_step_s,
                time.composition.maximum_step_s,
            );
        }
        let (
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
        ) = self.balances(&field, controls, &mut cancelled)?;
        let transferred_mol = transferred
            .chunks(n)
            .map(|amounts| amounts.iter().map(Sum::value).collect())
            .collect::<super::super::Result<Vec<Vec<f64>>>>()?;
        let transferred_energy_j = energy_transferred
            .iter()
            .map(Sum::value)
            .collect::<super::super::Result<Vec<_>>>()?;
        poll(&mut cancelled)?;
        Ok(MulticomponentNonisothermalAdvanceCandidate {
            field,
            advance: MulticomponentNonisothermalAdvance {
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
            },
        })
    }
}

fn count_step(
    step: &MulticomponentNonisothermalStep,
    steps: &mut usize,
    iterations: &mut u64,
) -> Result<()> {
    count_completed_step(step.caloric_iterations, steps, iterations)
}

pub(in crate::miscibility::multicomponent) fn count_completed_step(
    caloric_iterations: u64,
    steps: &mut usize,
    iterations: &mut u64,
) -> Result<()> {
    *steps = steps
        .checked_add(1)
        .ok_or(MulticomponentMiscibilityError::Arithmetic(
            "thermal completed-step count",
        ))?;
    *iterations = iterations.checked_add(caloric_iterations).ok_or(
        MulticomponentMiscibilityError::Arithmetic("thermal completed-caloric count"),
    )?;
    Ok(())
}
pub(in crate::miscibility::multicomponent) fn count_evaluations(
    rates: &solve::Rates,
    evaluations: &mut usize,
) -> Result<()> {
    *evaluations = evaluations.checked_add(rates.evaluations()).ok_or(
        MulticomponentMiscibilityError::Arithmetic("thermal constitutive-work count"),
    )?;
    Ok(())
}
pub(in crate::miscibility::multicomponent) fn admit_time(
    duration: f64,
    time: MulticomponentNonisothermalTimeControl,
) -> Result<()> {
    composition::admit_time(duration, time.composition)?;
    if !nonnegative(time.absolute_temperature_tolerance_k)
        || !nonnegative(time.relative_temperature_tolerance)
        || !positive(time.absolute_temperature_tolerance_k + time.relative_temperature_tolerance)
    {
        return Err(MulticomponentNonisothermalError::Invalid(
            "finite nonnegative Kelvin temporal allowances with positive combined accuracy required",
        ));
    }
    Ok(())
}
pub(in crate::miscibility::multicomponent) fn retryable(error: &MulticomponentNonisothermalError) -> bool {
    matches!(
        error,
        MulticomponentNonisothermalError::Spatial(
            MulticomponentMiscibilityError::Positivity { .. }
        ) | MulticomponentNonisothermalError::Energy { .. }
            | MulticomponentNonisothermalError::Entropy { .. }
            | MulticomponentNonisothermalError::Caloric(MulticomponentCaloricErrorV1::Caloric(
                BinaryCaloricErrorV1::HeatOutsideBracket { .. }
                    | BinaryCaloricErrorV1::IterationLimit
            ))
    )
}
pub(in crate::miscibility::multicomponent) fn error_ratio(
    a: &MulticomponentNonisothermalField,
    b: &MulticomponentNonisothermalField,
    time: MulticomponentNonisothermalTimeControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut worst = composition::error_ratio(&a.state, &b.state, time.composition, cancelled)?;
    for (&left, &right) in a.temperatures.iter().zip(&b.temperatures) {
        poll(cancelled)?;
        let tolerance = finite(
            time.absolute_temperature_tolerance_k
                + time.relative_temperature_tolerance * left.value().abs().max(right.value().abs()),
            "thermal temporal temperature allowance",
        )?;
        let difference = (left.value() - right.value()).abs();
        let ratio = if difference == 0.0 {
            0.0
        } else if tolerance == 0.0 {
            f64::INFINITY
        } else {
            difference / tolerance
        };
        if ratio.is_nan() {
            return Err(MulticomponentNonisothermalError::Invalid(
                "undefined thermal temporal error estimate",
            ));
        }
        worst = worst.max(ratio);
    }
    Ok(worst)
}

#[cfg(test)]
mod tests;
