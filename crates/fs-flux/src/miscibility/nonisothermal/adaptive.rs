//! Local composition AND temperature error control for the first-order sweep.
//! Two accepted half sweeps retain the original material/caloric owners; a
//! rejected coarse/fine trial never mutates the caller's accepted field.

use super::*;
use crate::miscibility::adaptive::BinaryMiscibilityTimeControl;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryNonisothermalTimeControl {
    pub composition: BinaryMiscibilityTimeControl,
    pub absolute_temperature_tolerance_k: f64,
    pub relative_temperature_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BinaryNonisothermalAdvanceError {
    Solver(BinaryNonisothermalError),
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
impl From<BinaryNonisothermalError> for BinaryNonisothermalAdvanceError {
    fn from(error: BinaryNonisothermalError) -> Self {
        Self::Solver(error)
    }
}
impl core::fmt::Display for BinaryNonisothermalAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive nonisothermal binary diffusion: {self:?}")
    }
}
impl std::error::Error for BinaryNonisothermalAdvanceError {}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalAdvance {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error_ratio: f64,
    pub suggested_next_step_s: f64,
    pub root_iterations: u64,
    pub caloric_iterations: u64,
    pub transferred_a_mol: Vec<f64>,
    pub transferred_energy_j: Vec<f64>,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalAdvanceCandidate {
    field: BinaryNonisothermalField,
    advance: BinaryNonisothermalAdvance,
}
impl BinaryNonisothermalAdvanceCandidate {
    pub const fn field(&self) -> &BinaryNonisothermalField {
        &self.field
    }
    pub const fn advance(&self) -> &BinaryNonisothermalAdvance {
        &self.advance
    }
    pub fn into_parts(self) -> (BinaryNonisothermalField, BinaryNonisothermalAdvance) {
        (self.field, self.advance)
    }
}

impl BinaryNonisothermalField {
    /// Advance a closed interval with independent mole-fraction and Kelvin
    /// tolerances. Refine rejected numerical/thermal brackets; original source
    /// errors and cancellation remain explicit immediate refusals. Whole-window
    /// conservation and entropy use one original allowance, not summed budgets.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: BinaryNonisothermalTimeControl,
        controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<BinaryNonisothermalAdvanceCandidate, BinaryNonisothermalAdvanceError>
    {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        let t = time.composition;
        admit_time(duration_s, time)?;
        let target = finite(self.elapsed_s() + duration_s, "nonisothermal target time")
            .map_err(BinaryNonisothermalError::from)?;
        if target <= self.elapsed_s() {
            return Err(BinaryNonisothermalError::Invalid(
                "requested physical interval does not advance time",
            )
            .into());
        }
        let mut field = self.clone();
        let mut completed = 0.0;
        let mut dt = t.initial_step_s.min(duration_s);
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut maximum_ratio = 0.0_f64;
        let mut root_iterations = 0u64;
        let mut caloric_iterations = 0u64;
        let mut amount = vec![Sum::default(); self.owner.faces.len()];
        let mut energy = vec![Sum::default(); self.owner.faces.len()];
        while completed < duration_s {
            poll(&mut cancelled)?;
            if attempts == t.maximum_attempts {
                return Err(BinaryNonisothermalAdvanceError::WorkLimit {
                    attempted_intervals: attempts,
                    accepted_intervals: accepted,
                    elapsed_s: completed,
                });
            }
            let remaining = duration_s - completed;
            dt = dt.min(remaining);
            let floor = t.minimum_step_s.min(remaining);
            let half = 0.5 * dt;
            if !positive(half) || completed + dt <= completed {
                return Err(BinaryNonisothermalError::Invalid(
                    "adaptive interval cannot be represented",
                )
                .into());
            }
            attempts += 1;
            let trial = (|| -> Result<_> {
                let coarse = field.step_candidate(dt, controls, &mut cancelled)?;
                let first = field.step_candidate(half, controls, &mut cancelled)?;
                let second = first
                    .field()
                    .step_candidate(dt - half, controls, &mut cancelled)?;
                let ratio = error_ratio(coarse.field(), second.field(), time, &mut cancelled)?;
                Ok((coarse, first, second, ratio))
            })();
            let (coarse, first, second, ratio) = match trial {
                Ok(result) => result,
                Err(error) if retryable(&error) && dt > floor => {
                    rejected += 1;
                    dt = (0.5 * dt).max(floor);
                    continue;
                }
                Err(error) => return Err(error.into()),
            };
            root_iterations += coarse.step().root_iterations
                + first.step().root_iterations
                + second.step().root_iterations;
            caloric_iterations += coarse.step().caloric_iterations
                + first.step().caloric_iterations
                + second.step().caloric_iterations;
            if ratio > 1.0 {
                rejected += 1;
                if dt <= floor {
                    return Err(BinaryNonisothermalAdvanceError::TemporalAccuracy {
                        elapsed_s: completed,
                        duration_s: dt,
                        error_ratio: ratio,
                    });
                }
                dt = (dt * factor(ratio)).max(floor);
                continue;
            }
            for report in [first.step(), second.step()] {
                for face in 0..amount.len() {
                    amount[face]
                        .add(report.transferred_a_mol[face])
                        .map_err(BinaryNonisothermalError::from)?;
                    energy[face]
                        .add(report.transferred_energy_j[face])
                        .map_err(BinaryNonisothermalError::from)?;
                }
            }
            field = second.into_parts().0;
            accepted += 1;
            maximum_ratio = maximum_ratio.max(ratio);
            completed = if dt == remaining {
                duration_s
            } else {
                completed + dt
            };
            dt = (dt * factor(ratio)).clamp(t.minimum_step_s, t.maximum_step_s);
        }
        field.state.elapsed_s = target;
        let (
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
        ) = self.balances(&field, controls, &mut cancelled)?;
        let transferred_a_mol = amount
            .into_iter()
            .map(|s| s.value().map_err(BinaryNonisothermalError::from))
            .collect::<Result<Vec<_>>>()?;
        let transferred_energy_j = energy
            .into_iter()
            .map(|s| s.value().map_err(BinaryNonisothermalError::from))
            .collect::<Result<Vec<_>>>()?;
        poll(&mut cancelled)?;
        Ok(BinaryNonisothermalAdvanceCandidate {
            field,
            advance: BinaryNonisothermalAdvance {
                duration_s,
                attempted_intervals: attempts,
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                maximum_accepted_error_ratio: maximum_ratio,
                suggested_next_step_s: dt,
                root_iterations,
                caloric_iterations,
                transferred_a_mol,
                transferred_energy_j,
                component_residual_mol,
                maximum_capacity_residual_mol,
                energy_residual_j,
                entropy_change_j_per_k,
            },
        })
    }
}

pub(super) fn admit_time(duration_s: f64, time: BinaryNonisothermalTimeControl) -> Result<()> {
    let t = time.composition;
    if !positive(duration_s)
        || !positive(t.initial_step_s)
        || !positive(t.minimum_step_s)
        || !positive(t.maximum_step_s)
        || t.minimum_step_s > t.initial_step_s
        || t.initial_step_s > t.maximum_step_s
        || !(1..=1_000_000).contains(&t.maximum_attempts)
        || !nonnegative(t.absolute_mole_fraction_tolerance)
        || !nonnegative(t.relative_mole_fraction_tolerance)
        || t.absolute_mole_fraction_tolerance + t.relative_mole_fraction_tolerance == 0.0
        || !nonnegative(time.absolute_temperature_tolerance_k)
        || !nonnegative(time.relative_temperature_tolerance)
        || time.absolute_temperature_tolerance_k + time.relative_temperature_tolerance == 0.0
    {
        return Err(BinaryNonisothermalError::Invalid(
            "finite ordered steps, positive work and composition/temperature accuracy required",
        ));
    }
    Ok(())
}

pub(super) fn retryable(error: &BinaryNonisothermalError) -> bool {
    matches!(
        error,
        BinaryNonisothermalError::PairIterationLimit
            | BinaryNonisothermalError::PairOutsideTemperatureDomain
            | BinaryNonisothermalError::Energy { .. }
            | BinaryNonisothermalError::Entropy { .. }
            | BinaryNonisothermalError::Caloric(BinaryCaloricErrorV1::HeatOutsideBracket { .. })
            | BinaryNonisothermalError::Caloric(BinaryCaloricErrorV1::IterationLimit)
    )
}
pub(super) fn factor(error: f64) -> f64 {
    if error == 0.0 {
        2.0
    } else {
        (0.9 / error.sqrt()).clamp(0.25, 2.0)
    }
}
pub(super) fn error_ratio(
    a: &BinaryNonisothermalField,
    b: &BinaryNonisothermalField,
    time: BinaryNonisothermalTimeControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut worst = 0.0_f64;
    for cell in 0..a.temperatures.len() {
        poll(cancelled)?;
        for s in 0..2 {
            let [x, y] = [a.state.fractions[cell][s], b.state.fractions[cell][s]];
            let tolerance = time.composition.absolute_mole_fraction_tolerance
                + time.composition.relative_mole_fraction_tolerance * x.abs().max(y.abs());
            // Include independent trace inventories and carried sub-ULP changes
            // even when the dominant mole fraction rounds to one.
            let mut difference = Sum::default();
            difference.add(a.state.amounts[cell][s])?;
            difference.add(-b.state.amounts[cell][s])?;
            difference.add(a.state.amount_remainders[cell][s])?;
            difference.add(-b.state.amount_remainders[cell][s])?;
            let amount_error = difference.value()?.abs();
            // Test a zero tolerance before normalizing: a representable amount
            // difference must not disappear when its fraction underflows.
            let fraction_error = if tolerance == 0.0 {
                amount_error
            } else {
                amount_error / a.network().owner.capacities_mol[cell]
            };
            worst = worst.max(scaled_difference(fraction_error, 0.0, tolerance)?);
        }
        let [x, y] = [a.temperatures[cell].value(), b.temperatures[cell].value()];
        let tolerance = time.absolute_temperature_tolerance_k
            + time.relative_temperature_tolerance * x.abs().max(y.abs());
        worst = worst.max(scaled_difference(x, y, tolerance)?);
    }
    if worst.is_nan() {
        return Err(BinaryNonisothermalError::Invalid(
            "undefined coupled temporal error estimate",
        ));
    }
    Ok(worst)
}

fn scaled_difference(x: f64, y: f64, tolerance: f64) -> Result<f64> {
    if !tolerance.is_finite() {
        return Err(BinaryNonisothermalError::Invalid(
            "coupled temporal tolerance overflow",
        ));
    }
    let difference = (x - y).abs();
    // A relative trace tolerance may underflow. Exact agreement passes; a
    // representable disagreement must refine or refuse, never become NaN and
    // disappear from the maximum error estimate.
    Ok(if difference == 0.0 {
        0.0
    } else if tolerance == 0.0 {
        f64::INFINITY
    } else {
        difference / tolerance
    })
}
