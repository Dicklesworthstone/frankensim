//! Conservative adaptive intervals over the existing source-owned midpoint step.
//!
//! Step doubling estimates local error in position, gas temperature and gas
//! pressure. Only the two half-step states and transfers are accepted; the
//! coarse trial is never published or Richardson-extrapolated. A successful
//! result advances the complete requested interval. Tube ends still require
//! their physical event/boundary model and never become a clipped success.

use super::{
    SealedCapillaryError, SealedCapillaryInput, SealedCapillarySolve, SealedCapillaryStep,
    SealedGasCapillary,
};
use crate::capillary::ContactLineMotion;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_time::PiController;

/// Local-error and bounded-work controls for implicit-midpoint step doubling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SealedCapillaryAdvanceControls {
    pub initial_step_s: f64,
    /// Smallest trial interval; its two accepted half steps may be half this size.
    /// The final interval may be shorter to reach the requested time exactly.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// Counts accepted and rejected doubling trials, each with at most three steps.
    pub maximum_attempts: usize,
    pub relative_tolerance: f64,
    pub length_absolute_tolerance_m: f64,
    pub temperature_absolute_tolerance_k: f64,
    pub pressure_absolute_tolerance_pa: f64,
}

/// Failure leaves the original owner unchanged, including all staged substeps.
#[derive(Clone, Debug, PartialEq)]
pub enum SealedCapillaryAdvanceError {
    Invalid(&'static str),
    Cancelled,
    Step {
        attempted_step_s: f64,
        error: SealedCapillaryError,
    },
    StepBudget {
        attempts: usize,
        accepted_intervals: usize,
    },
    Tolerance {
        normalized_error: f64,
        step_s: f64,
    },
    StepUnderflow,
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for SealedCapillaryAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive sealed capillary: {self:?}")
    }
}
impl std::error::Error for SealedCapillaryAdvanceError {}

/// Sums of actual accepted substep transfers, never transfers from rejected trials.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SealedCapillaryAdvanceStep {
    pub initial_wet_length_m: f64,
    pub wet_length_m: f64,
    pub elapsed_time_s: f64,
    pub transferred_liquid_volume_m3: f64,
    pub transferred_liquid_mass_kg: f64,
    pub inlet_pressure_work_j: f64,
    pub wall_surface_energy_change_j: f64,
    pub gravitational_energy_change_j: f64,
    pub liquid_viscous_dissipation_j: f64,
    pub contact_line_viscous_dissipation_j: f64,
    pub pinning_dissipation_j: f64,
    pub gas_supplied_heat_j: f64,
    pub gas_work_by_gas_j: f64,
    pub gas_stored_energy_change_j: f64,
    /// Thermodynamic entropy change, not a claim of nonnegative entropy production.
    pub gas_entropy_change_j_k: f64,
    pub energy_balance_residual_j: f64,
    pub attempts: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub step_evaluations: usize,
    /// Largest accepted local estimate; it is not a global error enclosure.
    pub maximum_normalized_local_error: f64,
}

#[derive(Clone, Debug)]
pub struct SealedCapillaryAdvanceCandidate<'a> {
    endpoint: SealedGasCapillary<'a>,
    step: SealedCapillaryAdvanceStep,
}
impl<'a> SealedCapillaryAdvanceCandidate<'a> {
    #[must_use]
    pub const fn endpoint(&self) -> &SealedGasCapillary<'a> {
        &self.endpoint
    }
    #[must_use]
    pub const fn step(&self) -> &SealedCapillaryAdvanceStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(self) -> (SealedGasCapillary<'a>, SealedCapillaryAdvanceStep) {
        (self.endpoint, self.step)
    }
}

impl<'a> SealedGasCapillary<'a> {
    /// Advance a full constant-boundary interval with controlled local error.
    ///
    /// Each accepted interval uses two conservative midpoint half steps. Local
    /// estimates compare their position, temperature and pressure against one
    /// coarse step and divide by three on a common smooth motion branch. Across
    /// a detected pinning or reversal transition the full difference is used,
    /// since the smooth second-order error cancellation no longer applies.
    /// A controller rejection consumes no physical time or heat. Boundary changes must be
    /// split into separate calls; adaptivity cannot discover an unsampled pulse.
    /// The original source domains and single-step energy gates remain active.
    pub fn advance_candidate(
        &self,
        input: SealedCapillaryInput,
        solve: SealedCapillarySolve,
        controls: SealedCapillaryAdvanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SealedCapillaryAdvanceCandidate<'a>, SealedCapillaryAdvanceError> {
        use SealedCapillaryAdvanceError as Error;
        poll(&mut cancelled)?;
        validate(input, solve, controls)?;
        // Existing PI settings, with order-three local-error exponents for
        // midpoint step doubling rather than RK45's order-five defaults.
        let pi = PiController {
            k_p: 0.7 / 3.0,
            k_i: 0.4 / 3.0,
            ..PiController::default()
        };
        let mut endpoint: Option<Self> = None;
        let mut elapsed = 0.0;
        let mut proposal = controls.initial_step_s;
        let mut previous_error: f64 = 1.0;
        let mut previous_branch = None;
        let mut attempts = 0;
        let mut accepted = 0;
        let mut rejected = 0;
        let mut evaluations = 0;
        let mut maximum_error: f64 = 0.0;
        let mut sums = [Sum::default(); 12];
        while elapsed < input.duration_s {
            poll(&mut cancelled)?;
            if attempts == controls.maximum_attempts {
                return Err(Error::StepBudget {
                    attempts,
                    accepted_intervals: accepted,
                });
            }
            attempts += 1;
            let remaining = input.duration_s - elapsed;
            let chosen = proposal.min(remaining);
            let end = if chosen == remaining {
                input.duration_s
            } else {
                elapsed + chosen
            };
            let h = end - elapsed;
            let half = 0.5 * h;
            if end <= elapsed || !h.is_finite() || half <= 0.0 {
                return Err(Error::StepUnderflow);
            }
            let current = endpoint.as_ref().unwrap_or(self);
            let trial = (|| {
                evaluations += 1;
                let coarse = current.step_candidate(
                    SealedCapillaryInput {
                        duration_s: h,
                        ..input
                    },
                    solve,
                    &mut cancelled,
                )?;
                evaluations += 1;
                let first = current.step_candidate(
                    SealedCapillaryInput {
                        duration_s: half,
                        ..input
                    },
                    solve,
                    &mut cancelled,
                )?;
                evaluations += 1;
                let second = first.endpoint().step_candidate(
                    SealedCapillaryInput {
                        duration_s: h - half,
                        ..input
                    },
                    solve,
                    &mut cancelled,
                )?;
                let first_branch = branch(first.step());
                let second_branch = branch(second.step());
                let changed = branch(coarse.step()) != first_branch
                    || first_branch != second_branch
                    || previous_branch.is_some_and(|b| b != first_branch);
                let error = error_ratio(
                    current,
                    coarse.endpoint(),
                    second.endpoint(),
                    controls,
                    if changed { 1.0 } else { 3.0 },
                );
                Ok::<_, SealedCapillaryError>((first, second, error, changed, second_branch))
            })();
            let (first, second, error, changed, second_branch) = match trial {
                Ok(value) => value,
                Err(SealedCapillaryError::Cancelled) => return Err(Error::Cancelled),
                Err(error) if refinable(&error) && h > controls.minimum_step_s => {
                    rejected += 1;
                    proposal = shrink(h, elapsed, controls.minimum_step_s, 0.5)?;
                    continue;
                }
                Err(error) => {
                    return Err(Error::Step {
                        attempted_step_s: h,
                        error,
                    });
                }
            };
            if error.is_nan() {
                return Err(Error::Invalid("nonfinite local error estimate"));
            }
            if changed {
                previous_error = 1.0;
            }
            if error > 1.0 {
                rejected += 1;
                if h <= controls.minimum_step_s {
                    return Err(Error::Tolerance {
                        normalized_error: error,
                        step_s: h,
                    });
                }
                let factor = (pi.safety * error.powf(-1.0 / 3.0)).clamp(pi.max_shrink, 0.9);
                proposal = shrink(h, elapsed, controls.minimum_step_s, factor)?;
                continue;
            }
            poll(&mut cancelled)?;
            accumulate(&mut sums, first.step());
            accumulate(&mut sums, second.step());
            let (next, _) = second.into_parts();
            endpoint = Some(next);
            elapsed = end;
            accepted += 1;
            maximum_error = maximum_error.max(error);
            let safe_error = error.max(1e-300);
            let factor = pi.safety * safe_error.powf(-pi.k_p) * previous_error.powf(pi.k_i);
            if !factor.is_finite() {
                return Err(Error::Invalid("nonfinite midpoint PI proposal"));
            }
            // A stick/slip or reversal transition is nonsmooth. Discard old PI
            // history and do not immediately grow across that transition.
            proposal = (h * factor.clamp(pi.max_shrink, if changed { 1.0 } else { pi.max_growth }))
                .clamp(controls.minimum_step_s, controls.maximum_step_s);
            previous_error = if changed { 1.0 } else { safe_error };
            previous_branch = Some(second_branch);
        }
        let endpoint = endpoint.unwrap_or_else(|| self.clone());
        let v = sums.map(Sum::value);
        if v.iter().any(|x| !x.is_finite()) {
            return Err(Error::Invalid("aggregate physical transfers overflowed"));
        }
        let residual = ((v[2] + v[8]) - v[10]) - v[3] - v[4] - v[5] - v[6] - v[7];
        let scale =
            v[2].abs() + v[8].abs() + v[10].abs() + v[3].abs() + v[4].abs() + v[5] + v[6] + v[7];
        let tolerance = solve.energy_absolute_tolerance_j + solve.energy_relative_tolerance * scale;
        if !residual.is_finite() || !tolerance.is_finite() {
            return Err(Error::Invalid(
                "aggregate energy balance is unrepresentable",
            ));
        }
        if residual.abs() > tolerance {
            return Err(Error::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        poll(&mut cancelled)?;
        Ok(SealedCapillaryAdvanceCandidate {
            step: SealedCapillaryAdvanceStep {
                initial_wet_length_m: self.geometry().wet_length_m,
                wet_length_m: endpoint.geometry().wet_length_m,
                elapsed_time_s: elapsed,
                transferred_liquid_volume_m3: v[0],
                transferred_liquid_mass_kg: v[1],
                inlet_pressure_work_j: v[2],
                wall_surface_energy_change_j: v[3],
                gravitational_energy_change_j: v[4],
                liquid_viscous_dissipation_j: v[5],
                contact_line_viscous_dissipation_j: v[6],
                pinning_dissipation_j: v[7],
                gas_supplied_heat_j: v[8],
                gas_work_by_gas_j: v[9],
                gas_stored_energy_change_j: v[10],
                gas_entropy_change_j_k: v[11],
                energy_balance_residual_j: residual,
                attempts,
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                step_evaluations: evaluations,
                maximum_normalized_local_error: maximum_error,
            },
            endpoint,
        })
    }
}

fn validate(
    input: SealedCapillaryInput,
    solve: SealedCapillarySolve,
    c: SealedCapillaryAdvanceControls,
) -> Result<(), SealedCapillaryAdvanceError> {
    let invalid = SealedCapillaryAdvanceError::Invalid;
    if !input.duration_s.is_finite()
        || input.duration_s < 0.0
        || !input.inlet_pressure_pa.is_finite()
        || input.inlet_pressure_pa < 0.0
        || !input.gas_heat_power_w.is_finite()
        || !input.opposing_gravity_m_s2.is_finite()
        || ![
            c.initial_step_s,
            c.minimum_step_s,
            c.maximum_step_s,
            c.length_absolute_tolerance_m,
            c.temperature_absolute_tolerance_k,
            c.pressure_absolute_tolerance_pa,
        ]
        .into_iter()
        .all(|x| x.is_finite() && x > 0.0)
        || c.minimum_step_s > c.initial_step_s
        || c.initial_step_s > c.maximum_step_s
        || !c.relative_tolerance.is_finite()
        || !(0.0..=1.0).contains(&c.relative_tolerance)
        || c.maximum_attempts == 0
        || c.maximum_attempts > usize::MAX / 3
    {
        return Err(invalid(
            "finite boundary values, forward interval and positive adaptive scales required",
        ));
    }
    let heat = input.gas_heat_power_w * input.duration_s;
    if !heat.is_finite() || (input.gas_heat_power_w != 0.0 && input.duration_s > 0.0 && heat == 0.0)
    {
        return Err(invalid("requested interval heat is unrepresentable"));
    }
    solve
        .gas
        .validate()
        .map_err(|error| SealedCapillaryAdvanceError::Step {
            attempted_step_s: 0.0,
            error: SealedCapillaryError::Gas(error),
        })?;
    if !(1..=128).contains(&solve.maximum_iterations)
        || ![
            solve.pressure_absolute_tolerance_pa,
            solve.length_absolute_tolerance_m,
            solve.energy_absolute_tolerance_j,
        ]
        .into_iter()
        .all(|x| x.is_finite() && x > 0.0)
        || !solve.energy_relative_tolerance.is_finite()
        || solve.energy_relative_tolerance < 0.0
    {
        return Err(invalid("valid single-step solve controls required"));
    }
    Ok(())
}

fn error_ratio(
    initial: &SealedGasCapillary<'_>,
    coarse: &SealedGasCapillary<'_>,
    fine: &SealedGasCapillary<'_>,
    c: SealedCapillaryAdvanceControls,
    divisor: f64,
) -> f64 {
    let gas0 = initial.gas().state().state();
    let gas1 = coarse.gas().state().state();
    let gas2 = fine.gas().state().state();
    let mut error = 0.0_f64;
    for (old, large, small, atol) in [
        (
            initial.geometry().wet_length_m,
            coarse.geometry().wet_length_m,
            fine.geometry().wet_length_m,
            c.length_absolute_tolerance_m,
        ),
        (
            gas0.temperature,
            gas1.temperature,
            gas2.temperature,
            c.temperature_absolute_tolerance_k,
        ),
        (
            gas0.pressure,
            gas1.pressure,
            gas2.pressure,
            c.pressure_absolute_tolerance_pa,
        ),
    ] {
        let scale = atol + c.relative_tolerance * old.abs().max(small.abs());
        let delta = (small - large).abs() / divisor;
        if !scale.is_finite() || scale <= 0.0 || !delta.is_finite() {
            return f64::NAN;
        }
        error = error.max(delta / scale);
    }
    error
}
fn refinable(error: &SealedCapillaryError) -> bool {
    matches!(
        error,
        SealedCapillaryError::StepTooLarge(_)
            | SealedCapillaryError::NotConverged { .. }
            | SealedCapillaryError::EnergyBalance { .. }
            | SealedCapillaryError::Gas(
                GasEnergyError::OutsideBracket { .. } | GasEnergyError::NotConverged { .. }
            )
    )
}
fn shrink(
    h: f64,
    time: f64,
    minimum: f64,
    factor: f64,
) -> Result<f64, SealedCapillaryAdvanceError> {
    let next = (factor * h).max(minimum);
    // A nominally smaller proposal can round back to the rejected interval,
    // particularly at the minimum step. Never repeat that same physical trial.
    if next >= h || time + next == time || time + next >= time + h {
        Err(SealedCapillaryAdvanceError::StepUnderflow)
    } else {
        Ok(next)
    }
}
fn branch(step: &SealedCapillaryStep) -> ContactLineMotion {
    step.contact_line_motion.unwrap_or_else(|| {
        if step.transferred_liquid_volume_m3 > 0.0 {
            ContactLineMotion::Advancing
        } else if step.transferred_liquid_volume_m3 < 0.0 {
            ContactLineMotion::Receding
        } else {
            ContactLineMotion::Pinned
        }
    })
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SealedCapillaryAdvanceError> {
    if cancelled() {
        Err(SealedCapillaryAdvanceError::Cancelled)
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    total: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.total + value;
        self.correction += if self.total.abs() >= value.abs() {
            (self.total - next) + value
        } else {
            (value - next) + self.total
        };
        self.total = next;
    }
    fn value(self) -> f64 {
        self.total + self.correction
    }
}
fn accumulate(sums: &mut [Sum; 12], s: &SealedCapillaryStep) {
    let values = [
        s.transferred_liquid_volume_m3,
        s.transferred_liquid_mass_kg,
        s.inlet_pressure_work_j,
        s.wall_surface_energy_change_j,
        s.gravitational_energy_change_j,
        s.liquid_viscous_dissipation_j,
        s.contact_line_viscous_dissipation_j,
        s.pinning_dissipation_j,
        s.gas.supplied_heat_j,
        s.gas.work_by_gas_j,
        s.gas.stored_energy_change_j,
        s.gas.entropy_change_j_k,
    ];
    for (sum, value) in sums.iter_mut().zip(values) {
        sum.add(value);
    }
}

#[cfg(test)]
mod tests {
    use super::{SealedCapillaryAdvanceError, shrink};

    #[test]
    fn rounded_minimum_step_cannot_repeat_the_rejected_interval() {
        let time = 1.0;
        let minimum = 0.1;
        let interval = (time + minimum) - time;
        assert!(interval > minimum);
        assert_eq!(
            shrink(interval, time, minimum, 0.5),
            Err(SealedCapillaryAdvanceError::StepUnderflow),
        );
    }
}
