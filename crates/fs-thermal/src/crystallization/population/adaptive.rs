//! Local time refinement without collapsing a distribution of crystal ages.

use super::{CHECKPOINT_STRIDE, CrystallizationPopulation, exposure, invalid};
use crate::crystallization::{
    CrystallizationAdaptivity, CrystallizationError, CrystallizationSolve, check_cancelled,
    check_duration, check_solve, finite, positive,
};
use fs_math::clock::clock_endpoint;

/// Complete accepted population interval and its local time-resolution evidence.
///
/// Discrepancies are one full midpoint step versus two retained half steps.
/// Their maxima include accepted trials only, without Richardson extrapolation.
/// The fraction tolerance applies to EVERY cohort as well as the aggregate;
/// a low-mass but unresolved history cannot disappear in a weighted mean.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveCrystallizationPopulationStep {
    pub duration_s: f64,
    pub specific_heat_j_kg: f64,
    pub crystallized_mass_kg: f64,
    /// Original mixture enthalpy minus final stored enthalpy, J/kg.
    pub caloric_residual_j_kg: f64,
    pub accepted_substeps: usize,
    pub rejected_trials: usize,
    pub attempted_trials: usize,
    pub maximum_temperature_discrepancy_k: f64,
    /// Maximum over individual cohorts AND aggregate crystalline mass fraction.
    pub maximum_crystalline_mass_fraction_discrepancy: f64,
    /// Difference in common exposure increments, unless all cohorts saturate.
    pub maximum_exposure_discrepancy: f64,
}

impl CrystallizationPopulation<'_> {
    /// Resolve the original irreversible population on a linear enthalpy path.
    ///
    /// Every candidate delegates to `advance_specific_enthalpy`. Immutable
    /// cohorts remain shared, with no per-step cohort allocation, age averaging,
    /// remelting approximation or extra latent heat. Both retained half steps
    /// follow the same authored path and the final enthalpy is kept exactly.
    /// Temperature, every cohort's fraction, and the shared exposure increment
    /// have separate local error tests using the existing phase controls.
    ///
    /// Stage overshoots and local nonconvergence request bisection. Source and
    /// state failures propagate unchanged. Cancellation, exhausted trial/time
    /// budgets, or a domain failure publish no partially advanced population.
    /// These local discrepancies do not certify global error, outer transport,
    /// or phase-event timing. The maximum timestep must resolve source-rate
    /// features; finite samples cannot discover an arbitrary hidden narrow peak.
    pub fn advance_specific_enthalpy_adaptive(
        &self,
        final_h_j_kg: f64,
        duration_s: f64,
        solve: CrystallizationSolve,
        adaptivity: CrystallizationAdaptivity,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, AdaptiveCrystallizationPopulationStep), CrystallizationError> {
        check_cancelled(&mut cancelled)?;
        adaptivity.validate()?;
        check_duration(duration_s)?;
        check_solve(solve)?;
        if solve.temperature_tolerance_k > adaptivity.temperature_tolerance_k {
            return Err(invalid(
                "adaptive temperature tolerance undercuts caloric inversion tolerance",
            ));
        }
        self.validate(solve, &mut cancelled)?;
        let heat = finite(final_h_j_kg - self.specific_enthalpy_j_kg)?;
        // A valid accepted boundary state can carry either sign of caloric
        // residual. An unchanged account needs no new exact-range admission;
        // each original midpoint step still validates its physical trajectory.
        if heat != 0.0 {
            self.admit_adaptive_enthalpy(final_h_j_kg)?;
        }
        let mut state = self.clone();
        let mut report = AdaptiveCrystallizationPopulationStep {
            duration_s,
            specific_heat_j_kg: heat,
            crystallized_mass_kg: 0.0,
            caloric_residual_j_kg: 0.0,
            accepted_substeps: 0,
            rejected_trials: 0,
            attempted_trials: 0,
            maximum_temperature_discrepancy_k: 0.0,
            maximum_crystalline_mass_fraction_discrepancy: 0.0,
            maximum_exposure_discrepancy: 0.0,
        };
        let h_at = |time: f64| {
            let weight = time / duration_s;
            if time == duration_s {
                final_h_j_kg
            } else if weight <= 0.5 {
                self.specific_enthalpy_j_kg + weight * heat
            } else {
                final_h_j_kg - (1.0 - weight) * heat
            }
        };
        let mut elapsed = 0.0;
        let mut requested_dt = adaptivity.maximum_dt_s.min(duration_s);
        while elapsed < duration_s {
            check_cancelled(&mut cancelled)?;
            if report.attempted_trials == adaptivity.maximum_trials {
                return Err(CrystallizationError::NonConverged);
            }
            // Each accepted trial is one rounded clock addition; a leftover
            // within their rounding bound joins this trial (`fs_math::clock`).
            let endpoint = clock_endpoint(
                0.0,
                elapsed,
                requested_dt,
                duration_s,
                report.attempted_trials - report.rejected_trials,
            );
            let dt = endpoint - elapsed;
            let midpoint = elapsed + 0.5 * dt;
            if !positive(dt) || midpoint <= elapsed || midpoint >= endpoint {
                return Err(CrystallizationError::NonConverged);
            }
            report.attempted_trials += 1;
            let trial = (|| {
                let (full, _) =
                    state.advance_specific_enthalpy(h_at(endpoint), dt, solve, &mut cancelled)?;
                let (first, _) = state.advance_specific_enthalpy(
                    h_at(midpoint),
                    midpoint - elapsed,
                    solve,
                    &mut cancelled,
                )?;
                let (refined, step) = first.advance_specific_enthalpy(
                    h_at(endpoint),
                    endpoint - midpoint,
                    solve,
                    &mut cancelled,
                )?;
                // Cooling can enter a growth window after every midpoint
                // sample. Agreement at zero growth then hides unresolved onset.
                if refined.exposure_increment == state.exposure_increment
                    && !refined.fully_crystallized
                    && refined.temperature_k < self.law.config().melting_temperature_k
                    && (self.law.rate_s_inv(state.temperature_k)? > 0.0
                        || self.law.rate_s_inv(refined.temperature_k)? > 0.0)
                {
                    return Err(CrystallizationError::NonConverged);
                }
                let temperature_error = (full.temperature_k - refined.temperature_k).abs();
                let fraction_error = full.fraction_discrepancy(&refined, &mut cancelled)?;
                let exposure_error = if full.fully_crystallized && refined.fully_crystallized {
                    0.0
                } else {
                    (full.exposure_increment - refined.exposure_increment).abs()
                };
                Ok((
                    refined,
                    step,
                    temperature_error,
                    fraction_error,
                    exposure_error,
                ))
            })();
            let (refined, step, temperature_error, fraction_error, exposure_error) = match trial {
                Ok(result) => result,
                Err(
                    error @ (CrystallizationError::NonConverged
                    | CrystallizationError::EnthalpyOutOfRange { .. }),
                ) => {
                    report.rejected_trials += 1;
                    requested_dt = halve(dt, adaptivity.minimum_dt_s, error)?;
                    continue;
                }
                Err(error) => return Err(error),
            };
            let relative_error = (temperature_error / adaptivity.temperature_tolerance_k)
                .max(fraction_error / adaptivity.crystalline_mass_fraction_tolerance)
                .max(exposure_error / adaptivity.exposure_tolerance);
            if relative_error > 1.0 {
                report.rejected_trials += 1;
                requested_dt = halve(
                    dt,
                    adaptivity.minimum_dt_s,
                    CrystallizationError::NonConverged,
                )?;
                continue;
            }
            state = refined;
            report.caloric_residual_j_kg = step.caloric_residual_j_kg;
            report.accepted_substeps += 2;
            report.maximum_temperature_discrepancy_k = report
                .maximum_temperature_discrepancy_k
                .max(temperature_error);
            report.maximum_crystalline_mass_fraction_discrepancy = report
                .maximum_crystalline_mass_fraction_discrepancy
                .max(fraction_error);
            report.maximum_exposure_discrepancy =
                report.maximum_exposure_discrepancy.max(exposure_error);
            elapsed = endpoint;
            requested_dt = if relative_error <= 0.125 {
                if dt >= 0.5 * adaptivity.maximum_dt_s {
                    adaptivity.maximum_dt_s
                } else {
                    2.0 * dt
                }
            } else {
                dt.min(adaptivity.maximum_dt_s)
            };
        }
        report.crystallized_mass_kg = finite(
            self.mass_kg * (state.crystalline_mass_fraction - self.crystalline_mass_fraction),
        )?;
        check_cancelled(&mut cancelled)?;
        Ok((state, report))
    }

    fn admit_adaptive_enthalpy(&self, final_h: f64) -> Result<(), CrystallizationError> {
        let config = self.law.config();
        let [lower, upper] = config.temperature_range_k;
        let largest_fraction = if lower <= config.melting_temperature_k {
            config.maximum_crystalline_mass_fraction
        } else {
            0.0
        };
        let maximum_temperature = if self.crystalline_mass_fraction > 0.0 {
            upper.min(config.melting_temperature_k)
        } else {
            upper
        };
        let minimum = self.law.specific_enthalpy_j_kg(lower, largest_fraction)?;
        let maximum = self
            .law
            .specific_enthalpy_j_kg(maximum_temperature, self.crystalline_mass_fraction)?;
        if final_h < minimum || final_h > maximum {
            return Err(CrystallizationError::EnthalpyOutOfRange {
                specific_enthalpy_j_kg: final_h,
                minimum_j_kg: minimum,
                maximum_j_kg: maximum,
            });
        }
        Ok(())
    }

    fn fraction_discrepancy(
        &self,
        other: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, CrystallizationError> {
        let mut maximum = (self.crystalline_mass_fraction - other.crystalline_mass_fraction).abs();
        // Candidates share precisely the same original cohort buffer. Query
        // only the forward fraction, without allocating states or caloric data.
        for (index, cohort) in self.cohorts.iter().enumerate() {
            if index % CHECKPOINT_STRIDE == 0 {
                check_cancelled(cancelled)?;
            }
            let full = self
                .law
                .crystalline_mass_fraction(exposure(cohort, self.exposure_increment)?)?;
            let refined = self
                .law
                .crystalline_mass_fraction(exposure(cohort, other.exposure_increment)?)?;
            maximum = maximum.max((full - refined).abs());
        }
        Ok(maximum)
    }
}

fn halve(dt: f64, minimum: f64, error: CrystallizationError) -> Result<f64, CrystallizationError> {
    let half = 0.5 * dt;
    if !positive(half) || half < minimum {
        Err(error)
    } else {
        Ok(half)
    }
}
