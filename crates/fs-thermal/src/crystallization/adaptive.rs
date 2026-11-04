//! Bounded step doubling on the original crystallization and remelting kernels.
//!
//! Every accepted interval uses the two half steps without extrapolating phase
//! fraction, age or enthalpy. The whole interval and its halves follow the same
//! prescribed linear enthalpy path. Their differences are local time-resolution
//! estimates, not certified global errors or event-time/order guarantees.
//! The caller's maximum timestep must resolve the original kinetic source: no
//! finite sample set can detect every narrow peak of an arbitrary rate closure.

use super::{
    CrystallizationError, CrystallizationLaw, CrystallizationSolve, CrystallizationState,
    CrystallizationStep, EquivalentAgeRemelting, check_cancelled, check_duration, check_solve,
    finite, positive,
};

/// Explicit local time-accuracy and work controls, independent of caloric inversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationAdaptivity {
    /// Upper bound on each trial interval. This is also the caller's resolution
    /// bound for features of the original temperature-dependent kinetic source.
    pub maximum_dt_s: f64,
    /// A rejected interval may not be bisected below this duration. A shorter
    /// final remainder can still be accepted if its first trial satisfies all tests.
    pub minimum_dt_s: f64,
    /// Counts accepted AND rejected trials; each invokes at most three original
    /// midpoint steps, each with the separate bounded caloric solve controls.
    pub maximum_trials: usize,
    pub temperature_tolerance_k: f64,
    /// Absolute tolerance on crystalline MASS fraction, in (0,1].
    pub crystalline_mass_fraction_tolerance: f64,
    /// Dimensionless history tolerance. It resolves age even when the forward
    /// Avrami fraction is initially very small. Two exactly saturated states
    /// need no age comparison: their unused excess exposure cannot affect X.
    pub exposure_tolerance: f64,
}

impl Default for CrystallizationAdaptivity {
    fn default() -> Self {
        Self {
            maximum_dt_s: 1.0,
            minimum_dt_s: 1e-10,
            maximum_trials: 10_000,
            temperature_tolerance_k: 1e-4,
            crystalline_mass_fraction_tolerance: 1e-7,
            exposure_tolerance: 1e-7,
        }
    }
}

impl CrystallizationAdaptivity {
    /// Validate without advancing or evaluating a material source.
    pub fn validate(self) -> Result<(), CrystallizationError> {
        if !positive(self.maximum_dt_s)
            || !positive(self.minimum_dt_s)
            || self.minimum_dt_s > self.maximum_dt_s
            || self.maximum_trials == 0
            || self.maximum_trials > usize::MAX / 2
            || !positive(self.temperature_tolerance_k)
            || !positive(self.crystalline_mass_fraction_tolerance)
            || self.crystalline_mass_fraction_tolerance > 1.0
            || !positive(self.exposure_tolerance)
        {
            return Err(CrystallizationError::Invalid(
                "invalid bounded crystallization adaptivity controls",
            ));
        }
        Ok(())
    }
}

/// The complete accepted enthalpy interval and its local resolution evidence.
///
/// Discrepancy maxima include accepted trials only. They are direct differences
/// between one step and two halves, with no smooth-regime Richardson factor.
/// The supplied caloric tolerances remain independent error contributions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveCrystallizationStep {
    pub initial: CrystallizationState,
    pub state: CrystallizationState,
    pub duration_s: f64,
    pub specific_heat_j_kg: f64,
    pub caloric_residual_j_kg: f64,
    /// The two retained half steps of each accepted trial.
    pub accepted_substeps: usize,
    pub rejected_trials: usize,
    pub attempted_trials: usize,
    pub maximum_temperature_discrepancy_k: f64,
    pub maximum_crystalline_mass_fraction_discrepancy: f64,
    pub maximum_exposure_discrepancy: f64,
}

impl CrystallizationLaw<'_> {
    /// Follow a prescribed linear enthalpy path with bounded local step doubling.
    ///
    /// `None` retains the original irreversible phase law. `Some` selects the
    /// already admitted equilibrium-remelting and equivalent-age model. All
    /// candidate steps use those same original kernels and immutable sources.
    /// No phase, temperature or enthalpy is clipped, and no latent heater is added.
    ///
    /// Local `NonConverged` and enthalpy-stage overshoots request bisection.
    /// Invalid source data, temperature-domain errors, state mismatch and
    /// cancellation propagate immediately. Exhausting `maximum_trials`, the
    /// minimum timestep, or represented time resolution refuses the entire
    /// interval; no partial state is published. The latter two retain a stage's
    /// specific enthalpy-domain error when one is available.
    ///
    /// The caller must resolve caloric error separately. An adaptive temperature
    /// tolerance tighter than the caloric inversion tolerance is refused. These
    /// local estimates do not certify the outer heat-transfer discretization,
    /// global time error, phase-event accuracy or experimental model validity.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_specific_enthalpy_adaptive(
        &self,
        current: CrystallizationState,
        final_h_j_kg: f64,
        duration_s: f64,
        solve: CrystallizationSolve,
        adaptivity: CrystallizationAdaptivity,
        remelting: Option<EquivalentAgeRemelting>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveCrystallizationStep, CrystallizationError> {
        check_cancelled(&mut cancelled)?;
        adaptivity.validate()?;
        check_duration(duration_s)?;
        check_solve(solve)?;
        if solve.temperature_tolerance_k > adaptivity.temperature_tolerance_k {
            return Err(CrystallizationError::Invalid(
                "adaptive temperature tolerance undercuts caloric inversion tolerance",
            ));
        }
        self.validate_state(current, solve)?;
        if let Some(model) = remelting {
            if self.equivalent_age_remelting()? != model {
                return Err(CrystallizationError::StateMismatch);
            }
        }
        let heat = finite(final_h_j_kg - current.specific_enthalpy_j_kg())?;
        self.check_adaptive_enthalpy_range(current, final_h_j_kg, remelting.is_some())?;
        let mut result = AdaptiveCrystallizationStep {
            initial: current,
            state: current,
            duration_s,
            specific_heat_j_kg: heat,
            caloric_residual_j_kg: 0.0,
            accepted_substeps: 0,
            rejected_trials: 0,
            attempted_trials: 0,
            maximum_temperature_discrepancy_k: 0.0,
            maximum_crystalline_mass_fraction_discrepancy: 0.0,
            maximum_exposure_discrepancy: 0.0,
        };
        let mut elapsed = 0.0;
        let mut requested_dt = adaptivity.maximum_dt_s.min(duration_s);
        while elapsed < duration_s {
            check_cancelled(&mut cancelled)?;
            if result.attempted_trials == adaptivity.maximum_trials {
                return Err(CrystallizationError::NonConverged);
            }
            let remaining = duration_s - elapsed;
            let endpoint = if requested_dt >= remaining {
                duration_s
            } else {
                elapsed + requested_dt
            };
            let dt = endpoint - elapsed;
            let midpoint = elapsed + 0.5 * dt;
            if !positive(dt) || midpoint <= elapsed || midpoint >= endpoint {
                return Err(CrystallizationError::NonConverged);
            }
            let h_at = |time: f64| {
                let weight = time / duration_s;
                // Retain the exact authored endpoint and avoid repeatedly adding
                // rounded local heats or multiplying a large formation datum.
                if time == duration_s {
                    final_h_j_kg
                } else if weight <= 0.5 {
                    current.specific_enthalpy_j_kg() + weight * heat
                } else {
                    final_h_j_kg - (1.0 - weight) * heat
                }
            };
            result.attempted_trials += 1;
            let trial = (|| {
                let full = self.adaptive_midpoint_step(
                    result.state,
                    h_at(endpoint),
                    dt,
                    solve,
                    remelting,
                    &mut cancelled,
                )?;
                let first = self.adaptive_midpoint_step(
                    result.state,
                    h_at(midpoint),
                    midpoint - elapsed,
                    solve,
                    remelting,
                    &mut cancelled,
                )?;
                let refined = self.adaptive_midpoint_step(
                    first.state,
                    h_at(endpoint),
                    endpoint - midpoint,
                    solve,
                    remelting,
                    &mut cancelled,
                )?;
                // During cooling or adiabatic recalescence, both trials may
                // overgrow and project to the very same plateau state. Their
                // zero discrepancy then hides an unresolved kinetic approach
                // to Tm. Resolve the remaining state excursion before accepting
                // such a projection. Positive heating can reach Tm directly;
                // it is not subjected to this non-heating approach guard.
                if remelting.is_some()
                    && h_at(endpoint) <= result.state.specific_enthalpy_j_kg()
                    && result.state.temperature_k() < self.config.melting_temperature_k
                    && refined.state.temperature_k() == self.config.melting_temperature_k
                    && ((refined.state.temperature_k() - result.state.temperature_k()).abs()
                        > adaptivity.temperature_tolerance_k
                        || (refined.state.crystalline_mass_fraction()
                            - result.state.crystalline_mass_fraction())
                        .abs()
                            > adaptivity.crystalline_mass_fraction_tolerance
                        || (refined.state.exposure() - result.state.exposure()).abs()
                            > adaptivity.exposure_tolerance)
                {
                    return Err(CrystallizationError::NonConverged);
                }
                // A cooling interval can enter growth after all three midpoint
                // samples. The two zero-age solutions then agree spuriously.
                // Positive endpoint kinetics (or lost undercooled initial growth)
                // require refinement, even when the numerical discrepancy is zero.
                if refined.state.exposure() == result.state.exposure()
                    && refined.state.crystalline_mass_fraction()
                        < self.config.maximum_crystalline_mass_fraction
                    && refined.state.temperature_k() < self.config.melting_temperature_k
                    && (self.rate_s_inv(result.state.temperature_k())? > 0.0
                        || self.rate_s_inv(refined.state.temperature_k())? > 0.0)
                {
                    return Err(CrystallizationError::NonConverged);
                }
                Ok((full, refined))
            })();
            let (full, refined) = match trial {
                Ok(pair) => pair,
                Err(
                    error @ (CrystallizationError::NonConverged
                    | CrystallizationError::EnthalpyOutOfRange { .. }),
                ) => {
                    result.rejected_trials += 1;
                    requested_dt = halve(dt, adaptivity.minimum_dt_s, error)?;
                    continue;
                }
                Err(error) => return Err(error),
            };
            let temperature_error =
                (full.state.temperature_k() - refined.state.temperature_k()).abs();
            let fraction_error = (full.state.crystalline_mass_fraction()
                - refined.state.crystalline_mass_fraction())
            .abs();
            let exposure_error = if full.state.crystalline_mass_fraction()
                == self.config.maximum_crystalline_mass_fraction
                && refined.state.crystalline_mass_fraction()
                    == self.config.maximum_crystalline_mass_fraction
            {
                0.0
            } else {
                (full.state.exposure() - refined.state.exposure()).abs()
            };
            let relative_error = (temperature_error / adaptivity.temperature_tolerance_k)
                .max(fraction_error / adaptivity.crystalline_mass_fraction_tolerance)
                .max(exposure_error / adaptivity.exposure_tolerance);
            if relative_error > 1.0 {
                result.rejected_trials += 1;
                requested_dt = halve(
                    dt,
                    adaptivity.minimum_dt_s,
                    CrystallizationError::NonConverged,
                )?;
                continue;
            }
            result.state = refined.state;
            result.caloric_residual_j_kg = refined.caloric_residual_j_kg;
            result.accepted_substeps += 2;
            result.maximum_temperature_discrepancy_k = result
                .maximum_temperature_discrepancy_k
                .max(temperature_error);
            result.maximum_crystalline_mass_fraction_discrepancy = result
                .maximum_crystalline_mass_fraction_discrepancy
                .max(fraction_error);
            result.maximum_exposure_discrepancy =
                result.maximum_exposure_discrepancy.max(exposure_error);
            elapsed = endpoint;
            requested_dt = if relative_error <= 0.125 {
                // Dyadic changes are deterministic and avoid an extra power
                // evaluation. Compare before doubling to avoid f64 overflow.
                if dt >= 0.5 * adaptivity.maximum_dt_s {
                    adaptivity.maximum_dt_s
                } else {
                    2.0 * dt
                }
            } else {
                dt.min(adaptivity.maximum_dt_s)
            };
        }
        check_cancelled(&mut cancelled)?;
        Ok(result)
    }

    fn adaptive_midpoint_step(
        &self,
        current: CrystallizationState,
        final_h_j_kg: f64,
        dt: f64,
        solve: CrystallizationSolve,
        remelting: Option<EquivalentAgeRemelting>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CrystallizationStep, CrystallizationError> {
        match remelting {
            Some(model) => self.advance_specific_enthalpy_with_remelting(
                current,
                final_h_j_kg,
                dt,
                solve,
                model,
                cancelled,
            ),
            None => self.advance_specific_enthalpy(current, final_h_j_kg, dt, solve, cancelled),
        }
    }

    fn check_adaptive_enthalpy_range(
        &self,
        current: CrystallizationState,
        final_h: f64,
        remelting: bool,
    ) -> Result<(), CrystallizationError> {
        let [lower, upper] = self.config.temperature_range_k;
        let largest_fraction = if lower <= self.config.melting_temperature_k {
            self.config.maximum_crystalline_mass_fraction
        } else {
            0.0
        };
        let smallest_fraction = if remelting {
            0.0
        } else {
            current.crystalline_mass_fraction()
        };
        let maximum_temperature = if smallest_fraction > 0.0 {
            upper.min(self.config.melting_temperature_k)
        } else {
            upper
        };
        let minimum = self.specific_enthalpy_j_kg(lower, largest_fraction)?;
        let maximum = self.specific_enthalpy_j_kg(maximum_temperature, smallest_fraction)?;
        if final_h < minimum || final_h > maximum {
            return Err(CrystallizationError::EnthalpyOutOfRange {
                specific_enthalpy_j_kg: final_h,
                minimum_j_kg: minimum,
                maximum_j_kg: maximum,
            });
        }
        Ok(())
    }
}

fn halve(
    duration_s: f64,
    minimum_dt_s: f64,
    error: CrystallizationError,
) -> Result<f64, CrystallizationError> {
    let half = 0.5 * duration_s;
    if half < minimum_dt_s || !positive(half) {
        Err(error)
    } else {
        Ok(half)
    }
}

#[cfg(test)]
#[path = "adaptive/tests.rs"]
mod tests;
