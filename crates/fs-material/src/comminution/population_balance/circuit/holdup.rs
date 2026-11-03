//! Close the continuous classifier/recycle circuit at a specified dry hold-up.
//!
//! For fixed selection frequencies, daughter fractions and classifier splits,
//! each class inventory has the triangular form
//! `H_i = tau/(escape_i + tau*S_i) * (fresh_i + sum_j b_ij*S_j*H_j)`.
//! Every factor `tau/(escape + tau*S)` is nondecreasing for nonnegative
//! coefficients, so the total dry inventory is nondecreasing. A fed circuit
//! must have an outlet; the original circuit owner refuses trapped material.
//! This does not assume that every positive inventory is reachable.
//!
//! Bisection uses an explicit residence bracket and a finite iteration budget.
//! Both endpoints and every trial run the original circuit, including its
//! actual circulating feed, classifier discharge, source domains and balances.
//! Any refusal stops the solve; no bracket expansion or source-domain fallback
//! occurs. The returned inventory is the original model result, never replaced
//! by the requested target. The mass tolerance is a binary64 stopping condition,
//! not a certified enclosure or a calibration claim.

use super::{ContinuousCircuitError, ContinuousCircuitSolution, ContinuousGrindingCircuit};
use crate::comminution::circuit::CircuitSettings;
use crate::separation::{ParticleStream, Partition};
use fs_matdb::QueryPoint;
use fs_qty::{Dimensionless, Mass, Time};

/// Caller-declared bracket, dry-mass accuracy and maximum number of midpoint trials.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HoldupSolveSettings {
    /// Strictly positive lower residence, evaluated through the original owners.
    pub minimum_residence: Time,
    /// Finite upper residence strictly greater than the lower residence.
    pub maximum_residence: Time,
    /// Strictly positive absolute tolerance on dry inventory, in kg.
    pub absolute_mass_tolerance: Mass,
    /// Strictly positive relative tolerance, less than one, relative to the target.
    pub relative_mass_tolerance: Dimensionless,
    /// Positive limit on midpoint trials, in addition to the two endpoint trials.
    pub maximum_iterations: usize,
}

/// Most recent endpoint evaluations retained by the bracketed solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HoldupBracket {
    /// Evaluated lower residence.
    pub lower_residence: Time,
    /// Original dry hold-up at the lower residence.
    pub lower_holdup: Mass,
    /// Evaluated upper residence.
    pub upper_residence: Time,
    /// Original dry hold-up at the upper residence.
    pub upper_holdup: Mass,
}

/// A rejected target, source/circuit refusal or exhausted numerical solve.
#[derive(Clone, Debug, PartialEq)]
pub enum HoldupError {
    /// Invalid target, bracket, tolerances or iteration budget.
    Invalid(&'static str),
    /// An original circuit owner refused an endpoint or midpoint trial.
    Circuit {
        /// Actual residence passed to the original circuit.
        residence: Time,
        /// Unmodified original refusal.
        error: ContinuousCircuitError,
    },
    /// The requested mass lies outside the evaluated endpoint masses.
    Unbracketed {
        /// Requested dry inventory.
        target: Mass,
        /// Both original endpoint evaluations.
        bracket: HoldupBracket,
    },
    /// The caller's midpoint budget ended without meeting the mass tolerance.
    IterationExhausted {
        /// Requested dry inventory.
        target: Mass,
        /// Final evaluated bounds.
        bracket: HoldupBracket,
        /// Number of midpoint evaluations performed.
        iterations: usize,
    },
    /// Binary64 cannot refine the bracket or preserve its monotone mass ordering.
    NumericalResolution {
        /// Concrete arithmetic condition that prevented further resolution.
        reason: &'static str,
        /// Requested dry inventory.
        target: Mass,
        /// Last evaluated bounds, with no synthetic interpolation.
        bracket: HoldupBracket,
        /// Number of midpoint evaluations performed.
        iterations: usize,
    },
}
impl core::fmt::Display for HoldupError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous circuit dry hold-up: {self:?}")
    }
}
impl std::error::Error for HoldupError {}

/// An original admitted circuit solution satisfying the declared dry-mass tolerance.
#[derive(Clone, Debug, PartialEq)]
pub struct HoldupSolution {
    /// Unmodified original circuit response at the accepted residence.
    pub solution: ContinuousCircuitSolution,
    /// Requested dry inventory; carrier water is not included.
    pub target_dry_holdup: Mass,
    /// Signed original achieved dry inventory minus the requested inventory.
    pub mass_residual: Mass,
    /// Finite `absolute + relative * target` stopping tolerance in kg.
    pub allowed_mass_error: Mass,
    /// Last evaluated endpoint bracket containing the requested inventory.
    pub bracket: HoldupBracket,
    /// Midpoint evaluations performed, excluding the two endpoint evaluations.
    pub iterations: usize,
    /// Original caller-declared numerical settings.
    pub settings: HoldupSolveSettings,
}

impl ContinuousGrindingCircuit {
    /// Solve residence for a prescribed dry inventory of the same closed circuit.
    ///
    /// The target must lie inside the original evaluated endpoint masses, even
    /// if an outside target would be within the declared tolerance of an endpoint.
    /// Both endpoints must be admitted before an endpoint can be returned.
    /// Midpoint success depends on dry-mass residual, not a narrow time interval.
    /// All inputs remain untouched if any trial or numerical condition refuses.
    pub fn at_dry_holdup(
        &self,
        fresh: &ParticleStream,
        target: Mass,
        mill_conditions: &QueryPoint,
        classifier_conditions: &QueryPoint,
        circuit_settings: CircuitSettings,
        underflow_recycle: Partition,
        settings: HoldupSolveSettings,
    ) -> Result<HoldupSolution, HoldupError> {
        let allowed_error = settings.validate(target)?;
        let evaluate = |residence| {
            self.at_residence_time(
                fresh,
                residence,
                mill_conditions,
                classifier_conditions,
                circuit_settings,
                underflow_recycle,
            )
            .map_err(|error| HoldupError::Circuit { residence, error })
        };
        let lower = evaluate(settings.minimum_residence)?;
        let upper = evaluate(settings.maximum_residence)?;
        let mut bracket = HoldupBracket {
            lower_residence: settings.minimum_residence,
            lower_holdup: lower.mill.dry_solids_holdup,
            upper_residence: settings.maximum_residence,
            upper_holdup: upper.mill.dry_solids_holdup,
        };
        let resolution = |reason, bracket, iterations| HoldupError::NumericalResolution {
            reason,
            target,
            bracket,
            iterations,
        };
        if bracket.lower_holdup > bracket.upper_holdup {
            return Err(resolution(
                "endpoint dry masses are not monotone",
                bracket,
                0,
            ));
        }
        if target < bracket.lower_holdup || target > bracket.upper_holdup {
            return Err(HoldupError::Unbracketed { target, bracket });
        }
        let accepts = |actual: Mass| (actual.value() - target.value()).abs() <= allowed_error;
        let finish = |solution: ContinuousCircuitSolution, bracket, iterations| HoldupSolution {
            mass_residual: Mass::new(solution.mill.dry_solids_holdup.value() - target.value()),
            solution,
            target_dry_holdup: target,
            allowed_mass_error: Mass::new(allowed_error),
            bracket,
            iterations,
            settings,
        };
        if accepts(bracket.lower_holdup) {
            return Ok(finish(lower, bracket, 0));
        }
        if accepts(bracket.upper_holdup) {
            return Ok(finish(upper, bracket, 0));
        }
        for iteration in 1..=settings.maximum_iterations {
            let lo = bracket.lower_residence.value();
            let hi = bracket.upper_residence.value();
            let middle = lo + (hi - lo) * 0.5;
            if middle <= lo || middle >= hi {
                return Err(resolution(
                    "no representable interior residence remains",
                    bracket,
                    iteration - 1,
                ));
            }
            let residence = Time::new(middle);
            let solution = evaluate(residence)?;
            let mass = solution.mill.dry_solids_holdup;
            if mass < bracket.lower_holdup || mass > bracket.upper_holdup {
                return Err(resolution(
                    "trial dry mass is outside its monotone endpoint bounds",
                    bracket,
                    iteration,
                ));
            }
            if accepts(mass) {
                return Ok(finish(solution, bracket, iteration));
            }
            if mass < target {
                bracket.lower_residence = residence;
                bracket.lower_holdup = mass;
            } else {
                bracket.upper_residence = residence;
                bracket.upper_holdup = mass;
            }
        }
        Err(HoldupError::IterationExhausted {
            target,
            bracket,
            iterations: settings.maximum_iterations,
        })
    }
}

impl HoldupSolveSettings {
    fn validate(self, target: Mass) -> Result<f64, HoldupError> {
        if !target.value().is_finite() || target.value() <= 0.0 {
            return Err(HoldupError::Invalid(
                "target dry hold-up must be finite and positive",
            ));
        }
        let lo = self.minimum_residence.value();
        let hi = self.maximum_residence.value();
        if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || hi <= lo {
            return Err(HoldupError::Invalid(
                "residence bracket must be finite, positive and increasing",
            ));
        }
        let absolute = self.absolute_mass_tolerance.value();
        let relative = self.relative_mass_tolerance.value();
        if !absolute.is_finite()
            || absolute <= 0.0
            || !relative.is_finite()
            || relative <= 0.0
            || relative >= 1.0
        {
            return Err(HoldupError::Invalid(
                "mass tolerances must be finite and positive, with relative tolerance below one",
            ));
        }
        if self.maximum_iterations == 0 {
            return Err(HoldupError::Invalid(
                "midpoint iteration budget must be positive",
            ));
        }
        let allowed_error = relative.mul_add(target.value(), absolute);
        if !allowed_error.is_finite() {
            return Err(HoldupError::Invalid(
                "combined dry-mass tolerance overflows",
            ));
        }
        Ok(allowed_error)
    }
}

#[cfg(test)]
#[path = "holdup/tests.rs"]
mod tests;
