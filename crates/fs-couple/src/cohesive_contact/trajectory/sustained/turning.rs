//! Angular resolution of original two-dimensional sliding responses.
//!
//! These helpers classify actual endpoint velocities and choose a strictly
//! shorter duration. They do not supply motion, contact force, a stopping event,
//! or a continuous-flow error estimate. Every refined mechanical trial must
//! restart from the same accepted states and pass the unchanged physical audit.

/// Explicit resolution and work bounds for original sliding-direction changes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TurningControls {
    /// Smallest duration the original mechanical owner may evaluate, in s.
    pub minimum_step_s: f64,
    /// Maximum admitted resolved change in sliding direction, in radians.
    /// Strictly between zero and pi/2; no change to the Coulomb coefficient.
    pub maximum_angle_rad: f64,
    /// Number of new shortened mechanical trials, sharing the owner's budgets.
    pub maximum_refinements: usize,
}

/// Resolution receipt for an actual admitted endpoint, not a physical stop.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TurningReport {
    /// Original contact which first required angular refinement.
    pub trigger_contact_id: u64,
    /// Actual represented duration of the accepted original native trial.
    pub elapsed_s: f64,
    /// Largest measured endpoint angle among the relevant loaded contacts whose
    /// two endpoint speeds exceed the physical velocity tolerance.
    pub maximum_angle_rad: f64,
    /// All shortened trials, including ones rejected by another contact.
    pub refinements: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum TurningError {
    Invalid(&'static str),
    NonFinite(&'static str),
    MinimumStep {
        elapsed_s: f64,
        minimum_step_s: f64,
    },
    Refinement {
        refinements: usize,
        maximum: usize,
    },
    Unrepresentable {
        previous_s: f64,
        requested_s: f64,
        actual_s: f64,
    },
}

impl core::fmt::Display for TurningError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "invalid sliding-turn controls: {reason}"),
            Self::NonFinite(quantity) => write!(f, "non-finite sliding-turn {quantity}"),
            Self::MinimumStep {
                elapsed_s,
                minimum_step_s,
            } => write!(
                f,
                "sliding turn at {elapsed_s} s cannot refine below {minimum_step_s} s"
            ),
            Self::Refinement {
                refinements,
                maximum,
            } => write!(
                f,
                "sliding-turn refinement budget exhausted after {refinements}/{maximum} trials"
            ),
            Self::Unrepresentable {
                previous_s,
                requested_s,
                actual_s,
            } => write!(
                f,
                "sliding turn cannot shorten {previous_s} s to {requested_s} s (actual {actual_s} s)"
            ),
        }
    }
}
impl std::error::Error for TurningError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct DirectionChange {
    pub angle_rad: f64,
    /// False for a collinear reversal and for transverse components within the
    /// two endpoints' physical velocity bands. It never establishes a stop.
    pub transverse_resolved: bool,
}

pub(super) fn validate(controls: TurningControls) -> Result<(), TurningError> {
    if !positive(controls.minimum_step_s)
        || !positive(controls.maximum_angle_rad)
        || controls.maximum_angle_rad >= core::f64::consts::FRAC_PI_2
        || controls.maximum_refinements == 0
    {
        return Err(TurningError::Invalid(
            "positive minimum duration, angle below pi/2 and nonzero refinement budget required",
        ));
    }
    Ok(())
}

/// Measure two actual relative tangent velocities in their unchanged frame.
///
/// A near-rest endpoint has no resolved direction and returns `None`; the
/// original stopping and static-contact owners still decide its physical mode.
/// Normalization before products avoids overflowing finite speed products.
/// Both physical velocity bands contribute to the transverse-resolution test.
pub(super) fn measure(
    before: [f64; 2],
    after: [f64; 2],
    velocity_tolerance_m_s: f64,
) -> Result<Option<DirectionChange>, TurningError> {
    if !positive(velocity_tolerance_m_s) {
        return Err(TurningError::Invalid(
            "positive finite physical velocity tolerance required",
        ));
    }
    if before.iter().chain(&after).any(|v| !v.is_finite()) {
        return Err(TurningError::NonFinite("endpoint velocity"));
    }
    let speed0 = before[0].hypot(before[1]);
    let speed1 = after[0].hypot(after[1]);
    if !speed0.is_finite() || !speed1.is_finite() {
        return Err(TurningError::NonFinite("endpoint speed"));
    }
    if speed0 <= velocity_tolerance_m_s || speed1 <= velocity_tolerance_m_s {
        return Ok(None);
    }
    let u = [before[0] / speed0, before[1] / speed0];
    let v = [after[0] / speed1, after[1] / speed1];
    let cosine = u[0].mul_add(v[0], u[1] * v[1]).clamp(-1.0, 1.0);
    // Compensate the rounded second product. A lone FMA can otherwise leave
    // a false transverse residue even for exactly opposite normalized vectors
    // when the physical velocity band is smaller than machine relative error.
    let product = u[1] * v[0];
    let product_error = u[1].mul_add(v[0], -product);
    let sine = (u[0].mul_add(v[1], -product) - product_error).abs();
    let angle_rad = sine.atan2(cosine);
    let transverse_band = velocity_tolerance_m_s / speed0 + velocity_tolerance_m_s / speed1;
    Ok(Some(DirectionChange {
        angle_rad,
        transverse_resolved: sine > transverse_band,
    }))
}

/// Choose a bounded halving from the duration of an actual original response.
/// The receiving owner counts and evaluates the returned trial, then checks its
/// actual represented duration; rounding must not leave it at the old endpoint.
pub(super) fn shorter_duration(
    elapsed_s: f64,
    controls: TurningControls,
    completed_refinements: usize,
) -> Result<f64, TurningError> {
    validate(controls)?;
    if !positive(elapsed_s) {
        return Err(TurningError::Invalid(
            "positive finite original duration required",
        ));
    }
    if completed_refinements >= controls.maximum_refinements {
        return Err(TurningError::Refinement {
            refinements: completed_refinements,
            maximum: controls.maximum_refinements,
        });
    }
    if elapsed_s <= controls.minimum_step_s {
        return Err(TurningError::MinimumStep {
            elapsed_s,
            minimum_step_s: controls.minimum_step_s,
        });
    }
    let requested_s = (0.5 * elapsed_s).max(controls.minimum_step_s);
    if !(requested_s > 0.0 && requested_s < elapsed_s) {
        return Err(TurningError::Unrepresentable {
            previous_s: elapsed_s,
            requested_s,
            actual_s: requested_s,
        });
    }
    Ok(requested_s)
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

#[cfg(test)]
#[path = "turning/tests.rs"]
mod tests;
