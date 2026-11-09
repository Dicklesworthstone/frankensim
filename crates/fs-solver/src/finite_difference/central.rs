//! Bounded central directional actions on an original anchored source.
//!
//! Sources receive offsets from an immutable reference, not a new physical
//! state representation. Two successive central stencils must agree before
//! publication. Their disagreement is a numerical diagnostic, NOT a certified
//! derivative error or a proof of smoothness; the source owns domain admission.

use fs_math::dd::Dd;

/// Scales and stopping controls in the caller's explicit input/output basis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CentralDifferenceControls {
    /// Initial largest absolute offset component.
    pub initial_step: f64,
    /// Allowed componentwise disagreement between successive derivative actions.
    pub absolute_tolerance: f64,
    /// Additional allowance relative to the larger of the two actions.
    pub relative_tolerance: f64,
    /// Successive halvings in 1..=16. Each stencil queries both signs.
    pub maximum_refinements: usize,
    /// All attempted source calls, including a refused call.
    pub maximum_evaluations: usize,
}
impl Default for CentralDifferenceControls {
    fn default() -> Self {
        Self { initial_step: 1e-3, absolute_tolerance: 1e-7, relative_tolerance: 1e-6,
            maximum_refinements: 8, maximum_evaluations: 18 }
    }
}

/// Original refusal or an unresolved numerical stencil. No partial result exists.
#[derive(Clone, Debug, PartialEq)]
pub enum CentralDifferenceError<E> {
    Invalid(&'static str),
    Source(E),
    EvaluationBudget { used: usize, maximum: usize },
    Unresolved { refinements: usize, maximum_disagreement: f64 },
}
impl<E: core::fmt::Debug> core::fmt::Display for CentralDifferenceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "central directional difference: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for CentralDifferenceError<E> {}

/// Refined action and its actual two-stencil numerical disagreement.
#[derive(Clone, Debug, PartialEq)]
pub struct CentralDifference {
    pub action: Vec<f64>,
    pub disagreement: Vec<f64>,
    pub evaluations: usize,
    pub refinements: usize,
    pub final_step: f64,
}

/// Evaluate d/ds F(s*direction) at s=0 in a caller-owned reference chart.
///
/// The largest offset starts at initial_step, independently of the magnitude
/// of the requested direction. Components are normalized by division, retaining
/// subnormal directions without an overflowing reciprocal. Unrepresentable
/// nonzero offsets refuse rather than silently dropping a coordinate.
/// Source calls must be repeatable, include their own physical/chart checks and
/// preserve accepted history. No clipping, one-sided fallback, extrapolation or
/// source-error replacement is performed. Zero direction makes no source calls
/// but still validates controls and polls cancellation before publication.
pub fn differentiate<E>(
    direction: &[f64],
    output_dimension: usize,
    controls: CentralDifferenceControls,
    mut source: impl FnMut(&[f64]) -> Result<Vec<f64>, E>,
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<CentralDifference, CentralDifferenceError<E>> {
    checkpoint().map_err(CentralDifferenceError::Source)?;
    if direction.is_empty() || direction.len() > 1_048_576
        || output_dimension == 0 || output_dimension > 1_048_576
        || direction.iter().any(|v| !v.is_finite())
        || !controls.initial_step.is_finite() || controls.initial_step <= 0.0
        || !controls.absolute_tolerance.is_finite() || controls.absolute_tolerance < 0.0
        || !controls.relative_tolerance.is_finite() || controls.relative_tolerance < 0.0
        || !(1..=16).contains(&controls.maximum_refinements)
        || !(4..=128).contains(&controls.maximum_evaluations)
    {
        return Err(CentralDifferenceError::Invalid("finite bounded dimensions, scales and controls required"));
    }
    let scale = direction.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
    if scale == 0.0 {
        checkpoint().map_err(CentralDifferenceError::Source)?;
        return Ok(CentralDifference { action: vec![0.0; output_dimension],
            disagreement: vec![0.0; output_dimension], evaluations: 0, refinements: 0,
            final_step: 0.0 });
    }
    let mut used = 0;
    let mut offset = vec![0.0; direction.len()];
    let mut query = |step: f64| -> Result<Vec<f64>, CentralDifferenceError<E>> {
        checkpoint().map_err(CentralDifferenceError::Source)?;
        if used >= controls.maximum_evaluations {
            return Err(CentralDifferenceError::EvaluationBudget {
                used, maximum: controls.maximum_evaluations,
            });
        }
        for (value, &d) in offset.iter_mut().zip(direction) {
            *value = (d / scale) * step;
            if !value.is_finite() || (d != 0.0 && *value == 0.0) {
                return Err(CentralDifferenceError::Invalid("nonzero offset is not representable"));
            }
        }
        used += 1;
        let result = source(&offset).map_err(CentralDifferenceError::Source)?;
        if result.len() != output_dimension || result.iter().any(|v| !v.is_finite()) {
            return Err(CentralDifferenceError::Invalid("source returned incomplete or nonfinite output"));
        }
        Ok(result)
    };
    let mut previous: Option<Vec<f64>> = None;
    let mut step = controls.initial_step;
    let mut maximum_disagreement = f64::INFINITY;
    for refinement in 0..=controls.maximum_refinements {
        let plus = query(step)?;
        let minus = query(-step)?;
        let mut action = Vec::with_capacity(output_dimension);
        for (&a, &b) in plus.iter().zip(&minus) {
            let value = ((Dd::from_f64(a) - Dd::from_f64(b))
                / Dd::from_f64(step) * Dd::from_f64(0.5) * Dd::from_f64(scale)).to_f64();
            if !value.is_finite() || (a != b && value == 0.0) {
                return Err(CentralDifferenceError::Invalid("derivative action is not representable"));
            }
            action.push(value);
        }
        if let Some(old) = previous {
            let mut accepted = true;
            let mut disagreement = Vec::with_capacity(output_dimension);
            maximum_disagreement = 0.0;
            for (&a, &b) in action.iter().zip(&old) {
                let error = (a-b).abs();
                maximum_disagreement = maximum_disagreement.max(error);
                let allowance = controls.absolute_tolerance
                    + controls.relative_tolerance * a.abs().max(b.abs());
                if !allowance.is_finite() || !error.is_finite() || error > allowance {
                    accepted = false;
                }
                disagreement.push(error);
            }
            if accepted {
                drop(query);
                checkpoint().map_err(CentralDifferenceError::Source)?;
                return Ok(CentralDifference { action, disagreement, evaluations: used,
                    refinements: refinement, final_step: step });
            }
        }
        previous = Some(action);
        let next = 0.5 * step;
        if next == 0.0 || next == step {
            return Err(CentralDifferenceError::Invalid("refined stencil is not representable"));
        }
        step = next;
    }
    Err(CentralDifferenceError::Unresolved {
        refinements: controls.maximum_refinements, maximum_disagreement,
    })
}

#[cfg(test)]
#[path = "central_tests.rs"]
mod tests;
