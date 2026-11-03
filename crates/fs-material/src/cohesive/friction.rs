//! Common tangent-disk algebra for damaged-area rigid Coulomb contact.
//!
//! The cohesive owner supplies its original path integral of damage times
//! compressive pressure. Receiving mechanical owners supply sourced static or
//! kinetic coefficients and solve `t = project(t + k delta_g, mu P)` jointly
//! with their native motion. The positive numerical scale k is not a spring:
//! the zero of this map is the unregularized rigid-friction inclusion.
//! No additional storage, slip state or material coefficient is inferred.

use super::CohesiveError;

/// Euclidean projection onto the closed two-dimensional friction disk.
/// Input and capacity must use the same units (traction or integrated force).
/// Preserves an interior input exactly; the zero-capacity disk returns zero.
pub fn project_tangent_disk(input: [f64; 2], capacity: f64) -> Result<[f64; 2], CohesiveError> {
    let length = input[0].hypot(input[1]);
    if !length.is_finite() || !capacity.is_finite() || capacity < 0.0 {
        return Err(CohesiveError::Invalid(
            "finite tangent effort and nonnegative Coulomb capacity required",
        ));
    }
    Ok(if length > capacity {
        input.map(|value| capacity * (value / length))
    } else {
        input
    })
}
