//! Scale-safe inversion of one positive affine conductivity segment.

/// Invert a positive affine segment for its traversed temperature fraction.
/// Scaling both conductivity and energy prevents k² overflow, and the
/// rationalized positive root avoids subtracting nearly equal roots.
pub(super) fn inverse_fraction(
    amount: f64,
    capacity: f64,
    k_start: f64,
    k_end: f64,
) -> Result<f64, &'static str> {
    let scale = k_start.max(k_end);
    let a = k_start / scale;
    let b = k_end / scale;
    let ratio = amount / capacity;
    if !ratio.is_finite() || ratio <= 0.0 || ratio > 1.0 {
        return Err("scaled segment energy is not finite and positive");
    }
    // The convex form cannot suffer cancellation near the end of a decreasing
    // segment. Avoid halving ratio: the smallest positive energy fraction can
    // still have a representable square-root-sized temperature fraction.
    let discriminant = (1.0 - ratio) * (a * a) + ratio * (b * b);
    if !discriminant.is_finite() || discriminant < 0.0 {
        return Err("conductivity inverse discriminant is not representable");
    }
    let denominator = a + fs_math::det::sqrt(discriminant);
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err("conductivity inverse denominator is not finite and positive");
    }
    let fraction = (ratio * (a + b)) / denominator;
    if !fraction.is_finite() || fraction <= 0.0 || fraction > 1.0 {
        return Err("conductivity inverse left its admitted segment");
    }
    Ok(fraction)
}
