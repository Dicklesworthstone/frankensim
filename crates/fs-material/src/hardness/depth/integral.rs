//! Exact segment integration of a prepared affine hardness profile.
//!
//! At fixed pressure, differential Archard removal satisfies
//! `integral(H(z) dz) = k * p * ds` [Pa m]. This integral concerns exposed
//! material depth, not indentation depth or frictional heat. Pressure feedback
//! belongs to the contact integrator. No source query or allocation occurs here.

use fs_math::{dd::Dd, det::sqrt};

use super::{DepthHardnessError, ResolvedDepthHardness, require};

/// Retained exposed depth and the original profile's mean along its increment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DepthHardnessAdvance {
    /// New absolute exposed material depth [m].
    pub depth_m: f64,
    /// Actual representable difference from the supplied original depth [m].
    pub increment_m: f64,
    /// `integral(H dz) / increment` [Pa], or initial H for zero removal.
    /// This is an interval quadrature mean, not a new material property claim.
    pub path_mean_hardness_pa: f64,
}

impl ResolvedDepthHardness {
    /// Integrate the original piecewise-affine hardness between admitted depths.
    ///
    /// Units are Pa m (force per length), not energy or removed volume. The
    /// integral uses local segment coordinates and compensated accumulation;
    /// subtracting two large reference-surface primitives would lose small wear
    /// increments late in a history. Work is bounded by the prepared source's
    /// native 4096-knot ceiling. Equal endpoints return zero without removal.
    pub fn hardness_integral_between(
        &self,
        start_m: f64,
        end_m: f64,
    ) -> Result<f64, DepthHardnessError> {
        let mut hardness = self.hardness_pa_at(start_m)?;
        let final_hardness = self.hardness_pa_at(end_m)?;
        require(
            end_m >= start_m,
            "hardness integration requires ordered exposed depths",
        )?;
        if end_m == start_m {
            return Ok(0.0);
        }
        let mut cursor = start_m;
        let mut right = self.knots.partition_point(|knot| knot.0 <= cursor);
        let mut total = Dd::ZERO;
        while cursor < end_m {
            let next = self.knots[right].0.min(end_m);
            let next_hardness = if next == end_m {
                final_hardness
            } else {
                self.knots[right].1
            };
            total = total + segment_integral(next - cursor, hardness, next_hardness)?;
            require(
                finite_dd(total),
                "accumulated hardness integral is not representable",
            )?;
            cursor = next;
            hardness = next_hardness;
            right += 1;
        }
        let value = total.to_f64();
        require(
            value > 0.0,
            "positive hardness integral is below arithmetic resolution",
        )?;
        Ok(value)
    }

    /// Advance exposed material depth by a nonnegative hardness integral [Pa m].
    ///
    /// For an Archard interval the supplied measure is `k * p * ds`. Positive
    /// source hardness makes its integral strictly increasing. Only crossed
    /// original affine segments are visited, and the final segment is inverted
    /// analytically using a scaled, cancellation-resistant quadratic root.
    ///
    /// Source limits are never extrapolated. A positive measure that cannot
    /// produce a representable depth increment is refused. The returned mean
    /// comes from the ORIGINAL profile over the actual retained increment, so a
    /// consumer can retain its existing Archard volume and roundoff accounting.
    pub fn advance_by_integral(
        &self,
        depth_m: f64,
        hardness_integral_pa_m: f64,
    ) -> Result<DepthHardnessAdvance, DepthHardnessError> {
        let initial_hardness = self.hardness_pa_at(depth_m)?;
        require(
            hardness_integral_pa_m.is_finite() && hardness_integral_pa_m >= 0.0,
            "hardness-integral increment must be finite and nonnegative",
        )?;
        if hardness_integral_pa_m == 0.0 {
            return Ok(DepthHardnessAdvance {
                depth_m,
                increment_m: 0.0,
                path_mean_hardness_pa: initial_hardness,
            });
        }
        let target = Dd::from_f64(hardness_integral_pa_m);
        let mut consumed = Dd::ZERO;
        let mut cursor = depth_m;
        let mut start_hardness = initial_hardness;
        let mut right = self.knots.partition_point(|knot| knot.0 <= cursor);
        while right < self.knots.len() {
            let (end, end_hardness) = self.knots[right];
            let span = end - cursor;
            let capacity = segment_integral(span, start_hardness, end_hardness)?;
            let remaining = (target - consumed).to_f64();
            if remaining <= capacity.to_f64() {
                let candidate = if remaining == capacity.to_f64() {
                    end
                } else {
                    let increment = inverse_segment(
                        span,
                        start_hardness,
                        end_hardness,
                        remaining,
                        capacity.to_f64(),
                    )?;
                    // The input measure was admitted inside this segment. This
                    // bracket is for rounding of its root, never extrapolation.
                    (cursor + increment).min(end)
                };
                let increment_m = candidate - depth_m;
                require(
                    candidate.is_finite() && candidate > depth_m && increment_m > 0.0,
                    "positive wear is below retained exposed-depth resolution",
                )?;
                let integral = if candidate == cursor {
                    consumed
                } else {
                    consumed
                        + segment_integral(
                            candidate - cursor,
                            start_hardness,
                            self.hardness_pa_at(candidate)?,
                        )?
                }
                .to_f64();
                let mean = integral / increment_m;
                require(
                    mean.is_finite() && mean > 0.0,
                    "path-mean hardness is not representable",
                )?;
                return Ok(DepthHardnessAdvance {
                    depth_m: candidate,
                    increment_m,
                    path_mean_hardness_pa: mean,
                });
            }
            consumed = consumed + capacity;
            require(
                finite_dd(consumed),
                "accumulated hardness integral is not representable",
            )?;
            cursor = end;
            start_hardness = end_hardness;
            right += 1;
        }
        Err(DepthHardnessError::Invalid(
            "hardness-integral increment exceeds the prepared exposed-depth support",
        ))
    }
}

fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

fn segment_integral(span: f64, a: f64, b: f64) -> Result<Dd, DepthHardnessError> {
    // This convex mean avoids overflow in a+b and retains equal subnormal H.
    let mean = a.min(b) + 0.5 * (a.max(b) - a.min(b));
    let value = Dd::from_f64(mean) * Dd::from_f64(span);
    require(
        finite_dd(value) && value.to_f64() > 0.0,
        "positive affine hardness integral is not representable",
    )?;
    Ok(value)
}

fn inverse_segment(
    span: f64,
    h0: f64,
    h1: f64,
    measure: f64,
    capacity: f64,
) -> Result<f64, DepthHardnessError> {
    let depth = if h0 == h1 {
        measure / h0
    } else {
        let scale = h0.max(h1);
        let a = h0 / scale;
        let b = h1 / scale;
        let fraction = measure / capacity;
        // For an admitted fraction f, the discriminant is the convex sum
        // (1-f)*a^2 + f*b^2. This avoids subtracting nearly equal terms for
        // softening and avoids squaring dimensional hardness values.
        let root = sqrt((1.0 - fraction) * (a * a) + fraction * (b * b));
        (measure / scale) / (0.5 * a + 0.5 * root)
    };
    require(
        depth.is_finite() && depth > 0.0 && span > 0.0,
        "affine hardness inverse is below arithmetic resolution",
    )?;
    Ok(depth)
}
