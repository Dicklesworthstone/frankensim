//! Derivatives of the original scalar contact work, without numerical probes.
use super::{
    DmtError, FiniteGapResponseCurve, FiniteGapResponseCurveError, NormalLawSource,
    NormalLawSourceError, NormalPatchError, SourceDmtSphereContact, unilateral_path,
};

impl NormalLawSource {
    /// Derivative of `mean_force_between(from, to)` with respect to `to`, N/m.
    /// The initial approach is held fixed. Unilateral paths retain their full
    /// signed length, including the inactive fraction before contact.
    ///
    /// At coincident piecewise-law corners the original positive-approach
    /// tangent supplies a one-sided Newton representative (the final table
    /// node uses its admitted left segment). This is not a smooth tangent at
    /// the corner. Unequal paths integrate both sides of each crossed corner.
    /// Coincident JKR and Maugis folds remain singular. Maugis integrates its
    /// original analytic branch derivatives within the source's work budget;
    /// no finite-difference probe or fitted law is substituted.
    pub fn mean_force_end_derivative(
        &self,
        from_m: f64,
        to_m: f64,
    ) -> Result<Option<f64>, NormalLawSourceError> {
        match self {
            Self::JkrAttached(source) => source
                .attached_path_end_derivative(from_m, to_m)
                .map(Some)
                .map_err(NormalLawSourceError::Jkr),
            Self::MaugisAttached(source) => source
                .path_end_derivative(from_m, to_m)
                .map(Some)
                .map_err(NormalLawSourceError::Maugis),
            Self::HertzPoint(source) => {
                let (a, b, fraction) = unilateral_path(from_m, to_m).ok_or(
                    NormalLawSourceError::Point(NormalPatchError::InvalidInput {
                        field: "finite signed normal derivative path",
                    }),
                )?;
                let path = source
                    .conservative_path_average(a, b)
                    .map_err(NormalLawSourceError::Point)?;
                let derivative = if from_m >= 0.0 && to_m >= 0.0 {
                    path.endpoint_derivative_n_per_m
                } else if a == 0.0 && b == 0.0 {
                    0.0
                } else if to_m > 0.0 {
                    // The active path starts at zero: its endpoint derivative
                    // is 3 K sqrt(b)/5. Differentiate its active time fraction
                    // analytically, without dividing a tiny force by b.
                    path.endpoint_derivative_n_per_m
                        * fraction
                        * (1.0 + (2.0 / 3.0) * (1.0 - fraction))
                } else {
                    // The active path ends at zero: its endpoint derivative
                    // is 2 K sqrt(a)/5, though the clipped endpoint is fixed.
                    path.endpoint_derivative_n_per_m * fraction * fraction
                };
                if !derivative.is_finite() {
                    return Err(NormalLawSourceError::Point(NormalPatchError::Overflow {
                        field: "mechanical path derivative",
                    }));
                }
                Ok(Some(derivative))
            }
            Self::DmtSphere(source) => dmt_derivative(source, from_m, to_m)
                .map(Some)
                .map_err(NormalLawSourceError::Dmt),
            Self::FiniteGap(source) => curve_derivative(source, from_m, to_m)
                .map(Some)
                .map_err(NormalLawSourceError::FiniteGap),
        }
    }
}

fn dmt_derivative(
    source: &SourceDmtSphereContact,
    from_m: f64,
    to_m: f64,
) -> Result<f64, DmtError> {
    let from = source.evaluate(from_m)?;
    let to = source.evaluate(to_m)?;
    if from_m == to_m {
        return Ok(0.5 * from.one_sided_tangents_n_per_m[1]);
    }
    let lo = from_m.min(to_m);
    let hi = from_m.max(to_m);
    let span = hi - lo;
    if !span.is_finite() {
        return Err(DmtError::Unrepresentable {
            field: "DMT derivative path span",
        });
    }
    let peak_tangent = if hi <= 0.0 {
        0.0
    } else if from_m > to_m {
        from.one_sided_tangents_n_per_m[1]
    } else {
        to.one_sided_tangents_n_per_m[1]
    };
    let hertz = if lo >= 0.0 {
        let a = (from_m / hi).sqrt();
        let b = (to_m / hi).sqrt();
        let a2 = a * a;
        let b2 = b * b;
        (peak_tangent * (2.0 / 15.0)) * (3.0 * b2 * b + 6.0 * a * b2 + 4.0 * a2 * b + 2.0 * a2 * a)
            / ((a + b) * (a + b))
    } else if hi > 0.0 {
        let fraction = hi / span;
        if to_m > 0.0 {
            (peak_tangent * (2.0 / 3.0)) * fraction * (1.0 - 0.4 * fraction)
        } else {
            (peak_tangent * (4.0 / 15.0)) * fraction * fraction
        }
    } else {
        0.0
    };
    // The original attraction has constant slope only on (-z0, 0).
    // Integrate t F'(from+t*(to-from)) over that phase exactly, retaining
    // negative stiffness and avoiding subtraction of nearly equal forces.
    let left = lo.max(-source.interaction_range_m());
    let right = hi.min(0.0);
    let tail = if right > left {
        let moment = linear_path_moment(from_m, to_m, left, right);
        source.evaluate(0.0)?.one_sided_tangents_n_per_m[0] * moment
    } else {
        0.0
    };
    let derivative = hertz + tail;
    if !derivative.is_finite() {
        return Err(DmtError::Unrepresentable {
            field: "DMT mean-force derivative",
        });
    }
    Ok(derivative)
}

fn curve_derivative(
    source: &FiniteGapResponseCurve,
    from_m: f64,
    to_m: f64,
) -> Result<f64, FiniteGapResponseCurveError> {
    let (from, to, fraction) =
        unilateral_path(from_m, to_m).ok_or(FiniteGapResponseCurveError::InvalidInput {
            field: "finite signed normal derivative path",
        })?;
    // Reuse the source's bounded segment and stored-energy continuity checks.
    // These public nodes cannot supply a derivative of inconsistent work.
    let active_mean = source.mean_force_between(from, to)?;
    if from_m < 0.0 && to_m <= 0.0 || to_m < 0.0 && from_m <= 0.0 {
        return Ok(0.0);
    }
    let active_derivative = if from == to {
        0.5 * source.evaluate(from)?.normal_tangent_n_per_m
    } else {
        let lo = from.min(to);
        let hi = from.max(to);
        let mut index = source.nodes.partition_point(|node| node.approach_m <= lo);
        index = index.clamp(1, source.nodes.len() - 1);
        let mut position = lo;
        let mut derivative = 0.0;
        let mut correction = 0.0;
        loop {
            let left = &source.nodes[index - 1];
            let right = &source.nodes[index];
            let slope =
                (right.normal_force_n - left.normal_force_n) / (right.approach_m - left.approach_m);
            let end = hi.min(right.approach_m);
            let part = slope * linear_path_moment(from, to, position, end);
            let (sum, low) = fs_math::eft::two_sum(derivative, part);
            derivative = sum;
            correction += low;
            if end == hi {
                break;
            }
            position = end;
            index += 1;
        }
        derivative + correction
    };
    let derivative = if from_m >= 0.0 && to_m >= 0.0 {
        active_derivative
    } else if to_m > 0.0 {
        fraction * active_derivative + (active_mean / (to_m - from_m)) * (1.0 - fraction)
    } else {
        (active_mean / (from_m - to_m)) * fraction
    };
    if !derivative.is_finite() || derivative < 0.0 {
        return Err(FiniteGapResponseCurveError::InvalidInput {
            field: "response curve path derivative",
        });
    }
    Ok(derivative)
}

// Integral of t dt over one physical subinterval of the signed affine path.
// Factor t1^2-t0^2 so tiny segments near t=1 retain their represented width.
fn linear_path_moment(from: f64, to: f64, left: f64, right: f64) -> f64 {
    let span = (to - from).abs();
    let (first, last) = if to > from {
        ((left - from) / span, (right - from) / span)
    } else {
        ((from - right) / span, (from - left) / span)
    };
    ((right - left) / span) * (0.5 * first + 0.5 * last)
}
