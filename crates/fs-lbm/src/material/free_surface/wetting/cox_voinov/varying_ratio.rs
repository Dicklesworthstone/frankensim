//! Invert the original Cox integral only over the actual angle change.
//!
//! A changing thermal/composition source can give every wall cell a different
//! viscosity ratio. Rebuilding a full angular table at each such sample spends
//! almost all its work on unused angles. Here a safeguarded scalar Newton solve
//! integrates g' between the microscopic angle and its candidate, with bounded
//! adaptive Gauss quadrature and a fixed stack. Scaling by max(1,M) preserves
//! the same relation at large viscosity contrast. Quadrature differences are
//! numerical error estimates, not rigorous enclosures or model-error bounds.

use super::{CancelGate, CoxVoinovControls, Error, checkpoint, scaled_angular_slope};

const ANGLE_TOLERANCE: f64 = 2e-13;
const MAX_ITERATIONS: usize = 64;
const MAX_SLOPE_EVALUATIONS: usize = 32_768;
const MAX_DEPTH: usize = 24;

pub(super) struct Answer {
    pub angle: f64,
    #[cfg(test)]
    pub evaluations: usize,
}

pub(super) fn apparent_angle(
    controls: &CoxVoinovControls,
    logarithmic_scale_ratio: f64,
    microscopic: f64,
    capillary: f64,
    ratio: f64,
    gate: &CancelGate,
) -> Result<Answer, Error> {
    checkpoint(gate)?;
    let direction = capillary.signum();
    let boundary = if capillary > 0.0 {
        controls.maximum_angle_radians
    } else {
        controls.minimum_angle_radians
    };
    let extent = (boundary - microscopic).abs();
    let target = (capillary * ratio.max(1.0)).abs() * logarithmic_scale_ratio;
    if !target.is_finite() || extent == 0.0 {
        return Err(Error::Invalid(
            "Cox moving contact line has no angle in its admitted domain",
        ));
    }
    let initial_slope = scaled_angular_slope(microscopic, ratio);
    let mut distance = (target / initial_slope).min(extent);
    let (mut lower, mut upper) = (0.0, extent);
    let mut evaluations = 1;
    for _ in 0..MAX_ITERATIONS {
        checkpoint(gate)?;
        // An exact endpoint stays inside the original physical domain even
        // when adding a large increment would round beyond that endpoint.
        let angle = if distance == extent {
            boundary
        } else {
            microscopic + direction * distance
        };
        if evaluations == MAX_SLOPE_EVALUATIONS {
            return Err(Error::Invalid(
                "Cox local angular quadrature exhausted its budget",
            ));
        }
        let slope = scaled_angular_slope(angle, ratio);
        evaluations += 1;
        if !slope.is_finite() || slope <= 0.0 {
            return Err(Error::Invalid(
                "Cox local angular derivative is unrepresentable",
            ));
        }
        let tolerance = ANGLE_TOLERANCE * angle.abs().max(1.0);
        // A remote trial (especially an angular endpoint) can have a much
        // smaller slope than the eventual root. First resolve its residual;
        // request final angular accuracy only when uncertainty reaches that
        // residual. Requiring final precision at every trial can exhaust the
        // quadrature at a perfectly admissible, well-conditioned root.
        let angular_tolerance = tolerance * slope;
        let mut integral = integrate(
            microscopic.min(angle),
            microscopic.max(angle),
            ratio,
            (0.125 * angular_tolerance).max(1e-8 * target),
            &mut evaluations,
            gate,
        )?;
        let mut residual = integral.value - target;
        if integral.error > 0.125 * angular_tolerance && residual.abs() <= 8.0 * integral.error {
            integral = integrate(
                microscopic.min(angle),
                microscopic.max(angle),
                ratio,
                0.125 * angular_tolerance,
                &mut evaluations,
                gate,
            )?;
            residual = integral.value - target;
        }
        if residual.abs() + integral.error <= angular_tolerance {
            return Ok(Answer {
                angle,
                #[cfg(test)]
                evaluations,
            });
        }
        if distance == extent && residual + integral.error < 0.0 {
            return Err(Error::Invalid(
                "Cox moving contact line has no angle in its admitted domain",
            ));
        }
        if residual > 0.0 {
            upper = distance;
        } else {
            lower = distance;
        }
        let newton = distance - residual / slope;
        let next = if newton > lower && newton < upper {
            newton
        } else {
            f64::midpoint(lower, upper)
        };
        if next == distance {
            break;
        }
        distance = next;
    }
    Err(Error::Invalid(
        "Cox local angular inverse did not resolve its accuracy",
    ))
}

#[derive(Clone, Copy, Default)]
struct Integral {
    value: f64,
    error: f64,
}

// Positive integrands need no cancellation-prone primitive subtraction. Split
// only intervals whose independent 4/8-point rules disagree, and keep both
// quadrature discrepancy and accumulation roundoff in the returned estimate.
fn integrate(
    from: f64,
    to: f64,
    ratio: f64,
    tolerance: f64,
    evaluations: &mut usize,
    gate: &CancelGate,
) -> Result<Integral, Error> {
    if from == to {
        return Ok(Integral::default());
    }
    let mut stack = [(0.0, 0.0, 0.0, 0_usize); MAX_DEPTH + 1];
    stack[0] = (from, to, tolerance, 0);
    let mut count = 1;
    let mut result = Integral::default();
    let mut correction = 0.0;
    while count != 0 {
        checkpoint(gate)?;
        count -= 1;
        let (left, right, local_tolerance, depth) = stack[count];
        if *evaluations + 12 > MAX_SLOPE_EVALUATIONS {
            return Err(Error::Invalid(
                "Cox local angular quadrature exhausted its budget",
            ));
        }
        let midpoint = f64::midpoint(left, right);
        let half = 0.5 * (right - left);
        let quadrature = |nodes: &[(f64, f64)]| {
            let mut sum = 0.0;
            for &(x, w) in nodes {
                sum += w
                    * (scaled_angular_slope(midpoint - half * x, ratio)
                        + scaled_angular_slope(midpoint + half * x, ratio));
            }
            half * sum
        };
        let coarse = quadrature(&GAUSS4);
        let fine = quadrature(&GAUSS8);
        *evaluations += 12;
        let discrepancy = (fine - coarse).abs();
        let roundoff = 16.0 * f64::EPSILON * fine;
        let error = discrepancy + roundoff;
        if !fine.is_finite() || fine <= 0.0 || !error.is_finite() {
            return Err(Error::Invalid(
                "Cox local angular integral is unrepresentable",
            ));
        }
        // Subdivision cannot reduce a relative floating-point noise floor.
        // Carry it into the final angular gate instead of recursively asking
        // each smaller interval to satisfy an impossible local tolerance.
        if discrepancy <= local_tolerance || discrepancy <= roundoff {
            let adjusted = fine - correction;
            let next = result.value + adjusted;
            correction = (next - result.value) - adjusted;
            result.value = next;
            result.error += error;
        } else {
            if depth == MAX_DEPTH || midpoint == left || midpoint == right {
                return Err(Error::Invalid(
                    "Cox local angular quadrature did not resolve its accuracy",
                ));
            }
            stack[count] = (midpoint, right, 0.5 * local_tolerance, depth + 1);
            stack[count + 1] = (left, midpoint, 0.5 * local_tolerance, depth + 1);
            count += 2;
        }
    }
    result.error += 4.0 * f64::EPSILON * result.value;
    Ok(result)
}

const GAUSS4: [(f64, f64); 2] = [
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
];
const GAUSS8: [(f64, f64); 4] = [
    (0.183_434_642_495_649_8, 0.362_683_783_378_362),
    (0.525_532_409_916_329, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];
