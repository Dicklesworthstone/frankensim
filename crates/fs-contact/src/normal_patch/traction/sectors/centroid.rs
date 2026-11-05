//! Positive radial centroids of the original Hertz, JKR and Maugis fields.
//!
//! The caller supplies exact annular force, so a first moment is force times
//! this bounded radius, never pressure times an avoidably overflowing a³.
//! This production numerical leaf has no project dependencies; its unit tests
//! can also run directly with rustc when the full runtime is unavailable.

use core::f64::consts::FRAC_PI_2;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum CentroidError {
    Invalid,
    Unresolved,
    Budget { evaluations: usize },
    Cancelled,
}

#[derive(Clone, Copy)]
pub(crate) struct Controls {
    pub absolute_fraction_tolerance: f64,
    pub relative_tolerance: f64,
    pub maximum_evaluations: usize,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Centroid {
    pub radius_m: f64,
    pub error_estimate_m: f64,
    pub evaluations: usize,
}

#[derive(Clone, Copy)]
struct Geometry {
    u: f64,
    v: f64,
    x: f64,
    y: f64,
    width: f64,
}

fn geometry(a: f64, lo: f64, hi: f64) -> Option<Geometry> {
    if !a.is_finite()
        || !lo.is_finite()
        || !hi.is_finite()
        || a <= 0.0
        || lo < 0.0
        || hi <= lo
        || hi > a
    {
        return None;
    }
    let u = lo / a;
    let v = hi / a;
    let x = (((a - lo) / a) * (1.0 + u)).sqrt();
    let y = (((a - hi) / a) * (1.0 + v)).sqrt();
    // sin(delta theta) from a factored difference of squares. Subtracting
    // asin endpoints or v*x-u*y loses one-ULP interior and central rings.
    let cross = ((hi - lo) / a) * ((u + v) / (v * x + u * y));
    let width = cross.atan2(u * v + x * y);
    (width.is_finite() && width > 0.0 && v > 0.0).then_some(Geometry { u, v, x, y, width })
}

/// Mean radius for constant pressure with radial area measure.
pub(crate) fn uniform(lo: f64, hi: f64) -> f64 {
    let r = lo / hi;
    hi * ((2.0 / 3.0) * (1.0 + r + r * r) / (1.0 + r))
}

// sinc(x), and (1-sinc(x))/x². The latter remains well-conditioned at zero.
fn sinc_defect(x: f64) -> (f64, f64) {
    let x2 = x * x;
    if x.abs() < 0.5 {
        let defect = 1.0 / 6.0
            + x2 * (-1.0 / 120.0
                + x2 * (1.0 / 5040.0
                    + x2 * (-1.0 / 362880.0
                        + x2 * (1.0 / 39916800.0
                            + x2 * (-1.0 / 6227020800.0 + x2 / 1307674368000.0)))));
        (1.0 - x2 * defect, defect)
    } else {
        let sinc = x.sin() / x;
        (sinc, (1.0 - sinc) / x2)
    }
}

fn support_radius(value: f64, lo: f64, hi: f64) -> Option<f64> {
    let roundoff = 64.0 * f64::EPSILON * hi;
    if !value.is_finite() || value < lo - roundoff || value > hi + roundoff {
        return None;
    }
    // Round only within the known centroid support, after the arithmetic
    // check. No pressure, tensile edge, force or physical branch is clipped.
    Some(value.clamp(lo, hi))
}

/// Exact continuum centroid formulas for Hertz and tensile JKR contributions.
pub(crate) fn jkr(a: f64, lo: f64, hi: f64) -> Option<[f64; 2]> {
    let g = geometry(a, lo, hi)?;
    let norm = (g.u + g.v).hypot(g.x + g.y);
    let sm = (g.u + g.v) / norm;
    let cm = (g.x + g.y) / norm;
    let (s1, d1) = sinc_defect(g.width);
    let (s2, d2) = sinc_defect(2.0 * g.width);
    let (sh, _) = sinc_defect(0.5 * g.width);
    let ratio = g.width / sm;
    // Integrals of sin²(theta), and sin²(theta)cos²(theta), divided by
    // their respective force integrals. Factored expressions avoid small
    // center powers and cancellation at the singular contact edge.
    let tensile = (sm * s1 + 0.5 * g.width * ratio * d1) / sh;
    let hertz = (3.0 / 8.0) * (8.0 * sm * cm * cm * s2 + 4.0 * g.width * ratio * d2)
        / (sh * (g.x * g.x + g.x * g.y + g.y * g.y));
    Some([
        support_radius(a * hertz, lo, hi)?,
        support_radius(a * tensile, lo, hi)?,
    ])
}

const G4: [(f64, f64); 2] = [
    (0.3399810435848563, 0.6521451548625461),
    (0.8611363115940526, 0.3478548451374538),
];
const G8: [(f64, f64); 4] = [
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
];
const G16: [(f64, f64); 8] = [
    (0.09501250983763744, 0.1894506104550685),
    (0.2816035507792589, 0.1826034150449236),
    (0.4580167776572274, 0.1691565193950025),
    (0.6178762444026438, 0.1495959888165767),
    (0.7554044083550030, 0.1246289712555339),
    (0.8656312023878318, 0.09515851168249278),
    (0.9445750230732326, 0.06225352393864789),
    (0.9894009349916499, 0.02715245941175409),
];

fn quadrature(f: &impl Fn(f64) -> [f64; 2], lo: f64, hi: f64, rule: &[(f64, f64)]) -> [f64; 2] {
    let mid = 0.5 * lo + 0.5 * hi;
    let half = 0.5 * (hi - lo);
    let mut sum = [0.0; 2];
    for &(node, weight) in rule {
        let left = f(mid - half * node);
        let right = f(mid + half * node);
        for k in 0..2 {
            sum[k] += weight * (left[k] + right[k]);
        }
    }
    [sum[0] * half, sum[1] * half]
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let y = value - *correction;
    let next = *sum + y;
    *correction = (next - *sum) - y;
    *sum = next;
}

/// Maugis centroid against the original analytic mean cohesive angle. The
/// 4/8/16-point ladder and source-normalization check report an estimate, not
/// an interval certificate. Geometric edge bands resolve arbitrarily thin
/// cohesive layers subject to the explicit work budget.
#[allow(clippy::too_many_arguments)]
pub(crate) fn maugis(
    a: f64,
    l: f64,
    lo: f64,
    hi: f64,
    mean_angle: f64,
    controls: Controls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Centroid, CentroidError> {
    let g = geometry(a, lo, hi).ok_or(CentroidError::Invalid)?;
    let q = (l / a).min(1.0);
    let absolute = controls.absolute_fraction_tolerance;
    let relative = controls.relative_tolerance;
    if !l.is_finite()
        || l <= 0.0
        || q == 0.0
        || !mean_angle.is_finite()
        || mean_angle <= 0.0
        || mean_angle > FRAC_PI_2
        || !absolute.is_finite()
        || absolute < 0.0
        || !relative.is_finite()
        || relative < 0.0
    {
        return Err(CentroidError::Invalid);
    }
    let normalization = (((hi - lo) / a) / g.width) * (0.5 * (1.0 + g.u / g.v)) * (mean_angle / q);
    let requested = absolute + relative * (uniform(lo, hi) / hi);
    if !normalization.is_finite()
        || normalization <= 0.0
        || !requested.is_finite()
        || requested < 256.0 * f64::EPSILON
    {
        return Err(CentroidError::Invalid);
    }
    let tolerance = 0.25 * normalization * requested;
    let f = |t: f64| {
        // Offset from the outer endpoint, not a difference of nearly equal
        // absolute angles. These sums remain accurate at both disk endpoints.
        let offset = g.width * t;
        let remaining = g.width * (1.0 - t);
        let cosine = g.y * offset.cos() + g.v * offset.sin();
        let sine = g.u * remaining.cos() + g.x * remaining.sin();
        let r = sine / g.v;
        let root = a * cosine;
        let weight = if l <= a {
            if root >= l {
                let z = l / root;
                if z == 0.0 { 1.0 } else { z.atan() / z }
            } else {
                (root / l) * l.atan2(root)
            }
        } else {
            cosine * l.atan2(root)
        };
        [r * weight, r * r * weight]
    };
    let mut sums = [0.0; 4];
    let mut corrections = [0.0; 4];
    let mut evaluations = 0usize;
    let alpha = g.y.atan2(g.v);
    let mut band = if l < a {
        (l / a).max(alpha)
    } else {
        f64::INFINITY
    };
    if band <= alpha {
        band *= 2.0;
    }
    let mut left = 0.0;
    loop {
        let right = ((band - alpha) / g.width).min(1.0);
        if right <= left {
            band *= 2.0;
            continue;
        }
        let mut stack = [(0.0, 0.0, 0u8); 34];
        stack[0] = (left, right, 0);
        let mut count = 1;
        while count > 0 {
            if cancelled() {
                return Err(CentroidError::Cancelled);
            }
            if controls.maximum_evaluations.saturating_sub(evaluations) < 28 {
                return Err(CentroidError::Budget { evaluations });
            }
            count -= 1;
            let (start, end, depth) = stack[count];
            let coarse = quadrature(&f, start, end, &G4);
            let medium = quadrature(&f, start, end, &G8);
            let fine = quadrature(&f, start, end, &G16);
            evaluations += 28;
            let mut error = [0.0; 2];
            let mut settled = true;
            for k in 0..2 {
                let floor = 32.0 * f64::EPSILON * fine[k].abs();
                let difference = (fine[k] - medium[k]).abs();
                error[k] = 4.0 * difference + floor;
                settled &= fine[k].is_finite()
                    && fine[k] >= 0.0
                    && error[k] <= tolerance * (end - start)
                    && (difference <= 0.5 * (medium[k] - coarse[k]).abs() || difference <= floor);
            }
            if settled {
                for k in 0..2 {
                    add(&mut sums[k], &mut corrections[k], fine[k]);
                    add(&mut sums[k + 2], &mut corrections[k + 2], error[k]);
                }
            } else {
                let mid = 0.5 * start + 0.5 * end;
                if depth >= 32 || mid == start || mid == end {
                    return Err(CentroidError::Unresolved);
                }
                stack[count] = (mid, end, depth + 1);
                stack[count + 1] = (start, mid, depth + 1);
                count += 2;
            }
        }
        if right == 1.0 {
            break;
        }
        left = right;
        band *= 2.0;
    }
    let fraction = sums[1] / normalization;
    let error_fraction = (sums[2] + sums[3] + (sums[0] - normalization).abs()) / normalization;
    let radius_m = support_radius(hi * fraction, lo, hi).ok_or(CentroidError::Unresolved)?;
    if !error_fraction.is_finite() || error_fraction > absolute + relative * fraction {
        return Err(CentroidError::Unresolved);
    }
    Ok(Centroid {
        radius_m,
        error_estimate_m: hi * error_fraction,
        evaluations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    fn close(x: f64, y: f64, scale: f64, tolerance: f64) {
        assert!((x - y).abs() <= scale * tolerance, "{x:e} != {y:e}");
    }

    // Independent full-disk pressure integral, with a large-zone series.
    fn full_angle(b: f64) -> f64 {
        if b <= 1.0 {
            b + (1.0 + b * b) * b.atan() - FRAC_PI_2 * b * b
        } else {
            let x = 1.0 / b;
            let mut term = x;
            let mut sum = 0.0;
            for n in 0..40 {
                sum += 2.0 * term / ((2 * n + 1) * (2 * n + 3)) as f64;
                term *= -x * x;
            }
            FRAC_PI_2 - sum
        }
    }

    fn controls() -> Controls {
        Controls {
            absolute_fraction_tolerance: 0.0,
            relative_tolerance: 1e-11,
            maximum_evaluations: 16384,
        }
    }

    // Independent antiderivative of s atan(b/s) for well-conditioned intervals.
    fn annular_angle(b: f64, lo: f64, hi: f64) -> f64 {
        let primitive = |s: f64| 0.5 * (s * s * b.atan2(s) + b * s - b * b * (s / b).atan());
        let x = ((1.0 - lo) * (1.0 + lo)).sqrt();
        let y = ((1.0 - hi) * (1.0 + hi)).sqrt();
        2.0 * (primitive(x) - primitive(y)) / ((hi - lo) * (hi + lo))
    }

    #[test]
    fn maugis_partial_annuli_preserve_first_moment_under_subdivision() {
        for b in [1e-12, 1e-8, 0.02, 0.2, 1.0] {
            let edges = [0.13, 0.52, 0.86, 1.0];
            let mut moment = 0.0;
            for pair in edges.windows(2) {
                let (lo, hi) = (pair[0], pair[1]);
                let angle = annular_angle(b, lo, hi);
                let result = maugis(1.0, b, lo, hi, angle, controls(), &mut || false).unwrap();
                assert!(result.radius_m >= lo && result.radius_m <= hi);
                assert!(result.evaluations <= controls().maximum_evaluations);
                moment += 0.5 * (hi - lo) * (hi + lo) * angle * result.radius_m;
            }
            let angle = annular_angle(b, edges[0], 1.0);
            let total = maugis(1.0, b, edges[0], 1.0, angle, controls(), &mut || false).unwrap();
            close(
                moment,
                0.5 * (1.0 - edges[0] * edges[0]) * angle * total.radius_m,
                moment,
                4e-12,
            );
        }
        // Separate Simpson integration of the original pressure moment. It
        // neither uses the production quadrature rules nor a fitted centroid.
        let (b, lo, hi) = (0.2_f64, 0.13_f64, 0.86_f64);
        let h = (hi.asin() - lo.asin()) / 8192.0;
        let mut reference = 0.0;
        for i in 0..=8192 {
            let t = lo.asin() + i as f64 * h;
            let value = b.atan2(t.cos()) * t.sin().powi(2) * t.cos();
            let weight = if i == 0 || i == 8192 {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            reference += weight * value;
        }
        reference *= h / 3.0;
        let angle = annular_angle(b, lo, hi);
        let result = maugis(1.0, b, lo, hi, angle, controls(), &mut || false).unwrap();
        close(
            0.5 * (hi - lo) * (hi + lo) * angle * result.radius_m,
            reference,
            reference,
            2e-12,
        );
        for (b, lo, hi, angle) in [
            (
                1e-12,
                1.0 - 1e-9,
                1.0,
                annular_angle(1e-12, 1.0 - 1e-9, 1.0),
            ),
            (
                0.2,
                0.4,
                0.4_f64.next_up(),
                0.2_f64.atan2((1.0_f64 - 0.16).sqrt()),
            ),
        ] {
            let result = maugis(1.0, b, lo, hi, angle, controls(), &mut || false).unwrap();
            assert!(result.radius_m >= lo && result.radius_m <= hi);
        }
    }

    #[test]
    fn jkr_centroids_retain_exact_disk_values_and_narrow_edge_support() {
        let [h, t] = jkr(1.0, 0.0, 1.0).unwrap();
        close(h, 3.0 * PI / 16.0, 1.0, 5e-15);
        close(t, PI / 4.0, 1.0, 5e-15);
        for (lo, hi) in [
            (0.0, 1e-140),
            (0.4, 0.4_f64.next_up()),
            (1.0_f64.next_down(), 1.0),
        ] {
            let radii = jkr(1.0, lo, hi).unwrap();
            for r in radii {
                assert!(r >= lo && r <= hi && r.is_finite());
            }
        }
        let [h, t] = jkr(1.0, 0.0, 1e-140).unwrap();
        close(h, (2.0 / 3.0) * 1e-140, 1e-140, 5e-15);
        close(t, h, 1e-140, 5e-15);
    }

    #[test]
    fn maugis_uses_original_force_normalization_and_resolves_thin_zones() {
        for b in [1e-12, 1e-8, 0.02, 0.2, 1.0, 5.0, 1000.0, 1e8] {
            let result =
                maugis(1.0, b, 0.0, 1.0, full_angle(b), controls(), &mut || false).unwrap();
            assert!(result.radius_m >= 2.0 / 3.0 && result.radius_m <= PI / 4.0);
            assert!(result.error_estimate_m <= 1e-11 * result.radius_m);
            if b <= 1e-8 {
                close(result.radius_m, PI / 4.0, 1.0, 2e-8);
            }
            if b >= 1e8 {
                close(result.radius_m, 2.0 / 3.0, 1.0, 2e-8);
            }
            if b == 1.0 {
                let exact = PI / 3.0 * (3.5 - 2.0_f64.powf(1.5));
                close(result.radius_m, exact, 1.0, 2e-12);
            }
        }
    }

    #[test]
    fn moment_work_and_cancellation_budgets_are_actual_refusals() {
        let mut tiny = controls();
        tiny.maximum_evaluations = 27;
        assert_eq!(
            maugis(1.0, 1.0, 0.0, 1.0, 1.0, tiny, &mut || false).unwrap_err(),
            CentroidError::Budget { evaluations: 0 }
        );
        assert_eq!(
            maugis(1.0, 1.0, 0.0, 1.0, 1.0, controls(), &mut || true).unwrap_err(),
            CentroidError::Cancelled
        );
        assert!(maugis(1.0, 1.0, 0.0, 1.0, 0.9, controls(), &mut || false).is_err());
        assert!(jkr(1.0, 0.5, 0.5).is_none());
    }
}
