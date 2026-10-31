//! Inverse-angle enclosures using only outward basic arithmetic and IEEE sqrt.
//!
//! For 0 <= r <= 1, atan(r) = 2 atan(r/(1 + sqrt(1+r*r))). The reduced
//! argument is at most sqrt(2)-1 < 1/2. Twenty-four alternating series terms
//! end in a negative term, so the remainder is in [0, z^49/49]. See NIST
//! DLMF 4.24.3: <https://dlmf.nist.gov/4.24.E3>. No atan/atan2 ULP assumption
//! is used. Quadrant and rectangle extrema are treated geometrically.

use crate::Interval;

/// Version of the real-angle domain, series and endpoint policies.
pub const ANGLE_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;

/// Enclose the real polar angle of `(x, y)` in the closed upper half-plane.
///
/// `y` must have finite nonnegative endpoints; `x` must be finite. The result
/// is in `[0, pi]`, in radians, and covers every NONZERO vector in the box.
/// A box consisting only of the origin, a negative `y`, or nonfinite input
/// returns `None`. A box containing the origin and other vectors retains all
/// their possible directions; it is not assigned an arbitrary axis angle.
/// Signed zero has its real-number meaning, not an atan2 branch-cut tag.
///
/// Angle extrema over a rectangle in this half-plane occur at its nonzero
/// corners (or their limits). There are at most four bounded series evaluations,
/// no allocation or iteration to a tolerance. For a periodic application,
/// split upper/lower half-planes before wrapping: hulling across a seam first
/// would unnecessarily include the opposite hemisphere.
#[must_use]
pub fn atan2_upper_half_plane(y: Interval, x: Interval) -> Option<Interval> {
    if y.lo() < 0.0
        || [y.lo(), y.hi(), x.lo(), x.hi()]
            .into_iter()
            .any(|v| !v.is_finite())
    {
        return None;
    }
    let mut result: Option<Interval> = None;
    for yy in [y.lo(), y.hi()] {
        for xx in [x.lo(), x.hi()] {
            if let Some(angle) = endpoint(yy, xx) {
                result = Some(result.map_or(angle, |old| old.hull(angle)));
            }
        }
    }
    result
}

fn pi() -> Interval {
    Interval::new(
        fs_math::next_down(core::f64::consts::PI),
        fs_math::next_up(core::f64::consts::PI),
    )
}

fn endpoint(y: f64, x: f64) -> Option<Interval> {
    if y == 0.0 {
        return if x > 0.0 {
            Some(Interval::point(0.0))
        } else if x < 0.0 {
            Some(pi())
        } else {
            None
        };
    }
    let half_pi = pi() * Interval::point(0.5);
    if x == 0.0 {
        return Some(half_pi);
    }
    let a = x.abs();
    // Divide the smaller by the larger to avoid overflow, even for opposing
    // extreme exponents. Outward underflow remains a nonzero upper bound.
    let first_quadrant = if y <= a {
        atan_unit(unit(Interval::point(y) / Interval::point(a)))
    } else {
        half_pi - atan_unit(unit(Interval::point(a) / Interval::point(y)))
    };
    let angle = if x < 0.0 { pi() - first_quadrant } else { first_quadrant };
    Some(Interval::new(angle.lo().max(0.0), angle.hi().min(pi().hi())))
}

fn atan_unit(r: Interval) -> Interval {
    if r.hi() == 0.0 {
        return Interval::point(0.0);
    }
    let one = Interval::point(1.0);
    let reduced = r / (one + (one + r * r).sqrt());
    // The exact half-angle identity gives z <= sqrt(2)-1 < 1/2.
    let z = Interval::new(reduced.lo().max(0.0), reduced.hi().min(0.5));
    let square = nonnegative(z * z);
    let mut power = z;
    let mut sum = z;
    for n in 1..24 {
        power = nonnegative(power * square);
        let term = power / Interval::point(f64::from(2 * n + 1));
        sum = if n % 2 == 0 { sum + term } else { sum - term };
    }
    let next = nonnegative(power * square) / Interval::point(49.0);
    let enclosed = sum + Interval::new(0.0, next.hi());
    nonnegative(Interval::point(2.0) * enclosed)
}

fn nonnegative(v: Interval) -> Interval {
    Interval::new(v.lo().max(0.0), v.hi().max(0.0))
}

fn unit(v: Interval) -> Interval {
    Interval::new(v.lo().max(0.0), v.hi().min(1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn independently_rounded_angles_and_quadrants_are_enclosed() {
        // Rounded independently from high-precision atan ratios, not fs-math.
        for (y, x, reference) in [
            (1.0, 1.0, 0.785_398_163_397_448_3_f64),
            (1.0, 2.0, 0.463_647_609_000_806_1),
            (2.0, 1.0, 1.107_148_717_794_090_4),
            (1.0, -1.0, 2.356_194_490_192_345),
            (1.0, -2.0, 2.677_945_044_588_987),
            (2.0, -1.0, 2.034_443_935_795_702_7),
        ] {
            let value = atan2_upper_half_plane(Interval::point(y), Interval::point(x)).unwrap();
            assert!(value.lo() <= fs_math::next_down(reference));
            assert!(value.hi() >= fs_math::next_up(reference));
            assert!(value.width() < 5e-14, "{value:?}");
        }
    }

    #[test]
    fn axes_and_origin_do_not_invent_a_direction() {
        let zero = Interval::point(0.0);
        assert_eq!(atan2_upper_half_plane(zero, zero), None);
        assert_eq!(atan2_upper_half_plane(zero, Interval::new(1.0, 3.0)), Some(zero));
        let upper = atan2_upper_half_plane(Interval::new(0.0, 1.0), Interval::new(-1.0, 1.0)).unwrap();
        assert_eq!(upper.lo(), 0.0);
        assert!(upper.encloses(pi()));
        let vertical = atan2_upper_half_plane(Interval::new(0.0, 1.0), zero).unwrap();
        assert!(vertical.contains(core::f64::consts::FRAC_PI_2));
        assert!(vertical.width() < 3e-15);
    }

    #[test]
    fn tiny_ratios_and_extreme_scales_keep_positive_angle_support() {
        for (y, x) in [(f64::from_bits(1), f64::MAX), (f64::from_bits(1), 1.0)] {
            let angle = atan2_upper_half_plane(Interval::point(y), Interval::point(x)).unwrap();
            assert_eq!(angle.lo(), 0.0);
            assert!(angle.hi() > 0.0 && angle.hi() < f64::MIN_POSITIVE);
        }
        for magnitude in [f64::from_bits(1), f64::MIN_POSITIVE, 1.0, f64::MAX] {
            let angle = atan2_upper_half_plane(Interval::point(magnitude), Interval::point(magnitude)).unwrap();
            assert!(angle.contains(core::f64::consts::FRAC_PI_4));
            assert!(angle.width() < 5e-14);
        }
    }

    #[test]
    fn rectangle_extrema_are_not_only_the_midpoint_angle() {
        let right = atan2_upper_half_plane(Interval::new(1.0, 2.0), Interval::new(1.0, 2.0)).unwrap();
        assert!(right.contains(0.463_647_609_000_806_1));
        assert!(right.contains(1.107_148_717_794_090_4));
        let left = atan2_upper_half_plane(Interval::new(1.0, 2.0), Interval::new(-2.0, -1.0)).unwrap();
        assert!(left.contains(2.034_443_935_795_702_7));
        assert!(left.contains(2.677_945_044_588_987));
    }

    #[test]
    fn invalid_domain_refuses() {
        for (y, x) in [
            (Interval::new(-f64::from_bits(1), 1.0), Interval::point(1.0)),
            (Interval::point(f64::INFINITY), Interval::point(1.0)),
            (Interval::point(1.0), Interval::WHOLE),
        ] {
            assert!(atan2_upper_half_plane(y, x).is_none());
        }
    }
}
