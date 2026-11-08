//! Exact interaction of two uniform spherical mass distributions.
//!
//! Integrate the larger sphere's potential over the smaller sphere. In partial
//! overlap, start with the quadratic interior potential averaged over the whole
//! smaller ball, then integrate its exterior-cap correction. With a >= b,
//! q=b/a, x=d/a and h=x+q-1, the positive correction to -U*a/(G*m1*m2) is
//! h^4 [30q - 6(1-q)h - h^2] / (160 q^3 x).
//! The code factors this in bounded cap coordinates; no runtime quadrature,
//! radius cubes, or subtraction of singular point fields is required.

use super::super::{binary_parts, force_component, pair_magnitudes, scaled_nonnegative};
use crate::Vec3;
use crate::gravity::{GravityError, finite, inverse_radius_potential_increment};
use fs_math::dd::Dd;

/// Force on the first center toward the second, and one pair potential [J].
/// The adapter admits positive finite G, native masses and source radii.
pub(super) fn load(
    g: f64,
    masses: [f64; 2],
    radii: [f64; 2],
    displacement: Vec3,
) -> Result<(Vec3, f64), GravityError> {
    let distance = displacement.stable_norm("uniform sphere pair separation")?;
    let a = radii[0].max(radii[1]);
    let b = radii[0].min(radii[1]);
    if exterior(distance, a, b) {
        // Exact original arithmetic and mass-order symmetry outside both balls.
        let (potential, force) = pair_magnitudes(g, masses[0], masses[1], distance)?;
        return Ok((
            Vec3::new(
                force_component(force, displacement.x, distance)?,
                force_component(force, displacement.y, distance)?,
                force_component(force, displacement.z, distance)?,
            ),
            -potential,
        ));
    }
    let (potential_shape, force_shape) = shape(distance, a, b);
    let potential = -pair_product(g, masses, a, 1, potential_shape, Dd::ONE)?.to_f64();
    let component = |value| -> Result<f64, GravityError> {
        Ok(pair_product(g, masses, a, 3, force_shape, Dd::from_f64(value))?.to_f64())
    };
    Ok((
        Vec3::new(
            component(displacement.x)?,
            component(displacement.y)?,
            component(displacement.z)?,
        ),
        finite(potential)?,
    ))
}

/// U_after - U_before from both original center trajectories. In particular,
/// a transverse drift can do representable work while both rounded radii and
/// both absolute binding energies remain identical.
pub(super) fn potential_change(
    g: f64,
    masses: [f64; 2],
    radii: [f64; 2],
    before: [Vec3; 2],
    after: [Vec3; 2],
) -> Result<Dd, GravityError> {
    let rb = before[1]
        .sub(before[0])
        .stable_norm("uniform sphere pair separation")?;
    let ra = after[1]
        .sub(after[0])
        .stable_norm("uniform sphere pair separation")?;
    let a = radii[0].max(radii[1]);
    let b = radii[0].min(radii[1]);
    if exterior(rb, a, b) && exterior(ra, a, b) {
        return exterior_change(g, masses, before, after, [rb, ra]);
    }

    let rb_precise = precise_radius(before, rb)?;
    let ra_precise = precise_radius(after, ra)?;
    let before_region = region(rb_precise, a, b);
    let after_region = region(ra_precise, a, b);
    let (start, end, first_region, last_region, sign) = if ra_precise.lt(rb_precise) {
        (ra_precise, rb_precise, after_region, before_region, -1.0)
    } else {
        (rb_precise, ra_precise, before_region, after_region, 1.0)
    };
    if first_region == last_region {
        let mean = match first_region {
            0 => Dd::ONE,
            1 => overlap_mean(a, b, start, end),
            _ => return exterior_change(g, masses, before, after, [rb, ra]),
        };
        return coordinate_change(g, masses, a, mean, before, after);
    }

    // Split crossings at the exact physical interfaces. All outward pieces
    // have the same sign, so no large absolute binding energies cancel.
    let mut change = Dd::ZERO;
    let mut first = start;
    for index in first_region..=last_region {
        let last = if index == last_region {
            end
        } else if index == 0 {
            Dd::from_f64(a) - Dd::from_f64(b)
        } else {
            // An admitted exterior endpoint guarantees this sum is finite.
            Dd::from_f64(a) + Dd::from_f64(b)
        };
        let piece = if index == 2 {
            exterior_interval(g, masses, first, last)?
        } else {
            let mean = if index == 0 {
                Dd::ONE
            } else {
                overlap_mean(a, b, first, last)
            };
            let travel = last - first;
            let weight = mean * Dd::from_f64(0.5);
            weighted_pair_product(g, masses, a, 3, travel, first, weight)?
                + weighted_pair_product(g, masses, a, 3, travel, last, weight)?
        };
        change = change + piece;
        first = last;
    }
    finite(change.to_f64())?;
    Ok(change * Dd::from_f64(sign))
}

fn exterior_change(
    g: f64,
    masses: [f64; 2],
    before: [Vec3; 2],
    after: [Vec3; 2],
    radii: [f64; 2],
) -> Result<Dd, GravityError> {
    // Keep the original mutual monopole arithmetic outside both supports.
    let (magnitude_before, _) = pair_magnitudes(g, masses[0], masses[1], radii[0])?;
    if radii[0] < 0.5 * radii[1] || radii[1] < 0.5 * radii[0] {
        let (magnitude_after, _) = pair_magnitudes(g, masses[0], masses[1], radii[1])?;
        return Ok(Dd::from_f64(magnitude_before) - Dd::from_f64(magnitude_after));
    }
    Ok(inverse_radius_potential_increment(
        before,
        after,
        radii,
        magnitude_before,
    ))
}

fn precise_radius(centers: [Vec3; 2], radius: f64) -> Result<Dd, GravityError> {
    if radius == 0.0 {
        return Ok(Dd::ZERO);
    }
    let mut square = Dd::ZERO;
    for (first, second) in [
        (centers[0].x, centers[1].x),
        (centers[0].y, centers[1].y),
        (centers[0].z, centers[1].z),
    ] {
        let q = (Dd::from_f64(second) - Dd::from_f64(first)) / Dd::from_f64(radius);
        square = square + q * q;
    }
    let (radius, exponent) = binary_parts(radius);
    scale_dd(Dd::from_f64(radius) * square.sqrt(), exponent)
}

// 0: containment, 1: partial overlap, 2: exterior. No a+b overflow.
fn region(distance: Dd, a: f64, b: f64) -> usize {
    let a = Dd::from_f64(a);
    let b = Dd::from_f64(b);
    if !(a - b).lt(distance) {
        0
    } else if distance.lt(a) || (distance - a).lt(b) {
        1
    } else {
        2
    }
}

fn coordinate_change(
    g: f64,
    masses: [f64; 2],
    a: f64,
    mean: Dd,
    before: [Vec3; 2],
    after: [Vec3; 2],
) -> Result<Dd, GravityError> {
    let mut change = Dd::ZERO;
    for (a0, a1, b0, b1) in [
        (after[0].x, after[1].x, before[0].x, before[1].x),
        (after[0].y, after[1].y, before[0].y, before[1].y),
        (after[0].z, after[1].z, before[0].z, before[1].z),
    ] {
        let [a0, a1, b0, b1] = [a0, a1, b0, b1].map(Dd::from_f64);
        let relative_after = a1 - a0;
        let relative_before = b1 - b0;
        let mut travel = (a1 - b1) - (a0 - b0);
        let mut weight = mean * Dd::from_f64(0.5);
        if !travel.hi.is_finite() {
            // Common translations can overflow individual trajectories even
            // while both endpoint separations and the work remain admitted.
            travel = relative_after * Dd::from_f64(0.5) - relative_before * Dd::from_f64(0.5);
            weight = mean;
        }
        // Keep the endpoint products separate, and scale each physical factor
        // before multiplying. Neither a radius cube nor travel*radius nor the
        // standalone G*m1*m2/a^3 coefficient needs to fit in f64.
        change = change
            + weighted_pair_product(g, masses, a, 3, travel, relative_after, weight)?
            + weighted_pair_product(g, masses, a, 3, travel, relative_before, weight)?;
    }
    finite(change.to_f64())?;
    Ok(change)
}

// Mean K with respect to d^2 on one outward partial-overlap interval. The
// correction C=q^2*n(u)/(160*x), n(u)=u^4*(30-6*(1-q)*u-q*u^2).
// Factor its divided difference with geometric sums, rather than subtracting
// nearby potentials. r=u0/u1, t=cap/d, u<=2, and u/x<=2 are all bounded.
fn overlap_mean(a: f64, b: f64, first: Dd, last: Dd) -> Dd {
    let a = Dd::from_f64(a);
    let b = Dd::from_f64(b);
    let x1 = last / a;
    let u0 = ((first - a) + b) / b;
    let u1 = ((last - a) + b) / b;
    if u1 == Dd::ZERO || x1 == Dd::ZERO {
        return Dd::ONE;
    }
    let q = b / a;
    let c = Dd::ONE - q;
    let t0 = if first == Dd::ZERO {
        Dd::ZERO
    } else {
        ((first - a) + b) / first
    };
    let t1 = ((last - a) + b) / last;
    let r = u0 / u1;
    let r2 = r * r;
    let r3 = r2 * r;
    let s4 = Dd::ONE + r + r2 + r3;
    let s5 = s4 + r3 * r;
    let s6 = s5 + r3 * r2;
    let bracket0 = Dd::from_f64(30.0) - Dd::from_f64(6.0) * c * u0 - q * u0 * u0;
    let divided = Dd::from_f64(30.0) * s4
        - Dd::from_f64(6.0) * c * u1 * s5
        - q * u1 * u1 * s6
        - t0 * r3 * bracket0;
    Dd::ONE - t1 * u1 * (u1 / x1) * divided / (Dd::from_f64(80.0) * (Dd::ONE + first / last))
}

fn exterior_interval(g: f64, masses: [f64; 2], first: Dd, last: Dd) -> Result<Dd, GravityError> {
    if first == last {
        return Ok(Dd::ZERO);
    }
    if first.lt(last * Dd::from_f64(0.5)) {
        // Well-separated radii have a large relative potential difference.
        let first_binding = pair_product(
            g,
            masses,
            first.hi,
            1,
            Dd::from_f64(first.hi) / first,
            Dd::ONE,
        )?;
        let last_binding =
            pair_product(g, masses, last.hi, 1, Dd::from_f64(last.hi) / last, Dd::ONE)?;
        return Ok(first_binding - last_binding);
    }
    // Rationalize 1/first - 1/last, retaining a tiny surface crossing.
    pair_product(
        g,
        masses,
        first.hi,
        2,
        last - first,
        (Dd::from_f64(first.hi) / first) * (Dd::from_f64(first.hi) / last),
    )
}

fn exterior(distance: f64, a: f64, b: f64) -> bool {
    // Never form a+b: the sum can overflow, or round away a small radius.
    distance >= a && !(Dd::from_f64(distance) - Dd::from_f64(a)).lt(Dd::from_f64(b))
}

// Return positive P=-U*a/(G*m1*m2), K=F/(G*m1*m2*d/a^3).
// In overlap K stays between 1/8 and 1, including the equal-radius d->0 limit.
fn shape(distance: f64, a: f64, b: f64) -> (Dd, Dd) {
    let a = Dd::from_f64(a);
    let b = Dd::from_f64(b);
    let d = Dd::from_f64(distance);
    let q = b / a;
    let x = d / a;
    let interior = Dd::from_f64(1.5)
        - x * x * Dd::from_f64(0.5)
        - q * q * (Dd::from_f64(3.0) / Dd::from_f64(10.0));
    if !(a - b).lt(d) || x == Dd::ZERO {
        return (interior, Dd::ONE);
    }

    // Physical-unit subtraction retains a very small cap when a-b or a+b
    // would round. u is in [0,2], t in [0,1], and u/x in [0,2]. For equal
    // radii these reduce directly to u=x and t=1, without a 1/d singularity.
    let cap = (d - a) + b;
    let u = cap / b;
    let t = cap / d;
    let c = Dd::ONE - q;
    let u2 = u * u;
    let potential_bracket = Dd::from_f64(30.0) - Dd::from_f64(6.0) * c * u - q * u2;
    let correction = q * t * u * u2 * potential_bracket / Dd::from_f64(160.0);
    let force_bracket = t
        * (Dd::from_f64(90.0) - Dd::from_f64(24.0) * c * u - Dd::from_f64(5.0) * q * u2)
        + (Dd::ONE - t)
            * (Dd::from_f64(120.0) - Dd::from_f64(30.0) * c * u - Dd::from_f64(6.0) * q * u2);
    let coefficient = Dd::ONE - t * u * (u / x) * force_bracket / Dd::from_f64(160.0);
    (interior + correction, coefficient)
}

// G*m1*m2*f1*f2/a^power with sorted mass factors and bounded mantissas. The
// scale itself need not be representable when the final potential/force is.
fn pair_product(
    g: f64,
    masses: [f64; 2],
    radius: f64,
    power: i32,
    first: Dd,
    second: Dd,
) -> Result<Dd, GravityError> {
    weighted_pair_product(g, masses, radius, power, first, second, Dd::ONE)
}

fn weighted_pair_product(
    g: f64,
    masses: [f64; 2],
    radius: f64,
    power: i32,
    first: Dd,
    second: Dd,
    weight: Dd,
) -> Result<Dd, GravityError> {
    if first == Dd::ZERO || second == Dd::ZERO || weight == Dd::ZERO {
        return Ok(Dd::ZERO);
    }
    let (g, eg) = binary_parts(g);
    let (small, es) = binary_parts(masses[0].min(masses[1]));
    let (large, el) = binary_parts(masses[0].max(masses[1]));
    let (radius, er) = binary_parts(radius);
    let (first, e1) = normalized(first)?;
    let (second, e2) = normalized(second)?;
    let (weight, ew) = normalized(weight)?;
    let mut coefficient =
        Dd::from_f64(g) * Dd::from_f64(small) * Dd::from_f64(large) * first * second * weight;
    for _ in 0..power {
        coefficient = coefficient / Dd::from_f64(radius);
    }
    let exponent = eg + es + el + e1 + e2 + ew - power * er;
    scale_dd(coefficient, exponent)
}

fn scale_dd(coefficient: Dd, exponent: i32) -> Result<Dd, GravityError> {
    let value = Dd::from_pair(
        signed_scale(coefficient.hi, exponent)?,
        signed_scale(coefficient.lo, exponent)?,
    );
    finite(value.to_f64())?;
    Ok(value)
}

fn normalized(value: Dd) -> Result<(Dd, i32), GravityError> {
    finite(value.hi)?;
    finite(value.lo)?;
    let (mantissa, exponent) = binary_parts(value.hi.abs());
    Ok((
        Dd::from_pair(
            mantissa.copysign(value.hi),
            signed_scale(value.lo, -exponent)?,
        ),
        exponent,
    ))
}

fn signed_scale(value: f64, exponent: i32) -> Result<f64, GravityError> {
    Ok(scaled_nonnegative(value.abs(), exponent)?.copysign(value))
}
