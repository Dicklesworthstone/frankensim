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
use crate::gravity::{GravityError, finite};
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
    if first == Dd::ZERO || second == Dd::ZERO {
        return Ok(Dd::ZERO);
    }
    let (g, eg) = binary_parts(g);
    let (small, es) = binary_parts(masses[0].min(masses[1]));
    let (large, el) = binary_parts(masses[0].max(masses[1]));
    let (radius, er) = binary_parts(radius);
    let (first, e1) = normalized(first)?;
    let (second, e2) = normalized(second)?;
    let mut coefficient =
        Dd::from_f64(g) * Dd::from_f64(small) * Dd::from_f64(large) * first * second;
    for _ in 0..power {
        coefficient = coefficient / Dd::from_f64(radius);
    }
    let exponent = eg + es + el + e1 + e2 - power * er;
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
