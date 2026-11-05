//! Compensated angular momentum through the original rigid Hamiltonian split.
//!
//! The median-inertia spherical Hamiltonian commutes with the two residual
//! principal-axis Hamiltonians. Their symmetric composition is the same map
//! as the native scalar Lie owner. Here the momentum, subflow phases and
//! quaternion products retain DD tails until the accepted attitude is formed.
//! A sub-ULP half-impulse can therefore affect both spin and finite orientation.
//! This module does not own another public quaternion or inertia representation.

use super::LieStepError;
use fs_ga::{Quat, So3, Vec3};
use fs_math::{dd::Dd, det};

type Vector = [Dd; 3];
type Product = [Dd; 4];

const IDENTITY: Product = [Dd::ONE, Dd::ZERO, Dd::ZERO, Dd::ZERO];

/// Free rotation with retained angular momentum in principal body coordinates.
///
/// `momentum` contains normalized or unnormalized finite hi/lo pairs; admission
/// renormalizes their representation without rounding them to a single scalar.
/// For median inertia `I0`, the existing split is
/// `H = |L|²/(2 I0) + (1/Ia-1/I0) La²/2 + (1/Ib-1/I0) Lb²/2`.
/// Each exact residual flow pairs `R <- R F` with `L <- F^-1 L`. The spherical
/// flow rotates about the retained `L` and leaves that momentum unchanged.
/// The symmetric composition is second order and reversible; spherical and
/// axisymmetric motion remain exact up to arithmetic error.
///
/// The retained phase drives the existing deterministic DD trigonometry.
/// Quaternion products remain private fixed-size arithmetic until the final
/// canonical `So3` admission; the stored attitude is never renormalized.
/// No allocation, eigensolve, nonlinear iteration or momentum projection occurs.
/// The caller may apply its endpoint impulse to the returned DD momentum before
/// rounding its accepted native state. This is not an arbitrary-range precision
/// guarantee or an exact-energy claim for asymmetric bodies at finite step size.
///
/// # Errors
/// Refuses nonfinite/zero duration, nonpositive/nonfinite principal moments,
/// nonfinite or unrepresentable momentum/phase arithmetic, and failed `So3`
/// admission. Negative durations are admitted for the inverse map.
pub fn rigid_body_momentum_step(
    rotation: So3,
    momentum: [Dd; 3],
    inertia: Vec3,
    h: f64,
) -> Result<(So3, [Dd; 3]), LieStepError> {
    let momentum = admit(momentum, inertia, h)?;
    let (increment, next) = principal_step(momentum, inertia, h)?;
    Ok((finish_attitude(rotation, increment)?, next))
}

/// The same compensated map for a full SPD inertia in fixed body coordinates.
///
/// `principal_to_body` is the admitted, fixed principal frame and `inertia`
/// its positive moments, as retained by `gyrostat::InertiaTensor3`. Both input
/// and output momentum use the caller's body coordinates. Only the increment
/// is conjugated through the fixed frame, and its inverse transports the
/// original body momentum. No rounded full-state basis roundtrip occurs.
///
/// The hot path has fixed storage and work and performs no factorization.
/// The frame must remain fixed over the step; changing inertia, internal rotors
/// and deforming geometry require their own momentum transport law.
///
/// # Errors
/// The same finite-parameter, derived-arithmetic and group refusals as
/// [`rigid_body_momentum_step`].
pub fn rigid_body_momentum_step_in_frame(
    rotation: So3,
    momentum_body: [Dd; 3],
    inertia: Vec3,
    principal_to_body: So3,
    h: f64,
) -> Result<(So3, [Dd; 3]), LieStepError> {
    if principal_to_body == So3::identity() || (inertia.x == inertia.y && inertia.y == inertia.z) {
        // Spherical inertia has no preferred principal frame. Preserve its
        // analytic unchanged body momentum instead of rotating it roundtrip.
        return rigid_body_momentum_step(rotation, momentum_body, inertia, h);
    }
    let original = admit(momentum_body, inertia, h)?;
    let basis = product(principal_to_body);
    let principal = rotate(conjugate(basis), original)?;
    let (increment, _) = principal_step(principal, inertia, h)?;
    let vector = rotate(basis, [increment[1], increment[2], increment[3]])?;
    let body_increment = [increment[0], vector[0], vector[1], vector[2]];
    let next = rotate(conjugate(body_increment), original)?;
    Ok((finish_attitude(rotation, body_increment)?, next))
}

fn admit(momentum: Vector, inertia: Vec3, h: f64) -> Result<Vector, LieStepError> {
    if !h.is_finite() || h == 0.0 {
        return Err(invalid("compensated rigid-body step size"));
    }
    if ![inertia.x, inertia.y, inertia.z]
        .into_iter()
        .all(|value| value.is_finite() && value > 0.0)
    {
        return Err(invalid("compensated rigid-body principal inertia"));
    }
    let mut normalized = [Dd::ZERO; 3];
    for (target, value) in normalized.iter_mut().zip(momentum) {
        checked(value, "compensated rigid-body angular momentum")?;
        *target = checked(
            Dd::from_pair(value.hi, value.lo),
            "compensated rigid-body angular momentum",
        )?;
    }
    Ok(normalized)
}

fn principal_step(
    mut momentum: Vector,
    inertia: Vec3,
    h: f64,
) -> Result<(Product, Vector), LieStepError> {
    let moments = [inertia.x, inertia.y, inertia.z];
    let mut axes = [0, 1, 2];
    axes.sort_unstable_by(|&a, &b| moments[a].total_cmp(&moments[b]));
    let reference = moments[axes[1]];
    let a = axes[0];
    let b = axes[2];
    let mut increment = IDENTITY;
    let mut axis_step = |axis: usize, fraction: f64| -> Result<(), LieStepError> {
        let component = momentum[axis];
        let rate = divide(component, Dd::from_f64(moments[axis]))?
            - divide(component, Dd::from_f64(reference))?;
        let half_angle = phase(rate, h, 0.5 * fraction)?;
        if half_angle.hi == 0.0 && half_angle.lo == 0.0 {
            return Ok(());
        }
        let (sine, cosine) = det::sin_cos_dd(half_angle);
        let mut flow = [Dd::ZERO; 4];
        flow[0] = checked(cosine, "compensated rigid-body rotation")?;
        flow[axis + 1] = checked(sine, "compensated rigid-body rotation")?;
        momentum = rotate(conjugate(flow), momentum)?;
        increment = multiply(increment, flow);
        Ok(())
    };
    if moments[a] == reference {
        axis_step(b, 1.0)?;
    } else if moments[b] == reference {
        axis_step(a, 1.0)?;
    } else {
        axis_step(a, 0.5)?;
        axis_step(b, 1.0)?;
        axis_step(a, 0.5)?;
    }
    let mut half_tangent = [Dd::ZERO; 3];
    for (target, value) in half_tangent.iter_mut().zip(momentum) {
        *target = phase(divide(value, Dd::from_f64(reference))?, h, 0.5)?;
    }
    // This exact subflow rotates about L itself. Preserve L analytically,
    // avoiding a spurious small transport of a potentially large momentum.
    increment = multiply(increment, exponential(half_tangent)?);
    Ok((increment, momentum))
}

// Evaluate the quaternion exponential from its half rotation vector. Scaling
// avoids squaring large components; normalization stays in DD so a small
// component is not lost from the rotation axis or its angle.
fn exponential(half_tangent: Vector) -> Result<Product, LieStepError> {
    let scale = half_tangent
        .iter()
        .map(|value| value.hi.abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Ok(IDENTITY);
    }
    let mut unit_scale = [Dd::ZERO; 3];
    for (target, value) in unit_scale.iter_mut().zip(half_tangent) {
        *target = divide(value, Dd::from_f64(scale))?;
    }
    let length = checked(
        dot(unit_scale, unit_scale).sqrt(),
        "compensated rotation norm",
    )?;
    let angle = checked(Dd::from_f64(scale) * length, "compensated rotation phase")?;
    let (sine, cosine) = det::sin_cos_dd(angle);
    let mut result = [Dd::ZERO; 4];
    result[0] = checked(cosine, "compensated rigid-body rotation")?;
    for d in 0..3 {
        result[d + 1] = checked(
            divide(unit_scale[d], length)? * sine,
            "compensated rigid-body rotation",
        )?;
    }
    Ok(result)
}

// Halve the larger operand first. This avoids both an overflowing full phase
// whose half is finite and an underflowed half-duration with a resolved phase.
pub(super) fn phase(rate: Dd, h: f64, fraction: f64) -> Result<Dd, LieStepError> {
    let time = Dd::from_f64(h);
    let fraction = Dd::from_f64(fraction);
    let value = if rate.hi.abs() >= h.abs() {
        (rate * fraction) * time
    } else {
        rate * (time * fraction)
    };
    checked(value, "compensated rotation phase")
}

// Quotient refinement uses fused remainders, not an independently rounded
// denominator*q. The latter can overflow at a finite quotient (MAX/3).
pub(super) fn divide(numerator: Dd, denominator: Dd) -> Result<Dd, LieStepError> {
    let remainder = |value: Dd, quotient: f64| {
        Dd::from_f64((-quotient).mul_add(denominator.hi, value.hi)) + Dd::from_f64(value.lo)
            - Dd::from_f64(quotient) * Dd::from_f64(denominator.lo)
    };
    let first = numerator.hi / denominator.hi;
    if !first.is_finite() || denominator.hi <= 0.0 {
        return Err(invalid("compensated rotation quotient"));
    }
    let residual = checked(remainder(numerator, first), "compensated rotation quotient")?;
    let second = residual.hi / denominator.hi;
    let residual = checked(remainder(residual, second), "compensated rotation quotient")?;
    let third = residual.hi / denominator.hi;
    checked(
        Dd::from_f64(first) + Dd::from_f64(second) + Dd::from_f64(third),
        "compensated rotation quotient",
    )
}

fn product(rotation: So3) -> Product {
    let q = rotation.as_quat();
    [q.w, q.x, q.y, q.z].map(Dd::from_f64)
}

pub(super) fn spatial_vector_to_body(
    rotation: So3,
    value: [Dd; 3],
) -> Result<[Dd; 3], LieStepError> {
    rotate(conjugate(product(rotation)), value)
}

fn conjugate(q: Product) -> Product {
    [q[0], -q[1], -q[2], -q[3]]
}

fn multiply(a: Product, b: Product) -> Product {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}

fn rotate(q: Product, value: Vector) -> Result<Vector, LieStepError> {
    let vector = [q[1], q[2], q[3]];
    let norm = q[0] * q[0] + dot(vector, vector);
    let twice_inverse_norm = divide(Dd::from_f64(2.0), norm)?;
    let t = cross(vector, value).map(|v| v * twice_inverse_norm);
    let correction = cross(vector, t);
    let mut result = [Dd::ZERO; 3];
    for d in 0..3 {
        result[d] = checked(
            value[d] + q[0] * t[d] + correction[d],
            "compensated rigid-body momentum transport",
        )?;
    }
    Ok(result)
}

fn dot(a: Vector, b: Vector) -> Dd {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn finish_attitude(rotation: So3, increment: Product) -> Result<So3, LieStepError> {
    let product = multiply(product(rotation), increment);
    let mut q = [0.0; 4];
    for (target, value) in q.iter_mut().zip(product) {
        *target = checked(value, "compensated rigid-body attitude")?.to_f64();
    }
    Ok(So3::try_from_quat(Quat {
        w: q[0],
        x: q[1],
        y: q[2],
        z: q[3],
    })?)
}

pub(super) fn checked(value: Dd, context: &'static str) -> Result<Dd, LieStepError> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(invalid(context))
    }
}

fn invalid(context: &'static str) -> LieStepError {
    LieStepError::InvalidParameter { context }
}
