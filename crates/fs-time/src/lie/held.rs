//! Held torque consumers of the compensated free rigid-body flow.

use super::{LieStepError, compensated};
use fs_ga::{So3, Vec3};
use fs_math::dd::Dd;

/// Second-order rotation under torque held constant in the body frame.
///
/// Both `omega_body` and `torque_body` use the principal-axis frame, and the
/// torque is about the center of mass. Symmetric momentum half-kicks enclose
/// the free Hamiltonian drift in [`compensated::rigid_body_momentum_step`].
/// Momentum products and both kicks retain DD tails through that drift and
/// the endpoint velocity conversion. The retained midpoint impulse drives
/// attitude as well as spin; only the accepted public state rounds to f64.
/// The second kick acts in the accepted body's frame, not its initial
/// spatial frame. No nonlinear solve or allocation is needed; input state
/// is never mutated on refusal.
/// Exactly zero torque retains the established scalar free-drift result
/// bit for bit, including its admission/refusal behavior.
///
/// This map is time reversible for a held load. Sampling a state-dependent
/// torque once per step does not preserve second-order accuracy of that
/// coupled model. A constant *spatial* torque instead requires
/// [`rigid_body_spatial_torque_step`].
///
/// # Errors
/// Refuses nonfinite/zero duration, nonfinite angular velocity or torque,
/// nonpositive/nonfinite inertia, nonrepresentable momentum/velocity or
/// free-rotation arithmetic, and failed canonical `So3` admission.
pub fn rigid_body_torque_step(
    rotation: So3,
    omega_body: Vec3,
    inertia: Vec3,
    torque_body: Vec3,
    h: f64,
) -> Result<(So3, Vec3), LieStepError> {
    if torque_body == Vec3::default() {
        return super::rigid_body_step(rotation, omega_body, inertia, h);
    }
    let momentum = momentum_from_velocity(omega_body, inertia)?;
    let torque = finite_vector(torque_body, "rigid-body torque")?;
    let kicked = half_kick(momentum, torque, h)?;
    let (next_rotation, drifted) =
        compensated::rigid_body_momentum_step(rotation, kicked, inertia, h)?;
    let next_momentum = half_kick(drifted, torque, h)?;
    Ok((
        next_rotation,
        velocity_from_momentum(next_momentum, inertia)?,
    ))
}

/// Second-order rotation under torque held constant in the spatial frame.
///
/// `torque_space` is about the center of mass; `omega_body` and `inertia`
/// use the body's principal axes. Each momentum half-kick rotates this
/// torque into the body frame at that endpoint. Together with the matched
/// free drift, this gives `R_next I omega_next - R I omega = h torque_space`
/// to roundoff for any admitted step. Reusing the initial body-frame torque
/// at both endpoints would instead apply a rotating, body-attached load.
///
/// The impulse tails, frame transforms and momentum products remain in DD
/// through the free attitude evolution and final velocity conversion.
/// This symmetric, allocation-free held-load map is time reversible. It
/// makes no second-order claim when a time- or state-dependent torque is
/// frozen once per step.
/// Exactly zero torque retains the established scalar free-drift result
/// bit for bit, including its admission/refusal behavior.
///
/// # Errors
/// The same finite-parameter, derived-arithmetic and group refusals as
/// [`rigid_body_torque_step`].
pub fn rigid_body_spatial_torque_step(
    rotation: So3,
    omega_body: Vec3,
    inertia: Vec3,
    torque_space: Vec3,
    h: f64,
) -> Result<(So3, Vec3), LieStepError> {
    if torque_space == Vec3::default() {
        return super::rigid_body_step(rotation, omega_body, inertia, h);
    }
    let momentum = momentum_from_velocity(omega_body, inertia)?;
    let torque = finite_vector(torque_space, "rigid-body spatial torque")?;
    let torque_start = compensated::spatial_vector_to_body(rotation, torque)?;
    let kicked = half_kick(momentum, torque_start, h)?;
    let (next_rotation, drifted) =
        compensated::rigid_body_momentum_step(rotation, kicked, inertia, h)?;
    let torque_end = compensated::spatial_vector_to_body(next_rotation, torque)?;
    let next_momentum = half_kick(drifted, torque_end, h)?;
    Ok((
        next_rotation,
        velocity_from_momentum(next_momentum, inertia)?,
    ))
}

fn finite_vector(vector: Vec3, context: &'static str) -> Result<[Dd; 3], LieStepError> {
    let mut result = [Dd::ZERO; 3];
    for (target, value) in result.iter_mut().zip([vector.x, vector.y, vector.z]) {
        *target = compensated::checked(Dd::from_f64(value), context)?;
    }
    Ok(result)
}

fn momentum_from_velocity(omega: Vec3, inertia: Vec3) -> Result<[Dd; 3], LieStepError> {
    let omega = finite_vector(omega, "rigid-body angular velocity")?;
    let moments = [inertia.x, inertia.y, inertia.z];
    if !moments
        .into_iter()
        .all(|value| value.is_finite() && value > 0.0)
    {
        return Err(LieStepError::InvalidParameter {
            context: "rigid-body principal inertia",
        });
    }
    let mut result = [Dd::ZERO; 3];
    for d in 0..3 {
        result[d] = compensated::checked(
            Dd::from_f64(moments[d]) * omega[d],
            "rigid-body angular momentum",
        )?;
    }
    Ok(result)
}

fn velocity_from_momentum(momentum: [Dd; 3], inertia: Vec3) -> Result<Vec3, LieStepError> {
    let moments = [inertia.x, inertia.y, inertia.z];
    let mut result = [0.0; 3];
    for d in 0..3 {
        result[d] = compensated::divide(momentum[d], Dd::from_f64(moments[d]))?.to_f64();
    }
    Ok(Vec3::new(result[0], result[1], result[2]))
}

fn half_kick(momentum: [Dd; 3], torque: [Dd; 3], h: f64) -> Result<[Dd; 3], LieStepError> {
    let mut result = [Dd::ZERO; 3];
    for d in 0..3 {
        result[d] = compensated::checked(
            momentum[d] + compensated::phase(torque[d], h, 0.5)?,
            "rigid-body kicked angular momentum",
        )?;
    }
    Ok(result)
}
