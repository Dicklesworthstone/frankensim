//! Rigid carriers with internal axial rotor momentum.
//!
//! The state momentum is the **total** body-frame momentum `L = I omega + h`.
//! `I` is the reduced inertia: retain the rotor transverse inertia but exclude
//! its axial spin inertia. Each component of `h` is built from an absolute
//! rotor axial momentum `p = J (relative_spin + axis.dot(omega))`, not from
//! `J * relative_spin`. This is the standard reaction-wheel reduction (e.g.
//! Basilisk's `mrpSteering` equations (2)--(6)). It prevents axial inertia and
//! motor reaction from being counted twice.
//!
//! With fixed `h`, the reduced Hamiltonian is
//! `H(L) = sum_i (L_i - h_i)^2 / (2 I_i)`. Its three axis subflows are exact:
//! `F = Exp(e_i dt (L_i-h_i)/I_i)`, `R <- R F`, `L <- F^-1 L`.
//! The symmetric x/2, y/2, z, y/2, x/2 composition is second order. Every
//! subflow preserves `R L`, so total spatial angular momentum is preserved
//! to roundoff even when carrier and rotors exchange angular momentum.
//! This is a finite-step approximation, not an exact kinetic-energy claim.

use super::{
    LieStepError, momentum_from_velocity, rigid_body_momentum_step, torque_kick,
    velocity_from_momentum,
};
use fs_ga::{So3, So3Tangent, Vec3};

/// Frame in which the externally applied torque is held over the step.
/// Internal motor torques belong in the rotor momentum change, not here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HeldTorque {
    /// A torque fixed in the rotating principal-axis frame.
    Body(Vec3),
    /// A torque fixed in the inertial frame.
    Space(Vec3),
}

fn finite(vector: Vec3) -> bool {
    [vector.x, vector.y, vector.z]
        .iter()
        .all(|value| value.is_finite())
}

fn admit_rotor(momentum: Vec3) -> Result<(), LieStepError> {
    if !finite(momentum) {
        return Err(LieStepError::InvalidParameter {
            context: "gyrostat rotor momentum",
        });
    }
    Ok(())
}

/// Free gyrostat drift for total momentum in the body's principal frame.
///
/// `rotor_momentum` is held constant in the body frame. With zero stored
/// momentum, this uses the existing rigid-body map exactly, including its
/// exact symmetric-top special case. Nonzero rotors require at most five
/// elementary rotations, no allocation and no iterative solve. A carrier
/// at rest (`momentum == rotor_momentum`) stays at rest exactly.
///
/// # Errors
/// Refuses non-finite momenta, zero/non-finite step size, non-positive or
/// non-finite reduced inertia, and non-representable group operations.
pub fn momentum_step(
    mut rotation: So3,
    mut momentum: Vec3,
    inertia: Vec3,
    rotor_momentum: Vec3,
    h: f64,
) -> Result<(So3, Vec3), LieStepError> {
    admit_rotor(rotor_momentum)?;
    if rotor_momentum == Vec3::new(0.0, 0.0, 0.0) {
        return rigid_body_momentum_step(rotation, momentum, inertia, h);
    }
    if !h.is_finite() || h == 0.0 {
        return Err(LieStepError::InvalidParameter {
            context: "gyrostat step size",
        });
    }
    if ![inertia.x, inertia.y, inertia.z]
        .iter()
        .all(|value| value.is_finite() && *value > 0.0)
    {
        return Err(LieStepError::InvalidParameter {
            context: "gyrostat reduced inertia",
        });
    }
    if !finite(momentum) {
        return Err(LieStepError::InvalidParameter {
            context: "gyrostat total momentum",
        });
    }
    let moments = [inertia.x, inertia.y, inertia.z];
    let stored = [rotor_momentum.x, rotor_momentum.y, rotor_momentum.z];
    for (axis, duration) in [
        (0, h * 0.5),
        (1, h * 0.5),
        (2, h),
        (1, h * 0.5),
        (0, h * 0.5),
    ] {
        let component = [momentum.x, momentum.y, momentum.z][axis];
        let angle = duration * ((component - stored[axis]) / moments[axis]);
        if angle == 0.0 {
            continue;
        }
        let mut tangent = [0.0; 3];
        tangent[axis] = angle;
        let increment = So3::exp(So3Tangent::new(Vec3::new(
            tangent[0], tangent[1], tangent[2],
        )))?;
        momentum = increment.inverse().rotate(momentum)?;
        rotation = rotation.compose(increment)?;
    }
    Ok((rotation, momentum))
}

/// Advance the carrier while rotor momentum changes linearly over a step.
///
/// Half of the rotor momentum change precedes the gyrostat drift, half
/// follows it. These internal transfers leave total momentum untouched;
/// carrier recoil follows from `omega = I^-1 (L - h)`. External torque
/// half-kicks are evaluated in their declared endpoint frames. For a
/// spatially held torque, `R1 (I omega1 + h1) - R0 (I omega0 + h0)` equals
/// `dt * torque` to roundoff. A body-held torque has the trapezoidal spatial
/// impulse `dt/2 * (R0 torque + R1 torque)`.
///
/// The map is symmetric and second order for these held loads/linear rotor
/// momenta. State-dependent motor and aerodynamic torques must be supplied
/// from an appropriate coupled stage to retain that order. No acceleration,
/// rotor-speed, or energy clamp is applied.
///
/// # Errors
/// The refusals of [`momentum_step`], plus non-finite angular velocity,
/// external torque, or resulting angular velocity. Inputs are immutable.
pub fn driven_step(
    rotation: So3,
    omega_body: Vec3,
    inertia: Vec3,
    rotor_before: Vec3,
    rotor_after: Vec3,
    torque: HeldTorque,
    h: f64,
) -> Result<(So3, Vec3), LieStepError> {
    admit_rotor(rotor_before)?;
    admit_rotor(rotor_after)?;
    let midpoint_rotor = rotor_before.scale(0.5) + rotor_after.scale(0.5);
    let initial_momentum = momentum_from_velocity(omega_body, inertia)? + rotor_before;
    let torque_at = |orientation: So3| -> Result<Vec3, LieStepError> {
        match torque {
            HeldTorque::Body(value) => Ok(value),
            HeldTorque::Space(value) => Ok(orientation.inverse().rotate(value)?),
        }
    };
    let kicked = torque_kick(initial_momentum, torque_at(rotation)?, h * 0.5)?;
    let (next_rotation, drifted) = momentum_step(rotation, kicked, inertia, midpoint_rotor, h)?;
    let accepted = torque_kick(drifted, torque_at(next_rotation)?, h * 0.5)?;
    Ok((
        next_rotation,
        velocity_from_momentum(accepted - rotor_after, inertia)?,
    ))
}
