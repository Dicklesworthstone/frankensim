//! Lie-group integrators on `SO(3)`.
//!
//! Group and tangent ownership lives in `fs-ga`; this module owns only time
//! integration algorithms. Body angular velocity advances by right
//! multiplication, while space angular velocity advances by left
//! multiplication. Both conventions are named in the API rather than inferred
//! from an array multiplication order.

use fs_ga::{GaError, So3, So3Tangent, Vec3};

#[path = "lie/gyrostat.rs"]
pub mod gyrostat;

/// Typed failures for `SO(3)` integration steps.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LieStepError {
    /// A group/tangent value failed `fs-ga` validation.
    Group(GaError),
    /// Step size or principal inertia was non-finite, zero, or non-positive as
    /// required by the operation.
    InvalidParameter {
        /// Parameter family that was invalid.
        context: &'static str,
    },
}

impl core::fmt::Display for LieStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Group(error) => write!(f, "SO(3) group refusal: {error}"),
            Self::InvalidParameter { context } => write!(f, "invalid {context}"),
        }
    }
}

impl std::error::Error for LieStepError {}

impl From<GaError> for LieStepError {
    fn from(value: GaError) -> Self {
        Self::Group(value)
    }
}

/// One exponential-map attitude step for body-frame angular velocity:
/// `R_next = R * Exp(h * omega_body)`.
///
/// # Errors
/// Refuses non-finite state, angular velocity, or step size through the
/// canonical `fs-ga` validation boundary.
pub fn so3_body_exp_step(rotation: So3, omega_body: Vec3, h: f64) -> Result<So3, LieStepError> {
    if !h.is_finite() {
        return Err(LieStepError::InvalidParameter {
            context: "SO(3) step size",
        });
    }
    Ok(rotation.body_plus(So3Tangent::new(omega_body.scale(h)))?)
}

/// One exponential-map attitude step for space-frame angular velocity:
/// `R_next = Exp(h * omega_space) * R`.
///
/// # Errors
/// Refuses non-finite state, angular velocity, or step size through the
/// canonical `fs-ga` validation boundary.
pub fn so3_space_exp_step(rotation: So3, omega_space: Vec3, h: f64) -> Result<So3, LieStepError> {
    if !h.is_finite() {
        return Err(LieStepError::InvalidParameter {
            context: "SO(3) step size",
        });
    }
    Ok(rotation.space_plus(So3Tangent::new(omega_space.scale(h)))?)
}

/// Second-order, symmetric free rotation with diagonal principal inertia.
///
/// Exact Hamiltonian subflows transport body momentum by the inverse of the
/// same rotation that advances the attitude. Spatial angular momentum is
/// therefore preserved to roundoff, with bounded symplectic energy error.
/// Spherical and axisymmetric free rotation are exact up to roundoff. There
/// are no allocations, nonlinear solves, or quaternion renormalizations.
///
/// # Errors
/// Refuses non-finite values, a zero step, or non-positive inertia.
pub fn rigid_body_step(
    rotation: So3,
    omega_body: Vec3,
    inertia: Vec3,
    h: f64,
) -> Result<(So3, Vec3), LieStepError> {
    let momentum = momentum_from_velocity(omega_body, inertia)?;
    let (next_rotation, next_momentum) = rigid_body_momentum_step(rotation, momentum, inertia, h)?;
    Ok((
        next_rotation,
        velocity_from_momentum(next_momentum, inertia)?,
    ))
}

fn momentum_from_velocity(omega_body: Vec3, inertia: Vec3) -> Result<Vec3, LieStepError> {
    if ![omega_body.x, omega_body.y, omega_body.z]
        .iter()
        .all(|value| value.is_finite())
    {
        return Err(LieStepError::InvalidParameter {
            context: "rigid-body angular velocity",
        });
    }
    Ok(Vec3::new(
        inertia.x * omega_body.x,
        inertia.y * omega_body.y,
        inertia.z * omega_body.z,
    ))
}

fn velocity_from_momentum(momentum: Vec3, inertia: Vec3) -> Result<Vec3, LieStepError> {
    let next_omega = Vec3::new(
        momentum.x / inertia.x,
        momentum.y / inertia.y,
        momentum.z / inertia.z,
    );
    if ![next_omega.x, next_omega.y, next_omega.z]
        .iter()
        .all(|value| value.is_finite())
    {
        return Err(LieStepError::InvalidParameter {
            context: "rigid-body angular velocity",
        });
    }
    Ok(next_omega)
}

fn torque_kick(momentum: Vec3, torque: Vec3, h: f64) -> Result<Vec3, LieStepError> {
    if ![torque.x, torque.y, torque.z]
        .iter()
        .all(|value| value.is_finite())
    {
        return Err(LieStepError::InvalidParameter {
            context: "rigid-body torque",
        });
    }
    Ok(Vec3::new(
        h.mul_add(torque.x, momentum.x),
        h.mul_add(torque.y, momentum.y),
        h.mul_add(torque.z, momentum.z),
    ))
}

/// Second-order rotation under torque held constant in the body frame.
///
/// Both `omega_body` and `torque_body` use the principal-axis frame, and the
/// torque is about the center of mass. Symmetric momentum half-kicks enclose
/// the free Hamiltonian drift in [`rigid_body_momentum_step`]. The second
/// kick acts in the accepted body's frame, not its initial spatial frame.
/// The map is time reversible for this held load. No nonlinear solve or
/// allocation is needed; input state is never mutated on refusal.
///
/// This is a held-load integrator: sampling a state-dependent torque once
/// per step does not preserve second-order accuracy of that coupled model.
/// A constant *spatial* torque instead requires
/// [`rigid_body_spatial_torque_step`].
///
/// # Errors
/// Refuses the same invalid states/steps as [`rigid_body_step`], non-finite
/// torque, and non-representable kicked momentum or angular velocity.
pub fn rigid_body_torque_step(
    rotation: So3,
    omega_body: Vec3,
    inertia: Vec3,
    torque_body: Vec3,
    h: f64,
) -> Result<(So3, Vec3), LieStepError> {
    let momentum = momentum_from_velocity(omega_body, inertia)?;
    let half = 0.5 * h;
    let kicked = torque_kick(momentum, torque_body, half)?;
    let (next_rotation, drifted) = rigid_body_momentum_step(rotation, kicked, inertia, h)?;
    let next_momentum = torque_kick(drifted, torque_body, half)?;
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
/// This symmetric, allocation-free held-load map is time reversible. It
/// makes no second-order claim when a time- or state-dependent torque is
/// frozen once per step.
///
/// # Errors
/// Refuses the same invalid states/steps as [`rigid_body_step`], non-finite
/// torque, and non-representable kicked momentum or angular velocity.
pub fn rigid_body_spatial_torque_step(
    rotation: So3,
    omega_body: Vec3,
    inertia: Vec3,
    torque_space: Vec3,
    h: f64,
) -> Result<(So3, Vec3), LieStepError> {
    let momentum = momentum_from_velocity(omega_body, inertia)?;
    let half = 0.5 * h;
    let torque_start = rotation.inverse().rotate(torque_space)?;
    let kicked = torque_kick(momentum, torque_start, half)?;
    let (next_rotation, drifted) = rigid_body_momentum_step(rotation, kicked, inertia, h)?;
    let torque_end = next_rotation.inverse().rotate(torque_space)?;
    let next_momentum = torque_kick(drifted, torque_end, half)?;
    Ok((
        next_rotation,
        velocity_from_momentum(next_momentum, inertia)?,
    ))
}

/// Momentum-state version of [`rigid_body_step`].
///
/// The input and output momentum are about the center of mass in the body's
/// principal-axis frame. Keeping momentum as the state avoids repeated
/// velocity-to-momentum conversions in long conservative trajectories.
///
/// For the median principal inertia `I0`, split
/// `H = |L|²/(2 I0) + (1/Ia - 1/I0) La²/2 + (1/Ib - 1/I0) Lb²/2`.
/// The spherical part commutes with both residual flows. A symmetric
/// `a(h/2), b(h), a(h/2)` composition followed by the spherical flow is
/// second order and time reversible. Each residual pairs `R <- R F` with
/// `L <- F⁻¹ L`; the spherical flow rotates about `L` itself. This is the
/// rigid-body Hamiltonian splitting of Dullweber, Leimkuhler and McLachlan
/// (J. Chem. Phys. 107, 5840, 1997), with the commuting spherical part
/// extracted to integrate symmetric tops exactly.
///
/// # Errors
/// Refuses non-finite momentum, zero/non-finite step size, non-positive
/// principal inertia, and non-representable rotation or momentum results.
pub fn rigid_body_momentum_step(
    mut rotation: So3,
    mut momentum: Vec3,
    inertia: Vec3,
    h: f64,
) -> Result<(So3, Vec3), LieStepError> {
    if !h.is_finite() || h == 0.0 {
        return Err(LieStepError::InvalidParameter {
            context: "rigid-body step size",
        });
    }
    if ![inertia.x, inertia.y, inertia.z]
        .iter()
        .all(|value| value.is_finite() && *value > 0.0)
    {
        return Err(LieStepError::InvalidParameter {
            context: "rigid-body principal inertia",
        });
    }
    if ![momentum.x, momentum.y, momentum.z]
        .iter()
        .all(|value| value.is_finite())
    {
        return Err(LieStepError::InvalidParameter {
            context: "rigid-body angular momentum",
        });
    }

    let moments = [inertia.x, inertia.y, inertia.z];
    let mut axes = [0, 1, 2];
    axes.sort_unstable_by(|&a, &b| moments[a].total_cmp(&moments[b]));
    let reference = moments[axes[1]];
    let a = axes[0];
    let b = axes[2];
    let mut axis_step = |axis: usize, step: f64| -> Result<(), LieStepError> {
        let component = [momentum.x, momentum.y, momentum.z][axis];
        let angle = step * (component / moments[axis] - component / reference);
        if angle == 0.0 {
            return Ok(());
        }
        let mut angular = [0.0; 3];
        angular[axis] = angle;
        let increment = So3::exp(So3Tangent::new(Vec3::new(
            angular[0], angular[1], angular[2],
        )))?;
        momentum = increment.inverse().rotate(momentum)?;
        rotation = rotation.compose(increment)?;
        Ok(())
    };
    if moments[a] == reference {
        axis_step(b, h)?;
    } else if moments[b] == reference {
        axis_step(a, h)?;
    } else {
        axis_step(a, 0.5 * h)?;
        axis_step(b, h)?;
        axis_step(a, 0.5 * h)?;
    }
    rotation = so3_body_exp_step(
        rotation,
        Vec3::new(
            momentum.x / reference,
            momentum.y / reference,
            momentum.z / reference,
        ),
        h,
    )?;
    if ![momentum.x, momentum.y, momentum.z]
        .iter()
        .all(|value| value.is_finite())
    {
        return Err(LieStepError::InvalidParameter {
            context: "rigid-body angular momentum",
        });
    }
    Ok((rotation, momentum))
}
