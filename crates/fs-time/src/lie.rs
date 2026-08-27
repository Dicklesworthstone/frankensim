//! Lie-group integrators on `SO(3)`.
//!
//! Group and tangent ownership lives in `fs-ga`; this module owns only time
//! integration algorithms. Body angular velocity advances by right
//! multiplication, while space angular velocity advances by left
//! multiplication. Both conventions are named in the API rather than inferred
//! from an array multiplication order.

use fs_ga::{GaError, So3, So3Tangent, Vec3};

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

/// One commutator-free CG2 midpoint step of a free rigid body with diagonal
/// principal inertia. Euler's equations advance body angular velocity and the
/// attitude uses the same midpoint velocity in a right/body group step.
///
/// # Errors
/// Refuses non-finite values, a zero step, or non-positive inertia.
pub fn rigid_body_step(
    rotation: So3,
    omega_body: Vec3,
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
    if ![omega_body.x, omega_body.y, omega_body.z]
        .iter()
        .all(|value| value.is_finite())
    {
        return Err(LieStepError::InvalidParameter {
            context: "rigid-body angular velocity",
        });
    }

    let torque_free = |omega: Vec3| -> Vec3 {
        let momentum = Vec3::new(
            inertia.x * omega.x,
            inertia.y * omega.y,
            inertia.z * omega.z,
        );
        Vec3::new(
            momentum.y.mul_add(omega.z, -(momentum.z * omega.y)) / inertia.x,
            momentum.z.mul_add(omega.x, -(momentum.x * omega.z)) / inertia.y,
            momentum.x.mul_add(omega.y, -(momentum.y * omega.x)) / inertia.z,
        )
    };
    let first_slope = torque_free(omega_body);
    let midpoint = Vec3::new(
        (0.5 * h).mul_add(first_slope.x, omega_body.x),
        (0.5 * h).mul_add(first_slope.y, omega_body.y),
        (0.5 * h).mul_add(first_slope.z, omega_body.z),
    );
    let second_slope = torque_free(midpoint);
    let next_omega = Vec3::new(
        h.mul_add(second_slope.x, omega_body.x),
        h.mul_add(second_slope.y, omega_body.y),
        h.mul_add(second_slope.z, omega_body.z),
    );
    let next_rotation = so3_body_exp_step(rotation, midpoint, h)?;
    Ok((next_rotation, next_omega))
}
