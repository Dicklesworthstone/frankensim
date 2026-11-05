//! Lie-group integrators on `SO(3)`.
//!
//! Group and tangent ownership lives in `fs-ga`; this module owns only time
//! integration algorithms. Body angular velocity advances by right
//! multiplication, while space angular velocity advances by left
//! multiplication. Both conventions are named in the API rather than inferred
//! from an array multiplication order.

use fs_ga::Vec3;

#[path = "lie/error.rs"]
mod error;
pub use error::LieStepError;

#[path = "lie/free.rs"]
mod free;
use free::{momentum_from_velocity, velocity_from_momentum};
pub use free::{rigid_body_momentum_step, rigid_body_step, so3_body_exp_step, so3_space_exp_step};

/// Free rigid rotation retaining impulse tails through momentum and attitude.
#[path = "lie/compensated.rs"]
pub mod compensated;

#[path = "lie/held.rs"]
mod held;
pub use held::{rigid_body_spatial_torque_step, rigid_body_torque_step};

#[path = "lie/gyrostat.rs"]
pub mod gyrostat;

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
