//! A maintained rigid ambient flow, with no rate of strain.
//!
//! Linearity of the Stokes equations permits subtracting the ambient rigid
//! translation and rotation before applying the ellipsoid resistance tensors.
//! See the relative-vorticity construction in JFM 949 (2022), A22, section 2.1,
//! <https://doi.org/10.1017/jfm.2022.704>. The angular velocity here is half the
//! curl of the fluid velocity, not the curl itself. Ambient strain would need
//! the corresponding stresslet and work channel and is deliberately excluded.

use super::{Error, finite, norm, positive, vector};
use fs_mbd::Vec3;

/// Prescribed steady velocity `u(x) = u_ref + Omega × (x - x_ref)` in world axes.
///
/// This is an inertial-frame fluid field, not a rotating body coordinate system.
/// No Coriolis, centrifugal pressure, ambient strain or resolved finite vessel
/// is inferred. The Stokes model needs low ambient-rotation Reynolds number as
/// well as low relative translation/spin Reynolds numbers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAmbientFlow {
    reference_position_world_m: Vec3,
    velocity_at_reference_world_m_s: Vec3,
    angular_velocity_world_rad_s: Vec3,
    angular_speed_rad_s: f64,
}

impl RigidAmbientFlow {
    pub fn new(
        reference_position_world_m: Vec3,
        velocity_at_reference_world_m_s: Vec3,
        angular_velocity_world_rad_s: Vec3,
    ) -> Result<Self, Error> {
        vector(reference_position_world_m)?;
        vector(velocity_at_reference_world_m_s)?;
        vector(angular_velocity_world_rad_s)?;
        Ok(Self {
            reference_position_world_m,
            velocity_at_reference_world_m_s,
            angular_velocity_world_rad_s,
            angular_speed_rad_s: norm(angular_velocity_world_rad_s)?,
        })
    }

    pub const fn reference_position_world_m(self) -> Vec3 {
        self.reference_position_world_m
    }
    pub const fn velocity_at_reference_world_m_s(self) -> Vec3 {
        self.velocity_at_reference_world_m_s
    }
    pub const fn angular_velocity_world_rad_s(self) -> Vec3 {
        self.angular_velocity_world_rad_s
    }
    pub(super) const fn angular_speed_rad_s(self) -> f64 {
        self.angular_speed_rad_s
    }

    /// Sample the same maintained fluid field at the body's actual COM.
    pub fn velocity_at(self, position_world_m: Vec3) -> Result<Vec3, Error> {
        vector(position_world_m)?;
        if self.angular_velocity_world_rad_s == Vec3::ZERO {
            // Keep uniform-flow arithmetic independent of the reference point.
            return Ok(self.velocity_at_reference_world_m_s);
        }
        let offset = vector(position_world_m.sub(self.reference_position_world_m))?;
        vector(
            self.velocity_at_reference_world_m_s
                .add(self.angular_velocity_world_rad_s.cross(offset)),
        )
    }
}

/// Eliminate native midpoint translation in the rotating ambient field.
///
/// In shape axes let `K = h*C/m`, `b = h*Omega`, `h = dt/2`, and
/// `rhs = v_before + h*a_external`. Native position satisfies
/// `x_mid = x_before + h*v_mid`, so the ambient sample is
/// `u_mid = u_before + b × v_mid`. Therefore
/// `(I + K - K*[b]x) v_mid = rhs + K*u_before`.
///
/// Row scaling gives `(I - F*[b]x)v = q`, with diagonal `F=K/(1+K)`.
/// Its inverse applied to q is
/// `(q + F*(b×q) + b*((adj(F)*b)·q)) / (1+b·adj(F)*b)`.
/// The denominator is at least one. This is a closed fixed-size elimination
/// of this physical skew coupling, with no matrix inverse, allocation or
/// iterative trajectory owner. The surrounding caller admits `K_i<=1/2`
/// and the ambient rotation increment before this function is evaluated.
pub(super) fn midpoint_velocity_axes(
    relaxation: [f64; 3],
    half_rotation: Vec3,
    free_velocity: Vec3,
    ambient_before: Vec3,
) -> Result<Vec3, Error> {
    let free = [free_velocity.x, free_velocity.y, free_velocity.z];
    let ambient = [ambient_before.x, ambient_before.y, ambient_before.z];
    let mut fraction = [0.0; 3];
    let mut uncoupled = [0.0; 3];
    for i in 0..3 {
        let k = positive(relaxation[i])?;
        let denominator = positive(1.0 + k)?;
        fraction[i] = finite(k / denominator)?;
        uncoupled[i] = finite((free[i] + k * ambient[i]) / denominator)?;
    }
    let q = Vec3::new(uncoupled[0], uncoupled[1], uncoupled[2]);
    let f = fraction;
    let b = vector(half_rotation)?;
    let adjugate_b = vector(Vec3::new(
        f[1] * f[2] * b.x,
        f[0] * f[2] * b.y,
        f[0] * f[1] * b.z,
    ))?;
    let cross = b.cross(q);
    let weighted_cross = vector(Vec3::new(f[0] * cross.x, f[1] * cross.y, f[2] * cross.z))?;
    let denominator = positive(1.0 + b.dot(adjugate_b))?;
    vector(
        q.add(weighted_cross)
            .add(b.scale(finite(adjugate_b.dot(q))?))
            .scale(1.0 / denominator),
    )
}
