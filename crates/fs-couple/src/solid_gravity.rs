//! Conservative rigid-body gravity fields on the original deforming solid.
//!
//! This adapter binds a field's explicit inertial frame to hyper3's quadrature
//! callback. The shared field owns acceleration, potential and their analytic
//! spatial derivative; the solid owns current positions, reference material
//! mass, quadrature, force assembly, Newton stiffness and support reactions.
//! There is no second gravity law, mass inventory or equilibrium solver here.

use fs_mbd::Vec3;
use fs_mbd::gravity::{GravityError, PointGravityField};
use fs_solid::HyperTetError;
use fs_solid::hyper3::gravity::GravitySample3;

/// Bind one stationary conservative field to the solid's world coordinates.
///
/// Pass the returned callback to `HyperTetProblem::solve_with_gravity`,
/// `solve_with_gravity_and_pressure`, or a prepared problem's `gravity_load`.
/// Solid positions must be expressed in this same frame and in meters. The
/// field's frame and domain checks remain authoritative. No rotating-frame
/// inertial force, moving-source work or solid self-gravity is inferred.
///
/// The binding borrows the source and allocates nothing. Each solid quadrature
/// sample requests all three conservative quantities from that source once.
pub fn bind_field<F: PointGravityField + ?Sized>(
    field: &F,
    frame: u64,
) -> impl Fn([f64; 3]) -> Result<GravitySample3, HyperTetError> + '_ {
    move |position| {
        let sample = field
            .sample_point(frame, Vec3::new(position[0], position[1], position[2]))
            .and_then(|sample| sample.validate())
            .map_err(|error| match error {
                GravityError::Cancelled => HyperTetError::Cancelled,
                error => HyperTetError::InvalidInput {
                    what: error.to_string(),
                },
            })?;
        Ok(GravitySample3 {
            potential_j_kg: sample.specific_potential_j_kg,
            acceleration_m_s2: [
                sample.acceleration_world_m_s2.x,
                sample.acceleration_world_m_s2.y,
                sample.acceleration_world_m_s2.z,
            ],
            acceleration_gradient_s2: std::array::from_fn(|row| {
                std::array::from_fn(|column| sample.acceleration_gradient_s2.m[3 * row + column])
            }),
        })
    }
}
