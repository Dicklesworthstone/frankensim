//! Prescribed vessel acceleration through the source-owned sparse VOF trajectory.

use super::*;
use crate::d3q19::RigidFrame3;

/// Prescribed motion of a rigid vessel, with all vectors in its moving axes.
/// The material grid has physical origin zero and cell centres `(i+0.5)*dx`.
/// This descriptor changes no geometry or initial velocity representation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MaterialRigidFrame3 {
    /// Rotation origin in the vessel's physical coordinates [m].
    pub origin_m: [f64; 3],
    /// Right-handed angular velocity [rad/s].
    pub angular_velocity_rad_s: [f64; 3],
    /// Right-handed angular acceleration [rad/s²].
    pub angular_acceleration_rad_s2: [f64; 3],
    /// Inertial acceleration of the origin, resolved in vessel axes [m/s²].
    pub origin_acceleration_m_s2: [f64; 3],
}

impl MaterialFreeSurface3<'_> {
    /// Advance the original liquid by one SI timestep in an accelerating vessel.
    ///
    /// `body_acceleration_m_s2` includes gravity, already rotated into vessel
    /// axes. Supply prescribed motion coefficients at the interval midpoint.
    /// Relative fluid velocity, material coefficients, capillarity, wetting and
    /// the original pressure reconstruction retain their existing owners.
    ///
    /// The complete held force is replaced, so initial gravity is not counted
    /// twice. Translation, centrifugal/Euler acceleration and the shared
    /// energy-neutral midpoint Coriolis source use each wet cell's CURRENT
    /// numerical density. This differs from the constructor's frozen rho0*g
    /// reference force. No thermal/buoyancy owner is implicitly composed here.
    ///
    /// Physical units always come from this trajectory's original source and
    /// dx/dt. Conversion, force preparation and native evolution take place in
    /// one private candidate: failure or cancellation leaves accepted liquid,
    /// forces, VOF mass and topology unchanged. The caller owns the external
    /// clock. This does not infer finite-inertia vessel or gas dynamics, prove
    /// full-flow energy conservation, or relax the existing resolution limits.
    pub fn step_in_rigid_frame_si(
        &mut self,
        frame: MaterialRigidFrame3,
        body_acceleration_m_s2: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), MaterialLbmError> {
        if gate.is_requested() {
            return Err(crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled.into());
        }
        let scale = &self.scaling;
        let mut native = RigidFrame3::default();
        let mut acceleration = [0.0; 3];
        for axis in 0..3 {
            native.origin[axis] = converted(
                frame.origin_m[axis],
                frame.origin_m[axis] / scale.cell_size_m,
            )?;
            native.angular_velocity[axis] = converted(
                frame.angular_velocity_rad_s[axis],
                frame.angular_velocity_rad_s[axis] * scale.timestep_s,
            )?;
            native.angular_acceleration[axis] = converted(
                frame.angular_acceleration_rad_s2[axis],
                (frame.angular_acceleration_rad_s2[axis] * scale.timestep_s) * scale.timestep_s,
            )?;
            native.origin_acceleration[axis] = converted(
                frame.origin_acceleration_m_s2[axis],
                scale.gravity_lattice(frame.origin_acceleration_m_s2[axis])?,
            )?;
            acceleration[axis] = converted(
                body_acceleration_m_s2[axis],
                scale.gravity_lattice(body_acceleration_m_s2[axis])?,
            )?;
        }
        native.step_free_surface(&mut self.lattice, acceleration, gate)?;
        Ok(())
    }
}

fn converted(source: f64, lattice: f64) -> Result<f64, MaterialLbmError> {
    if !source.is_finite() || !lattice.is_finite() || (source != 0.0 && lattice == 0.0) {
        return Err(MaterialLbmError::Invalid(
            "rigid-frame motion is not representable in the original lattice units",
        ));
    }
    Ok(lattice)
}
