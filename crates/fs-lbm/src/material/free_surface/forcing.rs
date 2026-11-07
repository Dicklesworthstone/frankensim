//! Prescribed acceleration fields and vessel motion on the sourced VOF trajectory.

use super::*;
use crate::d3q19::{Coriolis3, RigidFrame3};
use crate::d3q19::freesurface3::BodyAccelerationField3;

/// Prescribed affine body acceleration in the fixed physical grid axes.
///
/// `a(x) = acceleration_m_s2 + gradient_per_s2 * (x - origin_m)`.
/// The grid origin is zero and native cell centres are `(i + 1/2) * dx`.
/// This is a position-dependent load, not a complete rotating frame, Coriolis
/// model, source-density buoyancy law or thermodynamic expansion model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialAffineAcceleration3 {
    /// Fixed reference position in physical grid coordinates [m].
    pub origin_m: [f64; 3],
    /// Acceleration at that position [m/s²].
    pub acceleration_m_s2: [f64; 3],
    /// Row-major spatial derivative [1/s²], component then coordinate.
    pub gradient_per_s2: [[f64; 3]; 3],
}

impl MaterialAffineAcceleration3 {
    fn lattice_field(
        self,
        scale: &MaterialFreeSurfaceScaling<'_>,
    ) -> Result<BodyAccelerationField3, MaterialLbmError> {
        let mut origin = [0.0; 3];
        let mut acceleration = [0.0; 3];
        let mut gradient = [[0.0; 3]; 3];
        for axis in 0..3 {
            origin[axis] = converted(self.origin_m[axis], self.origin_m[axis] / scale.cell_size_m)?;
            acceleration[axis] = converted(
                self.acceleration_m_s2[axis],
                scale.gravity_lattice(self.acceleration_m_s2[axis])?,
            )?;
            for coordinate in 0..3 {
                let value = self.gradient_per_s2[axis][coordinate];
                gradient[axis][coordinate] =
                    converted(value, (value * scale.timestep_s) * scale.timestep_s)?;
            }
        }
        Ok(BodyAccelerationField3::Affine {
            origin,
            acceleration,
            gradient,
        })
    }
}

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

impl MaterialRigidFrame3 {
    // Frozen and evolving-viscosity owners retain one SI conversion convention.
    pub(super) fn lattice_forcing(
        self,
        scale: &MaterialFreeSurfaceScaling<'_>,
        body_acceleration_m_s2: [f64; 3],
    ) -> Result<(RigidFrame3, [f64; 3]), MaterialLbmError>{
        let mut native = RigidFrame3::default();
        let mut acceleration = [0.0; 3];
        for axis in 0..3 {
            native.origin[axis] =
                converted(self.origin_m[axis], self.origin_m[axis] / scale.cell_size_m)?;
            native.angular_velocity[axis] = converted(
                self.angular_velocity_rad_s[axis],
                self.angular_velocity_rad_s[axis] * scale.timestep_s,
            )?;
            native.angular_acceleration[axis] = converted(
                self.angular_acceleration_rad_s2[axis],
                (self.angular_acceleration_rad_s2[axis] * scale.timestep_s) * scale.timestep_s,
            )?;
            native.origin_acceleration[axis] = converted(
                self.origin_acceleration_m_s2[axis],
                scale.gravity_lattice(self.origin_acceleration_m_s2[axis])?,
            )?;
            acceleration[axis] = converted(
                body_acceleration_m_s2[axis],
                scale.gravity_lattice(body_acceleration_m_s2[axis])?,
            )?;
        }
        Ok((native, acceleration))
    }
}

/// Uniform local Coriolis model in the original material trajectory's SI units.
///
/// Effective gravity already includes the background centrifugal potential.
/// This is not complete vessel kinematics and does not transform initial data.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MaterialCoriolis3 {
    /// Right-handed angular velocity [rad/s] in the grid's velocity axes.
    pub angular_velocity_rad_s: [f64; 3],
}

impl MaterialCoriolis3 {
    pub(super) fn lattice_forcing(
        self,
        scale: &MaterialFreeSurfaceScaling<'_>,
        effective_acceleration_m_s2: [f64; 3],
    ) -> Result<(Coriolis3, [f64; 3]), MaterialLbmError> {
        let mut native = Coriolis3::default();
        let mut acceleration = [0.0; 3];
        for axis in 0..3 {
            native.angular_velocity[axis] = converted(
                self.angular_velocity_rad_s[axis],
                self.angular_velocity_rad_s[axis] * scale.timestep_s,
            )?;
            acceleration[axis] = converted(
                effective_acceleration_m_s2[axis],
                scale.gravity_lattice(effective_acceleration_m_s2[axis])?,
            )?;
        }
        Ok((native, acceleration))
    }
}

impl MaterialFreeSurface3<'_> {
    /// Advance one SI timestep with local Coriolis and effective acceleration.
    ///
    /// The native skew midpoint source uses `omega_SI * dt` and
    /// `a_SI * dt²/dx`, with the same original material/scaling identity.
    /// Preparation and the native VOF step publish as one transaction; invalid
    /// conversion, cell admission or cancellation leaves accepted flow intact.
    /// Refresh the source every tick. This replaces initial held gravity and
    /// adds no thermal source or independent body-acceleration owner.
    pub fn step_with_coriolis_si(
        &mut self,
        rotation: MaterialCoriolis3,
        effective_acceleration_m_s2: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), MaterialLbmError> {
        if gate.is_requested() {
            return Err(crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled.into());
        }
        let (native, acceleration) =
            rotation.lattice_forcing(&self.scaling, effective_acceleration_m_s2)?;
        native.step_free_surface(&mut self.lattice, acceleration, gate)?;
        Ok(())
    }

    /// Bind a fixed-grid affine acceleration without advancing time or applying
    /// an impulse. Density, mass, physical velocity and force-corrected stress
    /// retain the same source/Guo conventions as [`Self::rebind_gravity`].
    /// The affine coefficients are converted once; native cells and future
    /// newly wet cells evaluate the same field at their actual positions.
    pub fn rebind_affine_acceleration(
        &mut self,
        field: MaterialAffineAcceleration3,
        gate: &CancelGate,
    ) -> Result<(), MaterialLbmError> {
        let field = field.lattice_field(&self.scaling)?;
        self.lattice.rebind_body_acceleration_field(field, gate)?;
        self.scaling.density_dependent_gravity = true;
        Ok(())
    }

    /// Advance to the next prescribed affine field using both actual-density
    /// force endpoints. `a0` scales by dt²/dx and the spatial gradient by dt².
    /// The field stays in the original physical grid coordinates, including
    /// through piston relocation; no per-cell material query or new field
    /// buffer is needed. The complete native candidate publishes atomically.
    pub fn advance_affine_acceleration(
        &mut self,
        next: MaterialAffineAcceleration3,
        gate: &CancelGate,
    ) -> Result<(), MaterialLbmError> {
        let next = next.lattice_field(&self.scaling)?;
        self.lattice.advance_body_acceleration_field(next, gate)?;
        Ok(())
    }

    /// Compose the affine endpoint with absolute gas-link pressure [Pa] about
    /// the original fixed pressure reference. The actual next liquid density
    /// multiplies the field at each interface, including new interfaces. The
    /// same pressure, source-frame and cancellation admission owns publication.
    pub fn advance_affine_acceleration_with_gas_link_pressure(
        &mut self,
        next: MaterialAffineAcceleration3,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn(MaterialGasPressureLink3) -> f64,
    ) -> Result<(), MaterialLbmError> {
        let next = next.lattice_field(&self.scaling)?;
        let reference = self
            .scaling
            .reference_pressure_pa
            .ok_or(MaterialLbmError::Invalid(
                "absolute gas loading requires an explicit fixed pressure reference",
            ))?;
        self.lattice
            .advance_body_acceleration_field_with_gas_link_pressure(
                next,
                gate,
                maximum_density_excursion,
                super::pressure::lattice_link_pressure_field(
                    reference,
                    self.scaling.cell_size_m,
                    self.scaling.pressure_unit_pa,
                    pressure_pa,
                ),
            )?;
        Ok(())
    }

/// Advance one physical timestep from the accepted gravity/acceleration
    /// to `next_gravity_m_s2`, retaining both actual-density force endpoints.
    ///
    /// Material-bound construction already selects density-aware gravity;
    /// [`Self::rebind_gravity`] can replace it without advancing time.
    /// For homogeneous constant-density interior flow the velocity
    /// increment is `dt * (gravity_before + gravity_after) / 2` at any admitted
    /// sourced density and fixed numerical pressure-density level. When density
    /// evolves, the momentum source uses each endpoint's own `rho * gravity`;
    /// advection and pressure also contribute to the resulting velocity.
    ///
    /// SI conversion and the complete native candidate must succeed before
    /// either liquid or the next acceleration publishes. Material coefficients,
    /// pressure reference, dx/dt and the original speed limits stay attached.
    /// The vectors are in the same fixed axes; rotating-frame inertia belongs
    /// to [`Self::step_in_rigid_frame_si`]. This is prescribed loading, with no
    /// added vessel evolution, thermal feedback or mechanical-work heating.
    pub fn advance_gravity(
        &mut self,
        next_gravity_m_s2: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), MaterialLbmError> {
        let next = self.scaling.gravity_vector_lattice(next_gravity_m_s2)?;
        self.lattice.advance_body_acceleration(next, gate)?;
        Ok(())
    }

    /// Evolve the same gravity endpoints under absolute gas-link pressure [Pa].
    ///
    /// The original fixed pressure reference and physical link coordinates are
    /// those of [`Self::step_with_gas_link_pressure`]. Pressure reconstruction
    /// uses the actual next liquid density and acceleration; it cannot retain
    /// the old gravity force merely because pressure changed the carrier mass.
    /// Load conversion, pressure admission, native flow and cancellation are
    /// one transaction. The callback describes one frozen pressure field.
    pub fn advance_gravity_with_gas_link_pressure(
        &mut self,
        next_gravity_m_s2: [f64; 3],
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn(MaterialGasPressureLink3) -> f64,
    ) -> Result<(), MaterialLbmError> {
        let next = self.scaling.gravity_vector_lattice(next_gravity_m_s2)?;
        let reference = self
            .scaling
            .reference_pressure_pa
            .ok_or(MaterialLbmError::Invalid(
                "absolute gas loading requires an explicit fixed pressure reference",
            ))?;
        self.lattice
            .advance_body_acceleration_with_gas_link_pressure(
                next,
                gate,
                maximum_density_excursion,
                super::pressure::lattice_link_pressure_field(
                    reference,
                    self.scaling.cell_size_m,
                    self.scaling.pressure_unit_pa,
                    pressure_pa,
                ),
            )?;
        Ok(())
    }

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
    /// numerical density, as does the constructor's gravity. No thermal/buoyancy
    /// owner is implicitly composed here.
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
        let (native, acceleration) =
            frame.lattice_forcing(&self.scaling, body_acceleration_m_s2)?;
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
