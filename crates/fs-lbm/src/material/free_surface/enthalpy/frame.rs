//! Prescribed vessel acceleration on the original thermal mass/enthalpy tick.

use super::*;
use crate::d3q19::RigidFrame3;
use crate::material::free_surface::MaterialRigidFrame3;

#[derive(Clone, Copy)]
pub(super) struct ThermalRigidFrame {
    native: RigidFrame3,
    acceleration: [f64; 3],
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Advance heat, sourced viscosity and moving liquid in a prescribed frame.
    ///
    /// Frame coefficients are interval-midpoint SI values in the moving axes,
    /// as in [`MaterialFreeSurface3::step_in_rigid_frame_si`]. Body acceleration
    /// includes gravity. Initial fluid/wall velocities are already relative to
    /// this frame; the method does not transform an inertial initial state.
    ///
    /// Heat refreshes the original material/rheology law before ONE collision.
    /// The original midpoint Coriolis solve includes all retained surface/body
    /// forces in that same source step. Its actual native mass routes transport
    /// extensive enthalpy once, including interface conversion. Only one private
    /// sparse fluid candidate is cloned. All source, speed, transport, mass,
    /// enthalpy and cancellation checks precede joint state/clock acceptance.
    ///
    /// Constant-density thermal/rheological flow admits arbitrary-axis rotation.
    /// A bound Boussinesq source admits translation only: its original density
    /// departure multiplies `gravity - origin_acceleration`, so free fall does
    /// not create buoyancy. Current numerical density converts this acceleration
    /// to force, matching the isothermal frame owner; reference caloric mass,
    /// viscosity scaling and the source pressure datum remain unchanged.
    /// Rotating Boussinesq flow, including angular acceleration at zero speed,
    /// refuses: its omitted density/inertia couplings need a separate model.
    /// See Journal of Fluid Mechanics (2021), <https://doi.org/10.1017/jfm.2021.640>.
    ///
    /// Like the original additional-body-force entries, the accepted endpoint
    /// restores its thermal Guo force representation without another impulse.
    /// Supply the prescribed frame again on EVERY frame-driven interval.
    /// Mechanical frame work is not heat. Thermal storage remains the existing
    /// first-order source/upwind split; this adds neither compressive heating,
    /// finite-inertia vessel motion nor a total mechanical-plus-thermal audit.
    /// Bound electrical/surface-heat owners still require their full timelines.
    pub fn step_in_rigid_frame_si(
        &mut self,
        frame: MaterialRigidFrame3,
        body_acceleration_m_s2: [f64; 3],
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceEnthalpyReport, Error> {
        poll(gate)?;
        let frame = self.thermal_frame(frame, body_acceleration_m_s2)?;
        let (heated, external, absolute_heat) = self.prepare_heat(heat, gate)?;
        let tick = self.prepare_enthalpy_tick_with_frame(
            heated,
            external,
            0.0,
            0.0,
            absolute_heat,
            None,
            None,
            None,
            None,
            None,
            Some(frame),
            gate,
        )?;
        Ok(self.accept_enthalpy_tick(tick).0)
    }

    /// The same prescribed-frame tick with original molecular strain heating.
    ///
    /// The force-corrected initial strain and its actual source viscosity supply
    /// `2 mu D:D` once, before thermal feedback and flow. Uniform acceleration,
    /// Coriolis work, boundary work and numerical kinetic loss are not converted
    /// to heat. Authored heat must exclude this explicitly selected molecular
    /// source. Strong source feedback needs timestep refinement, as for
    /// [`Self::step_with_viscous_heating`]. All state publishes atomically.
    pub fn step_in_rigid_frame_with_viscous_heating_si(
        &mut self,
        frame: MaterialRigidFrame3,
        body_acceleration_m_s2: [f64; 3],
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousReport, Error> {
        poll(gate)?;
        let frame = self.thermal_frame(frame, body_acceleration_m_s2)?;
        let prepared = self.prepare_viscous_heat(heat, gate)?;
        let tick = self.prepare_enthalpy_tick_with_frame(
            prepared.heated,
            prepared.external,
            prepared.viscous,
            0.0,
            prepared.absolute_heat,
            None,
            None,
            None,
            None,
            None,
            Some(frame),
            gate,
        )?;
        Ok(viscous::FreeSurfaceViscousReport {
            flow: self.accept_enthalpy_tick(tick).0,
            applied: prepared.applied,
        })
    }

    fn thermal_frame(
        &self,
        frame: MaterialRigidFrame3,
        body_acceleration_m_s2: [f64; 3],
    ) -> Result<ThermalRigidFrame, Error> {
        let (native, acceleration) =
            frame.lattice_forcing(&self.flow.scaling, body_acceleration_m_s2)?;
        if self.buoyancy.is_some()
            && native
                .angular_velocity
                .into_iter()
                .chain(native.angular_acceleration)
                .any(|value| value != 0.0)
        {
            return Err(Error::Invalid(
                "rotating thermal buoyancy requires density-dependent inertial coupling",
            ));
        }
        Ok(ThermalRigidFrame {
            native,
            acceleration,
        })
    }

    pub(super) fn prepare_frame_forces(
        &self,
        lattice: &mut FreeSurface3,
        resolved: &ResolvedFlowMaterials,
        frame: ThermalRigidFrame,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        let reference = self.flow.lattice.reference_force();
        frame.native.prepare_free_surface_step_with_acceleration(
            lattice,
            frame.acceleration,
            |key, lane, density| {
                if !density.is_finite() || density <= 0.0 {
                    return Err(Error::Invalid(
                        "frame forcing requires positive native density",
                    ));
                }
                let state =
                    resolved
                        .states
                        .get(&MassLocation3::Cell(key, lane))
                        .ok_or(Error::Invalid(
                            "frame forcing requires the original heated source",
                        ))?;
                // Remove exactly the existing owner's complete gravity load,
                // including its density departure. Surface stress stays as a
                // force density, never as another gravity or inertial term.
                let (gravity, departure) = match self.buoyancy {
                    Some(binding) => (
                        binding.force(&self.flow, state)?,
                        binding.density_departure(state),
                    ),
                    None => (reference, 0.0),
                };
                let total = match &resolved.forces {
                    Some(forces) => *forces.get(&(key, lane)).ok_or(Error::Invalid(
                        "thermal force requires every accepted wet cell",
                    ))?,
                    None => gravity,
                };
                let acceleration = core::array::from_fn(|axis| {
                    let effective_gravity =
                        frame.acceleration[axis] - frame.native.origin_acceleration[axis];
                    frame.acceleration[axis]
                        + departure * effective_gravity
                        + (total[axis] - gravity[axis]) / density
                });
                if !acceleration.into_iter().all(f64::is_finite) {
                    return Err(Error::Invalid("combined thermal frame force is not finite"));
                }
                // Frame inertia is subtracted once by the original midpoint
                // kernel. Its force reframe preserves physical strain, so the
                // already-resolved heated viscosity remains the collision tau.
                Ok(acceleration)
            },
            gate,
        )
    }
}
