//! Prescribed frame and local Coriolis sources on one thermal mass/enthalpy tick.

use super::*;
use crate::d3q19::{Coriolis3, RigidFrame3};
use crate::material::free_surface::{MaterialCoriolis3, MaterialRigidFrame3};

#[derive(Clone, Copy)]
enum ThermalFrameSource {
    Rigid(RigidFrame3),
    Coriolis(Coriolis3),
}

#[derive(Clone, Copy)]
pub(super) struct ThermalRigidFrame {
    native: ThermalFrameSource,
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
    /// Full rigid-frame Boussinesq rotation, including angular acceleration at
    /// zero speed, refuses: omitted density/inertia couplings need a separate
    /// model. Use [`Self::step_with_coriolis_si`] for an explicit local model.
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

    /// Advance one thermal/VOF tick with the local Coriolis model in SI units.
    ///
    /// The effective gravity already includes the background centrifugal
    /// potential. With bound Boussinesq buoyancy, only this gravity is multiplied
    /// by `1 + (rho(T)-rho0)/rho0`; the Coriolis acceleration `-2 Omega cross u`
    /// uses the original reference inertia and is never multiplied by the
    /// thermal-density departure. Already-resolved nongravitational forces keep
    /// their force-density convention. The original material bounds, caloric
    /// mass, pressure datum and heated viscosity remain attached.
    ///
    /// This explicit local model admits arbitrary-axis rotation with bounded
    /// Boussinesq liquid. It does not admit full accelerating-vessel buoyancy:
    /// the refusal in [`Self::step_in_rigid_frame_si`] remains in force.
    /// The reduced momentum law is the rotating Boussinesq model described in
    /// the Oceananigans governing equations, section 3.1.1:
    /// <https://oceans.mit.edu/JohnMarshall/wp-content/uploads/2025/02/Oceananigans.pdf>.
    ///
    /// Heat is resolved before one native collision and its actual mass routes
    /// transport enthalpy once. Coriolis supplies no heat. One private candidate
    /// contains forcing, flow, thermal transport and endpoint admission; failure
    /// changes neither accepted state nor clock. The endpoint restores its
    /// original thermal force representation without an extra impulse.
    /// Supply rotation and gravity again on every interval. The original
    /// first-order thermal split and complete-boundary ownership limits remain;
    /// no global mechanical-energy or spherical-geometry claim is added.
    pub fn step_with_coriolis_si(
        &mut self,
        rotation: MaterialCoriolis3,
        effective_gravity_m_s2: [f64; 3],
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceEnthalpyReport, Error> {
        poll(gate)?;
        let frame = self.thermal_coriolis(rotation, effective_gravity_m_s2)?;
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

    /// The same local Coriolis tick with original molecular strain heating.
    ///
    /// Only the force-corrected `2 mu D:D` source is added to authored heat;
    /// neither Coriolis work, gravity work nor numerical kinetic loss becomes
    /// heat. The original source/strain admission and single thermal transaction
    /// are those of [`Self::step_in_rigid_frame_with_viscous_heating_si`].
    pub fn step_with_coriolis_and_viscous_heating_si(
        &mut self,
        rotation: MaterialCoriolis3,
        effective_gravity_m_s2: [f64; 3],
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousReport, Error> {
        poll(gate)?;
        let frame = self.thermal_coriolis(rotation, effective_gravity_m_s2)?;
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
            native: ThermalFrameSource::Rigid(native),
            acceleration,
        })
    }

    fn thermal_coriolis(
        &self,
        rotation: MaterialCoriolis3,
        effective_gravity_m_s2: [f64; 3],
    ) -> Result<ThermalRigidFrame, Error> {
        let (native, acceleration) =
            rotation.lattice_forcing(&self.flow.scaling, effective_gravity_m_s2)?;
        Ok(ThermalRigidFrame {
            native: ThermalFrameSource::Coriolis(native),
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
        let reference = lattice.material_gravity().ok_or(Error::Invalid(
            "thermal frame forcing requires the material-owned uniform gravity",
        ))?;
        let origin_acceleration = match frame.native {
            ThermalFrameSource::Rigid(native) => native.origin_acceleration,
            ThermalFrameSource::Coriolis(_) => [0.0; 3],
        };
        let acceleration_at = |key: u64, lane: usize, density: f64| {
            if !density.is_finite() || density <= 0.0 {
                return Err(Error::Invalid(
                    "frame forcing requires positive native density",
                ));
            }
            let state = resolved
                .states
                .get(&MassLocation3::Cell(key, lane))
                .ok_or(Error::Invalid(
                    "frame forcing requires the original heated source",
                ))?;
            // Remove exactly the existing owner's complete gravity load,
            // including its density departure. Surface stress stays as a
            // force density, never as another gravity or inertial term.
            let actual_gravity = reference.map(|acceleration| density * acceleration);
            let (gravity, departure) = match self.buoyancy {
                Some(binding) => (
                    binding.force(state, actual_gravity)?,
                    binding.density_departure(state),
                ),
                None => (actual_gravity, 0.0),
            };
            let total = match &resolved.forces {
                Some(forces) => *forces.get(&(key, lane)).ok_or(Error::Invalid(
                    "thermal force requires every accepted wet cell",
                ))?,
                None => gravity,
            };
            let acceleration = core::array::from_fn(|axis| {
                let effective_gravity = frame.acceleration[axis] - origin_acceleration[axis];
                frame.acceleration[axis]
                    + departure * effective_gravity
                    + (total[axis] - gravity[axis]) / density
            });
            if !acceleration.into_iter().all(f64::is_finite) {
                return Err(Error::Invalid("combined thermal frame force is not finite"));
            }
            // Only gravity receives the source-density departure.
            // Coriolis uses reference inertia and no thermal-density factor.
            // Reframing preserves physical strain and the heated viscosity.
            Ok(acceleration)
        };
        match frame.native {
            ThermalFrameSource::Rigid(native) => native.prepare_free_surface_step_with_acceleration(
                lattice,
                frame.acceleration,
                acceleration_at,
                gate,
            ),
            ThermalFrameSource::Coriolis(native) => native.prepare_free_surface_step_with_acceleration(
                lattice,
                frame.acceleration,
                acceleration_at,
                gate,
            ),
        }
    }
}
