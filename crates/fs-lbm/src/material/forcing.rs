//! Material-bound SI momentum and heat sources for the existing thermal solvers.

mod heating;
pub use heating::transport;
pub mod electrothermal;
mod thermal_walls;
pub use thermal_walls::{MaterialThermalWall, MaterialWallExchange};
mod moving_walls;
pub use moving_walls::MaterialPlateExchange;

use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use crate::ThermalLbm;
use crate::d3q19::ThermalLbm3;
use crate::thermal::ThermalForcingError;
use super::{MaterialLbmError, MaterialScaling, MaterialSlab};

// Capture the slab's original units, so stepping cannot accidentally use a
// different material's density or a different plan's dx/dt conversion.
#[derive(Clone, Copy, PartialEq)]
pub(super) struct ForceUnits {
    transport: FluidTransportState,
    cell_size_m: f64,
    timestep_s: f64,
    velocity_scale: f64,
}
impl ForceUnits {
    pub(super) fn from_scale<F: FluidTransport + ?Sized>(scale: &MaterialScaling<'_, F>) -> Self {
        Self {
            transport: scale.transport,
            cell_size_m: scale.cell_size_m,
            timestep_s: scale.timestep_s,
            velocity_scale: scale.velocity_scale,
        }
    }

    fn angular_velocity(self, angular_velocity_rad_s: f64) -> Result<f64, MaterialLbmError> {
        let converted = angular_velocity_rad_s * self.timestep_s;
        if !angular_velocity_rad_s.is_finite() || !converted.is_finite()
            || (angular_velocity_rad_s != 0.0 && converted == 0.0)
        {
            return Err(MaterialLbmError::Invalid(
                "physical angular velocity is not representable in lattice units",
            ));
        }
        Ok(converted)
    }

    fn component(self, force_n_m3: f64) -> Result<f64, MaterialLbmError> {
        // Lattice force is momentum per lattice volume per lattice step,
        // NOT acceleration: F_lb = F_SI/rho_ref * dt²/dx.
        let value = ((force_n_m3 / self.transport.density_kg_m3()) * self.timestep_s) / self.velocity_scale;
        if !force_n_m3.is_finite() || !value.is_finite() || (force_n_m3 != 0.0 && value == 0.0) {
            return Err(MaterialLbmError::Invalid("physical body-force density is not representable in lattice units"));
        }
        Ok(value)
    }

    fn field<const D: usize>(self, forces: &[[f64; D]], cells: usize) -> Result<Vec<[f64; D]>, MaterialLbmError> {
        if forces.len() != cells {
            return Err(ThermalForcingError::Length { expected: cells, actual: forces.len() }.into());
        }
        forces.iter().map(|force| {
            let mut converted = [0.0; D];
            for (target, &source) in converted.iter_mut().zip(force) {
                *target = self.component(source)?;
            }
            Ok(converted)
        }).collect()
    }
}

impl<F: FluidTransport + ?Sized> MaterialScaling<'_, F> {
    /// A signed body-force density component [N/m³] to Guo lattice force.
    /// Unlike acceleration, this already includes the physical mass density.
    /// Use for volumetric mechanical or electromagnetic loads supplied by
    /// another model. This conversion does not supply that model's force law.
    ///
    /// # Errors
    /// Refuses nonfinite input, overflow and disappearance of a nonzero force.
    pub fn force_density_lattice(&self, force_n_m3: f64) -> Result<f64, MaterialLbmError> {
        ForceUnits::from_scale(self).component(force_n_m3)
    }
}

impl MaterialSlab<ThermalLbm> {
    /// Configure persistent local Coriolis rotation about +z [rad/s].
    ///
    /// Uses this slab's retained timestep: omega_lattice = omega_SI * dt.
    /// All ordinary, prescribed-heat, material transport and rheology steps
    /// then compose rotation with their original selected force owner.
    /// Boussinesq density departures affect effective gravity only; rotation
    /// supplies no heat. Gravity must already include centrifugal potential.
    /// Initial velocities use the moving axes, as in the native setter.
    ///
    /// Zero disables rotation on the next step. Invalid/unrepresentable input
    /// leaves the setting, physical populations and thermal history unchanged.
    pub fn set_coriolis_si(&mut self, angular_velocity_rad_s: f64) -> Result<(), MaterialLbmError> {
        let angular_velocity = self.force_units.angular_velocity(angular_velocity_rad_s)?;
        self.lattice.set_coriolis(crate::core2::frame::Coriolis2 { angular_velocity })?;
        Ok(())
    }

    /// Advance 2-D heat and flow with external body forces [N/m³] PLUS buoyancy.
    /// Supply one vector per grid cell, in its existing linear order, including
    /// wall cells (whose forces are ignored). No force persists implicitly:
    /// submit the external field for each step; never resubmit total `fext`.
    /// The slab retains its original material density and dx/dt conversion.
    ///
    /// # Errors
    /// Refuses invalid shape, nonfinite/unrepresentable forces or overflow of
    /// their sum with buoyancy before mutating the solver. This does not turn
    /// an unstable or otherwise invalid underlying solver state into a valid one.
    pub fn step_with_force_density(&mut self, forces_n_m3: &[[f64; 2]]) -> Result<(), MaterialLbmError> {
        let forces = self.force_units.field(forces_n_m3, self.lattice.gpop.len())?;
        self.lattice.step_with_external_forces(&forces)?;
        Ok(())
    }
}

impl MaterialSlab<ThermalLbm3> {
    /// Configure persistent arbitrary-axis local Coriolis rotation [rad/s].
    ///
    /// Shares the two-dimensional SI/force/heat/rollback contract. Components
    /// are resolved in the slab axes and converted with its retained timestep.
    /// Use the local effective-gravity model, not full vessel kinematics.
    pub fn set_coriolis_si(
        &mut self,
        angular_velocity_rad_s: [f64; 3],
    ) -> Result<(), MaterialLbmError> {
        let mut angular_velocity = [0.0; 3];
        for axis in 0..3 {
            angular_velocity[axis] = self.force_units.angular_velocity(angular_velocity_rad_s[axis])?;
        }
        self.lattice.set_coriolis(crate::d3q19::Coriolis3 { angular_velocity })?;
        Ok(())
    }

    /// Advance 3-D heat and flow with external body forces [N/m³] PLUS buoyancy.
    /// Uses the same material-bound conversion and refusal contract as the 2-D
    /// method. The field includes both wall planes and is ordered x-fastest.
    pub fn step_with_force_density(&mut self, forces_n_m3: &[[f64; 3]]) -> Result<(), MaterialLbmError> {
        let forces = self.force_units.field(forces_n_m3, self.lattice.gpop.len())?;
        self.lattice.step_with_external_forces(&forces)?;
        Ok(())
    }
}
