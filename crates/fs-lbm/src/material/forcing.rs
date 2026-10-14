//! Material-bound SI body forces for both existing thermal flow solvers.

use fs_material::fluid_transport::FluidTransport;
use crate::ThermalLbm;
use crate::d3q19::ThermalLbm3;
use crate::thermal::ThermalForcingError;
use super::{MaterialLbmError, MaterialScaling, MaterialSlab};

// Capture the slab's original units, so stepping cannot accidentally use a
// different material's density or a different plan's dx/dt conversion.
#[derive(Clone, Copy)]
pub(super) struct ForceUnits {
    reference_density: f64,
    timestep_s: f64,
    velocity_scale: f64,
}
impl ForceUnits {
    pub(super) fn from_scale<F: FluidTransport + ?Sized>(scale: &MaterialScaling<'_, F>) -> Self {
        Self {
            reference_density: scale.transport.density_kg_m3(),
            timestep_s: scale.timestep_s,
            velocity_scale: scale.velocity_scale,
        }
    }

    fn component(self, force_n_m3: f64) -> Result<f64, MaterialLbmError> {
        // Lattice force is momentum per lattice volume per lattice step,
        // NOT acceleration: F_lb = F_SI/rho_ref * dt²/dx.
        let value = ((force_n_m3 / self.reference_density) * self.timestep_s) / self.velocity_scale;
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
    /// Advance 3-D heat and flow with external body forces [N/m³] PLUS buoyancy.
    /// Uses the same material-bound conversion and refusal contract as the 2-D
    /// method. The field includes both wall planes and is ordered x-fastest.
    pub fn step_with_force_density(&mut self, forces_n_m3: &[[f64; 3]]) -> Result<(), MaterialLbmError> {
        let forces = self.force_units.field(forces_n_m3, self.lattice.gpop.len())?;
        self.lattice.step_with_external_forces(&forces)?;
        Ok(())
    }
}
