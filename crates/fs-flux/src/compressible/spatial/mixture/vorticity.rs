//! Vorticity production from the original evolving multicomponent field.
//!
//! The shared spatial budget retains stretching, compressible dilation and
//! baroclinicity. Every cell's pressure is recovered from its own species
//! densities and referenced total energy, including temperature-dependent heat
//! capacities; no frozen-composition or constant-gamma surrogate is introduced.
//! Composition gradients influence the pressure/density misalignment through
//! that EOS. They are not an additional, separately added momentum source.
//!
//! This Estimated smooth-field observation has the same physical derivative,
//! boundary and resolution restrictions as the frozen-gas budget. It excludes
//! viscous, imposed-force and modeled stress contributions and never advances
//! accepted state, time or loads.

use super::{Error, MixtureEuler3d, Primitive, allocated, cell, stepping};
pub use crate::compressible::spatial::vorticity::{VorticityBudget3d, VorticityWorkspace3d};
use crate::compressible::spatial::vorticity::{admit_axes, budget};

impl MixtureEuler3d {
    /// Numerical bytes for source primitives and one vorticity budget per cell.
    pub fn required_vorticity_workspace_bytes(&self) -> Result<usize, Error> {
        self.grid
            .len()
            .checked_mul(size_of::<Primitive>() + size_of::<VorticityBudget3d>())
            .ok_or(Error::Invalid("mixture vorticity workspace size overflow"))
    }

    /// Allocate reusable observation scratch within an explicit byte cap.
    ///
    /// Resolved axes need at least three cells; singleton periodic axes denote
    /// homogeneous extrusion. The workspace is independent of species count.
    pub fn vorticity_workspace(&self, maximum_bytes: usize) -> Result<VorticityWorkspace3d, Error> {
        admit_axes::<Error>(&self.grid, self.boundaries)?;
        if self.required_vorticity_workspace_bytes()? > maximum_bytes {
            return Err(Error::Invalid(
                "mixture vorticity workspace exceeds its byte budget",
            ));
        }
        let work = VorticityWorkspace3d {
            grid: self.grid,
            primitive: allocated(self.grid.len(), Primitive::default())?,
            output: allocated(self.grid.len(), VorticityBudget3d::default())?,
        };
        if work.workspace_bytes() > maximum_bytes {
            return Err(Error::Invalid(
                "mixture vorticity workspace capacity exceeds its byte budget",
            ));
        }
        Ok(work)
    }

    /// Observe accepted cells in x-fast order through their original mixture EOS.
    ///
    /// Decode each cell once, allocating no fields. Source recovery, gradient
    /// traversal and final publication poll cancellation; failed evaluations
    /// expose no partial results and leave this owner unchanged.
    pub fn vorticity<'a>(
        &self,
        work: &'a mut VorticityWorkspace3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&'a [VorticityBudget3d], Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        admit_axes::<Error>(&self.grid, self.boundaries)?;
        if work.grid != self.grid {
            return Err(Error::Invalid(
                "mixture vorticity workspace and physical grid differ",
            ));
        }
        let species = self.gas.species_count();
        for (index, (primitive, state)) in work
            .primitive
            .iter_mut()
            .zip(self.states.chunks_exact(species + 4))
            .enumerate()
        {
            if index % 64 == 0 {
                stepping::poll::<Error>(&mut cancelled)?;
            }
            *primitive = cell(state, species).primitive(&self.gas, &mut cancelled)?;
        }
        for (index, output) in work.output.iter_mut().enumerate() {
            if index % 64 == 0 {
                stepping::poll::<Error>(&mut cancelled)?;
            }
            *output = budget::<Error>(&self.grid, self.boundaries, &work.primitive, index)?;
        }
        stepping::poll::<Error>(&mut cancelled)?;
        Ok(&work.output)
    }
}
