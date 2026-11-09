//! Common turbulent species mixing from the original momentum closure.
//!
//! Gradient diffusion uses j_k = -rho (nu_t / Sc_t) grad(Y_k). The same
//! diffusivity for every constituent makes the sum of mass fluxes zero when
//! their fractions partition unity. See NASA/TM-2022-0008776, VULCAN-CFD
//! Theory Manual, equation (46):
//! https://vulcan-cfd.larc.nasa.gov/WebPage/Documentation/manual_theory.pdf.
//! The generic inventory entry is isothermal. The thermochemical entry below
//! delegates original constituent enthalpy and calorics to the material owner.
//! Neither entry infers species variance, reaction rates or molecular diffusion.

use super::FvFlow;
use crate::conjugate::{
    ChtError, MAX_MATERIAL_SPECIES, MaterialAdvectionStep, SpeciesDiffusionControls,
    SpeciesDiffusionStep, VoxelDomain, poll,
};
use fs_exec::CancelGate;

mod thermochemical;
pub use thermochemical::{TurbulentGasSpeciesAdvance, TurbulentGasSpeciesStep};

fn refuse(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.turbulent_species_diffusion",
        reason: reason.into(),
    }
}

impl FvFlow {
    /// Mix constituent masses using this flow's actual momentum eddy viscosity.
    ///
    /// `turbulent_schmidt` is an explicitly authored positive dimensionless
    /// number, independent of the thermal Prandtl number. Each fluid cell uses
    /// `D_t = eddy_viscosity_m2_s / turbulent_schmidt`, held fixed for this
    /// interval. Source-resolved molecular viscosity and the selected LVEL or
    /// transient SGS closure remain owned by the original flow solve. The heat
    /// wall-law secant is not a species wall model and is not used here.
    ///
    /// `carrier` supplies the accepted immutable cell masses. Supply constituent
    /// inventories from its `advect_species` endpoint or a previous mixing
    /// endpoint. The original graded grid, fluid mask, species order and the
    /// flow's physical epoch must match. Dimensions and the carrier's original
    /// mask are checked; matching the flow's mask, physical widths and epoch
    /// remains the caller's obligation.
    ///
    /// All components use one conservative face exchange. Solid and exterior
    /// faces are impermeable to this mixing operator; advective inlet/backflow
    /// composition stays with `advect_species`. No fraction is clipped or
    /// normalized. Zero eddy viscosity gives the exact unchanged inventories.
    /// Per-species conservation and carrier partition pass the original
    /// diffusion audit before publication, with no input or flow mutation.
    ///
    /// The storage budget includes the temporary cell diffusivity field as well
    /// as the common diffusion workspace. Coefficients and face exponentials
    /// are prepared once; sweeps perform O(N * faces * substeps) scalar work.
    /// Mesh and timestep refinement remain necessary. This is a frozen,
    /// isothermal gradient-diffusion component, not a complete reacting or
    /// nonisothermal mixture model or a qualification of the turbulence closure.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_turbulent_species<const N: usize>(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species_mass_kg: &[[f64; N]],
        turbulent_schmidt: f64,
        duration_s: f64,
        mut controls: SpeciesDiffusionControls<N>,
        gate: &CancelGate,
    ) -> Result<SpeciesDiffusionStep<N>, ChtError> {
        poll(gate)?;
        if !turbulent_schmidt.is_finite() || turbulent_schmidt <= 0.0 {
            return Err(refuse(
                "turbulent Schmidt number must be finite and positive",
            ));
        }
        let cells = domain.cell_count();
        if N == 0
            || N > MAX_MATERIAL_SPECIES
            || self.field.dims() != domain.dims()
            || self.eddy_viscosity_m2_s.len() != cells
            || species_mass_kg.len() != cells
            || carrier.mass_kg.len() != cells
        {
            return Err(refuse("flow, carrier and species dimensions must match"));
        }
        let bytes = cells
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or_else(|| refuse("diffusivity storage extent overflows"))?;
        controls.maximum_storage_bytes = controls
            .maximum_storage_bytes
            .checked_sub(bytes)
            .ok_or_else(|| refuse("turbulent species storage budget exhausted"))?;
        let mut diffusivity = Vec::new();
        diffusivity
            .try_reserve_exact(cells)
            .map_err(|_| refuse("unable to allocate turbulent species diffusivity"))?;
        for (cell, &eddy) in self.eddy_viscosity_m2_s.iter().enumerate() {
            if cell % 256 == 0 {
                poll(gate)?;
            }
            if !eddy.is_finite() || eddy < 0.0 || (!domain.is_fluid(cell) && eddy != 0.0) {
                return Err(refuse(
                    "eddy viscosity must be nonnegative, finite and zero in solids",
                ));
            }
            let value = eddy / turbulent_schmidt;
            if !value.is_finite() || (eddy > 0.0 && value == 0.0) {
                return Err(refuse(
                    "positive turbulent species diffusivity is unrepresentable",
                ));
            }
            diffusivity.push(value);
        }
        carrier.diffuse_species(
            domain,
            species_mass_kg,
            &diffusivity,
            duration_s,
            controls,
            gate,
        )
    }
}

#[cfg(test)]
mod tests;
