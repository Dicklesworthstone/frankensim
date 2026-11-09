//! Common turbulent species mixing from the original momentum closure.
//!
//! Gradient diffusion uses j_k = -rho (nu_t / Sc_t) grad(Y_k). The same
//! diffusivity for every constituent makes the sum of mass fluxes zero when
//! their fractions partition unity. See NASA/TM-2022-0008776, VULCAN-CFD
//! Theory Manual, equation (46):
//! https://vulcan-cfd.larc.nasa.gov/WebPage/Documentation/manual_theory.pdf.
//! The inventory entries are isothermal. The combined entry consumes explicit
//! common molecular diffusivity. The thermochemical entry below
//! delegates original constituent enthalpy and calorics to the material owner.
//! No entry infers species variance, reaction rates or molecular coefficients.

use super::FvFlow;
use crate::conjugate::{
    AdaptiveSpeciesDiffusionControls, AdaptiveSpeciesDiffusionStep, ChtError, MAX_MATERIAL_SPECIES,
    MaterialAdvectionStep, SpeciesDiffusionControls, SpeciesDiffusionStep, VoxelDomain, poll,
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

// Cellwise common-D constitutive combination; geometry and solid-mask admission
// remain with the native field owner. The result retains physical m²/s units.
fn combined_species_diffusivity(
    molecular: f64,
    eddy: f64,
    turbulent_schmidt: f64,
) -> Result<f64, &'static str> {
    if !molecular.is_finite() || molecular < 0.0 {
        return Err("molecular diffusivity must be finite and nonnegative");
    }
    if !eddy.is_finite() || eddy < 0.0 {
        return Err("eddy viscosity must be finite and nonnegative");
    }
    if !turbulent_schmidt.is_finite() || turbulent_schmidt <= 0.0 {
        return Err("turbulent Schmidt number must be finite and positive");
    }
    let turbulent = eddy / turbulent_schmidt;
    if !turbulent.is_finite() || (eddy > 0.0 && turbulent == 0.0) {
        return Err("positive turbulent species diffusivity is unrepresentable");
    }
    let combined = molecular + turbulent;
    if !combined.is_finite() {
        return Err("combined species diffusivity is unrepresentable");
    }
    Ok(combined)
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
        let diffusivity = self.species_mixing_diffusivity(
            domain,
            carrier,
            species_mass_kg,
            None,
            turbulent_schmidt,
            &mut controls.maximum_storage_bytes,
            gate,
        )?;
        carrier.diffuse_species(
            domain,
            species_mass_kg,
            &diffusivity,
            duration_s,
            controls,
            gate,
        )
    }

    /// Resolve native turbulent species mixing to explicit temporal tolerances.
    ///
    /// Uses the same momentum eddy viscosity and independently authored Schmidt
    /// number as [`Self::diffuse_turbulent_species`]. The original adaptive
    /// conservative solver chooses substeps and reports its per-species
    /// mass-fraction estimate. Every rejected resolution counts against the
    /// work budget. The temporary D field and live coarse comparison both count
    /// against storage. No flow, carrier or species input is modified.
    ///
    /// This controls frozen-coefficient temporal splitting only. Mesh, closure,
    /// physical epoch and original geometry obligations are unchanged. For
    /// nonzero molecular diffusion in a laminar region, use
    /// [`Self::diffuse_molecular_and_turbulent_species_adaptive`].
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_turbulent_species_adaptive<const N: usize>(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species_mass_kg: &[[f64; N]],
        turbulent_schmidt: f64,
        duration_s: f64,
        mut controls: AdaptiveSpeciesDiffusionControls<N>,
        gate: &CancelGate,
    ) -> Result<AdaptiveSpeciesDiffusionStep<N>, ChtError> {
        let diffusivity = self.species_mixing_diffusivity(
            domain,
            carrier,
            species_mass_kg,
            None,
            turbulent_schmidt,
            &mut controls.solver.maximum_storage_bytes,
            gate,
        )?;
        carrier.diffuse_species_adaptive(
            domain,
            species_mass_kg,
            &diffusivity,
            duration_s,
            controls,
            gate,
        )
    }

    /// Resolve additive molecular and turbulent mixing on the native carrier.
    ///
    /// Every fluid cell supplies one explicit nonnegative molecular D [m²/s].
    /// The common mass-frame law is
    /// `j_k = -rho (D_molecular + nu_t / Sc_t) grad(Y_k)`, with nu_t from this
    /// flow's original MOMENTUM closure, not its thermal wall model. Molecular
    /// diffusion therefore survives where the eddy viscosity vanishes. The
    /// total cell coefficient enters the original two half-cell resistances;
    /// independently harmonizing molecular and eddy face coefficients would
    /// describe a different spatial interface law.
    ///
    /// This admits binary Fick diffusion, or a molecular model whose COMMON
    /// diffusivity across all N species is independently justified by the
    /// caller. It does not invent unequal-species molecular coefficients,
    /// Maxwell–Stefan transport, thermodynamic factors or partial enthalpies.
    /// Coefficients are frozen for this interval; the same temporal-estimate,
    /// geometry, epoch, conservation and cancellation boundaries as
    /// [`Self::diffuse_turbulent_species_adaptive`] apply. No physical diffusion
    /// is inferred from molecular viscosity or the thermal Prandtl number.
    ///
    /// Solids require exact zero molecular and eddy coefficients. Zero total
    /// diffusivity blocks an interface. Negative, nonfinite, overflowing or
    /// positively underflowed coefficients refuse without floors or clipping.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_molecular_and_turbulent_species_adaptive<const N: usize>(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species_mass_kg: &[[f64; N]],
        molecular_diffusivity_m2_s: &[f64],
        turbulent_schmidt: f64,
        duration_s: f64,
        mut controls: AdaptiveSpeciesDiffusionControls<N>,
        gate: &CancelGate,
    ) -> Result<AdaptiveSpeciesDiffusionStep<N>, ChtError> {
        let diffusivity = self.species_mixing_diffusivity(
            domain,
            carrier,
            species_mass_kg,
            Some(molecular_diffusivity_m2_s),
            turbulent_schmidt,
            &mut controls.solver.maximum_storage_bytes,
            gate,
        )?;
        carrier.diffuse_species_adaptive(
            domain,
            species_mass_kg,
            &diffusivity,
            duration_s,
            controls,
            gate,
        )
    }

    // One owned coefficient field serves fixed and adaptive native mixing.
    // The molecular input stays borrowed, while the combined copy is budgeted.
    #[allow(clippy::too_many_arguments)]
    fn species_mixing_diffusivity<const N: usize>(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species_mass_kg: &[[f64; N]],
        molecular_diffusivity_m2_s: Option<&[f64]>,
        turbulent_schmidt: f64,
        remaining_storage_bytes: &mut usize,
        gate: &CancelGate,
    ) -> Result<Vec<f64>, ChtError> {
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
            || molecular_diffusivity_m2_s.is_some_and(|d| d.len() != cells)
        {
            return Err(refuse("flow, carrier and species dimensions must match"));
        }
        let bytes = cells
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or_else(|| refuse("diffusivity storage extent overflows"))?;
        *remaining_storage_bytes = remaining_storage_bytes
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
            let molecular = molecular_diffusivity_m2_s.map_or(0.0, |d| d[cell]);
            if !domain.is_fluid(cell) && (molecular != 0.0 || eddy != 0.0) {
                return Err(refuse(
                    "molecular diffusivity and eddy viscosity must be zero in solids",
                ));
            }
            diffusivity.push(
                combined_species_diffusivity(molecular, eddy, turbulent_schmidt).map_err(refuse)?,
            );
        }
        Ok(diffusivity)
    }
}

#[cfg(test)]
mod tests;
