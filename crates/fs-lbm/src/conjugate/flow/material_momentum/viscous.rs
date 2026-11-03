//! Source-viscosity diffusion and physical heating before the original advance.

use super::super::density_projection::face_dual_mass;
use super::*;
use crate::conjugate::simple::FvBoundary;
pub use crate::conjugate::simple::viscous_heat::material_momentum::{
    MaterialViscosityControls, MaterialViscosityReport,
};
use crate::conjugate::simple::viscous_heat::material_momentum::{
    NativeMaterialViscosityInput, kick,
};
use fs_material::fluid_transport::FluidTransport;
use fs_math::dd::Dd;

/// Borrow original transport evaluations; no viscosity or EOS is inferred.
#[derive(Clone, Copy)]
pub struct MaterialMomentumViscosityInput<'a> {
    /// One source for each fluid cell, None for each fixed solid cell. This
    /// delegates to the original FluidTransport owner (including resolved
    /// liquid/gas sources); a raw FluidTransportState retains raw authority.
    /// Only dynamic viscosity is consumed. The caller evaluates it at the
    /// current material state; actual M/V, not source reference rho, is inertia.
    pub sources: &'a [Option<&'a dyn FluidTransport>],
    /// Closed no-slip walls, optionally moving tangentially, or symmetry.
    /// Solid voxels remain stationary no-slip. Openings are not admitted.
    pub boundaries: [FvBoundary; 6],
}

/// One atomic molecular-viscosity/force/pressure/material-advection candidate.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscousMaterialMomentumAdvance {
    pub advance: MaterialMomentumAdvance,
    /// Original reconstruction sampled by viscosity. This sampling does not
    /// separately apply its held-face reactions to cell momentum; the normal
    /// force/pressure advance owns that constraint once on its current input.
    pub viscosity_reconstruction: DensityMomentumReconstructionReport,
    pub viscosity: MaterialViscosityReport,
    /// Actual sum M*(h_heated-h_original) before material advection.
    pub represented_molecular_heat_j: f64,
    /// Non-cancelling local heat-insertion roundoff, included in the viscosity
    /// energy gate. Numerical time/scatter/projection losses never enter h.
    pub absolute_heat_insertion_defect_j: f64,
}

impl FlowField {
    /// Extend the original material-momentum advance with physical molecular
    /// stress and its reciprocal heat, reusing the existing MAC operator.
    ///
    /// Original source evaluations supply frozen dynamic viscosity. One native
    /// backward-Euler diffusion predictor uses the same actual M/V face inertia
    /// as pressure and the complete symmetric stress, including harmonic
    /// material-jump transmission, moving walls and stationary solid obstacles.
    /// Its native coupled rows are solved without convection or pressure.
    /// Half-cell mass shares transfer the viscous impulse to the ORIGINAL cell
    /// momenta. Only duration times physical molecular dissipation is inserted
    /// into h; both backward-Euler and scatter losses remain numerical.
    ///
    /// The existing advance then applies force/pressure once and transports
    /// the resulting M, M*h and P on its accepted mass routes. Its retained
    /// flow is still the transport-stage field, not endpoint P/M. Subsequent
    /// calls use current inventory and freshly evaluated material sources.
    ///
    /// This is a first-order split on closed fixed geometry. Source mu is
    /// frozen over the step; viscosity feedback, time/space convergence and
    /// consistent evaluation of current material state remain caller-owned.
    /// There is no added EOS, phase mobility, bulk-viscosity law, turbulence,
    /// pressure-opening work model or complete Navier-Stokes accuracy claim.
    ///
    /// # Errors
    /// Native admission, source evaluation, solve, physical-balance, caloric
    /// representation, transport or cancellation refusal publishes no state.
    #[allow(clippy::too_many_lines)] // one atomic composition of original owners
    pub fn advance_material_momentum_with_viscosity(
        &self,
        domain: &VoxelDomain,
        input: MaterialMomentumAdvanceInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: MaterialMomentumAdvanceControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousMaterialMomentumAdvance, ChtError> {
        poll(gate)?;
        let n = domain.cell_count();
        if n > viscosity_controls.max_cells
            || n > controls.mechanics.reconstruction.max_cells
            || n > controls.mechanics.forcing.max_cells
            || n > controls.material.max_cells
            || input.material.mass_kg.len() != n
            || input.material.specific_enthalpy_j_kg.len() != n
            || viscosity.sources.len() != n
            || controls.mechanics.projection != controls.material.projection
        {
            return Err(refuse(
                "viscous advance needs bounded original inventory and matching clocks",
            ));
        }
        let duration = controls.mechanics.projection.duration_s;
        finite_positive("flow.viscosity_duration", duration)?;
        let reconstructed = self.reconstruct_material_momentum(
            domain,
            [FlowFace::Wall; 6],
            input.material.mass_kg,
            input.momentum_kg_m_s,
            controls.mechanics.reconstruction,
            gate,
        )?;
        let mut molecular = vec![0.0; n];
        let mut first = None;
        for c in 0..n {
            poll_every(gate, c)?;
            match (domain.is_fluid(c), viscosity.sources[c]) {
                (true, Some(source)) => {
                    let state =
                        source
                            .transport_state()
                            .map_err(|error| ChtError::InvalidInput {
                                field: "flow.viscosity_source",
                                reason: format!("cell {c}: {error}"),
                            })?;
                    molecular[c] = state.dynamic_viscosity_pa_s();
                    first.get_or_insert(molecular[c]);
                }
                (false, None) => {}
                _ => {
                    return Err(refuse(
                        "viscosity sources must cover exactly the mobile cells",
                    ));
                }
            }
        }
        let reference = first.ok_or_else(|| refuse("viscosity requires mobile material"))?;
        // The old native coefficient vector has a positive slot for every
        // voxel. Solid slots are never sampled or given any fluid dissipation.
        for (c, mu) in molecular.iter_mut().enumerate() {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                *mu = reference;
            }
        }
        let mut face_mass: [Vec<f64>; 3] =
            std::array::from_fn(|a| vec![0.0; reconstructed.field.face_arrays()[a].len()]);
        for c in 0..n {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                continue;
            }
            for (a, masses) in face_mass.iter_mut().enumerate() {
                let face = 2 * a + 1;
                if let Some(neighbor) = domain.neighbor(c, face).filter(|&p| domain.is_fluid(p)) {
                    let (slot, _) = reconstructed.field.slot(domain, c, face);
                    masses[slot] = face_dual_mass(
                        domain,
                        &reconstructed.density_kg_m3,
                        c,
                        face,
                        Some(neighbor),
                    );
                }
            }
        }
        let viscous = kick(
            domain,
            NativeMaterialViscosityInput {
                field: &reconstructed.field,
                mass_kg: input.material.mass_kg,
                momentum_kg_m_s: input.momentum_kg_m_s,
                viscosity_pa_s: &molecular,
                face_mass_kg: face_mass,
                boundaries: viscosity.boundaries,
                duration_s: duration,
            },
            viscosity_controls,
            gate,
        )?;
        let mut enthalpy = input.material.specific_enthalpy_j_kg.to_vec();
        let mut represented = Dd::ZERO;
        let mut insertion_defect = Dd::ZERO;
        for c in 0..n {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                continue;
            }
            finite("flow.viscosity_original_enthalpy", enthalpy[c])?;
            let heat = viscous.report.cell_molecular_heat_j[c];
            if heat != 0.0 {
                enthalpy[c] += heat / input.material.mass_kg[c];
            }
            finite("flow.viscosity_heated_enthalpy", enthalpy[c])?;
            let actual = Dd::from_f64(input.material.mass_kg[c])
                * (Dd::from_f64(enthalpy[c])
                    - Dd::from_f64(input.material.specific_enthalpy_j_kg[c]));
            represented = represented + actual;
            insertion_defect =
                insertion_defect + Dd::from_f64((actual - Dd::from_f64(heat)).to_f64().abs());
        }
        let defect = insertion_defect.to_f64();
        finite("flow.viscosity_heat_insertion_defect", defect)?;
        if (insertion_defect + Dd::from_f64(viscous.report.absolute_energy_defect_j)).to_f64()
            > viscosity_controls.maximum_energy_defect_j
        {
            return Err(refuse(
                "actual molecular heating exceeds the viscosity energy budget",
            ));
        }
        let represented = represented.to_f64();
        finite("flow.viscosity_represented_heat", represented)?;
        let advance = self.advance_material_momentum(
            domain,
            [FlowFace::Wall; 6],
            MaterialMomentumAdvanceInput {
                material: MaterialAdvectionInput {
                    specific_enthalpy_j_kg: &enthalpy,
                    ..input.material
                },
                momentum_kg_m_s: &viscous.momentum_kg_m_s,
                ..input
            },
            controls,
            gate,
        )?;
        poll(gate)?;
        Ok(ViscousMaterialMomentumAdvance {
            advance,
            viscosity_reconstruction: reconstructed.report,
            viscosity: viscous.report,
            represented_molecular_heat_j: represented,
            absolute_heat_insertion_defect_j: defect,
        })
    }
}

#[cfg(test)]
mod tests;
