//! Source-viscosity diffusion and physical heating before the original advance.

use super::super::density_projection::face_dual_mass;
use super::super::{DensityMomentumReconstructionControls, StaggeredBoundaryPressure};
use super::*;
use crate::conjugate::simple::FvBoundary;
pub use crate::conjugate::simple::viscous_heat::material_momentum::{
    MaterialViscosityControls, MaterialViscosityReport,
};
use crate::conjugate::simple::viscous_heat::material_momentum::{
    NativeMaterialViscosityInput, kick,
};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
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
    /// The pressure-profile entry additionally admits one opposed `Outlet`
    /// pair: normal traction with zero tangential trace, not legacy static
    /// pressure/zero-gradient openings. Their pressure comes only from the
    /// separate profile. Solid voxels remain stationary no-slip.
    pub boundaries: [FvBoundary; 6],
}

/// One privately prepared original-source molecular kick. Caloric consumers
/// decide how its physical heat enters their original material state owner.
pub(in crate::conjugate::flow) struct PreparedMaterialViscosity {
    pub momentum_kg_m_s: Vec<[f64; 3]>,
    pub report: MaterialViscosityReport,
    pub reconstruction: DensityMomentumReconstructionReport,
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
    pub fn advance_material_momentum_with_viscosity(
        &self,
        domain: &VoxelDomain,
        input: MaterialMomentumAdvanceInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: MaterialMomentumAdvanceControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousMaterialMomentumAdvance, ChtError> {
        self.advance_viscous_material(
            domain,
            [FlowFace::Wall; 6],
            None,
            None,
            input,
            viscosity,
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// Evolve actual material inertia through one opposed reservoir pair,
    /// retaining native molecular stress, heat, pressure work and mass routes.
    ///
    /// Exactly two opposed `Free` planes must correspond to `Outlet` in the
    /// viscosity input; remaining planes are `Wall` with no-slip or symmetry
    /// viscosity conditions. Each mouth prescribes normal traction and zero
    /// tangential velocity trace. Nonuniform mechanical pressure is supplied
    /// by the original staggered profile; it does not set thermodynamic rho.
    /// The entering material's mass, enthalpy and velocity remain explicit in
    /// `input`, including flow reversal. The supplied incoming velocity must
    /// have zero tangential components at either reservoir mouth.
    ///
    /// The native viscosity predictor uses actual boundary half-cell inertia
    /// with zero normal traction. The following force/pressure stage applies
    /// the reservoir load exactly once and retains its original work account.
    /// Only physical molecular stress loss enters h; pressure, projection,
    /// temporal, scatter and advective numerical losses do not become heat.
    ///
    /// This is the same first-order frozen-viscosity split on fixed geometry.
    /// Arbitrary opening configurations, nonzero tangential reservoir traces,
    /// an EOS and a complete Navier-Stokes accuracy claim are not supplied.
    ///
    /// # Errors
    /// Unsupported boundary combinations, conflicting pressure authority or
    /// any original source, work, transport or cancellation gate returns no
    /// candidate and leaves all inputs unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_viscosity_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        input: MaterialMomentumAdvanceInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: MaterialMomentumAdvanceControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousMaterialMomentumAdvance, ChtError> {
        poll(gate)?;
        let axis = Self::admit_material_viscosity_pressure_mouth(
            faces,
            viscosity.boundaries,
            input.boundary_velocity_m_s,
        )?;
        self.advance_viscous_material(
            domain,
            faces,
            Some(profile),
            Some(axis),
            input,
            viscosity,
            controls,
            viscosity_controls,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_viscous_material(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        pressure_mouth_axis: Option<usize>,
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
        let viscous = self.prepare_material_viscosity(
            domain,
            faces,
            pressure_mouth_axis,
            input.material.mass_kg,
            input.momentum_kg_m_s,
            viscosity,
            controls.mechanics.reconstruction,
            duration,
            viscosity_controls,
            gate,
            |_, _| Ok(()),
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
        let heated = MaterialMomentumAdvanceInput {
            material: MaterialAdvectionInput {
                specific_enthalpy_j_kg: &enthalpy,
                ..input.material
            },
            momentum_kg_m_s: &viscous.momentum_kg_m_s,
            ..input
        };
        let advance = match profile {
            Some(profile) => self.advance_material_momentum_with_pressure_profile(
                domain, faces, profile, heated, controls, gate,
            )?,
            None => self.advance_material_momentum(domain, faces, heated, controls, gate)?,
        };
        poll(gate)?;
        Ok(ViscousMaterialMomentumAdvance {
            advance,
            viscosity_reconstruction: viscous.reconstruction,
            viscosity: viscous.report,
            represented_molecular_heat_j: represented,
            absolute_heat_insertion_defect_j: defect,
        })
    }

    /// Share the original opposed-mouth admission with other caloric owners.
    pub(in crate::conjugate::flow) fn admit_material_viscosity_pressure_mouth(
        faces: [FlowFace; 6],
        boundaries: [FvBoundary; 6],
        boundary_velocity_m_s: [Option<[f64; 3]>; 6],
    ) -> Result<usize, ChtError> {
        let axis = (0..3)
            .find(|&a| faces[2 * a] == FlowFace::Free && faces[2 * a + 1] == FlowFace::Free)
            .ok_or_else(|| refuse("viscous pressure flow needs one opposed reservoir pair"))?;
        for f in 0..6 {
            if f / 2 == axis {
                if !matches!(boundaries[f], FvBoundary::Outlet)
                    || boundary_velocity_m_s[f]
                        .is_some_and(|u| (0..3).any(|a| a != axis && u[a] != 0.0))
                {
                    return Err(refuse(
                        "reservoir pressure belongs to the profile and tangential trace must be zero",
                    ));
                }
            } else if faces[f] != FlowFace::Wall
                || !matches!(
                    boundaries[f],
                    FvBoundary::Wall { .. } | FvBoundary::Symmetry
                )
            {
                return Err(refuse(
                    "viscous pressure flow admits one opposed pair and closed sides",
                ));
            }
        }
        Ok(axis)
    }

    /// Prepare molecular impulse and heat once on the original actual M/V
    /// inertia. The private source callback admits the SAME sampled transport
    /// state for a stricter caloric owner without reevaluating its source.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::conjugate::flow) fn prepare_material_viscosity(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        pressure_mouth_axis: Option<usize>,
        mass_kg: &[f64],
        momentum_kg_m_s: &[[f64; 3]],
        viscosity: MaterialMomentumViscosityInput<'_>,
        reconstruction_controls: DensityMomentumReconstructionControls,
        duration_s: f64,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        mut admit_source: impl FnMut(usize, &FluidTransportState) -> Result<(), ChtError>,
    ) -> Result<PreparedMaterialViscosity, ChtError> {
        poll(gate)?;
        let n = domain.cell_count();
        if n > viscosity_controls.max_cells
            || n > reconstruction_controls.max_cells
            || mass_kg.len() != n
            || viscosity.sources.len() != n
        {
            return Err(refuse("viscosity needs bounded original inventory"));
        }
        finite_positive("flow.viscosity_duration", duration_s)?;
        let reconstructed = self.reconstruct_material_momentum(
            domain,
            faces,
            mass_kg,
            momentum_kg_m_s,
            reconstruction_controls,
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
                    admit_source(c, &state)?;
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
                if pressure_mouth_axis == Some(a) {
                    for side in 0..2 {
                        let f = 2 * a + side;
                        if domain.neighbor(c, f).is_none() {
                            let (slot, _) = reconstructed.field.slot(domain, c, f);
                            masses[slot] =
                                face_dual_mass(domain, &reconstructed.density_kg_m3, c, f, None);
                        }
                    }
                }
            }
        }
        let viscous = kick(
            domain,
            NativeMaterialViscosityInput {
                field: &reconstructed.field,
                mass_kg,
                momentum_kg_m_s,
                viscosity_pa_s: &molecular,
                face_mass_kg: face_mass,
                boundaries: viscosity.boundaries,
                pressure_mouth_axis,
                duration_s,
            },
            viscosity_controls,
            gate,
        )?;
        Ok(PreparedMaterialViscosity {
            momentum_kg_m_s: viscous.momentum_kg_m_s,
            report: viscous.report,
            reconstruction: reconstructed.report,
        })
    }
}

#[cfg(test)]
mod tests;
