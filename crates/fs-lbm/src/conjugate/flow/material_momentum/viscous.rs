//! Source-viscosity diffusion and physical heating before the original advance.
//! Existing entries retain 2 mu D; the bulk opt-in admits independent zeta.
//! The SGS opt-in adds deviatoric modeled stress with a separate energy ledger.

use super::super::density_projection::face_dual_mass;
use super::super::{
    DensityMomentumReconstruction, DensityMomentumReconstructionControls, StaggeredBoundaryPressure,
};
use super::*;
use crate::conjugate::FvSgsConfig;
use crate::conjugate::simple::FvBoundary;
pub use crate::conjugate::simple::viscous_heat::material_momentum::{
    MaterialViscosityControls, MaterialViscosityReport,
};
use crate::conjugate::simple::viscous_heat::material_momentum::{
    NativeMaterialRheologyInput, NativeMaterialRheologyStep, NativeMaterialViscosityInput, kick,
    kick_with_sgs, kick_with_source, kick_with_source_and_sgs,
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
    /// Original accepted cell SGS kinematic viscosity [m^2/s], transferred
    /// from the native stress solve for thermal transport. Empty means the
    /// original or zero-model path; active quiescent SGS retains zero slots.
    #[cfg(feature = "implicit-enthalpy-coupling")]
    pub cell_eddy_viscosity_m2_s: Vec<f64>,
    pub report: MaterialViscosityReport,
    pub reconstruction: DensityMomentumReconstructionReport,
}

/// One endpoint-rheology solve on original native material inertia.
pub(in crate::conjugate::flow) struct PreparedMaterialRheology {
    pub viscosity: PreparedMaterialViscosity,
    /// Native physical rates at the accepted molecular endpoint, before
    /// scattering to cells or applying pressure/material transport.
    pub endpoint_shear_rate_per_s: Vec<f64>,
    /// Includes initial and all post-sweep source resolutions.
    pub source_updates: usize,
}

/// One atomic viscosity/force/pressure/material-advection candidate.
/// Modeled SGS transfer, when selected, remains distinct from molecular heat.
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
    /// energy gate. Modeled SGS transfer and numerical time/scatter/projection
    /// losses never enter h.
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
    /// This original entry retains `2 mu D`; the bulk-viscosity entry below
    /// admits an independently supplied volumetric response. There is no added
    /// EOS, phase mobility, turbulence,
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
            None,
            None,
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// Add an explicit resolved-velocity SGS closure to original-source
    /// molecular stress on the same actual-density material-momentum step.
    ///
    /// The existing reconstruction supplies actual cell `rho = M/V`, native
    /// face inertia and the velocity gradient for the selected Smagorinsky,
    /// WALE or Vreman model. Each nonlinear sweep refreshes the closure on
    /// its current native velocity; cell and reciprocal shear-edge viscosity
    /// caps use the original source mu and actual adjacent material density.
    /// Source reference density never replaces material inertia or converts
    /// the modeled kinematic viscosity into dynamic stress.
    ///
    /// Molecular stress retains the original `2 mu D` response. The added
    /// stress is `2 rho nu_t dev(D)`; it supplies no bulk viscosity. Only
    /// molecular dissipation enters specific enthalpy. Modeled transfer is
    /// retained separately in `viscosity.cell_sgs_transfer_j` and
    /// `viscosity.sgs_transfer_j`, alongside native temporal/scatter losses.
    /// The original force, pressure and same-mass-route transport then run
    /// once on this candidate.
    ///
    /// This is an exploratory, first-order split on closed fixed geometry
    /// with walls or symmetry and frozen source mu. The original filter
    /// policy and applicability limits remain explicit. A zero coefficient
    /// is the exact existing molecular-only path; acceptance renders refuse.
    ///
    /// # Errors
    /// Invalid model, unsupported boundaries, source, cap, nonlinear solve,
    /// physical-energy, caloric or cancellation failure publishes no state.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_viscosity_and_sgs(
        &self,
        domain: &VoxelDomain,
        input: MaterialMomentumAdvanceInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        sgs: FvSgsConfig,
        controls: MaterialMomentumAdvanceControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousMaterialMomentumAdvance, ChtError> {
        poll(gate)?;
        sgs.validate()?;
        self.advance_viscous_material(
            domain,
            [FlowFace::Wall; 6],
            None,
            None,
            input,
            viscosity,
            None,
            Some(sgs),
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// Advance with independently supplied bulk viscosity [Pa s] in addition
    /// to the original transport source's dynamic shear viscosity.
    ///
    /// Each mobile cell uses `2 mu dev(D) + zeta div(u) I` and deposits only
    /// `2 mu dev(D):dev(D) + zeta div(u)^2` as physical heat. Bulk viscosity
    /// must be finite and nonnegative; choosing zero is explicit. The slice
    /// has one slot per voxel; fixed-solid slots are ignored. The caller
    /// resolves zeta at the same material state as the original shear source.
    /// No bulk law is inferred from shear viscosity or thermodynamic pressure.
    /// Both coefficients remain frozen during the original first-order split.
    ///
    /// # Errors
    /// Invalid coefficients or any original source, solve, work, heat-insertion
    /// or cancellation refusal returns no candidate and changes no input.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_bulk_viscosity(
        &self,
        domain: &VoxelDomain,
        input: MaterialMomentumAdvanceInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: &[f64],
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
            Some(bulk_viscosity_pa_s),
            None,
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// Combine an independently supplied bulk viscosity with explicit SGS on
    /// the original closed material-momentum update.
    ///
    /// Molecular stress is `2 mu dev(D) + zeta div(u) I`, where the original
    /// transport sources supply mu and `bulk_viscosity_pa_s` supplies zeta at
    /// the same state. Only its molecular dissipation enters enthalpy. The
    /// selected SGS adds `2 rho nu_t dev(D)` using actual `rho = M/V`; its
    /// modeled transfer remains in the separate SGS ledger. The independent
    /// bulk coefficient never enters the closure or its molecular shear cap.
    ///
    /// The original nonlinear stress solve, reconstruction, force/projection
    /// and accepted mass-route transport are shared with the other entries.
    /// This retains the same first-order frozen-source split, exploratory
    /// policy, wall/symmetry scope and exact zero-model molecular control.
    ///
    /// # Errors
    /// Invalid or unrepresentable bulk/shear coefficients, source, model,
    /// applicability, native-work, caloric or cancellation failure publishes
    /// no candidate and leaves all inputs unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_bulk_viscosity_and_sgs(
        &self,
        domain: &VoxelDomain,
        input: MaterialMomentumAdvanceInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: &[f64],
        sgs: FvSgsConfig,
        controls: MaterialMomentumAdvanceControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousMaterialMomentumAdvance, ChtError> {
        poll(gate)?;
        sgs.validate()?;
        self.advance_viscous_material(
            domain,
            [FlowFace::Wall; 6],
            None,
            None,
            input,
            viscosity,
            Some(bulk_viscosity_pa_s),
            Some(sgs),
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
            None,
            None,
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// Add resolved-velocity SGS to the original-source reservoir advance.
    ///
    /// The opposed mouths retain normal traction and zero tangential velocity
    /// trace. Cell gradients and reciprocal mouth half-edge gradients use
    /// those same boundary kinematics; the selected Smagorinsky, WALE or
    /// Vreman closure is refreshed before each native residual check. Actual
    /// cell `rho = M/V` supplies dynamic SGS stress and each material's cap,
    /// while the original physical half-cell masses supply mouth inertia.
    ///
    /// Molecular stress retains `2 mu D`. Its heat alone enters enthalpy;
    /// the added `2 rho nu_t dev(D)` transfers resolved kinetic energy into
    /// the separate SGS ledger. Both shares use the same harmonic total
    /// stress at material jumps, including two-cell mouth half-edges. The
    /// following force/projection stage applies the original spatial pressure
    /// profile exactly once and retains its independent work account.
    ///
    /// This shares the existing first-order frozen-source split, exploratory
    /// SGS policy and opposed-reservoir scope. It supplies no incoming
    /// turbulence spectrum, evolved subgrid-energy store or wall function.
    /// A zero model coefficient is the exact molecular reservoir path.
    ///
    /// # Errors
    /// Boundary/profile authority, source, closure, local cap, nonlinear
    /// stress, energy, caloric, transport or cancellation refusal publishes
    /// no state and leaves the original inputs unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_viscosity_sgs_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        input: MaterialMomentumAdvanceInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        sgs: FvSgsConfig,
        controls: MaterialMomentumAdvanceControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousMaterialMomentumAdvance, ChtError> {
        poll(gate)?;
        sgs.validate()?;
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
            None,
            Some(sgs),
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// The independent shear/bulk response with the same opposed reservoir
    /// normal traction and zero tangential trace as the existing pressure entry.
    /// The viscous stage uses zero total normal traction, including bulk stress;
    /// the following pressure stage owns the authored reservoir load once.
    ///
    /// # Errors
    /// The bulk, boundary, source and physical-work gates are identical to the
    /// closed bulk entry and the original reservoir entry; refusal is atomic.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_bulk_viscosity_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        input: MaterialMomentumAdvanceInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: &[f64],
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
            Some(bulk_viscosity_pa_s),
            None,
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// Compose independently supplied molecular bulk viscosity, resolved SGS
    /// and the original spatial reservoir-pressure profile on one material step.
    ///
    /// Original transport sources supply mu and `bulk_viscosity_pa_s` supplies
    /// zeta at the same material state. Molecular stress and physical heating
    /// are `2 mu dev(D) + zeta div(u) I` and
    /// `2 mu dev(D):dev(D) + zeta div(u)^2`. SGS adds only
    /// `2 rho nu_t dev(D)` using actual `rho = M/V`; its resolved-energy
    /// transfer remains separate from molecular heat. The independent bulk
    /// coefficient never enters the closure or relaxes its molecular shear cap.
    ///
    /// The native predictor retains physical mouth half-cell inertia, zero
    /// total normal traction and zero tangential reservoir trace. The original
    /// pressure stage then applies the authored profile once. The same
    /// nonlinear stress solve, total-stress harmonic mouth transmission and
    /// accepted material routes are shared with the other entries.
    ///
    /// This retains the first-order frozen-source split, opposed-reservoir
    /// scope and exploratory SGS policy. Zero model coefficients take the
    /// exact original molecular shear/bulk reservoir path.
    ///
    /// # Errors
    /// Invalid bulk/shear coefficients, boundary/profile authority, source,
    /// model, local cap, native-work, caloric, transport or cancellation
    /// failure publishes no candidate and leaves all inputs unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_bulk_viscosity_sgs_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        input: MaterialMomentumAdvanceInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: &[f64],
        sgs: FvSgsConfig,
        controls: MaterialMomentumAdvanceControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousMaterialMomentumAdvance, ChtError> {
        poll(gate)?;
        sgs.validate()?;
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
            Some(bulk_viscosity_pa_s),
            Some(sgs),
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
        bulk_viscosity_pa_s: Option<&[f64]>,
        sgs: Option<FvSgsConfig>,
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
        let viscous = self.prepare_material_viscosity_mode(
            domain,
            faces,
            pressure_mouth_axis,
            input.material.mass_kg,
            input.momentum_kg_m_s,
            viscosity,
            bulk_viscosity_pa_s,
            sgs,
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
        bulk_viscosity_pa_s: Option<&[f64]>,
        reconstruction_controls: DensityMomentumReconstructionControls,
        duration_s: f64,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        admit_source: impl FnMut(usize, &FluidTransportState) -> Result<(), ChtError>,
    ) -> Result<PreparedMaterialViscosity, ChtError> {
        self.prepare_material_viscosity_mode(
            domain,
            faces,
            pressure_mouth_axis,
            mass_kg,
            momentum_kg_m_s,
            viscosity,
            bulk_viscosity_pa_s,
            None,
            reconstruction_controls,
            duration_s,
            viscosity_controls,
            gate,
            admit_source,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn prepare_material_viscosity_mode(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        pressure_mouth_axis: Option<usize>,
        mass_kg: &[f64],
        momentum_kg_m_s: &[[f64; 3]],
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: Option<&[f64]>,
        sgs: Option<FvSgsConfig>,
        reconstruction_controls: DensityMomentumReconstructionControls,
        duration_s: f64,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        mut admit_source: impl FnMut(usize, &FluidTransportState) -> Result<(), ChtError>,
    ) -> Result<PreparedMaterialViscosity, ChtError> {
        poll(gate)?;
        let n = domain.cell_count();
        if viscosity.sources.len() != n || bulk_viscosity_pa_s.is_some_and(|bulk| bulk.len() != n) {
            return Err(refuse("viscosity needs bounded original inventory"));
        }
        self.prepare_material_viscosity_with_source_mode(
            domain,
            faces,
            pressure_mouth_axis,
            mass_kg,
            momentum_kg_m_s,
            viscosity.boundaries,
            bulk_viscosity_pa_s,
            sgs,
            reconstruction_controls,
            duration_s,
            viscosity_controls,
            gate,
            |_| {
                let mut molecular = vec![0.0; n];
                for c in 0..n {
                    poll_every(gate, c)?;
                    match (domain.is_fluid(c), viscosity.sources[c]) {
                        (true, Some(source)) => {
                            let state = source.transport_state().map_err(|error| {
                                ChtError::InvalidInput {
                                    field: "flow.viscosity_source",
                                    reason: format!("cell {c}: {error}"),
                                }
                            })?;
                            admit_source(c, &state)?;
                            molecular[c] = state.dynamic_viscosity_pa_s();
                        }
                        (false, None) => {}
                        _ => {
                            return Err(refuse(
                                "viscosity sources must cover exactly the mobile cells",
                            ));
                        }
                    }
                }
                Ok(molecular)
            },
        )
    }

    /// Resolve one molecular coefficient vector from the SAME original M/P
    /// reconstruction that receives the native viscosity kick. A source owner
    /// may observe physical shear on this field before evaluating its material
    /// state, without another reconstruction or momentum solve. Fixed-solid
    /// coefficient slots are placeholders and receive no fluid dissipation.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::conjugate::flow) fn prepare_material_viscosity_with_source(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        pressure_mouth_axis: Option<usize>,
        mass_kg: &[f64],
        momentum_kg_m_s: &[[f64; 3]],
        boundaries: [FvBoundary; 6],
        bulk_viscosity_pa_s: Option<&[f64]>,
        reconstruction_controls: DensityMomentumReconstructionControls,
        duration_s: f64,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        molecular_at_field: impl FnOnce(&FlowField) -> Result<Vec<f64>, ChtError>,
    ) -> Result<PreparedMaterialViscosity, ChtError> {
        self.prepare_material_viscosity_with_source_mode(
            domain,
            faces,
            pressure_mouth_axis,
            mass_kg,
            momentum_kg_m_s,
            boundaries,
            bulk_viscosity_pa_s,
            None,
            reconstruction_controls,
            duration_s,
            viscosity_controls,
            gate,
            molecular_at_field,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn prepare_material_viscosity_with_source_mode(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        pressure_mouth_axis: Option<usize>,
        mass_kg: &[f64],
        momentum_kg_m_s: &[[f64; 3]],
        boundaries: [FvBoundary; 6],
        bulk_viscosity_pa_s: Option<&[f64]>,
        sgs: Option<FvSgsConfig>,
        reconstruction_controls: DensityMomentumReconstructionControls,
        duration_s: f64,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        molecular_at_field: impl FnOnce(&FlowField) -> Result<Vec<f64>, ChtError>,
    ) -> Result<PreparedMaterialViscosity, ChtError> {
        poll(gate)?;
        let n = domain.cell_count();
        if n > viscosity_controls.max_cells
            || n > reconstruction_controls.max_cells
            || mass_kg.len() != n
            || bulk_viscosity_pa_s.is_some_and(|bulk| bulk.len() != n)
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
        let mut molecular = molecular_at_field(&reconstructed.field)?;
        if molecular.len() != n {
            return Err(refuse("molecular sources must match the original grid"));
        }
        let mut first = None;
        for (c, &mu) in molecular.iter().enumerate() {
            poll_every(gate, c)?;
            if domain.is_fluid(c) {
                finite_positive("flow.viscosity_source_mu", mu)?;
                first.get_or_insert(mu);
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
        let face_mass =
            material_viscosity_face_masses(domain, &reconstructed, pressure_mouth_axis, gate)?;
        let native = NativeMaterialViscosityInput {
            field: &reconstructed.field,
            mass_kg,
            momentum_kg_m_s,
            viscosity_pa_s: &molecular,
            bulk_viscosity_pa_s,
            face_mass_kg: face_mass,
            boundaries,
            pressure_mouth_axis,
            duration_s,
        };
        let viscous = match sgs {
            Some(sgs) => kick_with_sgs(domain, native, sgs, viscosity_controls, gate)?,
            None => kick(domain, native, viscosity_controls, gate)?,
        };
        Ok(PreparedMaterialViscosity {
            momentum_kg_m_s: viscous.momentum_kg_m_s,
            #[cfg(feature = "implicit-enthalpy-coupling")]
            cell_eddy_viscosity_m2_s: viscous.cell_eddy_viscosity_m2_s,
            report: viscous.report,
            reconstruction: reconstructed.report,
        })
    }

    /// Resolve generalized-Newtonian shear and independent bulk response on
    /// the exact retained native endpoint velocity. Source owners validate
    /// their phase/chart/context and retain receipts in the callback.
    ///
    /// Both returned coefficient vectors have one slot per voxel; immobile
    /// slots are ignored. All native coupled sweeps, Krylov work, physical
    /// impulse, heat and numerical losses belong to this ONE original M/P
    /// solve. No scattered-cell velocity is inverted to obtain source rates.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::conjugate::flow) fn prepare_material_rheology_with_source(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        pressure_mouth_axis: Option<usize>,
        mass_kg: &[f64],
        momentum_kg_m_s: &[[f64; 3]],
        boundaries: [FvBoundary; 6],
        reconstruction_controls: DensityMomentumReconstructionControls,
        duration_s: f64,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        resolve: impl FnMut(&[f64]) -> Result<(Vec<f64>, Vec<f64>), ChtError>,
    ) -> Result<PreparedMaterialRheology, ChtError> {
        self.prepare_material_rheology_with_source_mode(
            domain,
            faces,
            pressure_mouth_axis,
            mass_kg,
            momentum_kg_m_s,
            boundaries,
            reconstruction_controls,
            duration_s,
            None,
            viscosity_controls,
            gate,
            resolve,
        )
    }

    /// Refresh source rheology and actual-density SGS on the SAME native
    /// endpoint. The source observes the unit molecular strain quadrature,
    /// without eddy-viscosity weighting; both stresses then share one residual
    /// and retain separate physical heat and modeled-transfer accounts.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::conjugate::flow) fn prepare_material_rheology_with_source_and_sgs(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        pressure_mouth_axis: Option<usize>,
        mass_kg: &[f64],
        momentum_kg_m_s: &[[f64; 3]],
        boundaries: [FvBoundary; 6],
        reconstruction_controls: DensityMomentumReconstructionControls,
        duration_s: f64,
        sgs: FvSgsConfig,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        resolve: impl FnMut(&[f64]) -> Result<(Vec<f64>, Vec<f64>), ChtError>,
    ) -> Result<PreparedMaterialRheology, ChtError> {
        poll(gate)?;
        sgs.validate()?;
        self.prepare_material_rheology_with_source_mode(
            domain,
            faces,
            pressure_mouth_axis,
            mass_kg,
            momentum_kg_m_s,
            boundaries,
            reconstruction_controls,
            duration_s,
            Some(sgs),
            viscosity_controls,
            gate,
            resolve,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_material_rheology_with_source_mode(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        pressure_mouth_axis: Option<usize>,
        mass_kg: &[f64],
        momentum_kg_m_s: &[[f64; 3]],
        boundaries: [FvBoundary; 6],
        reconstruction_controls: DensityMomentumReconstructionControls,
        duration_s: f64,
        sgs: Option<FvSgsConfig>,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        resolve: impl FnMut(&[f64]) -> Result<(Vec<f64>, Vec<f64>), ChtError>,
    ) -> Result<PreparedMaterialRheology, ChtError> {
        poll(gate)?;
        let n = domain.cell_count();
        if n > viscosity_controls.max_cells
            || n > reconstruction_controls.max_cells
            || mass_kg.len() != n
        {
            return Err(refuse("rheology needs bounded original inventory"));
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
        let face_mass =
            material_viscosity_face_masses(domain, &reconstructed, pressure_mouth_axis, gate)?;
        let native = NativeMaterialRheologyInput {
            field: &reconstructed.field,
            mass_kg,
            momentum_kg_m_s,
            face_mass_kg: face_mass,
            boundaries,
            pressure_mouth_axis,
            duration_s,
        };
        let NativeMaterialRheologyStep {
            viscosity,
            endpoint_shear_rate_per_s,
            source_updates,
        } = match sgs {
            Some(sgs) => {
                kick_with_source_and_sgs(domain, native, sgs, viscosity_controls, gate, resolve)?
            }
            None => kick_with_source(domain, native, viscosity_controls, gate, resolve)?,
        };
        Ok(PreparedMaterialRheology {
            viscosity: PreparedMaterialViscosity {
                momentum_kg_m_s: viscosity.momentum_kg_m_s,
                #[cfg(feature = "implicit-enthalpy-coupling")]
                cell_eddy_viscosity_m2_s: viscosity.cell_eddy_viscosity_m2_s,
                report: viscosity.report,
                reconstruction: reconstructed.report,
            },
            endpoint_shear_rate_per_s,
            source_updates,
        })
    }
}

/// Original density-projection half-cell inertia, shared by frozen and
/// endpoint constitutive preparation without another M/P reconstruction.
fn material_viscosity_face_masses(
    domain: &VoxelDomain,
    reconstructed: &DensityMomentumReconstruction,
    pressure_mouth_axis: Option<usize>,
    gate: &CancelGate,
) -> Result<[Vec<f64>; 3], ChtError> {
    let n = domain.cell_count();
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
    Ok(face_mass)
}

#[cfg(test)]
mod tests;
