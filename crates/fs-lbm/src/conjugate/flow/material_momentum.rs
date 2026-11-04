//! One evolving-material momentum step on original native force and mass owners.

use super::density_transport::advect_projected_material;
use super::{
    CancelGate, ChtError, DensityForceReport, DensityMomentumReconstructionReport, FlowFace,
    FlowField, ForcedDensityProjection, ForcedDensityProjectionReport, MaterialAdvectionControls,
    MaterialAdvectionInput, MaterialAdvectionStep, MaterialDensityMomentumControls,
    MaterialDensityMomentumInput, MaterialDensityMomentumReport, MaterialDensityMomentumStep,
    MaterialMomentumControls, MaterialMomentumInput, MaterialMomentumStep, MaterialVolumeSource,
    ReservoirPressure, StaggeredBoundaryPressure, VoxelDomain, finite, finite_positive, poll,
    poll_every,
};

mod resistance;
pub use resistance::{
    PhaseResistanceMomentumControls, PhaseResistanceMomentumInput,
    ResistiveMaterialMomentumAdvance, ResolvedPhaseResistanceMomentumStep,
};
mod viscous;
pub(in crate::conjugate::flow) use viscous::PreparedMaterialViscosity;
pub use viscous::{
    MaterialMomentumViscosityInput, MaterialViscosityControls, MaterialViscosityReport,
    ViscousMaterialMomentumAdvance,
};

/// Original material inventory, extensive momentum and authored exchanges.
#[derive(Debug, Clone, Copy)]
pub struct MaterialMomentumAdvanceInput<'a> {
    /// Mass and specific enthalpy in original cells, with physical volume sources.
    pub material: MaterialAdvectionInput<'a>,
    /// Current extensive momentum in world axes, kg m/s; zero in solids.
    pub momentum_kg_m_s: &'a [[f64; 3]],
    /// Uniform physical acceleration, m/s^2, including gravity if present.
    pub acceleration_m_s2: [f64; 3],
    /// Additional body force, N/m^3; excludes pressure and is zero in solids.
    pub force_density_n_m3: &'a [[f64; 3]],
    /// Required for every local material inflow, including when total injected
    /// and extracted volume cancel. Extraction uses the actual donor velocity.
    pub source_velocity_m_s: &'a [Option<[f64; 3]>],
    /// Required on each actually incoming boundary; outflow uses its donor.
    pub boundary_velocity_m_s: [Option<[f64; 3]>; 6],
}

/// Existing independent stage budgets for one first-order split step.
#[derive(Debug, Clone, Copy)]
pub struct MaterialMomentumAdvanceControls {
    pub mechanics: MaterialDensityMomentumControls,
    pub material: MaterialAdvectionControls,
    pub momentum: MaterialMomentumControls,
}

/// Endpoint cell inventory and the original, separately audited stage reports.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialMomentumAdvance {
    /// Endpoint mass and specific enthalpy; extensive enthalpy is their product.
    /// The retained projection is the FORCE/PROJECT TRANSPORT-STAGE field which
    /// produced the accepted mass routes. It is not an endpoint reconstruction
    /// of the returned mass and momentum, and need not have their velocity.
    pub transport: MaterialAdvectionStep,
    /// Actual endpoint extensive momentum and its actual P/M velocity.
    pub momentum: MaterialMomentumStep,
    pub reconstruction: DensityMomentumReconstructionReport,
    pub forcing: DensityForceReport,
    pub force_projection: ForcedDensityProjectionReport,
    pub feedback: MaterialDensityMomentumReport,
}

impl FlowField {
    /// Advance material mass, enthalpy and momentum together with one fresh
    /// native force/pressure solve and the original immutable mass routes.
    ///
    /// Interior and free faces are reconstructed from the supplied CURRENT M/P
    /// on every call. `self` supplies held boundary fluxes only. The native
    /// half-cell force and pressure impulses reach the original cell momentum
    /// exactly once; the accepted mass routes then transport that momentum with
    /// zero additional impulse. No second projection or density solve is used
    /// to manufacture a separate momentum flux.
    ///
    /// The pressure volume target comes from each original material source:
    /// injection is positive, extraction negative. Opposed sources retain their
    /// separate mass, enthalpy and momentum exchanges even at zero net volume.
    /// Mechanics and material must contain exactly matching projection controls,
    /// including duration and acceptance gates; no stage silently overrides one.
    ///
    /// This is a first-order force/pressure/advection split, with explicit
    /// source and boundary momentum. It does not add viscosity, an EOS, thermal
    /// expansion or phase mobility. Native constraint/projection/scatter and
    /// donor-transport losses remain separate numerical accounts and never heat
    /// the material. This bounded composition is not a full Navier-Stokes step.
    ///
    /// For another step, supply the returned cell mass, enthalpy and momentum.
    /// The transport-stage field may supply the next held boundary values, but
    /// its interior/free velocities are not reused as endpoint momentum.
    ///
    /// # Errors
    /// A native stage refusal, mismatched clocks/gates, missing incoming data,
    /// invalid work budget, failed balance readback or cancellation returns no
    /// candidate. All supplied fields and inventories remain unchanged.
    pub fn advance_material_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        input: MaterialMomentumAdvanceInput<'_>,
        controls: MaterialMomentumAdvanceControls,
        gate: &CancelGate,
    ) -> Result<MaterialMomentumAdvance, ChtError> {
        self.advance_material_momentum_with_reservoirs(domain, faces, None, input, controls, gate)
    }

    /// Advance the same conserved material inventories under reservoir pressures
    /// on the actual staggered open faces. Each pressure impulse acts through
    /// the original current-density half-cell inertia, reaches cell momentum
    /// once, and drives the SAME mass routes carrying enthalpy and momentum.
    /// Spatial heads may balance gravity or drive inflow, outflow and backflow.
    ///
    /// The borrowed profile has the admission and gauge convention of
    /// [`Self::project_density_volume_source_with_pressure_profile`]. It supplies
    /// mechanical correction pressure, not a thermodynamic EOS state. Every
    /// actual incoming boundary still needs its explicit material and velocity.
    /// Reservoir work stays in the native mechanical reports; it does not heat
    /// the liquid or alter the supplied enthalpy reference. An enclosing total-
    /// energy owner must retain its existing pressure-flow/enthalpy convention.
    ///
    /// The first-order split, stage budgets and atomic refusal guarantees of
    /// [`Self::advance_material_momentum`] remain in force. No second pressure
    /// solve or rebuilt mass-route operator is introduced.
    ///
    /// # Errors
    /// An invalid pressure profile or any original mechanics, material,
    /// momentum, incoming-state or cancellation refusal returns no candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        input: MaterialMomentumAdvanceInput<'_>,
        controls: MaterialMomentumAdvanceControls,
        gate: &CancelGate,
    ) -> Result<MaterialMomentumAdvance, ChtError> {
        self.advance_material_momentum_with_reservoir_pressure(
            domain, faces, profile, input, controls, gate,
        )
    }

    /// Advance the original material owners with either six optional reservoir
    /// plane pressures in Pa or a borrowed spatial profile. Uniform reservoirs
    /// require no staggered pressure buffers. The existing profile method and
    /// this entry share one force/projection/feedback stage and one set of
    /// accepted mass routes, with no extra pressure or momentum solve.
    ///
    /// Opening pressure work stays in `transport.projection.report` and the
    /// composed kinetic accounts. Material advection still carries exactly
    /// `M*h`, with explicit incoming state and velocity for every actual inflow;
    /// pressure work is not added again to enthalpy. These remain separate
    /// mechanical and enthalpy accounts, without a new EOS or total-energy law.
    ///
    /// # Errors
    /// Nonfinite pressure, pressure on a held plane, invalid profile, or any
    /// refusal of [`Self::advance_material_momentum_with_pressure_profile`].
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_reservoir_pressure<'p>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: impl Into<ReservoirPressure<'p>>,
        input: MaterialMomentumAdvanceInput<'_>,
        controls: MaterialMomentumAdvanceControls,
        gate: &CancelGate,
    ) -> Result<MaterialMomentumAdvance, ChtError> {
        self.advance_material_momentum_with_reservoirs(
            domain,
            faces,
            Some(reservoir_pressure.into()),
            input,
            controls,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_material_momentum_with_reservoirs(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: Option<ReservoirPressure<'_>>,
        input: MaterialMomentumAdvanceInput<'_>,
        controls: MaterialMomentumAdvanceControls,
        gate: &CancelGate,
    ) -> Result<MaterialMomentumAdvance, ChtError> {
        poll(gate)?;
        if controls.mechanics.projection != controls.material.projection {
            return Err(refuse(
                "mechanics and material need identical projection controls",
            ));
        }
        let n = domain.cell_count();
        if n > controls.mechanics.reconstruction.max_cells
            || n > controls.mechanics.forcing.max_cells
            || n > controls.material.max_cells
            || input.material.mass_kg.len() != n
            || input.material.specific_enthalpy_j_kg.len() != n
            || input.material.volume_sources.len() != n
            || input.momentum_kg_m_s.len() != n
            || input.force_density_n_m3.len() != n
            || input.source_velocity_m_s.len() != n
        {
            return Err(refuse("inventory shape or native cell budget is invalid"));
        }
        for value in input.boundary_velocity_m_s.into_iter().flatten().flatten() {
            finite("flow.advance_incoming_velocity", value)?;
        }
        let mut target_outflow = vec![0.0; n];
        for (c, target) in target_outflow.iter_mut().enumerate() {
            poll_every(gate, c)?;
            for value in input.source_velocity_m_s[c].into_iter().flatten() {
                finite("flow.advance_source_velocity", value)?;
            }
            *target = match input.material.volume_sources[c] {
                MaterialVolumeSource::None => 0.0,
                MaterialVolumeSource::Inflow {
                    volume_rate_m3_s, ..
                } => {
                    finite_positive("flow.advance_injected_volume_rate", volume_rate_m3_s)?;
                    if input.source_velocity_m_s[c].is_none() {
                        return Err(refuse("incoming material source lacks velocity"));
                    }
                    volume_rate_m3_s
                }
                MaterialVolumeSource::Outflow { volume_rate_m3_s } => {
                    finite_positive("flow.advance_removed_volume_rate", volume_rate_m3_s)?;
                    -volume_rate_m3_s
                }
            };
        }
        let mechanical_input = MaterialDensityMomentumInput {
            mass_kg: input.material.mass_kg,
            momentum_kg_m_s: input.momentum_kg_m_s,
            acceleration_m_s2: input.acceleration_m_s2,
            force_density_n_m3: input.force_density_n_m3,
            target_outflow_m3_s: &target_outflow,
        };
        let mechanics = match reservoir_pressure {
            Some(reservoir_pressure) => self
                .force_project_material_momentum_with_reservoir_pressure(
                    domain,
                    faces,
                    reservoir_pressure,
                    mechanical_input,
                    controls.mechanics,
                    gate,
                )?,
            None => self.force_project_material_momentum(
                domain,
                faces,
                mechanical_input,
                controls.mechanics,
                gate,
            )?,
        };
        let MaterialDensityMomentumStep {
            momentum_kg_m_s: kicked_momentum,
            density_kg_m3: _,
            reconstruction,
            evolution:
                ForcedDensityProjection {
                    projection,
                    forcing,
                    report: force_projection,
                },
            report: feedback,
        } = mechanics;
        let transport =
            advect_projected_material(domain, projection, input.material, controls.material, gate)?;
        let no_additional_impulse = vec![[0.0; 3]; n];
        let momentum = transport.advect_momentum(
            MaterialMomentumInput {
                momentum_kg_m_s: &kicked_momentum,
                impulse_kg_m_s: &no_additional_impulse,
                source_velocity_m_s: input.source_velocity_m_s,
                boundary_velocity_m_s: input.boundary_velocity_m_s,
            },
            controls.momentum,
            gate,
        )?;
        poll(gate)?;
        Ok(MaterialMomentumAdvance {
            transport,
            momentum,
            reconstruction,
            forcing,
            force_projection,
            feedback,
        })
    }
}

fn refuse(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.material_momentum_advance",
        reason: reason.into(),
    }
}

#[cfg(test)]
mod tests;
