//! Native physical shear and original phase sources on evolving liquid density.

mod endpoint;
pub use endpoint::EndpointRheologyThermalVolumeMomentumStep;
#[cfg(feature = "implicit-enthalpy-coupling")]
pub(in crate::conjugate::flow) use endpoint::resolve_endpoint_rheology_sources;

use super::*;
use crate::conjugate::FvBoundary;
use crate::conjugate::flow::poll_every;
use crate::conjugate::simple::viscous_heat::material_momentum::physical_shear_rates;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::phase::transport::ResolvedPhaseTransportState;

/// A sourced molecular step with its original physical operating coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedViscousThermalVolumeMomentumStep {
    pub evolution: ViscousThermalVolumeMomentumStep,
    /// Cell RMS sqrt(2 dev(D):dev(D)), in 1/s, on the native stress quadrature.
    /// Fixed-solid cells are zero; these are not endpoint rates.
    pub original_shear_rate_per_s: Vec<f64>,
    /// One original phase/material snapshot per mobile cell, None in fixed
    /// solids. Source receipts keep their own authority and support domain.
    pub original_sources: Vec<Option<ResolvedPhaseTransportState>>,
}

impl FlowField {
    /// Query source rheology at the physical shear of the current mass/momentum
    /// inventory, then advance its molecular heat and changing liquid density.
    ///
    /// `resolve` receives each mobile cell's ORIGINAL phase state and native
    /// RMS sqrt(2 dev(D):dev(D)) exactly once. Bind these coordinates through
    /// `fs_material::phase::transport::resolve_phase_transport_state_with_shear_for_volume_change`
    /// to retain source-selected density, viscosity, conductivity and auxiliary
    /// caloric receipts. The returned source must contain the exact supplied
    /// phase state and current chart density. The resolver owns its fixed
    /// pressure/composition context and must consume the supplied rate.
    ///
    /// Shear is observed on the same M/P reconstruction that enters viscosity,
    /// including native walls and fixed obstacles. The unit molecular operator
    /// removes viscosity weighting from the observation. Source queries and
    /// viscosity run once, before bounded thermal-volume trials; no additional
    /// momentum solve or inverse velocity reconstruction is introduced.
    ///
    /// This source-driven entry admits fully liquid mobile material only.
    /// Original, declared incoming, calorically heated and transported states
    /// must remain liquid;
    /// a nonzero declared Darcy resistance or solid branch requires a separate
    /// phase-support owner. Sensible thermal expansion and its compatible
    /// volume flux retain the original chart. Viscosity is frozen at the
    /// original rate for this first-order step, so temporal refinement remains
    /// necessary for nonlinear rheology. Numerical losses remain separate from
    /// physical heat and thermodynamic pressure work stays caller-owned.
    ///
    /// # Errors
    /// Source, phase-support, native solve, conservation, budget or cancellation
    /// refusal returns no endpoint. A retry starts from unchanged inventories.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_source_rheology(
        &self,
        domain: &VoxelDomain,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedViscousThermalVolumeMomentumStep, ChtError> {
        self.advance_source_rheology_thermal_volume(
            domain,
            [FlowFace::Wall; 6],
            None,
            None,
            curve,
            input,
            boundaries,
            None,
            controls,
            viscosity_controls,
            gate,
            resolve,
        )
    }

    /// Resolve original source shear together with independent bulk viscosity.
    ///
    /// The bulk slice has one finite nonnegative coefficient [Pa s] per voxel;
    /// fixed-solid slots are ignored. The caller evaluates it at the same
    /// original material state as `resolve`. Native shear is still observed
    /// and resolved exactly once per mobile cell. The frozen predictor uses
    /// `2 mu dev(D) + zeta div(u) I` and deposits its matching physical heat
    /// into the same caloric/volume owner. No bulk law is inferred from shear.
    ///
    /// # Errors
    /// Invalid bulk response or any original source, phase-support, solve,
    /// thermal or cancellation refusal returns no endpoint or changed input.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_source_rheology_and_bulk_viscosity(
        &self,
        domain: &VoxelDomain,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        bulk_viscosity_pa_s: &[f64],
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedViscousThermalVolumeMomentumStep, ChtError> {
        self.advance_source_rheology_thermal_volume(
            domain,
            [FlowFace::Wall; 6],
            None,
            None,
            curve,
            input,
            boundaries,
            Some(bulk_viscosity_pa_s),
            controls,
            viscosity_controls,
            gate,
            resolve,
        )
    }

    /// The same source-driven liquid step with the original opposed pressure
    /// mouths. The physical rate includes their native zero tangential trace;
    /// pressure acts once through the accepted thermal-volume projection.
    /// Inflow phase and velocity remain explicit through either flow direction.
    ///
    /// # Errors
    /// Unsupported boundaries, any source refusal or an endpoint needing
    /// unrepresented phase support returns no partial state.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_source_rheology_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedViscousThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        let axis = Self::admit_material_viscosity_pressure_mouth(
            faces,
            boundaries,
            input.boundary_velocity_m_s,
        )?;
        self.advance_source_rheology_thermal_volume(
            domain,
            faces,
            Some(profile),
            Some(axis),
            curve,
            input,
            boundaries,
            None,
            controls,
            viscosity_controls,
            gate,
            resolve,
        )
    }

    /// Source-resolved shear and independent bulk viscosity with the original
    /// opposed pressure mouths and matching physical half-cell inertia.
    ///
    /// The viscous predictor imposes zero total normal traction, including
    /// bulk stress. The subsequent original pressure stage owns the reservoir
    /// load once. Bulk coefficients obey the closed bulk entry's state,
    /// units and frozen-coefficient convention.
    ///
    /// # Errors
    /// Source, bulk, boundary, phase-support, solve and cancellation refusals
    /// retain the original all-or-nothing endpoint publication.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_source_rheology_and_bulk_viscosity_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        bulk_viscosity_pa_s: &[f64],
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedViscousThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        let axis = Self::admit_material_viscosity_pressure_mouth(
            faces,
            boundaries,
            input.boundary_velocity_m_s,
        )?;
        self.advance_source_rheology_thermal_volume(
            domain,
            faces,
            Some(profile),
            Some(axis),
            curve,
            input,
            boundaries,
            Some(bulk_viscosity_pa_s),
            controls,
            viscosity_controls,
            gate,
            resolve,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_source_rheology_thermal_volume(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        pressure_mouth_axis: Option<usize>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        boundaries: [FvBoundary; 6],
        bulk_viscosity_pa_s: Option<&[f64]>,
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
            f64,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedViscousThermalVolumeMomentumStep, ChtError> {
        let heat_budget = admit_viscous_thermal_input(domain, input, controls, gate)?;
        let n = domain.cell_count();
        if input
            .volume
            .boundary_inflow
            .iter()
            .flatten()
            .any(|phase| phase.liquid_mass_fraction() != 1.0)
        {
            return Err(refuse(
                "source rheology requires fully liquid incoming material",
            ));
        }
        let mut original_shear_rate_per_s = Vec::new();
        let mut original_sources = Vec::new();
        let viscous = self.prepare_material_viscosity_with_source(
            domain,
            faces,
            pressure_mouth_axis,
            input.volume.mass_kg,
            input.momentum_kg_m_s,
            boundaries,
            bulk_viscosity_pa_s,
            controls.mechanics.reconstruction,
            controls.mechanics.projection.duration_s,
            viscosity_controls,
            gate,
            |field| {
                original_shear_rate_per_s = physical_shear_rates(
                    domain,
                    field,
                    boundaries,
                    pressure_mouth_axis,
                    viscosity_controls.max_cells,
                    gate,
                )?;
                original_sources = Vec::with_capacity(n);
                let mut molecular = vec![0.0; n];
                for cell in 0..n {
                    poll_every(gate, cell)?;
                    if !domain.is_fluid(cell) {
                        original_sources.push(None);
                        continue;
                    }
                    let phase = input.volume.phase_states[cell].ok_or_else(|| {
                        refuse("source rheology needs an original phase state in each mobile cell")
                    })?;
                    if phase.phase_curve_identity() != curve.identity()
                        || phase.material_card_identity() != curve.material_card_identity()
                        || phase.liquid_mass_fraction() != 1.0
                    {
                        return Err(refuse(
                            "source rheology needs the original fully liquid chart state",
                        ));
                    }
                    let source = resolve(cell, &phase, original_shear_rate_per_s[cell])?;
                    if *source.phase_state() != phase
                        || source.source_density_kg_m3().to_bits()
                            != phase.bulk_density_kg_m3().to_bits()
                    {
                        return Err(refuse(
                            "rheology source must match the exact original phase and density",
                        ));
                    }
                    if source.mobility().mobile().is_none_or(|state| {
                        state.inverse_permeability_per_m2 != 0.0 || state.volumes.liquid != 1.0
                    }) {
                        return Err(refuse(
                            "source rheology requires liquid with no unrepresented Darcy support",
                        ));
                    }
                    molecular[cell] = source
                        .mobile_transport()
                        .ok_or_else(|| refuse("liquid source has no mobile transport"))?
                        .dynamic_viscosity_pa_s();
                    original_sources.push(Some(source));
                }
                Ok(molecular)
            },
        )?;
        // The molecular and authored heat use the same addition order as the
        // original caloric owner. Advection cannot hide a temporary departure
        // into a phase whose support this liquid-only entry does not represent.
        for cell in 0..n {
            poll_every(gate, cell)?;
            if domain.is_fluid(cell) {
                let phase = input.volume.phase_states[cell]
                    .ok_or_else(|| refuse("missing original liquid state"))?;
                let heat = input.volume.net_enthalpy_increment_j[cell]
                    + viscous.report.cell_molecular_heat_j[cell];
                let h = phase.specific_enthalpy_j_kg() + heat / input.volume.mass_kg[cell];
                let heated = curve.values_at_specific_enthalpy(h).map_err(|error| {
                    ChtError::InvalidInput {
                        field: "flow.rheology_heated_phase",
                        reason: error.to_string(),
                    }
                })?;
                if heated.liquid_mass_fraction != 1.0 {
                    return Err(refuse("rheology heating leaves the admitted liquid branch"));
                }
            }
        }
        let evolution = self.finish_viscous_thermal_volume(
            domain,
            faces,
            profile,
            curve,
            input,
            controls,
            viscosity_controls,
            viscous,
            heat_budget,
            gate,
        )?;
        for (cell, phase) in evolution
            .evolution
            .thermal
            .volume
            .phase_states
            .iter()
            .enumerate()
        {
            poll_every(gate, cell)?;
            if phase.is_some_and(|phase| phase.liquid_mass_fraction() != 1.0) {
                return Err(refuse(
                    "rheology transport leaves the admitted liquid branch",
                ));
            }
        }
        poll(gate)?;
        Ok(ResolvedViscousThermalVolumeMomentumStep {
            evolution,
            original_shear_rate_per_s,
            original_sources,
        })
    }
}

#[cfg(test)]
mod tests;
