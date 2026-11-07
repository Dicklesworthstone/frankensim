//! Original shear and optional bulk stress heat the chart that owns material volume.

mod rheology;
pub use rheology::ResolvedViscousThermalVolumeMomentumStep;

use super::*;
use crate::conjugate::flow::material_momentum::PreparedMaterialViscosity;
use crate::conjugate::flow::{
    MaterialMomentumViscosityInput, MaterialViscosityControls, MaterialViscosityReport,
};
use fs_math::dd::Dd;

/// One admitted molecular-stress, caloric-density and momentum endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscousThermalVolumeMomentumStep {
    pub evolution: ThermalVolumeMomentumStep,
    /// Original M/P sampling for viscosity, before the thermal/pressure step.
    /// The later reconstruction applies its own held-face reaction once.
    pub viscosity_reconstruction: DensityMomentumReconstructionReport,
    /// Includes the physical cell heat and separate temporal/scatter losses.
    pub viscosity: MaterialViscosityReport,
    /// Absolute local roundoff from adding molecular heat to authored net heat.
    /// This is included in both the viscosity work and thermal enthalpy gates.
    pub absolute_heat_assembly_defect_j: f64,
    /// Assembly plus the accepted original caloric/transport defect, checked
    /// against the same native thermal-volume enthalpy budget.
    pub absolute_enthalpy_defect_j: f64,
}

pub(in crate::conjugate::flow) fn admit_viscous_thermal_input(
    domain: &VoxelDomain,
    input: ThermalVolumeMomentumInput<'_>,
    controls: ThermalVolumeMomentumControls,
    gate: &CancelGate,
) -> Result<f64, ChtError> {
    poll(gate)?;
    let n = domain.cell_count();
    let transport = controls.volume.thermal.volume.transport;
    if n > transport.max_cells
        || n > controls.mechanics.forcing.max_cells
        || controls.mechanics.projection != transport.projection
        || input.volume.mass_kg.len() != n
        || input.volume.phase_states.len() != n
        || input.volume.net_enthalpy_increment_j.len() != n
        || input.momentum_kg_m_s.len() != n
        || input.force_density_n_m3.len() != n
    {
        return Err(refuse(
            "viscous thermal evolution needs bounded original inventory and identical clocks",
        ));
    }
    let heat_budget = transport.maximum_enthalpy_defect_j;
    if !heat_budget.is_finite() || heat_budget < 0.0 {
        return Err(refuse("invalid combined molecular thermal enthalpy budget"));
    }
    Ok(heat_budget)
}

impl FlowField {
    /// Evolve molecular momentum diffusion and its density-changing heat on
    /// the original material chart and accepted conservative mass routes.
    ///
    /// The existing symmetric-stress owner samples each original viscosity
    /// source once and uses actual M/V face inertia. Each sampled density must
    /// exactly match the supplied original phase state. The caller evaluates
    /// the remaining source coordinates and admits phase mobility; a raw
    /// transport tuple retains raw authority. This is frozen source viscosity,
    /// not an implicit constitutive iteration or a solid-support model.
    ///
    /// Physical molecular dissipation is added once to NET caloric increments,
    /// retaining the original phase states as the initial chart state. The
    /// existing thermal owner then advances h, phase and density, derives the
    /// compatible volume source, and returns the accepted pressure impulse to
    /// viscously updated cell momentum. Its same mass routes carry M, M*h and P.
    /// Thus molecular heat can expand or shrink occupied volume without an
    /// invented material mass source or a stale post-heating density.
    ///
    /// The first-order viscosity/thermal/pressure/advection split retains
    /// moving-wall work and physical heat separately from all numerical losses.
    /// Required source-pressure work keeps its original caller-owned caloric
    /// convention; it is not silently added to or subtracted from heat. Fixed
    /// solid voxels remain stationary. A sealed geometry that cannot accommodate
    /// its chart's net expansion refuses instead of inventing compression.
    ///
    /// # Errors
    /// Source/chart disagreement, unsupported boundaries, mismatched clocks,
    /// heat representation, any native balance/budget or cancellation refusal
    /// returns no endpoint and leaves every input unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_viscosity(
        &self,
        domain: &VoxelDomain,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousThermalVolumeMomentumStep, ChtError> {
        self.advance_viscous_thermal_volume(
            domain,
            [FlowFace::Wall; 6],
            None,
            None,
            curve,
            input,
            viscosity,
            None,
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// Use independent nonnegative bulk viscosity [Pa s] in the original
    /// thermal-volume split. Shear remains owned by `viscosity.sources`;
    /// the caller resolves one bulk coefficient per voxel at the same state.
    /// Fixed-solid slots are ignored. Stress is `2 mu dev(D) + zeta div(u) I`.
    /// Its matching shear plus bulk heat enters the original caloric chart
    /// once, so density/occupied-volume feedback uses that same physical loss.
    /// The original entry retains its `2 mu D` constitutive assumption.
    ///
    /// # Errors
    /// Invalid bulk response or any original source/chart, work, volume,
    /// transport or cancellation refusal returns no endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_bulk_viscosity(
        &self,
        domain: &VoxelDomain,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: &[f64],
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousThermalVolumeMomentumStep, ChtError> {
        self.advance_viscous_thermal_volume(
            domain, [FlowFace::Wall; 6], None, None, curve, input, viscosity,
            Some(bulk_viscosity_pa_s), controls, viscosity_controls, gate,
        )
    }

    /// The same viscous thermal-density step through the existing opposed
    /// pressure-reservoir pair. Normal traction and zero tangential trace use
    /// the native boundary half-cell inertia. Viscosity applies zero pressure;
    /// the subsequent accepted thermal-volume projection applies the profile
    /// once. Incoming chart material and normal velocity remain explicit.
    ///
    /// # Errors
    /// Any original pressure-mouth or viscous thermal-step refusal returns no
    /// candidate. In particular, missing incoming momentum can refuse after
    /// heat/pressure acceptance without publishing a partial heated state.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_viscosity_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        let axis = Self::admit_material_viscosity_pressure_mouth(
            faces,
            viscosity.boundaries,
            input.boundary_velocity_m_s,
        )?;
        self.advance_viscous_thermal_volume(
            domain,
            faces,
            Some(profile),
            Some(axis),
            curve,
            input,
            viscosity,
            None,
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// Independent shear/bulk stress and caloric-volume feedback through the
    /// existing opposed reservoir pair. Zero total viscous normal traction
    /// and the original half-cell inertia include the bulk response; authored
    /// reservoir pressure remains owned by the subsequent pressure stage.
    ///
    /// # Errors
    /// The same bulk, source/chart, boundary and physical-balance gates as the
    /// closed bulk entry and original pressure entry apply atomically.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_bulk_viscosity_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: &[f64],
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        let axis = Self::admit_material_viscosity_pressure_mouth(
            faces, viscosity.boundaries, input.boundary_velocity_m_s,
        )?;
        self.advance_viscous_thermal_volume(
            domain, faces, Some(profile), Some(axis), curve, input, viscosity,
            Some(bulk_viscosity_pa_s), controls, viscosity_controls, gate,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_viscous_thermal_volume(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        pressure_mouth_axis: Option<usize>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: Option<&[f64]>,
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ViscousThermalVolumeMomentumStep, ChtError> {
        let (viscous, heat_budget) = self.prepare_thermal_volume_viscosity(
            domain,
            faces,
            pressure_mouth_axis,
            curve,
            input,
            viscosity,
            bulk_viscosity_pa_s,
            controls,
            viscosity_controls,
            gate,
        )?;
        self.finish_viscous_thermal_volume(
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
        )
    }

    /// Share the original frozen source sampling and molecular kick with
    /// caloric owners that must solve its physical heat before volume transport.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::conjugate::flow) fn prepare_thermal_volume_viscosity(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        pressure_mouth_axis: Option<usize>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        bulk_viscosity_pa_s: Option<&[f64]>,
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<(PreparedMaterialViscosity, f64), ChtError> {
        let heat_budget = admit_viscous_thermal_input(domain, input, controls, gate)?;
        let transport = controls.volume.thermal.volume.transport;
        let viscous = self.prepare_material_viscosity(
            domain, faces, pressure_mouth_axis, input.volume.mass_kg,
            input.momentum_kg_m_s, viscosity, bulk_viscosity_pa_s,
            controls.mechanics.reconstruction,
            transport.projection.duration_s, viscosity_controls, gate,
            |cell, source| {
                let phase = input.volume.phase_states[cell]
                    .ok_or_else(|| refuse("viscosity needs the original phase state in every mobile cell"))?;
                if phase.phase_curve_identity() != curve.identity()
                    || phase.material_card_identity() != curve.material_card_identity()
                    || source.density_kg_m3().to_bits() != phase.bulk_density_kg_m3().to_bits()
                {
                    return Err(refuse("viscosity source density must match the original phase chart at this state"));
                }
                Ok(())
            },
        )?;
        Ok((viscous, heat_budget))
    }

    /// Add physical molecular heat once, retaining the exact original local
    /// addition residual and joint molecular-work/thermal acceptance gates.
    pub(in crate::conjugate::flow) fn assemble_thermal_volume_viscous_heat(
        domain: &VoxelDomain,
        authored_heat_j: &[f64],
        viscosity: &MaterialViscosityReport,
        viscosity_controls: MaterialViscosityControls,
        heat_budget: f64,
        gate: &CancelGate,
    ) -> Result<(Vec<f64>, Dd), ChtError> {
        let n = domain.cell_count();
        let mut combined = Vec::with_capacity(n);
        let mut assembly = Dd::ZERO;
        for cell in 0..n {
            super::super::poll_every(gate, cell)?;
            let authored = authored_heat_j[cell];
            finite("flow.viscous_thermal_authored_heat", authored)?;
            let heat = viscosity.cell_molecular_heat_j[cell];
            let sum = authored + heat;
            finite("flow.viscous_thermal_combined_heat", sum)?;
            let residual = (Dd::from_f64(authored) + Dd::from_f64(heat) - Dd::from_f64(sum))
                .to_f64()
                .abs();
            finite("flow.viscous_thermal_heat_addition_residual", residual)?;
            assembly = assembly + Dd::from_f64(residual);
            combined.push(sum);
        }
        let assembly_defect = assembly.to_f64();
        finite("flow.viscous_thermal_heat_assembly_defect", assembly_defect)?;
        if assembly_defect > heat_budget
            || (assembly + Dd::from_f64(viscosity.absolute_energy_defect_j)).to_f64()
                > viscosity_controls.maximum_energy_defect_j
        {
            return Err(refuse(
                "molecular heat assembly exceeds the shared thermal or viscosity work budget",
            ));
        }
        Ok((combined, assembly))
    }

    /// Consume the prepared impulse and physical heat without reconstructing
    /// or sampling its original molecular source a second time.
    #[allow(clippy::too_many_arguments)]
    fn finish_viscous_thermal_volume(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        viscous: PreparedMaterialViscosity,
        heat_budget: f64,
        gate: &CancelGate,
    ) -> Result<ViscousThermalVolumeMomentumStep, ChtError> {
        let (combined, assembly) = Self::assemble_thermal_volume_viscous_heat(
            domain,
            input.volume.net_enthalpy_increment_j,
            &viscous.report,
            viscosity_controls,
            heat_budget,
            gate,
        )?;
        let assembly_defect = assembly.to_f64();
        let heated = ThermalVolumeMomentumInput {
            volume: CoexistenceVolumeInput {
                net_enthalpy_increment_j: &combined,
                ..input.volume
            },
            momentum_kg_m_s: &viscous.momentum_kg_m_s,
            ..input
        };
        let evolution = match profile {
            Some(profile) => self.advance_piecewise_thermal_volume_momentum_with_pressure_profile(
                domain, faces, profile, curve, heated, controls, gate,
            )?,
            None => self.advance_piecewise_thermal_volume_momentum(
                domain, faces, curve, heated, controls, gate,
            )?,
        };
        let absolute_enthalpy_defect_j = (assembly
            + Dd::from_f64(evolution.thermal.volume.report.absolute_enthalpy_defect_j))
        .to_f64();
        finite(
            "flow.viscous_thermal_total_enthalpy_defect",
            absolute_enthalpy_defect_j,
        )?;
        if absolute_enthalpy_defect_j > heat_budget {
            return Err(refuse(
                "molecular heat and native volume transport exceed the shared enthalpy budget",
            ));
        }
        poll(gate)?;
        Ok(ViscousThermalVolumeMomentumStep {
            evolution,
            viscosity_reconstruction: viscous.reconstruction,
            viscosity: viscous.report,
            absolute_heat_assembly_defect_j: assembly_defect,
            absolute_enthalpy_defect_j,
        })
    }
}

#[cfg(test)]
mod tests;
