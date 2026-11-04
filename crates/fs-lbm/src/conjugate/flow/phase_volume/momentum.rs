//! Original-chart volume change, pressure impulse and conservative cell momentum.

use fs_material::phase::EquilibriumEnthalpyPhaseCurve;

use super::super::density_projection::{
    BoundaryPressure, admit_feedback_controls, audit_combined, scatter_projected_material_momentum,
};
use super::super::{
    CancelGate, ChtError, DensityForceInput, DensityForceReport,
    DensityMomentumReconstructionReport, FlowFace, FlowField, ForcedDensityProjectionReport,
    MaterialDensityMomentumControls, MaterialDensityMomentumInput, MaterialDensityMomentumReport,
    MaterialMomentumControls, MaterialMomentumInput, MaterialMomentumStep,
    StaggeredBoundaryPressure, VoxelDomain, finite, poll,
};
use super::{
    ChartVolumeRange, CoexistenceVolumeInput, PiecewiseThermalVolumeControls, ThermalVolumeStep,
    refuse,
};

/// Authoritative thermal inventory and momentum before the same physical step.
#[derive(Clone, Copy, Debug)]
pub struct ThermalVolumeMomentumInput<'a> {
    /// Original chart states, masses and prescribed NET enthalpy increments.
    /// Thermal expansion is a volume source; it injects no material mass.
    pub volume: CoexistenceVolumeInput<'a>,
    /// Current extensive cell momentum in world axes, kg m/s; zero in solids.
    pub momentum_kg_m_s: &'a [[f64; 3]],
    /// Uniform physical acceleration, m/s^2, including gravity if present.
    pub acceleration_m_s2: [f64; 3],
    /// Additional body force, N/m^3; excludes pressure and is zero in solids.
    pub force_density_n_m3: &'a [[f64; 3]],
    /// Required on every actually incoming boundary. Outflow carries the
    /// actual donor momentum after the accepted native force/pressure impulse.
    pub boundary_velocity_m_s: [Option<[f64; 3]>; 6],
}

/// Original thermal, mechanical and momentum-route acceptance controls.
#[derive(Clone, Copy, Debug)]
pub struct ThermalVolumeMomentumControls {
    pub volume: PiecewiseThermalVolumeControls,
    /// Its projection controls must exactly equal the thermal-volume controls,
    /// including physical duration and all continuity/work acceptance budgets.
    pub mechanics: MaterialDensityMomentumControls,
    pub momentum: MaterialMomentumControls,
}

/// One accepted thermal, phase, mass and authoritative momentum endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalVolumeMomentumStep {
    /// Original chart states and conservative M/Mh transport. Its face field is
    /// the accepted transport-stage field, not an endpoint P/M reconstruction.
    pub thermal: ThermalVolumeStep,
    /// Endpoint extensive momentum, transported on those SAME mass routes.
    pub momentum: MaterialMomentumStep,
    pub reconstruction: DensityMomentumReconstructionReport,
    pub forcing: DensityForceReport,
    pub force_projection: ForcedDensityProjectionReport,
    pub feedback: MaterialDensityMomentumReport,
}

impl FlowField {
    /// Advance original-chart thermal volume and authoritative momentum together.
    /// Current M/P reconstructs interior/free face momentum once; `self` supplies
    /// held boundary fluxes only. Body forcing uses actual current M/V inertia.
    /// All bounded thermal-volume trials start from that same frozen predictor
    /// and the same once-applied net caloric increment. Only the accepted trial's
    /// pressure impulse is audited and returned to the original cell momenta.
    ///
    /// The original accepted mass routes then carry this kicked momentum with
    /// zero additional impulse. Thermal expansion/contraction and nonlinear
    /// enthalpy mixing supply continuity volume without creating material mass.
    /// There is no second pressure solve, projection clone or inverse velocity
    /// reconstruction. Subsequent steps consume the returned M, h, phase and P.
    ///
    /// The caller owns phase mobility and the net enthalpy/work convention.
    /// Required source-pressure work, held/reservoir work and each numerical
    /// projection/scatter/transport loss remain in their original separate
    /// reports. They never become heat automatically. This first-order split
    /// adds no solid stress, free surface, viscosity or pressure-dependent EOS.
    ///
    /// # Errors
    /// Any original chart, force, pressure, material, momentum or cancellation
    /// refusal returns no candidate. Mismatched stage controls and missing
    /// incoming material or velocity also refuse. Inputs remain unchanged.
    pub fn advance_piecewise_thermal_volume_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ThermalVolumeMomentumControls,
        gate: &CancelGate,
    ) -> Result<ThermalVolumeMomentumStep, ChtError> {
        self.advance_piecewise_thermal_volume_momentum_with_reservoirs(
            domain, faces, None, curve, input, controls, gate,
        )
    }

    /// Evolve the same thermal-volume/momentum step under the actual staggered
    /// reservoir pressure profile. Pressure acts once through the accepted
    /// chart-compatible native projection and its current-density inertia.
    /// Component-relative pressure is retained through the impulse/work audit,
    /// preserving small physical heads under a large common ambient reference.
    ///
    /// The profile has the admission and mechanical-pressure convention of
    /// [`Self::project_density_volume_source_with_pressure_profile`]. Actual
    /// incoming material and velocity remain explicit. Reservoir work is
    /// reported separately; it does not change the source chart's enthalpy
    /// reference or supply an unmodeled thermodynamic pressure response.
    ///
    /// # Errors
    /// Invalid pressure data or any original thermal, mechanical, transport or
    /// cancellation refusal returns no candidate and changes no supplied state.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ThermalVolumeMomentumControls,
        gate: &CancelGate,
    ) -> Result<ThermalVolumeMomentumStep, ChtError> {
        self.advance_piecewise_thermal_volume_momentum_with_reservoirs(
            domain,
            faces,
            Some(profile),
            curve,
            input,
            controls,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_piecewise_thermal_volume_momentum_with_reservoirs(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ThermalVolumeMomentumControls,
        gate: &CancelGate,
    ) -> Result<ThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        let transport_controls = controls.volume.thermal.volume.transport;
        if controls.mechanics.projection != transport_controls.projection {
            return Err(refuse(
                "thermal-volume and momentum need identical projection controls",
            ));
        }
        let n = domain.cell_count();
        if n > transport_controls.max_cells
            || n > controls.mechanics.reconstruction.max_cells
            || n > controls.mechanics.forcing.max_cells
            || input.volume.mass_kg.len() != n
            || input.volume.phase_states.len() != n
            || input.volume.net_enthalpy_increment_j.len() != n
            || input.momentum_kg_m_s.len() != n
            || input.force_density_n_m3.len() != n
        {
            return Err(refuse(
                "thermal-volume momentum inventory shape or cell budget is invalid",
            ));
        }
        admit_feedback_controls(controls.mechanics)?;
        for value in input.boundary_velocity_m_s.into_iter().flatten().flatten() {
            finite("flow.thermal_volume_incoming_velocity", value)?;
        }
        let reconstructed = self.reconstruct_material_momentum(
            domain,
            faces,
            input.volume.mass_kg,
            input.momentum_kg_m_s,
            controls.mechanics.reconstruction,
            gate,
        )?;
        let force_input = DensityForceInput {
            density_kg_m3: &reconstructed.density_kg_m3,
            acceleration_m_s2: input.acceleration_m_s2,
            force_density_n_m3: input.force_density_n_m3,
        };
        let predicted = reconstructed.field.predict_density_body_force(
            domain,
            faces,
            force_input,
            controls.mechanics.projection.duration_s,
            controls.mechanics.forcing,
            gate,
        )?;
        let (thermal, pressure_readback) = Self::advance_chart_volume_with_projection(
            domain,
            curve,
            ChartVolumeRange::Full {
                max_mixing_intervals: controls.volume.max_mixing_intervals,
            },
            input.volume,
            controls.volume.thermal,
            gate,
            |density, source| match profile {
                None => predicted
                    .field
                    .project_density_volume_source(
                        domain,
                        faces,
                        density,
                        source,
                        transport_controls.projection,
                        gate,
                    )
                    .map(|projection| (projection, None)),
                Some(profile) => predicted
                    .field
                    .project_density_volume_source_with_pressure_readback(
                        domain,
                        faces,
                        BoundaryPressure::Profile(profile),
                        density,
                        source,
                        transport_controls.projection,
                        gate,
                    )
                    .map(|(projection, readback)| (projection, Some(readback))),
            },
        )?;
        let projection = &thermal.volume.transport.projection;
        let force_projection = audit_combined(
            &reconstructed.field,
            &predicted,
            projection,
            domain,
            faces,
            profile.map(BoundaryPressure::Profile),
            pressure_readback.as_ref(),
            force_input,
            controls.mechanics.projection.duration_s,
            controls.mechanics.forcing,
            gate,
        )?;
        let (kicked_momentum, feedback) = scatter_projected_material_momentum(
            self,
            domain,
            faces,
            profile.map(BoundaryPressure::Profile),
            pressure_readback.as_ref(),
            MaterialDensityMomentumInput {
                mass_kg: input.volume.mass_kg,
                momentum_kg_m_s: input.momentum_kg_m_s,
                acceleration_m_s2: input.acceleration_m_s2,
                force_density_n_m3: input.force_density_n_m3,
                target_outflow_m3_s: &thermal.volume.volume_source_m3_s,
            },
            controls.mechanics,
            &reconstructed,
            projection,
            gate,
        )?;
        let no_additional_impulse = vec![[0.0; 3]; n];
        let no_material_source_velocity = vec![None; n];
        let momentum = thermal.volume.transport.advect_momentum(
            MaterialMomentumInput {
                momentum_kg_m_s: &kicked_momentum,
                impulse_kg_m_s: &no_additional_impulse,
                source_velocity_m_s: &no_material_source_velocity,
                boundary_velocity_m_s: input.boundary_velocity_m_s,
            },
            controls.momentum,
            gate,
        )?;
        poll(gate)?;
        Ok(ThermalVolumeMomentumStep {
            thermal,
            momentum,
            reconstruction: reconstructed.report,
            forcing: predicted.report,
            force_projection,
            feedback,
        })
    }
}

#[cfg(test)]
mod tests;
