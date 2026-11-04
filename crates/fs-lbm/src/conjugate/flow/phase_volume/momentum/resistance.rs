//! Coupled physical Darcy heating and original-chart density/volume evolution.

use std::borrow::Cow;

use fs_material::phase::EquilibriumPhaseState;
use fs_math::dd::Dd;

use super::super::super::density_projection::{
    admit_resistance_field, scatter_resistive_projected_material_momentum,
};
use super::super::super::{DensityResistanceReport, ResistiveDensityProjection, poll_every};
use super::*;

/// Bounded nonlinear heat/constitutive closure around the original volume solve.
#[derive(Clone, Copy, Debug)]
pub struct ResistiveThermalVolumeMomentumControls {
    pub thermal: ThermalVolumeMomentumControls,
    /// Maximum outer heat/resistance trials. Each retains the original bounded
    /// thermal-volume/mixing solve; rejected trials never advance physical time.
    pub max_coupling_iterations: usize,
    /// Relaxation of physical cell heat, strictly positive and at most one.
    pub heat_relaxation: f64,
    /// Maximum local physical-versus-caloric heat defect, including addition
    /// roundoff and endpoint constitutive change, J. Global defects also share
    /// the native energy/enthalpy budgets.
    pub maximum_cell_heat_defect_j: f64,
}

/// Independent coupling residuals in the accepted original-state trial.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalResistanceCouplingReport {
    pub coupling_iterations: usize,
    pub maximum_cell_heat_defect_j: f64,
    /// Sum of absolute physical-versus-used heat and caloric addition defects.
    pub absolute_heat_coupling_defect_j: f64,
    /// Non-cancelling change in drag impulse under endpoint constitutive data.
    /// Exactly zero for a prescribed frozen resistance field.
    pub absolute_constitutive_impulse_defect_kg_m_s: f64,
    /// Non-cancelling change in physical heat under endpoint constitutive data.
    pub absolute_constitutive_heat_defect_j: f64,
}

/// Accepted thermal, phase, extensive momentum and physical resistance endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ResistiveThermalVolumeMomentumStep {
    pub thermal: ThermalVolumeStep,
    pub momentum: MaterialMomentumStep,
    pub reconstruction: DensityMomentumReconstructionReport,
    pub resistance: DensityResistanceReport,
    pub feedback: MaterialDensityMomentumReport,
    /// Coefficients actually used by the accepted native mechanics, Pa s/m2.
    pub resistance_pa_s_m2: Vec<f64>,
    /// Physical heat supplied to the original caloric state in the accepted
    /// trial, J/cell. Its difference from the native physical heat is budgeted.
    pub caloric_resistance_heat_j: Vec<f64>,
    pub coupling: ThermalResistanceCouplingReport,
}

/// Constitutive data for one trial; a prescribed field borrows without copying.
pub(super) struct ResistanceResponse<'a, T> {
    pub(super) resistance_pa_s_m2: Cow<'a, [f64]>,
    pub(super) sources: T,
}

impl FlowField {
    /// Couple physical drag heating to phase density, continuity and momentum.
    ///
    /// Every trial starts from the ORIGINAL M/P/h. The original net external
    /// caloric increment plus the current physical-heat guess advances the same
    /// enthalpy chart. Its density change and nonlinear enthalpy-mixing volume
    /// determine the original continuity source. One native implicit resistance
    /// solve supplies pressure, body/drag impulses and endpoint physical heat.
    /// That heat must agree with the caloric input before anything is returned.
    ///
    /// Current M/P is reconstructed once. Only the accepted native impulse is
    /// scattered, and the accepted mass routes transport that kicked momentum.
    /// Neither trial heat nor force is accumulated across nonlinear iterations.
    /// Actual M/V supplies inertia; pressure mobility includes M + dt B.
    /// Numerical backward-Euler, reconstruction, scatter and transport losses
    /// remain separate and never enter the caloric source.
    ///
    /// This entry uses prescribed frozen nonnegative isotropic resistance.
    /// The original chart owns density and latent/sensible enthalpy. Required
    /// pressure-source and opening work retain the original NET-enthalpy
    /// convention and their explicit separate reports. No exact-solid locking,
    /// moving skeleton, pressure-dependent EOS or free surface is introduced.
    ///
    /// # Errors
    /// Invalid original data, any native refusal, an exhausted coupling budget,
    /// or cancellation returns no candidate and leaves caller fields unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_resistance(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        resistance_pa_s_m2: &[f64],
        controls: ResistiveThermalVolumeMomentumControls,
        gate: &CancelGate,
    ) -> Result<ResistiveThermalVolumeMomentumStep, ChtError> {
        self.advance_thermal_volume_resistance_with_response(
            domain,
            faces,
            None,
            curve,
            input,
            controls,
            gate,
            |_| {
                Ok(ResistanceResponse {
                    resistance_pa_s_m2: Cow::Borrowed(resistance_pa_s_m2),
                    sources: (),
                })
            },
        )
        .map(|(step, ())| step)
    }

    /// The same coupled thermal/resistance step with spatial reservoir pressure.
    /// Component-relative pressure is retained through the accepted impulse
    /// readback. Incoming material and boundary velocity remain explicit;
    /// mechanical reservoir pressure does not change the thermodynamic chart.
    ///
    /// # Errors
    /// Invalid profile or any ordinary coupled-step refusal returns no endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_piecewise_thermal_volume_momentum_with_resistance_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        resistance_pa_s_m2: &[f64],
        controls: ResistiveThermalVolumeMomentumControls,
        gate: &CancelGate,
    ) -> Result<ResistiveThermalVolumeMomentumStep, ChtError> {
        self.advance_thermal_volume_resistance_with_response(
            domain,
            faces,
            Some(BoundaryPressure::Profile(profile)),
            curve,
            input,
            controls,
            gate,
            |_| {
                Ok(ResistanceResponse {
                    resistance_pa_s_m2: Cow::Borrowed(resistance_pa_s_m2),
                    sources: (),
                })
            },
        )
        .map(|(step, ())| step)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn advance_thermal_volume_resistance_with_response<'a, T>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        boundary_pressure: Option<BoundaryPressure<'_>>,
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ResistiveThermalVolumeMomentumControls,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            &[Option<EquilibriumPhaseState>],
        ) -> Result<ResistanceResponse<'a, T>, ChtError>,
    ) -> Result<(ResistiveThermalVolumeMomentumStep, T), ChtError> {
        poll(gate)?;
        let mechanical = controls.thermal.mechanics;
        let volume = controls.thermal.volume;
        let transport_controls = volume.thermal.volume.transport;
        if mechanical.projection != transport_controls.projection {
            return Err(refuse(
                "resistance and thermal volume need identical projection controls",
            ));
        }
        if controls.max_coupling_iterations == 0
            || !controls.heat_relaxation.is_finite()
            || controls.heat_relaxation <= 0.0
            || controls.heat_relaxation > 1.0
            || !controls.maximum_cell_heat_defect_j.is_finite()
            || controls.maximum_cell_heat_defect_j < 0.0
        {
            return Err(refuse("invalid resistance/thermal coupling controls"));
        }
        admit_feedback_controls(mechanical)?;
        let n = domain.cell_count();
        if n > transport_controls.max_cells
            || n > mechanical.reconstruction.max_cells
            || n > mechanical.forcing.max_cells
            || input.volume.mass_kg.len() != n
            || input.volume.phase_states.len() != n
            || input.volume.net_enthalpy_increment_j.len() != n
            || input.momentum_kg_m_s.len() != n
            || input.force_density_n_m3.len() != n
        {
            return Err(refuse(
                "resistance/thermal inventory shape or cell budget is invalid",
            ));
        }
        for value in input.boundary_velocity_m_s.into_iter().flatten().flatten() {
            finite("flow.resistance_thermal_incoming_velocity", value)?;
        }
        let reconstructed = self.reconstruct_material_momentum(
            domain,
            faces,
            input.volume.mass_kg,
            input.momentum_kg_m_s,
            mechanical.reconstruction,
            gate,
        )?;
        let mut response = resolve(input.volume.phase_states)?;
        admit_resistance_field(domain, &response.resistance_pa_s_m2, gate)?;
        let mut heat = vec![0.0; n];
        let mut total_increment = vec![0.0; n];
        let mut addition_defect = vec![0.0; n];
        for coupling_iterations in 1..=controls.max_coupling_iterations {
            poll(gate)?;
            for c in 0..n {
                poll_every(gate, c)?;
                let external = input.volume.net_enthalpy_increment_j[c];
                finite("flow.resistance_thermal_external_increment", external)?;
                total_increment[c] = if heat[c] == 0.0 {
                    external
                } else {
                    external + heat[c]
                };
                finite(
                    "flow.resistance_thermal_total_increment",
                    total_increment[c],
                )?;
                addition_defect[c] = (Dd::from_f64(external) + Dd::from_f64(heat[c])
                    - Dd::from_f64(total_increment[c]))
                .to_f64()
                .abs();
            }
            // The existing owner performs all chart, material, incoming-state,
            // caloric representation, mixing and final occupied-volume checks.
            let (thermal, (native, pressure_readback)) =
                Self::advance_chart_volume_with_projection(
                    domain,
                    curve,
                    ChartVolumeRange::Full {
                        max_mixing_intervals: volume.max_mixing_intervals,
                    },
                    CoexistenceVolumeInput {
                        net_enthalpy_increment_j: &total_increment,
                        ..input.volume
                    },
                    volume.thermal,
                    gate,
                    |density, source| {
                        let (evolution, readback) = reconstructed
                            .field
                            .force_and_project_density_with_resistance_readback(
                                domain,
                                faces,
                                boundary_pressure,
                                DensityForceInput {
                                    density_kg_m3: density,
                                    acceleration_m_s2: input.acceleration_m_s2,
                                    force_density_n_m3: input.force_density_n_m3,
                                },
                                &response.resistance_pa_s_m2,
                                source,
                                mechanical.projection,
                                mechanical.forcing,
                                gate,
                            )?;
                        let ResistiveDensityProjection { projection, report } = evolution;
                        Ok((projection, (report, readback)))
                    },
                )?;
            let endpoint = resolve(&thermal.volume.phase_states)?;
            admit_resistance_field(domain, &endpoint.resistance_pa_s_m2, gate)?;
            let projection = &thermal.volume.transport.projection;
            let coupling = coupling_defect(
                domain,
                &projection.field,
                &response.resistance_pa_s_m2,
                &endpoint.resistance_pa_s_m2,
                mechanical.projection.duration_s,
                &native.cell_resistance_heat_j,
                &heat,
                &addition_defect,
                coupling_iterations,
                gate,
            )?;
            let heat_defect = Dd::from_f64(coupling.absolute_heat_coupling_defect_j);
            let coupling_defect =
                heat_defect + Dd::from_f64(coupling.absolute_constitutive_heat_defect_j);
            let caloric_defect =
                coupling_defect + Dd::from_f64(thermal.volume.report.absolute_enthalpy_defect_j);
            let energy_defect = caloric_defect + Dd::from_f64(native.absolute_energy_defect_j);
            let impulse_defect = Dd::from_f64(native.absolute_impulse_defect_kg_m_s)
                + Dd::from_f64(coupling.absolute_constitutive_impulse_defect_kg_m_s);
            for value in [heat_defect, caloric_defect, energy_defect, impulse_defect] {
                finite("flow.resistance_thermal_coupled_defect", value.to_f64())?;
            }
            let accepted = coupling.maximum_cell_heat_defect_j
                <= controls.maximum_cell_heat_defect_j
                && caloric_defect.to_f64() <= transport_controls.maximum_enthalpy_defect_j
                && energy_defect.to_f64() <= mechanical.forcing.maximum_energy_defect_j
                && impulse_defect.to_f64() <= mechanical.forcing.maximum_impulse_defect_kg_m_s;
            if accepted {
                let (kicked_momentum, feedback) = scatter_resistive_projected_material_momentum(
                    self,
                    domain,
                    faces,
                    boundary_pressure,
                    &pressure_readback,
                    MaterialDensityMomentumInput {
                        mass_kg: input.volume.mass_kg,
                        momentum_kg_m_s: input.momentum_kg_m_s,
                        acceleration_m_s2: input.acceleration_m_s2,
                        force_density_n_m3: input.force_density_n_m3,
                        target_outflow_m3_s: &thermal.volume.volume_source_m3_s,
                    },
                    mechanical,
                    &reconstructed,
                    projection,
                    &native,
                    &response.resistance_pa_s_m2,
                    gate,
                )?;
                let no_additional_impulse = vec![[0.0; 3]; n];
                let no_material_source = vec![None; n];
                let momentum = thermal.volume.transport.advect_momentum(
                    MaterialMomentumInput {
                        momentum_kg_m_s: &kicked_momentum,
                        impulse_kg_m_s: &no_additional_impulse,
                        source_velocity_m_s: &no_material_source,
                        boundary_velocity_m_s: input.boundary_velocity_m_s,
                    },
                    controls.thermal.momentum,
                    gate,
                )?;
                poll(gate)?;
                return Ok((
                    ResistiveThermalVolumeMomentumStep {
                        thermal,
                        momentum,
                        reconstruction: reconstructed.report,
                        resistance: native,
                        feedback,
                        resistance_pa_s_m2: response.resistance_pa_s_m2.into_owned(),
                        caloric_resistance_heat_j: heat,
                        coupling,
                    },
                    endpoint.sources,
                ));
            }
            if coupling_iterations == controls.max_coupling_iterations {
                return Err(refuse("resistance/heat/volume coupling budget exhausted"));
            }
            // Relax only a numerical heat guess. Every next thermal and native
            // mechanics trial still starts from the unchanged original state.
            for c in 0..n {
                poll_every(gate, c)?;
                heat[c] += controls.heat_relaxation * (native.cell_resistance_heat_j[c] - heat[c]);
                finite("flow.resistance_thermal_heat_guess", heat[c])?;
                if heat[c] < 0.0 {
                    return Err(refuse("physical resistance heat guess became negative"));
                }
            }
            response = endpoint;
        }
        unreachable!("positive bounded coupling loop returns an endpoint or refusal")
    }
}

/// A half-cell absolute sum bounds the face-integrated coefficient change
/// without allowing contributions from opposite sides to cancel. Held faces
/// contribute controller/heat differences too. Kinetic inertia is untouched.
#[allow(clippy::too_many_arguments)]
fn coupling_defect(
    domain: &VoxelDomain,
    field: &FlowField,
    used: &[f64],
    endpoint: &[f64],
    duration_s: f64,
    native_heat: &[f64],
    caloric_heat: &[f64],
    addition_defect: &[f64],
    coupling_iterations: usize,
    gate: &CancelGate,
) -> Result<ThermalResistanceCouplingReport, ChtError> {
    let mut impulse = Dd::ZERO;
    let mut heat = Dd::ZERO;
    let mut caloric = Dd::ZERO;
    let mut maximum_cell_heat = 0.0_f64;
    for c in 0..domain.cell_count() {
        poll_every(gate, c)?;
        if !domain.is_fluid(c) {
            continue;
        }
        let mut cell_constitutive_heat = Dd::ZERO;
        if used[c] != endpoint[c] {
            let change = (endpoint[c] - used[c]).abs();
            for face in 0..6 {
                let axis = face / 2;
                let area = domain.face_area(c, axis);
                let (slot, _) = field.slot(domain, c, face);
                let speed = (field.face(axis, slot) / area).abs();
                let coefficient = area * (0.5 * domain.widths(c)[axis] * change);
                let local_impulse = (duration_s * coefficient) * speed;
                let local_heat = local_impulse * speed;
                finite("flow.resistance_constitutive_impulse_defect", local_impulse)?;
                finite("flow.resistance_constitutive_heat_defect", local_heat)?;
                impulse = impulse + Dd::from_f64(local_impulse);
                cell_constitutive_heat = cell_constitutive_heat + Dd::from_f64(local_heat);
            }
        }
        heat = heat + cell_constitutive_heat;
        let residual = (Dd::from_f64(native_heat[c]) - Dd::from_f64(caloric_heat[c]))
            .to_f64()
            .abs();
        let local = Dd::from_f64(residual) + Dd::from_f64(addition_defect[c]);
        caloric = caloric + local;
        let total = (local + cell_constitutive_heat).to_f64();
        finite("flow.resistance_thermal_cell_heat_defect", total)?;
        maximum_cell_heat = maximum_cell_heat.max(total);
    }
    finite(
        "flow.resistance_total_constitutive_impulse_defect",
        impulse.to_f64(),
    )?;
    finite(
        "flow.resistance_total_constitutive_heat_defect",
        heat.to_f64(),
    )?;
    finite("flow.resistance_total_caloric_defect", caloric.to_f64())?;
    Ok(ThermalResistanceCouplingReport {
        coupling_iterations,
        maximum_cell_heat_defect_j: maximum_cell_heat,
        absolute_heat_coupling_defect_j: caloric.to_f64(),
        absolute_constitutive_impulse_defect_kg_m_s: impulse.to_f64(),
        absolute_constitutive_heat_defect_j: heat.to_f64(),
    })
}

#[cfg(test)]
mod tests;
