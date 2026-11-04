//! Implicit stationary resistance, physical heat and the original mass routes.

use super::super::{
    DensityResistanceReport, ResistiveDensityProjection, ResistiveMaterialDensityMomentumStep,
};
use super::*;
use fs_math::dd::Dd;

/// One atomic implicit-resistance/pressure/material-momentum candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct ResistiveMaterialMomentumAdvance {
    /// The same native transport-stage projection used for all accepted routes;
    /// endpoint M and P need not reconstruct its divergence-constrained field.
    pub transport: MaterialAdvectionStep,
    pub momentum: MaterialMomentumStep,
    pub reconstruction: DensityMomentumReconstructionReport,
    /// Actual-density native mechanics, including physical drag heat, pressure
    /// and support impulses. Pressure-source/opening work is in `transport`.
    pub resistance: DensityResistanceReport,
    pub feedback: MaterialDensityMomentumReport,
    /// Actual sum M*(h_heated-h_original) before the accepted material transport.
    pub represented_resistance_heat_j: f64,
    /// Sum of absolute local insertion defects. This shares the original
    /// native energy budget and the material enthalpy budget with their stages.
    pub absolute_heat_insertion_defect_j: f64,
}

impl FlowField {
    /// Advance current mass, enthalpy and momentum through frozen stationary
    /// Darcy resistance using one actual-density native pressure solve.
    ///
    /// Each original half-cell contributes its volume-integrated beta to the
    /// implicit face resistance. The accepted drag impulse reaches original
    /// cell P once through actual M/2 shares. Only dt B u_final^2 enters h, in
    /// the native report's cell partition, before the original mass routes
    /// transport M, M*h and kicked momentum. No numerical backward-Euler,
    /// pressure, reconstruction, scatter or advection loss is converted to heat.
    ///
    /// `resistance_pa_s_m2` is supplied physical beta at the current material
    /// state, nonnegative in fluids and zero in solids. It does not replace
    /// inertial density, lock a phase, or define an EOS. This first-order split
    /// uses fixed geometry and frozen resistance. Thermodynamic density/phase
    /// evolution remains with the existing material owner.
    ///
    /// Source and boundary material/velocity are explicit, including backflow.
    /// Every stage retains its original budget and cancellation contract; any
    /// failure returns no candidate and leaves all supplied state unchanged.
    ///
    /// # Errors
    /// Invalid resistance, input, clock or incoming state; an unrepresentable
    /// heat insertion; or any native mechanical/transport/cancellation refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_resistance(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        input: MaterialMomentumAdvanceInput<'_>,
        resistance_pa_s_m2: &[f64],
        controls: MaterialMomentumAdvanceControls,
        gate: &CancelGate,
    ) -> Result<ResistiveMaterialMomentumAdvance, ChtError> {
        self.advance_material_momentum_with_resistance_reservoirs(
            domain,
            faces,
            None,
            input,
            resistance_pa_s_m2,
            controls,
            gate,
        )
    }

    /// The same resistance-aware advance with uniform or spatial reservoir
    /// pressure. One native pressure impulse and one set of accepted mass
    /// routes carry all material and momentum exchanges. Reservoir work is
    /// mechanical and is not deposited again as enthalpy or resistance heat.
    ///
    /// # Errors
    /// Invalid pressure or any refusal of
    /// [`Self::advance_material_momentum_with_resistance`].
    #[allow(clippy::too_many_arguments)]
    pub fn advance_material_momentum_with_resistance_and_reservoir_pressure<'p>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: impl Into<ReservoirPressure<'p>>,
        input: MaterialMomentumAdvanceInput<'_>,
        resistance_pa_s_m2: &[f64],
        controls: MaterialMomentumAdvanceControls,
        gate: &CancelGate,
    ) -> Result<ResistiveMaterialMomentumAdvance, ChtError> {
        self.advance_material_momentum_with_resistance_reservoirs(
            domain,
            faces,
            Some(reservoir_pressure.into()),
            input,
            resistance_pa_s_m2,
            controls,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_material_momentum_with_resistance_reservoirs(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: Option<ReservoirPressure<'_>>,
        input: MaterialMomentumAdvanceInput<'_>,
        resistance_pa_s_m2: &[f64],
        controls: MaterialMomentumAdvanceControls,
        gate: &CancelGate,
    ) -> Result<ResistiveMaterialMomentumAdvance, ChtError> {
        poll(gate)?;
        if controls.mechanics.projection != controls.material.projection {
            return Err(refuse(
                "resistance mechanics and material need identical projection controls",
            ));
        }
        finite(
            "flow.resistance_enthalpy_budget",
            controls.material.maximum_enthalpy_defect_j,
        )?;
        if controls.material.maximum_enthalpy_defect_j < 0.0 {
            return Err(refuse("material enthalpy budget must be nonnegative"));
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
            || resistance_pa_s_m2.len() != n
        {
            return Err(refuse("inventory shape or native cell budget is invalid"));
        }
        for value in input.boundary_velocity_m_s.into_iter().flatten().flatten() {
            finite("flow.resistance_incoming_velocity", value)?;
        }
        let mut target_outflow = vec![0.0; n];
        for (c, target) in target_outflow.iter_mut().enumerate() {
            poll_every(gate, c)?;
            for value in input.source_velocity_m_s[c].into_iter().flatten() {
                finite("flow.resistance_source_velocity", value)?;
            }
            *target = match input.material.volume_sources[c] {
                MaterialVolumeSource::None => 0.0,
                MaterialVolumeSource::Inflow {
                    volume_rate_m3_s, ..
                } => {
                    finite_positive("flow.resistance_injected_volume_rate", volume_rate_m3_s)?;
                    if input.source_velocity_m_s[c].is_none() {
                        return Err(refuse("incoming material source lacks velocity"));
                    }
                    volume_rate_m3_s
                }
                MaterialVolumeSource::Outflow { volume_rate_m3_s } => {
                    finite_positive("flow.resistance_removed_volume_rate", volume_rate_m3_s)?;
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
            Some(pressure) => self
                .force_project_material_momentum_with_resistance_and_reservoir_pressure(
                    domain,
                    faces,
                    pressure,
                    mechanical_input,
                    resistance_pa_s_m2,
                    controls.mechanics,
                    gate,
                )?,
            None => self.force_project_material_momentum_with_resistance(
                domain,
                faces,
                mechanical_input,
                resistance_pa_s_m2,
                controls.mechanics,
                gate,
            )?,
        };
        let ResistiveMaterialDensityMomentumStep {
            momentum_kg_m_s: kicked_momentum,
            density_kg_m3: _,
            reconstruction,
            evolution:
                ResistiveDensityProjection {
                    projection,
                    report: resistance,
                },
            report: feedback,
        } = mechanics;
        let mut enthalpy = input.material.specific_enthalpy_j_kg.to_vec();
        let mut represented = Dd::ZERO;
        let mut insertion_defect = Dd::ZERO;
        for c in 0..n {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                continue;
            }
            finite("flow.resistance_original_enthalpy", enthalpy[c])?;
            let heat = resistance.cell_resistance_heat_j[c];
            if heat != 0.0 {
                enthalpy[c] += heat / input.material.mass_kg[c];
            }
            finite("flow.resistance_heated_enthalpy", enthalpy[c])?;
            let actual = Dd::from_f64(input.material.mass_kg[c])
                * (Dd::from_f64(enthalpy[c])
                    - Dd::from_f64(input.material.specific_enthalpy_j_kg[c]));
            represented = represented + actual;
            insertion_defect =
                insertion_defect + Dd::from_f64((actual - Dd::from_f64(heat)).to_f64().abs());
        }
        let defect = insertion_defect.to_f64();
        finite("flow.resistance_heat_insertion_defect", defect)?;
        let native_defect =
            (insertion_defect + Dd::from_f64(resistance.absolute_energy_defect_j)).to_f64();
        finite("flow.resistance_heated_energy_defect", native_defect)?;
        if native_defect > controls.mechanics.forcing.maximum_energy_defect_j {
            return Err(refuse(
                "actual resistance heating exceeds the native energy budget",
            ));
        }
        if defect > controls.material.maximum_enthalpy_defect_j {
            return Err(refuse(
                "actual resistance heating exceeds the material enthalpy budget",
            ));
        }
        let represented = represented.to_f64();
        finite("flow.resistance_represented_heat", represented)?;
        let transport = advect_projected_material(
            domain,
            projection,
            MaterialAdvectionInput {
                specific_enthalpy_j_kg: &enthalpy,
                ..input.material
            },
            controls.material,
            gate,
        )?;
        let caloric_defect =
            (insertion_defect + Dd::from_f64(transport.report.absolute_enthalpy_defect_j)).to_f64();
        finite("flow.resistance_advected_enthalpy_defect", caloric_defect)?;
        if caloric_defect > controls.material.maximum_enthalpy_defect_j {
            return Err(refuse(
                "heat insertion and transport exceed the shared material enthalpy budget",
            ));
        }
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
        Ok(ResistiveMaterialMomentumAdvance {
            transport,
            momentum,
            reconstruction,
            resistance,
            feedback,
            represented_resistance_heat_j: represented,
            absolute_heat_insertion_defect_j: defect,
        })
    }
}

#[cfg(test)]
mod tests;
