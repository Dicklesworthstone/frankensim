//! Source conduction and molecular heat on one original caloric-volume owner.

use super::*;
use crate::conjugate::flow::{
    MaterialMomentumViscosityInput, MaterialViscosityControls, ViscousThermalVolumeMomentumStep,
};

/// Original conductive heat and one accepted viscous thermal-density endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductiveViscousThermalVolumeMomentumStep {
    pub evolution: ViscousThermalVolumeMomentumStep,
    /// Original-state paired conductive heat, before all other caloric sources.
    pub conductive_enthalpy_increment_j: Vec<f64>,
    /// Its shared enthalpy defect includes BOTH heat-assembly stages and the
    /// original caloric/transport defect, without counting either stage twice.
    pub report: ConductiveThermalVolumeReport,
}

impl FlowField {
    /// Conduct source heat and diffuse molecular momentum before evolving the
    /// resulting chart density, compatible volume flux and material momentum.
    ///
    /// Conductivity and viscosity use the same original phase states. The
    /// original isotropic conductivity owner assembles paired face heat once;
    /// the original viscosity owner samples mu once, uses actual M/V inertia,
    /// and supplies only its physical molecular dissipation to caloric heat.
    /// Authored, conductive and molecular increments enter the existing phase
    /// chart together. Native phase-volume trials neither repeat conduction
    /// nor apply viscosity, body forcing or pressure impulse repeatedly.
    ///
    /// Both heat-addition residuals and the complete caloric/transport residual
    /// share ONE enthalpy budget. Wall work, required source-pressure work and
    /// numerical losses remain in their existing mechanical accounts. The
    /// caller owns phase mobility and the caloric/work convention. This is the
    /// existing first-order frozen-source split, with insulated exterior and
    /// fixed-solid thermal faces, and requires time/space refinement.
    ///
    /// # Errors
    /// Any source, heat-step, chart, viscosity, volume, momentum or cancellation
    /// refusal returns no endpoint and changes none of the supplied state.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_conductive_thermal_volume_momentum_with_viscosity(
        &self,
        domain: &VoxelDomain,
        conductivity: &PreparedPhaseConductivity<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ConductiveViscousThermalVolumeMomentumStep, ChtError> {
        self.conduct_viscous_thermal_momentum(
            domain,
            None,
            conductivity,
            input,
            viscosity,
            controls,
            viscosity_controls,
            gate,
        )
    }

    /// The same source-conductive/viscous density evolution with the existing
    /// opposed reservoir mouths and original staggered pressure profile.
    /// Pressure enters the accepted native projection once; incoming material
    /// and normal velocity remain explicit through either flow direction.
    ///
    /// # Errors
    /// Unsupported boundary combinations or any original coupled-stage refusal
    /// returns no partial endpoint. A retry starts from unchanged inventories.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_conductive_thermal_volume_momentum_with_viscosity_and_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        conductivity: &PreparedPhaseConductivity<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ConductiveViscousThermalVolumeMomentumStep, ChtError> {
        self.conduct_viscous_thermal_momentum(
            domain,
            Some((faces, profile)),
            conductivity,
            input,
            viscosity,
            controls,
            viscosity_controls,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn conduct_viscous_thermal_momentum(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        gate: &CancelGate,
    ) -> Result<ConductiveViscousThermalVolumeMomentumStep, ChtError> {
        let accepted = self.conduct_then_advance(
            domain,
            conductivity,
            input.volume,
            ConductiveThermalVolumeControls {
                volume: controls.evolution.volume,
                max_conductive_faces: controls.max_conductive_faces,
            },
            gate,
            |volume| {
                let heated = ThermalVolumeMomentumInput { volume, ..input };
                let evolution = match reservoir {
                    Some((faces, profile)) => self
                        .advance_piecewise_thermal_volume_momentum_with_viscosity_and_pressure_profile(
                            domain, faces, profile, conductivity.caloric_curve(), heated,
                            viscosity, controls.evolution, viscosity_controls, gate,
                        ),
                    None => self.advance_piecewise_thermal_volume_momentum_with_viscosity(
                        domain, conductivity.caloric_curve(), heated, viscosity,
                        controls.evolution, viscosity_controls, gate,
                    ),
                }?;
                // Molecular heat is added inside the viscous owner, after
                // conduction's original caloric admission. Check that actual
                // intermediate state too: later advection must not hide a
                // departure from this original conductivity source's support.
                let bounds = conductivity.enthalpy_bounds_j_kg();
                for (cell, phase) in input.volume.phase_states.iter().enumerate() {
                    poll_every(gate, cell)?;
                    if let Some(phase) = phase {
                        let heat = volume.net_enthalpy_increment_j[cell]
                            + evolution.viscosity.cell_molecular_heat_j[cell];
                        let h = phase.specific_enthalpy_j_kg()
                            + heat / input.volume.mass_kg[cell];
                        if !(bounds[0]..=bounds[1]).contains(&h) {
                            return Err(refuse(
                                "molecular heating leaves the original conductivity source",
                            ));
                        }
                    }
                }
                Ok(evolution)
            },
            |evolution| (
                &evolution.evolution.thermal,
                evolution.absolute_heat_assembly_defect_j,
            ),
        )?;
        Ok(ConductiveViscousThermalVolumeMomentumStep {
            evolution: accepted.evolution,
            conductive_enthalpy_increment_j: accepted.increments,
            report: accepted.report,
        })
    }
}
