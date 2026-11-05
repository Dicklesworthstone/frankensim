//! One original molecular kick supplies the same implicit caloric endpoint.

use super::*;
use crate::conjugate::flow::{
    MaterialMomentumViscosityInput, MaterialViscosityControls, StaggeredBoundaryPressure,
    ThermalVolumeMomentumInput, ViscousThermalVolumeMomentumStep,
};

#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConductiveViscousThermalVolumeMomentumStep {
    pub evolution: ViscousThermalVolumeMomentumStep,
    pub conductive_enthalpy_increment_j: Vec<f64>,
    /// The shared thermal defect also includes molecular heat-addition roundoff.
    pub report: ImplicitConductiveThermalVolumeReport,
}

impl FlowField {
    /// Prepare molecular momentum and physical heat once, solve that heat with
    /// source conduction and baths, then advance original volume and momentum.
    ///
    /// Viscosity retains its original frozen-phase sampling and actual M/V
    /// inertia. Authored and molecular heat enter the same implicit k(T,f)/latent
    /// endpoint; its accepted heat drives native flow using the prepared cell
    /// momentum. Pressure work and numerical viscous losses stay separate.
    /// This is the existing first-order split, without a thermal-viscosity
    /// constitutive iteration. The optional reservoir retains the original
    /// opposed mechanical mouths and pressure profile.
    ///
    /// # Errors
    /// Original source, viscosity, implicit heat, pressure, material, work-budget
    /// and cancellation refusals return no endpoint or changed input inventory.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conductive_thermal_volume_momentum_with_viscosity_and_baths(
        &self,
        domain: &VoxelDomain,
        reservoir: Option<([FlowFace; 6], StaggeredBoundaryPressure<'_>)>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        input: ThermalVolumeMomentumInput<'_>,
        viscosity: MaterialMomentumViscosityInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        viscosity_controls: MaterialViscosityControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConductiveViscousThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        let (faces, profile, mouth) = match reservoir {
            Some((faces, profile)) => (
                faces,
                Some(profile),
                Some(Self::admit_material_viscosity_pressure_mouth(
                    faces,
                    viscosity.boundaries,
                    input.boundary_velocity_m_s,
                )?),
            ),
            None => ([FlowFace::Wall; 6], None, None),
        };
        let (viscous, heat_budget) = self.prepare_thermal_volume_viscosity(
            domain,
            faces,
            mouth,
            conductivity.caloric_curve(),
            input,
            viscosity,
            controls.evolution,
            viscosity_controls,
            gate,
        )?;
        let (combined, molecular_assembly) = Self::assemble_thermal_volume_viscous_heat(
            domain,
            input.volume.net_enthalpy_increment_j,
            &viscous.report,
            viscosity_controls,
            heat_budget,
            gate,
        )?;
        let heated = ThermalVolumeMomentumInput {
            volume: CoexistenceVolumeInput {
                net_enthalpy_increment_j: &combined,
                ..input.volume
            },
            momentum_kg_m_s: &viscous.momentum_kg_m_s,
            ..input
        };
        let mut accepted = self.implicit_conduct_then_advance(
            domain,
            conductivity,
            baths,
            heated.volume,
            ConductiveThermalVolumeControls {
                volume: controls.evolution.volume,
                max_conductive_faces: controls.max_conductive_faces,
            },
            solve,
            workspace,
            gate,
            |volume| {
                let input = ThermalVolumeMomentumInput { volume, ..heated };
                match profile {
                    Some(profile) => self
                        .advance_piecewise_thermal_volume_momentum_with_pressure_profile(
                            domain,
                            faces,
                            profile,
                            conductivity.caloric_curve(),
                            input,
                            controls.evolution,
                            gate,
                        ),
                    None => self.advance_piecewise_thermal_volume_momentum(
                        domain,
                        faces,
                        conductivity.caloric_curve(),
                        input,
                        controls.evolution,
                        gate,
                    ),
                }
            },
            |evolution| &evolution.thermal,
        )?;
        let viscous_defect = (molecular_assembly
            + Dd::from_f64(
                accepted
                    .evolution
                    .thermal
                    .volume
                    .report
                    .absolute_enthalpy_defect_j,
            ))
        .to_f64();
        let total = (molecular_assembly + Dd::from_f64(accepted.report.absolute_enthalpy_defect_j))
            .to_f64();
        finite("flow.implicit_viscous_enthalpy_defect", viscous_defect)?;
        finite("flow.implicit_viscous_total_enthalpy_defect", total)?;
        if total > heat_budget || viscous_defect > heat_budget {
            return Err(refuse(
                "implicit conduction and molecular heat exceed shared enthalpy budget",
            ));
        }
        accepted.report.absolute_enthalpy_defect_j = total;
        poll(gate)?;
        Ok(ImplicitConductiveViscousThermalVolumeMomentumStep {
            evolution: ViscousThermalVolumeMomentumStep {
                evolution: accepted.evolution,
                viscosity_reconstruction: viscous.reconstruction,
                viscosity: viscous.report,
                absolute_heat_assembly_defect_j: molecular_assembly.to_f64(),
                absolute_enthalpy_defect_j: viscous_defect,
            },
            conductive_enthalpy_increment_j: accepted.increments,
            report: accepted.report,
        })
    }
}
