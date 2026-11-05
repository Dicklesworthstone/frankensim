//! The accepted implicit face heat drives the original material momentum once.

use super::*;
use crate::conjugate::flow::{
    StaggeredBoundaryPressure, ThermalVolumeMomentumInput, ThermalVolumeMomentumStep,
};

#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConductiveThermalVolumeMomentumStep {
    pub evolution: ThermalVolumeMomentumStep,
    pub conductive_enthalpy_increment_j: Vec<f64>,
    pub report: ImplicitConductiveThermalVolumeReport,
}

impl FlowField {
    /// Apply one accepted implicit source-conduction/bath interval to the
    /// original phase-density volume and material-momentum evolution.
    ///
    /// Heat and its limits follow
    /// [`Self::advance_isotropic_implicit_conductive_thermal_volume_with_baths`].
    /// The optional existing mechanical pressure profile, actual mass routes,
    /// caller-supplied inflow velocity and single accepted pressure impulse keep
    /// their original owner. Bath heat adds no mass, pressure work or momentum.
    ///
    /// # Errors
    /// Any original heat, source, nonlinear, volume, momentum or cancellation
    /// refusal returns no endpoint and leaves the supplied inventories unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conductive_thermal_volume_momentum_with_baths(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConductiveThermalVolumeMomentumStep, ChtError> {
        let accepted = self.implicit_conduct_then_advance(
            domain,
            conductivity,
            baths,
            input.volume,
            ConductiveThermalVolumeControls {
                volume: controls.evolution.volume,
                max_conductive_faces: controls.max_conductive_faces,
            },
            solve,
            workspace,
            gate,
            |volume| {
                let heated = ThermalVolumeMomentumInput { volume, ..input };
                match profile {
                    Some(profile) => self
                        .advance_piecewise_thermal_volume_momentum_with_pressure_profile(
                            domain,
                            faces,
                            profile,
                            conductivity.caloric_curve(),
                            heated,
                            controls.evolution,
                            gate,
                        ),
                    None => self.advance_piecewise_thermal_volume_momentum(
                        domain,
                        faces,
                        conductivity.caloric_curve(),
                        heated,
                        controls.evolution,
                        gate,
                    ),
                }
            },
            |evolution| &evolution.thermal,
        )?;
        Ok(ImplicitConductiveThermalVolumeMomentumStep {
            evolution: accepted.evolution,
            conductive_enthalpy_increment_j: accepted.increments,
            report: accepted.report,
        })
    }
}
