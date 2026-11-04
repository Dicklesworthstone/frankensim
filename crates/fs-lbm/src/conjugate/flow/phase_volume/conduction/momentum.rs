//! Original-source heat, thermal volume and actual material momentum together.

mod viscous;
pub use viscous::ConductiveViscousThermalVolumeMomentumStep;

use super::super::{
    ThermalVolumeMomentumControls, ThermalVolumeMomentumInput, ThermalVolumeMomentumStep,
};
use super::*;
use crate::conjugate::flow::StaggeredBoundaryPressure;

/// Existing native evolution budgets plus the original explicit conduction cap.
#[derive(Clone, Copy, Debug)]
pub struct ConductiveThermalVolumeMomentumControls {
    pub evolution: ThermalVolumeMomentumControls,
    pub max_conductive_faces: usize,
}

/// One admitted source-heat, phase-density, volume and momentum endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductiveThermalVolumeMomentumStep {
    /// Original phase/volume result and the momentum it actually transported.
    pub evolution: ThermalVolumeMomentumStep,
    /// Paired original-source heat, before authored net caloric increments [J].
    pub conductive_enthalpy_increment_j: Vec<f64>,
    /// Same source identity, heat, work cap and shared enthalpy gates as the
    /// original conduction owner. Mechanical source work remains separate.
    pub report: ConductiveThermalVolumeReport,
}

impl FlowField {
    /// Conduct original material heat, evolve its chart density and return the
    /// resulting pressure impulse to authoritative material momentum.
    ///
    /// The same original-source face heat and enthalpy-secant monotonicity gate
    /// as [`Self::advance_isotropic_conductive_thermal_volume`] run once. Their
    /// increments enter the original phase/volume trials. Only the accepted
    /// pressure trial supplies cell impulse; the same mass routes then carry
    /// that momentum, mass and mass*h. A trial is not another physical tick.
    ///
    /// The source remains explicitly isotropic and evaluated at the original
    /// phase states. Exterior and fixed-solid thermal faces are insulated.
    /// Pressure/source work and all numerical losses retain their existing
    /// mechanical accounts and are not deposited as heat. The caller supplies
    /// admitted phase mobility and the work associated with authored net
    /// caloric increments; this does not infer solid stress or free-surface
    /// motion. Original first-order split/refinement requirements remain.
    ///
    /// # Errors
    /// Any source, heat-step, chart, volume, momentum, work or cancellation
    /// refusal returns no candidate and changes none of the supplied state.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_conductive_thermal_volume_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        conductivity: &PreparedPhaseConductivity<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        gate: &CancelGate,
    ) -> Result<ConductiveThermalVolumeMomentumStep, ChtError> {
        self.conduct_thermal_momentum(domain, faces, None, conductivity, input, controls, gate)
    }

    /// The same coupled source-heat step with the existing borrowed mechanical
    /// pressure profile on actual open faces. Inflow/backflow phase states and
    /// velocities remain explicit. Pressure uses the native component gauge;
    /// it is not a thermodynamic EOS update or an additional heat increment.
    ///
    /// # Errors
    /// An invalid profile or any original conductive/thermal/momentum refusal
    /// returns no endpoint. Reservoir pressure acts in the accepted trial once.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_conductive_thermal_volume_momentum_with_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        conductivity: &PreparedPhaseConductivity<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        gate: &CancelGate,
    ) -> Result<ConductiveThermalVolumeMomentumStep, ChtError> {
        self.conduct_thermal_momentum(
            domain,
            faces,
            Some(profile),
            conductivity,
            input,
            controls,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn conduct_thermal_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        conductivity: &PreparedPhaseConductivity<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveThermalVolumeMomentumControls,
        gate: &CancelGate,
    ) -> Result<ConductiveThermalVolumeMomentumStep, ChtError> {
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
            |evolution| (&evolution.thermal, 0.0),
        )?;
        Ok(ConductiveThermalVolumeMomentumStep {
            evolution: accepted.evolution,
            conductive_enthalpy_increment_j: accepted.increments,
            report: accepted.report,
        })
    }
}
